"""
MinerU Content Parser for pyPaperFlow

Parse mineru >= 3.0 ``content_list_v2.json`` into canonical sectioned JSON
with metadata extraction and section aggregation.

Two classification backends are supported:

- **regex** (default): improved regex patterns with cursor-based ordering
  and context keyword scanning.  No API calls needed.
- **ai**: batch-classifies section titles via Claude / GPT API.

Usage::

    from pyPaperFlow.integrations.mineru_parser import (
        MinerUContentParser, RegexSectionClassifier, AISectionClassifier,
    )
    classifier = RegexSectionClassifier.from_config()
    parser = MinerUContentParser(classifier)
    result = parser.parse("path/to/content_list_v2.json")

For more details, please refer to the official documentation
https://opendatalab.github.io/MinerU/zh/reference/output_files/
"""
from __future__ import annotations
import json
import os
import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import *


#############################################################
#  1, Mineru Content Parser
#############################################################

# ── shared constants ────────────────────────────────────────────────

_NOISE_TYPES = frozenset({
    "page_header", "page_footer", "page_number",
    "page_aside_text", "page_footnote",
})

_YEAR_RE = re.compile(r"\b(19|20)\d{2}\b")
_DOI_RE = re.compile(r"\b10\.\d{4,}/[^\s]+\b")
_EMAIL_RE = re.compile(r"E-mail:\s*\S+@\S+", re.I)

_LEADING_NUM_RE = re.compile(r"^\s*(?:\d+[\.\)]\s*)+(.*)$")

_DEFAULT_CANONICAL_ORDER = [
    "abstract", "introduction", "results", "discussion",
    "methods", "conclusion", "supplementary", "availability",
    "funding", "acknowledgements", "author_contributions",
    "keywords", "conflicts", "references", "other",
]

_DEFAULT_DISPLAY_NAMES: dict[str, str] = {
    "abstract": "Abstract",
    "introduction": "Introduction",
    "results": "Results",
    "discussion": "Discussion",
    "methods": "Methods",
    "conclusion": "Conclusion",
    "supplementary": "Supplementary Material",
    "availability": "Data & Code Availability",
    "funding": "Funding",
    "acknowledgements": "Acknowledgements",
    "author_contributions": "Author Contributions",
    "keywords": "Keywords",
    "conflicts": "Competing Interests",
    "references": "References",
    "other": "Other",
}

def _strip_section_number(title: str) -> str:
    m = _LEADING_NUM_RE.match(title)
    return m.group(1) if m else title


# ── config loader ───────────────────────────────────────────────────

_DEFAULT_CONFIG_PATH = Path(__file__).parent / "mineru_config.yaml"


def _load_yaml_config(path: str | Path | None = None) -> dict:
    """Load a YAML config file.

    If *path* is ``None``, tries the default ``mineru_config.yaml`` shipped
    alongside this module.  Returns ``{}`` if YAML is not installed or the
    file is missing / unreadable.
    """
    try:
        import yaml
    except ImportError:
        return {}

    p = Path(path) if path else _DEFAULT_CONFIG_PATH
    if not p.exists():
        return {}
    try:
        with open(p, "r") as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


# ── SectionClassifier ───────────────────────────────────────────────

class SectionClassifier(ABC):
    """Base class for pluggable section-title classifiers."""

    def __init__(
        self,
        canonical_order: list[str] | None = None,
        display_names: dict[str, str] | None = None,
    ):
        self.canonical_order = canonical_order or list(_DEFAULT_CANONICAL_ORDER)
        self.display_names = display_names or dict(_DEFAULT_DISPLAY_NAMES)

    @abstractmethod
    def normalize(self, raw_title: str, context_text: str = "") -> tuple[str, str]:
        """Return (canonical_type, display_title)."""
        ...

    def batch_normalize(
        self, titles: list[dict[str, Any]]
    ) -> list[tuple[str, str]]:
        """Classify multiple titles. Default impl calls normalize() one-by-one."""
        return [self.normalize(t["raw"], t.get("context", "")) for t in titles]

    def display_for(self, canonical_type: str) -> str:
        return self.display_names.get(
            canonical_type,
            canonical_type.replace("_", " ").title(),
        )

    def order_key(self, canonical_type: str) -> int:
        try:
            return self.canonical_order.index(canonical_type)
        except ValueError:
            return len(self.canonical_order)


# ── alias builders (patterns come from YAML config only) ─────────────

def _build_strong_aliases(
    canonical_order: list[str], alias_cfg: dict
) -> dict[str, set[str]]:
    """Build strong (exact) aliases from config.  Config is the single source."""
    result: dict[str, set[str]] = {}
    for ctype in canonical_order:
        entry = alias_cfg.get(ctype, {})
        strong = {s.lower() for s in (entry.get("strong") or [])}
        if strong:
            result[ctype] = strong
    return result


def _build_weak_patterns(
    canonical_order: list[str], alias_cfg: dict
) -> dict[str, list[re.Pattern]]:
    """Build weak (regex) patterns from config.  Config is the single source."""
    result: dict[str, list[re.Pattern]] = {}
    for ctype in canonical_order:
        entry = alias_cfg.get(ctype, {})
        pat_strs = entry.get("weak") or []
        if pat_strs:
            result[ctype] = [re.compile(p, re.I) for p in pat_strs]
    return result


def _build_context_keywords(
    canonical_order: list[str], alias_cfg: dict
) -> dict[str, list[re.Pattern]]:
    """Build context keyword patterns from config.  Config is the single source."""
    result: dict[str, list[re.Pattern]] = {}
    for ctype in canonical_order:
        entry = alias_cfg.get(ctype, {})
        pat_strs = entry.get("context_keywords") or []
        if pat_strs:
            result[ctype] = [re.compile(p, re.I) for p in pat_strs]
    return result


# ── RegexSectionClassifier ──────────────────────────────────────────

class RegexSectionClassifier(SectionClassifier):
    """Classify section titles using configurable alias patterns + context.

    Match order:
    1. **strong** — exact lowercase string match (fast, best precision)
    2. **weak** — regex ``re.search`` against stripped title
    3. **context_keywords** — regex against first ~300 chars of body text
       following the title (fallback for ambiguous one-word titles like
       ``"Overview"``)

    A sliding cursor tracks the expected next position in ``canonical_order``.
    When a title matches multiple canonical types, the one closest to the
    cursor (in document order) wins.
    """

    def __init__(
        self,
        canonical_order: list[str] | None = None,
        display_names: dict[str, str] | None = None,
        strong_aliases: dict[str, set[str]] | None = None,
        weak_patterns: dict[str, list[re.Pattern]] | None = None,
        context_keywords: dict[str, list[re.Pattern]] | None = None,
    ):
        super().__init__(canonical_order, display_names)
        self._strong = strong_aliases or {}
        self._weak = weak_patterns or {}
        self._ctx_kw = context_keywords or {}
        self._cursor = 0  # index into canonical_order

    def reset_cursor(self) -> None:
        self._cursor = 0

    @classmethod
    def from_config(cls, config_path: str | Path | None = None) -> RegexSectionClassifier:
        """Build from a YAML config file, falling back to built-in defaults."""
        cfg = _load_yaml_config(config_path)

        canonical_order = cfg.get("canonical_order") or list(_DEFAULT_CANONICAL_ORDER)
        display_names = cfg.get("display_names") or dict(_DEFAULT_DISPLAY_NAMES)
        alias_cfg = cfg.get("aliases") or {}

        strong_aliases: dict[str, set[str]] = _build_strong_aliases(canonical_order, alias_cfg)
        weak_patterns: dict[str, list[re.Pattern]] = _build_weak_patterns(canonical_order, alias_cfg)
        context_keywords: dict[str, list[re.Pattern]] = _build_context_keywords(canonical_order, alias_cfg)

        return cls(
            canonical_order=canonical_order,
            display_names=display_names,
            strong_aliases=strong_aliases,
            weak_patterns=weak_patterns,
            context_keywords=context_keywords,
        )

    # ── public API ────────────────────────────────────────────

    def normalize(self, raw_title: str, context_text: str = "") -> tuple[str, str]:
        stripped = _strip_section_number(raw_title).strip()
        lower = stripped.lower().strip(" .:-—–\t\n\r")

        if not lower:
            return "other", raw_title.strip()

        # 1. strong match
        ctype = self._strong_match(lower)
        if ctype:
            self._advance_cursor(ctype)
            return ctype, self.display_for(ctype)

        # 2. weak regex match
        ctype = self._weak_match(lower)
        if ctype:
            self._advance_cursor(ctype)
            return ctype, self.display_for(ctype)

        # 3. context keyword fallback
        if context_text:
            ctype = self._context_match(context_text)
            if ctype:
                self._advance_cursor(ctype)
                return ctype, self.display_for(ctype)

        return "other", raw_title.strip()

    # ── matching helpers ──────────────────────────────────────

    def _strong_match(self, lower: str) -> str | None:
        candidates: list[tuple[int, str]] = []
        for ctype, aliases in self._strong.items():
            if lower == ctype or lower in aliases:
                candidates.append((self.order_key(ctype), ctype))
        return self._pick_best(candidates)

    def _weak_match(self, lower: str) -> str | None:
        candidates: list[tuple[int, str]] = []
        for ctype, patterns in self._weak.items():
            for pat in patterns:
                if pat.search(lower):
                    candidates.append((self.order_key(ctype), ctype))
                    break
        return self._pick_best(candidates)

    def _context_match(self, context_text: str) -> str | None:
        if not context_text:
            return None
        candidates: list[tuple[int, str]] = []
        for ctype, patterns in self._ctx_kw.items():
            for pat in patterns:
                if pat.search(context_text):
                    candidates.append((self.order_key(ctype), ctype))
                    break
        return self._pick_best(candidates)

    def _pick_best(self, candidates: list[tuple[int, str]]) -> str | None:
        """Pick the candidate closest to the cursor (lowest positive distance)."""
        if not candidates:
            return None
        if len(candidates) == 1:
            return candidates[0][1]
        # prefer the candidate whose canonical index is >= cursor and closest
        best = min(
            candidates,
            key=lambda x: (
                0 if x[0] >= self._cursor else 1,
                abs(x[0] - self._cursor),
            ),
        )
        return best[1]

    def _advance_cursor(self, ctype: str) -> None:
        idx = self.order_key(ctype)
        if idx >= self._cursor:
            self._cursor = idx + 1


# ── AISectionClassifier ─────────────────────────────────────────────

_AI_SYSTEM_PROMPT = """\
You are a scientific document section classifier. Given a list of section \
titles from an academic paper (biomedical / computational biology), classify \
each title into exactly one canonical type.

Canonical types:
- abstract: paper abstract or summary
- introduction: introduction, background, motivation
- results: results, findings, experimental outcomes
- discussion: discussion, interpretation, conclusion and future directions
- methods: methods, materials and methods, experimental section, \
  computational details, model architecture, training regime, algorithm details
- conclusion: conclusion, concluding remarks (standalone, not combined with discussion)
- supplementary: supplementary material, supporting information, \
  reporting summary, supplemental figures/tables
- availability: data availability, code availability, software availability
- funding: funding, financial support, grant information
- acknowledgements: acknowledgements, thanks
- author_contributions: author contributions, author statement
- keywords: keywords or key words
- conflicts: competing interests, conflict of interest, declaration of interests
- references: references, bibliography
- other: anything that does not fit the above (notes, article masthead, \
  online content, key resources table, diversity statement, etc.)

Return a JSON object:
{"classifications": [{"index": <int>, "canonical_type": "<type>"}, ...]}
Only return JSON, no other text."""

_AI_USER_PROMPT_TEMPLATE = """\
Classify the following section titles from a scientific paper.

Titles:
{titles_json}"""


class AISectionClassifier(SectionClassifier):
    """Classify section titles using an LLM API.

    Two modes, selected automatically:

    - **base_url set** → OpenAI chat-completions format against that endpoint.
      Covers DeepSeek, university proxies, self-hosted vLLM, etc.
    - **base_url not set** → Native Anthropic SDK (ANTHROPIC_API_KEY) or
      native OpenAI SDK (OPENAI_API_KEY), tried in that order.

    Sends all section titles from a paper in a single batch API call.
    """

    def __init__(
        self,
        canonical_order: list[str] | None = None,
        display_names: dict[str, str] | None = None,
        model: str = "claude-haiku-4-5",
        api_key: str | None = None,
        base_url: str | None = None,
    ):
        super().__init__(canonical_order, display_names)
        self.model = model
        self.api_key = api_key
        self.base_url = base_url

    def _get_api_key(self) -> str:
        if self.api_key:
            return self.api_key
        if self.base_url:
            return os.environ.get("OPENAI_API_KEY", "")
        return os.environ.get("ANTHROPIC_API_KEY", "") or os.environ.get("OPENAI_API_KEY", "")

    @classmethod
    def from_config(cls, config_path: str | Path | None = None) -> AISectionClassifier:
        cfg = _load_yaml_config(config_path)
        ai_cfg = cfg.get("ai", {}) if cfg else {}

        model = ai_cfg.get("model", "claude-haiku-4-5")
        api_key = ai_cfg.get("api_key") or None
        base_url = ai_cfg.get("base_url") or None

        canonical_order = cfg.get("canonical_order") or list(_DEFAULT_CANONICAL_ORDER)
        display_names = cfg.get("display_names") or dict(_DEFAULT_DISPLAY_NAMES)

        return cls(
            canonical_order=canonical_order,
            display_names=display_names,
            model=model,
            api_key=api_key,
            base_url=base_url,
        )

    def normalize(self, raw_title: str, context_text: str = "") -> tuple[str, str]:
        # Single-title fallback: delegate to batch
        results = self.batch_normalize([{"raw": raw_title, "context": context_text}])
        return results[0]

    def batch_normalize(
        self, titles: list[dict[str, Any]]
    ) -> list[tuple[str, str]]:
        if not titles:
            return []

        # Build the JSON payload
        items = []
        for i, t in enumerate(titles):
            item: dict = {"index": i, "title": t["raw"].strip()}
            ctx = t.get("context", "")
            if ctx:
                item["context_preview"] = ctx[:200]
            items.append(item)

        titles_json = json.dumps(items, ensure_ascii=False, indent=2)
        user_prompt = _AI_USER_PROMPT_TEMPLATE.format(titles_json=titles_json)

        try:
            raw_response = self._call_api(user_prompt)
            classifications = self._parse_response(raw_response)
        except Exception:
            # Fall back to "other" for all titles on any API failure
            return [("other", t["raw"].strip()) for t in titles]

        # Build result list preserving input order
        classified: dict[int, str] = {}
        for c in classifications:
            idx = c.get("index", -1)
            ctype = c.get("canonical_type", "other")
            if ctype not in set(self.canonical_order):
                ctype = "other"
            classified[idx] = ctype

        results: list[tuple[str, str]] = []
        for i, t in enumerate(titles):
            ctype = classified.get(i, "other")
            results.append((ctype, self.display_for(ctype)))
        return results

    def _call_api(self, user_prompt: str) -> str:
        if self.base_url:
            return self._call_openai_compatible(user_prompt)
        # No base_url: try Anthropic first, fall back to OpenAI
        if os.environ.get("ANTHROPIC_API_KEY") or self.api_key:
            try:
                return self._call_anthropic(user_prompt)
            except Exception:
                pass
        return self._call_openai(user_prompt)

    def _call_openai_compatible(self, user_prompt: str) -> str:
        """OpenAI chat-completions against a custom base_url (DeepSeek, proxy, etc.)."""
        from openai import OpenAI

        api_key = self._get_api_key()
        if not api_key:
            raise ValueError(
                "API key not set. Set OPENAI_API_KEY env var or pass --api-key."
            )

        client = OpenAI(api_key=api_key, base_url=self.base_url)
        response = client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": _AI_SYSTEM_PROMPT},
                {"role": "user", "content": user_prompt},
            ],
        )
        return response.choices[0].message.content or ""

    def _call_openai(self, user_prompt: str) -> str:
        from openai import OpenAI

        api_key = self.api_key or os.environ.get("OPENAI_API_KEY", "")
        if not api_key:
            raise ValueError("OPENAI_API_KEY not set")

        client = OpenAI(api_key=api_key)
        response = client.chat.completions.create(
            model=self.model,
            messages=[
                {"role": "system", "content": _AI_SYSTEM_PROMPT},
                {"role": "user", "content": user_prompt},
            ],
        )
        return response.choices[0].message.content or ""

    def _call_anthropic(self, user_prompt: str) -> str:
        import anthropic

        api_key = self.api_key or os.environ.get("ANTHROPIC_API_KEY", "")
        if not api_key:
            raise ValueError("ANTHROPIC_API_KEY not set")

        client = anthropic.Anthropic(api_key=api_key)
        message = client.messages.create(
            model=self.model,
            max_tokens=500,
            system=_AI_SYSTEM_PROMPT,
            messages=[{"role": "user", "content": user_prompt}],
        )
        return message.content[0].text

    def _parse_response(self, raw: str) -> list[dict]:
        # Strip markdown code fences if present
        text = raw.strip()
        if text.startswith("```"):
            text = re.sub(r"^```(?:json)?\s*", "", text)
            text = re.sub(r"```\s*$", "", text)
        data = json.loads(text)
        return data.get("classifications", [])


# ── Parser ──────────────────────────────────────────────────────────

class MinerUContentParser:
    """Parse MinerU ``content_list_v2.json`` into canonical sectioned JSON.

    Parameters
    ----------
    classifier:
        A ``SectionClassifier`` instance.  If *None*, a default
        ``RegexSectionClassifier`` is created from the built-in config.
    """

    def __init__(self, classifier: SectionClassifier | None = None):
        self.classifier = classifier or RegexSectionClassifier.from_config()

    def parse(self, json_path: str | Path) -> dict[str, Any]:
        json_path = Path(json_path)
        with open(json_path, "r") as f:
            pages: list[list[dict]] = json.load(f)

        blocks = self._flatten(pages)
        if not blocks:
            return self._empty(json_path)

        raw = self._collect_raw(pages)
        metadata = self._extract_metadata(raw, blocks)
        abstract = self._extract_abstract(blocks)
        sections = self._build_sections(blocks)
        sections = self._aggregate_sections(sections)

        # Ensure abstract is always the first section
        if abstract and not any(s["canonical_type"] == "abstract" for s in sections):
            sections.insert(0, {
                "canonical_type": "abstract",
                "raw_title": "Abstract",
                "display_title": self.classifier.display_for("abstract"),
                "level": 2,
                "paragraphs": [abstract],
            })

        backend = (
            "ai" if isinstance(self.classifier, AISectionClassifier) else "regex"
        )

        return {
            "source": "mineru",
            "file": json_path.name,
            "backend": backend,
            "metadata": metadata,
            "sections": sections,
        }

    # ── flatten / collect ─────────────────────────────────────

    def _flatten(self, pages: list[list[dict]]) -> list[dict]:
        out: list[dict] = []
        for pi, page in enumerate(pages):
            for block in page:
                if not isinstance(block, dict):
                    continue
                if block.get("type", "") in _NOISE_TYPES:
                    continue
                block["_page"] = pi
                out.append(block)
        return out

    def _collect_raw(self, pages: list[list[dict]]) -> list[dict]:
        out: list[dict] = []
        for pi, page in enumerate(pages):
            for block in page:
                if isinstance(block, dict):
                    block["_page"] = pi
                    out.append(block)
        return out

    # ── text extraction ───────────────────────────────────────

    def _extract_text(self, block: dict) -> str | None:
        btype = block.get("type", "")
        content = block.get("content", {})

        if btype == "title":
            items = content.get("title_content") or []
            return " ".join(
                i.get("content", "") for i in items if i.get("type") == "text"
            ).strip() or None

        if btype == "paragraph":
            items = content.get("paragraph_content") or []
            texts: list[str] = []
            for i in items:
                if i.get("type") == "text":
                    t = i.get("content", "").strip()
                    if t:
                        texts.append(t)
                elif i.get("type") == "equation_inline":
                    eq = i.get("content", "").strip()
                    if eq:
                        texts.append(f"${eq}$")
            return " ".join(texts).strip() or None

        if btype == "equation_interline":
            math = content.get("math_content", "").strip()
            return f"$${math}$$" if math else None

        if btype == "code":
            return content.get("code_content", "") or None

        if btype in ("list", "index"):
            items = content.get("list_items") or []
            lines = []
            for it in items:
                t = self._extract_text(
                    {"type": "paragraph", "content": {"paragraph_content": [it]}}
                )
                if t:
                    lines.append(t)
            return "\n".join(lines) or None

        if btype in ("image", "chart"):
            caption = content.get("image_caption") or content.get("chart_caption") or []
            text = " ".join(
                i.get("content", "") for i in caption if i.get("type") == "text"
            ).strip()
            return f"[Figure: {text}]" if text else None

        return None

    def _extract_text_raw(self, block: dict) -> str:
        t = self._extract_text(block)
        if t:
            return t
        content = block.get("content", {})
        for key in (
            "page_header_content",
            "page_footer_content",
            "page_footnote_content",
            "page_aside_text_content",
        ):
            items = content.get(key) or []
            text = " ".join(
                i.get("content", "") for i in items if i.get("type") == "text"
            ).strip()
            if text:
                return text
        return ""

    def _is_title(self, block: dict) -> bool:
        return block.get("type") == "title"

    def _title_level(self, block: dict) -> int:
        return (block.get("content") or {}).get("level", 2)

    # ── context snippet (for classifier) ──────────────────────

    def _collect_context(
        self, blocks: list[dict], title_idx: int, max_chars: int = 300
    ) -> str:
        """Collect text from the 2-3 paragraphs following a title block."""
        snippets: list[str] = []
        total = 0
        for i in range(title_idx + 1, min(title_idx + 4, len(blocks))):
            b = blocks[i]
            if b.get("type") in ("paragraph",):
                t = self._extract_text(b)
                if t:
                    snippets.append(t)
                    total += len(t)
                    if total >= max_chars:
                        break
        return " ".join(snippets)

    # ── metadata ──────────────────────────────────────────────

    def _extract_metadata(
        self, raw_blocks: list[dict], blocks: list[dict]
    ) -> dict:
        return {
            "title": self._extract_title(blocks),
            "authors": self._extract_authors(blocks),
            "year": self._extract_year(raw_blocks),
            "doi": self._extract_doi(raw_blocks),
            "journal": self._extract_journal(raw_blocks),
        }

    def _extract_title(self, blocks: list[dict]) -> str:
        for b in blocks:
            if self._is_title(b) and self._title_level(b) == 1:
                return self._extract_text(b) or ""
        return ""

    def _extract_authors(self, blocks: list[dict]) -> str:
        found_title = False
        for b in blocks:
            if self._is_title(b) and self._title_level(b) == 1:
                found_title = True
                continue
            if not found_title:
                continue
            if self._is_title(b) and self._title_level(b) >= 2:
                break
            if b.get("type") != "paragraph":
                continue
            text = self._extract_text(b)
            if not text:
                continue
            if re.search(r"^https?://", text.strip()):
                continue
            if re.search(
                r"^(?:Received|Accepted|Published|Open access|Check for)",
                text.strip(),
                re.I,
            ):
                continue
            if _EMAIL_RE.search(text) and text.count(",") < 3:
                continue
            if text.count(",") >= 2:
                lines = text.split("  ")
                clean = []
                for line in lines:
                    if _EMAIL_RE.search(line):
                        continue
                    if re.match(r"^\d", line.strip()):
                        continue
                    clean.append(line.strip())
                return " ".join(clean).strip()
        return ""

    def _extract_year(self, blocks: list[dict]) -> int | None:
        # 1. page_footer — journal citation lines (most reliable)
        for b in blocks:
            if b.get("type") != "page_footer":
                continue
            text = self._extract_text_raw(b)
            m = _YEAR_RE.search(text)
            if m:
                y = int(m.group())
                if 1900 <= y <= 2100:
                    return y

        # 2. page_aside_text — arXiv date stamp
        #    e.g. "arXiv:2409.02240v1 [physics.bio-ph] 3 Sep 2024"
        #    Only consider if the block looks like an arXiv header (has "arXiv" or a month name).
        for b in blocks:
            if b.get("type") != "page_aside_text":
                continue
            text = self._extract_text_raw(b)
            if "arXiv" not in text and "arxiv" not in text:
                continue
            for part in text.split():
                m = _YEAR_RE.search(part)
                if m:
                    y = int(m.group())
                    if 2000 <= y <= 2100:
                        return y

        # 3. fallback — scan first 3 pages first (abstract/body), then the rest
        for b in blocks:
            if b.get("_page", 999) > 2:
                continue
            text = self._extract_text_raw(b)
            m = _YEAR_RE.search(text)
            if m:
                y = int(m.group())
                if 2000 <= y <= 2100:
                    return y
        for b in blocks:
            text = self._extract_text_raw(b)
            m = _YEAR_RE.search(text)
            if m:
                y = int(m.group())
                if 2000 <= y <= 2100:
                    return y
        return None

    def _extract_doi(self, blocks: list[dict]) -> str:
        for b in blocks:
            text = self._extract_text_raw(b)
            m = _DOI_RE.search(text)
            if m:
                return m.group().rstrip(".,;")
        return ""

    def _extract_journal(self, blocks: list[dict]) -> str:
        candidates: list[tuple[int, str]] = []

        def _scan(b: dict) -> None:
            text = self._extract_text_raw(b)
            if not text:
                return
            if "bioRxiv" in text or "medRxiv" in text:
                m = re.search(r"(?:bioRxiv|medRxiv)", text, re.I)
                if m:
                    nonlocal journal
                    journal = m.group()
            upper_chars = sum(1 for c in text if c.isupper())
            if upper_chars > len(text) * 0.5 and 2 <= len(text.split()) <= 5:
                candidates.append((upper_chars, text))

        journal = ""
        for b in blocks:
            if b.get("type") == "page_header":
                _scan(b)
                if journal:
                    return journal
        for b in blocks:
            if b.get("type") == "page_footer":
                text = self._extract_text_raw(b)
                if not text:
                    continue
                m = re.match(
                    r"^([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\s*\|", text
                )
                if m:
                    return m.group(1)
        if candidates:
            candidates.sort(
                key=lambda x: (-(x[0] / max(len(x[1]), 1)), len(x[1]))
            )
            return candidates[0][1]
        return ""

    # ── abstract ──────────────────────────────────────────────

    def _extract_abstract(self, blocks: list[dict]) -> str:
        _META_SKIP_RES = [
            re.compile(r"^https?://", re.I),
            re.compile(r"^(?:Received|Accepted|Published)", re.I),
            re.compile(r"^(?:Open access|Check for)", re.I),
            re.compile(r"^\*?\s*(?:Corresponding|Lead|Co-)", re.I),
            re.compile(r"^Email:", re.I),
            re.compile(r"^\d+\s", re.I),
            re.compile(r"^bioRxiv\s+preprint", re.I),
            re.compile(r"^©\s", re.I),
            re.compile(r"^\+\s", re.I),
            re.compile(r"^Shruthi|^shruthiv", re.I),
        ]

        paras: list[str] = []
        hit_abstract_heading = False
        authors_parsed = False

        for b in blocks:
            if self._is_title(b):
                raw = self._extract_text(b) or ""
                canonical, _ = self.classifier.normalize(raw)
                if canonical == "abstract":
                    hit_abstract_heading = True
                    continue
                if hit_abstract_heading or self._title_level(b) >= 2:
                    break
                continue

            if b.get("type") != "paragraph":
                continue

            text = self._extract_text(b)
            if not text:
                continue

            if not authors_parsed:
                if text.count(",") >= 2:
                    authors_parsed = True
                continue

            if any(pat.search(text) for pat in _META_SKIP_RES):
                continue

            paras.append(text)
            if len(" ".join(paras)) > 3000:
                break

        return " ".join(paras) if paras else ""

    # ── section building ──────────────────────────────────────

    def _build_sections(self, blocks: list[dict]) -> list[dict]:
        sections: list[dict] = []
        current: dict | None = None
        sub_current: dict | None = None
        pending: list[str] = []

        # Reset cursor for regex classifier before processing document
        if hasattr(self.classifier, "reset_cursor"):
            self.classifier.reset_cursor()

        # Collect title blocks with context for AI batch classification
        if isinstance(self.classifier, AISectionClassifier):
            sections = self._build_sections_ai(blocks)
        else:
            sections = self._build_sections_regex(blocks)

        if not sections and pending:
            sections.append({
                "canonical_type": "other",
                "raw_title": "Other",
                "display_title": "Other",
                "level": 0,
                "paragraphs": pending,
            })

        return sections

    def _build_sections_regex(self, blocks: list[dict]) -> list[dict]:
        """Build sections using the regex classifier (per-title calls)."""
        sections: list[dict] = []
        current: dict | None = None
        sub_current: dict | None = None
        pending: list[str] = []

        def _flush_section() -> None:
            nonlocal current, sub_current
            if current is not None and (
                current.get("paragraphs") or current.get("subsections") or pending
            ):
                if pending:
                    current.setdefault("paragraphs", []).extend(pending)
                sections.append(current)
            current = None
            sub_current = None

        for bi, b in enumerate(blocks):
            if self._is_title(b):
                raw = self._extract_text(b) or ""
                level = self._title_level(b)

                if level <= 1:
                    continue

                is_sub = (
                    current is not None
                    and (
                        level >= 3
                        or bool(re.match(r"^\s*(?:\d+[\.\)]\s*){2,}", raw))
                    )
                )

                if is_sub and current is not None:
                    subs = current.setdefault("subsections", [])
                    sub_current = {
                        "raw_title": raw,
                        "display_title": raw,
                        "level": level,
                        "paragraphs": [],
                    }
                    subs.append(sub_current)
                    continue

                _flush_section()
                pending.clear()

                ctx = self._collect_context(blocks, bi)
                canonical, display = self.classifier.normalize(raw, ctx)
                current = {
                    "canonical_type": canonical,
                    "raw_title": raw,
                    "display_title": display,
                    "level": level,
                    "paragraphs": [],
                }
                sub_current = None
                continue

            text = self._extract_text(b)
            if not text:
                continue

            if sub_current is not None:
                sub_current["paragraphs"].append(text)
            elif current is not None:
                current["paragraphs"].append(text)
            else:
                pending.append(text)

        if pending:
            pending.clear()
        _flush_section()
        return sections

    def _build_sections_ai(self, blocks: list[dict]) -> list[dict]:
        """Build sections using AI batch classification (one API call)."""
        # First pass: collect all title blocks
        title_infos: list[dict] = []  # {bi, raw, level}
        for bi, b in enumerate(blocks):
            if self._is_title(b):
                raw = self._extract_text(b) or ""
                level = self._title_level(b)
                if level >= 2:
                    ctx = self._collect_context(blocks, bi)
                    title_infos.append({
                        "index": bi,
                        "raw": raw,
                        "level": level,
                        "context": ctx,
                    })

        # Batch classify
        if title_infos:
            results = self.classifier.batch_normalize(title_infos)
            classified: dict[int, tuple[str, str]] = {
                title_infos[i]["index"]: results[i] for i in range(len(title_infos))
            }
        else:
            classified = {}

        # Second pass: build sections
        sections: list[dict] = []
        current: dict | None = None
        sub_current: dict | None = None
        pending: list[str] = []

        def _flush_section() -> None:
            nonlocal current, sub_current
            if current is not None and (
                current.get("paragraphs") or current.get("subsections") or pending
            ):
                if pending:
                    current.setdefault("paragraphs", []).extend(pending)
                sections.append(current)
            current = None
            sub_current = None

        for bi, b in enumerate(blocks):
            if self._is_title(b):
                raw = self._extract_text(b) or ""
                level = self._title_level(b)

                if level <= 1:
                    continue

                is_sub = (
                    current is not None
                    and (
                        level >= 3
                        or bool(re.match(r"^\s*(?:\d+[\.\)]\s*){2,}", raw))
                    )
                )

                if is_sub and current is not None:
                    subs = current.setdefault("subsections", [])
                    sub_current = {
                        "raw_title": raw,
                        "display_title": raw,
                        "level": level,
                        "paragraphs": [],
                    }
                    subs.append(sub_current)
                    continue

                _flush_section()
                pending.clear()

                canonical, display = classified.get(
                    bi, self.classifier.normalize(raw)
                )
                current = {
                    "canonical_type": canonical,
                    "raw_title": raw,
                    "display_title": display,
                    "level": level,
                    "paragraphs": [],
                }
                sub_current = None
                continue

            text = self._extract_text(b)
            if not text:
                continue

            if sub_current is not None:
                sub_current["paragraphs"].append(text)
            elif current is not None:
                current["paragraphs"].append(text)
            else:
                pending.append(text)

        if pending:
            pending.clear()
        _flush_section()
        return sections

    # ── section aggregation ───────────────────────────────────

    def _aggregate_sections(self, sections: list[dict]) -> list[dict]:
        merged: list[dict] = []
        seen: dict[str, dict] = {}

        for sec in sections:
            ctype = sec["canonical_type"]
            if ctype != "other" and ctype in seen:
                parent = seen[ctype]
                parent["paragraphs"].extend(sec.get("paragraphs", []))
                subs = parent.setdefault("subsections", [])
                subs.append({
                    "raw_title": sec["raw_title"],
                    "display_title": sec["display_title"],
                    "level": sec["level"],
                    "paragraphs": sec.get("paragraphs", []),
                })
            elif ctype == "other":
                merged.append(sec)
            else:
                merged.append(sec)
                seen[ctype] = sec

        # Sort by canonical order
        order_map = {
            t: i for i, t in enumerate(self.classifier.canonical_order)
        }
        merged.sort(
            key=lambda s: order_map.get(s["canonical_type"], len(order_map))
        )
        return merged

    def _empty(self, json_path: Path) -> dict:
        return {
            "source": "mineru",
            "file": json_path.name,
            "backend": (
                "ai"
                if isinstance(self.classifier, AISectionClassifier)
                else "regex"
            ),
            "metadata": {},
            "sections": [],
        }



#############################################################
#  2, Mineru Markdown Export
#############################################################


# ── Markdown export ───────────────────────────────────────────────────

def _slugify(text: str) -> str:
    """Create a URL-friendly anchor slug from heading text."""
    value = re.sub(r"[^a-z0-9\s-]", "", text.lower())
    value = re.sub(r"\s+", "-", value.strip())
    value = re.sub(r"-+", "-", value)
    return value or "paper"


def export_mineru_json_to_md(
    input_json: str | Path,
    output_md: str | Path,
    yaml_cfg: str | Path | None = None,
) -> dict:
    """Export a structured mineru JSON to a Markdown file for LLM consumption.

    Parameters
    ----------
    input_json:
        Path to a JSON file produced by ``MinerUContentParser.parse()``,
        or a directory containing multiple such files.
    output_md:
        Output Markdown file path.
    yaml_cfg:
        Optional YAML config specifying which sections to include.
        Example::

            content_sections:
              - abstract
              - introduction
              - methods
              - results
              - discussion

        If not provided, ALL sections are included.

    Returns
    -------
    A dict with keys: ``total``, ``output``, ``sections_exported``.
    """
    input_path = Path(input_json)
    output_path = Path(output_md)

    # ── Collect JSON files ──
    json_files: list[Path] = []
    if input_path.is_dir():
        json_files = sorted(input_path.glob("*.json"))
    elif input_path.is_file():
        json_files = [input_path]
    else:
        raise FileNotFoundError(f"Input not found: {input_json}")

    # ── Load config ──
    cfg = _load_yaml_config(yaml_cfg)
    requested_sections: list[str] | None = None
    if cfg:
        requested = cfg.get("content_sections")
        if isinstance(requested, list) and requested:
            requested_sections = [s.strip().lower() for s in requested]

    # ── Load papers ──
    papers: list[dict] = []
    for jf in json_files:
        try:
            with open(jf, "r") as fh:
                papers.append(json.load(fh))
        except Exception:
            continue

    if not papers:
        raise ValueError("No valid JSON files found")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    sections_exported: dict[str, int] = {}

    # ── Build index entries (title → anchor) ──
    index_entries: list[tuple[str, str]] = []
    for paper in papers:
        meta = paper.get("metadata", {})
        title = meta.get("title", "Untitled").strip()
        index_entries.append((title, _slugify(title)))

    with open(output_path, "w") as f:
        # ── Index ──
        if len(index_entries) > 1:
            f.write("# Index\n\n")
            for title, anchor in index_entries:
                f.write(f"- [{title}](#{anchor})\n")
            f.write("\n---\n\n")

        # ── Papers ──
        for i, paper in enumerate(papers):
            meta = paper.get("metadata", {})
            title = meta.get("title", "Untitled").strip()
            authors = meta.get("authors", "")
            year = meta.get("year", "")
            doi = meta.get("doi", "")
            journal = meta.get("journal", "")
            anchor = _slugify(title)

            # anchor + paper title
            f.write(f'<a id="{anchor}"></a>\n\n')
            f.write(f"# {title}\n\n")

            # metadata block
            if authors:
                f.write(f"**Authors:** {authors}\n\n")
            meta_parts = []
            if year:
                meta_parts.append(f"**Year:** {year}")
            if journal:
                meta_parts.append(f"**Journal:** {journal}")
            if doi:
                meta_parts.append(f"**DOI:** {doi}")
            if meta_parts:
                f.write(" | ".join(meta_parts) + "\n\n")

            # ── sections ──
            sections = paper.get("sections", [])
            for sec in sections:
                ctype = sec.get("canonical_type", "other")
                if requested_sections is not None and ctype not in requested_sections:
                    continue

                sections_exported[ctype] = sections_exported.get(ctype, 0) + 1

                display = sec.get("display_title", ctype)
                f.write(f"## {display}\n\n")

                for para in sec.get("paragraphs", []):
                    if para.strip():
                        f.write(f"{para}\n\n")

                for sub in sec.get("subsections", []):
                    sub_title = sub.get("display_title", "")
                    f.write(f"### {sub_title}\n\n")
                    for para in sub.get("paragraphs", []):
                        if para.strip():
                            f.write(f"{para}\n\n")

            # paper separator
            if i < len(papers) - 1:
                f.write("\n---\n\n")

    return {
        "total": len(papers),
        "output": str(output_path),
        "sections_exported": sections_exported,
    }
