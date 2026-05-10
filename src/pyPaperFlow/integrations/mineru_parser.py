"""
MinerU Content Parser for pyPaperFlow

Parse mineru >= 3.0 ``content_list_v2.json`` into canonical sectioned JSON
with metadata extraction and section aggregation.

Usage::

    from pyPaperFlow.integrations.mineru_parser import MinerUContentParser
    parser = MinerUContentParser()
    result = parser.parse("path/to/content_list_v2.json")

For more details, please refer to the official documentation https://opendatalab.github.io/MinerU/zh/reference/output_files/
"""

import json
import re
from pathlib import Path
from typing import *

# ── Section normalisation ──────────────────────────────────────

CANONICAL_ORDER = [
    "abstract", "introduction", "results", "discussion",
    "methods", "conclusion", "supplementary", "availability",
    "funding", "acknowledgements", "author_contributions",
    "references", "other",
]


# here we define regex patterns for matching section titles to canonical types
# you can refer to 🌟 pubmed/pubmed_merger.py 🌟 for similar pattern definitions 
'''
SECTION_CANONICAL_ORDER = [
    'abstract',
    'introduction',
    'results',
    'discussion',
    'methods',
    'conclusion',
    'supplementary',
    'availability',
    'funding',
    'acknowledgements',
    'author_contributions',
    'other',
]
'''

_SECTION_PATTERNS: dict[str, list[re.Pattern]] = {
    # summary = bioRxiv abstract-style heading; explicit "abstract" heading also covered
    "abstract":       [re.compile(r"^\s*abstract\s*$", re.I),
                       re.compile(r"^\s*summary\s*$", re.I)],
    "introduction":   [re.compile(r"^\s*introduction\s*$", re.I),
                       re.compile(r"^\s*intro(?:duction)?\s*$", re.I),
                       re.compile(r"^\s*background\s*$", re.I)],
    "results":        [re.compile(r"^\s*results?(?:\s+and\s+discussion)?\s*$", re.I),
                       re.compile(r"^\s*findings?\s*$", re.I)],
    "discussion":     [re.compile(r"^\s*discussions?\s*$", re.I)],
    "methods":        [re.compile(r"^\s*methods?\s*$", re.I),
                       re.compile(r"^\s*materials?\s*(?:and|&)\s*methods?\s*$", re.I),
                       re.compile(r"^\s*methodology\s*$", re.I),
                       re.compile(r"^\s*experimental\s+(?:section|procedures?|methods?)\s*$", re.I),
                       re.compile(r"^\s*computational\s+(?:details|methods?)\s*$", re.I),
                       re.compile(r"^\s*model\s+(?:architecture|design|formulation)\s*$", re.I),
                       re.compile(r"^\s*star\s+methods\s*$", re.I)],
    "conclusion":     [re.compile(r"^\s*conclusions?\s*$", re.I),
                       re.compile(r"^\s*concluding\s+remarks?\s*$", re.I)],
    "supplementary":  [re.compile(r"^\s*supplementary\s+(?:material|information|data|note)\s*$", re.I),
                       re.compile(r"^\s*supporting\s+information\s*$", re.I)],
    "availability":   [re.compile(r"^\s*(?:data|code|software)(?:\s+and\s+(?:data|code|software))*\s+availability\s*(?:statement)?\s*$", re.I),
                       re.compile(r"^\s*availability(?:\s+and\s+implementation)?\s*$", re.I)],
    "funding":        [re.compile(r"^\s*funding\s*$", re.I),
                       re.compile(r"^\s*financial\s+support\s*$", re.I)],
    # "acknowledg?e?ments?" covers: acknowledgment, acknowledgments, acknowledgement, acknowledgements
    "acknowledgements":[re.compile(r"^\s*acknowledg?e?ments?\s*$", re.I)],
    "author_contributions": [re.compile(r"^\s*authors?'?\s*contributions?\s*$", re.I),
                             re.compile(r"^\s*author\s+contribution\s*$", re.I)],
    "references":     [re.compile(r"^\s*references?\s*$", re.I),
                       re.compile(r"^\s*bibliography\s*$", re.I)],
}

_LEADING_NUM_RE = re.compile(r"^\s*(?:\d+[\.\)]\s*)+(.*)$")

_NOISE_TYPES = frozenset({
    "page_header", "page_footer", "page_number",
    "page_aside_text", "page_footnote",
})

_YEAR_RE = re.compile(r"\b(19|20)\d{2}\b")
_DOI_RE = re.compile(r"\b10\.\d{4,}/[^\s]+\b")
_EMAIL_RE = re.compile(r"E-mail:\s*\S+@\S+", re.I)


def _strip_section_number(title: str) -> str:
    m = _LEADING_NUM_RE.match(title)
    return m.group(1) if m else title


def normalize_section_title(raw_title: str) -> tuple[str, str]:
    """(canonical_type, display_title)"""
    stripped = _strip_section_number(raw_title)
    for canonical in CANONICAL_ORDER:
        if canonical == "other":
            continue
        for pat in _SECTION_PATTERNS.get(canonical, []):
            if pat.search(stripped):
                return canonical, canonical.replace("_", " ").title()
    return "other", raw_title.strip()


# ── Parser ─────────────────────────────────────────────────────

class MinerUContentParser:

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
        tables = self._extract_tables(pages)

        return {
            "source": "mineru",
            "file": json_path.name,
            "metadata": metadata,
            "abstract": abstract,
            "sections": sections,
            "tables": tables,
        }

    # ── flatten / collect ─────────────────────────────────

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
        """Collect all blocks including noise for metadata scanning."""
        out: list[dict] = []
        for pi, page in enumerate(pages):
            for block in page:
                if isinstance(block, dict):
                    block["_page"] = pi
                    out.append(block)
        return out

    # ── text extraction ───────────────────────────────────

    def _extract_text(self, block: dict) -> str | None:
        btype = block.get("type", "")
        content = block.get("content", {})

        if btype == "title":
            items = content.get("title_content") or []
            return " ".join(i.get("content", "") for i in items if i.get("type") == "text").strip() or None

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
                t = self._extract_text({"type": "paragraph", "content": {"paragraph_content": [it]}})
                if t:
                    lines.append(t)
            return "\n".join(lines) or None

        return None

    def _extract_text_raw(self, block: dict) -> str:
        t = self._extract_text(block)
        if t:
            return t
        content = block.get("content", {})
        for key in ("page_header_content", "page_footer_content", "page_footnote_content"):
            items = content.get(key) or []
            text = " ".join(i.get("content", "") for i in items if i.get("type") == "text").strip()
            if text:
                return text
        return ""

    def _is_title(self, block: dict) -> bool:
        return block.get("type") == "title"

    def _title_level(self, block: dict) -> int:
        return (block.get("content") or {}).get("level", 2)

    # ── metadata ──────────────────────────────────────────

    def _extract_metadata(self, raw_blocks: list[dict], blocks: list[dict]) -> dict:
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
        """Locate the author paragraph after the paper title.

        Handles both short author lines (e.g. IDPFold) and very long author lists
        (e.g. AlphaFold3 with 30+ co-authors and superscript affiliation markers).
        Scans paragraphs after the title until a level-2 heading appears, then picks
        the FIRST paragraph with >= 2 commas (preferring proximity over count).
        """
        found_title = False
        for b in blocks:
            if self._is_title(b) and self._title_level(b) == 1:
                found_title = True
                continue
            if not found_title:
                continue
            # Stop at level-2 headings (we're past the title+author+abstract region)
            if self._is_title(b) and self._title_level(b) >= 2:
                break
            if b.get("type") != "paragraph":
                continue
            text = self._extract_text(b)
            if not text:
                continue
            # Skip metadata paragraphs
            if re.search(r"^https?://", text.strip()):
                continue
            if re.search(r"^(?:Received|Accepted|Published|Open access|Check for)", text.strip(), re.I):
                continue
            if _EMAIL_RE.search(text) and text.count(",") < 3:
                continue
            # Keep the FIRST paragraph with >= 2 commas
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
        """Extract publication year preferring footer/journal-citation lines.

        Page footers often contain "Nature | Vol 630 | 13 June 2024" which
        gives the correct publication year, whereas random paragraphs may
        contain misleading years (e.g. "Received: 19 December 2023").
        """
        # First pass: only page_footer blocks (journal citation lines)
        for b in blocks:
            if b.get("type") != "page_footer":
                continue
            text = self._extract_text_raw(b)
            m = _YEAR_RE.search(text)
            if m:
                y = int(m.group())
                if 1900 <= y <= 2100:
                    return y
        # Fallback: any block
        for b in blocks:
            text = self._extract_text_raw(b)
            m = _YEAR_RE.search(text)
            if m:
                y = int(m.group())
                if 1900 <= y <= 2100:
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
        """Extract journal name from page headers, falling back to page footers.

        Falls back to detecting bioRxiv / medRxiv preprint headers when no
        standard journal header is found.  For Nature-family journals the header
        often just says "Article", but the footer carries the full citation line
        ("Nature | Vol 630 | 13 June 2024").
        """
        candidates: list[tuple[int, str]] = []

        def _scan(b: dict) -> None:
            text = self._extract_text_raw(b)
            if not text:
                return
            if "bioRxiv" in text or "medRxiv" in text:
                m = re.search(r"(?:bioRxiv|medRxiv)", text, re.I)
                if m:
                    # return from outer function — early exit
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
        # Fallback: page_footer blocks (e.g. Nature-style "Journal | Vol | Date")
        for b in blocks:
            if b.get("type") == "page_footer":
                text = self._extract_text_raw(b)
                if not text:
                    continue
                # "Nature | Vol 630 | 13 June 2024" → extract "Nature"
                m = re.match(r"^([A-Z][a-z]+(?:\s+[A-Z][a-z]+)*)\s*\|", text)
                if m:
                    return m.group(1)

        if candidates:
            candidates.sort(key=lambda x: (-(x[0] / max(len(x[1]), 1)), len(x[1])))
            return candidates[0][1]
        return ""

    # ── abstract ──────────────────────────────────────────

    def _extract_abstract(self, blocks: list[dict]) -> str:
        """Collect paragraphs between the paper title and the first section heading.

        Skips metadata paragraphs (DOIs, dates, affiliation lines, email blocks,
        corresponding-author notes, header text) that MinerU often places after the
        title on the first page.
        """
        _META_SKIP_RES = [
            re.compile(r"^https?://", re.I),
            re.compile(r"^(?:Received|Accepted|Published)", re.I),
            re.compile(r"^(?:Open access|Check for)", re.I),
            re.compile(r"^\*?\s*(?:Corresponding|Lead|Co-)", re.I),
            re.compile(r"^Email:", re.I),
            re.compile(r"^\d+\s", re.I),
            re.compile(r"^bioRxiv\s+preprint", re.I),
            re.compile(r"^\u00a9\s", re.I),
            re.compile(r"^\+\s", re.I),            # "+ Lead contact"
            re.compile(r"^Shruthi|^shruthiv", re.I),  # specific disbind noise
        ]

        paras: list[str] = []
        hit_abstract_heading = False
        authors_parsed = False

        for b in blocks:
            if self._is_title(b):
                raw = self._extract_text(b) or ""
                canonical, _ = normalize_section_title(raw)
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

            # Phase 1: skip past the author paragraph (comma-rich)
            if not authors_parsed:
                if text.count(",") >= 2:
                    authors_parsed = True
                continue

            # Phase 2: skip known metadata / boilerplate lines
            if any(pat.search(text) for pat in _META_SKIP_RES):
                continue

            # Phase 3: collect real abstract paragraphs
            paras.append(text)
            if len(" ".join(paras)) > 3000:
                break

        return " ".join(paras) if paras else ""

    # ── section building ──────────────────────────────────

    def _build_sections(self, blocks: list[dict]) -> list[dict]:
        sections: list[dict] = []
        current: dict | None = None       # current main section
        sub_current: dict | None = None   # current subsection (within current)
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

        for b in blocks:
            if self._is_title(b):
                raw = self._extract_text(b) or ""
                level = self._title_level(b)

                # Skip paper title (level 1) — it goes to metadata, not sections
                if level <= 1:
                    continue

                # Subsection: level >= 3, or multi-part number like "2.1."
                is_sub = (
                    current is not None
                    and (level >= 3
                         or bool(re.match(r"^\s*(?:\d+[\.\)]\s*){2,}", raw)))
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

                # New main section — flush previous
                _flush_section()
                pending.clear()

                canonical, display = normalize_section_title(raw)
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

        if not sections and pending:
            sections.append({
                "canonical_type": "other",
                "raw_title": "Other",
                "display_title": "Other",
                "level": 0,
                "paragraphs": pending,
            })

        return sections

    # ── section aggregation ───────────────────────────────

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

        return merged

    # ── table extraction ──────────────────────────────────

    def _extract_tables(self, pages: list[list[dict]]) -> list[dict]:
        tables: list[dict] = []
        for page in pages:
            for block in page:
                if not isinstance(block, dict) or block.get("type") != "table":
                    continue
                content = block.get("content", {})
                caption_items = content.get("table_caption") or []
                caption = " ".join(
                    i.get("content", "") for i in caption_items if i.get("type") == "text"
                ).strip()
                html = content.get("html", "")
                tables.append({"caption": caption, "html": html})
        return tables

    def _empty(self, json_path: Path) -> dict:
        return {
            "source": "mineru",
            "file": json_path.name,
            "metadata": {},
            "abstract": "",
            "sections": [],
            "tables": [],
        }
