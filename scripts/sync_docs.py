#!/usr/bin/env python3
"""Regenerate mkdoc_site/*.md content pages from README_zh.md.

Splits the Chinese README by its H2/H3 headings. ``RULES`` gives known
sections stable English filenames; any H2 not listed there is auto-named from
its heading text, and any H3 not listed is folded into its parent H2 page. So
adding or removing a top-level section in the README needs no edits here —
the page appears/disappears automatically, and ``nav.yml`` follows the same
structure. Cross-page links and local asset/config paths are rewritten, and
heading levels are demoted so each page has a single H1. ``figs/`` is copied
into ``mkdoc_site/assets/`` so local image links resolve.

The homepage ``mkdoc_site/index.md`` is hand-written and left untouched; it is
the only tracked file under ``mkdoc_site/`` (everything else is gitignored).
``mkdoc_site/nav.yml`` is generated here and merged into ``mkdocs.yml`` via
``INHERIT``, so the nav tracks the README structure automatically.

Run from the repo root (also run by the CI docs workflow before ``mkdocs build``):

    python scripts/sync_docs.py
"""
from __future__ import annotations

import re
import pathlib
import shutil

ROOT = pathlib.Path(__file__).resolve().parent.parent
README = ROOT / "README_zh.md"
DOCS = ROOT / "mkdoc_site"
GH = "https://github.com/MaybeBio/pyPaperFlow"

# (heading_level, keyword, target_relpath). ``keyword`` is a substring of the
# heading text; first match wins. H2s not listed are auto-named from their
# heading text; H3s not listed are folded into their parent H2 page.
RULES = [
    (2, "简介", "overview.md"),
    (2, "功能特性", "features.md"),
    (2, "架构", "architecture.md"),
    (2, "安装", "installation.md"),
    (2, "使用方法", "usage/index.md"),
    (3, "模块概述", "usage/index.md"),
    (3, "研究起点", "usage/research-start.md"),
    (3, "文献检索", "usage/search.md"),
    (3, "文献获取", "usage/fetch.md"),
    (3, "内容提取", "usage/extraction.md"),
    (3, "其他文献数据平台", "usage/other-databases.md"),
    (3, "批判性阅读", "usage/downstream.md"),
    (3, "Reading与Coding", "usage/reading-coding.md"),
    (2, "组合利用方案", "combination.md"),
    (2, "测试示例", "test-cases.md"),
    (2, "文献调研示例", "full-survey.md"),
    (2, "维护待办", "todo.md"),
]

HEADING_RE = re.compile(r"^(#{1,6})\s+(.*)$")


def find_rule(level: int, text: str) -> str | None:
    for rule_level, keyword, target in RULES:
        if rule_level == level and keyword in text:
            return target
    return None


def slugify(text: str) -> str:
    """Heading text -> filename stem (keeps CJK/alnum, maps the rest to '-')."""
    stem = text.strip()
    stem = re.sub(r"[*_`]", "", stem)
    stem = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", stem)
    stem = re.sub(r"[^\w]+", "-", stem)
    return stem.strip("-") or "section"


def gh_url(path: str) -> str:
    """GitHub URL for a repo path: tree/ for directories, blob/ for files."""
    kind = "tree" if path.endswith("/") else "blob"
    return f"{GH}/{kind}/main/{path}"


def process(content: str, prefix: str) -> str:
    # 1) repo-file links (docs/config/test) -> GitHub URLs
    for top in ("docs", "config", "test"):
        content = re.sub(
            rf"\]\(\./{top}/([^)]*)\)",
            lambda m, t=top: f"]({gh_url(f'{t}/{m.group(1)}')})",
            content,
        )
    # 2) local figs -> assets
    content = re.sub(
        r"\]\(\./figs/([^)]+)\)",
        lambda m: f"]({prefix}assets/{m.group(1)})",
        content,
    )
    # 3) demote headings so the page's first heading becomes H1 (skip fences)
    in_fence = False
    first_level: int | None = None
    out: list[str] = []
    for line in content.splitlines(keepends=True):
        stripped = line.lstrip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_fence = not in_fence
            out.append(line)
            continue
        m = re.match(r"^(#{1,6})\s", line)
        if m and not in_fence:
            level = len(m.group(1))
            if first_level is None:
                first_level = level
            new_level = max(1, level - (first_level - 1))
            out.append("#" * new_level + line[level:])
        else:
            out.append(line)
    return "".join(out)


def parse_structure(lines: list[str]) -> list[dict]:
    """Walk the README into an ordered H2 tree (each node records its line index)."""
    sections: list[dict] = []
    current: dict | None = None
    in_fence = False
    for idx, line in enumerate(lines):
        stripped = line.lstrip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        m = HEADING_RE.match(line)
        if not m:
            continue
        level = len(m.group(1))
        title = m.group(2).strip()
        if level == 2:
            if title == "目录":
                current = None
                continue
            current = {"level": 2, "title": title, "idx": idx, "children": []}
            sections.append(current)
        elif level == 3 and current is not None:
            current["children"].append({"level": 3, "title": title, "idx": idx})
    return sections


def assign_targets(sections: list[dict]) -> None:
    for sec in sections:
        rule2 = find_rule(2, sec["title"])
        if rule2:
            sec["target"] = rule2
        elif sec["children"]:
            sec["target"] = f"{slugify(sec['title'])}/index.md"
        else:
            sec["target"] = f"{slugify(sec['title'])}.md"

        for child in sec["children"]:
            rule3 = find_rule(3, child["title"])
            child["target"] = rule3 if rule3 else sec["target"]


def build_buffers(lines: list[str], sections: list[dict]) -> tuple[dict[str, list[str]], list[str]]:
    boundaries: list[tuple[int, str]] = []
    for sec in sections:
        boundaries.append((sec["idx"], sec["target"]))
        for child in sec["children"]:
            boundaries.append((child["idx"], child["target"]))
    boundaries.sort()

    buffers: dict[str, list[str]] = {}
    order: list[str] = []
    for i, (start, target) in enumerate(boundaries):
        end = boundaries[i + 1][0] if i + 1 < len(boundaries) else len(lines)
        if target not in buffers:
            buffers[target] = []
            order.append(target)
        buffers[target].extend(lines[start:end])
    return buffers, order


def build_nav(sections: list[dict]) -> str:
    def q(text: str) -> str:
        return f'"{text}"'

    lines = ["nav:", "  - 首页: index.md"]
    for sec in sections:
        subs = [c for c in sec["children"] if c["target"] != sec["target"]]
        if subs:
            lines.append(f"  - {q(sec['title'])}:")
            for c in sec["children"]:
                lines.append(f"      - {q(c['title'])}: {c['target']}")
        else:
            lines.append(f"  - {q(sec['title'])}: {sec['target']}")
    return "\n".join(lines) + "\n"


def clean_docs() -> None:
    for p in DOCS.iterdir():
        if p.name == "index.md":
            continue
        if p.is_dir():
            shutil.rmtree(p, ignore_errors=True)
        else:
            p.unlink(missing_ok=True)


def main() -> None:
    clean_docs()
    shutil.copytree(ROOT / "figs", DOCS / "assets")

    lines = README.read_text(encoding="utf-8").splitlines(keepends=True)
    sections = parse_structure(lines)
    assign_targets(sections)
    buffers, order = build_buffers(lines, sections)

    for target in order:
        body = "".join(buffers[target])
        rendered = process(body, "../" * target.count("/"))
        out = DOCS / target
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(rendered, encoding="utf-8")
        print(f"wrote {target}")

    (DOCS / "nav.yml").write_text(build_nav(sections), encoding="utf-8")
    print("wrote nav.yml")
    print("done")


if __name__ == "__main__":
    main()
