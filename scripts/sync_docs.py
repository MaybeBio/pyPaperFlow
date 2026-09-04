#!/usr/bin/env python3
"""Regenerate mkdoc_site/*.md content pages from README_zh.md.

Splits the Chinese README by its section headings (H2/H3), rewrites
cross-page links and local asset/config paths, and demotes heading levels
so each generated page has a single H1. It also copies ``figs/`` into
``mkdoc_site/assets/`` so local image links resolve.

``mkdoc_site/index.md`` (the homepage) is left untouched. The hand-written
reference docs live in ``docs/`` (Design.md, Cases.md, Skills.md,
mineru_parse.md, undetected_fallback.md, PaperDB/) and are not part of the
MkDocs source (``docs_dir`` is ``mkdoc_site``), so they are never rendered.

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

# (heading_level, keyword, target_relpath). ``keyword`` is matched as a
# substring of the heading text (emoji included in the heading are ignored by
# this check); the first matching rule wins.
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


def match_heading(line: str) -> tuple[int, str] | None:
    m = HEADING_RE.match(line)
    if not m:
        return None
    return len(m.group(1)), m.group(2).strip()


def resolve_target(level: int, text: str) -> str | None:
    for rule_level, keyword, target in RULES:
        if rule_level == level and keyword in text:
            return target
    return None


def gh_url(path: str) -> str:
    """GitHub URL for a repo path: tree/ for directories, blob/ for files."""
    kind = "tree" if path.endswith("/") else "blob"
    return f"{GH}/{kind}/main/{path}"


def process(content: str, prefix: str) -> str:
    # 1) repo-file links (docs/config/test) -> GitHub URLs. These files are
    #    not rendered into the site, so link out to the repo instead.
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

    # 4) demote headings so the page's first heading becomes H1 (skip code fences)
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


def main() -> None:
    # refresh docs/assets from figs/ so image links resolve
    shutil.rmtree(DOCS / "assets", ignore_errors=True)
    shutil.copytree(ROOT / "figs", DOCS / "assets")

    lines = README.read_text(encoding="utf-8").splitlines(keepends=True)

    buffers: dict[str, list[str]] = {}
    order: list[str] = []
    current: str | None = None
    in_fence = False

    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("```") or stripped.startswith("~~~"):
            in_fence = not in_fence
            if current is not None:
                buffers[current].append(line)
            continue

        if not in_fence:
            heading = match_heading(line)
            if heading is not None:
                target = resolve_target(*heading)
                if target is not None:
                    current = target
                    if target not in buffers:
                        buffers[target] = []
                        order.append(target)
                    buffers[current].append(line)
                    continue

        if current is not None:
            buffers[current].append(line)

    for target in order:
        body = "".join(buffers[target])
        rendered = process(body, "../" * target.count("/"))
        out = DOCS / target
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(rendered, encoding="utf-8")
        print(f"wrote {target} ({len(buffers[target])} lines)")

    print("done")


if __name__ == "__main__":
    main()
