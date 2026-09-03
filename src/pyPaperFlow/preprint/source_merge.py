from __future__ import annotations

import re
from typing import List

from .source_models import SourcePaper
from .source_utils import normalize_text

_DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.IGNORECASE)


def extract_doi(text) -> str:
    """Extract a DOI from arbitrary text or a URL, or return "" if none is present."""
    if not text:
        return ""
    match = _DOI_RE.search(str(text))
    return match.group(0).rstrip(".,;)") if match else ""


def _paper_key(record: SourcePaper) -> str:
    doi = normalize_text(record.doi).lower()
    if doi:
        return f"doi:{doi}"
    title = normalize_text(record.title).lower()
    authors = " ".join(a.strip().lower() for a in record.authors)
    if title:
        return f"title:{title}|authors:{authors}"
    return f"id:{normalize_text(record.source_id)}"


def dedupe_papers(records: List[SourcePaper]) -> List[SourcePaper]:
    """Deduplicate records by DOI, falling back to title+authors and finally source_id.

    First occurrence wins; later records sharing the same key are dropped.
    """
    seen: set = set()
    out: List[SourcePaper] = []
    for record in records:
        key = _paper_key(record)
        if key in seen:
            continue
        seen.add(key)
        out.append(record)
    return out


def merge_papers(records: List[SourcePaper]) -> List[SourcePaper]:
    """Merge records from multiple sources: backfill missing DOIs, then deduplicate."""
    for record in records:
        if not normalize_text(record.doi):
            record.doi = extract_doi(record.landing_url or record.abstract or record.title)
    return dedupe_papers(records)
