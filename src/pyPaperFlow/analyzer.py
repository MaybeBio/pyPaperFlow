"""Statistical analysis utilities for collected PubMed papers."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Any


@dataclass
class AnalyzerResult:
    paper_count: int
    year_distribution: Dict[str, int]
    journal_distribution: Dict[str, int]
    top_authors: List[Dict[str, Any]]
    citation_summary: Dict[str, float]


class PaperAnalyzer:
    """Compute lightweight corpus-level statistics from merged paper records."""

    def analyze(self, papers: List[Dict[str, Any]], top_k_authors: int = 20) -> AnalyzerResult:
        year_counter: Counter[str] = Counter()
        journal_counter: Counter[str] = Counter()
        author_counter: Counter[str] = Counter()

        total_refs = 0
        total_review_links = 0

        for paper in papers:
            source = paper.get("source", {}) or {}
            contributors = paper.get("contributors", {}) or {}
            links = paper.get("links", {}) or {}

            pub_date = str(source.get("pub_date", "")).strip()
            year = pub_date[:4] if len(pub_date) >= 4 and pub_date[:4].isdigit() else "unknown"
            year_counter[year] += 1

            journal = str(source.get("journal_title", "unknown")).strip() or "unknown"
            journal_counter[journal] += 1

            medline = contributors.get("medline", {}) if isinstance(contributors, dict) else {}
            full_names = medline.get("full_names", []) if isinstance(medline, dict) else []
            for name in full_names:
                if name:
                    author_counter[str(name).strip()] += 1

            refs = links.get("refs", []) if isinstance(links, dict) else []
            review = links.get("review", []) if isinstance(links, dict) else []
            total_refs += len(refs) if isinstance(refs, list) else 0
            total_review_links += len(review) if isinstance(review, list) else 0

        paper_count = len(papers)
        avg_refs = total_refs / paper_count if paper_count else 0.0
        avg_review_links = total_review_links / paper_count if paper_count else 0.0

        top_authors = [
            {"author": name, "count": count}
            for name, count in author_counter.most_common(top_k_authors)
        ]

        return AnalyzerResult(
            paper_count=paper_count,
            year_distribution=dict(year_counter),
            journal_distribution=dict(journal_counter),
            top_authors=top_authors,
            citation_summary={
                "total_refs": float(total_refs),
                "total_review_links": float(total_review_links),
                "avg_refs_per_paper": avg_refs,
                "avg_review_links_per_paper": avg_review_links,
            },
        )
