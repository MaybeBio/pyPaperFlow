"""Search enhancement helpers for paper corpus retrieval."""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List


@dataclass
class SearchResult:
    paper_id: str
    score: float


class SearchEnhancer:
    """A lightweight lexical search baseline with query expansion."""

    DEFAULT_SYNONYMS = {
        "idr": ["intrinsically disordered region", "disordered protein"],
        "llm": ["large language model", "foundation model"],
        "ptm": ["post-translational modification", "post translational modification"],
    }

    def __init__(self, synonyms: Dict[str, List[str]] | None = None):
        self.synonyms = synonyms or self.DEFAULT_SYNONYMS

    def expand_query(self, query: str) -> List[str]:
        terms = [query.strip()]
        low = query.lower()
        for key, values in self.synonyms.items():
            if key in low:
                terms.extend(values)
        return [t for t in terms if t]

    def search(self, query: str, docs: Dict[str, str], top_k: int = 20) -> List[SearchResult]:
        expanded = self.expand_query(query)
        keywords = [re.escape(term.lower()) for term in expanded]
        if not keywords:
            return []

        pattern = re.compile("|".join(keywords), re.I)
        scored = []

        for paper_id, text in docs.items():
            hits = pattern.findall((text or "").lower())
            if hits:
                score = float(len(hits))
                scored.append(SearchResult(paper_id=paper_id, score=score))

        scored.sort(key=lambda x: x.score, reverse=True)
        return scored[:top_k]
