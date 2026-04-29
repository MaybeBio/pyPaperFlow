"""Semantic analysis helpers without heavy ML dependencies."""

from __future__ import annotations

import math
import re
from collections import Counter
from typing import Dict, List, Tuple


TOKEN_RE = re.compile(r"[A-Za-z][A-Za-z0-9_-]{1,}")
STOPWORDS = {
    "the", "and", "for", "with", "that", "this", "from", "are", "was", "were", "have",
    "has", "had", "into", "their", "they", "our", "not", "can", "may", "also", "using",
    "used", "than", "which", "such", "over", "between", "both", "within", "each", "more",
}


class SemanticAnalyzer:
    """Keyword extraction and pairwise similarity over paper text."""

    def tokenize(self, text: str) -> List[str]:
        return [t.lower() for t in TOKEN_RE.findall(text or "") if t.lower() not in STOPWORDS]

    def top_keywords(self, docs: List[str], top_k: int = 30) -> List[Tuple[str, float]]:
        tf = Counter()
        df = Counter()

        for doc in docs:
            tokens = self.tokenize(doc)
            tf.update(tokens)
            df.update(set(tokens))

        n_docs = max(len(docs), 1)
        scored = []
        for token, freq in tf.items():
            idf = math.log((1 + n_docs) / (1 + df[token])) + 1.0
            scored.append((token, float(freq) * idf))

        scored.sort(key=lambda x: x[1], reverse=True)
        return scored[:top_k]

    def jaccard_similarity(self, text_a: str, text_b: str) -> float:
        a = set(self.tokenize(text_a))
        b = set(self.tokenize(text_b))
        if not a and not b:
            return 1.0
        if not a or not b:
            return 0.0
        return len(a & b) / len(a | b)

    def pairwise_similarity(self, docs: Dict[str, str]) -> Dict[str, Dict[str, float]]:
        keys = list(docs.keys())
        sim: Dict[str, Dict[str, float]] = {k: {} for k in keys}
        for i, a in enumerate(keys):
            for j, b in enumerate(keys):
                if j < i:
                    sim[a][b] = sim[b][a]
                else:
                    sim[a][b] = self.jaccard_similarity(docs.get(a, ""), docs.get(b, ""))
        return sim
