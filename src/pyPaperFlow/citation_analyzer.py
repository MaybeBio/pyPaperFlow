"""Citation graph construction and basic network metrics."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Dict, List, Set, Tuple


@dataclass
class CitationGraph:
    nodes: Set[str]
    edges: List[Tuple[str, str]]
    in_degree: Dict[str, int]
    out_degree: Dict[str, int]


class CitationAnalyzer:
    """Build citation relationships from paper metadata links."""

    def build_direct_citation_graph(self, papers: List[Dict]) -> CitationGraph:
        nodes: Set[str] = set()
        edges: List[Tuple[str, str]] = []
        in_deg = Counter()
        out_deg = Counter()

        for paper in papers:
            pmid = str((paper.get("identity", {}) or {}).get("pmid") or paper.get("pmid") or "")
            if not pmid:
                continue

            nodes.add(pmid)
            links = paper.get("links", {}) or {}
            refs = links.get("refs", [])
            if not isinstance(refs, list):
                refs = []

            for ref in refs:
                ref_id = str(ref).strip()
                if not ref_id:
                    continue
                nodes.add(ref_id)
                edges.append((pmid, ref_id))
                out_deg[pmid] += 1
                in_deg[ref_id] += 1

        return CitationGraph(nodes=nodes, edges=edges, in_degree=dict(in_deg), out_degree=dict(out_deg))

    def co_citation_pairs(self, graph: CitationGraph, min_shared: int = 2) -> Dict[str, int]:
        citing_by_target: Dict[str, Set[str]] = defaultdict(set)
        for src, dst in graph.edges:
            citing_by_target[dst].add(src)

        targets = list(citing_by_target.keys())
        score: Dict[str, int] = {}

        for i in range(len(targets)):
            for j in range(i + 1, len(targets)):
                t1, t2 = targets[i], targets[j]
                shared = len(citing_by_target[t1] & citing_by_target[t2])
                if shared >= min_shared:
                    score[f"{t1}|{t2}"] = shared

        return score

    def core_papers(self, graph: CitationGraph, top_k: int = 20) -> List[Tuple[str, int]]:
        ranked = sorted(graph.in_degree.items(), key=lambda kv: kv[1], reverse=True)
        return ranked[:top_k]
