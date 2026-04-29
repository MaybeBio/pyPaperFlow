"""LLM-oriented prompt builders for downstream paper analysis."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List


@dataclass
class PromptPack:
    system: str
    user: str


class AIAnalyzer:
    """Build stable prompts for summarization, comparison, and extraction."""

    def build_multilevel_summary_prompt(self, paper_markdown: str, language: str = "zh") -> PromptPack:
        system = (
            "You are a scientific literature analyst. Produce factual, structured outputs "
            "without fabrication, and cite evidence snippets from the source text."
        )
        user = (
            f"Language: {language}\n"
            "Task:\n"
            "1) Write a 3-sentence ultra-short summary.\n"
            "2) Write a structured summary with sections: Problem, Method, Data, Results, Limitations.\n"
            "3) List 5 key claims with direct evidence quotes.\n\n"
            "Paper content:\n"
            f"{paper_markdown}"
        )
        return PromptPack(system=system, user=user)

    def build_comparison_prompt(self, papers: List[str], language: str = "zh") -> PromptPack:
        joined = "\n\n".join(f"[Paper {i + 1}]\n{p}" for i, p in enumerate(papers))
        system = "You compare scientific papers and return concise, evidence-grounded differences."
        user = (
            f"Language: {language}\n"
            "Compare these papers in a table with columns: method, dataset, key result, limitation, novelty.\n"
            "Then provide: (a) consensus points, (b) conflicts, (c) open questions.\n\n"
            f"{joined}"
        )
        return PromptPack(system=system, user=user)

    def build_method_tagging_prompt(self, paper_markdown: str, candidate_tags: List[str]) -> PromptPack:
        system = "You classify methods in scientific papers using only provided candidate tags."
        user = (
            "Select up to 5 tags from the candidate list and provide confidence scores in [0,1].\n"
            "Return JSON with fields: tags, rationale.\n"
            f"Candidate tags: {candidate_tags}\n\n"
            f"Paper content:\n{paper_markdown}"
        )
        return PromptPack(system=system, user=user)
