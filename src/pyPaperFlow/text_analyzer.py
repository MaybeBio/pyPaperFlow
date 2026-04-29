"""Text-level analyzers for parsed paper content."""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List


SECTION_PATTERNS = {
    "introduction": re.compile(r"^\s{0,3}(\d+\.?\s*)?(introduction)\b", re.I | re.M),
    "methods": re.compile(r"^\s{0,3}(\d+\.?\s*)?(materials and methods|methods)\b", re.I | re.M),
    "results": re.compile(r"^\s{0,3}(\d+\.?\s*)?(results)\b", re.I | re.M),
    "discussion": re.compile(r"^\s{0,3}(\d+\.?\s*)?(discussion)\b", re.I | re.M),
    "conclusion": re.compile(r"^\s{0,3}(\d+\.?\s*)?(conclusion|conclusions)\b", re.I | re.M),
    "references": re.compile(r"^\s{0,3}(references|bibliography)\b", re.I | re.M),
}


@dataclass
class SectionMatch:
    name: str
    start: int


class TextAnalyzer:
    """Parse section-level structure and extract lightweight signals from full text."""

    def detect_sections(self, text: str) -> Dict[str, int]:
        matches: Dict[str, int] = {}
        for name, pattern in SECTION_PATTERNS.items():
            found = pattern.search(text or "")
            if found:
                matches[name] = found.start()
        return dict(sorted(matches.items(), key=lambda kv: kv[1]))

    def extract_references_block(self, text: str) -> str:
        sec = self.detect_sections(text)
        if "references" not in sec:
            return ""
        start = sec["references"]
        return (text or "")[start:].strip()

    def extract_scientific_entities(self, text: str) -> Dict[str, List[str]]:
        content = text or ""

        proteins = sorted(set(re.findall(r"\b[A-Z0-9]{3,10}\b", content)))
        genes = sorted(set(re.findall(r"\b[A-Z]{2,8}\d{0,3}\b", content)))
        compounds = sorted(set(re.findall(r"\b[A-Z][a-z]{2,20}(?:ine|ate|ide|ol|one|acid)\b", content)))

        return {
            "proteins": proteins[:200],
            "genes": genes[:200],
            "compounds": compounds[:200],
        }
