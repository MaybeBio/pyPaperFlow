from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List


@dataclass
class SourcePaper:
    source: str
    source_id: str
    title: str = ""
    doi: str = ""
    abstract: str = ""
    authors: List[str] = field(default_factory=list)
    published_date: str = ""
    updated_date: str = ""
    journal: str = ""
    category: str = ""
    landing_url: str = ""
    pdf_url: str = ""
    pdf_path: str = ""
    pdf_downloaded: bool = False
    query: str = ""
    version: str = ""
    keywords: List[str] = field(default_factory=list)
    extra: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)