from __future__ import annotations

import re
import time
import xml.etree.ElementTree as ET
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from .source_models import SourcePaper
from .source_utils import (
    build_source_record_dir,
    download_binary,
    extract_year,
    normalize_text,
    safe_filename,
    save_json,
)


ARXIV_API_URL = "https://export.arxiv.org/api/query"
ARXIV_ATOM_NS = {
    "atom": "http://www.w3.org/2005/Atom",
    "arxiv": "http://arxiv.org/schemas/atom",
}


class ArxivFetcher:
    def __init__(
        self,
        root_dir: str,
        batch_size: int = 100,
        max_retries: int = 3,
        request_timeout: float = 60.0,
    ):
        self.root_dir = root_dir
        self.batch_size = max(1, int(batch_size))
        self.max_retries = max(1, int(max_retries))
        self.request_timeout = float(request_timeout)
        self.headers = {
            "User-Agent": "pyPaperFlow/0.1.0 (+https://github.com/MaybeBio/pyPaperFlow)",
            "Accept": "application/atom+xml,application/xml;q=0.9,*/*;q=0.8",
        }

    def build_query(
        self,
        query: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
    ) -> str:
        raw_query = normalize_text(query)
        if not raw_query:
            raise ValueError("query must be non-empty")

        if any(marker in raw_query for marker in (":", "submittedDate:", "all:", "ti:", "au:", "abs:")):
            search_query = raw_query
        else:
            tokens = re.findall(r'"[^"]+"|\'[^\']+\'|\S+', raw_query)
            normalized_tokens: List[str] = []
            for token in tokens:
                clean = token.strip('"\'')
                if not clean:
                    continue
                if clean.upper() in {"AND", "OR", "NOT"}:
                    continue
                if " " in clean:
                    normalized_tokens.append(f'all:"{clean}"')
                else:
                    normalized_tokens.append(f"all:{clean}")
            search_query = " AND ".join(normalized_tokens)

        if start_date or end_date:
            start_text = self._format_arxiv_date(start_date or "1991-01-01", suffix="0000")
            end_text = self._format_arxiv_date(end_date or datetime.utcnow().strftime("%Y-%m-%d"), suffix="2359")
            date_filter = f"submittedDate:[{start_text} TO {end_text}]"
            search_query = f"({search_query}) AND {date_filter}"

        return search_query

    def search(
        self,
        query: str,
        max_results: int = 100,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
    ) -> List[SourcePaper]:
        search_query = self.build_query(query, start_date=start_date, end_date=end_date)
        results: List[SourcePaper] = []

        start = 0
        remaining = max(1, int(max_results))
        while remaining > 0:
            page_size = min(self.batch_size, remaining)
            feed = self._request_feed(search_query, start=start, max_results=page_size)
            entries = feed.findall("atom:entry", ARXIV_ATOM_NS)
            if not entries:
                break

            for entry in entries:
                record = self._normalize_entry(entry, query=query)
                results.append(record)
                remaining -= 1
                if remaining <= 0:
                    break

            if len(entries) < page_size or remaining <= 0:
                break

            start += len(entries)

        return results

    def fetch_from_query(
        self,
        query: str,
        output_dir: Optional[str] = None,
        max_results: int = 100,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        download_pdf: bool = True,
    ) -> List[SourcePaper]:
        records = self.search(
            query=query,
            max_results=max_results,
            start_date=start_date,
            end_date=end_date,
        )

        for record in records:
            self._save_record(record, output_dir=output_dir, download_pdf=download_pdf)

        return records

    def _request_feed(self, search_query: str, start: int, max_results: int) -> ET.Element:
        last_error: Optional[Exception] = None
        params = {
            "search_query": search_query,
            "start": start,
            "max_results": max_results,
            "sortBy": "submittedDate",
            "sortOrder": "descending",
        }

        for attempt in range(self.max_retries):
            try:
                response = requests.get(
                    ARXIV_API_URL,
                    params=params,
                    headers=self.headers,
                    timeout=self.request_timeout,
                )
                response.raise_for_status()
                return ET.fromstring(response.text)
            except Exception as exc:
                last_error = exc
                if attempt + 1 < self.max_retries:
                    time.sleep(min(2.0, 0.5 * (attempt + 1)))

        if last_error is not None:
            raise last_error
        raise RuntimeError("Failed to query arXiv API")

    def _normalize_entry(self, entry: ET.Element, query: str) -> SourcePaper:
        entry_id = normalize_text(entry.findtext("atom:id", default="", namespaces=ARXIV_ATOM_NS))
        source_id = entry_id.rsplit("/", 1)[-1] if entry_id else ""
        title = normalize_text(entry.findtext("atom:title", default="", namespaces=ARXIV_ATOM_NS))
        summary = normalize_text(entry.findtext("atom:summary", default="", namespaces=ARXIV_ATOM_NS))
        published = normalize_text(entry.findtext("atom:published", default="", namespaces=ARXIV_ATOM_NS))
        updated = normalize_text(entry.findtext("atom:updated", default="", namespaces=ARXIV_ATOM_NS))

        authors = [
            normalize_text(author.findtext("atom:name", default="", namespaces=ARXIV_ATOM_NS))
            for author in entry.findall("atom:author", ARXIV_ATOM_NS)
        ]
        authors = [author for author in authors if author]

        categories = [
            normalize_text(category.attrib.get("term", ""))
            for category in entry.findall("atom:category", ARXIV_ATOM_NS)
        ]
        categories = [category for category in categories if category]

        doi = normalize_text(entry.findtext("arxiv:doi", default="", namespaces=ARXIV_ATOM_NS))
        journal_ref = normalize_text(entry.findtext("arxiv:journal_ref", default="", namespaces=ARXIV_ATOM_NS))
        landing_url = f"https://arxiv.org/abs/{source_id}" if source_id else ""

        pdf_url = ""
        for link in entry.findall("atom:link", ARXIV_ATOM_NS):
            href = normalize_text(link.attrib.get("href", ""))
            link_type = normalize_text(link.attrib.get("type", ""))
            link_title = normalize_text(link.attrib.get("title", ""))
            if href and (link_type == "application/pdf" or link_title.lower() == "pdf"):
                pdf_url = href
                break
        if not pdf_url and source_id:
            pdf_url = f"https://arxiv.org/pdf/{source_id}.pdf"

        return SourcePaper(
            source="arxiv",
            source_id=source_id or safe_filename(title or query),
            title=title,
            doi=doi,
            abstract=summary,
            authors=authors,
            published_date=published,
            updated_date=updated,
            journal=journal_ref,
            category=", ".join(categories),
            landing_url=landing_url,
            pdf_url=pdf_url,
            query=query,
            version=source_id.rsplit("v", 1)[-1] if source_id and "v" in source_id else "",
            keywords=categories,
            extra={
                "entry_id": entry_id,
                "raw_categories": categories,
                "journal_ref": journal_ref,
            },
        )

    def _save_record(self, record: SourcePaper, output_dir: Optional[str], download_pdf: bool) -> None:
        base_dir = output_dir or self.root_dir
        year = extract_year(record.published_date)
        record_dir = build_source_record_dir(base_dir, record.source, year, record.source_id)
        file_stem = safe_filename(record.source_id)

        if download_pdf and record.pdf_url:
            pdf_path = record_dir / f"{file_stem}.pdf"
            if download_binary(record.pdf_url, pdf_path, headers=self.headers, timeout=self.request_timeout):
                record.pdf_downloaded = True
                record.pdf_path = str(pdf_path)

        save_json(record_dir / f"{file_stem}.json", record.to_dict())

    def _format_arxiv_date(self, date_text: str, suffix: str) -> str:
        parsed = datetime.strptime(date_text, "%Y-%m-%d")
        return parsed.strftime("%Y%m%d") + suffix