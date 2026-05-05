from __future__ import annotations

import time
import re
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import requests
from bs4 import BeautifulSoup

from .source_models import SourcePaper
from .source_utils import (
    basic_boolean_text_match,
    build_source_record_dir,
    download_binary,
    ensure_directory,
    extract_year,
    normalize_text,
    safe_filename,
    save_json,
)


BIO_RXIV_API_BASE = "https://api.biorxiv.org/details/biorxiv"
BIO_RXIV_CROSSREF_API = "https://api.crossref.org/works"
BIO_RXIV_CROSSREF_PREFIX = "10.64898"
BIO_RXIV_LANDING_BASE = "https://www.biorxiv.org/content"
BIO_RXIV_LAUNCH_DATE = datetime(2013, 1, 1)


class BioRxivFetcher:
    def __init__(
        self,
        root_dir: str,
        window_days: int = 365,
        max_retries: int = 3,
        request_timeout: float = 60.0,
    ):
        self.root_dir = root_dir
        self.window_days = max(1, int(window_days))
        self.max_retries = max(1, int(max_retries))
        self.request_timeout = float(request_timeout)
        self.headers = {
            "User-Agent": "pyPaperFlow/0.1.0 (+https://github.com/MaybeBio/pyPaperFlow)",
            "Accept": "application/json,text/html;q=0.9,*/*;q=0.8",
        }

    def search(
        self,
        query: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[SourcePaper]:
        query_text = normalize_text(query)
        if not query_text:
            raise ValueError("query must be non-empty")

        records: List[SourcePaper] = []
        start_dt, end_dt = self._normalize_date_range(start_date, end_date)

        cursor = "*"
        while True:
            page_size = 1000
            if max_results is not None:
                remaining = max(1, int(max_results) - len(records))
                page_size = min(page_size, remaining)

            payload = self._request_crossref_page(
                query_text=query_text,
                cursor=cursor,
                page_size=page_size,
                start_dt=start_dt,
                end_dt=end_dt,
            )
            message = payload.get("message") or {}
            items = message.get("items") or []
            if not items:
                break

            for raw_record in items:
                if normalize_text(raw_record.get("publisher", "")).lower() != "openrxiv":
                    continue
                if not basic_boolean_text_match(self._search_text_crossref(raw_record), query_text):
                    continue
                record = self._normalize_crossref_record(raw_record, query=query_text)
                records.append(record)
                if max_results is not None and len(records) >= max_results:
                    return records

            if len(items) < page_size:
                break

            next_cursor = normalize_text(message.get("next-cursor", ""))
            if not next_cursor or next_cursor == cursor:
                break
            cursor = next_cursor

        return records

    def fetch_from_query(
        self,
        query: str,
        output_dir: Optional[str] = None,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        max_results: Optional[int] = None,
        download_pdf: bool = True,
    ) -> List[SourcePaper]:
        records = self.search(
            query=query,
            start_date=start_date,
            end_date=end_date,
            max_results=max_results,
        )

        for record in records:
            self._save_record(record, output_dir=output_dir, download_pdf=download_pdf)

        return records

    def _normalize_date_range(
        self,
        start_date: Optional[str],
        end_date: Optional[str],
    ) -> tuple[datetime, datetime]:
        if start_date:
            start_dt = datetime.strptime(start_date, "%Y-%m-%d")
        else:
            start_dt = BIO_RXIV_LAUNCH_DATE

        if end_date:
            end_dt = datetime.strptime(end_date, "%Y-%m-%d")
        else:
            end_dt = datetime.utcnow()

        if start_dt < BIO_RXIV_LAUNCH_DATE:
            start_dt = BIO_RXIV_LAUNCH_DATE
        if start_dt > end_dt:
            raise ValueError("start_date cannot be later than end_date")
        return start_dt, end_dt

    def _request_crossref_page(
        self,
        query_text: str,
        cursor: str,
        page_size: int,
        start_dt: datetime,
        end_dt: datetime,
    ) -> Dict[str, Any]:
        params: Dict[str, Any] = {
            "filter": f"prefix:{BIO_RXIV_CROSSREF_PREFIX},type:posted-content",
            "query.bibliographic": query_text,
            "rows": page_size,
            "cursor": cursor,
            "sort": "relevance",
        }
        if start_dt:
            params["filter"] += f",from-pub-date:{start_dt.strftime('%Y-%m-%d')}"
        if end_dt:
            params["filter"] += f",until-pub-date:{end_dt.strftime('%Y-%m-%d')}"

        last_error: Optional[Exception] = None

        for attempt in range(self.max_retries):
            try:
                response = requests.get(
                    BIO_RXIV_CROSSREF_API,
                    params=params,
                    headers=self.headers,
                    timeout=self.request_timeout,
                )
                response.raise_for_status()
                return response.json()
            except Exception as exc:
                last_error = exc
                if attempt + 1 < self.max_retries:
                    time.sleep(min(2.0, 0.5 * (attempt + 1)))

        if last_error is not None:
            raise last_error
        raise RuntimeError("Failed to query bioRxiv Crossref search")

    def _search_text_crossref(self, record: Dict[str, Any]) -> str:
        pieces = [
            record.get("title", ""),
            record.get("doi", ""),
            record.get("publisher", ""),
            record.get("container-title", ""),
            record.get("subject", ""),
            record.get("abstract", ""),
        ]
        normalized_parts: List[str] = []
        for piece in pieces:
            if isinstance(piece, list):
                normalized_piece = " ".join(normalize_text(item) for item in piece if normalize_text(item))
            else:
                normalized_piece = normalize_text(piece)
            if normalized_piece:
                normalized_parts.append(normalized_piece)
        return " ".join(normalized_parts)

    def _normalize_crossref_record(self, record: Dict[str, Any], query: str) -> SourcePaper:
        title = normalize_text((record.get("title") or [""])[0])
        abstract = self._clean_crossref_abstract(record.get("abstract", ""))
        published_date = self._extract_crossref_date(record)
        updated_date = normalize_text(record.get("created", {}).get("date-time", ""))

        authors = self._normalize_crossref_authors(record.get("author", []))
        doi = self._normalize_crossref_doi(record.get("DOI", ""))
        landing_url = self._landing_url(doi, "")

        source_id = doi or safe_filename(f"biorxiv_{published_date}_{title}")
        keywords = self._normalize_crossref_keywords(record.get("subject", []))
        pdf_url = f"{BIO_RXIV_LANDING_BASE}/{doi}.full.pdf" if doi else ""

        return SourcePaper(
            source="biorxiv",
            source_id=source_id,
            title=title,
            doi=doi,
            abstract=abstract,
            authors=authors,
            published_date=published_date,
            updated_date=updated_date,
            journal="bioRxiv",
            category=", ".join(keywords),
            landing_url=landing_url,
            pdf_url=pdf_url,
            query=query,
            version="",
            keywords=keywords,
            extra={
                "publisher": record.get("publisher", ""),
                "prefix": record.get("prefix", ""),
                "type": record.get("type", ""),
                "raw_record": record,
            },
        )

    def _normalize_crossref_authors(self, authors_value: Any) -> List[str]:
        if not isinstance(authors_value, list):
            return []

        authors: List[str] = []
        for author in authors_value:
            if not isinstance(author, dict):
                continue
            name = normalize_text(
                author.get("name")
                or author.get("given")
                or " ".join(
                    part
                    for part in [author.get("given", ""), author.get("family", "")]
                    if normalize_text(part)
                )
            )
            if name:
                authors.append(name)
        return authors

    def _normalize_crossref_keywords(self, keywords_value: Any) -> List[str]:
        if isinstance(keywords_value, list):
            return [normalize_text(item) for item in keywords_value if normalize_text(item)]
        if isinstance(keywords_value, str) and keywords_value.strip():
            return [part.strip() for part in keywords_value.split(",") if part.strip()]
        return []

    def _normalize_crossref_doi(self, doi: Any) -> str:
        return normalize_text(doi)

    def _extract_crossref_date(self, record: Dict[str, Any]) -> str:
        for key in ("published-online", "published-print", "issued", "created"):
            value = record.get(key)
            if isinstance(value, dict):
                date_parts = value.get("date-parts") or []
                if date_parts and date_parts[0]:
                    parts = date_parts[0]
                    year = parts[0] if len(parts) > 0 else None
                    month = parts[1] if len(parts) > 1 else 1
                    day = parts[2] if len(parts) > 2 else 1
                    if year:
                        try:
                            return datetime(int(year), int(month), int(day)).strftime("%Y-%m-%d")
                        except Exception:
                            return normalize_text(year)
        return ""

    def _clean_crossref_abstract(self, abstract: Any) -> str:
        text = normalize_text(abstract)
        if not text:
            return ""
        soup = BeautifulSoup(text, "html.parser")
        cleaned = soup.get_text(separator=" ")
        cleaned = re.sub(r"\s+", " ", cleaned).strip()
        return cleaned

    def _normalize_authors(self, authors_value: Any) -> List[str]:
        if isinstance(authors_value, list):
            authors: List[str] = []
            for author in authors_value:
                if isinstance(author, dict):
                    name = normalize_text(
                        author.get("name")
                        or author.get("author_name")
                        or author.get("full_name")
                        or " ".join(
                            part
                            for part in [author.get("given", ""), author.get("family", "")]
                            if normalize_text(part)
                        )
                    )
                else:
                    name = normalize_text(author)
                if name:
                    authors.append(name)
            return authors

        if isinstance(authors_value, str):
            parts = [part.strip() for part in authors_value.split(";") if part.strip()]
            return parts or [authors_value]

        return []

    def _normalize_keywords(self, category_value: Any) -> List[str]:
        if isinstance(category_value, list):
            return [normalize_text(item) for item in category_value if normalize_text(item)]
        if isinstance(category_value, str) and category_value.strip():
            return [part.strip() for part in category_value.split(",") if part.strip()]
        return []

    def _normalize_doi(self, doi: Any, version: Any) -> str:
        doi_text = normalize_text(doi)
        version_text = normalize_text(version)
        if doi_text and version_text and not doi_text.lower().endswith(f"v{version_text.lower()}"):
            return f"{doi_text}v{version_text}"
        return doi_text

    def _landing_url(self, doi: str, version: str) -> str:
        if not doi:
            return ""
        return f"{BIO_RXIV_LANDING_BASE}/{doi}"

    def _candidate_pdf_urls(self, record: Dict[str, Any], doi: str, version: str) -> List[str]:
        candidates: List[str] = []
        for key in ("pdf_url", "full_text_url", "fulltext_url", "pdf"):
            value = normalize_text(record.get(key, ""))
            if value:
                candidates.append(value)

        if doi:
            candidates.append(f"{BIO_RXIV_LANDING_BASE}/{doi}.full.pdf")
            if version and not doi.lower().endswith(f"v{version.lower()}"):
                candidates.insert(0, f"{BIO_RXIV_LANDING_BASE}/{doi}v{version}.full.pdf")

        landing_url = self._landing_url(doi, version)
        if landing_url:
            candidates.append(landing_url)

        deduped: List[str] = []
        seen = set()
        for candidate in candidates:
            if candidate not in seen:
                seen.add(candidate)
                deduped.append(candidate)
        return deduped

    def _save_record(self, record: SourcePaper, output_dir: Optional[str], download_pdf: bool) -> None:
        base_dir = output_dir or self.root_dir
        year = extract_year(record.published_date)
        record_dir = build_source_record_dir(base_dir, record.source, year, record.source_id)
        file_stem = safe_filename(record.source_id)

        if download_pdf and record.pdf_url:
            pdf_path = record_dir / f"{file_stem}.pdf"
            if self._download_pdf(record, pdf_path):
                record.pdf_downloaded = True
                record.pdf_path = str(pdf_path)

        save_json(record_dir / f"{file_stem}.json", record.to_dict())

    def _download_pdf(self, record: SourcePaper, pdf_path: Path) -> bool:
        candidates = [candidate for candidate in [record.pdf_url, record.landing_url] if candidate]

        for candidate in candidates:
            if candidate.lower().endswith(".pdf"):
                if download_binary(candidate, pdf_path, headers=self.headers, timeout=self.request_timeout):
                    return True
                continue

            try:
                response = requests.get(candidate, headers=self.headers, timeout=self.request_timeout)
                response.raise_for_status()
                soup = BeautifulSoup(response.text, "html.parser")
                meta = soup.find("meta", {"name": "citation_pdf_url"})
                if meta and meta.get("content"):
                    pdf_url = normalize_text(meta.get("content"))
                    if download_binary(pdf_url, pdf_path, headers=self.headers, timeout=self.request_timeout):
                        return True
            except Exception:
                continue

        return False