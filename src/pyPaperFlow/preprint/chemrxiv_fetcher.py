from __future__ import annotations

import re
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
from urllib.parse import quote

import httpx
from bs4 import BeautifulSoup

from ..integrations.cloak_fallback import cloak_fetch_pdf, is_cloak_enabled
from ..integrations.undetected_fallback import is_undetected_enabled, undetected_fetch_pdf

from .source_models import SourcePaper
from .source_utils import (
    basic_boolean_text_match,
    build_source_record_dir,
    download_binary,
    extract_year,
    normalize_text,
    safe_filename,
    save_json,
)


CHEM_RXIV_CROSSREF_API = "https://api.crossref.org/works"
CHEM_RXIV_CROSSREF_PREFIX = "10.26434"
CHEM_RXIV_JOURNAL = "ChemRxiv"
CHEM_RXIV_LAUNCH_DATE = datetime(2017, 8, 1)
CHEM_RXIV_PDF_BASE = "https://chemrxiv.org/doi/pdf"

DOI_RE = re.compile(r"^10\.\d{4,9}/[^\s]+$")


class ChemRxivFetcher:
    """Search ChemRxiv preprints via Crossref and download their PDFs.

    ChemRxiv metadata is deposited to Crossref under prefix 10.26434
    (publisher "American Chemical Society (ACS)", type ``posted-content``), so
    search mirrors the bioRxiv Crossref path. The ChemRxiv public API is
    Cloudflare-walled; the PDF route (``/doi/pdf/{doi}``) usually downloads
    directly, but falls back to the same opt-in browser routes (CloakBrowser /
    undetected_chromedriver) as bioRxiv/medRxiv when it does not.
    """

    def __init__(
        self,
        root_dir: str,
        max_retries: int = 3,
        request_timeout: float = 60.0,
    ):
        self.root_dir = root_dir
        self.max_retries = max(1, int(max_retries))
        self.request_timeout = float(request_timeout)
        self.headers = {
            "User-Agent": "pyPaperFlow/0.1.0 (+https://github.com/MaybeBio/pyPaperFlow)",
            "Accept": "application/json,text/html;q=0.9,*/*;q=0.8",
        }
        self._http_client: Optional[httpx.Client] = None

    def close(self) -> None:
        if self._http_client is not None:
            self._http_client.close()
            self._http_client = None

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _get_http_client(self) -> httpx.Client:
        if self._http_client is None:
            self._http_client = httpx.Client(
                headers=self.headers,
                timeout=self.request_timeout,
                follow_redirects=True,
            )
        return self._http_client

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

        if DOI_RE.match(query_text):
            return self._search_by_doi(query_text)

        return self._search_crossref(
            query_text=query_text,
            start_date=start_date,
            end_date=end_date,
            max_results=max_results,
        )

    def _search_crossref(
        self,
        query_text: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[SourcePaper]:
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
                if self._is_chemrxiv_record(raw_record):
                    if basic_boolean_text_match(self._search_text_crossref(raw_record), query_text):
                        records.append(self._normalize_crossref_record(raw_record, query=query_text))
                        if max_results is not None and len(records) >= max_results:
                            return records

            if len(items) < page_size:
                break

            next_cursor = normalize_text(message.get("next-cursor", ""))
            if not next_cursor or next_cursor == cursor:
                break
            cursor = next_cursor

        return records

    def _search_by_doi(self, doi: str) -> List[SourcePaper]:
        record = self._fetch_crossref_work(doi)
        if not record or not self._is_chemrxiv_record(record):
            return []
        return [self._normalize_crossref_record(record, query=doi)]

    def _fetch_crossref_work(self, doi: str) -> Dict[str, Any]:
        url = f"{CHEM_RXIV_CROSSREF_API}/{quote(doi, safe='')}"
        last_error: Optional[Exception] = None
        for attempt in range(self.max_retries):
            try:
                response = self._get_http_client().get(url)
                if response.status_code == 404:
                    return {}
                response.raise_for_status()
                return response.json().get("message") or {}
            except Exception as exc:
                last_error = exc
                if attempt + 1 < self.max_retries:
                    time.sleep(min(2.0, 0.5 * (attempt + 1)))
        if last_error is not None:
            raise last_error
        raise RuntimeError(f"Failed to fetch Crossref work for DOI {doi}")

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

    def fetch_from_dois(
        self,
        dois: Iterable[str],
        output_dir: Optional[str] = None,
        download_pdf: bool = True,
    ) -> List[SourcePaper]:
        records: List[SourcePaper] = []
        for raw_doi in dois:
            doi = normalize_text(raw_doi)
            if not doi:
                continue
            matched = self._search_by_doi(doi)
            if not matched:
                continue
            record = matched[0]
            self._save_record(record, output_dir=output_dir, download_pdf=download_pdf)
            records.append(record)
        return records

    def _normalize_date_range(
        self,
        start_date: Optional[str],
        end_date: Optional[str],
    ) -> tuple[datetime, datetime]:
        if start_date:
            start_dt = datetime.strptime(start_date, "%Y-%m-%d")
        else:
            start_dt = CHEM_RXIV_LAUNCH_DATE

        if end_date:
            end_dt = datetime.strptime(end_date, "%Y-%m-%d")
        else:
            end_dt = datetime.utcnow()

        if start_dt < CHEM_RXIV_LAUNCH_DATE:
            start_dt = CHEM_RXIV_LAUNCH_DATE
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
            "filter": f"prefix:{CHEM_RXIV_CROSSREF_PREFIX},type:posted-content",
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
                response = self._get_http_client().get(CHEM_RXIV_CROSSREF_API, params=params)
                response.raise_for_status()
                return response.json()
            except Exception as exc:
                last_error = exc
                if attempt + 1 < self.max_retries:
                    time.sleep(min(2.0, 0.5 * (attempt + 1)))

        if last_error is not None:
            raise last_error
        raise RuntimeError("Failed to query ChemRxiv Crossref search")

    def _is_chemrxiv_record(self, record: Dict[str, Any]) -> bool:
        doi = normalize_text(record.get("DOI", ""))
        return doi.startswith(CHEM_RXIV_CROSSREF_PREFIX + "/")

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
        doi = normalize_text(record.get("DOI", ""))
        landing_url = f"https://doi.org/{doi}" if doi else ""

        source_id = doi or safe_filename(f"chemrxiv_{published_date}_{title}")
        keywords = self._normalize_crossref_keywords(record.get("subject", []))
        if not keywords:
            group_title = normalize_text(record.get("group-title", ""))
            keywords = [group_title] if group_title else []
        pdf_url = self._crossref_pdf_url(record) or (f"{CHEM_RXIV_PDF_BASE}/{doi}" if doi else "")

        return SourcePaper(
            source="chemrxiv",
            source_id=source_id,
            title=title,
            doi=doi,
            abstract=abstract,
            authors=authors,
            published_date=published_date,
            updated_date=updated_date,
            journal=CHEM_RXIV_JOURNAL,
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
                "group_title": record.get("group-title", ""),
                "raw_record": record,
            },
        )

    def _crossref_pdf_url(self, record: Dict[str, Any]) -> str:
        links = record.get("link") or []
        if not isinstance(links, list):
            return ""
        for link in links:
            if not isinstance(link, dict):
                continue
            url = normalize_text(link.get("URL", ""))
            if "/doi/pdf/" in url or url.rstrip("/").endswith("/pdf"):
                return url
        return ""

    def _normalize_crossref_authors(self, authors_value: Any) -> List[str]:
        if not isinstance(authors_value, list):
            return []

        authors: List[str] = []
        for author in authors_value:
            if not isinstance(author, dict):
                continue
            given = normalize_text(author.get("given", ""))
            family = normalize_text(author.get("family", ""))
            full_name = " ".join(part for part in [given, family] if part)
            name = normalize_text(author.get("name") or full_name or given)
            if name:
                authors.append(name)
        return authors

    def _normalize_crossref_keywords(self, keywords_value: Any) -> List[str]:
        if isinstance(keywords_value, list):
            return [normalize_text(item) for item in keywords_value if normalize_text(item)]
        if isinstance(keywords_value, str) and keywords_value.strip():
            return [part.strip() for part in keywords_value.split(",") if part.strip()]
        return []

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
        cleaned = re.sub(r"^(?:Abstract|Summary)\b[\s:：\-–—]*", "", cleaned, flags=re.IGNORECASE)
        return cleaned

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
        candidates = [
            candidate
            for candidate in [record.pdf_url]
            if candidate
        ]

        for candidate in candidates:
            if download_binary(candidate, pdf_path, headers=self.headers, timeout=self.request_timeout):
                return True

        # When chemrxiv.org serves a Cloudflare challenge to non-browser
        # clients (HTTP 403), two opt-in browser fallbacks can clear it,
        # mirroring the bioRxiv/medRxiv path.
        if is_cloak_enabled():
            for candidate in candidates:
                data = cloak_fetch_pdf(candidate, timeout=int(self.request_timeout))
                if data and data[:5].startswith(b"%PDF"):
                    try:
                        pdf_path.parent.mkdir(parents=True, exist_ok=True)
                        pdf_path.write_bytes(data)
                        return True
                    except OSError:
                        return False

        if is_undetected_enabled():
            for candidate in candidates:
                data = undetected_fetch_pdf(candidate, timeout=int(self.request_timeout))
                if data and data[:5].startswith(b"%PDF"):
                    try:
                        pdf_path.parent.mkdir(parents=True, exist_ok=True)
                        pdf_path.write_bytes(data)
                        return True
                    except OSError:
                        return False

        return False
