from __future__ import annotations

import sys
import time
import re
from datetime import datetime, timedelta
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
    detect_platform_from_doi,
    download_binary,
    ensure_directory,
    extract_year,
    normalize_text,
    safe_filename,
    save_json,
)


BIO_RXIV_CROSSREF_API = "https://api.crossref.org/works"
BIO_RXIV_CROSSREF_PREFIX = "10.64898"
BIO_RXIV_LAUNCH_DATE = datetime(2013, 1, 1)
MED_RXIV_LAUNCH_DATE = datetime(2019, 6, 1)

DOI_RE = re.compile(r"^10\.\d{4,9}/[^\s]+$")

PLATFORM_CONFIG = {
    "biorxiv": {
        "journal": "bioRxiv",
        "landing_base": "https://www.biorxiv.org/content",
        "launch_date": BIO_RXIV_LAUNCH_DATE,
    },
    "medrxiv": {
        "journal": "medRxiv",
        "landing_base": "https://www.medrxiv.org/content",
        "launch_date": MED_RXIV_LAUNCH_DATE,
    },
}


class BioRxivFetcher:
    def __init__(
        self,
        root_dir: str,
        platform: str = "biorxiv",
        window_days: int = 365,
        max_retries: int = 3,
        request_timeout: float = 60.0,
    ):
        if platform not in PLATFORM_CONFIG:
            raise ValueError(f"platform must be one of {sorted(PLATFORM_CONFIG)}, got {platform!r}")
        self.root_dir = root_dir
        self.platform = platform
        self.journal_name = PLATFORM_CONFIG[platform]["journal"]
        self.landing_base = PLATFORM_CONFIG[platform]["landing_base"]
        self.launch_date = PLATFORM_CONFIG[platform]["launch_date"]
        self.window_days = max(1, int(window_days))
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
        use_europepmc: bool = True,
    ) -> List[SourcePaper]:
        query_text = normalize_text(query)
        if not query_text:
            raise ValueError("query must be non-empty")

        if DOI_RE.match(query_text):
            return self._search_by_doi(query_text)

        records = self._search_crossref(
            query_text=query_text,
            start_date=start_date,
            end_date=end_date,
            max_results=max_results,
        )

        if use_europepmc:
            europepmc_records = self._search_europepmc(
                query_text=query_text,
                start_date=start_date,
                end_date=end_date,
                max_results=max_results,
            )
            records = self._union_records(records, europepmc_records)
            if max_results is not None:
                records = records[: int(max_results)]

        return records

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
                if normalize_text(raw_record.get("publisher", "")).lower() != "openrxiv":
                    continue
                if self._platform_from_record(raw_record) != self.platform:
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

    def _search_europepmc(
        self,
        query_text: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[SourcePaper]:
        try:
            from .europepmc_fetcher import EuropePMCSearch
        except Exception:
            return []

        start_dt, end_dt = self._normalize_date_range(start_date, end_date)
        start_str = start_dt.strftime("%Y-%m-%d")
        end_str = end_dt.strftime("%Y-%m-%d")

        searcher = EuropePMCSearch(
            request_timeout=self.request_timeout,
            max_retries=self.max_retries,
        )
        try:
            raw_records = searcher.search(
                query=query_text,
                start_date=start_str,
                end_date=end_str,
                max_results=max_results,
            )
        except Exception as exc:
            print(
                f"[biorxiv] Europe PMC search failed ({exc}); returning Crossref-only results.",
                file=sys.stderr,
            )
            return []
        finally:
            searcher.close()

        results: List[SourcePaper] = []
        for raw_record in raw_records:
            doi = normalize_text(raw_record.get("doi", ""))
            if not doi:
                continue
            if detect_platform_from_doi(doi) != self.platform:
                continue
            results.append(self._normalize_europepmc_record(raw_record, query=query_text))
        return results

    def _union_records(
        self,
        crossref_records: List[SourcePaper],
        europepmc_records: List[SourcePaper],
    ) -> List[SourcePaper]:
        merged: List[SourcePaper] = []
        seen: set = set()
        for record in list(crossref_records) + list(europepmc_records):
            key = normalize_text(record.doi or record.source_id).lower()
            if key in seen:
                continue
            seen.add(key)
            merged.append(record)
        return merged

    def _normalize_europepmc_record(self, record: Dict[str, Any], query: str) -> SourcePaper:
        title = self._clean_html_text(record.get("title", ""))
        doi = normalize_text(record.get("doi", ""))
        abstract = normalize_text(record.get("abstractText", ""))
        published_date = (
            normalize_text(record.get("firstPublicationDate", ""))
            or normalize_text(record.get("firstIndexDate", ""))
        )
        authors = self._normalize_europepmc_authors(
            record.get("authorList"), record.get("authorString", "")
        )
        landing_url = f"{self.landing_base}/{doi}" if doi else ""
        pdf_url = f"{self.landing_base}/{doi}.full.pdf" if doi else ""

        return SourcePaper(
            source=self.platform,
            source_id=doi or safe_filename(f"{self.platform}_{published_date}_{title}"),
            title=title,
            doi=doi,
            abstract=abstract,
            authors=authors,
            published_date=published_date,
            updated_date="",
            journal=self.journal_name,
            category="",
            landing_url=landing_url,
            pdf_url=pdf_url,
            query=query,
            version="",
            keywords=[],
            extra={
                "provider": "europepmc",
                "epmc_id": record.get("id", ""),
                "raw_record": record,
            },
        )

    def _normalize_europepmc_authors(self, author_list: Any, author_string: Any) -> List[str]:
        if isinstance(author_list, dict):
            authors = author_list.get("author") or []
            if isinstance(authors, list):
                names: List[str] = []
                for author in authors:
                    if not isinstance(author, dict):
                        continue
                    full = normalize_text(author.get("fullName", ""))
                    if full:
                        names.append(full)
                        continue
                    given = normalize_text(author.get("firstName", ""))
                    family = normalize_text(author.get("lastName", ""))
                    name = " ".join(part for part in [given, family] if part)
                    if name:
                        names.append(name)
                if names:
                    return names

        text = normalize_text(author_string)
        if text:
            return [part.strip() for part in re.split(r"[,;]\s*", text) if part.strip()]
        return []

    def _clean_html_text(self, text: Any) -> str:
        cleaned = normalize_text(text)
        if not cleaned:
            return ""
        soup = BeautifulSoup(cleaned, "html.parser")
        cleaned = soup.get_text(separator=" ")
        cleaned = re.sub(r"\s+", " ", cleaned).strip()
        return cleaned

    def _search_by_doi(self, doi: str) -> List[SourcePaper]:
        record = self._fetch_crossref_work(doi)
        if not record:
            return []
        if normalize_text(record.get("publisher", "")).lower() != "openrxiv":
            return []
        if self._platform_from_record(record) != self.platform:
            return []
        return [self._normalize_crossref_record(record, query=doi)]

    def _fetch_crossref_work(self, doi: str) -> Dict[str, Any]:
        url = f"{BIO_RXIV_CROSSREF_API}/{quote(doi, safe='')}"
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
        use_europepmc: bool = True,
    ) -> List[SourcePaper]:
        records = self.search(
            query=query,
            start_date=start_date,
            end_date=end_date,
            max_results=max_results,
            use_europepmc=use_europepmc,
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
            start_dt = self.launch_date

        if end_date:
            end_dt = datetime.strptime(end_date, "%Y-%m-%d")
        else:
            end_dt = datetime.utcnow()

        if start_dt < self.launch_date:
            start_dt = self.launch_date
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
                response = self._get_http_client().get(BIO_RXIV_CROSSREF_API, params=params)
                response.raise_for_status()
                return response.json()
            except Exception as exc:
                last_error = exc
                if attempt + 1 < self.max_retries:
                    time.sleep(min(2.0, 0.5 * (attempt + 1)))

        if last_error is not None:
            raise last_error
        raise RuntimeError("Failed to query bioRxiv Crossref search")

    def _platform_from_record(self, record: Dict[str, Any]) -> str:
        # The authoritative signal: openRxiv deposits a per-server primary URL.
        primary_url = normalize_text((record.get("resource") or {}).get("primary", {}).get("URL", ""))
        if "medrxiv.org" in primary_url:
            return "medrxiv"
        if "biorxiv.org" in primary_url:
            return "biorxiv"

        # Fallback: medRxiv accessions are 8 digits, bioRxiv accessions are 6.
        return detect_platform_from_doi(record.get("DOI", ""))

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
        landing_url = f"{self.landing_base}/{doi}" if doi else ""

        source_id = doi or safe_filename(f"{self.platform}_{published_date}_{title}")
        keywords = self._normalize_crossref_keywords(record.get("subject", []))
        if not keywords:
            group_title = normalize_text(record.get("group-title", ""))
            keywords = [group_title] if group_title else []
        pdf_url = f"{self.landing_base}/{doi}.full.pdf" if doi else ""

        return SourcePaper(
            source=self.platform,
            source_id=source_id,
            title=title,
            doi=doi,
            abstract=abstract,
            authors=authors,
            published_date=published_date,
            updated_date=updated_date,
            journal=self.journal_name,
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
        pdf_candidates = [
            candidate
            for candidate in [record.pdf_url, self._api_pdf_url(record), self._early_pdf_url(record)]
            if candidate
        ]

        for candidate in pdf_candidates:
            if download_binary(candidate, pdf_path, headers=self.headers, timeout=self.request_timeout):
                return True

        citation_pdf = self._scrape_citation_pdf(record.landing_url)
        if citation_pdf and download_binary(citation_pdf, pdf_path, headers=self.headers, timeout=self.request_timeout):
            return True

        # bioRxiv/medRxiv serve a Cloudflare challenge to non-browser clients
        # (HTTP 403). Two opt-in browser fallbacks can clear it, in increasing
        # order of strength: CloakBrowser (PAPER_FETCH_CLOAK), then
        # undetected_chromedriver (PAPER_FETCH_UNDETECTED) for hosts where
        # CloakBrowser's stealth still stalls at "Just a moment...".
        if is_cloak_enabled():
            for candidate in pdf_candidates:
                data = cloak_fetch_pdf(candidate, timeout=int(self.request_timeout))
                if data and data[:5].startswith(b"%PDF"):
                    try:
                        pdf_path.parent.mkdir(parents=True, exist_ok=True)
                        pdf_path.write_bytes(data)
                        return True
                    except OSError:
                        return False

        if is_undetected_enabled():
            for candidate in pdf_candidates:
                data = undetected_fetch_pdf(candidate, timeout=int(self.request_timeout))
                if data and data[:5].startswith(b"%PDF"):
                    try:
                        pdf_path.parent.mkdir(parents=True, exist_ok=True)
                        pdf_path.write_bytes(data)
                        return True
                    except OSError:
                        return False

        return False

    def _scrape_citation_pdf(self, landing_url: str) -> str:
        if not landing_url:
            return ""
        try:
            response = self._get_http_client().get(landing_url)
            response.raise_for_status()
            soup = BeautifulSoup(response.text, "html.parser")
            meta = soup.find("meta", {"name": "citation_pdf_url"})
            if meta and meta.get("content"):
                return normalize_text(meta.get("content"))
        except Exception:
            return ""
        return ""

    def _api_pdf_url(self, record: SourcePaper) -> str:
        # The bioRxiv/medRxiv details API returns the exact latest version, which
        # lets us build the versioned {doi}v{version}.full.pdf URL instead of
        # guessing. This endpoint is not behind the Cloudflare bot wall.
        if not record.doi:
            return ""
        try:
            response = self._get_http_client().get(
                f"https://api.biorxiv.org/details/{self.platform}/{record.doi}"
            )
            response.raise_for_status()
            collection = response.json().get("collection") or []
            if not collection:
                return ""
            latest = collection[-1]
            version = normalize_text(latest.get("version", "")) or "1"
            return f"{self.landing_base}/{record.doi}v{version}.full.pdf"
        except Exception:
            return ""

    def _early_pdf_url(self, record: SourcePaper) -> str:
        # bioRxiv/medRxiv serve preprints under a HighWire "early" path as an
        # alternative to the versioned {doi}.full.pdf route.
        published_date = normalize_text(record.published_date)
        match = re.search(r"(\d{4})-(\d{2})-(\d{2})", published_date)
        if not match or not record.doi:
            return ""
        year, month, day = match.groups()
        accession = record.doi.split("/", 1)[-1]
        return f"{self.landing_base}/{self.platform}/early/{year}/{month}/{day}/{accession}.full.pdf"
