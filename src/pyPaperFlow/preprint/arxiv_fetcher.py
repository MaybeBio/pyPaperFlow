from __future__ import annotations

import re
import importlib
import time
import xml.etree.ElementTree as ET
from datetime import date, datetime, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import httpx

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

PAPERSCRAPER_FIELDS = [
    "title",
    "authors",
    "date",
    "abstract",
    "journal",
    "doi",
    "entry_id",
    "pdf_url",
    "updated",
    "primary_category",
    "comment",
]


def _safe_lower(text: Any) -> str:
    return normalize_text(text).lower()


class ArxivFetcher:
    def __init__(
        self,
        root_dir: str,
        backend: str = "native",
        batch_size: int = 100,
        max_retries: int = 3,
        request_timeout: float = 60.0,
    ):
        self.root_dir = root_dir
        if backend not in {"native", "paperscraper"}:
            raise ValueError("backend must be 'native' or 'paperscraper'")
        self.backend = backend
        self.batch_size = max(1, int(batch_size))
        self.max_retries = max(1, int(max_retries))
        self.request_timeout = float(request_timeout)
        self.headers = {
            "User-Agent": "pyPaperFlow/0.1.0 (+https://github.com/MaybeBio/pyPaperFlow)",
            "Accept": "application/atom+xml,application/xml;q=0.9,*/*;q=0.8",
        }
        self._http_client: Optional[httpx.Client] = None

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
            start_bound, end_bound = self._normalize_date_bounds(start_date, end_date)
            start_text = self._format_arxiv_date(start_bound.isoformat() if start_bound else "1991-01-01", suffix="0000")
            end_text = self._format_arxiv_date(end_bound.isoformat() if end_bound else datetime.now(timezone.utc).strftime("%Y-%m-%d"), suffix="2359")
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
        if self.backend == "paperscraper":
            return self._search_with_paperscraper(
                search_query=search_query,
                query=query,
                max_results=max_results,
            )
        return self._search_native(
            search_query=search_query,
            query=query,
            max_results=max_results,
            start_date=start_date,
            end_date=end_date,
        )

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

    def fetch_from_ids(
        self,
        ids: Iterable[str],
        output_dir: Optional[str] = None,
        download_pdf: bool = True,
    ) -> List[SourcePaper]:
        records: List[SourcePaper] = []
        for raw_id in ids:
            arxiv_id = normalize_text(raw_id)
            if not arxiv_id:
                continue
            # id_list does not accept a version suffix; drop a trailing vN.
            clean_id = re.sub(r"v\d+$", "", arxiv_id)
            feed = self._request_id_feed(clean_id)
            entries = feed.findall("atom:entry", ARXIV_ATOM_NS)
            if not entries:
                continue
            record = self._normalize_entry(entries[0], query=arxiv_id)
            self._save_record(record, output_dir=output_dir, download_pdf=download_pdf)
            records.append(record)
        return records

    def close(self) -> None:
        if self._http_client is not None:
            self._http_client.close()
            self._http_client = None

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def _search_native(
        self,
        search_query: str,
        query: str,
        max_results: int,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
    ) -> List[SourcePaper]:
        results: List[SourcePaper] = []
        start = 0
        remaining = max(1, int(max_results))
        tried_fallback = False

        while remaining > 0:
            page_size = min(self.batch_size, remaining)
            feed = self._request_feed(search_query, start=start, max_results=page_size)
            entries = feed.findall("atom:entry", ARXIV_ATOM_NS)
            # If there are no entries on the first page, try some fallback query forms
            if not entries and start == 0 and not tried_fallback:
                tried_fallback = True
                raw_query = normalize_text(query)
                date_filter = ""
                if ") AND " in search_query and "submittedDate:" in search_query:
                    try:
                        _, date_filter = search_query.split(") AND ", 1)
                        date_filter = " AND " + date_filter
                    except Exception:
                        date_filter = ""

                stopwords = {"a", "an", "the", "for", "of", "in", "on", "and", "or", "to", "from", "by", "with"}
                tokens = [t for t in re.findall(r'"[^"]+"|\'[^\']+\'|\S+', raw_query) if t]

                alt_queries = []
                alt_queries.append(raw_query)
                alt_queries.append(f'all:"{raw_query}"')
                # remove common stopwords and re-build tokenized form
                cleaned = []
                for t in tokens:
                    tclean = t.strip('"\'')
                    if not tclean:
                        continue
                    low = tclean.lower()
                    if low in stopwords:
                        continue
                    if " " in tclean:
                        cleaned.append(f'all:"{tclean}"')
                    else:
                        cleaned.append(f'all:{tclean}')
                if cleaned:
                    alt_queries.append(" AND ".join(cleaned))

                found = False
                for alt in alt_queries:
                    try_query = alt + date_filter
                    feed = self._request_feed(try_query, start=0, max_results=page_size)
                    entries = feed.findall("atom:entry", ARXIV_ATOM_NS)
                    if entries:
                        found = True
                        break
                if not found:
                    break

            for entry in entries:
                results.append(self._normalize_entry(entry, query=query))
                remaining -= 1
                if remaining <= 0:
                    break

            if len(entries) < page_size or remaining <= 0:
                break

            start += len(entries)

        return results

    def _search_with_paperscraper(self, search_query: str, query: str, max_results: int) -> List[SourcePaper]:
        api = self._load_paperscraper_api()
        result = api(
            search_query,
            fields=PAPERSCRAPER_FIELDS,
            max_results=max(1, int(max_results)),
            verbose=False,
        )
        rows = self._records_from_paperscraper_result(result)

        records: List[SourcePaper] = []
        for row in rows[: max(1, int(max_results))]:
            records.append(self._normalize_paperscraper_row(row, query=query))
        return records

    def _records_from_paperscraper_result(self, result: Any) -> List[Dict[str, Any]]:
        if result is None:
            return []
        if isinstance(result, list):
            return [row for row in result if isinstance(row, dict)]
        if hasattr(result, "to_dict"):
            try:
                records = result.to_dict(orient="records")
            except TypeError:
                records = result.to_dict()
            if isinstance(records, list):
                return [row for row in records if isinstance(row, dict)]
        return []

    def _request_feed(self, search_query: str, start: int, max_results: int) -> ET.Element:
        params = {
            "search_query": search_query,
            "start": start,
            "max_results": max_results,
            "sortBy": "submittedDate",
            "sortOrder": "descending",
        }
        return self._request_feed_params(
            params,
            description=f"query={search_query!r} start={start} max_results={max_results}",
        )

    def _request_id_feed(self, arxiv_id: str) -> ET.Element:
        params = {
            "id_list": arxiv_id,
            "start": 0,
            "max_results": 1,
        }
        return self._request_feed_params(params, description=f"id={arxiv_id!r}")

    def _request_feed_params(self, params: Dict[str, Any], description: str) -> ET.Element:
        last_error: Optional[Exception] = None

        for attempt in range(self.max_retries):
            response: Optional[httpx.Response] = None
            try:
                response = self._get_http_client().get(
                    ARXIV_API_URL,
                    params=params,
                    headers=self.headers,
                    timeout=self.request_timeout,
                )
                if response.status_code == 429:
                    last_error = RuntimeError(f"arXiv API rate limited request for {description}")
                    if attempt + 1 < self.max_retries:
                        self._sleep_before_retry(response, attempt)
                        continue
                    break
                response.raise_for_status()
                try:
                    return ET.fromstring(response.text)
                except ET.ParseError as exc:
                    last_error = exc
                    if attempt + 1 < self.max_retries:
                        self._sleep_before_retry(response, attempt)
                        continue
                    break
            except (httpx.HTTPStatusError, httpx.TimeoutException, httpx.TransportError) as exc:
                last_error = exc
                if attempt + 1 < self.max_retries:
                    self._sleep_before_retry(response, attempt)
                    continue
                break
            except Exception as exc:
                last_error = exc
                if attempt + 1 < self.max_retries:
                    self._sleep_before_retry(response, attempt)
                    continue
                break

        if last_error is not None:
            raise last_error
        raise RuntimeError("Failed to query arXiv API")

    def _get_http_client(self) -> httpx.Client:
        if self._http_client is not None:
            return self._http_client

        client_kwargs = {
            "headers": self.headers,
            "timeout": self.request_timeout,
            "follow_redirects": True,
        }
        try:
            self._http_client = httpx.Client(http2=True, **client_kwargs)
        except ImportError:
            self._http_client = httpx.Client(http2=False, **client_kwargs)
        return self._http_client

    def _sleep_before_retry(self, response: Optional[httpx.Response], attempt: int) -> None:
        retry_after = self._retry_after_seconds(response)
        delay = retry_after if retry_after is not None else min(30.0, 1.5 * (2**attempt))
        time.sleep(max(0.0, delay))

    def _retry_after_seconds(self, response: Optional[httpx.Response]) -> Optional[float]:
        if response is None:
            return None

        raw_retry_after = normalize_text(response.headers.get("Retry-After", ""))
        if not raw_retry_after:
            return None

        if raw_retry_after.isdigit():
            return float(raw_retry_after)

        try:
            retry_after_dt = parsedate_to_datetime(raw_retry_after)
        except (TypeError, ValueError, IndexError):
            return None

        if retry_after_dt.tzinfo is None:
            retry_after_dt = retry_after_dt.replace(tzinfo=timezone.utc)
        now = datetime.now(retry_after_dt.tzinfo)
        return max(0.0, (retry_after_dt - now).total_seconds())

    def _normalize_date_bounds(
        self,
        start_date: Optional[str],
        end_date: Optional[str],
    ) -> tuple[Optional[date], Optional[date]]:
        start_bound = self._parse_date(start_date) if start_date else None
        end_bound = self._parse_date(end_date) if end_date else None
        today = datetime.now(timezone.utc).date()

        if end_bound and end_bound > today:
            end_bound = today
        if start_bound and start_bound > today and end_bound is None:
            end_bound = today
        if start_bound and end_bound and start_bound > end_bound:
            raise ValueError("start_date must not be after end_date after clamping future dates")

        return start_bound, end_bound

    def _parse_date(self, date_text: str) -> date:
        return datetime.strptime(date_text, "%Y-%m-%d").date()

    def _load_paperscraper_api(self):
        try:
            module = importlib.import_module("paperscraper.arxiv.arxiv")
        except ImportError as exc:
            raise ModuleNotFoundError(
                "paperscraper backend requested but the 'paperscraper' package is not installed"
            ) from exc
        return module.get_arxiv_papers_api

    def _normalize_paperscraper_row(self, row: Dict[str, Any], query: str) -> SourcePaper:
        title = normalize_text(row.get("title", ""))
        abstract = normalize_text(row.get("abstract", ""))
        journal = normalize_text(row.get("journal", ""))
        doi = normalize_text(row.get("doi", ""))
        entry_id = normalize_text(row.get("entry_id", ""))
        pdf_url = normalize_text(row.get("pdf_url", ""))
        source_id = self._extract_source_id(entry_id=entry_id, doi=doi, pdf_url=pdf_url, title=title or query)
        authors = self._normalize_authors(row.get("authors"))
        published_date = self._normalize_date_value(row.get("date"))
        updated_date = self._normalize_date_value(row.get("updated"))
        category = normalize_text(row.get("primary_category", "") or row.get("category", ""))

        if not pdf_url and source_id:
            pdf_url = f"https://arxiv.org/pdf/{source_id}.pdf"

        if not entry_id and source_id:
            entry_id = f"https://arxiv.org/abs/{source_id}"

        return SourcePaper(
            source="arxiv",
            source_id=source_id,
            title=title,
            doi=doi,
            abstract=abstract,
            authors=authors,
            published_date=published_date,
            updated_date=updated_date,
            journal=journal,
            category=category,
            landing_url=entry_id,
            pdf_url=pdf_url,
            query=query,
            version=source_id.rsplit("v", 1)[-1] if source_id and "v" in source_id else "",
            keywords=[item for item in [category] if item],
            extra={
                "backend": "paperscraper",
                "entry_id": entry_id,
                "primary_category": category,
                "comment": normalize_text(row.get("comment", "")),
            },
        )

    def _normalize_authors(self, value: Any) -> List[str]:
        if isinstance(value, list):
            authors: List[str] = []
            for item in value:
                if isinstance(item, str):
                    name = normalize_text(item)
                else:
                    name = normalize_text(getattr(item, "name", item))
                if name:
                    authors.append(name)
            return authors
        text = normalize_text(value)
        if not text:
            return []
        return [item.strip() for item in re.split(r"\s*,\s*", text) if item.strip()]

    def _normalize_date_value(self, value: Any) -> str:
        text = normalize_text(value)
        if not text:
            return ""
        if len(text) >= 10:
            return text[:10]
        return text

    def _extract_source_id(self, entry_id: str, doi: str, pdf_url: str, title: str) -> str:
        if entry_id:
            tail = entry_id.rstrip("/").rsplit("/", 1)[-1]
            if tail:
                return tail

        if doi:
            doi_match = re.search(r"10\.48550/arXiv\.(?P<source_id>[^\s]+)", doi, flags=re.IGNORECASE)
            if doi_match:
                return doi_match.group("source_id")

        if pdf_url:
            pdf_match = re.search(r"/pdf/(?P<source_id>[^/?#]+?)(?:\.pdf)?(?:[?#].*)?$", pdf_url, flags=re.IGNORECASE)
            if pdf_match:
                return pdf_match.group("source_id")

        return safe_filename(title or "arxiv-paper")

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