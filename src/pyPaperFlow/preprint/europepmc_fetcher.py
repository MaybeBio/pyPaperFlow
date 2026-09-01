from __future__ import annotations

import re
import time
from typing import Any, Dict, List, Optional

import httpx

EUROPE_PMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPE_PMC_PREPRINT_FILTER = "SRC:PPR"

_BOOLEAN_TOKEN_RE = re.compile(r'"[^"]+"|\'[^\']+\'|\S+')
_BOOLEAN_OPERATOR_RE = re.compile(r"\b(AND|OR|NOT)\b", re.IGNORECASE)


class EuropePMCSearch:
    """Query Europe PMC's REST API for preprint records.

    Europe PMC indexes bioRxiv/medRxiv preprints (source ``PPR``) with full
    text and text-mined terms, so it can run strict boolean AND searches that
    Crossref's metadata-only relevance search cannot.

    The client bypasses the local HTTP proxy via ``trust_env=False`` because
    the proxy commonly times out (HTTP 504) on ``ebi.ac.uk``.
    """

    def __init__(
        self,
        request_timeout: float = 60.0,
        max_retries: int = 3,
        page_size: int = 100,
    ):
        self.request_timeout = float(request_timeout)
        self.max_retries = max(1, int(max_retries))
        self.page_size = max(1, min(1000, int(page_size)))
        self.headers = {
            "User-Agent": "pyPaperFlow/0.1.0 (+https://github.com/MaybeBio/pyPaperFlow)",
            "Accept": "application/json",
        }
        self._client = httpx.Client(
            headers=self.headers,
            timeout=self.request_timeout,
            follow_redirects=True,
            trust_env=False,
        )

    def close(self) -> None:
        self._client.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def build_query(
        self,
        query: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
    ) -> str:
        raw = (query or "").strip()
        if not raw:
            raise ValueError("query must be non-empty")

        if _BOOLEAN_OPERATOR_RE.search(raw):
            # User wrote explicit AND/OR/NOT; pass through unchanged.
            term_query = raw
        else:
            # Bare query: AND the terms explicitly (do not rely on the
            # server's default operator).
            tokens = [t for t in _BOOLEAN_TOKEN_RE.findall(raw) if t]
            terms = [t for t in tokens if t.upper() not in {"AND", "OR", "NOT"}]
            term_query = " AND ".join(terms)

        parts = [EUROPE_PMC_PREPRINT_FILTER]
        if term_query:
            parts.append(f"({term_query})")

        if start_date or end_date:
            start = start_date or "0000-00-00"
            end = end_date or "9999-12-31"
            parts.append(f"FIRST_PDATE:[{start} TO {end}]")

        return " AND ".join(parts)

    def search(
        self,
        query: str,
        start_date: Optional[str] = None,
        end_date: Optional[str] = None,
        max_results: Optional[int] = None,
    ) -> List[Dict[str, Any]]:
        search_query = self.build_query(query, start_date=start_date, end_date=end_date)
        records: List[Dict[str, Any]] = []
        cursor_mark = "*"
        page_size = self.page_size

        while True:
            if max_results is not None:
                remaining = max(1, int(max_results) - len(records))
                page_size = min(self.page_size, remaining)

            payload = self._request_page(search_query, cursor_mark, page_size)
            result_list = payload.get("resultList") or {}
            items = result_list.get("result") or []
            records.extend(items)

            if max_results is not None and len(records) >= int(max_results):
                break

            next_cursor = payload.get("nextCursorMark") or ""
            if not items or not next_cursor or next_cursor == cursor_mark:
                break
            cursor_mark = next_cursor

        return records

    def _request_page(self, query: str, cursor_mark: str, page_size: int) -> Dict[str, Any]:
        params = {
            "query": query,
            "format": "json",
            "pageSize": page_size,
            "cursorMark": cursor_mark,
            "resultType": "core",
        }
        last_error: Optional[Exception] = None
        for attempt in range(self.max_retries):
            try:
                response = self._client.get(EUROPE_PMC_SEARCH_URL, params=params)
                response.raise_for_status()
                return response.json()
            except Exception as exc:
                last_error = exc
                if attempt + 1 < self.max_retries:
                    time.sleep(min(2.0, 0.5 * (attempt + 1)))
        if last_error is not None:
            raise last_error
        raise RuntimeError("Failed to query Europe PMC")
