from __future__ import annotations

import json
import re
import time
from datetime import datetime, timedelta, timezone
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import httpx


BOOLEAN_OR_SPLIT_RE = re.compile(r"\s+OR\s+", re.IGNORECASE)
BOOLEAN_AND_SPLIT_RE = re.compile(r"\s+AND\s+", re.IGNORECASE)
TOKEN_RE = re.compile(r'"([^"]+)"|\'([^\']+)\'|(\S+)')


def normalize_text(value: Any) -> str:
    if value is None:
        return ""
    normalized = str(value).replace("\n", " ").strip()
    normalized = re.sub(r"\s+", " ", normalized)
    return normalized


def safe_filename(value: Any, fallback: str = "record") -> str:
    text = normalize_text(value)
    if not text:
        text = fallback
    text = text.replace("/", "_")
    text = text.replace("\\", "_")
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("._-")
    return text or fallback


def detect_platform_from_doi(doi: Any) -> str:
    """Infer 'biorxiv' vs 'medrxiv' from a bioRxiv/medRxiv DOI accession.

    medRxiv accessions are 8 digits; bioRxiv accessions are 6 digits.
    Returns '' when the DOI does not look like a bioRxiv/medRxiv accession.
    """
    doi_text = normalize_text(doi)
    if not doi_text:
        return ""
    suffix = doi_text.split("/", 1)[-1]
    match = re.search(r"(?:\d{4}\.\d{2}\.\d{2}\.)?(\d+)$", suffix)
    if not match:
        return ""
    digits = match.group(1)
    if len(digits) == 8:
        return "medrxiv"
    if len(digits) == 6:
        return "biorxiv"
    return ""


def extract_year(date_text: Any) -> str:
    text = normalize_text(date_text)
    if not text:
        return "unknown"

    for fmt in ("%Y-%m-%d", "%Y/%m/%d", "%Y-%m", "%Y/%m"):
        try:
            return datetime.strptime(text, fmt).strftime("%Y")
        except Exception:
            continue

    match = re.search(r"(19|20)\d{2}", text)
    if match:
        return match.group(0)
    return "unknown"


def ensure_directory(path: Path | str) -> Path:
    directory = Path(path)
    directory.mkdir(parents=True, exist_ok=True)
    return directory


def build_source_record_dir(base_dir: str | Path, source: str, year: str, source_id: str) -> Path:
    directory = Path(base_dir) / safe_filename(source) / safe_filename(year, fallback="unknown") / safe_filename(source_id)
    return ensure_directory(directory)


def save_json(path: Path | str, payload: Dict[str, Any]) -> None:
    output_path = Path(path)
    ensure_directory(output_path.parent)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)


def download_binary(url: str, output_path: Path | str, headers: Optional[Dict[str, str]] = None, timeout: float = 60.0) -> bool:
    try:
        response = httpx.get(url, headers=headers, timeout=timeout, follow_redirects=True)
        response.raise_for_status()
        content = response.content or b""
        if not content.startswith(b"%PDF"):
            return False

        output_file = Path(output_path)
        ensure_directory(output_file.parent)
        with output_file.open("wb") as handle:
            handle.write(content)
        return True
    except Exception:
        return False


def retry_after_seconds(response: Optional[httpx.Response]) -> Optional[float]:
    """Seconds to wait per the Retry-After header, or None when unusable."""
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


# Default retry budget for every preprint fetcher. This is the interactive/tool
# default: fail fast (~4.5s of backoff) so a human gets a quick answer plus the
# degradation notice, instead of a ~22s silent stall. Unattended callers (e.g.
# monitor.py) pass a larger max_retries explicitly.
DEFAULT_MAX_RETRIES = 3


def sleep_before_retry(response: Optional[httpx.Response], attempt: int) -> None:
    """Back off between attempts, capped at 30s.

    A 503 usually asks for a pause via Retry-After; otherwise back off
    exponentially. A short linear cap cannot outlast a transient outage, which
    pushes callers into whatever lossy fallback they have.
    """
    retry_after = retry_after_seconds(response)
    delay = retry_after if retry_after is not None else min(30.0, 1.5 * (2**attempt))
    time.sleep(max(0.0, delay))


def parse_boolean_query(query: str) -> List[List[str]]:
    text = normalize_text(query)
    if not text:
        return []

    clauses: List[List[str]] = []
    for or_clause in BOOLEAN_OR_SPLIT_RE.split(text):
        and_parts = BOOLEAN_AND_SPLIT_RE.split(or_clause)
        terms: List[str] = []
        for part in and_parts:
            for token_match in TOKEN_RE.findall(part):
                token = next((group for group in token_match if group), "")
                token = normalize_text(token).strip('"\'')
                if not token:
                    continue
                if token.upper() in {"AND", "OR", "NOT"}:
                    continue
                terms.append(token.lower())
        if terms:
            clauses.append(terms)
    return clauses


def basic_boolean_text_match(text: Any, query: str) -> bool:
    haystack = normalize_text(text).lower()
    clauses = parse_boolean_query(query)
    if not clauses:
        return True

    for clause in clauses:
        if all(term in haystack for term in clause):
            return True
    return False


def iter_date_windows(start_date: datetime, end_date: datetime, window_days: int) -> Iterable[tuple[datetime, datetime]]:
    if window_days < 1:
        raise ValueError(f"window_days must be >= 1, got {window_days}")

    current_start = start_date
    while current_start <= end_date:
        current_end = min(current_start + timedelta(days=window_days - 1), end_date)
        yield current_start, current_end
        current_start = current_end + timedelta(days=1)
