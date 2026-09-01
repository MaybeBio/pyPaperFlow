from __future__ import annotations

import json
import re
from datetime import datetime, timedelta
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
