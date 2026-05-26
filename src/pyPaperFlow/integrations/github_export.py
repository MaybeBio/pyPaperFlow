import csv
import json
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import *
from urllib.parse import urlsplit, urlunsplit

import httpx


def _load_merged_papers_any(input_path: str) -> List[Dict[str, Any]]:
    """Load merged papers from JSON dict or JSONL list file."""
    path = Path(input_path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    suffix = path.suffix.lower()
    if suffix == ".jsonl":
        papers: List[Dict[str, Any]] = []
        with path.open("r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                obj = json.loads(line)
                if isinstance(obj, dict):
                    papers.append(obj)
        return papers

    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    if isinstance(data, dict):
        return [paper for paper in data.values() if isinstance(paper, dict)]
    if isinstance(data, list):
        return [paper for paper in data if isinstance(paper, dict)]
    raise ValueError("Unsupported merged file format. Expect JSON dict/list or JSONL.")


def _get_pmid(paper: Dict[str, Any]) -> str:
    pmid = paper.get("pmid")
    if pmid:
        return str(pmid)
    pmid = paper.get("meta", {}).get("identity", {}).get("pmid")
    return str(pmid) if pmid else "N/A"


def _extract_github_entries(papers: List[Dict[str, Any]]) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for paper in papers:
        pmid = _get_pmid(paper)
        links = paper.get("meta", {}).get("links", {})
        text_mined = links.get("text_mined", [])
        if not isinstance(text_mined, list):
            continue
        for item in text_mined:
            if not isinstance(item, dict):
                continue
            category = str(item.get("category", "")).strip().lower()
            if category != "github":
                continue
            url = str(item.get("url", "")).strip()
            if not url:
                continue
            rows.append({"pmid": pmid, "original_url": url})
    return rows


def _ensure_https(url: str) -> str:
    """Prepend ``https://`` to protocol-less GitHub URLs."""
    url = url.strip()
    if url.startswith("github.com/"):
        url = "https://" + url
    return url


def _normalize_url(url: str) -> str:
    # Strip wrapping punctuation first so _ensure_https can match bare
    # hostnames like '"github.com/a/b"'.
    url = url.strip().strip("\"'[]()<>{}")
    url = _ensure_https(url)
    while url and url[-1] in ".,;:!?":
        url = url[:-1]
    if " " in url:
        return url
    try:
        parsed = urlsplit(url)
    except Exception:
        return url
    host = parsed.netloc.lower()
    # Unicode dashes / curly quotes in the hostname signal data corruption.
    if host and any(ch in host for ch in "‐‑‒–—―‘’“”"):
        return url
    return url


def _remove_trailing_number_noise(url: str) -> str:
    """Remove terminal numeric noise from last path segment.

    Strips trailing digits and any preceding ``_`` / ``-`` separator.
    Used as a fallback -- the original URL is always tried first.
    """
    try:
        parsed = urlsplit(url)
    except Exception:
        return url

    segments = [seg for seg in parsed.path.split("/") if seg]
    if not segments:
        return url

    last = segments[-1]
    m = re.match(r"^(.*?\D)(\d+)$", last)
    if not m:
        return url

    cleaned = m.group(1).rstrip("_-")
    if not cleaned:
        segments.pop()
    else:
        segments[-1] = cleaned
    if not segments:
        return url
    new_path = "/" + "/".join(segments)
    return urlunsplit((parsed.scheme, parsed.netloc, new_path, parsed.query, parsed.fragment))


def _is_accessible_url(url: str, timeout: float = 10.0, retries: int = 1) -> Tuple[bool, Optional[int], str]:
    headers = {
        "User-Agent": "pyPaperFlow-github-export/0.1",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
    }
    last_status: Optional[int] = None
    last_effective = url
    for _ in range(max(1, retries + 1)):
        try:
            resp = httpx.head(url, headers=headers, timeout=timeout, follow_redirects=True)
            last_status = resp.status_code
            last_effective = str(resp.url) or url
            if resp.status_code < 400:
                return True, resp.status_code, last_effective
        except Exception:
            continue
    return False, last_status, last_effective


def _extract_repo_slug(url: str) -> Optional[str]:
    try:
        parsed = urlsplit(url)
    except Exception:
        return None
    host = parsed.netloc.lower()
    if "github.com" not in host:
        return None

    segments = [seg for seg in parsed.path.split("/") if seg]
    if len(segments) < 2:
        return None
    owner, repo = segments[0], segments[1]
    if repo.endswith(".git"):
        repo = repo[:-4]
    if not owner or not repo:
        return None
    return f"{owner}/{repo}"


def _write_report_csv(report_path: str, rows: List[Dict[str, Any]]) -> None:
    Path(report_path).parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pmid",
        "original_url",
        "modified_url",
        "effective_url",
        "used_modified",
        "accessible",
        "http_status",
        "repo_slug",
    ]
    with open(report_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def _append_repo_section(output_md: str, separator: str, slug: str, source_url: str, content: str) -> None:
    Path(output_md).parent.mkdir(parents=True, exist_ok=True)
    header = (
        f"\n{separator}\n"
        f"### REPO: {slug}\n"
        f"- source_url: {source_url}\n"
    )
    with open(output_md, "a", encoding="utf-8") as fh:
        fh.write(header)
        fh.write(content.strip())
        fh.write("\n")


def run_github_export(
    input_json: str,
    output_md: str,
    report: Optional[str] = None,
    separator: str = "<<<PY_PAPERFLOW_REPO_BOUNDARY>>>",
    timeout: float = 10.0,
    retries: int = 3,
    sleep_sec: float = 0.2,
    max_repos: Optional[int] = None,
    strict_ghresearcher: bool = False,
    overwrite: bool = False,
    include_tree: bool = True,
    echo: Optional[Callable[[str], None]] = None,
) -> Dict[str, Any]:
    log = echo or (lambda msg: None)

    papers = _load_merged_papers_any(input_json)
    entries = _extract_github_entries(papers)
    if not entries:
        return {
            "total_url_entries": 0,
            "unique_urls": 0,
            "accessible_entries": 0,
            "accessible_ratio": 0.0,
            "report": report or f"{os.path.splitext(output_md)[0]}_github_report.csv",
            "selected_repos": 0,
            "exported": 0,
            "failed": 0,
            "output": output_md,
            "ghresearcher_skipped": True,
        }

    log(f"Found {len(entries)} GitHub URL entries from {len(papers)} papers.")
    unique_urls = sorted({_normalize_url(e["original_url"]) for e in entries})
    validation_cache: Dict[str, Dict[str, Any]] = {}

    for raw_url in unique_urls:
        ok, status, effective = _is_accessible_url(raw_url, timeout=timeout, retries=retries)
        modified_url = ""
        used_modified = False
        check_url = raw_url

        if not ok:
            fallback = _remove_trailing_number_noise(raw_url)
            if fallback != raw_url:
                ok2, status2, effective2 = _is_accessible_url(fallback, timeout=timeout, retries=retries)
                modified_url = fallback
                if ok2:
                    ok = True
                    status = status2
                    effective = effective2
                    check_url = fallback
                    used_modified = True

        repo_slug = _extract_repo_slug(check_url if ok else raw_url)
        validation_cache[raw_url] = {
            "modified_url": modified_url,
            "used_modified": used_modified,
            "accessible": ok,
            "http_status": status,
            "effective_url": effective,
            "repo_slug": repo_slug,
        }

    report_rows: List[Dict[str, Any]] = []
    for entry in entries:
        normalized = _normalize_url(entry["original_url"])
        meta = validation_cache.get(normalized, {})
        report_rows.append(
            {
                "pmid": entry["pmid"],
                "original_url": entry["original_url"],
                "modified_url": meta.get("modified_url", ""),
                "effective_url": meta.get("effective_url", ""),
                "used_modified": meta.get("used_modified", False),
                "accessible": meta.get("accessible", False),
                "http_status": meta.get("http_status", ""),
                "repo_slug": meta.get("repo_slug", ""),
            }
        )

    accessible_rows = [r for r in report_rows if r["accessible"]]
    total = len(report_rows)
    accessible_count = len(accessible_rows)
    ratio = (accessible_count / total) if total else 0.0

    report_path = report or f"{os.path.splitext(output_md)[0]}_github_report.csv"
    _write_report_csv(report_path, report_rows)

    slug_to_source: Dict[str, str] = {}
    for row in accessible_rows:
        slug = row.get("repo_slug")
        if not slug:
            continue
        if slug not in slug_to_source:
            source = row.get("modified_url") if row.get("used_modified") else row.get("original_url")
            slug_to_source[slug] = str(source)

    slugs = sorted(slug_to_source.keys())
    if max_repos is not None:
        slugs = slugs[: max(0, max_repos)]

    ghresearcher_skipped = False
    if slugs:
        if include_tree:
            # ---- ghresearcher path (with file tree) ----
            ghresearcher_path = shutil.which("ghresearcher")
            if not ghresearcher_path:
                if strict_ghresearcher:
                    raise RuntimeError("ghresearcher command not found in PATH")
                ghresearcher_skipped = True
                slugs = []
            else:
                if overwrite:
                    Path(output_md).parent.mkdir(parents=True, exist_ok=True)
                    with open(output_md, "w", encoding="utf-8") as fh:
                        fh.write("")

                success = 0
                failed = 0
                for idx, slug in enumerate(slugs, start=1):
                    log(f"[{idx}/{len(slugs)}] ghresearcher parse {slug} --view --clear")
                    proc = subprocess.run(
                        [ghresearcher_path, "parse", slug, "--view", "--clear"],
                        capture_output=True,
                        text=True,
                    )
                    if proc.returncode != 0:
                        failed += 1
                        if strict_ghresearcher:
                            err = proc.stderr.strip() or "unknown error"
                            raise RuntimeError(f"ghresearcher failed for {slug}: {err}")
                        continue
                    _append_repo_section(
                        output_md=output_md,
                        separator=separator,
                        slug=slug,
                        source_url=slug_to_source.get(slug, ""),
                        content=proc.stdout,
                    )
                    success += 1
                    if sleep_sec > 0:
                        time.sleep(sleep_sec)

                return {
                    "total_url_entries": total,
                    "unique_urls": len(unique_urls),
                    "accessible_entries": accessible_count,
                    "accessible_ratio": ratio,
                    "report": report_path,
                    "selected_repos": len(slugs),
                    "exported": success,
                    "failed": failed,
                    "output": output_md,
                    "ghresearcher_skipped": False,
                }
        else:
            # ---- gh repo view path (no file tree) ----
            gh_path = shutil.which("gh")
            if not gh_path:
                if strict_ghresearcher:
                    raise RuntimeError("gh (GitHub CLI) command not found in PATH")
                ghresearcher_skipped = True
                slugs = []
            else:
                if overwrite:
                    Path(output_md).parent.mkdir(parents=True, exist_ok=True)
                    with open(output_md, "w", encoding="utf-8") as fh:
                        fh.write("")

                success = 0
                failed = 0
                for idx, slug in enumerate(slugs, start=1):
                    log(f"[{idx}/{len(slugs)}] gh repo view {slug}")
                    proc = subprocess.run(
                        [gh_path, "repo", "view", slug],
                        capture_output=True,
                        text=True,
                    )
                    if proc.returncode != 0:
                        failed += 1
                        if strict_ghresearcher:
                            err = proc.stderr.strip() or "unknown error"
                            raise RuntimeError(f"gh repo view failed for {slug}: {err}")
                        continue
                    _append_repo_section(
                        output_md=output_md,
                        separator=separator,
                        slug=slug,
                        source_url=slug_to_source.get(slug, ""),
                        content=proc.stdout,
                    )
                    success += 1
                    if sleep_sec > 0:
                        time.sleep(sleep_sec)

                return {
                    "total_url_entries": total,
                    "unique_urls": len(unique_urls),
                    "accessible_entries": accessible_count,
                    "accessible_ratio": ratio,
                    "report": report_path,
                    "selected_repos": len(slugs),
                    "exported": success,
                    "failed": failed,
                    "output": output_md,
                    "ghresearcher_skipped": False,
                }

    return {
        "total_url_entries": total,
        "unique_urls": len(unique_urls),
        "accessible_entries": accessible_count,
        "accessible_ratio": ratio,
        "report": report_path,
        "selected_repos": 0,
        "exported": 0,
        "failed": 0,
        "output": output_md,
        "ghresearcher_skipped": ghresearcher_skipped,
    }
