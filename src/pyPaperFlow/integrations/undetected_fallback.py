from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def is_undetected_enabled() -> bool:
    """True iff the operator opted into the undetected_chromedriver fallback.

    Controlled by PAPER_FETCH_UNDETECTED=1. This is the last-resort route for
    bioRxiv/medRxiv PDFs that are Cloudflare-walled even against CloakBrowser:
    undetected_chromedriver + Xvfb headed Chrome clears the challenge where
    headless/stealth browsers are fingerprint-detected.
    """
    return bool(os.environ.get("PAPER_FETCH_UNDETECTED"))


def resolve_undetected_python() -> str | None:
    """Locate a Python interpreter that can import undetected_chromedriver.

    Order: UNDETECTED_PYTHON env -> the current interpreter. Returns None when
    no candidate can import undetected_chromedriver (the helper then fails closed).
    """
    candidates = [
        os.environ.get("UNDETECTED_PYTHON", "").strip(),
        sys.executable,
    ]
    for candidate in candidates:
        if not candidate:
            continue
        if not (os.path.isfile(candidate) or shutil.which(candidate)):
            continue
        try:
            result = subprocess.run(
                [candidate, "-c", "import undetected_chromedriver"],
                capture_output=True,
                timeout=30,
            )
            if result.returncode == 0:
                return candidate
        except Exception:
            continue
    return None


def undetected_fetch_pdf(url: str, *, timeout: int = 60) -> bytes | None:
    """Fetch PDF bytes through undetected_chromedriver. Returns bytes, or None on failure.

    Shells out to the companion undetected_pdf.py via an
    undetected_chromedriver-importable Python. Fails closed: a missing
    dependency, missing Chrome/Xvfb, or any error returns None so the caller
    falls through (or gives up) as if the fallback were absent.
    """
    if not is_undetected_enabled():
        return None
    python = resolve_undetected_python()
    if not python:
        return None
    helper = Path(__file__).with_name("undetected_pdf.py")
    if not helper.exists():
        return None
    try:
        # Browser launch + challenge solve is slow (and Xvfb may need to spin
        # up); give it generous headroom over the per-request timeout.
        result = subprocess.run(
            [python, str(helper), url, str(timeout)],
            capture_output=True,
            timeout=timeout + 120,
        )
    except Exception:
        return None
    if result.returncode != 0 or not result.stdout:
        return None
    return result.stdout
