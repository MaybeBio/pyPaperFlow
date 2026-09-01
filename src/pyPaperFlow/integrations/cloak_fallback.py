from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path


def is_cloak_enabled() -> bool:
    """True iff the operator opted into the CloakBrowser fallback (PAPER_FETCH_CLOAK=1)."""
    return bool(os.environ.get("PAPER_FETCH_CLOAK"))


def resolve_cloak_python() -> str | None:
    """Locate a Python interpreter that can import cloakbrowser.

    Order: CLOAKBROWSER_PYTHON env -> ~/github/CloakBrowser/.venv/bin/python ->
    the current interpreter. Mirrors pdf_fetch.py's resolution.
    Returns None when no candidate can import cloakbrowser.
    """
    candidates = [
        os.environ.get("CLOAKBROWSER_PYTHON", "").strip(),
        str(Path.home() / "github" / "CloakBrowser" / ".venv" / "bin" / "python"),
        sys.executable,
    ]
    for candidate in candidates:
        if not candidate:
            continue
        if not (os.path.isfile(candidate) or shutil.which(candidate)):
            continue
        try:
            result = subprocess.run(
                [candidate, "-c", "import cloakbrowser"],
                capture_output=True,
                timeout=30,
            )
            if result.returncode == 0:
                return candidate
        except Exception:
            continue
    return None


def cloak_fetch_pdf(url: str, *, timeout: int = 60) -> bytes | None:
    """Fetch PDF bytes through CloakBrowser. Returns bytes, or None on failure.

    Shells out to the companion cloak_pdf.py via a cloakbrowser-importable
    Python. Fails closed: a missing dependency or any error returns None so the
    caller falls through to its next candidate.
    """
    if not is_cloak_enabled():
        return None
    python = resolve_cloak_python()
    if not python:
        return None
    helper = Path(__file__).with_name("cloak_pdf.py")
    if not helper.exists():
        return None
    try:
        # Browser launch + challenge solve is slow; give it headroom over the
        # per-request timeout.
        result = subprocess.run(
            [python, str(helper), url, str(timeout)],
            capture_output=True,
            timeout=timeout + 90,
        )
    except Exception:
        return None
    if result.returncode != 0 or not result.stdout:
        return None
    return result.stdout
