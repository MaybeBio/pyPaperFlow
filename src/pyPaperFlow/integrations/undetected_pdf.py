#!/usr/bin/env python3

"""Fetch a Cloudflare-gated PDF (bioRxiv/medRxiv) via undetected_chromedriver.

Usage: undetected_pdf.py <url> [timeout_seconds]

Stdout: raw PDF bytes (binary).
Stderr: progress / error messages (safe to /dev/null).
Exit:   0 on success (bytes written), 1 on any failure.

This is the last-resort companion to biorxiv_fetcher._download_pdf(). bioRxiv
serves `*.full.pdf` behind Cloudflare bot management: plain httpx/curl get a
403/429 JS challenge, headless Chrome is fingerprint-detected ("Attention
Required"), and CloakBrowser's Playwright stealth still stalls at "Just a
moment...". The only route found to work (verified 2026-09-02) is
undetected_chromedriver running Chrome *headed* under Xvfb with a Linux X11
user-agent that matches the Chrome major version. The browser solves the
challenge on the origin, then a direct navigation to the PDF URL triggers a
download (via `plugins.always_open_pdf_externally` + CDP
`Page.setDownloadBehavior`); the helper reads the downloaded file and writes
its bytes to stdout. FlareSolverr's HTTP API cannot be used here because it
only returns `driver.page_source` (HTML), not the PDF bytes.

Environment:
  DISPLAY                 If unset, this helper starts Xvfb on a free display and
                          tears it down on exit.
  UNDETECTED_CHROME_PATH  Optional path to a Chrome/Chromium binary (used when it
                          is not on PATH or in a standard location).
  UNDETECTED_DRIVER_PATH  Optional path to a pre-installed chromedriver binary.
                          When unset, undetected_chromedriver downloads/caches a
                          matching driver itself.
"""

import ipaddress
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from urllib.parse import urlparse

# Match the fetch path's cap so an oversized response is rejected before we
# hand bytes back to the caller.
MAX_PDF_SIZE = 50 * 1024 * 1024  # 50 MB

# Loopback aliases + cloud metadata hostnames.
# KEEP IN SYNC with the identical _BLOCKED_HOSTS set in cloak_pdf.py — this
# companion runs as its own process and re-implements the SSRF gate rather than
# importing it.
_BLOCKED_HOSTS = {
    "localhost",
    "localhost.localdomain",
    "ip6-localhost",
    "ip6-loopback",
    "metadata.google.internal",
    "metadata.aws.internal",
    "metadata",
}

# Linux X11 platform token; the Chrome major version is filled in dynamically
# so the UA matches the installed browser (a Windows UA on a Linux browser, or a
# version mismatch, is itself a fingerprint signal Cloudflare keys on).
_BASE_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/{version}.0.0.0 Safari/537.36"
)

_xvfb_proc: "subprocess.Popen | None" = None


def _err(msg: str) -> None:
    print(f"[undetected] {msg}", file=sys.stderr)


def _url_is_safe(url: str) -> "tuple[bool, str]":
    """SSRF gate mirroring cloak_pdf.py's check.

    undetected_pdf runs as its own process and is never imported by the fetcher,
    so the guard is duplicated here to keep this companion independently safe
    when invoked directly: block non-http(s) schemes, non-80/443 ports, private
    / loopback / metadata hosts, and hostnames that resolve into private space.
    """
    try:
        p = urlparse(url)
    except Exception:
        return False, "malformed_url"
    if p.scheme not in ("http", "https"):
        return False, "scheme_not_allowed"
    if p.port is not None and p.port not in (80, 443):
        return False, "port_not_allowed"
    host = (p.hostname or "").lower()
    if not host:
        return False, "empty_host"
    if host in _BLOCKED_HOSTS:
        return False, "blocked_host"
    try:
        ip = ipaddress.ip_address(host)
        if ip.is_private or ip.is_loopback or ip.is_link_local or ip.is_reserved or ip.is_multicast or ip.is_unspecified:
            return False, "private_ip"
        return True, ""
    except ValueError:
        pass  # hostname is a name, not a literal — resolve it below
    try:
        infos = socket.getaddrinfo(host, None)
    except OSError:
        return False, "dns_error"
    for info in infos:
        try:
            ip = ipaddress.ip_address(info[4][0])
        except ValueError:
            continue
        if ip.is_private or ip.is_loopback or ip.is_link_local or ip.is_reserved or ip.is_multicast or ip.is_unspecified:
            return False, "private_ip"
    return True, ""


def _find_chrome() -> "str | None":
    """Locate a Chrome/Chromium binary so uc.Chrome does not guess."""
    for name in (
        "google-chrome",
        "google-chrome-stable",
        "chromium",
        "chromium-browser",
        "chrome",
    ):
        path = shutil.which(name)
        if path:
            return path
    for path in (
        "/usr/bin/google-chrome",
        "/usr/bin/google-chrome-stable",
        "/usr/bin/chromium",
        "/usr/bin/chromium-browser",
        "/opt/google/chrome/chrome",
        str(Path.home() / ".local" / "chrome" / "opt" / "google" / "chrome" / "chrome"),
    ):
        if os.path.isfile(path):
            return path
    return None


def _chrome_major_version(chrome: "str | None") -> "int | None":
    if not chrome:
        return None
    try:
        out = subprocess.run(
            [chrome, "--version"], capture_output=True, timeout=10
        )
        text = (out.stdout or b"") + (out.stderr or b"")
        m = re.search(rb"(\d+)\.", text)
        if m:
            return int(m.group(1))
    except Exception:
        pass
    return None


def _start_xvfb() -> None:
    """Start Xvfb on a free display when DISPLAY is unset; teardown on exit."""
    global _xvfb_proc
    if os.environ.get("DISPLAY"):
        return
    if not shutil.which("Xvfb"):
        _err("no DISPLAY and no Xvfb binary found; cannot run headed Chrome")
        return
    for display_num in range(99, 140):
        display = f":{display_num}"
        try:
            proc = subprocess.Popen(
                ["Xvfb", display, "-screen", "0", "1920x1080x24"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
        except Exception:
            return
        time.sleep(0.5)
        if proc.poll() is None:
            _xvfb_proc = proc
            os.environ["DISPLAY"] = display
            _err(f"started Xvfb on {display}")
            return


def _stop_xvfb() -> None:
    global _xvfb_proc
    if _xvfb_proc is not None:
        try:
            _xvfb_proc.terminate()
            _xvfb_proc.wait(timeout=5)
        except Exception:
            try:
                _xvfb_proc.kill()
            except Exception:
                pass
        _xvfb_proc = None


def main() -> int:
    if not (2 <= len(sys.argv) <= 3):
        _err("usage: undetected_pdf.py <url> [timeout_seconds]")
        return 1

    url = sys.argv[1]
    ok, reason = _url_is_safe(url)
    if not ok:
        _err(f"refusing unsafe url ({reason}): {url}")
        return 1
    timeout_s = int(sys.argv[2]) if len(sys.argv) == 3 else 60

    # A Clash-style proxy (http_proxy/https_proxy -> 127.0.0.1:PORT) inherited
    # from the parent routes selenium's localhost chromedriver NEW_SESSION
    # handshake through the proxy, which answers "502 Bad Gateway". Chrome
    # itself does not read these env vars (it goes direct or via
    # --proxy-server), so dropping them only un-proxies the Python-side
    # handshake and never the browser's own page fetches.
    for _k in (
        "http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY",
        "all_proxy", "ALL_PROXY", "no_proxy", "NO_PROXY",
    ):
        os.environ.pop(_k, None)

    _start_xvfb()

    try:
        import undetected_chromedriver as uc
    except ImportError as e:
        _err(f"undetected_chromedriver import failed: {e}")
        _err("install via: pip install undetected-chromedriver")
        _stop_xvfb()
        return 1

    dl_dir = Path(tempfile.mkdtemp(prefix="undetected_pdf_"))
    driver = None
    try:
        chrome = os.environ.get("UNDETECTED_CHROME_PATH", "").strip() or _find_chrome()
        major = _chrome_major_version(chrome)

        options = uc.ChromeOptions()
        if chrome:
            options.binary_location = chrome
        # FlareSolverr's hardened flag set, minus `--headless` (headless is
        # fingerprint-detected). Run headed under Xvfb instead.
        for flag in (
            "--no-sandbox",
            "--disable-setuid-sandbox",
            "--disable-dev-shm-usage",
            "--no-zygote",
            "--ignore-certificate-errors",
            "--ignore-ssl-errors",
            "--window-size=1920,1080",
            "--disable-search-engine-choice-screen",
        ):
            options.add_argument(flag)
        if major:
            options.add_argument(f"--user-agent={_BASE_UA.format(version=major)}")
        options.add_experimental_option(
            "prefs",
            {
                "download.default_directory": str(dl_dir),
                "download.prompt_for_download": False,
                "download.directory_upgrade": True,
                "plugins.always_open_pdf_externally": True,
            },
        )

        kwargs = {"options": options, "headless": False}
        if chrome:
            kwargs["browser_executable_path"] = chrome
        if major:
            kwargs["version_main"] = major
        driver_path = os.environ.get("UNDETECTED_DRIVER_PATH", "").strip()
        if driver_path:
            kwargs["driver_executable_path"] = driver_path

        _err("launching headed undetected_chromedriver under Xvfb")
        driver = uc.Chrome(**kwargs)
        driver.execute_cdp_cmd(
            "Page.setDownloadBehavior",
            {"behavior": "allow", "downloadPath": str(dl_dir)},
        )

        # Solve the Cloudflare challenge on the origin first so the
        # cf_clearance cookie is set for the domain before hitting the PDF URL.
        parsed = urlparse(url)
        origin = f"{parsed.scheme}://{parsed.netloc}/"
        _err(f"clearing challenge at {origin}")
        driver.get(origin)

        deadline = time.time() + min(timeout_s, 60)
        cleared = False
        while time.time() < deadline:
            try:
                title = driver.title or ""
            except Exception:
                title = ""
            if title and "Just a moment" not in title and not title.startswith("Loading"):
                cleared = True
                break
            time.sleep(1)
        if not cleared:
            _err("challenge did not clear within timeout")
            return 1
        _err(f"origin title: {(driver.title or '')[:60]!r}")

        _err(f"downloading {url}")
        before = set(dl_dir.glob("*.pdf"))
        try:
            driver.get(url)
        except Exception as e:
            _err(f"pdf navigation warning: {e}")

        # Poll for a fully-downloaded PDF (ignore .crdownload partials).
        pdf_file = None
        deadline = time.time() + timeout_s
        while time.time() < deadline:
            for candidate in dl_dir.glob("*.pdf"):
                if candidate in before:
                    continue
                if candidate.stat().st_size > 0:
                    pdf_file = candidate
                    break
            if pdf_file:
                break
            time.sleep(0.5)

        if not pdf_file:
            _err("no PDF file appeared in download directory")
            return 1

        body = pdf_file.read_bytes()
        if not body.startswith(b"%PDF"):
            _err("downloaded file is not a PDF")
            return 1
        if len(body) > MAX_PDF_SIZE:
            _err(f"response exceeds {MAX_PDF_SIZE} byte cap")
            return 1
        sys.stdout.buffer.write(body)
        sys.stdout.buffer.flush()
        _err(f"done, {len(body)} bytes")
        return 0
    except Exception as e:
        _err(f"failed: {e}")
        return 1
    finally:
        if driver is not None:
            try:
                driver.quit()
            except Exception:
                pass
        try:
            shutil.rmtree(dl_dir, ignore_errors=True)
        except Exception:
            pass
        _stop_xvfb()


if __name__ == "__main__":
    sys.exit(main())
