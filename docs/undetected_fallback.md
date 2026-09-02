# bioRxiv/medRxiv PDF — undetected_chromedriver fallback

Last resort for downloading `www.biorxiv.org/content/{doi}.full.pdf` (and the
medRxiv equivalent). These PDF endpoints sit behind Cloudflare bot management,
which defeats every ordinary client. This document records *why*, *what finally
worked*, and *how to reproduce/debug it* on a fresh host.

## Why this exists (the Cloudflare wall)

| Approach | Result |
| --- | --- |
| Plain `httpx` / `curl` (direct) | HTTP 429 JS challenge |
| Plain `httpx` / `curl` (via Clash proxy) | HTTP 403 |
| Headless Chrome (`--headless`) | fingerprint-detected → "Attention Required" |
| CloakBrowser (Playwright stealth), headless *and* headed under Xvfb | stalls at "Just a moment..." |
| Reuse `cf_clearance` cookie via `httpx` | Cloudflare binds the cookie to the browser's TLS/JA3 fingerprint |
| `uc.Chrome(headless=True)` | headless is still fingerprint-detected |
| **`undetected_chromedriver` + Xvfb *headed* mode** | **works** (verified 2026-09-02) |

The one working recipe (borrowed from FlareSolverr's engine) is:

- `undetected_chromedriver` launching **real Chrome headed** under Xvfb
  (`headless=False`, NOT `--headless` — headless is itself a fingerprint signal),
- a Linux X11 `--user-agent` whose Chrome major version **matches the installed
  browser** (a Windows UA or a version mismatch is another signal Cloudflare
  keys on),
- FlareSolverr's hardening flags (`--no-sandbox --disable-dev-shm-usage
  --no-zygote --ignore-certificate-errors` etc.).

The browser solves the Cloudflare challenge on the origin, then a direct
navigation to the PDF URL triggers a download (via
`plugins.always_open_pdf_externally` + CDP `Page.setDownloadBehavior`); the
helper reads the downloaded file and writes its bytes to stdout.

> **FlareSolverr itself cannot be used as the download service.** Its HTTP API
> returns only `driver.page_source` (HTML), never the PDF bytes. It can only be
> used as a pre-built "Chrome + uc + Xvfb" *environment*, which is exactly the
> stack reproduced here natively — so running it natively is strictly better.

## Architecture

Three pieces, mirroring the existing CloakBrowser fallback:

- `src/pyPaperFlow/integrations/undetected_pdf.py` — standalone helper script,
  run as its **own process** (never imported). Takes `<url> [timeout_seconds]`,
  writes raw PDF bytes to stdout. Re-implements the SSRF gate from
  `cloak_pdf.py` and strips proxy env vars (see "Proxy leak" below).
- `src/pyPaperFlow/integrations/undetected_fallback.py` — the in-process
  wrapper: `is_undetected_enabled()`, `resolve_undetected_python()`, and
  `undetected_fetch_pdf(url, *, timeout)`. Fails closed to `None` on any
  missing dependency.
- `src/pyPaperFlow/preprint/biorxiv_fetcher.py` → `_download_pdf()` — the
  download chain is: ordinary binary download → CloakBrowser (if
  `PAPER_FETCH_CLOAK=1`) → undetected_chromedriver (if
  `PAPER_FETCH_UNDETECTED=1`) → give up.

## Install (host, no sudo)

Tested with Python 3.14 + selenium 4.43 + `undetected-chromedriver` 3.5.5.

```bash
# 1. Chrome — extract a .deb into a user dir (no root needed)
mkdir -p ~/.local/chrome
wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
dpkg -x google-chrome-stable_current_amd64.deb ~/.local/chrome
# binary: ~/.local/chrome/opt/google/chrome/chrome

# 2. chromedriver — MUST match the Chrome major version.
#    Download from Chrome for Testing:
#    https://googlechromelabs.github.io/chrome-for-testing/
#    (known-good pair: Chrome 152.0.7977.75 + chromedriver 152.0.7977.75)
unzip chromedriver-linux64.zip -d ~/.local/bin/
chmod +x ~/.local/bin/chromedriver

# 3. undetected_chromedriver (in the same env that runs paperflow)
pip install undetected-chromedriver

# 4. Xvfb (Ubuntu/Debian; often already present)
sudo apt install xvfb
```

## Environment variables

| Variable | Meaning | Required |
| --- | --- | --- |
| `PAPER_FETCH_UNDETECTED=1` | Enable the fallback (off by default) | **yes** |
| `UNDETECTED_CHROME_PATH` | Path to the Chrome binary | optional (auto-detected) |
| `UNDETECTED_DRIVER_PATH` | Path to the chromedriver binary | optional (uc auto-downloads) |
| `UNDETECTED_PYTHON` | Interpreter that can `import undetected_chromedriver` | optional (defaults to `sys.executable`) |

Minimal persistent setup (append to `~/.zshrc`):

```bash
export PAPER_FETCH_UNDETECTED=1
export UNDETECTED_CHROME_PATH="$HOME/.local/chrome/opt/google/chrome/chrome"
export UNDETECTED_DRIVER_PATH="$HOME/.local/bin/chromedriver"
```

## Verify it works

```bash
PAPER_FETCH_UNDETECTED=1 \
UNDETECTED_CHROME_PATH=~/.local/chrome/opt/google/chrome/chrome \
UNDETECTED_DRIVER_PATH=~/.local/bin/chromedriver \
  paperflow biorxiv-fetch --doi 10.1101/2021.10.04.463034 -o /tmp/e2e

# expect: /tmp/e2e/biorxiv/<year>/<doi>/<doi>.pdf exists and starts with %PDF
head -c 8 /tmp/e2e/biorxiv/*/*/*.pdf
```

First launch includes a ~25 s+ Chrome cold start; each subsequent paper is
~19 s.

## Debugging

### `[undetected] failed: Message: Bad Gateway` at launch

**Root cause: proxy env leak.** A Clash-style proxy
(`http_proxy`/`https_proxy` → `127.0.0.1:7892`) is inherited by the
`undetected_pdf.py` subprocess and routes selenium's `127.0.0.1:<port>/session`
NEW_SESSION handshake through the proxy, which answers `502 Bad Gateway`.

Chrome itself **does not read** `http_proxy`/`https_proxy` (it goes direct or
via `--proxy-server`), so dropping those vars only un-proxies the Python-side
handshake — never the browser's page fetches. The fix is already applied at the
top of `main()` in `undetected_pdf.py` (`os.environ.pop(...)` for all proxy
vars).

Quick confirmation that this is the problem:

```bash
env -u http_proxy -u https_proxy -u HTTP_PROXY -u HTTPS_PROXY \
  python src/pyPaperFlow/integrations/undetected_pdf.py "<pdf url>" 120
```

### `Missing X server or $DISPLAY`

`undetected_pdf.py` starts its own Xvfb on a free display when `DISPLAY` is
unset. If you see this from chromedriver, either Xvfb is missing (`apt install
xvfb`) or a stale Xvfb/chromedriver from a previous run is holding the display —
kill strays and retry.

### Challenge never clears / "Just a moment..." forever

- Confirm the Chrome major version matches the chromedriver major version
  (mismatch = the browser never actually launches).
- Confirm headed mode is in effect — `headless=True` is fingerprint-detected.
- Raise the helper's `timeout_seconds` (Chrome 152 cold start is slow; the
  default `request_timeout` of 60 s is fine, but very cold hosts may need more).

### Related proxy note (Europe PMC)

The same Clash proxy times out on `ebi.ac.uk` (HTTP 504), which is why the
Europe PMC client uses `httpx.Client(trust_env=False)`. If Europe PMC search
"isn't working", check the proxy before touching the code.
