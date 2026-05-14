# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

pyPaperFlow is an automated paper reading platform that fetches and processes scientific literature from multiple sources (PubMed, arXiv, bioRxiv). It implements a 7-stage workflow for literature retrieval, processing, and knowledge extraction.

## Common Commands

### Installation
```bash
pip install -e .
```

### Running Tests

No Python test files currently exist; the `test/` directory contains only data fixtures.

### CLI Usage (paperflow command)

**PubMed Operations:**
```bash
# Search PubMed for PMIDs only
paperflow pubmed-search "query terms" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Fetch metadata from query
paperflow pubmed-meta --query "query terms" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Fetch metadata from PMID list file
paperflow pubmed-meta --file pmids.txt --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Download full text (PMC) for PMIDs
paperflow pubmed-content --file pmids.txt --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Fetch both metadata and full text
paperflow pubmed-all --query "query terms" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Merge papers to canonical JSON/JSONL
paperflow pubmed-merge-json --input ./papers_dir --output ./output.jsonl

# Export merged JSON to single Markdown file
paperflow pubmed-export-md --input ./merged.jsonl --output ./output.md --config config.yaml
```

**arXiv Operations:**
```bash
paperflow arxiv-search "query terms" --max-results 10 --output-dir ./output
paperflow arxiv-fetch "query terms" --max-results 10 --download-pdf --output-dir ./output
paperflow arxiv-fetch "query terms" --start-date 2024-01-01 --end-date 2024-12-31 --output-dir ./output
paperflow arxiv-fetch "query terms" --backend paperscraper --output-dir ./output
```

**bioRxiv Operations:**
```bash
paperflow biorxiv-search "query terms" --max-results 10 --output-dir ./output
paperflow biorxiv-fetch "query terms" --max-results 10 --download-pdf --output-dir ./output
```

**Paper Fetch (DOI → PDF):**
```bash
# Fetch PDF by DOI via Unpaywall, Semantic Scholar, arXiv, PMC, bioRxiv, Sci-Hub
paperflow paper-fetch 10.1038/s41586-020-2649-2 -o ./papers
paperflow paper-fetch --batch dois.txt -o ./papers --format text
paperflow paper-fetch --title "AlphaFold" -o ./papers
# Set UNPAYWALL_EMAIL env var for best results
```

**MinerU Integration:**
```bash
# Parse a PDF using MinerU engine
paperflow pdf-parse -i paper.pdf -o ./output --clear

# Parse MinerU content_list_v2.json into canonical sectioned JSON
paperflow mineru-parse -i content_list_v2.json -o paper.json
```

## Architecture

### Core Data Flow
```
Raw Data (PubMed/arXiv/bioRxiv) -> Fetcher -> Paper Object -> JSON/Database -> Merge -> Export
PDF -> MinerU -> content_list_v2.json -> MinerUContentParser -> Structured JSON
DOI/Title -> paper-fetch engine -> PDF file on disk
```

### Storage Structure
```
output_dir/
├── pubmed/ (or arxiv/, biorxiv/)
│   └── {year}/
│       └── {source_id}/
│           ├── {source_id}.json           # Metadata
│           ├── {source_id}.xml            # Full text XML (PubMed only)
│           ├── {source_id}_parsed.json    # Parsed full text structure
│           ├── {source_id}_parsed.md      # Markdown full text
│           └── {source_id}.pdf            # Downloaded PDF (arXiv/bioRxiv)
```

### Key Modules

**Fetchers:**
- `src/pyPaperFlow/pubmed/pubmed_fetcher.py` — PubMed metadata and full text via NCBI Entrez API. Defines the PubMed paper data model (PaperIdentity, PaperContent, PaperContributors, PaperSource, PaperLinks, PaperMetadata, Paper_MetaData, Paper_TextData). Uses Biopython's Entrez; requires email, optionally api-key for 10 req/s rate limit. Fetches both Medline and XML formats, performs batch ELink queries for linked data.
- `src/pyPaperFlow/arxiv_fetcher.py` — arXiv metadata and PDF fetching (native httpx or paperscraper backend). Handles rate limiting with Retry-After support.
- `src/pyPaperFlow/biorxiv_fetcher.py` — bioRxiv metadata and PDF via Crossref API (filters publisher="openrxiv", prefix="10.64898"). Cursor-based pagination.

**Merger:**
- `src/pyPaperFlow/pubmed/pubmed_merger.py` — Two-stage pipeline: (1) `merge_json_from_directory` creates a canonical merged JSON/JSONL per paper, (2) `export_md_from_merged_json` produces a single Markdown file from the merged JSON, driven by an optional YAML config that selects metadata fields and content sections.

**Integrations:**
- `src/pyPaperFlow/integrations/pdf_fetch.py` — Paper-fetch engine that resolves DOIs/titles to PDFs via a chain: Unpaywall → Semantic Scholar → arXiv → Europe PMC → PMC → bioRxiv/medRxiv → publisher-direct (institutional mode) → Sci-Hub (last resort). Supports batch, dry-run, idempotency, and streaming output. Exit codes discriminate validation, transport, and unresolved failures.
- `src/pyPaperFlow/integrations/mineru_parser.py` — Parses MinerU's `content_list_v2.json` into canonical sectioned JSON. Extracts metadata (title, authors, year, DOI, journal), builds hierarchical sections matched to canonical types (abstract, introduction, methods, results, discussion, etc.) via regex patterns, and preserves tables as HTML. See `CANONICAL_ORDER` list for all section types.

**Data Models:**
- `src/pyPaperFlow/source_models.py` — Generic `SourcePaper` model for arXiv/bioRxiv papers.
- `src/pyPaperFlow/source_utils.py` — Common utilities (normalize_text, safe_filename, build_source_record_dir, download_binary).
- `src/pyPaperFlow/utils.py` — PubMed-specific utilities (URL extraction from text).

**CLI:**
- `src/pyPaperFlow/cli.py` — Typer-based CLI. Commands: `pubmed-search`, `pubmed-meta`, `pubmed-content`, `pubmed-all`, `pubmed-merge-json`, `pubmed-export-md`, `arxiv-search`, `arxiv-fetch`, `biorxiv-search`, `biorxiv-fetch`, `paper-fetch` (passthrough to pdf_fetch engine), `pdf-parse` (MinerU wrapper), `mineru-parse` (content_list_v2 → structured JSON).

### Important Implementation Details

**PubMed Fetcher (`PubmedFetcher`):**
- Uses Biopython's Entrez API for all PubMed operations
- Requires `--email` and optionally `--api-key` (increases rate limit to 10 req/s)
- Fetches both Medline and XML formats for comprehensive metadata
- Performs batch ELink queries for linked data (citations, references, similar articles, PMC links)
- Parses PMC XML into hierarchical JSON structure and Markdown

**arXiv Fetcher (`ArxivFetcher`):**
- Two backends: `native` (httpx) or `paperscraper` (optional dependency)
- Builds normalized queries with date filters
- Handles rate limiting with Retry-After header support

**bioRxiv Fetcher (`BioRxivFetcher`):**
- Uses Crossref API for server-side query (not the legacy bioRxiv details API)
- `--window-days` is a compatibility-only option and is ignored

**Paper Merger (`PubmedMerger`):**
- Two-stage pipeline: merge JSON first, then export to Markdown
- Merge stage produces canonical JSON per paper with metadata + parsed full text
- Export stage driven by YAML config selecting metadata fields and content sections

**Paper-Fetch Engine (`pdf_fetch`):**
- DOI resolution chain: Unpaywall → Semantic Scholar → arXiv → Europe PMC → PMC → bioRxiv/medRxiv → publisher-direct → Sci-Hub
- Title resolution via Crossref → Semantic Scholar fallback chain
- SSRF protection (blocks private IPs, non-http schemes, non-80/443 ports, cloud metadata hosts)
- Machine-readable NDJSON progress on stderr when `--format json`; prose when `--format text`
- Exit codes: 0=success, 1=unresolved, 3=validation error, 4=transport error

**MinerU Parser (`MinerUContentParser`):**
- Extracts metadata by scanning page headers/footers for journal name, year, DOI
- Author paragraph detection heuristics (comma-rich paragraph after level-1 title)
- Section matching: strips leading section numbers, then regex-matches against canonical types
- Handles subsections (level >= 3 or multi-part numbering like "2.1.")
- Aggregates multiple sections of the same canonical type into one, storing later occurrences as subsections

### NCBI API Requirements
- **Email**: Required by NCBI Entrez API (use `--email` parameter)
- **API Key**: Optional but recommended (increases rate limit from 3 to 10 req/s)
- **Rate Limiting**: All fetchers implement retry logic with exponential backoff
- **Batch Processing**: Default batch size is 50-100 records per request

### Date Handling
- PubMed: Uses DP (Date of Publication) field, extracts year (first 4 characters)
- arXiv: Uses submittedDate filter with YYYYMMDDHHmm format
- bioRxiv: Uses Crossref from-pub-date/until-pub-date filters (YYYY-MM-DD)

### Link Extraction
- PubMed automatically extracts URLs from abstract text during metadata parsing
- URLs are categorized (GitHub, GitLab, Zenodo, Figshare, HuggingFace, Google Drive)
- Additional URLs can be extracted from full text and merged into metadata

# CLAUDE.md

Behavioral guidelines to reduce common LLM coding mistakes. Merge with project-specific instructions as needed.

**Tradeoff:** These guidelines bias toward caution over speed. For trivial tasks, use judgment.

## 1. Think Before Coding

**Don't assume. Don't hide confusion. Surface tradeoffs.**

Before implementing:
- State your assumptions explicitly. If uncertain, ask.
- If multiple interpretations exist, present them - don't pick silently.
- If a simpler approach exists, say so. Push back when warranted.
- If something is unclear, stop. Name what's confusing. Ask.

## 2. Simplicity First

**Minimum code that solves the problem. Nothing speculative.**

- No features beyond what was asked.
- No abstractions for single-use code.
- No "flexibility" or "configurability" that wasn't requested.
- No error handling for impossible scenarios.
- If you write 200 lines and it could be 50, rewrite it.

Ask yourself: "Would a senior engineer say this is overcomplicated?" If yes, simplify.

## 3. Surgical Changes

**Touch only what you must. Clean up only your own mess.**

When editing existing code:
- Don't "improve" adjacent code, comments, or formatting.
- Don't refactor things that aren't broken.
- Match existing style, even if you'd do it differently.
- If you notice unrelated dead code, mention it - don't delete it.

When your changes create orphans:
- Remove imports/variables/functions that YOUR changes made unused.
- Don't remove pre-existing dead code unless asked.

The test: Every changed line should trace directly to the user's request.

## 4. Goal-Driven Execution

**Define success criteria. Loop until verified.**

Transform tasks into verifiable goals:
- "Add validation" → "Write tests for invalid inputs, then make them pass"
- "Fix the bug" → "Write a test that reproduces it, then make it pass"
- "Refactor X" → "Ensure tests pass before and after"

For multi-step tasks, state a brief plan:
```
1. [Step] → verify: [check]
2. [Step] → verify: [check]
3. [Step] → verify: [check]
```

Strong success criteria let you loop independently. Weak criteria ("make it work") require constant clarification.

---

**These guidelines are working if:** fewer unnecessary changes in diffs, fewer rewrites due to overcomplication, and clarifying questions come before implementation rather than after mistakes.
