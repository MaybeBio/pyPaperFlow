# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

pyPaperFlow is an automated paper reading platform that fetches and processes scientific literature from multiple sources (PubMed, arXiv, bioRxiv). It implements a 7-stage workflow for literature retrieval, processing, and knowledge extraction, with current implementation covering Stages 1-2 and parts of Stages 4-5.

## Common Commands

### Installation
```bash
pip install -e .
```

### Running Tests
```bash
# Run all tests
python -m pytest test/

# Run specific test file
python -m pytest test/test_arxiv_fetcher.py
python -m pytest test/test_biorxiv_fetcher.py
python -m pytest test/test_merger.py

# Run with verbose output
python -m pytest test/ -v
```

### CLI Usage (paperflow command)

**PubMed Operations:**
```bash
# Search PubMed for PMIDs only
paperflow pubmed-search "query terms" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Fetch metadata from query
paperflow pubmed-fetch --query "query terms" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Fetch metadata from PMID list file
paperflow pubmed-fetch --file pmids.txt --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Download full text (PMC) for PMIDs
paperflow download-fulltext --file pmids.txt --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output

# Fetch both metadata and full text
paperflow fetch-full --query "query terms" --email YOUR_EMAIL --api-key YOUR_NCBI_API_KEY --output-dir ./output
```

**arXiv Operations:**
```bash
# Search arXiv for IDs
paperflow arxiv-search "query terms" --max-results 10 --output-dir ./output

# Fetch arXiv papers with optional PDF download
paperflow arxiv-fetch "query terms" --max-results 10 --download-pdf --output-dir ./output

# With date filtering
paperflow arxiv-fetch "query terms" --start-date 2024-01-01 --end-date 2024-12-31 --output-dir ./output

# Using paperscraper backend (requires optional dependency)
paperflow arxiv-fetch "query terms" --backend paperscraper --output-dir ./output
```

**bioRxiv Operations:**
```bash
# Search bioRxiv for IDs (uses Crossref server-side query)
paperflow biorxiv-search "query terms" --max-results 10 --output-dir ./output

# Fetch bioRxiv papers with optional PDF download
paperflow biorxiv-fetch "query terms" --max-results 10 --download-pdf --output-dir ./output
```

**Merge Operations (PubMed only):**
```bash
# Merge papers to Markdown (metadata + full text)
paperflow merge ./papers_dir ./output.md --mode full --format md

# Merge with specific PMIDs from file
paperflow merge ./papers_dir ./output.md --pmid-file pmids.txt --mode full

# Merge to JSONL with LLM profile
paperflow merge ./papers_dir ./output.jsonl --mode full --format jsonl --profile llm

# Merge specific sections only
paperflow merge ./papers_dir ./output.md --mode full --include-sections abstract,introduction,results
```

## Architecture

### Core Data Flow
```
Raw Data (PubMed/arXiv/bioRxiv) -> Fetcher -> Paper Object -> JSON/Database -> AI Analyzer
```

### Storage Structure
Papers are stored in a hierarchical structure:
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
- `src/pyPaperFlow/pubmed/pubmed_fetcher.py` - PubMed metadata and full text fetching via NCBI Entrez API
- `src/pyPaperFlow/arxiv_fetcher.py` - arXiv metadata and PDF fetching (native or paperscraper backend)
- `src/pyPaperFlow/biorxiv_fetcher.py` - bioRxiv metadata and PDF fetching via Crossref API

**Merger:**
- `src/pyPaperFlow/pubmed/pubmed_merger.py` - Merges PubMed papers into unified formats (Markdown, JSONL, plain text) for AI analysis

**Data Models:**
- `src/pyPaperFlow/pubmed/pubmed_fetcher.py` defines PubMed paper structure:
  - `PaperIdentity` - PMID, DOI, title
  - `PaperContent` - abstract, keywords, MeSH terms, publication types
  - `PaperContributors` - authors and affiliations (medline + xml formats)
  - `PaperSource` - journal and publication details
  - `PaperLinks` - citations, references, PMC links, external links
  - `PaperMetadata` - fetch timestamps
  - `Paper_MetaData` - combines identity, content, contributors, source, metadata, links
  - `Paper_TextData` - full text XML, parsed JSON, parsed Markdown

- `src/pyPaperFlow/source_models.py` defines generic source paper structure:
  - `SourcePaper` - unified model for arXiv/bioRxiv papers

**Utilities:**
- `src/pyPaperFlow/source_utils.py` - Common utilities (normalize_text, safe_filename, build_source_record_dir, download_binary)
- `src/pyPaperFlow/utils.py` - PubMed-specific utilities (URL extraction from text)

**CLI:**
- `src/pyPaperFlow/cli.py` - Typer-based CLI with commands for all operations

### Important Implementation Details

**PubMed Fetcher (`PubmedFetcher`):**
- Uses Biopython's Entrez API for all PubMed operations
- Requires `--email` and optionally `--api-key` (increases rate limit to 10 req/s)
- Fetches both Medline and XML formats for comprehensive metadata
- Performs batch ELink queries for linked data (citations, references, similar articles, PMC links)
- Supports both query-based and PMID-list-based fetching
- Parses PMC XML into hierarchical JSON structure and Markdown

**arXiv Fetcher (`ArxivFetcher`):**
- Two backends: `native` (httpx) or `paperscraper` (optional dependency)
- Builds normalized queries with date filters
- Handles rate limiting (429 responses) with Retry-After header support
- Downloads PDFs when available

**bioRxiv Fetcher (`BioRxivFetcher`):**
- Uses Crossref API for server-side query over openRxiv records (not the legacy bioRxiv details API)
- Filters by publisher="openrxiv" and prefix="10.64898"
- Supports cursor-based pagination for large result sets

**Paper Merger (`PaperMerger`):**
- Three modes: metadata only, metadata + content, full
- Three output formats: Markdown, JSONL, plain text
- Two output profiles: analysis (full metadata), LLM (compact for AI consumption)
- Text source priority: parsed JSON first or Markdown first
- Handles both directory scan and PMID file list input

**Testing Patterns:**
- Uses `unittest` framework
- Tests use `tempfile.TemporaryDirectory()` for isolated file operations
- Mock HTTP clients for testing fetchers without real API calls
- Test fixtures build realistic directory structures with JSON/parsed files

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
