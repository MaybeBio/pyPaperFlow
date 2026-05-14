import subprocess
import typer
import os
import json
from pathlib import Path
from typing import *
from .pubmed.pubmed_fetcher import PubmedFetcher
from .preprint.arxiv_fetcher import ArxivFetcher
from .preprint.biorxiv_fetcher import BioRxivFetcher
from .pubmed.pubmed_merger import PubmedMerger
from .integrations import pdf_fetch
from datetime import datetime

app = typer.Typer(help="pyPaperFlow CLI", no_args_is_help=True)

# Common options
opt_storage = typer.Option("./Papers", "--storage-dir", "-s", help="Directory in Repository-level to store paper data for Initialization.") # Note: this is a repository-level default path
opt_email = typer.Option(..., "--email", help="Entrez Email.")
opt_api_key = typer.Option(None, "--api-key", help="NCBI API Key (recommended).")
opt_max_retries = typer.Option(3, "--max-retries", help="Maximum number of retries for Entrez API calls.")
opt_batch_size = typer.Option(50, "--batch-size", "-b", help="Batch size for fetching.")
opt_arxiv_backend = typer.Option("native", "--backend", help="arXiv backend: 'native' or 'paperscraper'.")


def _save_id_list(output_dir: str, filename: str, values: List[str]) -> str:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, filename)
    with open(output_file, "w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")
    return output_file


#############################################################
#  1, For Pubmed Parser
#############################################################

@app.command("pubmed-search")
def search_cmd(
    query: str = typer.Argument(..., help="PubMed search query."), # Argument means that this parameter is different from option, it must be provided and does not need flag --query like option!
    retmax: int = typer.Option(500, "--retmax", "-n", help="Max number of PMIDs to return every batch, must less than 10000."),
    email: str = opt_email,
    api_key: Optional[str] = opt_api_key,
    storage_dir: str = opt_storage, # Repository-level default path, used for initialing fetcher
    output_dir: Optional[str] = typer.Option(None, "--output-dir", "-o", help="Directory in result-level to store output IDs."), # User-specified output path in result-level
    max_retries: int = opt_max_retries
):
    """
    Search PubMed using Your customized query and return PMIDs.

    \b
    Notes:
    - 1, This command only searches and returns PMIDs, it does not fetch paper metadata.
    - 2, This command will print the found PMIDs and also save them to 'pubmed_searched_ids.txt' in the specified output directory. 
    If --output-dir is not specified, it will default to the storage directory.
    - 3, Note that storage_dir is used to initialize the fetcher for consistency, while output_dir is where the PMIDs are saved. They are different parameters!

    \b
    Example usage:
    - 1. Search for papers related to "machine learning" and return up to 500 PMIDs/per batch:
    paperflow pubmed-search "machine learning" --retmax 500 --output-dir ./MyPapers --email "YOUR_EMAIL@example.com" --api-key "YOUR_NCBI_API_KEY"
    """
    # we initialize fetcher with storage_dir for consistency
    fetcher = PubmedFetcher(root_dir=storage_dir, entrez_email=email, api_key=api_key or "", max_retries=max_retries)
    
    # 1. Search
    query_meta = fetcher.query_search(query)
    
    # 2. Get PMIDs
    pmids = fetcher.get_pubmedIDs_from_query(query_meta, retmax=retmax)
    
    typer.echo(f"Found {len(pmids)} PMIDs.")
    typer.echo(pmids)

    # 3. Optionally, save PMIDs to a file
    save_dir = output_dir if output_dir else storage_dir
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
    # note: search log should be saved with timestamp
    # pmid_file = os.path.join(save_dir, "pubmed_searched_ids.txt") 
    
    pmid_file = f"{save_dir}/pubmed_searched_ids_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.txt"
    with open(pmid_file, 'w') as f:
        for pmid in pmids:
            f.write(f"{pmid}\n")
    typer.echo(f"PMIDs saved to {pmid_file}.")

@app.command("pubmed-meta")
def fetch_cmd(
    query: Optional[str] = typer.Option(None, "--query", "-q", help="PubMed search query."),
    file: Optional[str] = typer.Option(None, "--file", "-f", help="Text file containing PMIDs (one per line), -q and -f are mutually exclusive."),
    batch_size: int = opt_batch_size,
    email: str = opt_email,
    api_key: Optional[str] = opt_api_key,
    storage_dir: str = opt_storage,
    max_retries: int = opt_max_retries,
    output_dir: Optional[str] = typer.Option(".", "--output-dir", "-o", help="Directory in result-level to store output papers, default is current directory. If not specified, will be set to root directory of the repository-level which is storage_dir. 🌟 We will create a '/pubmed' subfolder under the output directory to save all pubmed related data"),
):
    """
    Fetch paper metadata from PubMed using Your customized query, pmid list file and save to storage.

    \b
    Notes:
    - 1, You must provide one of --query, or --file to specify which papers to fetch. Note that they are mutually exclusive.
    - 2, -f can be used to fetch one or more PMIDs listed in a text file (one PMID per line).

    \b
    Example usage:
    - 1. Fetch papers for a query and save to storage:
      paperflow pubmed-fetch --query "machine learning" --output-dir ./MyPapers --email "YOUR_EMAIL@example.com" --api-key "YOUR_NCBI_API_KEY"
    - 2. Fetch papers from a list of PMIDs in a file:
      paperflow pubmed-fetch --file ./pmid_list.txt --output-dir ./MyPapers --email "YOUR_EMAIL@example.com" --api-key "YOUR_NCBI_API_KEY" 
    """
    fetcher = PubmedFetcher(root_dir=storage_dir, entrez_email=email, api_key=api_key or "", batch_size=batch_size, max_retries=max_retries)
    
    papers = []
    
    # for meta data or full paper data from different sources, we will save them to corresponding paper source database folders
    output_dir = f"{output_dir}/pubmed" if output_dir else f"{storage_dir}/pubmed"
    os.makedirs(output_dir, exist_ok=True)
    
    if query:
        typer.echo(f"Fetching papers for query: {query}")
        query_meta = fetcher.query_search(query)
        papers = fetcher.fetch_from_query(query_meta, output_dir=output_dir)
        
    elif file:
        if not os.path.exists(file):
            typer.echo(f"Error: File {file} not found.")
            raise typer.Exit(code=1)
        with open(file, 'r') as f:
            pmid_list = [line.strip() for line in f if line.strip()]
        typer.echo(f"Fetching {len(pmid_list)} papers from file {os.path.abspath(file)}.")
        papers = fetcher.fetch_from_pmid_list(pmid_list, output_dir=output_dir)
        
    else:
        typer.echo("Error: Must provide --query or --file.")
        raise typer.Exit(code=1)


@app.command("pubmed-content")
def download_fulltext_cmd(
    file: Optional[str] = typer.Option(None, "--file", "-f", help="File containing PMIDs (one per line)."),
    email: str = opt_email,
    api_key: Optional[str] = opt_api_key,
    storage_dir: str = opt_storage,
    max_retries: int = opt_max_retries,
    output_dir: Optional[str] = typer.Option(".", "--output-dir", "-o", help="Directory in result-level to store output full texts, default is current directory. If not specified, will be set to root directory of the repository-level which is storage_dir. 🌟 We will create a '/pubmed' subfolder under the output directory to save all pubmed related data"),
    pmid: Optional[List[str]] = typer.Option(None, "--pmid", "-p", help="Single PMID to download full text for, can be repeated."),
):
    """
    Download full text (PMC) for given PMIDs if the paper has a PMC ID.


    \b
    Notes: 
    - 1, This currently only supports PMC full text fetching if the paper has a PMC ID.


    \b
    Example usage:
    - 1. Download full text for PMIDs listed in a file:
      paperflow download-fulltext --file ./pmid_list.txt --email "YOUR_EMAIL@example" --api-key "YOUR_NCBI_API_KEY" --output-dir ./MyPapers
      
    """
    fetcher = PubmedFetcher(root_dir=storage_dir, entrez_email=email, api_key=api_key or "", max_retries=max_retries)
    
    # create pubmed subfolder in output directory to save full text data
    output_dir = f"{output_dir}/pubmed" if output_dir else f"{storage_dir}/pubmed"
    os.makedirs(output_dir, exist_ok=True)
    
    target_pmids = []
    if file:
        with open(file, 'r') as f:
            target_pmids = [line.strip() for line in f if line.strip()]
    elif pmid:
        target_pmids = pmid
    else:
        typer.echo("Error: Must provide --file or --pmid.")
        raise typer.Exit(code=1)
        
    typer.echo(f"Downloading full texts for {len(target_pmids)} PMIDs from file {os.path.abspath(file) if file else 'provided PMIDs'}.")
    fetcher.fetch_pmc_full_text(target_pmids, output_dir=output_dir)

@app.command("pubmed-all")
def fetch_full_cmd(
    query: Optional[str] = typer.Option(None, "--query", "-q", help="PubMed search query."),
    file: Optional[str] = typer.Option(None, "--file", "-f", help="Text file containing PMIDs (one per line), -q and -f are mutually exclusive."),
    pmid: Optional[List[str]] = typer.Option(None, "--pmid", "-p", help="Single PMID to download full text for, can be repeated."),
    batch_size: int = opt_batch_size,
    max_retries: int = opt_max_retries,
    email: str = opt_email,
    api_key: Optional[str] = opt_api_key,
    storage_dir: str = opt_storage,
    output_dir: Optional[str] = typer.Option(None, "--output-dir", "-o", help="Directory in result-level to store output papers. If not specified, defaults to storage-dir."),
):
    """
    Fetch BOTH metadata and full text (if available) for papers.
    Also extracts URLs from full text and updates metadata links.

    \b
    Example usage:
    - 1. Fetch full papers for a query:
      paperflow pubmed-all --query "machine learning" --output-dir ./MyPapers --email "YOUR_EMAIL"
    """
    fetcher = PubmedFetcher(root_dir=storage_dir, entrez_email=email, api_key=api_key or "", batch_size=batch_size, max_retries=max_retries)
    
    # create pubmed subfolder in output directory to save all pubmed related data (metadata + full text)
    output_dir = f"{output_dir}/pubmed" if output_dir else f"{storage_dir}/pubmed"
    os.makedirs(output_dir, exist_ok=True)
    
    pmid_list = []
    if file:
        if os.path.exists(file):
             with open(file, 'r') as f:
                pmid_list = [line.strip() for line in f if line.strip()]
    if pmid:
        pmid_list.extend(pmid)

    if not query and not pmid_list:
         typer.echo("Error: Must provide --query, --file, or --pmid.")
         raise typer.Exit(code=1)

    fetcher.fetch_and_save_full_papers(query=query, pmid_list=pmid_list if pmid_list else None, output_dir=output_dir)


@app.command("pubmed-merge-json")
def merge_json_cmd(
    paper_dir: str = typer.Option(..., "--input", "-i" , help="Directory containing paper data ({INPUT_PAPER_DIR_HERE}/pubmed/year/pmid/structure)."),
    output: str = typer.Option(..., "--output", "-o" , help="Output directory or file path. If a directory or path without extension is given, the merged file is auto-named as <input-directory-base-name>_<datetime>.json/.jsonl."),
    pmid_file: Optional[str] = typer.Option(None, "--pmid-file", "-p", help="File containing PMIDs to merge (one per line)."),
    jsonl: bool = typer.Option(False, "--jsonl", help="Write output as JSONL, one JSON per line."),
    stats_path: Optional[str] = typer.Option(".", "--stats-path", "-s", help="Optional path to save merge statistics file, defaults to current directory.")
):
    """
    Create a merged JSON (or JSONL) file from PubMed paper directories.

    This produces a canonical merged JSON representation per paper and is
    intended as the first stage in a two-stage pipeline (merge-json -> export-md).

    \b
    Example usage:
    - 1. Merge JSON files for all papers in a directory:
      paperflow pubmed-merge-json --input ./MyPapers --output ./MyPapers 
    - 2. Merge JSON files for PMIDs listed in a file:
      paperflow pubmed-merge-json --input ./MyPapers --output ./MyPapers --pmid-file ./pmid_list.txt --jsonl --stats-path ./MyPapers/stats
    """
    merger = PubmedMerger()

    try:
        stats = merger.merge_json_from_directory(paper_dir, output, pmid_file=pmid_file, jsonl=jsonl)
        # Save stats if path provided
        with open(f"{stats_path}/{os.path.basename(os.path.normpath(paper_dir))}_stats_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.json", "w") as f:
            json.dump(stats, f, indent=2)
        
        typer.echo(f"✅ Please check the merged pubmed JSON/JSONL file at {output} and the merge statistics file at {stats_path}. \
            Also, a JSON file per paper is created within the PMID subfolders.")
    except Exception as e:
        typer.echo(f"Error during merge-json: {e}")
        raise typer.Exit(code=1)


@app.command("pubmed-export-md")
def export_md_cmd(
    merged_json: str = typer.Option(...,"--input", "-i", help="Path to merged JSON or JSONL produced by pubmed-merge-json."),
    output_md: str = typer.Option(...,"--output", "-o", help="Output Markdown file path."),
    yaml_cfg: Optional[str] = typer.Option(None, "--config", "-c", help="YAML config file specifying metadata_fields and content_sections. If not provided, defaults to basic metadata and FULL content."),
    pmid_file: Optional[str] = typer.Option(None, "--pmid-file", "-p", help="Optional PMID file to filter exported papers."),
):
    """
    Export a single Markdown view from a merged JSON file using optional YAML config.

    \b
    Notes:
    - 1, The input merged JSON/JSONL should be produced by the pubmed-merge-json command, which creates a canonical representation of paper metadata and content.
    - 2, The optional YAML config can specify which metadata fields and content sections to include in the Markdown output. If not provided, it defaults to including basic metadata and the FULL content.

    \b
    Example usage:
    - 1. Export Markdown for all papers in a merged JSON:
    paperflow pubmed-export-md --input ./MyPapers/merged.jsonl --output ./MyPapers/exported.md --config ./config.yaml
    - 2. Export Markdown for PMIDs listed in a file:
    paperflow pubmed-export-md --input ./MyPapers/merged.jsonl --output ./MyPapers/exported.md --config ./config.yaml --pmid-file ./pmid_list.txt
    
    
    """
    merger = PubmedMerger()

    try:
        stats = merger.export_md_from_merged_json(merged_json, output_md, yaml_cfg=yaml_cfg, pmid_file=pmid_file)
        typer.secho(f"Successfully exported {stats.get('total', 0)} papers to {stats.get('output')}", fg=typer.colors.GREEN)
    except Exception as e:
        typer.echo(f"Error during export-md: {e}")
        raise typer.Exit(code=1)



#############################################################
#  2, For BioRxiv Parser
#############################################################


@app.command("arxiv-search")
def arxiv_search_cmd(
    query: str = typer.Argument(..., help="arXiv search query."),
    max_results: int = typer.Option(100, "--max-results", "-n", help="Maximum number of arXiv results to return."),
    storage_dir: str = opt_storage,
    output_dir: Optional[str] = typer.Option(None, "--output-dir", "-o", help="Directory to save searched arXiv IDs."),
    start_date: Optional[str] = typer.Option(None, "--start-date", help="Optional start date in YYYY-MM-DD."),
    end_date: Optional[str] = typer.Option(None, "--end-date", help="Optional end date in YYYY-MM-DD."),
    backend: str = opt_arxiv_backend,
):
    """Search arXiv and write matching IDs to a text file."""
    fetcher = ArxivFetcher(root_dir=storage_dir, backend=backend)
    records = fetcher.search(query=query, max_results=max_results, start_date=start_date, end_date=end_date)
    typer.echo(f"Found {len(records)} arXiv papers.")
    for record in records:
        typer.echo(record.source_id)

    save_dir = output_dir if output_dir else storage_dir
    output_file = _save_id_list(save_dir, "searched_arxiv_ids.txt", [record.source_id for record in records])
    typer.echo(f"arXiv IDs saved to {output_file}.")


@app.command("arxiv-fetch")
def arxiv_fetch_cmd(
    query: str = typer.Argument(..., help="arXiv search query."),
    max_results: int = typer.Option(100, "--max-results", "-n", help="Maximum number of arXiv records to fetch."),
    storage_dir: str = opt_storage,
    output_dir: Optional[str] = typer.Option(None, "--output-dir", "-o", help="Directory to save fetched arXiv papers."),
    start_date: Optional[str] = typer.Option(None, "--start-date", help="Optional start date in YYYY-MM-DD."),
    end_date: Optional[str] = typer.Option(None, "--end-date", help="Optional end date in YYYY-MM-DD."),
    download_pdf: bool = typer.Option(True, "--download-pdf/--no-download-pdf", help="Download PDFs when available."),
    backend: str = opt_arxiv_backend,
):
    """Fetch arXiv metadata and attempt to download PDFs."""
    fetcher = ArxivFetcher(root_dir=storage_dir, backend=backend)
    records = fetcher.fetch_from_query(
        query=query,
        output_dir=output_dir if output_dir else storage_dir,
        max_results=max_results,
        start_date=start_date,
        end_date=end_date,
        download_pdf=download_pdf,
    )
    typer.echo(f"Fetched {len(records)} arXiv papers.")


@app.command("biorxiv-search")
def biorxiv_search_cmd(
    query: str = typer.Argument(..., help="bioRxiv search query."),
    max_results: int = typer.Option(100, "--max-results", "-n", help="Maximum number of bioRxiv results to return."),
    storage_dir: str = opt_storage,
    output_dir: Optional[str] = typer.Option(None, "--output-dir", "-o", help="Directory to save searched bioRxiv IDs."),
    start_date: Optional[str] = typer.Option(None, "--start-date", help="Optional start date in YYYY-MM-DD."),
    end_date: Optional[str] = typer.Option(None, "--end-date", help="Optional end date in YYYY-MM-DD."),
    window_days: int = typer.Option(365, "--window-days", help="Compatibility-only option. Retained for older scripts; not used by current Crossref-backed direct query path."),
):
    """Search bioRxiv and write matching IDs to a text file.

    The current implementation uses Crossref server-side query over openRxiv records
    instead of date-window paging over the legacy bioRxiv details API.
    """
    if window_days != 365:
        typer.secho(
            "Note: --window-days is a compatibility-only option and is ignored by the current Crossref-backed direct query path.",
            fg=typer.colors.YELLOW,
        )

    fetcher = BioRxivFetcher(root_dir=storage_dir, window_days=window_days)
    records = fetcher.search(query=query, start_date=start_date, end_date=end_date, max_results=max_results)
    typer.echo(f"Found {len(records)} bioRxiv papers.")
    for record in records:
        typer.echo(record.source_id)

    save_dir = output_dir if output_dir else storage_dir
    output_file = _save_id_list(save_dir, "searched_biorxiv_ids.txt", [record.source_id for record in records])
    typer.echo(f"bioRxiv IDs saved to {output_file}.")


@app.command("biorxiv-fetch")
def biorxiv_fetch_cmd(
    query: str = typer.Argument(..., help="bioRxiv search query."),
    max_results: int = typer.Option(100, "--max-results", "-n", help="Maximum number of bioRxiv records to fetch."),
    storage_dir: str = opt_storage,
    output_dir: Optional[str] = typer.Option(None, "--output-dir", "-o", help="Directory to save fetched bioRxiv papers."),
    start_date: Optional[str] = typer.Option(None, "--start-date", help="Optional start date in YYYY-MM-DD."),
    end_date: Optional[str] = typer.Option(None, "--end-date", help="Optional end date in YYYY-MM-DD."),
    window_days: int = typer.Option(365, "--window-days", help="Compatibility-only option. Retained for older scripts; not used by current Crossref-backed direct query path."),
    download_pdf: bool = typer.Option(True, "--download-pdf/--no-download-pdf", help="Download PDFs when available."),
):
    """Fetch bioRxiv metadata and attempt to download PDFs.

    Metadata retrieval uses Crossref server-side query over openRxiv records.
    """
    if window_days != 365:
        typer.secho(
            "Note: --window-days is a compatibility-only option and is ignored by the current Crossref-backed direct query path.",
            fg=typer.colors.YELLOW,
        )

    fetcher = BioRxivFetcher(root_dir=storage_dir, window_days=window_days)
    records = fetcher.fetch_from_query(
        query=query,
        output_dir=output_dir if output_dir else storage_dir,
        start_date=start_date,
        end_date=end_date,
        max_results=max_results,
        download_pdf=download_pdf,
    )
    typer.echo(f"Fetched {len(records)} bioRxiv papers.")



#############################################################
#  3, For Third-Party integrations 
#############################################################

# 3.1 paper fetch related commands

@app.command(
    "paper-fetch",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
    add_help_option=False,  # Disable typer's --help so argparse handles it
)
def paper_fetch_cmd(ctx: typer.Context):
    """
    Fetch PDFs by DOI — passes through to the paper-fetch engine.

    Use ``paperflow paper-fetch --help`` to see the full parameter list.

    \b
    Notes:
    - 1, This command is a thin wrapper around the paper-fetch engine, which is a powerful tool for fetching papers by DOI, title, or batch lists. 
    - 2, ⚠️ Remember to set Unpaywall email in environment variables for best performance when fetching by DOI.

    \b
    Example usage:
      paperflow paper-fetch 10.1038/s41586-020-2649-2 -o ./papers
      paperflow paper-fetch --batch dois.txt -o ./papers --format text
      paperflow paper-fetch --title "AlphaFold" -o ./papers
    """

    try:
        pdf_fetch.run(["paper-fetch"] + ctx.args)
    except SystemExit as e:
        if e.code != 0:
            raise typer.Exit(code=e.code)


# 3.2 mineru related commands
# including pdf parse, post-parse structuring, and final markdown conversion

@app.command("pdf-parse")
def pdf_parse_cmd( 
    input_path: str = typer.Option(..., "--input", "-i", help="Input PDF file path."),
    output_dir: str = typer.Option(..., "--output", "-o", help="Output directory for parsed output."),
    clear: bool = typer.Option(False, "--clear", help="After conversion, keep only the .md files and necessary .json files(_content_list_v2.json/_content_list.json)."),
):
    """
    Parse a PDF file using MinerU engine, and clean up the output directory.

    \b
    Notes:
    - 1, MinerU generates a subfolder /auto under --output with .md, .json, .pdf, and images/.  Use --clear to strip anything unnecessary, 
    note that we only use .md files and _content_list_v2.json/_content_list.json files for further processing like structuring.
    - 2, ⚠️  Remember to switch to domestic mirror source when you can not access huggingface.

    \b
    Example usage:
      paperflow pdf-parse -i paper.pdf -o ./output
    """
    input_p = Path(input_path)
    if not input_p.exists():
        typer.echo(f"Error: Input path not found: {input_path}")
        raise typer.Exit(code=1)

    output_p = Path(output_dir)
    output_p.mkdir(parents=True, exist_ok=True)

    # Snapshot existing top-level entries before mineru creates its subfolder
    _before: set[Path] = {p for p in output_p.iterdir()} if output_p.exists() else set()

    cmd = ["mineru", "-p", str(input_p.resolve()), "-o", str(output_p.resolve()), "-b", "pipeline"]
    typer.echo(f"Running: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        typer.secho("Done.", fg=typer.colors.GREEN)
    except subprocess.CalledProcessError as e:
        typer.secho(f"MinerU failed with exit code {e.returncode}", fg=typer.colors.RED)
        raise typer.Exit(code=e.returncode)
    except FileNotFoundError:
        typer.secho("Error: mineru not found. Please install MinerU first.", fg=typer.colors.RED)
        raise typer.Exit(code=1)


    def _clean_mineru_dirs(dirs: list[Path]) -> None:
        """
        Clean up MinerU output: only keep .md files and necessary .json files(_content_list_v2.json/_content_list.json), 
        in case memory deficient when batch processing many PDFs with MinerU. 
        """
        import shutil
        removed = 0
        # subfolders under output_p are the ones created by mineru
        for top_dir in dirs:
            # generally, only one subfolder is created by mineru(only one top_dir)
            for pdf_file in top_dir.rglob("*.pdf"):
                pdf_file.unlink()  # delete all the pdf files 
                removed += 1
            for target_json in top_dir.rglob("*.json"):
                if target_json.name.endswith(("_content_list_v2.json", "_content_list.json")):
                    continue  # keep necessary json files for further processing
                target_json.unlink()  # delete all the unwanted json files 
                removed += 1
        if removed:
            typer.echo(f"✅ Removed {removed} source files. Only .md and necessary .json files are kept in the output directory {output_p}.")


    if clear:
        _new = [p for p in output_p.iterdir() if p not in _before and p.is_dir()]
        _clean_mineru_dirs(_new)



@app.command("mineru-parse")
def mineru_parse_cmd(
    input_json: str = typer.Option(..., "--input", "-i",
        help="Path to mineru content_list_v2.json."),
    output_json: str = typer.Option(..., "--output", "-o",
        help="Output path for the structured JSON file."),
    backend: str = typer.Option("regex", "--backend", "-b",
        help="Section classification backend: 'regex' (default, no API needed) or 'ai'."),
    config: Optional[str] = typer.Option(None, "--config", "-c",
        help="Path to YAML config file for canonical types, aliases, and AI settings."),
    api_key: Optional[str] = typer.Option(None, "--api-key",
        help="API key for AI backend. Overrides config file and env var."),
    model: Optional[str] = typer.Option(None, "--model",
        help="Override AI model (e.g. 'deepseek-v4-pro', 'claude-haiku-4-5', 'gpt-4o-mini')."),
    base_url: Optional[str] = typer.Option(None, "--base-url",
        help="Custom API base URL for OpenAI-compatible endpoints (e.g. 'https://api.deepseek.com')."),
):
    """
    Parse mineru output content_list_v2.json into canonical sectioned JSON.

    Extracts metadata (title, authors, year, DOI, journal),
    and sections normalised to canonical types (abstract, introduction, results,
    discussion, methods, etc.). Tables are preserved as HTML.

    \b
    Notes:
    - 1, Two backends: 'regex' (pattern + context, no API) and 'ai' (LLM batch classification).
    - 2, AI backend supports Anthropic native, OpenAI native, and any OpenAI-compatible
    endpoint via --base-url (DeepSeek, university proxies, self-hosted, etc.).
    - 3, Set the appropriate API key env var (ANTHROPIC_API_KEY, OPENAI_API_KEY,
    DEEPSEEK_API_KEY) or pass --api-key.
    - 4, Configure provider/model via --model, --base-url, or a YAML config file.

    \b
    Examples:
      paperflow mineru-parse -i content_list_v2.json -o paper.json
      paperflow mineru-parse -i content_list_v2.json -o paper.json --backend ai
      paperflow mineru-parse -i content_list_v2.json -o paper.json --backend ai \\
          --base-url https://api.deepseek.com --model deepseek-v4-pro --api-key sk-xxx
      paperflow mineru-parse -i content_list_v2.json -o paper.json --backend ai \\
          --base-url https://models.sjtu.edu.cn/api/v1 --model deepseek-chat
      paperflow mineru-parse -i content_list_v2.json -o paper.json --backend regex --config custom.yaml
    """
    import json as _json
    from .integrations.mineru_parser import (
        MinerUContentParser,
        RegexSectionClassifier,
        AISectionClassifier,
    )

    if backend == "ai":
        classifier = AISectionClassifier.from_config(config)
        if api_key:
            classifier.api_key = api_key
        if model:
            classifier.model = model
        if base_url:
            classifier.base_url = base_url
        if base_url:
            typer.echo(f"Using AI backend: {classifier.model} @ {classifier.base_url}")
        else:
            typer.echo(f"Using AI backend: {classifier.model}")
    else:
        classifier = RegexSectionClassifier.from_config(config)
        typer.echo("Using regex backend with configurable aliases")

    parser = MinerUContentParser(classifier)
    result = parser.parse(input_json)
    Path(output_json).parent.mkdir(parents=True, exist_ok=True)
    with open(output_json, "w") as f:
        _json.dump(result, f, ensure_ascii=False, indent=2)

    section_summary = ", ".join(
        f"{s['canonical_type']}({s.get('display_title', '?')})"
        for s in result["sections"]
    )
    typer.echo(
        f"Parsed {len(result['sections'])} sections -> {output_json}"
    )
    typer.echo(f"  Sections: {section_summary}")


@app.command("mineru-export-md")
def mineru_export_md_cmd(
    input_json: str = typer.Option(..., "--input", "-i",
        help="Path to structured JSON file (from mineru-parse), or a directory of such files."),
    output_md: str = typer.Option(..., "--output", "-o",
        help="Output Markdown file path."),
    yaml_cfg: Optional[str] = typer.Option(None, "--config", "-c",
        help="YAML config specifying content_sections to include. If not provided, all sections are included."),
):
    """
    Export structured mineru JSON to a clean Markdown file for LLM processing.

    Reads one or more JSON files produced by ``mineru-parse`` and writes a
    single Markdown file.  Metadata (title, authors, year, DOI, journal) is
    always included.  Content sections are included based on the optional
    YAML config.

    \b
    YAML config format:
      content_sections:
        - abstract
        - introduction
        - methods
        - results
        - discussion
        - conclusion

    \b
    Examples:
      paperflow mineru-export-md -i paper.json -o paper.md
      paperflow mineru-export-md -i paper.json -o paper.md --config extract.yaml
      paperflow mineru-export-md -i ./parsed_dir -o all_papers.md
    """
    from .integrations.mineru_parser import export_mineru_json_to_md

    try:
        stats = export_mineru_json_to_md(input_json, output_md, yaml_cfg)
        typer.secho(
            f"Exported {stats['total']} papers to {stats['output']}",
            fg=typer.colors.GREEN,
        )
        if stats.get("sections_exported"):
            sec_str = ", ".join(
                f"{k}({v})" for k, v in stats["sections_exported"].items()
            )
            typer.echo(f"  Sections exported: {sec_str}")
    except Exception as e:
        typer.echo(f"Error: {e}")
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
