import typer
import os
import json
from typing import *
from .fetcher import PubmedFetcher
from .arxiv_fetcher import ArxivFetcher
from .biorxiv_fetcher import BioRxivFetcher
from .merger import (
    PaperMerger,
    MergeConfig,
    MergeMode,
    OutputFormat,
    OutputProfile,
    TextSourcePriority,
)

app = typer.Typer(help="pyPaperFlow CLI", no_args_is_help=True)

# Common options
opt_storage = typer.Option("./Papers", "--storage-dir", "-s", help="Directory in Repository-level to store paper data for Initialization.") # Note: this is a repository-level default path
opt_email = typer.Option(..., "--email", help="Entrez Email.")
opt_api_key = typer.Option(None, "--api-key", help="NCBI API Key (recommended).")
opt_max_retries = typer.Option(3, "--max-retries", help="Maximum number of retries for Entrez API calls.")
opt_batch_size = typer.Option(50, "--batch-size", "-b", help="Batch size for fetching.")


def _save_id_list(output_dir: str, filename: str, values: List[str]) -> str:
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_file = os.path.join(output_dir, filename)
    with open(output_file, "w", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")
    return output_file

@app.command("search")
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
    - 2, This command will print the found PMIDs and also save them to 'searched_pmids.txt' in the specified output directory. 
    If --output-dir is not specified, it will default to the storage directory.
    - 3, Note that storage_dir is used to initialize the fetcher for consistency, while output_dir is where the PMIDs are saved. They are different parameters!

    \b
    Example usage:
    - 1. Search for papers related to "machine learning" and return up to 500 PMIDs/per batch:
    paperflow search "machine learning" --retmax 500 --output-dir ./MyPapers --email "YOUR_EMAIL@example.com" --api-key "YOUR_NCBI_API_KEY"
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
    pmid_file = os.path.join(save_dir, "searched_pmids.txt")
    with open(pmid_file, 'w') as f:
        for pmid in pmids:
            f.write(f"{pmid}\n")
    typer.echo(f"PMIDs saved to {pmid_file}.")

@app.command("fetch")
def fetch_cmd(
    query: Optional[str] = typer.Option(None, "--query", "-q", help="PubMed search query."),
    file: Optional[str] = typer.Option(None, "--file", "-f", help="Text file containing PMIDs (one per line), -q and -f are mutually exclusive."),
    batch_size: int = opt_batch_size,
    email: str = opt_email,
    api_key: Optional[str] = opt_api_key,
    storage_dir: str = opt_storage,
    max_retries: int = opt_max_retries,
    output_dir: Optional[str] = typer.Option(".", "--output-dir", "-o", help="Directory in result-level to store output papers, default is current directory. If not specified, will be set to root directory of the repository-level which is storage_dir."),
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
      paperflow fetch --query "machine learning" --output-dir ./MyPapers --email "YOUR_EMAIL@example.com" --api-key "YOUR_NCBI_API_KEY"
    - 2. Fetch papers from a list of PMIDs in a file:
      paperflow fetch --file ./pmid_list.txt --output-dir ./MyPapers --email "YOUR_EMAIL@example.com" --api-key "YOUR_NCBI_API_KEY" 
    """
    fetcher = PubmedFetcher(root_dir=storage_dir, entrez_email=email, api_key=api_key or "", batch_size=batch_size, max_retries=max_retries)
    storage = PaperStorage(storage_dir)
    
    papers = []
    
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
        
    # Save to storage
    # typer.echo(f"Saving {len(papers)} papers to storage...")
    # for paper in papers:
    #    storage.add_paper(paper)
    # typer.echo("Done.")

@app.command("download-fulltext")
def download_fulltext_cmd(
    file: Optional[str] = typer.Option(None, "--file", "-f", help="File containing PMIDs (one per line)."),
    email: str = opt_email,
    api_key: Optional[str] = opt_api_key,
    storage_dir: str = opt_storage,
    max_retries: int = opt_max_retries,
    output_dir: Optional[str] = typer.Option(".", "--output-dir", "-o", help="Directory in result-level to store output full texts, default is current directory. If not specified, will be set to root directory of the repository-level which is storage_dir."),
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
    storage = PaperStorage(storage_dir)
    
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

@app.command("fetch-full")
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
      paperflow fetch-full --query "machine learning" --output-dir ./MyPapers --email "YOUR_EMAIL"
    """
    fetcher = PubmedFetcher(root_dir=storage_dir, entrez_email=email, api_key=api_key or "", batch_size=batch_size, max_retries=max_retries)
    
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

@app.command("tag")
def tag_cmd(
    pmid: str = typer.Argument(..., help="PMID to tag."),
    tags: List[str] = typer.Option([], "--tag", "-t", help="Tag(s) to add. Can be repeated."),
    remove: List[str] = typer.Option([], "--remove", "-r", help="Tag(s) to remove. Can be repeated."),
    clear: bool = typer.Option(False, "--clear", help="Clear all existing tags before applying --tag/--remove."),
    storage_dir: str = opt_storage,
):
    """
    Add/remove tags for a paper.

    Examples:
    - Add tags:
      paperflow tag 12345678 -t "蛋白" -t "AI"
    - Remove tags:
      paperflow tag 12345678 -r "ai"
    - Replace all tags:
      paperflow tag 12345678 --clear -t "结构" -t "protein"
    """
    storage = PaperStorage(storage_dir)

    if clear:
        # Clear everything then add desired tags.
        storage.set_tags(pmid, [])

    if tags:
        storage.add_tags(pmid, tags)

    if remove:
        storage.remove_tags(pmid, remove)

    current = storage.list_tags(pmid)
    typer.echo(f"PMID {pmid} tags ({len(current)}): {current}")


@app.command("tag-set")
def tag_set_cmd(
    pmid: str = typer.Argument(..., help="PMID to tag."),
    tag: str = typer.Argument(..., help="Tag name."),
    value: int = typer.Argument(..., help="Tag value (0 or 1)."),
    storage_dir: str = opt_storage,
):
    """Set a tag explicitly to 0/1 (backward-compatible)."""
    storage = PaperStorage(storage_dir)
    storage.update_tags(pmid, {tag: value})
    typer.echo(f"Set tag '{tag}' to {value} for PMID {pmid}.")

@app.command("query")
def query_cmd(
    tags: List[str] = typer.Option([], "--tag", "-t", help="Tags to filter by (format: name=value)."),
    storage_dir: str = opt_storage
):
    """
    Query papers by tags.
    """
    storage = PaperStorage(storage_dir)
    
    query_dict = {}
    for t in tags:
        try:
            k, v = t.split("=")
            query_dict[k] = int(v)
        except ValueError:
            typer.echo(f"Invalid tag format: {t}. Use name=value.")
            raise typer.Exit(code=1)
            
    pmids = storage.query_papers(query_dict)
    typer.echo(f"Found {len(pmids)} matching papers:")
    for pmid in pmids:
        typer.echo(pmid)

@app.command("get")
def get_cmd(
    pmid: str = typer.Argument(..., help="PMID to retrieve."),
    storage_dir: str = opt_storage
):
    """
    Get paper details.
    """
    storage = PaperStorage(storage_dir)
    paper = storage.get_paper(pmid)
    
    if paper:
        typer.echo(json.dumps(paper.to_dict(), indent=2, ensure_ascii=False))
        
        # Also show tags
        tags = storage.get_feature_vector(pmid)
        typer.echo("\nTags:")
        typer.echo(json.dumps(tags, indent=2))
    else:
        typer.echo(f"Paper {pmid} not found.")

@app.command("build-query")
def build_query_cmd(
    openai_api_key: Optional[str] = typer.Option(None, "--openai-key", help="OpenAI API Key for AI assistance."),
    model: str = typer.Option("gpt-4o", "--model", help="AI Model to use (e.g., gpt-4o, gpt-3.5-turbo)."),
    api_base: str = typer.Option("https://api.openai.com/v1", "--api-base", help="Base URL for AI API.")
):
    """
    Interactive wizard to build complex PubMed queries with AI assistance.
    """
    ai_assistant = None
    if openai_api_key:
        ai_assistant = AIQueryAssistant(api_key=openai_api_key, api_base=api_base, model=model)
    elif os.environ.get("OPENAI_API_KEY"):
        ai_assistant = AIQueryAssistant(api_key=os.environ.get("OPENAI_API_KEY"), api_base=api_base, model=model)
    else:
        typer.secho("Warning: No OpenAI API Key provided. AI features will be disabled.", fg=typer.colors.YELLOW)
        typer.echo("You can provide it via --openai-key or OPENAI_API_KEY environment variable.")

    builder = QueryBuilder(ai_assistant)
    builder.run()


@app.command("arxiv-search")
def arxiv_search_cmd(
    query: str = typer.Argument(..., help="arXiv search query."),
    max_results: int = typer.Option(100, "--max-results", "-n", help="Maximum number of arXiv results to return."),
    storage_dir: str = opt_storage,
    output_dir: Optional[str] = typer.Option(None, "--output-dir", "-o", help="Directory to save searched arXiv IDs."),
    start_date: Optional[str] = typer.Option(None, "--start-date", help="Optional start date in YYYY-MM-DD."),
    end_date: Optional[str] = typer.Option(None, "--end-date", help="Optional end date in YYYY-MM-DD."),
):
    """Search arXiv and write matching IDs to a text file."""
    fetcher = ArxivFetcher(root_dir=storage_dir)
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
):
    """Fetch arXiv metadata and attempt to download PDFs."""
    fetcher = ArxivFetcher(root_dir=storage_dir)
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
    window_days: int = typer.Option(365, "--window-days", help="Date window size for bioRxiv API paging."),
):
    """Search bioRxiv and write matching IDs to a text file."""
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
    window_days: int = typer.Option(365, "--window-days", help="Date window size for bioRxiv API paging."),
    download_pdf: bool = typer.Option(True, "--download-pdf/--no-download-pdf", help="Download PDFs when available."),
):
    """Fetch bioRxiv metadata and attempt to download PDFs."""
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

@app.command("merge")
def merge_cmd(
    paper_dir: str = typer.Argument(..., help="Directory containing paper data (year/pmid/ structure)."),
    output: str = typer.Argument(..., help="Output file path for merged data."),
    pmid_file: Optional[str] = typer.Option(None, "--pmid-file", "-p", help="File containing PMIDs to merge (one per line). If not specified, merge all papers in directory."),
    mode: str = typer.Option("full", "--mode", "-m", help="Merge mode: 'meta' (metadata only), 'full' (metadata + content)."),
    format: str = typer.Option("md", "--format", "-f", help="Output format: 'md' (Markdown), 'jsonl' (JSON Lines), 'txt' (plain text)."),
    profile: str = typer.Option("analysis", "--profile", help="Output profile: 'analysis' (full metadata for modules), 'llm' (LLM-focused compact output)."),
    include_sections: Optional[str] = typer.Option(None, "--include-sections", help="Comma-separated section names to include from parsed JSON, e.g. 'abstract,introduction,results'."),
    metadata_fields: Optional[str] = typer.Option(None, "--metadata-fields", help="Comma-separated metadata field paths, e.g. 'identity,content.keywords,links'. Use 'all' for full metadata."),
    include_links: Optional[bool] = typer.Option(None, "--include-links/--exclude-links", help="Whether to include links in llm profile output. Ignored in analysis profile."),
    text_source: str = typer.Option("parsed-json-first", "--text-source", help="Text source priority: 'parsed-json-first' or 'parsed-md-first'."),
):
    """
    Merge PubMed paper data into unified format for AI analysis.

    This command merges paper metadata and/or full text content into a single
    file optimized for AI processing and analysis.

    Examples:
      - Merge all papers with content to Markdown:
        paperflow merge ./papers_dir ./merged_papers.md --mode full --format md

      - Merge papers from a PMID list file:
        paperflow merge ./papers_dir ./selected_papers.md --pmid-file pmids.txt --mode full
    """
    # Map mode string to enum
    mode_map = {
        'meta': MergeMode.METADATA_ONLY,
                'full': MergeMode.METADATA_CONTENT,
    }

    if mode not in mode_map:
        typer.echo(f"Error: Invalid mode '{mode}'. Use: meta, full")
        raise typer.Exit(code=1)

    merge_mode = mode_map[mode]

    # Map format string to enum
    format_map = {
        'md': OutputFormat.MARKDOWN,
        'jsonl': OutputFormat.JSONL,
        'txt': OutputFormat.PLAIN_TEXT
    }

    if format not in format_map:
        typer.echo(f"Error: Invalid format '{format}'. Use: md, jsonl, txt")
        raise typer.Exit(code=1)

    output_format = format_map[format]

    # Map profile string to enum
    profile_map = {
        'analysis': OutputProfile.ANALYSIS,
        'llm': OutputProfile.LLM,
    }

    if profile not in profile_map:
        typer.echo(f"Error: Invalid profile '{profile}'. Use: analysis, llm")
        raise typer.Exit(code=1)

    output_profile = profile_map[profile]

    text_source_map = {
        'parsed-json-first': TextSourcePriority.PARSED_JSON_FIRST,
        'parsed-md-first': TextSourcePriority.PARSED_MD_FIRST,
    }

    if text_source not in text_source_map:
        typer.echo("Error: Invalid text-source. Use: parsed-json-first, parsed-md-first")
        raise typer.Exit(code=1)

    include_sections_list = [item.strip() for item in include_sections.split(',')] if include_sections else None
    metadata_fields_list = [item.strip() for item in metadata_fields.split(',')] if metadata_fields else None

    # Create merge config
    config = MergeConfig(
        mode=merge_mode,
        output_format=output_format,
        output_profile=output_profile,
        include_sections=include_sections_list,
        include_metadata_fields=metadata_fields_list,
        include_links_in_llm=(include_links if include_links is not None else False),
        text_source_priority=text_source_map[text_source],
    )

    # Create merger and execute
    merger = PaperMerger(config)

    typer.echo(f"Merging papers from: {paper_dir}")

    try:
        if pmid_file:
            # Merge from PMID list file
            stats = merger.merge_from_file_list(paper_dir, pmid_file, output)
        else:
            # Merge all papers in directory
            stats = merger.merge_from_directory(paper_dir, output)

        merger.print_statistics()

        if stats['successful'] > 0:
            typer.secho(f"Successfully merged {stats['successful']} papers.", fg=typer.colors.GREEN)

        if stats['failed'] > 0:
            typer.secho(f"{stats['failed']} papers failed to merge.", fg=typer.colors.YELLOW)

    except Exception as e:
        typer.echo(f"Error during merge: {e}")
        raise typer.Exit(code=1)

if __name__ == "__main__":
    app()
