import typer
import os
import json
from typing import List, Optional
from .fetcher import PubmedFetcher
from .storage import PaperStorage
from .query_builder import QueryBuilder, AIQueryAssistant

app = typer.Typer(help="pyPaperFlow CLI", no_args_is_help=True)

# Common options
opt_storage = typer.Option("./Papers", "--storage-dir", "-s", help="Directory in Repository-level to store paper data for Initialization.") # Note: this is a repository-level default path
opt_email = typer.Option(..., "--email", help="Entrez Email.")
opt_api_key = typer.Option(None, "--api-key", help="NCBI API Key (recommended).")
opt_max_retries = typer.Option(3, "--max-retries", help="Maximum number of retries for Entrez API calls.")

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
    batch_size: int = typer.Option(50, "--batch-size", "-b", help="Batch size for fetching."),
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
    fetcher = PubmedFetcher(root_dir=storage_dir, entrez_email=email, api_key=api_key or "", max_retries=3)
    
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

if __name__ == "__main__":
    app()
