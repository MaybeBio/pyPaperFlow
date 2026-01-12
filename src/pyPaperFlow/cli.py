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
opt_max_retries = typer.Option(5, "--max-retries", help="Maximum number of retries for Entrez API calls.")

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
    max_retries: int = opt_max_retries
):
    """
    Download full text (PMC) for given PMIDs if the paper has a PMC ID.


    \b
    Notes: 
    - 1, This currently only supports PMC full text fetching if the paper has a PMC ID.


    \b
    Example usage:
    - 1. Download full text for PMIDs listed in a file:
      paperflow download-fulltext --file ./pmid_list.txt --email "YOUR_EMAIL@example" --api-key "YOUR_NCBI_API_KEY"
      
    """
    fetcher = PubmedFetcher(root_dir=storage_dir, entrez_email=email, api_key=api_key or "", max_retries=max_retries)
    storage = PaperStorage(storage_dir)
    
    target_pmids = []
    if file:
        with open(file, 'r') as f:
            target_pmids = [line.strip() for line in f if line.strip()]
    elif pmids:
        target_pmids = pmids
    else:
        typer.echo("Error: Must provide --file or --pmid.")
        raise typer.Exit(code=1)
        
    # We need to find PMC IDs for these PMIDs.
    # The Paper object might have links, or we can try to convert.
    # For now, let's assume the user might provide PMC IDs or we look up in our storage if we have the paper.
    
    pmc_ids = []
    for pmid in target_pmids:
        paper = storage.get_paper(pmid)
        if paper:
            # Check links for PMC
            # In fetcher.py, we didn't explicitly store PMC ID in PaperIdentity, but it might be in links.
            # However, fetch_pmc_full_text expects PMC IDs.
            # We can use Entrez to convert PMID to PMC ID.
            # Or we can check if the user provided PMC IDs directly (unlikely if they say PMIDs).
            pass
            
    # Since converting PMID to PMC ID is not explicitly in fetcher methods (except via links),
    # we might need to add a converter or rely on what we have.
    # For this MVP, let's just try to fetch using the fetcher's method which expects PMC IDs.
    # But wait, the user asked "Given pmid list, get... PMC...".
    # So we need a conversion step.
    
    # Let's add a quick conversion using Entrez.link if not already done.
    # Actually, fetcher.fetch_linked_data_for_single_pmid does fetch 'pubmed_pmc' links!
    # Let's use that.
    
    typer.echo(f"Resolving PMC IDs for {len(target_pmids)} PMIDs...")
    resolved_pmc_ids = []
    
    for pmid in target_pmids:
        # Try to get from storage first
        paper = storage.get_paper(pmid)
        links = None
        if paper and paper.links:
            links = paper.links
        else:
            # Fetch links live
            links_map = fetcher.fetch_linked_data_for_single_pmid(pmid)
            if pmid in links_map:
                links = links_map[pmid]
        
        if links:
            # Check for PMC link
            # In fetcher.py: pubmed_pmc linkname
            # But wait, fetcher.py stores 'pubmed_pmc' in links.entrez dictionary?
            # Let's check fetcher.py:
            # if db_link['linkname'] == "pubmed_pubmed_citedin": ...
            # else: links_map[pmid].entrez[db_link['linkname']] = db_link['links']
            
            # So we look for 'pubmed_pmc' in links.entrez
            pmc_link = links.entrez.get('pubmed_pmc', [])
            if pmc_link:
                # These are usually PMC IDs (numbers or PMC+numbers)
                # Entrez usually returns just the ID part for some links, but for PMC it might be the PMC ID.
                # Let's assume it is.
                for pid in pmc_link:
                    resolved_pmc_ids.append(f"PMC{pid}" if not pid.startswith("PMC") else pid)
            else:
                typer.echo(f"No PMC link found for PMID {pmid}")
        else:
             typer.echo(f"Could not fetch links for PMID {pmid}")

    if resolved_pmc_ids:
        typer.echo(f"Found {len(resolved_pmc_ids)} PMC IDs. Downloading...")
        fetcher.fetch_pmc_full_text(resolved_pmc_ids)
    else:
        typer.echo("No PMC IDs found.")

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
