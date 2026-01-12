from typing import List, Dict, Literal
import re


def extract_urls_from_text(text: str, source_tag: Literal["abstract", "full_text"]) -> List[Dict[str, str]]:
    """
    Description
    -----------
    Extract URLs from the given text, and attempt to categorize them.

    Args
    ----
        text (str): The input text from which to extract URLs.
        source (Literal["abstract", "full_text"]): The source of the text, either "abstract" or "full_text".

    Returns
    -------
        List[Dict[str, str]]: A list of dictionaries, each containing:
            - "url": The extracted URL.
            - "source": The source of the URL extraction (e.g., "abstract", "full_text").
            - "category": A simple category based on URL patterns (e.g., "GitHub", "Zenodo", etc.). 
    
    """
    if not text:
        return []

    # Match URLs starting with http://, https://, ftp://, or www.
    url_pattern = r'(https?://[^\s,;>)]+|www\.[^\s,;>)]+|ftp://[^\s,;>)]+)'
    
    found_urls = re.findall(url_pattern, text)
    
    results = []
    seen = set() # deduplicate URLs

    for url in found_urls:
        # clearing trailing punctuation like . , ; ) >
        url = url.rstrip('.,;)>')
        
        if url in seen:
            continue
        seen.add(url)

        # Simple categorization based on URL patterns
        category = "General"
        if "github.com" in url:
            category = "GitHub"
        elif "gitlab.com" in url:
            category = "GitLab"
        elif "zenodo.org" in url:
            category = "Zenodo"
        elif "figshare.com" in url:
            category = "Figshare"
        elif "huggingface.co" in url:
            category = "HuggingFace"
        elif "drive.google.com" in url:
            category = "Google Drive"

        results.append({
            "url": url,
            "source": source_tag, # for now we just label as abstract 
            "category": category
        })
    
    return results


