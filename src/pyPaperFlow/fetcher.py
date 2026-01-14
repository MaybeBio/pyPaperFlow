import json
import os
import time
from typing import List, Dict, Any, Optional, Tuple, Union
from Bio import Entrez, Medline
import traceback
from bs4 import BeautifulSoup
from dataclasses import dataclass, field, asdict
from .utils import extract_urls_from_text


#############################################################
#  1, Some pre-defined Data Classes for Paper Structure
#############################################################

@dataclass(frozen=True)
class PaperIdentity:
    """
    Description
    ----------
        Identity of the paper, storing basic identification information
    
    Args
    ----------
        pmid (str): pubmed ID of the paper
        doi (str): digital object identifier of the paper, default is an empty string
        title (str): title of the paper, default is an empty string
    """
    pmid: str 
    doi: str = ""
    title: str = ""

@dataclass(frozen=True)
class PaperContent:
    """
    Description
    ----------
        Content of the paper, storing main textual information, important for NLP tasks like ontology construction.

    Args
    ----------
        abstract (str): abstract of the paper
        keywords (List[str]): list of keywords associated with the paper
        mesh_terms (List[str]): list of MeSH terms associated with the paper
        pub_types (List[str]): list of publication types of the paper, e.g., Journal Article, Review, etc.
    
    """
    abstract: str
    keywords: List[str] = field(default_factory=list)
    mesh_terms: List[str] = field(default_factory=list)
    pub_types: List[str] = field(default_factory=list)

# we do not use frozen here cause we need to update xml field once the instance is created for medline 
@dataclass
class PaperContributors:
    """
    Description
    ----------
        Contributors of the paper, storing information about authors and affiliations.
        Data is separated by source format (Medline vs XML) to preserve structure.

    Args
    ----------
        medline (Dict[str, List[str]]): Flat lists extracted from Medline format (or flattened from XML).
                                        Keys: 'full_names', 'short_names', 'auids', 'affiliations'.
        xml (List[Dict[str, Any]]): Structured author list with detailed fields (from XML).
    """
    medline: Dict[str, List[str]] = field(default_factory=dict)
    xml: List[Dict[str, Any]] = field(default_factory=list)

@dataclass(frozen=True)
class PaperSource:
    """
    Description
    ----------
        Source information of the paper, storing journal and publication details.
    
    Args
    ----------
        journal_title (str): title of the journal
        journal_abbrev (str): abbreviation of the journal
        pub_date (str): publication date of the paper
        pub_year (str): publication year of the paper
        pub_types (List[str]): list of publication types of the paper, e.g., Journal Article, Review, etc. ⚠️ SAME as Pub_types in PaperContent

    """
    # journal: Dict[str, str] = field(default_factory=dict)
    journal_title: str = ""
    journal_abbrev: str = ""
    pub_date: str = ""
    pub_year: str = ""
    pub_types: List[str] = field(default_factory=list)

@dataclass
class PaperLinks:
    """
    Description
    ----------
        Stores linked data from other NCBI databases (retrieved via ELink).
        This is crucial for the "Deep Ontology" platform.
    
    Args
    ----------
        cites (List[str]): List of PMIDs that cite this paper (pubmed_pubmed_citedin).
        refs (List[str]): List of PMIDs that this paper cites (pubmed_pubmed_refs).
        similar (List[str]): List of related articles (pubmed_pubmed).Similar PubMed articles, obtained by matching text and MeSH terms
        review (List[str]): List of review articles (pubmed_pubmed_reviews).
        pmc (List[str]): List of PMC IDs linked to this paper (pubmed_pmc), we store free full-text links here.
        entrez (Dict[str, List[str]]): Dictionary of other internal Entrez Cross-database links, keys are linkname, values are lists of linked UIDs.
        external (List[Dict[str, str]]): List of external database links
        text_mined (List[Dict[str, str]]): List of text-mined links, store urls mined from abstract or full text
    
    Notes
    ----------
    - 1, Since links may change over time, we store the fetch timestamp for reference.
    That means we need to update link data(e.g. cites) periodically for each paper in our database, and update the fetched article 
    """

    cites: List[str] = field(default_factory=list)
    refs: List[str] = field(default_factory=list)
    similar: List[str] = field(default_factory=list)
    review: List[str] = field(default_factory=list)
    pmc: List[str] = field(default_factory=list)
    entrez: Dict[str, List[str]] = field(default_factory=dict)
    external: List[Dict[str, str]] = field(default_factory=list)
    text_mined: List[Dict[str, str]] = field(default_factory=list)

@dataclass(frozen=True)
class PaperMetadata:
    """
    Description
    ----------
        Metadata of the paper, storing additional information about the fetching process.
    
    Args
    ----------
        entrez_date (str): Entrez date of the paper, used for incremental updates
        fetched_at (str): timestamp when the paper was fetched
    
    Notes
    ----------
    - 1, entrez_date is crucial for incremental updates, as it indicates when the paper was added to PubMed.
    - 2, ⚠️ All information in this class is related to the fetching process, not the paper content itself.
    """
    entrez_date: str = ""
    fetched_at: str = ""

@dataclass
class Paper_MetaData:
    identity: PaperIdentity
    content: PaperContent
    contributors: PaperContributors
    source: PaperSource
    metadata: PaperMetadata
    links: PaperLinks = field(default_factory=PaperLinks)

    def to_dict(self):
        return asdict(self)

@dataclass
class Paper_TextData:
    """
    Description
    ----------
        Text data of the paper, storing the full text content.
    
    Args
    ----------
        pub_year (str): publication year of the paper, parsed from XML <pub-date>
        pmid (str): pubmed ID of the paper
        pmcid (str): pmc ID of the paper
        xml (str): raw XML content of the paper, can be exported to txt format or xml format (xml format is better here)
        parsed_json (Dict[str, Any]): parsed JSON structure of the paper, exported to json format
        parsed_text (str): parsed text content of the paper, can be exported to Markdown or plain text format (Md is better here, we choose Md format)

    Notes
    ----------
    - 1, There are mainly 3 ways to get full text:
        - via PMC Open Access Subset (best quality, XML or PDF converted to text), we can fetch via EFetch API if we have PMC IDs
        - via publisher website (if accessible, may require subscription), we can get urls via ELink such as llinks (test successfully in some cases, but it is hard to fetch automatically due to paywalls and different formats)
        - via Other Mature Github Projects (e.g., SciHub, PaperScraper, etc.), we can integrate their APIs or methods to fetch full text given DOI or URL
    We implement the first method here, and may consider the other two methods in future versions.
    """
    pub_year: str = "" # Parsed from XML <pub-date>
    pmid: str = ""
    pmcid: str = ""
    xml: str = ""
    parsed_json: Dict[str, Any] = field(default_factory=dict)
    parsed_text: str = ""

@dataclass
class Paper:
    Meta: Paper_MetaData
    Text: Paper_TextData


#############################################################
#  2, Main Class: PubmedFetcher
#############################################################


class PubmedFetcher:
    def __init__(self, root_dir: str, entrez_email: str,api_key: str, batch_size: int = 50, max_retries: int = 3):
        """
        Description
        -----------
        Initializes the PubmedFetcher with output directory, Entrez email, and batch size.
        It automatically sets up the internal directory structure:
        - root_dir/ (the main data repository, can be output directory)
            - papers/  (or any other name you prefer)
                - {pub_year}/
                    - {pmid}/ 
                        - metadata.json
                        - fulltext.txt
                        - others
            - LookUP tables 
            - Others
        
        Args
        -----
        root_dir (str): Directory to save the fetched JSON files.
        entrez_email (str): Email address for NCBI Entrez. It is required by NCBI and should be set to a valid email.
        api_key (str): NCBI API Key for higher rate limits (10 req/sec).
        batch_size (int): Number of articles to fetch per batch, 50~100 is recommended.
        max_retries (int): Maximum number of retries for Entrez API calls. Default is 3.
        
        """
        self.root_dir = root_dir # This is the ROOT of our data repository
        self.batch_size = batch_size
        self.entrez_email = entrez_email
        self.api_key = api_key
        self.max_retries = max_retries

        # Global setting for Entrez email (necessary for Biopython Entrez)
        Entrez.email = entrez_email
        
        if api_key:
            Entrez.api_key = api_key
            print(f"✅ NCBI API Key set successfully. Rate limit increased to 10 req/s.")
        
        # designed for central library, but not used currently
        # make sure that root directory exists
        # if not os.path.exists(self.root_dir):
        #    os.makedirs(self.root_dir)

    #############################################################
    #  2.1, Query search and fetch PMIDs
    #############################################################

    def query_search(self, query: str) -> Dict[str, Any]:
        """
        Description
        -----------
        Search for the first time in PubMed using the given query to get the count list of results and WebEnv info and QueryKey.

        Args
        -----
        query (str): The search query string.

        Returns
        -------
        Dict[str, Any]: A dictionary containing the count list of results, WebEnv, and QueryKey.

        Notes
        -----
        - 1, why ues retmax=0? Because we only want the count and WebEnv, not the actual IDs ——> we do not need to retrieve the actual IDs at this first stage.
        - 2, usehistory="y" is crucial as it tells the server to remember the search results for later retrieval using WebEnv.
        """
        print(f"Now searching PubMed with query [{query}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
        
        for attempt in range(self.max_retries):
            try:
                # set Entrez.email before making requests
                # Entrez.email = self.entrez_email

                # usehistory="y" to enable WebEnv
                handle = Entrez.esearch(db="pubmed", term=query, retmax=0, usehistory="y")
                results = Entrez.read(handle)
                handle.close()
                
                count = int(results["Count"])
                print(f"found {count} related articles about [{query}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                
                return {
                    "count": count,
                    "query": query,
                    "webenv": results["WebEnv"],
                    "query_key": results["QueryKey"]
                }
            except Exception as e:
                if attempt < self.max_retries - 1:
                    print(f"Search PubMed with query [{query}] failed (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 2s...")
                    time.sleep(2)
                else:
                    print(f"Search PubMed with query [{query}] failed after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                    traceback.print_exc()
                    return {"count": 0}
    

    def get_pubmedIDs_from_query(self, query_meta: Dict[str, Any], retmax: int = 500) -> List[str]:
        """
        Description
        -----------
        Get all PubMed IDs from a query search using the history.

        Args
        -----
        query_meta (str): Metadata from the query search containing count, WebEnv, and QueryKey and query. Just the return 
        value of query_search function.
        retmax (int): Batch size for fetching PubMed IDs to retrieve. Default is 500.

        Returns
        -------
        List[str]: A list of PubMed IDs retrieved from the query.

        Notes
        -----
        - 1, esearch function in Biopython's Entrez module forces us to use term parameter when we want to use history server,
        so we need to pass query_meta["query"] here again even we do not really use it. Note that we actually do not need term parameter as long as we have
        WebEnv and QueryKey. It is just a design quirk of Biopython Entrez.esearch function and we just need a placeholder(占位符).
        """

        count = query_meta.get("count", 0)
        webenv = query_meta.get("webenv", "")
        query_key = query_meta.get("query_key", "")
        query = query_meta.get("query", "")

        if count == 0 or not webenv or not query_key:
            # the reason why we do not check query is that we just need a placeholder for Biopython Entrez.esearch function
            print(f"No PMIDs to fetch due to failure in step query_search at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
            return []

        print(f"Retrieving {count} PMIDs from history server at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")

        all_pmids = []

        # Fetch in batches according to retmax
        for start in range(0, count, retmax):
            end = min(count, start + retmax)
            print(f"Fetching PMIDs {start + 1} to {end} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
            
            for attempt in range(self.max_retries):
                try:
                    handle = Entrez.esearch(
                        db="pubmed",
                        term=query,
                        retstart=start,
                        retmax=retmax,
                        webenv=webenv,
                        query_key=query_key
                    )
                    results = Entrez.read(handle)
                    handle.close()
                    
                    batch_pmids = results.get("IdList", [])
                    all_pmids.extend(batch_pmids)
                    
                    print(f"  -> Retrieved {len(batch_pmids)} PMIDs in this batch.")
                    
                    # Polite delay to avoid overwhelming NCBI servers
                    time.sleep(1)
                    break # Success, exit retry loop
                    
                except Exception as e:
                    if attempt < self.max_retries - 1:
                        print(f"Fetching PMIDs {start + 1} to {end} failed (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 2s...")
                        time.sleep(2)
                    else:
                        print(f"Fetching PMIDs {start + 1} to {end} failed after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                        traceback.print_exc()
                        # Continue to next batch even if this one failed
                        continue
        
        # Final count check
        print(f"Total PMIDs retrieved: {len(all_pmids)} out of {count} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
        return all_pmids


    #############################################################
    #  2.2, Fetch articles and parse metadata: paper_metadata only
    #############################################################


    def fetch_from_query(self, query_meta: Dict[str, Any], output_dir: str = None) -> List[Paper_MetaData]:
        """
        Description
        -----------
        Use the query search metadata to fetch articles in batches and save them as JSON files.

        Args
        -----
        query_meta (Dict[str, Any]): Metadata from the query search containing count, WebEnv, and QueryKey. Just the return \
        value of query_search function.

        Returns
        -----   
        List[Paper_MetaData]: A list of Paper_MetaData objects

        Notes
        -----
        - 1, The best strategy to collect large number of articles: medline parsing + xml parsing

        """
        count = query_meta.get("count", 0)
        webenv = query_meta.get("webenv", "")
        query_key = query_meta.get("query_key", "")
        
        if count == 0 or not webenv or not query_key:
            print(f"No articles to fetch due to failure in step query_search at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
            return []

        all_parsed_articles = []

        # Fetch in batches according to batch_size
        for start in range(0, count, self.batch_size):
            end = min(count, start + self.batch_size)
            print(f"Fetching articles {start + 1} to {end} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
            
            for attempt in range(self.max_retries):
                try:
                    # Use efetch to get detailed information in Medline format
                    # AS for rettype and retmode, please refer to: 
                    # https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly
                    fetch_handle_medline = Entrez.efetch(
                        db="pubmed",
                        rettype="medline",
                        retmode="text",
                        retstart=start,
                        retmax=self.batch_size,
                        webenv=webenv,
                        query_key=query_key
                    )
                    # and we also fetch xml format data
                    fetch_handle_xml = Entrez.efetch(
                        db="pubmed",
                        retmode="xml",
                        retstart=start,
                        retmax=self.batch_size,
                        webenv=webenv,
                        query_key=query_key
                    )

                    # Parse Medline format data and xml format data
                    data_medline = Medline.parse(fetch_handle_medline)
                    # data_Medline now is a generator, convert to list will get all data
                    # Now records_medline is a list of batch_size medline records
                    records_medline = list(data_medline)

                    # ⚠️ For xml, We read instead of parsing here, so we must be careful about the batch size and memory usage!
                    data_xml = Entrez.read(fetch_handle_xml)
                    # data_xml now is a huge dict, it only contains 2 keys:'PubmedBookArticle'(none) and 'PubmedArticle'(useful), we use latter
                    articles_list = data_xml['PubmedArticle']

                    fetch_handle_medline.close()
                    fetch_handle_xml.close() 
                    
                    print(f"  -> Retrieved {len(records_medline)} Medline records and {len(articles_list)} Xml articles. Please check whether they equal and the efetch number here with esearch count.")
                    
                    # Now we start to parse each record, note that each record is a PubMed article in Medline format
                    parsed_articles: List[Paper_MetaData] = []

                    # we parse medline and xml at the same time
                    for record,article in zip(records_medline, articles_list):
                        # one single record
                        parsed_medline = self.parse_medline_record(record)
                        parsed_xml = self.parse_article_xml(article)

                        if parsed_medline:
                            # merge medline and xml parsed results
                            # ⚠️ Note that parsed_medline is the main body, we only update contributors from parsed_xml
                            parsed_medline.contributors.xml = parsed_xml
                            parsed_articles.append(parsed_medline)

                    
                    # --- Fetch Linked Data ---
                    if parsed_articles:
                        current_batch_pmids = [p.identity.pmid for p in parsed_articles if p.identity.pmid] # list of pmids
                        
                        # Batch fetch all links for the current batch of PMIDs
                        links_map = self.fetch_linked_data_for_batch_pmid(current_batch_pmids) # Dict[str, PaperLinks]

                        # Map back to Paper objects
                        for p in parsed_articles:
                            if p.identity.pmid in links_map:
                                p.links = links_map[p.identity.pmid]
                    
                    # Save this batch as individual JSON files per paper
                    # We do not save batch files anymore since it is hard to manage and update single paper
                    if parsed_articles:
                        # Change to save single paper to json for better indexing and retrieval
                        for paper_meta in parsed_articles:
                            self.save_single_paper_to_json(paper_meta, output_dir=output_dir)
                        all_parsed_articles.extend(parsed_articles)
                    
                    # Polite delay to avoid overwhelming NCBI servers
                    time.sleep(1)
                    break # Success, exit retry loop
                    
                except Exception as e:
                    if attempt < self.max_retries - 1:
                        print(f"Fetching articles {start + 1} to {end} failed (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 1s...")
                        time.sleep(1)
                    else:
                        print(f"Fetching articles {start + 1} to {end} failed after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                        traceback.print_exc()
                        # Continue to next batch
                        pass
        
        return all_parsed_articles

    def fetch_from_pmid_list(self, pmid_list: List[str], output_dir: str = None) -> List[Paper_MetaData]:
        """
        Description
        -----------
        Fetch articles by a specific list of PMIDs in batches and save them as JSON files.

        Args
        -----
        pmid_list (List[str]): List of PubMed IDs to fetch.

        Returns
        -------
        List[Paper_MetaData]: A list of Paper_MetaData objects retrieved from the given PMIDs.

        Notes
        -----
        - 1, This function is similar to fetch_from_query, but it fetches articles based on a provided list of PMIDs.
        - 2, It can be used to fetch batch articles or a single article by providing a list with one PMID.
        """
        count = len([pmid_list] if isinstance(pmid_list, str) else pmid_list)
        if count == 0:
            print("No PMIDs to fetch.")
            return []

        if isinstance(pmid_list, str):
            pmid_list = [pmid_list]

        print(f"Total PMIDs to fetch: {count} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
        all_parsed_articles = []

        # Fetch in batches according to batch_size, same as fetch_from_query above
        for start in range(0, count, self.batch_size):
            end = min(count, start + self.batch_size)
            batch_pmids = pmid_list[start:end]
            print(f"Fetching articles {start + 1} to {end} (PMID: {batch_pmids}) at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
            
            for attempt in range(self.max_retries):
                try:
                    # same as fetch_from_query, we fetch both medline and xml formats
                    fetch_handle_medline = Entrez.efetch(
                        db="pubmed",
                        rettype="medline",
                        retmode="text",
                        id=batch_pmids
                    )
                    # and we also fetch xml format data
                    fetch_handle_xml = Entrez.efetch(
                        db="pubmed",
                        retmode="xml",
                        id=batch_pmids
                    )

                    # Parse Medline format data and xml format data
                    data_medline = Medline.parse(fetch_handle_medline)
                    records_medline = list(data_medline)

                    data_xml = Entrez.read(fetch_handle_xml)
                    articles_list = data_xml['PubmedArticle']

                    fetch_handle_medline.close()
                    fetch_handle_xml.close()

                    print(f"  -> Retrieved {len(records_medline)} Medline records and {len(articles_list)} Xml articles. Please check whether they equal and whether they match the number of this batch.")
                    
                    # Now we start to parse each record, note that each record is a PubMed article in Medline format
                    parsed_articles: List[Paper_MetaData] = []

                    # we parse medline and xml at the same time
                    for record,article in zip(records_medline, articles_list):
                        # one single record
                        parsed_medline = self.parse_medline_record(record)
                        parsed_xml = self.parse_article_xml(article)

                        if parsed_medline:
                            # merge medline and xml parsed results
                            # ⚠️ Note that parsed_medline is the main body, we only update contributors from parsed_xml
                            parsed_medline.contributors.xml = parsed_xml
                            parsed_articles.append(parsed_medline)

                    # --- Fetch Linked Data ---
                    if parsed_articles:

                        links_map = self.fetch_linked_data_for_batch_pmid(batch_pmids)

                        for p in parsed_articles:
                            if p.identity.pmid in links_map:
                                p.links = links_map[p.identity.pmid]

                    # Save this batch as individual JSON files per paper
                    if parsed_articles:
                        # Change to save single paper to json for better indexing and retrieval
                        for paper_meta in parsed_articles:
                            self.save_single_paper_to_json(paper_meta, output_dir=output_dir)
                        all_parsed_articles.extend(parsed_articles)
                    # Polite delay to avoid overwhelming NCBI servers
                    time.sleep(1)
                    break # Success, exit retry loop

                except Exception as e:
                    if attempt < self.max_retries - 1:
                        print(f"Fetching articles {start + 1} to {end} failed (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 1s...")
                        time.sleep(1)
                    else:
                        print(f"Fetching articles {start + 1} to {end} failed after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                        traceback.print_exc()
                        # Continue to next batch
                        pass
        
        return all_parsed_articles

    def parse_medline_record(self, medline_record: Dict[str, Any]) -> Paper_MetaData:
        """
        Description
        -----------
        Parse one single Medline record into a structured dictionary to get relevant fields for construction of Ontology of this article.
        
        Args
        ----
        medline_record: Dict[str, Any]
            A single Medline record as parsed by Bio.Medline.

        Returns
        -------
        Paper_MetaData
       
            A Paper_MetaData object containing structured information extracted from the Medline record.
        
        Notes
        -----
        - 1, For medline record field, refer to http://zatoka.icm.edu.pl/OVIDWEB/fldguide/medline.htm, https://medialab.github.io/sciencescape/medline_utils/
        - 2, the INPUT medline_record should be a single record of medline parsing result:
            handle = Entrez.efetch(db="pubmed",rettype="medline",retmode="text")
            records = Medline.parse(handle)
            medline_records = list(records)
            for medline_record in medline_records: where each medline_record is the input here.
        """

        try:
            # --- 1, Metadata extraction --- 
            pmid = medline_record.get("PMID", "")
            title = medline_record.get("TI", "")

            # ⚠️ get DOI from AID field if available (more reliable than LID)
            # AID and LID fields could be str or list
            doi = ""
            aids = medline_record.get("AID", [])
            if isinstance(aids, str):
                aids = [aids]
            for aid in aids:
                if "[doi]" in aid:
                    doi = aid.replace("[doi]", "").strip()
                    break
            
            # alternative DOI extraction from LID field if AID not available
            if not doi and 'LID' in medline_record:
                lids = medline_record.get("LID", [])
                if isinstance(lids, str):
                    lids = [lids]
                for lid in lids:
                    if '[doi]' in lid:
                        doi = lid.replace("[doi]", "").strip()
                        break

            # --- 2, Content extraction --- 
            abstract = medline_record.get("AB", "")

            # Mine URLs from abstract text
            mined_urls = extract_urls_from_text(abstract, source_tag="abstract")

            # keywords: OT(Other Terms) field, always converted to list
            keywords = medline_record.get("OT", [])
            if isinstance(keywords, str):
                keywords = [keywords]
            
            # MeSH terms: MH field, always converted to list
            # ⚠️ A Question: what's the difference between OT and MH fields?
            mesh_terms = medline_record.get("MH", [])
            if isinstance(mesh_terms, str):
                mesh_terms = [mesh_terms]

            # Publication Type: PT field, e.g., Journal Article, Review, etc.
            pub_types = medline_record.get("PT", [])
            if isinstance(pub_types, str):
                pub_types = [pub_types]

            # 3, ---- Entity: Authors & Affiliation extraction ---
            # Authors: FAU (Full Author Name) field
            full_author_names = medline_record.get("FAU", [])
            short_author_names = medline_record.get("AU", [])

            # Author Identifiers (e.g., ORCID): AUID field
            auids = medline_record.get("AUID", [])
            # Affiliations: AD field
            affiliations = medline_record.get("AD", [])

            # Ensure all author-related fields are lists
            for i in [full_author_names, short_author_names, auids, affiliations]:
                if isinstance(i, str):
                    i = [i]


            # Construct author entities/dictionaries
            # ⚠️ A Question: How to correctly map authors to their affiliations and identifiers?
            # we can not do this perfectly without more information
            # let's try xml format output

            # 4, ---- Entity: Journal & Publication extraction ---
            journal_title = medline_record.get("JT", ""), # Journal Title
            journal_abbrev = medline_record.get("TA", ""), # Journal Abbreviation

            # Date of Publication: DP field
            pub_date = medline_record.get("DP", "")
            pub_year = pub_date[:4] if len(pub_date) >= 4 else ""

            # Entrez Date: EDAT field, used for incremental updates (录入日期, 用于增量更新)
            entrez_date = medline_record.get("EDAT", "")
                 
            
            # ⚠️ Further fields can be extracted as needed

            return Paper_MetaData(
                identity=PaperIdentity(pmid=pmid, doi=doi, title=title),
                content=PaperContent(abstract=abstract, 
                                     keywords=keywords, 
                                     mesh_terms=mesh_terms, 
                                     pub_types=pub_types),
                contributors=PaperContributors(
                    medline={
                        "full_names": full_author_names,
                        "short_names": short_author_names,
                        "auids": auids,
                        "affiliations": affiliations,
                    }
                ),
                source=PaperSource(
                    journal_title=journal_title,
                    journal_abbrev=journal_abbrev,
                    pub_date=pub_date,
                    pub_year=pub_year,
                    pub_types=pub_types,
                ),
                metadata=PaperMetadata(
                    entrez_date=entrez_date,
                    fetched_at=time.strftime('%Y-%m-%d %H:%M:%S')
                ),
                # Initialize links with mined URLs
                links=PaperLinks(text_mined=mined_urls)
            )
    

        except Exception as e:
            print(f"Error parsing Medline record PMID {medline_record.get('PMID', 'Unknown')}: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
            traceback.print_exc()
            return None

    def parse_article_xml(self, article: Any) -> Optional[List[Dict[str, Any]]]:
        """
        Description
        -----------
        Parse one single article in XML format into a structured dictionary to get relevant fields for construction of Ontology of this article.

        Args
        ----
        article: Any
            A single article as parsed by Entrez.read() from XML format.
        
        Returns
        -------
        List[Dict[str, Any]]
            A list of dictionaries containing structured author information extracted from the article XML.
        
        Notes
        -----
        - 1, the INPUT article should be a single article of xml parsing result:
            handle = Entrez.efetch(db="pubmed",retmode="xml")
            data = Entrez.read(handle)
            articles_list = data['PubmedArticle']
            for article in articles_list: where each article is the input here.
        - 2, Note that we only extract author information here, other fields can be added as needed but medline format already covers most.
        """
        try:
            medline = article.get('MedlineCitation', {})
            article_data = medline.get('Article', {})
            pmid = str(medline.get('PMID', 'Unknown')) # defensive coding: in case to trace back errors

            # --- Authors & Affiliations ---
            authors_structured = []
            full_names_list = []
            short_names_list = []
            auids_list = []
            affiliations_list = []

            if 'AuthorList' in article_data:
                for author in article_data['AuthorList']:
                    # 1. Name Extraction
                    last_name = author.get('LastName', '')
                    fore_name = author.get('ForeName', '')
                    initials = author.get('Initials', '')
                    
                    # Construct Full Name (Medline FAU format: LastName, ForeName)
                    # else if only one of them exists, use that
                    full_name = f"{last_name}, {fore_name}" if last_name and fore_name else last_name or fore_name
                    if full_name: full_names_list.append(full_name)
                    
                    # Construct Short Name (Medline AU format: LastName Initials)
                    short_name = f"{last_name} {initials}" if last_name and initials else last_name or initials
                    if short_name: short_names_list.append(short_name)

                    # 2. Identifier Extraction (e.g. ORCID)
                    current_auids = []
                    if 'Identifier' in author:
                        # Note that author['Identifier'] could be a empty list or a list of identifiers
                        identifiers = author['Identifier']
                        # Defensive coding: handle both list and single item
                        if len(identifiers) > 0:
                            # it means that this author has at least one identifier
                            current_auids.append(str(identifiers[0]))

                    # 3. Affiliation Extraction
                    current_affiliations = []
                    if 'AffiliationInfo' in author:
                        aff_info = author['AffiliationInfo']
                        for aff in aff_info:
                            aff_text = aff.get('Affiliation', '')
                            if aff_text:
                                current_affiliations.append(aff_text)
                    
                    authors_structured.append({
                        "full_name": full_name,
                        "short_name": short_name,
                        "identifiers": current_auids,
                        "affiliations": current_affiliations,
                    })
            
            return authors_structured

        except Exception as e:
            print(f"Error parsing article XML PMID {pmid}: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
            traceback.print_exc()
            return None

    def fetch_linked_data_for_single_pmid(self, pmid: str) -> Dict[str, PaperLinks]:
        """
        Description
        -----------
        Uses ELink to fetch all connected data for a single PMID.
        
        Args
        ----
        pmid: str
            A single PubMed ID to fetch linked data for.

        Returns
        -------
        Dict[str, PaperLinks]
            A dictionary mapping the PMID to its corresponding PaperLinks object containing linked data.

        Notes
        -----
        - 1, ELink Strategy:
            1. 'acheck' to discover available internal Entrez links (Gene, Protein, etc.).
            2. 'neighbor' to fetch the actual IDs for those internal links. (1 and 2 are combined here)
            3. 'llinks' to fetch external URLs (LinkOuts) for datasets, full text, etc. (3 is separate)
        """
        if not pmid: 
            return {}
        
        # Map PMID -> PaperLinks
        links_map = {pmid: PaperLinks()}
        
        # --- Part A: Internal Entrez Links (Discovery & Fetching) ---
        # 1. Discovery (acheck)
        # Note: acheck returns LinkSetDbHistory
        acheck_results = []
        for attempt in range(self.max_retries):
            try:
                handle = Entrez.elink(dbfrom="pubmed", id=pmid, cmd="acheck")
                acheck_results = Entrez.read(handle)
                handle.close()
                break
            except Exception as e:
                if attempt < self.max_retries - 1:
                    print(f"Error in ELink acheck for {pmid} (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 2s...")
                    time.sleep(2)
                else:
                    print(f"Error in ELink acheck for {pmid} after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                    # Continue even if acheck fails
                    acheck_results = [] 

        # Collect all unique LinkNames to query
        # tuple() like (db, linkname)
        links_to_fetch : set[Tuple[str, str]] = set()
        
        # ⚠️ The following links are always useful to fetch
        # cited_by /被他引
        links_to_fetch.add(("pubmed", "pubmed_pubmed_citedin"))
        # references /参考文献
        links_to_fetch.add(("pubmed", "pubmed_pubmed_refs"))
        # similar articles /相似文章
        links_to_fetch.add(("pubmed", "pubmed_pubmed"))
        # related reviews / 相关综述文章
        # ⚠️ But actually, We recommend you to unite all the articles about your PMID papers, and then efetch them using medline format, 
        # and finally filter them by PT field to get all reviews, which is more reliable.
        # So we do not guarantee that this link is complete about all reviews.
        links_to_fetch.add(("pubmed", "pubmed_pubmed_reviews"))
        # pmc free full text articles / PMC 免费全文
        links_to_fetch.add(("pmc", "pubmed_pmc"))

        if acheck_results:
            try:
                # safe access to link info
                id_check_list = acheck_results[0].get('IdCheckList', {})
                if not id_check_list:
                    raise ValueError("IdCheckList missing in acheck results")
                id_link_set = id_check_list.get('IdLinkSet', [])
                if not id_link_set:
                    raise ValueError("IdLinkSet missing in acheck results IdCheckList")
                link_list = id_link_set[0].get('LinkInfo', [])
                if link_list:
                    for link in link_list:
                        db = link.get('DbTo', '')
                        linkname = link.get('LinkName', '')

                        # Filter out useless UI links
                        # ⚠️ You may customize this filtering based on your needs
                        # After testing, we find that some links are not useful for ontology construction, such as:
                        # pubmed_pubmed_five (part of pubmed_pubmed), pubmed_pubmed_reviews_five (part of pubmed_pubmed_reviews)
                        # ExternalLink (link to external resources BUT returns nothing ?))
                        filter_list = ["pubmed_pubmed_five", "pubmed_pubmed_reviews_five", "ExternalLink"]
                        if linkname in filter_list:
                            continue 

                        if db and linkname:
                            links_to_fetch.add((db, linkname))
            except Exception as e:
                traceback.print_exc()
                print(f"    Warning: Failed to parse acheck results for {pmid}: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
        
        print(f"  -> Deep mining {len(links_to_fetch)} types of internal connections for {pmid} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")

        # 2. Fetching (neighbor)
        for link in links_to_fetch:
            if link[0] == "LinkOut":
                # ("LinkOut", "ExternalLink") is handled again
                continue

            # for each link type, we initialize the db_link structure
            db_link = {"id": pmid, "db": link[0], "linkname": link[1], "links": []}
            
            for attempt in range(self.max_retries):
                try:
                    print(f"     Fetching {link[1]} from {link[0]} for {pmid} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                    handle = Entrez.elink(dbfrom="pubmed", id=pmid, db=link[0], linkname=link[1])
                    results = Entrez.read(handle)
                    handle.close()
                    
                    if not results: 
                        break # Success (empty result is valid)
                    
                    # Note that link[1] is the linkname, link[0] is the db
                    # ⚠️ For different db and linkname, THE retured results structure may vary slightly

                    # For db="pubmed" or "pmc", we can directly extract the linked PMIDs/PMCIDs, their uids are in results[0]['LinkSetDb'][0]['Link']
                    link_set_db = results[0].get('LinkSetDb', [])
                    if not link_set_db:
                        # raise ValueError("LinkSetDb missing in elink results")
                        # If missing, maybe just no links
                        break
                    
                    link_data = link_set_db[0].get('Link', [])
                    if not link_data:
                        # raise ValueError("Link missing in LinkSetDb in elink results")
                        break
                        
                    for uid in link_data:
                        if  'Id' in uid:
                            db_link['links'].append(uid['Id'])

                    # Now we have fetched all linked IDs for this link type
                    if db_link['linkname'] == "pubmed_pubmed_citedin":
                        links_map[pmid].cites = db_link['links']
                    elif db_link['linkname'] == "pubmed_pubmed_refs":
                        links_map[pmid].refs = db_link['links']
                    elif db_link['linkname'] == "pubmed_pubmed":
                        links_map[pmid].similar = db_link['links']
                    elif db_link['linkname'] == "pubmed_pubmed_reviews":
                        links_map[pmid].review = db_link['links']
                    elif db_link['linkname'] == "pubmed_pmc":
                        links_map[pmid].pmc = db_link['links']
                    else:
                        # other internal links will be stored in entrez dict
                        links_map[pmid].entrez[db_link['linkname']] = db_link['links']
                    
                    break # Success, exit retry loop

                except Exception as e:
                    if attempt < self.max_retries - 1:
                        print(f"    Warning: Failed to fetch {link[1]} from {link[0]} for {pmid} (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 2s...")
                        time.sleep(2)
                    else:
                        # since we force fetch some useful links above, some of them may not be available for certain PMIDs
                        # so we may fail to fetch them at db_link['links'] update, just warn and continue
                        print(f"    Warning: Failed to fetch {link[1]} from {link[0]} for {pmid} after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                        # If fetching fails in db_link['links'] update, we just continue, cause links_map[pmid].this_link is still an available empty list, we still can access it later
                        traceback.print_exc()
                        pass

                
        # --- Part B: External LinkOuts (llinks) ---
        print(f"  -> Fetching external LinkOuts (Datasets, Full Text, etc.) for {pmid} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
        
        for attempt in range(self.max_retries):
            try:
                handle = Entrez.elink(dbfrom="pubmed", id=pmid, cmd="llinks")
                llinks_results = Entrez.read(handle, validate=False)
                handle.close()
                
                if llinks_results:
                    # safe access to IdUrlList
                    id_url_list = llinks_results[0].get('IdUrlList', {})
                    if not id_url_list:
                        # raise ValueError("IdUrlList missing in llinks results")
                        break
                    id_url_set = id_url_list.get('IdUrlSet', [])
                    if not id_url_set:
                        # raise ValueError("IdUrlSet missing in IdUrlList in llinks results")
                        break
                    urls_list = id_url_set[0].get('ObjUrl', [])
                    if urls_list:
                        links_map[pmid].external = urls_list
                
                break # Success, exit retry loop
                
            except Exception as e:
                if attempt < self.max_retries - 1:
                    print(f"    Warning: Failed to fetch external LinkOuts for {pmid} (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 2s...")
                    time.sleep(2)
                else:
                    print(f"    Warning: Failed to fetch external LinkOuts for {pmid} after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                    traceback.print_exc()    

        return links_map

    def fetch_linked_data_for_batch_pmid(self, pmids: List[str]) -> Dict[str, PaperLinks]:
        """
        Description
        -----------
        Uses ELink to fetch all connected data for a batch of PMIDs.
        
        Args
        ----
        pmids: List[str]
            A list of PubMed IDs to fetch linked data for.

        Returns
        -------
        Dict[str, PaperLinks]
            A dictionary mapping each PMID to its corresponding PaperLinks object containing linked data.

        Notes
        -----
        - 1, This function internally calls fetch_linked_data_single for each PMID in the list.
        - 2, This batch processing helps to reduce the number of individual requests to NCBI servers.
        """

        if not pmids:
            return {}
        
        # Map PMID -> PaperLinks
        links_map = {pmid: PaperLinks() for pmid in pmids}

        # --- Part A: Internal Entrez Links (Discovery & Fetching) ---
        # 1. Discovery (acheck)
        # Note: acheck returns LinkSetDbHistory
        acheck_results = []
        for attempt in range(self.max_retries):
            try:
                handle = Entrez.elink(dbfrom="pubmed", id=pmids, cmd="acheck")
                acheck_results = Entrez.read(handle)
                handle.close()
                break
            except Exception as e:
                if attempt < self.max_retries - 1:
                    print(f"Error in batch Elink acheck for {len(pmids)} PMIDs (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 1s...")
                    time.sleep(1)
                else:
                    # we do not print traceback here to avoid flooding the logs
                    print(f"    [Warning] Batch Elink acheck failed for {len(pmids)} PMIDs after {self.max_retries} attempts at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...  Skipping discovery step (using default links). Error: {e}")
                    # Continue even if acheck fails
                    acheck_results = []


        # Collect all unique LinkNames to query
        # tuple() like (db, linkname)
        links_to_fetch : set[Tuple[str, str]] = set()

        # ⚠️ The following links are always useful to fetch
        # cited_by /被他引
        links_to_fetch.add(("pubmed", "pubmed_pubmed_citedin"))
        # references /参考文献
        links_to_fetch.add(("pubmed", "pubmed_pubmed_refs"))
        # similar articles /相似文章
        links_to_fetch.add(("pubmed", "pubmed_pubmed"))
        # related reviews / 相关综述文章
        links_to_fetch.add(("pubmed", "pubmed_pubmed_reviews"))
        # pmc free full text articles / PMC 免费全文
        links_to_fetch.add(("pmc", "pubmed_pmc"))

        if acheck_results:
            for linkset in acheck_results:
                try:
                    # safe access to LinkInfo
                    id_check_list = linkset.get('IdCheckList', {})
                    if not id_check_list:
                        continue
                    id_link_set = id_check_list.get('IdLinkSet', [])
                    if not id_link_set:
                        continue
                    link_list = id_link_set[0].get('LinkInfo', [])
                    if link_list:
                        for link in link_list:
                            db = link.get('DbTo', '')
                            linkname = link.get('LinkName', '')

                            # Filter out useless UI links
                            filter_list = ["pubmed_pubmed_five", "pubmed_pubmed_reviews_five", "ExternalLink"]
                            if linkname in filter_list:
                                continue 
                            if "combined" in linkname:
                                continue

                            if db and linkname:
                                links_to_fetch.add((db, linkname))
            
                except Exception as e:
                    print(f"    Warning: Failed to parse batch acheck results for {len(pmids)} PMIDs: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                    traceback.print_exc()

        print(f"  -> Deep mining {len(links_to_fetch)} types of internal connections for {len(pmids)} PMIDs at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")

        # 2. Fetching (neighbor)
        for link in links_to_fetch:
            if link[0] == "LinkOut":
                # ("LinkOut", "ExternalLink") is handled again
                continue

            # for each link type, we initialize the db_link structure
            for attempt in range(self.max_retries):
                try:
                    print(f"     Fetching {link[1]} from {link[0]} for {len(pmids)} PMIDs at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                    handle = Entrez.elink(dbfrom="pubmed", id=pmids, db=link[0], linkname=link[1])
                    results = Entrez.read(handle)
                    handle.close()
                    
                    if not results: 
                        break

                    # For each PMID in results
                    for linkset in results:
                        source_id = "Unknown"
                        try:
                            source_id = str(linkset['IdList'][0]) # the source PMID
                            
                            if source_id not in links_map:
                                continue
                            
                            # for each link type, we initialize the db_link structure
                            db_link = {"id": source_id, "db": link[0], "linkname": link[1], "links": []}

                            # For db="pubmed" or "pmc", we can directly extract the linked PMIDs/PMCIDs, their uids are in linkset['LinkSetDb'][0]['Link']
                            link_set_db = linkset.get('LinkSetDb', [])
                            # defensive coding
                            if not link_set_db:
                                continue
                            if len(link_set_db) == 0:
                                continue

                            links_list = link_set_db[0].get('Link', [])
                            if not links_list:
                                continue
                            
                            for uid in links_list:
                                if 'Id' in uid:
                                    db_link["links"].append(uid['Id'])

                            # Now we have fetched all linked IDs for this link type, for this PMID
                            if db_link["linkname"] == "pubmed_pubmed_citedin":
                                links_map[source_id].cites = db_link['links']
                            elif db_link["linkname"] == "pubmed_pubmed_refs":
                                links_map[source_id].refs = db_link['links']
                            elif db_link["linkname"] == "pubmed_pubmed":
                                links_map[source_id].similar = db_link['links']
                            elif db_link["linkname"] == "pubmed_pubmed_reviews":
                                links_map[source_id].review = db_link['links']
                            elif db_link["linkname"] == "pubmed_pmc":
                                links_map[source_id].pmc = db_link['links']
                            else:
                                # other internal links will be stored in entrez dict
                                links_map[source_id].entrez[db_link['linkname']] = db_link['links']

                        except Exception as e:
                            print(f"    Warning: Failed to parse fetched links for {link[1]} from {link[0]} for {source_id}: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                            traceback.print_exc()
                            continue
                    
                    break # Success, exit retry loop

                except Exception as e:
                    if attempt < self.max_retries - 1:
                        print(f"    Warning: Failed to fetch {link[1]} from {link[0]} for {len(pmids)} PMIDs (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 1s...")
                        time.sleep(1)
                    else:
                        print(f"    Warning: Failed to fetch {link[1]} from {link[0]} for {len(pmids)} PMIDs after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                        traceback.print_exc()
                        pass

        # --- Part B: External LinkOuts (llinks) ---
        print(f"  -> Fetching external LinkOuts (Datasets, Full Text, etc.) for {len(pmids)} PMIDs at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
        
        for attempt in range(self.max_retries):
            try:
                handle = Entrez.elink(dbfrom="pubmed", id=pmids, cmd="llinks")
                llinks_results = Entrez.read(handle, validate=False)
                handle.close()
                
                if llinks_results:
                        for linkset in llinks_results:
                            source_id = "Unknown"
                            try:
                                # safe access to IdUrlList
                                id_url_list = linkset.get('IdUrlList', {})
                                if not id_url_list:
                                    continue
                                id_url_set = id_url_list.get('IdUrlSet', [])
                                if not id_url_set:
                                    continue
                                source_id = str(id_url_set[0].get('Id', 'Unknown')) # the source PMID
                                if source_id not in links_map:
                                    continue
                                
                                # safe access to ObjUrl
                                urls_list = id_url_set[0].get('ObjUrl', [])
                                if urls_list:
                                    # we need clean urls here
                                    clean_urls = []

                                    for item in urls_list:
                                        # 1. Handle Category (it's a list, take first item)
                                        cat_list = item.get('Category', [])
                                        category = str(cat_list[0]) if isinstance(cat_list, list) and cat_list else ""

                                        # 2. Handle Attribute (it's a list, take first item)
                                        attr_list = item.get('Attribute', [])
                                        attribute = str(attr_list[0]) if isinstance(attr_list, list) and attr_list else ""
                                        
                                        # 3. Handle LinkName (optional string, defensive)
                                        linkname = str(item.get('LinkName', ''))

                                        # 4. Construct clean url entry
                                        url = str(item.get('Url', ''))
                                    
                                        # 5. handle Provider Name
                                        provider_name = str(item.get('Provider', {}).get('Name', "Unknown"))

                                        clean_urls.append({
                                            "url": url,
                                            "provider": provider_name,
                                            "category": category,
                                            "attribute": attribute,
                                            "linkname": linkname
                                        })
                                    links_map[source_id].external = clean_urls

                            except Exception as e:
                                print(f"    Warning: Failed to parse external LinkOuts for {source_id}: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                                continue
                break # Success, exit retry loop
                
            except Exception as e:
                if attempt < self.max_retries - 1:
                    print(f"    Warning: Failed to fetch external LinkOuts for {len(pmids)} PMIDs (Attempt {attempt+1}/{self.max_retries}): [{e}]. Retrying in 1s...")
                    time.sleep(1)
                else:
                    print(f"    Warning: Failed to fetch external LinkOuts for {len(pmids)} PMIDs after {self.max_retries} attempts: [{e}] at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                    traceback.print_exc()

        return links_map


    # ⚠️ For batch saving, we deprecate it for now, please use single paper saving instead.
    # def save_batch_to_json(self, articles: List[Paper], start_index: int)

    def save_single_paper_to_json(self, paper_meta: Paper_MetaData, output_dir: Optional[str] = None):
        """
        Description
        -----------
        Save a single paper (metadata only) to JSON.
        
        Args
        ----
        paper_meta: Paper_MetaData
            One single Paper_MetaData object to save.
        output_dir: Optional[str]
            Directory to save the JSON file. Defaults to self.root_dir (the library).
        """
        pmid = paper_meta.identity.pmid if paper_meta.identity.pmid else "unknown"
        filename = f"{pmid}.json"
        
        # Determine base directory: specific output_dir or default repo root_dir
        base_dir = output_dir if output_dir else self.root_dir
        
        # Structure: base_dir / Year / PMID / pmid.json
        pub_year = paper_meta.source.pub_year if paper_meta.source.pub_year else "Unknown_Year"
        
        # Create year directory
        year_path = os.path.join(base_dir, pub_year)
        
        # Create PMID directory inside year directory
        pmid_path = os.path.join(year_path, pmid)

        if not os.path.exists(pmid_path):
            os.makedirs(pmid_path, exist_ok=True)
            
        final_filepath = os.path.join(pmid_path, filename)
        
        with open(final_filepath, 'w') as f:
            json.dump(paper_meta.to_dict(), f, ensure_ascii=False, sort_keys=True, indent=4)
            
        print(f"  -> Saved {pmid} metadata to {final_filepath}")


    #############################################################
    #  2.3, Fetch Full Text from PMC: paper_text_data only
    #############################################################

    def fetch_pmc_full_text(self, pmid_list: Union[str, List[str]], output_dir: str = None) -> List[Paper_TextData]:
        """
        Description
        -----------
        Fetch full text from PMC for given PMIDs using detailed error handling and batch processing.
        
        Reasons for Try-Except Blocks:
        1. ELink/EFetch: Network calls are unstable. We use retry loops with backoff.
        2. Parsing: XML structure is unpredictable. We wrap individual article parsing in try-except to ensures 
           that one malformed article doesn't crash the entire batch.

        Args
        ---- 
        pmid_list: Union[str, List[str]]
            A single PMID (e.g., "12345678") or a list of PMIDs.
            pmids will be mapped to PMC IDs (e.g., "PMC8328303" or "8328303")
        output_dir: str
            Directory to save the full text XML and parsed content (markdown and JSON). Defaults to self.root_dir. 
            
        Returns
        -------
        List[Paper_TextData]
            A list of Paper_TextData objects containing the full text content for each PMID, including raw XML content, parsed json content, and parsed text content.
        """

        # 1. first check
        count = len([pmid_list] if isinstance(pmid_list, str) else pmid_list)
        if count == 0:
            print("No PMIDs provided for PMC full text fetching.")
            return []

        if isinstance(pmid_list, str):
            pmid_list = [pmid_list]
            
        all_paper_text_data = []

        print(f"Fetching full text for {count} Pubmed articles at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
        
        # 2. Batch Processing Loop (Outer Loop)
        # We process PMIDs in batches to respect API limits and manage memory.
        for start in range(0, count, self.batch_size):
            end = min(count, start + self.batch_size)
            batch_pmids = pmid_list[start:end]
            print(f" -> Converting Pubmed articles {start+1} to {end} (PMID : {batch_pmids}) to PMC IDs at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")

            # PMID -> PMCID Mapping (with Retry)
            pmid_to_pmcid = {}
            valid_pmcids = []
            
            # --- Retry Block for ELink ---
            for attempt in range(self.max_retries):
                try:
                    # Map PMIDs to PMC IDs using ELink
                    handles = Entrez.elink(dbfrom="pubmed", id=batch_pmids, db="pmc", linkname="pubmed_pmc")
                    pmc_links = Entrez.read(handles)
                    handles.close()

                    for linkset in pmc_links:
                        # Defensive coding: Check if IdList exists (input ID)
                        if not linkset.get('IdList'): 
                            continue
                        source_id = str(linkset['IdList'][0]) 
                        
                        link_set_db = linkset.get('LinkSetDb', [])
                        if link_set_db:
                            pmc_links_list = link_set_db[0].get('Link', [])
                            if pmc_links_list:
                                pmc_id = str(pmc_links_list[0]['Id']) 
                                valid_pmcids.append(pmc_id)
                                pmid_to_pmcid[pmc_id] = source_id # Map PMC ID back to PMID for later identification
                    
                    break # Success: Break the retry loop immediately
                    
                except Exception as e:
                    if attempt < self.max_retries - 1:
                        print(f"     [Warning] Map PMIDs to PMC IDs failed (Attempt {attempt+1}/{self.max_retries}): {e}, Retrying in 1s...")
                        time.sleep(1) # Backoff strategies: wait before retry
                    else:
                        print(f"     [Error] Map PMIDs to PMC IDs failed after {self.max_retries} attempts: {e} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                        traceback.print_exc()
                        # skip to next batch if mapping fails completely
                        continue
            
            # Reason for Continue: If mapping fails entirely, we cannot fetch anything for this batch.
            # We skip to the next batch of PMIDs.
            if not valid_pmcids:
                print(f"  -> No valid PMC IDs found for current batch of PMIDs: {batch_pmids}. Skipping full text fetching for this batch at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                continue

            # Fectch Full Text for valid PMC IDs only
            print(f"  -> Mapped {len(valid_pmcids)} out of {len(batch_pmids)} PMIDs to valid PMC IDs. Downloading full text XML for these PMC IDs at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                
            for attempt in range(self.max_retries):
                try:
                    # db="pmc", retmode="xml" gets the full structured text
                    handles = Entrez.efetch(db="pmc", id=valid_pmcids, retmode="xml")
                    pmc_full_xml = handles.read() # Read raw string
                    handles.close()
                    
                    break # Success: Break the retry loop immediately

                except Exception as e:
                    if attempt < self.max_retries - 1:
                        print(f"     [Warning] EFetch full text XML failed (Attempt {attempt+1}/{self.max_retries}): {e}, Retrying in 1s...")
                        time.sleep(1) # Backoff strategies: wait before retry
                    else:
                        print(f"     [Error] EFetch full text XML failed after {self.max_retries} attempts: {e} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                        traceback.print_exc()
                        # skip to next batch if fetching fails completely
                        continue

            # Parse Big XML using BeautifulSoup
            # This part handles the "Big XML" containing multiple articles.
            try:
                soup = BeautifulSoup(pmc_full_xml, 'xml')
                articles = soup.find_all('article')
                
                # Iterate over each article found in the batch XML
                for article in articles:
                    # Reason for Try-Except inside loop:
                    # Isolation. If one article has malformed XML or unexpected structure, 
                    # we log the error and CONTINUE to the next article, rather than crashing the whole batch.
                    try:
                        # 1. Identify the paper (Double Check)
                        # We need to know which PMID this XML belongs to.
                        current_pmid = None
                        
                        # Strategy 1: Look for explicit PMID in metadata
                        pmid_node = article.find("article-id", {"pub-id-type": "pmid"})
                        if pmid_node:
                            current_pmid = pmid_node.text.strip()
                        
                        # Strategy 2: Look for PMC ID and map back using our pmid_to_pmcid map
                        if not current_pmid:
                            pmc_node = article.find("article-id", {"pub-id-type": "pmcid"})
                            if pmc_node:
                                found_pmc = pmc_node.text.strip()
                                # The map keys are usually raw IDs, check compatibility
                                if found_pmc in pmid_to_pmcid:
                                    current_pmid = pmid_to_pmcid[found_pmc]
                                    
                        if not current_pmid:
                            # If we can't identify the paper, we can't safely store it.
                            continue
                            
                        # 2. Extract Full Text (Hierarchical JSON structure)
                        # Utilizing the BeautifulSoup logic
                        # Json file is hierarchical representation of the article structure, suitable for structured analysis
                        # For example, we can extract the same section between different papers easily and compare them
                        parsed_json = self._parse_soup_to_json(article)
                        
                        # 3. Flatten for text view (for NLP or reading)
                        # Plain text is a flattened version suitable for reading or NLP tasks, we use markdown-like formatting
                        # Just markdown file 
                        # Markdown file is suitable for human reading and simple text-based NLP tasks, and easy to be used in AI models
                        parsed_text = self._flatten_json_to_text(parsed_json)

                        # 4. Create Data Object
                        current_pmc = ""
                        pmc_node = article.find("article-id", {"pub-id-type": "pmcid"})
                        if pmc_node: 
                            current_pmc = pmc_node.text.strip()
                        
                        # 5. Extract Publication Year from <pub-date>
                        current_pub_year = "Unknown_Year"
                        pub_date = article.find("pub-date").find("year")
                        if pub_date:
                            current_pub_year = pub_date.text.strip()

                        paper_text_data = Paper_TextData(
                            pmid=current_pmid,
                            pmcid=current_pmc,
                            xml=str(article), # Store just this article's XML
                            parsed_json=parsed_json,
                            parsed_text=parsed_text,
                            pub_year=current_pub_year # ADDED from XML
                        )
                        
                        # save this paper's text data to output directory
                        self.save_single_paper_text(paper_text_data, output_dir=output_dir)

                        all_paper_text_data.append(paper_text_data)
                        
                    except Exception as inner_e:
                        print(f"     [Error] Failed to parse one article in batch: {inner_e} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                        traceback.print_exc()
                        continue
                        
            except Exception as e:
                # This catches errors if the entire XML batch is invalid (unlikely but possible)
                print(f"     [Error] Failed to parse XML batch with BeautifulSoup: {e} at [{time.strftime('%Y-%m-%d %H:%M:%S')}] ...")
                traceback.print_exc()
                continue
                
        return all_paper_text_data
    

    def _parse_soup_to_json(self, article_soup: Any) -> Dict[str, Any]:
        """
        Parses a single <article> soup object into hierarchical JSON.
        Robust strategy:
        1. Parse Title.
        2. Parse Abstract (Handle structured <sec> or flat text).
        3. Parse Body (Recursive <sec>).
        4. Parse Back (Recursive <sec> for Acknowledgements, etc.).
        """
        paper_structure = {
            "title": "N/A",
            "body": []
        }
        
        # 1. Title
        title_node = article_soup.find('article-title')
        if title_node:
            paper_structure["title"] = title_node.get_text().strip()

        # 2. Abstract
        abstract_node = article_soup.find('abstract')
        if abstract_node:
            # Check for structured sections in abstract
            # Note: recursive=False is safer to avoid finding sections inside other things erroneously
            abs_sections = abstract_node.find_all("sec", recursive=False)
            if abs_sections:
                for sec in abs_sections:
                    parsed_sec = self._parse_section_recursive(sec)
                    if parsed_sec["title"] == "N/A":
                        parsed_sec["title"] = "Abstract Section"
                    paper_structure["body"].append(parsed_sec)
            else:
                # If no top-level sections, it might be flat paragraphs or unstructured
                # Sometimes structured abstracts use <title> directly without <sec> (rare but possible)
                
                # Check for title
                abs_title_node = abstract_node.find('title')
                abs_title = abs_title_node.get_text().strip() if abs_title_node else "Abstract"
                
                # Collect paragraphs
                ps = abstract_node.find_all('p')
                if ps:
                    content = [p.get_text().strip() for p in ps]
                    paper_structure["body"].append({
                        "title": abs_title,
                        "content": content,
                        "subsections": []
                    })
                else:
                    # Fallback to full text if no p tags
                    text = abstract_node.get_text(separator=' ', strip=True)
                    if text:
                        paper_structure["body"].append({
                            "title": abs_title,
                            "content": [text], 
                            "subsections": []
                        })

        # 3. Body
        body = article_soup.find('body')
        if body:
            # Standard JATS: <body > <sec> ... </sec> </body>
            body_sections = body.find_all("sec", recursive=False)
            if body_sections:
                for sec in body_sections:
                     # Prevent parsing sub-sections as top-level
                    if sec.parent and sec.parent.name == 'sec':
                         continue
                    paper_structure["body"].append(self._parse_section_recursive(sec))
            else:
                # No sections? Look for direct paragraphs (e.g. Letters)
                direct_paragraphs = body.find_all("p", recursive=False)
                if direct_paragraphs:
                    content = [p.get_text().strip() for p in direct_paragraphs]
                    paper_structure["body"].append({
                        "title": "Main Text", 
                        "content": content,
                        "subsections": []
                    })

        # 4. Back Matter (Acknowledgements, etc.)
        back = article_soup.find('back')
        if back:
            # <ack> <sec> ...
            # Process standard sections in back
            back_sections = back.find_all("sec", recursive=False)
            for sec in back_sections:
                parsed = self._parse_section_recursive(sec)
                paper_structure["body"].append(parsed)
                
            # Specifically check for <ack>
            ack = back.find('ack')
            if ack:
                # ack usually contains title and p
                ack_title_node = ack.find('title')
                ack_title = ack_title_node.get_text().strip() if ack_title_node else "Acknowledgements"
                ack_ps = ack.find_all('p')
                if ack_ps:
                    content = [p.get_text().strip() for p in ack_ps]
                    paper_structure["body"].append({
                        "title": ack_title,
                        "content": content,
                        "subsections": []
                    })

        return paper_structure

    def _parse_section_recursive(self, sec_element: Any) -> Dict[str, Any]:
        """
        Description
        -----------
        Recursive helper for section parsing.
        Traverses <sec> -> <title>/<p>/<sec> structure.
        Every section is represented as a dictionary with title, content, and its subsections.

        Args
        ----
        sec_element: Any
            A BeautifulSoup element representing a <sec> node.

        Returns
        -------
        Dict[str, Any]
            A dictionary representing the section with title, content, and subsections.
        """
        section_data = {
            "title": "N/A",
            "content": [],
            "subsections": []
        }

        # Title
        title_node = sec_element.find("title", recursive=False)
        if title_node:
            section_data["title"] = title_node.get_text().strip()

        # Content (Paragraphs)
        # Only direct <p> children, not nested in subsections, so we use recursive=False here
        direct_paragraphs = sec_element.find_all("p", recursive=False)
        section_data["content"] = [p.get_text().strip() for p in direct_paragraphs]

        # Subsections (Recursive)
        sub_sections = sec_element.find_all("sec", recursive=False)
        for sub in sub_sections:
            child_data = self._parse_section_recursive(sub)
            section_data["subsections"].append(child_data)

        return section_data
        
    def _flatten_json_to_text(self, json_data: Dict[str, Any]) -> str:
        """
        Converts the hierarchical JSON into a flat text string for simple viewing (Markdown-like).
        """
        lines = []
        
        # Article Title (H1)
        # We explicitly mark it as Title
        article_title = json_data.get('title', 'N/A')
        lines.append(f"# {article_title}")
        lines.append("") # Blank line after title
        
        def recurse(sections, level):
            for sec in sections:
                # Section Title (Markdown Header)
                header_prefix = "#" * level
                title = sec.get('title', 'No Title')
                
                # Add a blank line before section title for better structure
                lines.append("") 
                lines.append(f"{header_prefix} {title}")
                
                # Content (Paragraphs)
                if sec['content']:
                    for para in sec['content']:
                        # Clean up newlines: replace newlines with spaces and strip extra spaces
                        # This fixes issues where citations like [1] cause line breaks
                        clean_para = str(para).replace('\n', ' ').strip()
                        
                        # Indent paragraphs slightly (2 spaces) to visually distinguish from headers
                        lines.append(f"  {clean_para}")
                        lines.append("") # Add blank line after paragraph for readability
                
                # Subsections (Recursive)
                recurse(sec['subsections'], level + 1)
        
        # Body sections start from Level 2 (##), assuming Paper Title is Level 1 (#)
        recurse(json_data.get('body', []), level=2)
        
        return "\n".join(lines)

    def save_single_paper_text(self, paper_text: Paper_TextData, output_dir: Optional[str] = None):
        """      
        Description
        -----------
        Save a single paper's text data (full text) to structured files (XML, Markdown, JSON).

        Args
        ----
        paper_text: Paper_TextData
            One single Paper_TextData object to save. It contains raw XML, parsed JSON, and parsed Markdown.
        output_dir: Optional[str]
            Directory to save the text files. Defaults to self.root_dir (the library).
        
        """

        # Determine base directory: specific output_dir or default repo root_dir
        base_dir = output_dir if output_dir else self.root_dir

        # structure: base_dir / Year / PMID / files
        pmid = paper_text.pmid if paper_text.pmid else "unknown"
        pub_year = paper_text.pub_year if paper_text.pub_year else "Unknown_Year"

        # create year directory
        year_path = os.path.join(base_dir, pub_year)
        # create PMID directory inside year directory
        pmid_path = os.path.join(year_path, pmid)

        if not os.path.exists(pmid_path):
            os.makedirs(pmid_path, exist_ok=True)

        # 1. Save Raw XML
        if paper_text.xml:
            xml_filename = f"{pmid}.xml"
            xml_filepath = os.path.join(pmid_path, xml_filename)
            with open(xml_filepath, 'w') as f:
                f.write(paper_text.xml)
            print(f"  -> Saved XML to {xml_filepath}")

        # 2. Save Parsed JSON
        if paper_text.parsed_json:
            json_filename = f"{pmid}_parsed.json"
            json_filepath = os.path.join(pmid_path, json_filename)
            with open(json_filepath, 'w') as f:
                json.dump(paper_text.parsed_json, f, ensure_ascii=False, indent=2)
            print(f"  -> Saved parsed JSON to {json_filepath}")

        # 3. Save Parsed Text (Markdown)
        if paper_text.parsed_text:
            md_filename = f"{pmid}_parsed.md"
            md_filepath = os.path.join(pmid_path, md_filename)
            with open(md_filepath, 'w') as f:
                f.write(paper_text.parsed_text)
            print(f"  -> Saved parsed text to {md_filepath}")


    #############################################################
    #  2.4, Fetch articles and parse metadata + full text together
    #############################################################
        
    def fetch_and_save_full_papers(self, query: str = None, pmid_list: List[str] = None, output_dir: str = None) -> List[Paper]:
        """
        Description
        -----------
        Integrated fetching of full papers (metadata + full text) based on either a search query or a list of PMIDs.
        we first fetch metadata, then fetch full text from PMC, and finally save both metadata (as JSON) and text (TXT, MARKDOWN, JSON) to structured directories.

        
        Args
        ----
        query: str
            A PubMed search query string. If provided, PMIDs will be fetched based on this query.
        pmid_list: List[str]
            A list of PubMed IDs to fetch. If provided, these PMIDs will be used directly.
        output_dir: str
            Directory to save the fetched papers. If not provided, defaults to self.root_dir.
        
        Returns
        -------
            List[Paper]: A list of Paper objects containing both metadata and text data.

        Notes
        -----
        - 1, we store the results in following structure:
            output_dir / Year / PMID / pmid.json (metadata)
                                    / pmid.txt (raw xml)
                                    / pmid_parsed.md (markdown text)
                                    / pmid_parsed.json (parsed json)

        """
        if not query and not pmid_list:
            raise ValueError("Must provide either query or pmid_list")
            
        # 1. Fetch Metadata
        print("=== Step 1: Fetching Metadata ===")
        meta_data_list: List[Paper_MetaData] = []
        if query:
            query_record = self.query_search(query)
            meta_data_list = self.fetch_from_query(query_record, output_dir) 
        else:
            meta_data_list = self.fetch_from_pmid_list(pmid_list, output_dir)
            
        if not meta_data_list:
            print("No papers found.")
            return []

        # Collect PMIDs for text fetching
        found_pmids = [p.identity.pmid for p in meta_data_list if p.identity.pmid]
        
        # 2. Fetch Full Text (returns List[Paper_TextData])
        # This function fetches, parses (using robust soup parser), and saves text data to Year/PMID/ code.
        print("\n=== Step 2: Fetching Full Text ===")
        text_data_list = self.fetch_pmc_full_text(found_pmids, output_dir=output_dir)
        
        # Convert list to map for easy lookup
        text_map = {td.pmid: td for td in text_data_list if td.pmid}
        
        # 3. Process and Compile
        print("\n=== Step 3: Processing and Saving Metadata ===")
        complete_papers = []
        
        for meta in meta_data_list:
            pmid = meta.identity.pmid
            pm_text = text_map.get(pmid, Paper_TextData())
            
            # Combine
            paper = Paper(Meta=meta, Text=pm_text)
            complete_papers.append(paper)
            
            # 3.1 Extract URLs from full text and update metadata
            if pm_text.parsed_text:
                extracted_urls = extract_urls_from_text(pm_text.parsed_text, "full_text")
                if extracted_urls:
                    print(f"  -> Extracted {len(extracted_urls)} URLs from full text for PMID {pmid}")
                    # Update metadata links
                    if not meta.links:
                        meta.links = PaperLinks()
                    
                    # Merge extracting URLs, avoiding duplicates if possible (though straightforward extend is okay here)
                    meta.links.text_mined.extend(extracted_urls)
            
            # 3.2 Save Metadata (JSON)
            # Text data (xml, parsed_json, parsed_md) is already saved by fetch_pmc_full_text
            # We save metadata NOW (save twice) because it might have been updated with mined links
            self.save_single_paper_to_json(meta, output_dir)
                
        return complete_papers
