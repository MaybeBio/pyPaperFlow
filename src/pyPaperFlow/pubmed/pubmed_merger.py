"""
Paper Merger Module for pyPaperFlow

This module provides functionality to merge PubMed paper metadata and content
into unified formats optimized for AI analysis.

⚠️ For Pubmed papers ONLY

This module provides a minimal, well-documented implementation that
performs exactly two tasks required by the user:

- Merge: scan a PubMed-style folder (pubmed/<year>/<pmid>/...) or a list
  of PMIDs and produce per-paper `{PMID}.json` sidecars plus a single
  merged JSON (or JSONL) file named `<source>_<timestamp>.json`.

- Export: read the merged JSON/JSONL and write a single Markdown file.

"""

import os
import json
import csv
from datetime import datetime
from typing import * 
import yaml


def _timestamp() -> str:
    return datetime.now().strftime('%Y-%m-%d_%H-%M-%S')


class PubmedMerger:
    """
    Description
    -----------
    Minimal merger with two public methods: merge_json and export_md.
    JSON/JSONL is for structured data storage and potential downstream use, 
    while Markdown is for human-readable output and LLM input(the most important part).

    Usage:
        merger = PubmedMerger()
        merger.merge_json_from_directory(paper_dir, output_path, pmid_file=None)
        merger.export_md_from_merged_json(merged_json, output_md, yaml_cfg=None)
    """

    def merge_json_from_directory(
        self,
        paper_dir: str,
        output_json: str,
        pmid_file: Optional[str] = None,
        jsonl: bool = False,
    ) -> Dict[str, Any]:
        """
        Description
        -----------
        Merge per-paper JSONs.

        - If `pmid_file` is provided, only PMIDs in that file are merged. 
        But the premise is that the pmid_file's PMIDs must be present in the directory structure. 
        If `pmid_file` is not provided, all discovered PMIDs in `paper_dir` are merged.
        - For each discovered paper we write `{PMID}.json` next to the
          paper's files (sidecar).
        - A single merged JSON (array) or JSONL is written to `output_json`.


        Args 
        ----
        paper_dir: Directory containing per-paper subdirectories (e.g. pubmed/2023/12345678/, withput pubmed/ containing)
        output_json: Path to the output merged JSON/JSONL file, or a directory to auto-name the output, default is current directory with auto-naming
        pmid_file: Optional path to a file containing a list of PMIDs to merge
        jsonl: If True, write output as JSONL instead of JSON 
  
        Returns
        -------
        A small stats dict.
          
        Example return value:
        {
            'total': 0,  of PMIDs found in directory (after optional filtering by pmid_file)
            'meta_missing': ['12345678'],  of papers that failed to merge in theory (missing meta JSON), we do not count the files failing at merge part
            'content_missing': ['12345678'],  of papers with missing content JSON
        }
        """
     
        pmids: Optional[List[str]] = None
        if pmid_file:
            # read PMIDs from file, if provided (filter), list like ['12345678', '23456789', ...]
            pmids = self._read_pmids_from_file(pmid_file)

        # discover PMIDs and their directories, list like [('12345678', '/path/to/paper_dir/pubmed/2023/12345678'), ...]
        pairs = self._iter_pmid_directories(paper_dir)

        # define stats
        # list for pmid records
        meta_missing = []
        content_missing = []
        # number of papers to process in total
        total_papers = len(pmids) if pmids else len(pairs)

        # build selected list filtered by pmid_file (if provided), and load paper data
        papers: List[Dict[str, Any]] = []
        for pmid, pmid_dir in pairs:
            if pmids is not None and pmid not in pmids:
                continue
            
            # load paper data from meta JSON and content JSON (if exists)
            # paper is like {
            #     'pmid': '12345678',
            #     'meta': {...}, # loaded from 12345678_meta.json
            #     'content': {...} # loaded from 12345678_content.json
            # }
            paper = self._load_paper_data(pmid_dir, pmid)
            # meta and content are both empty, skip (no useful info to merge)
            if not paper['meta']:
                meta_missing.append(pmid)
                continue
            elif not paper['content']:
                content_missing.append(pmid)
            
            # write sidecar {PMID}.json
            try:
                with open(os.path.join(pmid_dir, f"{pmid}.json"), 'w') as fh:
                    json.dump(paper, fh, ensure_ascii=False, indent=2)
            except Exception:
                # best-effort: do not fail merge for a single write error
                pass
            papers.append(paper)

        # resolve output path (if directory or no ext -> auto name)
        # output_json can be a file path or a directory
        out = output_json
        # directory existed or non-existent path without extension (treat as directory)
        # /data2/pyPaperFlow/test/full_paper_test or /data2/pyPaperFlow/test/Not-Exist-Folder/
        if os.path.isdir(output_json) or not os.path.splitext(output_json)[1]:
            # paper_dir is not pubmed/ containing
            source_name = os.path.basename(os.path.normpath(paper_dir)) 
            suffix = '.jsonl' if jsonl else '.json'   
            out = os.path.join(output_json, f"{source_name}_{_timestamp()}{suffix}")

        # write merged output
        # note that papers are list of dict
        try:
            # ensure output directory exists, default is current directory
            os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
            if jsonl:
                # list of dict -> JSONL with one JSON object per line, one paper per line
                with open(out, 'w') as fh:
                    for p in papers:
                        fh.write(json.dumps(p, ensure_ascii=False) + '\n')
            else:
                # json
                # here we create pmid: i in papers 
                with open(out, 'w') as fh:
                    papers_json_dict = {p.get('pmid'): p for p in papers}
                    json.dump(papers_json_dict, fh, ensure_ascii=False, indent=2)
        except Exception as e:
            raise RuntimeError(f"Failed to write merged JSON: {e}")

        # write the stats
        # actually we care more about the content missing than meta missing, cause the later is less likely to happen 
        # ⚠️ and we can use the pmids missing content to fetch the article again by DOI module 
        stats = {
            'total': total_papers,
            'meta_missing': meta_missing,
            'content_missing': content_missing
        }
        return stats 

    def export_md_from_merged_json(
        self,
        merged_json: str,
        output_md: str,
        yaml_cfg: Optional[str] = None,
        pmid_file: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Description
        -----------
        Export a single Markdown file from merged JSON/JSONL.

        Args
        ----
        merged_json: Path to the merged JSON or JSONL file produced by merge_json_from_directory
        output_md: Path to the output Markdown file
        yaml_cfg: Optional path to a YAML config file specifying which metadata fields and content sections to include in the Markdown output. If not provided, defaults to including 'identity.title' and 'identity.pmid' for metadata, and 'abstract' for content.

        `yaml_cfg` may specify:
            metadata_fields: ["identity.title", "identity.pmid"]
            content_sections: ["abstract", "introduction"]
        """
        papers = self._load_merged(merged_json)

        if pmid_file:
            filter_pmids = set(self._read_pmids_from_file(pmid_file))
            papers = [p for p in papers if str(p.get('pmid', '')) in filter_pmids]

        cfg = {}
        if yaml_cfg:
            with open(yaml_cfg, 'r', encoding='utf-8') as yf:
                cfg = yaml.safe_load(yf) or {}

        metadata_fields = cfg.get('metadata_fields', ['identity.title', 'identity.pmid'])
        content_sections = cfg.get('content_sections', ['abstract'])

        # write markdown
        os.makedirs(os.path.dirname(output_md) or '.', exist_ok=True)
        with open(output_md, 'w', encoding='utf-8') as f:
            f.write('# Merged Literature Corpus\n\n')
            f.write(f'Total papers: {len(papers)}\n\n')
            for i, p in enumerate(papers, 1):
                pmid = str(p.get('pmid', 'N/A'))
                title = (p.get('identity') or {}).get('title', 'N/A')
                f.write(f'## PMID {pmid} - {title}\n\n')

                # metadata
                for mf in metadata_fields:
                    val = self._get_by_path(p, mf)
                    if val is None:
                        continue
                    if isinstance(val, list):
                        val = ', '.join(str(x) for x in val)
                    f.write(f'**{mf}:** {val}  \n')
                f.write('\n')

                # content sections
                text_sections = p.get('text_sections') or []
                for sec in content_sections:
                    f.write(f'### {sec.title()}\n\n')
                    # find matching section by normalized title
                    found = None
                    for s in text_sections:
                        if str(s.get('title', '')).strip().lower() == sec.strip().lower():
                            found = s.get('content')
                            break
                    if not found:
                        # fallback to content.abstract
                        if sec == 'abstract':
                            found = (p.get('content') or {}).get('abstract')
                    f.write((found or 'None') + '\n\n')

                if i < len(papers):
                    f.write('\n<!-- PAPER_BREAK -->\n\n---\n\n')

        return {'total': len(papers), 'output': output_md}

    # ---- helpers ----
    def _load_merged(self, path: str) -> List[Dict[str, Any]]:
        if not os.path.exists(path):
            return []
        try:
            with open(path, 'r', encoding='utf-8') as fh:
                data = json.load(fh)
                if isinstance(data, list):
                    return data
        except Exception:
            # try JSONL
            out = []
            with open(path, 'r', encoding='utf-8') as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        out.append(json.loads(line))
                    except Exception:
                        continue
            return out
        return []

    def _get_by_path(self, data: Dict[str, Any], path: str) -> Any:
        parts = path.split('.')
        cur = data
        for p in parts:
            if not isinstance(cur, dict):
                return None
            cur = cur.get(p)
            if cur is None:
                return None
        return cur

    def _read_pmids_from_file(self, file_path: str) -> List[str]:
        """
        Description
        -----------
        Read PMIDs from a file. 
        
        Args
        ----
        file_path: Path to a text or csv file containing PMIDs. 
        
        Returns
        -------
        A list of PMIDs as strings. If the file does not exist, returns an empty list.
        """
        
        pmids: List[str] = []
        if not os.path.exists(file_path):
            return pmids
        # for csv file
        if file_path.endswith('.csv'):
            with open(file_path, 'r', encoding='utf-8') as f:
                for row in csv.reader(f):
                    if row:
                        pmids.append(str(row[0]).strip())
            return pmids
        # for tsv or plain text file
        with open(file_path, 'r', encoding='utf-8') as f:
            for line in f:
                v = line.strip()
                if v:
                    pmids.append(v)
        return pmids

    def _iter_pmid_directories(self, paper_dir: str) -> List[tuple[str, str]]:
        """
        Description
        -----------
        Yield (pmid, pmid_dir) pairs for two common layouts.
        
        Args
        ----
        paper_dir: Root directory to scan for PMIDs. Only accepts paper_dir(without pubmed)/pubmed/<year>/<pmid>/
        
        Returns
        -------
        A list of (pmid, pmid_dir) tuples for discovered PMIDs, like ('12345678', '/path/to/paper_dir/pubmed/2023/12345678'). If no PMIDs are found, returns an empty list.
        """
        
        pairs: List[tuple[str, str]] = []
        if not os.path.isdir(paper_dir):
            return pairs
        
        # paper_dir(without pubmed)/pubmed/<year>/<pmid>/
        # like /data2/pyPaperFlow/test/full_paper_test'
        if 'pubmed' in os.listdir(paper_dir):
            paper_dir = os.path.join(paper_dir, 'pubmed')

        # year subfolders within pubmed/
        entries = sorted(os.listdir(paper_dir))
        # year subfolders
        for e in entries:
            year_dir = os.path.join(paper_dir, e)
            if not os.path.isdir(year_dir):
                continue
            for pm in sorted(os.listdir(year_dir)):
                pm_dir = os.path.join(year_dir, pm)
                # pm is PMID
                if os.path.isdir(pm_dir) and pm.isdigit():
                    pairs.append((pm, pm_dir))
        return pairs

    def _resolve_paper_files(self, pmid_dir: str, pmid: str) -> Dict[str, Optional[str]]:
        """
        Description
        -----------
        Resolve file paths for a given PMID directory.
        
        Args
        ----
        pmid_dir: Directory containing the paper's files (e.g. pubmed/2023/12345678/)
        pmid: The PMID of the paper (e.g. '12345678')
        ('12345678', '/path/to/paper_dir/pubmed/2023/12345678')
        
        Returns
        -------
        A dict with keys 'meta_json', 'content_json', 'content_md' and values as the resolved file paths or None if not found.
        
        """
        
        # simple detection of common file names
        candidates = {
            'meta_json': f'{pmid}_meta.json',
            'content_json': f'{pmid}_content.json',
            'content_md': f'{pmid}_content.md',
        }
        resolved = {}
        for k, name in candidates.items():
            # in most cases, the meta JSON exits but content JSON/MD may be missing, so we do not fail if they are not found. We just set them to None in the resolved dict.
            resolved[k] = None
            # ('12345678', '/path/to/paper_dir/pubmed/2023/12345678') / '12345678_meta.json'
            p = os.path.join(pmid_dir, name)
            if os.path.exists(p):
                # {'meta_json': '/path/to/paper_dir/pubmed/2023/12345678/12345678_meta.json', ...}
                resolved[k] = p
        return resolved

    def _load_paper_data(self, pmid_dir: str, pmid: str) -> Optional[Dict[str, Any]]:
        """
        Description
        -----------
        Load paper data from a PMID directory.
        
        Args
        ----
        pmid_dir: Directory containing the paper's files (e.g. pubmed/2023/12345678/)
        pmid: The PMID of the paper (e.g. '12345678')
        ('12345678', '/path/to/paper_dir/pubmed/2023/12345678')
        
        Returns
        -------
        A dict containing the paper's data, with keys 'meta' and 'content'. 
        'meta' is loaded from the meta JSON file if it exists, otherwise an empty dict. 
        'content' is loaded from the content JSON file if it exists, otherwise an empty dict. 
        """
        
        # paths is a dict like {'meta_json': '/path/to/paper_dir/pubmed/2023/12345678/12345678_meta.json', 'content_json': None, 'content_md': None}
        paths = self._resolve_paper_files(pmid_dir, pmid)
        
        # meta and meta_json 
        meta = paths.get('meta_json', None)
        if meta:
            try:
                with open(meta, 'r') as fh:
                    meta_json = json.load(fh)
            except Exception:
                meta_json = {}
        else:
            meta_json = {}

        # content and content_json
        content = paths.get('content_json', None)
        if content:
            try:
                with open(content, 'r') as fh:
                    content_json = json.load(fh)
            except Exception:
                content_json = {}
        else:
            content_json = {}

        # Compose a small canonical dict
        paper = {
            'pmid': pmid,
            'meta': meta_json,
            'content': content_json
        }
        return paper
