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
import re
from collections import Counter
from datetime import datetime
from typing import * 
import yaml

# typical paper sections that we want to export 
SECTION_CANONICAL_ORDER = [
    # "title"
    # "year"
    # "authors"
    'abstract',
    'introduction',
    'results',
    'discussion',
    'methods',
    'conclusion',
    'supplementary',
    'availability',
    'funding',
    'acknowledgements',
    'author_contributions',
    'other',
]

# lower case aliases for section titles.
# Each canonical section uses two tiers:
# - strong: exact phrases that should map immediately
# - weak: regex patterns, most are "in" matches 
# The two tiers are OR'ed together during normalization.
SECTION_TITLE_ALIASES = {
    'abstract': {
        'strong': {'abstract'},
        'weak': {r'abstract', r'summary'},
    },
    'introduction': {
        'strong': {'introduction', 'intro', 'background'},
        'weak': {r'introduction', r'intro', r'background'},
    },
    'results': {
        'strong': {'results', 'result', 'findings', 'finding'},
        'weak': {r'result', r'finding'},
    },
    'discussion': {
        'strong': {'discussion', 'discussions'},
        'weak': {r'discussion'},
    },
    'methods': {
        'strong': {
            'methods',
            'method',
            'materials and methods',
            'material and methods',
            'materials & methods',
            'methodology',
        },
        'weak': {
            r'method',
            r'material',
            r'methodology',
            r'materials?\s*(?:and|&)\s*methods?',
        },
    },
    'conclusion': {
        'strong': {'conclusion', 'conclusions', 'concluding remarks', 'summary'},
        'weak': {r'conclusion', r'remark', r'summary'},
    },
    'supplementary': {
        'strong': {
            'supplementary material',
            'supplementary information',
            'supplementary data',
            'supporting information',
            'supplementary',
        },
        'weak': {r'supplementary', r'supporting', r'supplement'},
    },
    'availability': {
        'strong': {
            'data availability',
            'software and data availability',
            'data and code availability',
            'data availability statement',
            'data and code availability.',
            'availability and implementation',
            'availability',
        },
        'weak': {
            r'data',
            r'software',
            r'code',
            r'availability'        
        },
    },
    'funding': {
        'strong': {'funding'},
        'weak': {r'funding', r'financial'},
    },
    'acknowledgements': {
        'strong': {'acknowledgements', 'acknowledgments', 'acknowledgment'},
        'weak': {r'acknowledgement', r'thank', r'expression'},
    },
    'author_contributions': {
        'strong': {
            'author contributions',
            'authors contributions',
            'author contribution',
            'authors contribution',
        },
        'weak': {
            r'author',
            r'contribution',
            r'statement',
        },
    },
}

# display names for section titles (for better readability in the Markdown output)
# only shown in the Markdown output(final output)
SECTION_DISPLAY_NAMES = {
    'abstract': 'Abstract',
    'introduction': 'Introduction',
    'results': 'Results',
    'discussion': 'Discussion',
    'methods': 'Methods',
    'conclusion': 'Conclusion',
    'supplementary': 'Supplementary Material',
    'availability': 'Data Availability',
    'funding': 'Funding',
    'acknowledgements': 'Acknowledgements',
    'author_contributions': 'Author Contributions',
    'other': 'Other',
}

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

        pmid_file: Optional path to a file containing a list of PMIDs to filter the papers by. If not provided, defaults to including all papers. Note that the PMIDs in this file must be present in the merged JSON for them to be included in the output.

        Notes
        -----
        - ⚠️ 1. export_md_from_merged_json must be called after merge_json_from_directory, we design it like this to ensure the merged JSON/JSONL file is up-to-date.
        """

        # papers now is list of dict loaded from merged_json, 
        # each dict is like {'pmid': '12345678', 'meta': {...}, 'content': {...} }
        papers = self._load_merged(merged_json)

        if pmid_file:
            filter_pmids = set(self._read_pmids_from_file(pmid_file))
            # filter papers by pmid_file
            papers = [p for p in papers if str(p.get('pmid', '')) in filter_pmids]

        cfg = {}
        if yaml_cfg:
            with open(yaml_cfg, 'r') as yf:
                cfg = yaml.safe_load(yf) or {}

        use_yaml_layout = bool(yaml_cfg)
        # call config yaml first, if not exist, use default values for metadata only
        if use_yaml_layout:
            metadata_fields = cfg.get('metadata_fields', ['identity.title', 'identity.pmid', 'identity.doi', 'content.keywords', 'content.mesh_terms', 'content.pub_types', 'content.abstract'])
            content_sections = cfg.get('content_sections', ['abstract', 'introduction', 'methods', 'discussion', 'conclusion', 'availability'])
        else:
            # No YAML means raw meta + raw content tree. Keep metadata lean to avoid
            # duplicating abstract/content sections that will be rendered below.
            metadata_fields = cfg.get('metadata_fields', ['identity.title', 'identity.pmid', 'identity.doi', 'content.keywords', 'content.mesh_terms', 'content.pub_types'])
            content_sections = []

        # write markdown
        os.makedirs(os.path.dirname(output_md) or '.', exist_ok=True)
        section_summary = Counter()
        raw_title_summary = Counter()

        def slugify(text: str) -> str:
            """
            Description
            -----------
            create a slug for markdown anchor from the heading text, 
            like 'PMID 12345678 - Title of the Paper' -> 'pmid-12345678-title-of-the-paper'
            """
            value = re.sub(r'[^a-z0-9\s-]', '', text.lower())
            value = re.sub(r'\s+', '-', value.strip())
            value = re.sub(r'-+', '-', value)
            return value or 'paper'

        def resolve_field(paper: Dict[str, Any], path: str) -> Any:
            value = self._get_by_path(paper.get('meta', {}), path)
            if value is not None:
                return value

            if path.startswith('content.'):
                tail = path.split('.', 1)[1]
                content = paper.get('content')
                if isinstance(content, dict):
                    value = content.get(tail)
                    if value is not None:
                        return value

                    if tail == 'abstract':
                        body_nodes = content.get('body')
                        if isinstance(body_nodes, list):
                            for node in body_nodes:
                                if not isinstance(node, dict):
                                    continue
                                if self._normalize_section_title(node.get('title')) != 'abstract':
                                    continue
                                paragraphs = node.get('content') or []
                                if isinstance(paragraphs, str):
                                    return paragraphs.strip() or None
                                if isinstance(paragraphs, list):
                                    joined = '\n\n'.join(str(item).strip() for item in paragraphs if str(item).strip())
                                    if joined:
                                        return joined
                                break

            return None

        index_entries: List[tuple[str, str, str]] = []
        for p in papers:
            pmid = str(self._get_by_path(p['meta'], 'identity.pmid') or p.get('pmid', 'N/A'))
            title = str(self._get_by_path(p['meta'], 'identity.title') or 'N/A')
            heading_text = f'PMID {pmid} - {title}' if title != 'N/A' else f'PMID {pmid}'
            # anchor for markdown link
            index_entries.append((pmid, title, slugify(heading_text)))

        with open(output_md, 'w') as f:
            for pmid, title, anchor in index_entries:
                display = f'PMID {pmid} - {title}' if title != 'N/A' else f'PMID {pmid}'
                f.write(f'- [{display}](#{anchor})\n')

            if index_entries:
                f.write('\n---\n\n')

            for i, p in enumerate(papers, 1):
                # i for paper index, p for paper dict
                pmid = str(self._get_by_path(p['meta'], 'identity.pmid') or p.get('pmid', 'N/A'))
                title = str(self._get_by_path(p['meta'], 'identity.title') or 'N/A')
                heading_text = f'PMID {pmid} - {title}' if title != 'N/A' else f'PMID {pmid}'
                anchor = slugify(heading_text)
                content = p.get('content') if isinstance(p.get('content'), dict) else {}
                body_nodes = content.get('body') if isinstance(content, dict) else []
                metadata_fields_to_write = list(metadata_fields)
                if not use_yaml_layout and not body_nodes:
                    metadata_fields_to_write.append('content.abstract')

                f.write(f'<a id="{anchor}"></a>\n\n')
                f.write(f'# PMID {pmid} - {title}\n\n')

                # flat metadata block
                for mf in metadata_fields_to_write:
                    label = mf.split('.')[-1].replace('_', ' ').title()  # default label is the last part
                    val = resolve_field(p, mf)
                    if val is None:
                        continue
                    if isinstance(val, list):
                        val = ', '.join(str(x) for x in val)
                    f.write(f'## {label}\n\n{val}\n\n')

                if use_yaml_layout:
                    # YAML provided: canonicalize only first-level body nodes.
                    section_records = self._extract_section_records(p)
                    for record in section_records:
                        section_summary[record['canonical_type']] += 1
                        raw_title_summary[record['raw_title'] or 'N/A'] += 1
                    for record in section_records:
                        # skip canonical abstract (handled via metadata/resolve_field)
                        if record.get('canonical_type') == 'abstract':
                            continue
                        # only render sections requested in the YAML config
                        if record.get('canonical_type') not in content_sections:
                            continue
                        self._render_section_records(f, [record], base_level=2)
                else:
                    # No YAML provided: render raw meta + raw content tree with no canonical mapping.
                    section_records = self._extract_section_records(p)
                    for record in section_records:
                        section_summary[record['canonical_type']] += 1
                        raw_title_summary[record['raw_title'] or 'N/A'] += 1
                    self._render_raw_section_nodes(f, body_nodes, level=2)

                if i < len(papers):
                    f.write('\n<!-- PAPER_BREAK -->\n\n---\n\n')

        return {
            'total': len(papers),
            'output': output_md,
            'section_summary': dict(section_summary),
            'raw_title_summary': dict(raw_title_summary),
        }

    # ---- helpers ----
    def _load_merged(self, path: str) -> List[Dict[str, Any]]:
        """
        Description
        -----------
        Load merged JSON/JSONL file.

        Args
        ----
        path: Path to the merged JSON or JSONL file produced by merge_json_from_directory.

        Returns
        -------
        A list of paper dicts loaded from the merged file. If the file does not exist or cannot be parsed, returns an empty list.
        Like [ {'pmid': '12345678','meta': {...},'content': {...}},... ]
        """

        if not os.path.exists(path):
            return []
        try:
            if path.endswith('.json'):
                with open(path, 'r') as fh:
                    # data now is a dict of {pmid: paper}
                    data = json.load(fh)
                    if isinstance(data, dict):
                        return list(data.values())
            elif path.endswith('.jsonl'):
                out = []
                with open(path, 'r') as fh:
                    for line in fh:
                        line = line.strip()
                        if not line:
                            continue
                        try:
                            out.append(json.loads(line))
                        except Exception:
                            continue
                return out
            else:
                raise ValueError("Unsupported file format for merged_json. Only .json and .jsonl are supported.")
        except Exception:
            return []

    def _get_by_path(self, data: Dict[str, Any], path: str) -> Any:
        """
        Description
        -----------
        Safely get a nested value from a dict using a dot-separated path.

        Args
        ----
        data: The dict to search.
        path: The dot-separated path to the value to get.

        Returns
        -------
        The value at the specified path, or None if any part of the path is not found or if the path leads to a non-dict.
        """
        parts = path.split('.')
        cur = data
        for p in parts:
            if not isinstance(cur, dict):
                return None
            cur = cur.get(p)
            if cur is None:
                return None
        return cur

    def _normalize_section_title(self, title: Any) -> str:
        """ 
        Description
        -----------
        Normalize section title to a canonical form.

        Args
        ----
        title: The raw section title to normalize. (key of the node in the tree: title)
        
        Returns
        -------
        A canonical section type string, like 'abstract', 'introduction', 'methods', etc. If the title cannot be matched to any known section, returns 'other'.
        """

        canonical_type, _ = self._match_section_title_from_cursor(title, 0)
        return canonical_type or 'other'

    def _match_section_title_from_cursor(self, title: Any, start_index: int = 0) -> tuple[Optional[str], int]:
        """
        Description
        -----------
        Match a title against canonical sections starting from a given cursor.

        Returns
        -------
        A pair of (canonical_type, next_cursor). If no match is found, returns
        (None, start_index) so the caller can keep the cursor in place.
        """

        # normalize the title
        text = re.sub(r'\s+', ' ', str(title or '').strip())
        text = text.strip(' .:-—–\t\n\r')
        lower = text.lower()

        section_names = [name for name in SECTION_CANONICAL_ORDER if name != 'other']
        if start_index < 0:
            start_index = 0

        # 1) exact alias table first: cheapest and safest.
        for index in range(start_index, len(section_names)):
            canonical = section_names[index]
            aliases = SECTION_TITLE_ALIASES.get(canonical, {})
            strong_aliases = aliases.get('strong', set())
            weak_aliases = aliases.get('weak', set())
            if lower == canonical or lower in strong_aliases:
                return canonical, index + 1

            for pattern in weak_aliases:
                if re.search(pattern, lower, flags=re.IGNORECASE):
                    return canonical, index + 1

        return None, start_index

    def _display_section_title(self, canonical_type: str) -> str:
        return SECTION_DISPLAY_NAMES.get(canonical_type, canonical_type.replace('_', ' ').title())

    def _escape_md_table(self, text: Any) -> str:
        value = str(text or '')
        return value.replace('|', '\\|').replace('\n', ' ').strip()

    def _candidate_body_nodes(self, paper: Dict[str, Any]) -> List[Dict[str, Any]]:
        """
        Description
        -----------
        Find the body nodes in the paper dictionary.
        
        Args
        ----
        paper: The paper dictionary.
        
        Returns
        -------
        A list of body nodes. 
        """

        # papers now is list of dict loaded from merged_json, 
        # each dict is like {'pmid': '12345678', 'meta': {...}, 'content': {...} }
        candidates = [
            paper.get('content'),
            paper.get('meta'),
            (paper.get('meta') or {}).get('content') if isinstance(paper.get('meta'), dict) else None,
        ]

        for candidate in candidates:
            if isinstance(candidate, list) and candidate:
                return candidate
            if isinstance(candidate, dict):
                for key in ('body', 'sections', 'content'):
                    nested = candidate.get(key)
                    if isinstance(nested, list) and nested:
                        return nested

        return []

    def _extract_section_records(
        self,
        paper: Dict[str, Any],
    ) -> List[Dict[str, Any]]:
        body_nodes = self._candidate_body_nodes(paper)

        def copy_raw_tree(nodes: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
            copied: List[Dict[str, Any]] = []
            for node in nodes or []:
                if not isinstance(node, dict):
                    continue

                copied.append({
                    'title': node.get('title'),
                    'content': node.get('content'),
                    'subsections': copy_raw_tree(node.get('subsections') or []),
                })
            return copied

        records: List[Dict[str, Any]] = []
        cursor = 0
        for node in body_nodes or []:
            if not isinstance(node, dict):
                continue

            raw_title = str(node.get('title') or 'N/A').strip()
            canonical_type, next_cursor = self._match_section_title_from_cursor(raw_title, cursor)
            if canonical_type is None:
                canonical_type = 'other'
            else:
                cursor = next_cursor

            paragraphs = node.get('content') or []
            if isinstance(paragraphs, str):
                paragraphs = [paragraphs]
            paragraphs = [str(item).strip() for item in paragraphs if str(item).strip()]

            records.append({
                'raw_title': raw_title,
                'canonical_type': canonical_type,
                'display_title': self._display_section_title(canonical_type),
                'path': raw_title,
                'depth': 0,
                'paragraphs': paragraphs,
                'paragraph_count': len(paragraphs),
                'children': copy_raw_tree(node.get('subsections') or []),
            })

        return records


    def _flatten_section_records(self, records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        flattened: List[Dict[str, Any]] = []

        def recurse(items: List[Dict[str, Any]]) -> None:
            for item in items or []:
                flattened.append(item)
                recurse(item.get('children', []))

        recurse(records)
        return flattened

    def _render_raw_section_nodes(
        self,
        handle: Any,
        sections: List[Dict[str, Any]],
        level: int = 2,
    ) -> None:
        for sec in sections or []:
            if not isinstance(sec, dict):
                continue

            title = str(sec.get('title') or 'No Title').strip() or 'No Title'
            heading_level = '#' * max(2, min(level, 6))
            handle.write(f'\n{heading_level} {title}\n\n')

            content = sec.get('content') or []
            if isinstance(content, str):
                content = [content]
            for paragraph in content:
                clean_para = str(paragraph).replace('\n', ' ').strip()
                if clean_para:
                    handle.write(f'  {clean_para}\n\n')

            self._render_raw_section_nodes(handle, sec.get('subsections', []), level=level + 1)

    def _order_section_records(
        self,
        records: List[Dict[str, Any]],
        section_order: List[str],
    ) -> List[Dict[str, Any]]:
        order_index = {name: idx for idx, name in enumerate(section_order or SECTION_CANONICAL_ORDER)}

        def sort_key(item: Dict[str, Any]) -> tuple[int, int]:
            return (order_index.get(item.get('canonical_type', 'other'), len(order_index)), item.get('depth', 0))

        # Preserve intra-paper order as much as possible while keeping common sections first.
        return sorted(records, key=sort_key)

    def _aggregate_section_records(self, records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Aggregate multiple section records that share the same canonical_type into
        a single record. This preserves the first-seen ordering from `records`
        and concatenates paragraphs and children for non-'other' canonical types.
        'other' sections are kept separate to preserve their distinct titles.
        """
        aggregated: List[Dict[str, Any]] = []
        seen: Dict[str, Dict[str, Any]] = {}

        for rec in records or []:
            ctype = rec.get('canonical_type', 'other') or 'other'
            if ctype != 'other' and ctype in seen:
                target = seen[ctype]
                target_pars = target.setdefault('paragraphs', [])
                target_children = target.setdefault('children', [])
                target_pars.extend(rec.get('paragraphs', []) or [])
                target_children.extend(rec.get('children', []) or [])
                target['paragraph_count'] = target.get('paragraph_count', 0) + rec.get('paragraph_count', 0)
            else:
                # shallow copy to avoid mutating the original record list
                new_rec = {
                    'raw_title': rec.get('raw_title'),
                    'canonical_type': rec.get('canonical_type'),
                    'display_title': rec.get('display_title'),
                    'path': rec.get('path'),
                    'depth': rec.get('depth', 0),
                    'paragraphs': list(rec.get('paragraphs', []) or []),
                    'paragraph_count': rec.get('paragraph_count', 0),
                    'children': list(rec.get('children', []) or []),
                }
                aggregated.append(new_rec)
                if ctype != 'other':
                    seen[ctype] = new_rec

        return aggregated

    def _render_section_records(
        self,
        handle: Any,
        records: List[Dict[str, Any]],
        base_level: int = 2,
        prev_heading: Optional[str] = None
    ) -> None:
        for record in records or []:
            level = base_level + int(record.get('depth', 0))
            heading_level = '#' * max(2, min(level, 6))
            display_title = record.get('display_title') or self._display_section_title(record.get('canonical_type', 'other'))
            raw_title = record.get('raw_title') or ''
            if record.get('canonical_type') == 'other' and raw_title:
                heading = raw_title
            else:
                heading = display_title

            heading_line = f'{heading_level} {heading}'.strip()
            # skip writing the heading if it's identical to the previous written heading
            if heading_line != prev_heading:
                handle.write(f'{heading_line}\n\n')
                prev_heading = heading_line

            for paragraph in record.get('paragraphs', []):
                handle.write(f'{paragraph}\n\n')
            if record.get('children'):
                # children remain raw tree nodes and are appended under this top-level section.
                self._render_raw_section_nodes(handle, record['children'], level=level + 1)

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
