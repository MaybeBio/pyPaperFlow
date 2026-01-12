import os
import json
import pandas as pd
from typing import List, Dict, Optional, Union
from .fetcher import Paper

class PaperStorage:
    def __init__(self, storage_dir: str):
        """
        Initialize the PaperStorage.
        
        Args:
            storage_dir (str): The directory where the lookup table and paper data are stored.
        """
        self.storage_dir = storage_dir
        if not os.path.exists(self.storage_dir):
            os.makedirs(self.storage_dir)
            
        self.lookup_file = os.path.join(self.storage_dir, "paper_lookup.csv")
        self.papers_dir = os.path.join(self.storage_dir, "papers")
        if not os.path.exists(self.papers_dir):
            os.makedirs(self.papers_dir)
            
        self._load_lookup_table()

    def _normalize_tag(self, tag: str) -> str:
        """Normalize a tag into a stable column key.

        - Strips surrounding whitespace.
        - Case-folds English.
        - Collapses internal whitespace to single underscores.

        Notes:
        - Chinese tags are kept as-is except whitespace normalization.
        """
        if tag is None:
            return ""
        normalized = str(tag).strip()
        if not normalized:
            return ""
        normalized = normalized.casefold()
        normalized = "_".join(normalized.split())
        return normalized

    def _ensure_pmid_row(self, pmid: str, json_path: str = ""):
        """Ensure the lookup table contains a row for this PMID."""
        pmid = str(pmid).strip()
        if not pmid:
            raise ValueError("pmid must be non-empty")

        if pmid not in self.df['pmid'].astype(str).values:
            new_row = {'pmid': pmid, 'json_path': json_path}
            self.df = pd.concat([self.df, pd.DataFrame([new_row])], ignore_index=True)
            self._save_lookup_table()

    def _load_lookup_table(self):
        """Load the lookup table from CSV."""
        if os.path.exists(self.lookup_file):
            self.df = pd.read_csv(self.lookup_file, dtype={'pmid': str})
            if 'pmid' not in self.df.columns:
                 # Should not happen if initialized correctly, but safety check
                 self.df['pmid'] = []
        else:
            self.df = pd.DataFrame(columns=['pmid', 'json_path'])

    def _save_lookup_table(self):
        """Save the lookup table to CSV."""
        self.df.to_csv(self.lookup_file, index=False)

    def add_paper(self, paper: Paper, tags: Optional[Dict[str, int]] = None):
        """
        Add a paper to the storage.
        
        Args:
            paper (Paper): The paper object to add.
            tags (Dict[str, int], optional): Initial tags for the paper.
        """
        pmid = paper.identity.pmid
        if not pmid:
            print("Warning: Paper has no PMID, skipping storage.")
            return

        # Save individual paper JSON
        # We save each paper individually or we could reference the batch file.
        # For easier lookup, saving individually or by small groups in a hashed structure is better.
        # Here we save individually for simplicity as requested "locate paper then view metadata".
        paper_filename = f"{pmid}.json"
        paper_path = os.path.join(self.papers_dir, paper_filename)
        
        with open(paper_path, 'w', encoding='utf-8') as f:
            json.dump(paper.to_dict(), f, ensure_ascii=False, indent=2)
            
        # Update Lookup Table
        if pmid in self.df['pmid'].values:
            # Update existing entry
            idx = self.df.index[self.df['pmid'] == pmid][0]
            self.df.at[idx, 'json_path'] = paper_path
        else:
            # Add new entry
            new_row = {'pmid': pmid, 'json_path': paper_path}
            # Append using concat
            self.df = pd.concat([self.df, pd.DataFrame([new_row])], ignore_index=True)
            
        # Add tags if provided
        if tags:
            self.update_tags(pmid, tags)
        else:
            self._save_lookup_table()

    def update_tags(self, pmid: str, tags: Dict[str, int]):
        """
        Update tags for a specific paper.
        
        Args:
            pmid (str): The PMID of the paper.
            tags (Dict[str, int]): A dictionary of tag names and their values (0 or 1).
        """
        # Allow tagging even if the PMID hasn't been added via add_paper yet.
        self._ensure_pmid_row(pmid)

        idx = self.df.index[self.df['pmid'] == pmid][0]
        
        for tag, value in tags.items():
            norm_tag = self._normalize_tag(tag)
            if not norm_tag:
                continue
            if norm_tag not in self.df.columns:
                self.df[norm_tag] = 0  # Initialize new column with 0

            self.df.at[idx, norm_tag] = int(value)
            
        self._save_lookup_table()

    def add_tags(self, pmid: str, tags: List[str]):
        """Add (set to 1) a list of tags for a PMID."""
        self._ensure_pmid_row(pmid)
        updates: Dict[str, int] = {}
        for tag in tags or []:
            norm_tag = self._normalize_tag(tag)
            if norm_tag:
                updates[norm_tag] = 1
        if updates:
            self.update_tags(pmid, updates)

    def remove_tags(self, pmid: str, tags: List[str]):
        """Remove (set to 0) a list of tags for a PMID."""
        self._ensure_pmid_row(pmid)
        updates: Dict[str, int] = {}
        for tag in tags or []:
            norm_tag = self._normalize_tag(tag)
            if norm_tag:
                updates[norm_tag] = 0
        if updates:
            self.update_tags(pmid, updates)

    def set_tags(self, pmid: str, tags: List[str]):
        """Replace tags for a PMID (clear all existing tag columns, then set given tags to 1)."""
        self._ensure_pmid_row(pmid)
        idx = self.df.index[self.df['pmid'] == pmid][0]

        tag_cols = [c for c in self.df.columns if c not in ['pmid', 'json_path']]
        for col in tag_cols:
            self.df.at[idx, col] = 0

        for tag in tags or []:
            norm_tag = self._normalize_tag(tag)
            if not norm_tag:
                continue
            if norm_tag not in self.df.columns:
                self.df[norm_tag] = 0
            self.df.at[idx, norm_tag] = 1

        self._save_lookup_table()

    def list_tags(self, pmid: str) -> List[str]:
        """List enabled (value==1) tags for a PMID."""
        if pmid not in self.df['pmid'].astype(str).values:
            return []
        idx = self.df.index[self.df['pmid'] == pmid][0]
        row = self.df.iloc[idx]
        enabled: List[str] = []
        for col in self.df.columns:
            if col in ['pmid', 'json_path']:
                continue
            try:
                if int(row[col]) == 1:
                    enabled.append(col)
            except Exception:
                continue
        enabled.sort()
        return enabled

    def get_paper(self, pmid: str) -> Optional[Paper]:
        """
        Retrieve a paper object by PMID.
        
        Args:
            pmid (str): The PMID to retrieve.
            
        Returns:
            Paper: The paper object, or None if not found.
        """
        if pmid not in self.df['pmid'].values:
            return None
            
        idx = self.df.index[self.df['pmid'] == pmid][0]
        json_path = self.df.at[idx, 'json_path']
        
        if not os.path.exists(json_path):
            print(f"Error: JSON file for {pmid} not found at {json_path}")
            return None
            
        try:
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                # We need to reconstruct the Paper object from the dict
                # This requires a from_dict method or manual reconstruction
                # Since Paper is a dataclass, we can try to reconstruct it.
                # However, nested dataclasses need care.
                # For now, let's return the dict or implement a simple reconstructor if needed.
                # But the user asked for "Paper object".
                # Let's implement a simple from_dict in fetcher.py or here.
                # Ideally fetcher.py's Paper class should have a from_dict.
                # I'll assume I can just return the dict for now as the user might just want the data.
                # Wait, the user said "return Paper class data".
                # I'll add a helper to reconstruct.
                return self._reconstruct_paper(data)
        except Exception as e:
            print(f"Error reading paper {pmid}: {e}")
            return None

    def _reconstruct_paper(self, data: Dict) -> Paper:
        from .fetcher import PaperIdentity, PaperContent, PaperContributors, PaperSource, PaperMetadata, PaperLinks
        
        return Paper(
            identity=PaperIdentity(**data.get('identity', {})),
            content=PaperContent(**data.get('content', {})),
            contributors=PaperContributors(
                medline=data.get('contributors', {}).get('medline', {}),
                xml=data.get('contributors', {}).get('xml', [])
            ),
            source=PaperSource(**data.get('source', {})),
            metadata=PaperMetadata(**data.get('metadata', {})),
            links=PaperLinks(**data.get('links', {}))
        )

    def query_papers(self, query_tags: Dict[str, int]) -> List[str]:
        """
        Find PMIDs that match the given tags.
        
        Args:
            query_tags (Dict[str, int]): Tags to filter by (e.g., {'relevant': 1, 'reviewed': 0}).
            
        Returns:
            List[str]: List of matching PMIDs.
        """
        if self.df.empty:
            return []
            
        condition = pd.Series([True] * len(self.df))
        
        for tag, value in query_tags.items():
            if tag in self.df.columns:
                condition &= (self.df[tag] == value)
            else:
                # If tag doesn't exist, and we ask for 1, then no match.
                # If we ask for 0, then all match (assuming default is 0).
                if value == 1:
                    return []
        
        return self.df[condition]['pmid'].tolist()

    def get_feature_vector(self, pmid: str) -> Dict[str, int]:
        """
        Get the feature vector (tags) for a paper.
        """
        if pmid not in self.df['pmid'].values:
            return {}
            
        idx = self.df.index[self.df['pmid'] == pmid][0]
        row = self.df.iloc[idx]
        
        # Filter out non-tag columns
        tags = {}
        for col in self.df.columns:
            if col not in ['pmid', 'json_path']:
                tags[col] = int(row[col])
                
        return tags
