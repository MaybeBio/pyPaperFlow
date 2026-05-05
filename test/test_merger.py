import json
import os
import tempfile
import unittest
from pathlib import Path

from pyPaperFlow.pubmed.pubmed_merger import (
    MergeConfig,
    MergeMode,
    OutputFormat,
    OutputProfile,
    PaperMerger,
)


class TestPaperMerger(unittest.TestCase):
    def _build_fixture(self, root: str) -> str:
        base = Path(root) / "papers" / "2024" / "1001"
        base.mkdir(parents=True, exist_ok=True)

        meta = {
            "identity": {
                "pmid": "1001",
                "doi": "10.1/demo",
                "title": "Demo Paper",
            },
            "content": {
                "abstract": "A short abstract.",
                "keywords": ["demo", "ai"],
                "mesh_terms": ["Protein"],
            },
            "source": {
                "journal_title": "Demo Journal",
                "pub_date": "2024-01-01",
            },
            "contributors": {
                "medline": {
                    "full_names": ["Alice Smith", "Bob Lee"],
                }
            },
            "links": {
                "refs": ["2001", "2002"],
                "review": ["2003"],
            },
            "metadata": {
                "fetched_at": "2026-04-25 00:00:00",
            },
        }

        with open(base / "1001.json", "w", encoding="utf-8") as f:
            json.dump(meta, f)

        with open(base / "1001_parsed.md", "w", encoding="utf-8") as f:
            f.write("# Full text\n\nThis is full text.")

        parsed = {
            "title": "Demo Parsed",
            "body": [
                {
                    "title": "Abstract",
                    "content": ["A short abstract from parsed json."],
                    "subsections": [],
                },
                {
                    "title": "Introduction",
                    "content": ["Intro body."],
                    "subsections": [],
                },
                {
                    "title": "References",
                    "content": ["[1] Demo Ref"],
                    "subsections": [],
                },
            ],
        }

        with open(base / "1001_parsed.json", "w", encoding="utf-8") as f:
            json.dump(parsed, f)

        return str(Path(root) / "papers")

    def test_merge_metadata_only_markdown(self):
        with tempfile.TemporaryDirectory() as tmp:
            paper_dir = self._build_fixture(tmp)
            output = Path(tmp) / "merged.md"

            merger = PaperMerger(MergeConfig(mode=MergeMode.METADATA_ONLY, output_format=OutputFormat.MARKDOWN))
            stats = merger.merge_from_directory(paper_dir, str(output))

            self.assertEqual(stats["successful"], 1)
            content = output.read_text(encoding="utf-8")
            self.assertIn("Demo Paper", content)
            self.assertNotIn("### Full Text", content)

    def test_merge_with_pmid_file_jsonl(self):
        with tempfile.TemporaryDirectory() as tmp:
            paper_dir = self._build_fixture(tmp)
            output = Path(tmp) / "merged.jsonl"
            pmids = Path(tmp) / "pmids.txt"
            pmids.write_text("1001\n9999\n", encoding="utf-8")

            merger = PaperMerger(
                MergeConfig(
                    mode=MergeMode.METADATA_CONTENT,
                    output_format=OutputFormat.JSONL,
                    output_profile=OutputProfile.ANALYSIS,
                )
            )
            stats = merger.merge_from_file_list(paper_dir, str(pmids), str(output))

            self.assertEqual(stats["successful"], 1)
            self.assertGreaterEqual(stats["skipped"], 1)

            lines = [line for line in output.read_text(encoding="utf-8").splitlines() if line.strip()]
            self.assertEqual(len(lines), 1)
            row = json.loads(lines[0])
            self.assertEqual(row["pmid"], "1001")
            self.assertIsNotNone(row["full_text"])
            self.assertIn("contributors", row)
            self.assertIn("links", row)
            self.assertIn("metadata_selected", row)
            self.assertIn("identity", row["metadata_selected"])

    def test_merge_llm_profile_sections_and_no_links(self):
        with tempfile.TemporaryDirectory() as tmp:
            paper_dir = self._build_fixture(tmp)
            output = Path(tmp) / "merged_llm.jsonl"

            merger = PaperMerger(
                MergeConfig(
                    mode=MergeMode.METADATA_CONTENT,
                    output_format=OutputFormat.JSONL,
                    output_profile=OutputProfile.LLM,
                    include_sections=["abstract", "introduction"],
                    include_links_in_llm=False,
                )
            )
            stats = merger.merge_from_directory(paper_dir, str(output))
            self.assertEqual(stats["successful"], 1)

            row = json.loads(output.read_text(encoding="utf-8").strip())
            self.assertEqual(row["profile"], "llm")
            self.assertIn("Abstract", row["llm_text"])
            self.assertIn("Introduction", row["llm_text"])
            self.assertNotIn("References", row["llm_text"])
            self.assertEqual(row["links"], {})


if __name__ == "__main__":
    unittest.main()
