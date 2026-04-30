from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pyPaperFlow.arxiv_fetcher import ArxivFetcher
from pyPaperFlow.source_utils import safe_filename


class FakeResponse:
    def __init__(self, *, text: str = "", content: bytes = b"", json_data=None, status_code: int = 200):
        self.text = text
        self.content = content
        self._json_data = json_data
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")

    def json(self):
        if self._json_data is None:
            raise ValueError("No JSON payload")
        return self._json_data


ARXIV_FEED = """<?xml version='1.0' encoding='UTF-8'?>
<feed xmlns='http://www.w3.org/2005/Atom' xmlns:arxiv='http://arxiv.org/schemas/atom'>
  <entry>
    <id>http://arxiv.org/abs/2401.12345v1</id>
    <updated>2024-01-02T00:00:00Z</updated>
    <published>2024-01-01T00:00:00Z</published>
    <title>Deep Learning for Demo Papers</title>
    <summary>This paper studies deep learning for biology.</summary>
    <author><name>Alice Smith</name></author>
    <author><name>Bob Lee</name></author>
    <arxiv:doi>10.48550/arXiv.2401.12345</arxiv:doi>
    <arxiv:journal_ref>Demo Journal 2024</arxiv:journal_ref>
    <category term='cs.LG'/>
    <link rel='related' title='pdf' href='https://arxiv.org/pdf/2401.12345v1.pdf' type='application/pdf'/>
  </entry>
</feed>
"""


class TestArxivFetcher(unittest.TestCase):
    def test_build_query_from_free_text(self):
        fetcher = ArxivFetcher(root_dir="/tmp")
        query = fetcher.build_query("deep learning", start_date="2024-01-01", end_date="2024-12-31")

        self.assertIn("all:deep", query)
        self.assertIn("all:learning", query)
        self.assertIn("submittedDate:[202401010000 TO 202412312359]", query)

    def test_fetch_and_save_query(self):
        def fake_get(url, params=None, headers=None, timeout=None):
            if "export.arxiv.org" in url:
                return FakeResponse(text=ARXIV_FEED)
            if url.endswith(".pdf"):
                return FakeResponse(content=b"%PDF-1.4\n% arxiv demo pdf\n")
            raise AssertionError(f"Unexpected URL: {url}")

        with tempfile.TemporaryDirectory() as tmp:
            fetcher = ArxivFetcher(root_dir=tmp)
            with patch("pyPaperFlow.arxiv_fetcher.requests.get", side_effect=fake_get), patch(
                "pyPaperFlow.source_utils.requests.get", side_effect=fake_get
            ):
                records = fetcher.fetch_from_query("deep learning", output_dir=tmp, max_results=5)

            self.assertEqual(len(records), 1)
            record = records[0]
            self.assertEqual(record.source, "arxiv")
            self.assertEqual(record.source_id, "2401.12345v1")
            self.assertTrue(record.pdf_downloaded)

            paper_dir = Path(tmp) / "arxiv" / "2024" / "2401.12345v1"
            stem = safe_filename(record.source_id)
            self.assertTrue((paper_dir / f"{stem}.json").exists())
            self.assertTrue((paper_dir / f"{stem}.pdf").exists())


if __name__ == "__main__":
    unittest.main()