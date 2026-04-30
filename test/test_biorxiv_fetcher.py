from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from pyPaperFlow.biorxiv_fetcher import BioRxivFetcher
from pyPaperFlow.source_utils import basic_boolean_text_match, safe_filename


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


BIORXIV_PAYLOAD = {
    "messages": [{"status": "ok"}],
    "collection": [
        {
            "title": "AlphaFold improves protein modeling",
            "abstract": "AlphaFold changes the landscape for structure prediction.",
            "authors": [{"name": "Alice Smith"}, {"name": "Bob Lee"}],
            "date": "2026-01-05",
            "doi": "10.1101/2026.01.05.123456",
            "version": "1",
            "category": ["Bioinformatics"],
            "server": "biorxiv",
        }
    ],
}


class TestBioRxivFetcher(unittest.TestCase):
    def test_boolean_text_match_query(self):
        text = BIORXIV_PAYLOAD["collection"][0]["title"] + " " + BIORXIV_PAYLOAD["collection"][0]["abstract"]
        self.assertTrue(basic_boolean_text_match(text, "AlphaFold AND structure"))
        self.assertFalse(basic_boolean_text_match(text, "AlphaFold AND microscopy"))

    def test_search_paginates_past_first_page(self):
        first_page = {
            "messages": [{"status": "ok", "total": "31"}],
            "collection": [
                {
                    "title": f"Unrelated paper {index}",
                    "abstract": "This page does not mention the query terms.",
                    "authors": [{"name": "Example Author"}],
                    "date": "2026-01-05",
                    "doi": f"10.1101/2026.01.05.{index:06d}",
                    "version": "1",
                    "category": ["Bioinformatics"],
                    "server": "biorxiv",
                }
                for index in range(30)
            ],
        }
        second_page = {
            "messages": [{"status": "ok", "total": "31"}],
            "collection": [
                {
                    "title": "AlphaFold resolves structure questions",
                    "abstract": "A structure-focused analysis of AlphaFold.",
                    "authors": [{"name": "Alice Smith"}],
                    "date": "2026-01-06",
                    "doi": "10.1101/2026.01.06.123456",
                    "version": "1",
                    "category": ["Structural Biology"],
                    "server": "biorxiv",
                }
            ],
        }

        def fake_get(url, headers=None, timeout=None):
            if "api.biorxiv.org" in url:
                cursor = int(url.rsplit("/", 1)[-1])
                return FakeResponse(json_data=first_page if cursor == 0 else second_page)
            raise AssertionError(f"Unexpected URL: {url}")

        with patch("pyPaperFlow.biorxiv_fetcher.requests.get", side_effect=fake_get):
            fetcher = BioRxivFetcher(root_dir="/tmp", window_days=31)
            records = fetcher.search(
                "AlphaFold AND structure",
                start_date="2026-01-01",
                end_date="2026-01-31",
                max_results=5,
            )

        self.assertEqual(len(records), 1)
        self.assertEqual(records[0].title, "AlphaFold resolves structure questions")

    def test_fetch_and_save_query(self):
        def fake_get(url, headers=None, timeout=None):
            if "api.biorxiv.org" in url:
                cursor = int(url.rsplit("/", 1)[-1])
                if cursor == 0:
                    return FakeResponse(
                        json_data={
                            **BIORXIV_PAYLOAD,
                            "messages": [{"status": "ok", "total": "1"}],
                        }
                    )
                return FakeResponse(json_data={"messages": [{"status": "ok", "total": "1"}], "collection": []})
            if url.endswith(".pdf"):
                return FakeResponse(content=b"%PDF-1.4\n% biorxiv demo pdf\n")
            if "biorxiv.org/content" in url:
                html = """
                <html>
                  <head>
                    <meta name='citation_pdf_url' content='https://www.biorxiv.org/content/10.1101/2026.01.05.123456v1.full.pdf' />
                  </head>
                </html>
                """
                return FakeResponse(text=html)
            raise AssertionError(f"Unexpected URL: {url}")

        with tempfile.TemporaryDirectory() as tmp:
            fetcher = BioRxivFetcher(root_dir=tmp, window_days=31)
            with patch("pyPaperFlow.biorxiv_fetcher.requests.get", side_effect=fake_get), patch(
                "pyPaperFlow.source_utils.requests.get", side_effect=fake_get
            ):
                records = fetcher.fetch_from_query(
                    "AlphaFold AND structure",
                    output_dir=tmp,
                    start_date="2026-01-01",
                    end_date="2026-01-31",
                    max_results=5,
                    download_pdf=True,
                )

            self.assertEqual(len(records), 1)
            record = records[0]
            self.assertEqual(record.source, "biorxiv")
            self.assertTrue(record.source_id.startswith("10.1101/2026.01.05.123456"))
            self.assertTrue(record.pdf_downloaded)

            paper_dir = Path(tmp) / "biorxiv" / "2026" / record.source_id.replace("/", "_")
            stem = safe_filename(record.source_id)
            self.assertTrue((paper_dir / f"{stem}.json").exists())
            self.assertTrue((paper_dir / f"{stem}.pdf").exists())


if __name__ == "__main__":
    unittest.main()