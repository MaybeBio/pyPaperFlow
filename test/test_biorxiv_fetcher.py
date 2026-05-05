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
            "message": {
                "items": [
                    {
                        "DOI": f"10.1101/2025.01.01.00000{index}",
                        "title": [f"Unrelated paper {index}"],
                        "publisher": "openRxiv",
                        "type": "posted-content",
                        "prefix": "10.64898",
                        "author": [{"given": "Example", "family": "Author"}],
                        "issued": {"date-parts": [[2026, 1, 5]]},
                    }
                    for index in range(5)
                ],
                "next-cursor": "cursor-2",
            }
        }
        second_page = {
            "message": {
                "items": [
                    {
                        "DOI": "10.1101/2026.01.06.123456",
                        "title": ["AlphaFold resolves structure questions"],
                        "abstract": "<jats:p>A structure-focused analysis of AlphaFold.</jats:p>",
                        "publisher": "openRxiv",
                        "type": "posted-content",
                        "prefix": "10.64898",
                        "author": [{"given": "Alice", "family": "Smith"}],
                        "issued": {"date-parts": [[2026, 1, 6]]},
                    }
                ],
                "next-cursor": "cursor-3",
            }
        }

        def fake_get(url, params=None, headers=None, timeout=None):
            if "api.crossref.org" in url:
                cursor = (params or {}).get("cursor")
                return FakeResponse(json_data=first_page if cursor == "*" else second_page)
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
        def fake_get(url, params=None, headers=None, timeout=None):
            if "api.crossref.org" in url:
                return FakeResponse(
                    json_data={
                        "message": {
                            "items": [
                                {
                                    "DOI": "10.1101/2026.01.05.123456",
                                    "title": ["AlphaFold improves protein modeling"],
                                    "abstract": "<jats:p>AlphaFold changes the landscape for structure prediction.</jats:p>",
                                    "publisher": "openRxiv",
                                    "type": "posted-content",
                                    "prefix": "10.64898",
                                    "author": [{"given": "Alice", "family": "Smith"}, {"given": "Bob", "family": "Lee"}],
                                    "issued": {"date-parts": [[2026, 1, 5]]},
                                    "subject": ["Bioinformatics"],
                                }
                            ],
                            "next-cursor": "",
                        }
                    }
                )
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