from pyPaperFlow.preprint.arxiv_fetcher import ArxivFetcher, _html_to_text

HTML = """<html><body>
<h1 class="ltx_title">Title</h1>
<section class="ltx_section"><h2 class="ltx_title">Introduction</h2>
<p class="ltx_p">First paragraph.</p><p class="ltx_p">Second paragraph.</p></section>
<section class="ltx_section"><h2 class="ltx_title">Results</h2><p class="ltx_p">Result text.</p></section>
</body></html>"""


def test_html_to_text_extracts_sections():
    text = _html_to_text(HTML)
    assert "## Introduction" in text
    assert "First paragraph." in text
    assert "## Results" in text


def test_fetch_full_text_returns_body():
    fetcher = ArxivFetcher(root_dir="/tmp")

    class FakeResp:
        status_code = 200
        text = HTML

    class FakeClient:
        def get(self, url, follow_redirects=True):
            assert url == "https://ar5iv.labs.arxiv.org/html/2301.00001"
            return FakeResp()

    fetcher._http_client = FakeClient()
    text = fetcher.fetch_full_text("2301.00001v2")
    assert "First paragraph." in text


def test_fetch_full_text_404_returns_empty():
    fetcher = ArxivFetcher(root_dir="/tmp")

    class FakeResp:
        status_code = 404
        text = ""

    class FakeClient:
        def get(self, url, follow_redirects=True):
            return FakeResp()

    fetcher._http_client = FakeClient()
    assert fetcher.fetch_full_text("2301.00001") == ""
