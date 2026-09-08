import types

from pyPaperFlow.preprint.biorxiv_fetcher import BioRxivFetcher, _jats_xml_to_text
from pyPaperFlow.preprint.europepmc_fetcher import EuropePMCFullText

XML = """<article><body><sec><title>Abstract</title><p>Abs text.</p></sec>
<sec><title>Results</title><p>Result one.</p><p>Result two.</p></sec></body></article>"""


def test_jats_xml_to_text():
    text = _jats_xml_to_text(XML)
    assert "## Abstract" in text
    assert "Abs text." in text
    assert "## Results" in text
    assert "Result two." in text


def test_full_text_xml_resolves_then_fetches():
    ft = EuropePMCFullText()

    def search_resp():
        r = types.SimpleNamespace()
        r.raise_for_status = lambda: None
        r.json = lambda: {"resultList": {"result": [{"pmcid": "PMC123"}]}}
        return r

    def xml_resp():
        r = types.SimpleNamespace()
        r.raise_for_status = lambda: None
        r.text = XML
        return r

    class FakeClient:
        def __init__(self):
            self.urls = []

        def get(self, url, params=None):
            self.urls.append(url)
            if "fullTextXML" in url:
                return xml_resp()
            return search_resp()

        def close(self):
            pass

    ft._client = FakeClient()
    xml = ft.full_text_xml("10.1101/2026.01.01.123456")
    assert "Results" in xml
    assert any("PMC/PMC123/fullTextXML" in u for u in ft._client.urls)


def test_full_text_xml_no_result_empty():
    ft = EuropePMCFullText()
    r = types.SimpleNamespace()
    r.raise_for_status = lambda: None
    r.json = lambda: {"resultList": {"result": []}}

    class FakeClient:
        def get(self, url, params=None):
            return r

        def close(self):
            pass

    ft._client = FakeClient()
    assert ft.full_text_xml("10.1101/2026.01.01.999999") == ""


def test_biorxiv_fetch_full_text(monkeypatch):
    class FakeFT:
        def __init__(self, request_timeout, max_retries):
            pass

        def full_text_xml(self, doi):
            assert doi == "10.1101/2026.01.01.123456"
            return XML

        def close(self):
            pass

    monkeypatch.setattr("pyPaperFlow.preprint.europepmc_fetcher.EuropePMCFullText", FakeFT)
    fetcher = BioRxivFetcher(root_dir="/tmp", platform="biorxiv")
    text = fetcher.fetch_full_text("10.1101/2026.01.01.123456")
    assert "## Results" in text
    assert "Result one." in text


def test_biorxiv_fetch_full_text_empty(monkeypatch):
    class FakeFT:
        def __init__(self, request_timeout, max_retries):
            pass

        def full_text_xml(self, doi):
            return ""

        def close(self):
            pass

    monkeypatch.setattr("pyPaperFlow.preprint.europepmc_fetcher.EuropePMCFullText", FakeFT)
    fetcher = BioRxivFetcher(root_dir="/tmp", platform="biorxiv")
    assert fetcher.fetch_full_text("10.1101/2026.01.01.999999") == ""
