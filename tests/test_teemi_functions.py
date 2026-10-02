"""Tests for the vendored teemi helpers.

The two NEB helpers post to a live API, so the HTTP call is stubbed here; the
``integration`` marker covers the real endpoint elsewhere.
"""

import json

import pytest
from Bio.SeqRecord import SeqRecord

from streptocad import teemi_functions
from streptocad.teemi_functions import (
    primer_ta_neb,
    primer_tm_neb,
    read_fasta_files,
    read_genbank_files,
)

GENBANK = "tests/test_files/pOEX-PkasO.gb"


def test_read_genbank_files_returns_seqrecords():
    records = read_genbank_files(GENBANK)

    assert len(records) == 1
    assert isinstance(records[0], SeqRecord)
    assert records[0].name == "pOEX-PkasO"
    assert len(records[0].seq) == 5230


def test_read_fasta_files_returns_seqrecords(tmp_path):
    fasta = tmp_path / "two.fasta"
    fasta.write_text(">first\nATGC\n>second\nGGCC\n")

    records = read_fasta_files(str(fasta))

    assert [r.id for r in records] == ["first", "second"]
    assert str(records[0].seq) == "ATGC"


def test_read_fasta_files_on_a_genbank_file_finds_nothing(tmp_path):
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    assert read_fasta_files(str(empty)) == []


class FakeResponse:
    def __init__(self, payload):
        self.content = json.dumps(payload).encode()


@pytest.fixture
def fake_post(monkeypatch):
    """Capture the request and return a canned NEB response."""
    calls = []

    def post(url, data=None, headers=None):
        calls.append({"url": url, "payload": json.loads(data), "headers": headers})
        return FakeResponse(post.payload)

    post.payload = {"success": True, "data": [{"tm1": 64, "ta": 66}]}
    monkeypatch.setattr(teemi_functions.requests, "post", post)
    post.calls = calls
    return post


def test_primer_tm_neb_returns_the_reported_tm(fake_post):
    assert primer_tm_neb("ATGCATGCATGCATGCATGC") == 64

    call = fake_post.calls[0]
    assert call["url"] == "https://tmapi.neb.com/tm/batch"
    assert call["payload"]["seqpairs"] == [["ATGCATGCATGCATGCATGC"]]
    assert call["payload"]["conc"] == 0.5
    assert call["payload"]["prodcode"] == "q5-0"


def test_primer_ta_neb_returns_the_reported_ta(fake_post):
    assert primer_ta_neb("ATGC", "GGCC", conc=0.25, prodcode="q5-1") == 66

    payload = fake_post.calls[0]["payload"]
    assert payload["seqpairs"] == [["ATGC", "GGCC"]]
    assert payload["conc"] == 0.25
    assert payload["prodcode"] == "q5-1"


@pytest.mark.parametrize("call", [primer_tm_neb, primer_ta_neb])
def test_a_failed_neb_request_reports_the_error(fake_post, capsys, call):
    fake_post.payload = {"success": False, "error": ["bad sequence"]}

    args = ("ATGC",) if call is primer_tm_neb else ("ATGC", "GGCC")
    assert call(*args) is None

    out = capsys.readouterr().out
    assert "request failed" in out
    assert "bad sequence" in out
