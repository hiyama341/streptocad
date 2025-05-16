from streptocad.crispr.guideRNAcas3_9_12 import (
    SgRNAargs,
    extract_sgRNAs,
)

import pytest
import pandas as pd
from pydna.dseqrecord import Dseqrecord
from Bio import SeqIO


@pytest.fixture
def coelicolor_genbank_record():
    filepath = "tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb"
    record = Dseqrecord(SeqIO.read(filepath, "genbank"))
    return record


@pytest.fixture
def sgrna_args_cas9(coelicolor_genbank_record):
    return SgRNAargs(
        dseqrecord=coelicolor_genbank_record,
        locus_tag=["SCO5087"],
        cas_type="cas9",
        gc_upper=0.9999,
        gc_lower=0.0001,
        off_target_seed=13,
        off_target_upper=0,
        step=["find", "filter"],
    )


@pytest.fixture
def expected_sgrna_df():
    # load the “CRISPYweb…” CSV
    filepath = "tests/test_files/CRISPYweb_S_coelicolor_ SCO5087_0_mismatches.csv"
    df = pd.read_csv(filepath, index_col=0)
    df.reset_index(drop=True, inplace=True)
    return df


def test_extract_sgRNAs_output(sgrna_args_cas9, expected_sgrna_df):
    # run the function and normalize columns
    out = extract_sgRNAs(sgrna_args_cas9).copy()
    out.columns = out.columns.str.lower()

    # 0) length must match
    assert len(out) == len(expected_sgrna_df), (
        f"Row count mismatch: got {len(out)}, expected {len(expected_sgrna_df)}"
    )

    # 1) basic sanity
    for col in ("sgrna", "pam", "gc", "off_target_count", "sgrna_strand"):
        assert col in out.columns, f"Missing column: {col}"

    # 2) GC bounds
    assert (out["gc"] >= sgrna_args_cas9.gc_lower).all() and (
        out["gc"] <= sgrna_args_cas9.gc_upper
    ).all(), "GC content is out of bounds"

    # 3) off‐target sorted
    assert out["off_target_count"].is_monotonic_increasing, (
        "Off‐target counts must be ascending"
    )

    # 4) build lookup from expected CSV
    # note: expected has “Sequence”, “PAM”, “Strand”
    exp = expected_sgrna_df
    seq_to_pam = dict(zip(exp["Sequence"], exp["PAM"]))
    seq_to_strand = dict(zip(exp["Sequence"], exp["Strand"]))

    # 5) for each extracted guide, confirm seq, pam, strand
    for seq, pam, strand in zip(out["sgrna"], out["pam"], out["sgrna_strand"]):
        assert seq in seq_to_pam, f"Unexpected guide sequence: {seq}"
        # PAM match
        assert pam == seq_to_pam[seq], (
            f"PAM mismatch for {seq!r}: got {pam!r}, expected {seq_to_pam[seq]!r}"
        )
        # strand match
        assert strand == seq_to_strand[seq], (
            f"Strand mismatch for {seq!r}: got {strand}, expected {seq_to_strand[seq]}"
        )

    # 6) no duplicate guides
    assert out["sgrna"].is_unique, "Duplicate sgRNAs found in output"


@pytest.fixture
def sgrna_args_cas9_sco4543(coelicolor_genbank_record):
    # ← changed locus_tag to SCO4543
    return SgRNAargs(
        dseqrecord=coelicolor_genbank_record,
        locus_tag=["SCO4543"],
        cas_type="cas9",
        gc_upper=0.9999,
        gc_lower=0.0001,
        off_target_seed=13,
        off_target_upper=0,
        step=["find", "filter"],
    )


@pytest.fixture
def expected_sgrna_df_sco4543():
    # ← loading the new CSV
    filepath = (
        "tests/test_files/CRISPYweb_S_coelicolor_SCO4543_region_6815941-6816101.csv"
    )
    df = pd.read_csv(filepath, index_col=0)
    df.reset_index(drop=True, inplace=True)
    return df


def test_extract_sgRNAs_output_sco4543(
    sgrna_args_cas9_sco4543, expected_sgrna_df_sco4543
):
    # run & normalize
    out = extract_sgRNAs(sgrna_args_cas9_sco4543).copy()
    out.columns = out.columns.str.lower()

    # 0) row‐count
    assert len(out) == len(expected_sgrna_df_sco4543), (
        f"Row count mismatch: got {len(out)}, expected {len(expected_sgrna_df_sco4543)}"
    )

    # 1) required columns
    for col in ("sgrna", "pam", "gc", "off_target_count", "sgrna_strand"):
        assert col in out.columns, f"Missing column: {col}"

    # 2) GC‐bounds
    assert (out["gc"] >= sgrna_args_cas9_sco4543.gc_lower).all() and (
        out["gc"] <= sgrna_args_cas9_sco4543.gc_upper
    ).all(), "GC content is out of bounds"

    # 3) off‐target sorting
    assert out["off_target_count"].is_monotonic_increasing, (
        "Off‐target counts must be ascending"
    )

    # 4) build lookups from expected
    exp = expected_sgrna_df_sco4543
    seq2pam = dict(zip(exp["Sequence"], exp["PAM"]))
    seq2strand = dict(zip(exp["Sequence"], exp["Strand"]))

    # 5) validate each guide
    for seq, pam, strand in zip(out["sgrna"], out["pam"], out["sgrna_strand"]):
        assert seq in seq2pam, f"Unexpected guide sequence: {seq}"
        assert pam == seq2pam[seq], (
            f"PAM mismatch for {seq}: got {pam}, expected {seq2pam[seq]}"
        )
        assert strand == seq2strand[seq], (
            f"Strand mismatch for {seq}: got {strand}, expected {seq2strand[seq]}"
        )

    # 6) uniqueness
    assert out["sgrna"].is_unique, "Duplicate sgRNAs found in output"
