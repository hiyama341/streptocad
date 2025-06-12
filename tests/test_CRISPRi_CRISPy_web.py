from streptocad.crispr.guideRNAcas3_9_12 import (
    SgRNAargs,
    extract_sgRNAs,
)

from streptocad.crispr.guideRNA_crispri import (
    extract_sgRNAs_for_crispri,
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
        extension_to_promoter_region=200,
        upstream_tss=70,
        dwstream_tss=76,
        target_non_template_strand=True,
    )


@pytest.fixture
def expected_sgrna_df():
    # load the “CRISPYweb…” CSV
    filepath = "tests/test_files/CRISPy-web_CRISPRi_region_SCO5087_s_coelicolorA3_0mismatch.csv"
    df = pd.read_csv(filepath, index_col=0).sort_values(by=["Start"])
    df.reset_index(drop=True, inplace=True)
    return df


def test_extract_sgRNAs_output(sgrna_args_cas9, expected_sgrna_df):
    # run the function and normalize columns
    out = (
        extract_sgRNAs_for_crispri(sgrna_args_cas9)
        .sort_values(by="off_target_count", ascending=True)
        .copy()
    )
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

    # --- THE IMPORTANT PART: row-by-row matching ---
    # Normalize expected column names for matching
    exp = expected_sgrna_df.rename(
        columns={"Sequence": "sgrna", "PAM": "pam", "Strand": "sgrna_strand"}
    )

    # Merge to find matching rows (inner join)
    merged = out.merge(
        exp, on=["sgrna", "pam", "sgrna_strand"], how="left", indicator=True
    )

    # Rows in out that didn't find a match in expected will have _merge == 'left_only'
    not_found = merged[merged["_merge"] == "left_only"]

    assert not not_found.empty == False, (
        f"The following rows were not found in expected DataFrame:\n{not_found[['sgrna', 'pam', 'sgrna_strand']]}"
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
    filepath = "tests/test_files/CRISPYweb_S_coelicolor_SCO4543_region_6815941-6816101.csv"  # (4959583 - 4960435)
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
