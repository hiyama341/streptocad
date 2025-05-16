import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pydna.dseqrecord import Dseqrecord
import pandas as pd
from collections import Counter

# Import the relevant functions from crispr_best.py
from streptocad.crispr.guideRNA_crispri import (
    filter_crispri_guides,
    find_sgrna_hits_cas9_crispri,
    extract_sgRNAs_for_crispri,
)
from streptocad.crispr.guideRNAcas3_9_12 import SgRNAargs, revcomp


# Sample fixture for a GenBank record
@pytest.fixture
def coelicolor_genbank_record():
    filepath = "tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb"
    record = Dseqrecord(SeqIO.read(filepath, "genbank"))
    return record


# Sample fixture for SgRNAargs
@pytest.fixture
def sgrna_args_cas9_crispri(coelicolor_genbank_record):
    return SgRNAargs(
        dseqrecord=coelicolor_genbank_record,
        locus_tag=["SCO5087"],
        cas_type="cas9",
        gc_upper=0.72,
        gc_lower=0.2,
        off_target_seed=13,
        off_target_upper=0,
        step=["find", "filter"],
        extension_to_promoter_region=200,
        target_non_template_strand=True,
    )


# Fixture for expected output
@pytest.fixture
def expected_sgrna_df_crispri():
    filepath = "tests/test_files/crispri_guides_unfiltered_SCO5087.csv"
    df = pd.read_csv(filepath, index_col=0)
    df.reset_index(drop=True, inplace=True)
    return df


# Sample fixture for SgRNAargs
@pytest.fixture
def sgrna_args_cas9_crispri_CRISPy_web(coelicolor_genbank_record):
    return SgRNAargs(
        dseqrecord=coelicolor_genbank_record,
        locus_tag=["SCO5087"],
        step=["find", "filter"],
        gc_upper=0.9999,
        gc_lower=0.0001,
        off_target_seed=13,
        off_target_upper=0,
        cas_type="cas9",
        extension_to_promoter_region=200,
        target_non_template_strand=True,
        upstream_tss=70,
        dwstream_tss=76,
    )


# Fixture for expected output
@pytest.fixture
def expected_sgrna_df_crispri_from_CRISPy_web():
    filepath = "tests/test_files/CRISPy-web_CRISPRi_region_SCO5087_s_coelicolorA3_0mismatch.csv"
    df = pd.read_csv(filepath, index_col=0)
    df.reset_index(drop=True, inplace=True)
    return df


# Test the filter_crispri_guides function
def test_filter_crispri_guides(sgrna_args_cas9_crispri, expected_sgrna_df_crispri):
    # Create a sample hitframe DataFrame to filter
    hitframe = expected_sgrna_df_crispri.copy()

    # Apply the filtering function
    filtered_df = filter_crispri_guides(sgrna_args_cas9_crispri, hitframe)

    # Check that the DataFrame is not empty after filtering
    assert not filtered_df.empty, "The filtered DataFrame should not be empty."

    # Validate columns
    assert "strain_name" in filtered_df.columns, (
        "Column 'strain_name' is missing from the filtered DataFrame."
    )
    assert "gc" in filtered_df.columns, (
        "Column 'gc' is missing from the filtered DataFrame."
    )

    # Validate filtering by GC content
    assert (
        filtered_df["gc"]
        .between(sgrna_args_cas9_crispri.gc_lower, sgrna_args_cas9_crispri.gc_upper)
        .all()
    ), "All GC content should be within the specified range."


# Test the find_sgrna_hits_cas9_crispri function
def test_find_sgrna_hits_cas9_crispri(
    coelicolor_genbank_record, sgrna_args_cas9_crispri
):
    # Use a Counter to mock off-target results
    off_target_counter = Counter()
    off_target_counter["GCGCCGCGTCCAG"] = 5  # Mock off-target seed sequence count

    # Execute the function to find sgRNA hits
    sgrna_df = find_sgrna_hits_cas9_crispri(
        coelicolor_genbank_record,
        sgrna_args_cas9_crispri.strain_name,
        sgrna_args_cas9_crispri.locus_tag,
        off_target_counter,
        sgrna_args_cas9_crispri.off_target_seed,
        revcomp,
        sgrna_args_cas9_crispri.extension_to_promoter_region,
    )

    # Check that the DataFrame is not empty
    assert not sgrna_df.empty, "The resulting DataFrame should not be empty."

    # Validate expected columns
    expected_columns = [
        "strain_name",
        "locus_tag",
        "gene_loc",
        "gene_strand",
        "sgrna_strand",
        "sgrna_loc",
        "gc",
        "pam",
        "sgrna",
        "sgrna_seed_sequence",
        "off_target_count",
    ]
    assert list(sgrna_df.columns) == expected_columns, (
        "The columns do not match the expected output."
    )

    # Validate expected sgRNA sequences
    expected_sgrna = [
        "CGCGTGTTCGTCTGCGCGCG",
        "CCCGCAGGCTCGGTAAGGAG",
        "GCCCGCAGGCTCGGTAAGGA",
    ]
    missing_sgrnas = [
        sgrna for sgrna in expected_sgrna if sgrna not in sgrna_df["sgrna"].values
    ]

    assert not missing_sgrnas, (
        f"The following expected sgRNAs were not found in the results: {missing_sgrnas}"
    )


# Test the extract_sgRNAs_for_crispri function
def test_extract_sgRNAs_for_crispri(
    sgrna_args_cas9_crispri_CRISPy_web, expected_sgrna_df_crispri_from_CRISPy_web
):
    # Execute the extraction function
    out = extract_sgRNAs_for_crispri(sgrna_args_cas9_crispri_CRISPy_web)
    out.columns = out.columns.str.lower()

    # 0) row‐count
    assert len(out) == len(expected_sgrna_df_crispri_from_CRISPy_web), (
        f"Row count mismatch: got {len(out)}, expected {len(expected_sgrna_df_crispri_from_CRISPy_web)}"
    )

    # 1) required columns
    for col in ("sgrna", "pam", "gc", "off_target_count", "sgrna_strand"):
        assert col in out.columns, f"Missing column: {col}"

    # 2) GC‐bounds
    assert (out["gc"] >= sgrna_args_cas9_crispri_CRISPy_web.gc_lower).all() and (
        out["gc"] <= sgrna_args_cas9_crispri_CRISPy_web.gc_upper
    ).all(), "GC content is out of bounds"

    # 3) off‐target sorting
    assert out["off_target_count"].is_monotonic_increasing, (
        "Off‐target counts must be ascending"
    )

    # 4) build lookups from expected
    exp = expected_sgrna_df_crispri_from_CRISPy_web
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
