# tests/test_base_editing_filter.py

import pytest
import pandas as pd
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord

from streptocad.crispr.guideRNAcas3_9_12 import (
    SgRNAargs,
    extract_sgRNAs,
)
from streptocad.crispr.crispr_best import (
    identify_base_editing_sites,
    filter_sgrnas_for_base_editing,
    process_base_editing,
)

from streptocad.sequence_loading.sequence_loading import (
    process_specified_gene_sequences_from_record,
)


@pytest.fixture
def coelicolor_genbank_record():
    """Load the full S. coelicolor chromosome."""
    path = "tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb"
    rec = SeqIO.read(path, "genbank")
    return Dseqrecord(rec)


@pytest.fixture
def sgrna_args_cas9(coelicolor_genbank_record):
    """Cas9 guide‐RNA argument for locus SCO5087."""
    return SgRNAargs(
        dseqrecord=coelicolor_genbank_record,
        locus_tag=["SCO5087"],
        cas_type="cas9",
        gc_upper=0.9999,
        gc_lower=0.0001,
        off_target_seed=13,
        off_target_upper=10,
        step=["find", "filter"],
    )


def test_filter_sgrnas_for_base_editing(coelicolor_genbank_record, sgrna_args_cas9):
    # 1) Pick the genes you want to KO
    genes_to_KO = ["SCO5087"]
    genome = coelicolor_genbank_record

    # 2) Pull out just those CDS sequences
    gene_sequences = process_specified_gene_sequences_from_record(genome, genes_to_KO)
    # only keys in genes_to_KO survive
    assert set(gene_sequences) == set(genes_to_KO)
    # sanity: each sequence is non‐empty
    for seq in gene_sequences.values():
        assert len(seq) > 100, "CDS should be at least 100 nt long"

    # 3) Extract all Cas9 sgRNAs for SCO5087
    sgrna_df = extract_sgRNAs(sgrna_args_cas9).copy()
    sgrna_df.columns = sgrna_df.columns.str.lower()

    # 4) Annotate each guide with editable‐C positions
    editing_df = identify_base_editing_sites(sgrna_df)

    # 5) Keep only the ones that actually have at least one editable C
    filtered_df = filter_sgrnas_for_base_editing(editing_df)

    # ---- assertions on the filtered output ----
    # a) it’s a subset of the original
    assert len(filtered_df) <= len(editing_df)

    # b) every row has a non‐empty editable_cytosines string
    assert (filtered_df["editable_cytosines"].str.len() > 0).all()

    # c) columns from identify_base_editing_sites are present
    for col in ("editable_cytosines", "editing_context"):
        assert col in filtered_df.columns

    # d) GUIDES without any editable Cs were dropped
    zeros = editing_df["editable_cytosines"] == ""
    assert not any(filtered_df.index.isin(editing_df[zeros].index))


@pytest.fixture
def expected_only_stop_codons_df():
    """
    Load the validation CSV containing only those sgRNAs
    predicted to yield stop codons for SCO5087.
    """
    path = "tests/test_files/CRISPy-web_BEST_SCO5087_only_stop_codons_.csv"
    df = pd.read_csv(path, index_col=0)
    df.reset_index(drop=True, inplace=True)
    return df


def test_only_sgrna_and_mutations_match_validation(
    coelicolor_genbank_record,
    sgrna_args_cas9,
    expected_only_stop_codons_df,
):
    # 1) extract just the gene sequences we care about
    genes_to_KO = ["SCO5087"]
    gene_seqs = process_specified_gene_sequences_from_record(
        coelicolor_genbank_record, genes_to_KO
    )
    assert set(gene_seqs) == set(genes_to_KO)

    # 2) run the full pipeline
    sgrna_df = extract_sgRNAs(sgrna_args_cas9).copy()
    sgrna_df.columns = sgrna_df.columns.str.lower()

    editing_df = identify_base_editing_sites(sgrna_df)
    filtered_df = filter_sgrnas_for_base_editing(editing_df)
    result_df = process_base_editing(
        filtered_df,
        gene_sequences=gene_seqs,
        only_stop_codons=True,
        editing_context=False,
    ).reset_index(drop=True)

    # 3) pick only the two columns we care about
    actual = result_df[["sgrna", "mutations"]].copy()

    expected = expected_only_stop_codons_df[["Sequence", "C to T mutations"]].copy()
    expected.columns = ["sgrna", "mutations"]

    # 4) sort both by sgrna so order cannot trip us up
    actual = actual.sort_values("sgrna").reset_index(drop=True)
    expected = expected.sort_values("sgrna").reset_index(drop=True)

    # 5) the sets of guides must match exactly
    assert list(actual["sgrna"]) == list(expected["sgrna"]), (
        f"sgRNA lists differ:\n  actual:   {actual['sgrna'].tolist()}\n  expected: {expected['sgrna'].tolist()}"
    )

    # 6) compare each mutation list after splitting + stripping whitespace
    for sgrna, act_mut_str, exp_mut_str in zip(
        actual["sgrna"], actual["mutations"], expected["mutations"]
    ):
        act_list = [m.strip() for m in act_mut_str.split(",") if m.strip()]
        exp_list = [m.strip() for m in exp_mut_str.split(",") if m.strip()]

        assert act_list == exp_list, (
            f"For guide {sgrna!r}:\n"
            f"  actual mutations = {act_list}\n"
            f"  expected mutations = {exp_list}"
        )
