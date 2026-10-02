"""Tests for the experimental cas12a sgRNA module.

It carries its own copy of the sgRNA pipeline (cas9, cas12a and cas3 finders),
so each finder and each pipeline step is exercised here.
"""

from collections import Counter

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord

from streptocad.experimental.cas12a.guideRNAcas12a_ttn import (
    SgRNAargs,
    extract_sgRNAs,
    filter_guides,
    find_off_target_hits,
    find_sgrna_hits_cas9,
    find_sgrna_hits_cas12a,
    find_sgrna_hits_cas3,
    parse_genbank_record,
    revcomp,
)

GENOME = "tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb"


@pytest.fixture(scope="module")
def genome():
    return Dseqrecord(SeqIO.read(GENOME, "genbank"))


def args_for(genome, cas_type, **overrides):
    kwargs = dict(
        dseqrecord=genome,
        locus_tag=["SCO5087"],
        cas_type=cas_type,
        gc_upper=0.999,
        gc_lower=0.01,
        off_target_seed=13,
        off_target_upper=10,
        step=["find", "filter"],
    )
    kwargs.update(overrides)
    return SgRNAargs(**kwargs)


# --- SgRNAargs ---------------------------------------------------------


def test_sgrnaargs_keeps_every_setting():
    record = Dseqrecord(Seq("ATGCATGCATGCATGCATGC"), linear=True)
    record.id = "Mock_Strain"

    args = SgRNAargs(
        dseqrecord=record,
        locus_tag=["gene1", "gene2"],
        cas_type="cas12a",
        downstream=20,
        off_target_seed=10,
        gc_upper=0.8,
        gc_lower=0.2,
        pam_remove=["TTTA"],
        downstream_remove=["AAAA"],
        sgrna_remove=["GGGG"],
        off_target_upper=12,
        extension_to_promoter_region=5,
        upstream_tss=150,
        dwstream_tss=150,
        target_non_template_strand=False,
        protospacer_len=23,
    )

    assert args.dseqrecord is record
    assert args.locus_tag == ["gene1", "gene2"]
    assert args.cas_type == "cas12a"
    assert args.strain_name == "Mock_Strain"
    assert args.protospacer_len == 23
    assert args.target_non_template_strand is False
    assert args.pam_remove == ["TTTA"]
    assert args.extension_to_promoter_region == 5


def test_a_single_locus_tag_is_wrapped_in_a_list():
    record = Dseqrecord(Seq("ATGC"), linear=True)
    assert SgRNAargs(dseqrecord=record, locus_tag="only_one").locus_tag == ["only_one"]


def test_defaults():
    record = Dseqrecord(Seq("ATGC"), linear=True)
    args = SgRNAargs(dseqrecord=record, locus_tag=["g"])

    assert args.cas_type == "cas9"
    assert args.step == ["find", "filter"]
    assert args.downstream == 15
    assert args.off_target_seed == 9
    assert args.gc_upper == 1
    assert args.gc_lower == 0
    assert args.off_target_upper == 10
    assert args.pam_remove is None


def test_a_non_dseqrecord_is_rejected():
    with pytest.raises(ValueError, match="Dseqrecord"):
        SgRNAargs(dseqrecord="ATGC", locus_tag=["g"])


# --- revcomp -----------------------------------------------------------


@pytest.mark.parametrize(
    "sequence, expected",
    [
        ("ATGC", "GCAT"),
        ("AAAA", "TTTT"),
        ("", ""),
        ("N", "N"),
        ("atgc", "gcat"),
    ],
)
def test_revcomp(sequence, expected):
    assert revcomp(sequence) == expected


def test_revcomp_is_its_own_inverse():
    sequence = "ATGCATTTACGGCA"
    assert revcomp(revcomp(sequence)) == sequence


# --- parse_genbank_record ---------------------------------------------


def test_parse_genbank_record_returns_both_strands(genome):
    watson, crick = parse_genbank_record(genome)

    assert len(watson) == len(crick) == len(genome.seq)
    assert watson != crick


def test_parse_genbank_record_without_features():
    record = Dseqrecord(SeqRecord(Seq("ATGCATGCATGC"), id="tiny"))
    record.features = []

    watson, crick = parse_genbank_record(record)

    assert watson == "ATGCATGCATGC"
    assert crick == revcomp("ATGCATGCATGC")


# --- find_off_target_hits ---------------------------------------------


# Carries every PAM the finders look for: TTC (cas3), GG (cas9) and TTN (cas12a).
PAM_RICH = (
    "ATTCGGTTAACCGGTTCAACCGGTTAACCTTCGG",
    "GGTTCAACCGGTTAACCTTCGGAACCGGTTAACC",
)


@pytest.mark.parametrize("cas_type", ["cas9", "cas12a", "cas3"])
def test_find_off_target_hits_counts_seeds(cas_type):
    counter = find_off_target_hits(PAM_RICH, 6, cas_type=cas_type)

    assert isinstance(counter, Counter)
    assert counter, "no seeds were counted"
    assert all(count >= 1 for count in counter.values())


@pytest.mark.parametrize("cas_type", ["cas9", "cas12a", "cas3"])
def test_find_off_target_hits_on_a_pamless_sequence_counts_nothing(cas_type):
    assert find_off_target_hits(("AAAAAAAAAA",), 6, cas_type=cas_type) == Counter()


def test_a_longer_seed_yields_longer_keys():
    short = find_off_target_hits(PAM_RICH, 4, cas_type="cas3")
    long = find_off_target_hits(PAM_RICH, 8, cas_type="cas3")

    assert max(len(k) for k in short) == 4
    assert max(len(k) for k in long) == 8


# --- the three finders ------------------------------------------------


@pytest.fixture(scope="module")
def off_targets(genome):
    sequences = parse_genbank_record(genome)
    return find_off_target_hits(sequences, 13, cas_type="cas9")


EXPECTED_COLUMNS = {
    "locus_tag",
    "gc",
    "sgrna",
    "pam",
    "off_target_count",
}


@pytest.mark.parametrize(
    "finder, needs_revcomp_kwarg",
    [
        (find_sgrna_hits_cas9, True),
        (find_sgrna_hits_cas12a, False),
        (find_sgrna_hits_cas3, True),
    ],
)
def test_each_finder_returns_a_populated_frame(
    genome, off_targets, finder, needs_revcomp_kwarg
):
    args = (genome, genome.id, ["SCO5087"], off_targets, 13)
    frame = (
        finder(*args, revcomp=revcomp)
        if needs_revcomp_kwarg
        else finder(*args, revcomp)
    )

    assert isinstance(frame, pd.DataFrame)
    assert not frame.empty
    assert EXPECTED_COLUMNS <= set(frame.columns)
    assert set(frame["locus_tag"]) <= {"SCO5087"}
    assert ((frame["gc"] >= 0) & (frame["gc"] <= 1)).all()


def test_a_finder_on_an_unknown_locus_returns_no_rows(genome, off_targets):
    frame = find_sgrna_hits_cas9(
        genome, genome.id, ["NOT_A_GENE"], off_targets, 13, revcomp=revcomp
    )
    assert frame.empty


# --- filter_guides ----------------------------------------------------


@pytest.fixture
def hitframe():
    return pd.DataFrame(
        {
            "pam": ["TTTA", "TTTC", "TTTG"],
            "sgrna": ["AAAACCCC", "GGGGTTTT", "ACGTACGT"],
            "downstream": ["AAAA", "CCCC", "GGGG"],
            "gc": [0.1, 0.5, 0.9],
            "off_target_count": [1, 5, 50],
        }
    )


def test_filter_guides_applies_the_gc_window(genome, hitframe):
    args = args_for(genome, "cas9", gc_lower=0.4, gc_upper=0.6)
    assert list(filter_guides(args, hitframe)["gc"]) == [0.5]


def test_filter_guides_applies_the_off_target_ceiling(genome, hitframe):
    args = args_for(genome, "cas9", gc_lower=0, gc_upper=1, off_target_upper=5)
    assert list(filter_guides(args, hitframe)["off_target_count"]) == [1, 5]


@pytest.mark.parametrize(
    "field, patterns, surviving",
    [
        ("pam_remove", ["TTTA"], ["TTTC", "TTTG"]),
        ("sgrna_remove", ["AAAA"], ["TTTC", "TTTG"]),
        ("downstream_remove", ["GGGG"], ["TTTA", "TTTC"]),
    ],
)
def test_filter_guides_excludes_requested_patterns(
    genome, hitframe, field, patterns, surviving
):
    args = args_for(genome, "cas9", gc_lower=0, gc_upper=1, off_target_upper=100)
    setattr(args, field, patterns)

    assert list(filter_guides(args, hitframe)["pam"]) == surviving


def test_filter_guides_without_patterns_keeps_everything(genome, hitframe):
    args = args_for(genome, "cas9", gc_lower=0, gc_upper=1, off_target_upper=100)
    assert len(filter_guides(args, hitframe)) == 3


# --- extract_sgRNAs ---------------------------------------------------


@pytest.mark.parametrize("cas_type", ["cas9", "cas12a", "cas3"])
def test_extract_sgrnas_runs_the_whole_pipeline(genome, cas_type):
    seed = 8 if cas_type == "cas3" else 13
    frame = extract_sgRNAs(args_for(genome, cas_type, off_target_seed=seed))

    assert isinstance(frame, pd.DataFrame)
    assert not frame.empty
    assert EXPECTED_COLUMNS <= set(frame.columns)
    assert (frame["off_target_count"] <= 10).all()


def test_extract_sgrnas_find_only_skips_filtering(genome):
    unfiltered = extract_sgRNAs(args_for(genome, "cas9", step=["find"], gc_upper=0.4))
    filtered = extract_sgRNAs(
        args_for(genome, "cas9", step=["find", "filter"], gc_upper=0.4)
    )

    assert len(unfiltered) >= len(filtered)
    assert (filtered["gc"] <= 0.4).all()


def test_extract_sgrnas_is_sorted_by_off_target_count(genome):
    frame = extract_sgRNAs(args_for(genome, "cas9", step=["find"]))
    counts = list(frame["off_target_count"])
    assert counts == sorted(counts)


def test_extract_sgrnas_on_an_unknown_locus_returns_no_rows(genome):
    assert extract_sgRNAs(args_for(genome, "cas9", locus_tag=["NOT_A_GENE"])).empty


def test_the_predict_step_is_a_placeholder(genome, capsys):
    extract_sgRNAs(args_for(genome, "cas9", step=["find", "predict"]))
    assert "implement this part later" in capsys.readouterr().out
