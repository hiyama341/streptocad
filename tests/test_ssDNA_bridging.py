import pytest
import pandas as pd
from pydna.dseqrecord import Dseqrecord
from teemi.design.fetch_sequences import read_genbank_files
from Bio.Restriction import NcoI

# Source functions
from streptocad.cloning.ssDNA_bridging import assemble_plasmids_by_ssDNA_bridging, make_ssDNA_oligos


@pytest.fixture
def filtered_df():
    return pd.read_csv('tests/test_files/filtered_BEST_df.csv')


@pytest.fixture
def input_plasmid():
    plasmid = read_genbank_files('tests/test_files/pCRISPR-Cas9.gbk')[0]
    return Dseqrecord(plasmid, circular=True)


@pytest.fixture
def homology_sequences():
    up_homology = 'CGGTTGGTAGGATCGACGGC'
    dw_homology = 'GTTTTAGAGCTAGAAATAGC'
    return up_homology, dw_homology


@pytest.fixture
def list_of_ssDNAs(filtered_df, homology_sequences):
    up_homology, dw_homology = homology_sequences
    return make_ssDNA_oligos(filtered_df, upstream_ovh=up_homology, downstream_ovh=dw_homology)


@pytest.fixture
def linearized_input(input_plasmid):
    return input_plasmid.cut(NcoI)[0]


@pytest.fixture
def assembled_vectors(list_of_ssDNAs, linearized_input):
    return assemble_plasmids_by_ssDNA_bridging(list_of_ssDNAs, linearized_input)


# Test for assemble_plasmids_by_ssDNA_bridging
def test_assemble_plasmids_by_ssDNA_bridging(assembled_vectors):
    assert len(assembled_vectors) == 7
    assert len(assembled_vectors[0]) == 11279
    assert str(assembled_vectors[0].seq[0:20]) == 'GGTTTTAGAGCTAGAAATAG'

    # check correct features
    for feature in assembled_vectors[0].features:
        if feature.type == 'sgRNA':
            # Check if 'label' is a qualifier and if it has the specific value you're looking for
            sgRNA_label = feature.qualifiers.get('label', [None])
            assert sgRNA_label is not None, "sgRNA does not have a label qualifier"
            assert sgRNA_label == "sgRNA_SCO5090_loc_140", f"sgRNA label {sgRNA_label} does not match expected"

# Test for make_ssDNA_oligos
def test_make_ssDNA_oligos(filtered_df, homology_sequences):
    up_homology, dw_homology = homology_sequences
    result = make_ssDNA_oligos(filtered_df, upstream_ovh=up_homology, downstream_ovh=dw_homology)
    assert len(result) == 7
    assert all(isinstance(item, Dseqrecord) for item in result)