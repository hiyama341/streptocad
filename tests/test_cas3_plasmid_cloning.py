import pytest
import pandas as pd
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import Seq
from streptocad.cloning.cas3_plasmid_cloning import (
    generate_cas3_protospacer_primers, 
    cas3_plasmid_pcrs,
    assemble_cas3_plasmids,
)
from teemi.design.fetch_sequences import read_genbank_files

# Fixtures to set up the necessary components for the tests

@pytest.fixture
def pCRISPR_cas3_vector():
    path_to_plasmid = 'tests/test_files/pCRISPR_Cas3.gbk'
    plasmid = read_genbank_files(path_to_plasmid)[0]
    return Dseqrecord(plasmid, circular=True)
@pytest.fixture
def sample_spacer_table():
    data = {
        'strain_name': ['NC_003888.3', 'NC_003888.3', 'NC_003888.3', 'NC_003888.3'],
        'locus_tag': ['SCO5087', 'SCO5087', 'SCO5090', 'SCO5090'],
        'gene_loc': [5529801, 5529801, 5532706, 5532706],
        'gene_strand': [1, 1, 1, 1],
        'sgrna_strand': [1, -1, -1, -1],
        'sgrna_loc': [131, 163, 199, 505],
        'gc': [0.558824, 0.617647, 0.764706, 0.735294],
        'pam': ['TTC', 'TTC', 'TTC', 'TTC'],
        'sgrna': [
            'AATCCGAACAGCTCCTTCGATGCGTTCGTCCGGT',
            'GGATTGAAGCGCAGAGTCGTCATCACGGGCGTCG',
            'CGGATCTGGGCGCGCGTGGGCGGCCGGGTGAAGA',
            'ACGTTCGAGGACACGCTCCGGGTGCCGTCCGGGG'
        ],
        'sgrna_seed_sequence': [
            'AATCCGAA',
            'GGATTGAA',
            'CGGATCTG',
            'ACGTTCGA'
        ],
        'off_target_count': [0, 0, 3, 4],
    }
    return pd.DataFrame(data)

@pytest.fixture
def filtered_spacer_table_with_primers(sample_spacer_table):
    return generate_cas3_protospacer_primers(sample_spacer_table)


# Tests
def test_generate_cas3_primers(sample_spacer_table):
    df_with_primers = generate_cas3_protospacer_primers(sample_spacer_table)

    # Verify that the expected columns are present
    assert 'Fwd Primer' in df_with_primers.columns
    assert 'Rev Primer' in df_with_primers.columns

    # Check the content of the generated primers
    expected_fwd_primers = [
        'AATCCGAACAGCTCCTTCGATGCGTTCGTCCGGTGTCGCCCGGCAAAACCGG',
        'GGATTGAAGCGCAGAGTCGTCATCACGGGCGTCGGTCGCCCGGCAAAACCGG',
        'CGGATCTGGGCGCGCGTGGGCGGCCGGGTGAAGAGTCGCCCGGCAAAACCGG',
        'ACGTTCGAGGACACGCTCCGGGTGCCGTCCGGGGGTCGCCCGGCAAAACCGG'
    ]
    expected_rev_primers = [
        str(Seq('AATCCGAACAGCTCCTTCGATGCGTTCGTCCGGT').reverse_complement()) + "GTTTCAATCCACGCGCCCGT",
        str(Seq('GGATTGAAGCGCAGAGTCGTCATCACGGGCGTCG').reverse_complement()) + "GTTTCAATCCACGCGCCCGT",
        str(Seq('CGGATCTGGGCGCGCGTGGGCGGCCGGGTGAAGA').reverse_complement()) + "GTTTCAATCCACGCGCCCGT",
        str(Seq('ACGTTCGAGGACACGCTCCGGGTGCCGTCCGGGG').reverse_complement()) + "GTTTCAATCCACGCGCCCGT"
    ]
    assert df_with_primers['Fwd Primer'].tolist() == expected_fwd_primers
    assert df_with_primers['Rev Primer'].tolist() == expected_rev_primers
def test_cas3_plasmid_pcrs(pCRISPR_cas3_vector, filtered_spacer_table_with_primers):
    amplicons = cas3_plasmid_pcrs(pCRISPR_cas3_vector, filtered_spacer_table_with_primers)

    # Check that the correct number of amplicon pairs are generated
    assert len(amplicons) == len(filtered_spacer_table_with_primers)

    # Verify the content of each amplicon pair
    for amplicon_pair in amplicons:
        assert len(amplicon_pair) == 2  # Each pair should have two amplicons
        assert all(isinstance(amp, Dseqrecord) for amp in amplicon_pair)

def test_assemble_cas3_plasmids(pCRISPR_cas3_vector, filtered_spacer_table_with_primers):
    amplicons = cas3_plasmid_pcrs(pCRISPR_cas3_vector, filtered_spacer_table_with_primers)
    assembled_plasmids = assemble_cas3_plasmids(pCRISPR_cas3_vector, amplicons)

    # Check the number of assembled plasmids
    assert len(assembled_plasmids) == len(amplicons)

    # Verify that each plasmid is a circular Dseqrecord
    for plasmid in assembled_plasmids:
        assert isinstance(plasmid, Dseqrecord)
        assert plasmid.circular
