import pytest
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord
from pydna.amplicon import Amplicon
from pydna.primer import Primer
from teemi.build.PCR import primer_tm_neb
from src.teemi_functions.cloning.pcr_simulation import (
    perform_pcr_on_sequences, 
    make_amplicons,
)

# Fixtures for test setup

@pytest.fixture
def sample_primer_dataframe():
    """
    Fixture to create a sample dataframe for testing perform_pcr_on_sequences.
    """
    data = {
        'f_primer_sequences(5-3)': ['ATGCGTACGTAGC', 'TGCAGTGCAGTGC', 'GCTAGCTAGCTA'],
        'r_primer_sequences(5-3)': ['CGTACGTACGTA', 'TGCAGTCGATGC', 'TAGCTAGCTAGC'],
        'template': ['template_1', 'template_2', 'template_3']
    }
    return pd.DataFrame(data)

@pytest.fixture
def sample_clean_sequences():
    """
    Fixture to create a sample list of clean sequences for testing perform_pcr_on_sequences.
    """
    sequences = [
        Dseqrecord(Seq("ATGCGTACGTAGCCGTACGTACGTAGCTAGCTAGCTAGCTAGC"), id="template_1"),
        Dseqrecord(Seq("TGCAGTGCAGTGCTGCAGTCGATGCGTACGTACGTACGTACG"), id="template_2"),
        Dseqrecord(Seq("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"), id="template_3")
    ]
    return sequences

@pytest.fixture
def sample_amplicons():
    """
    Fixture to create a sample list of SeqRecords for testing make_amplicons.
    """
    return [
        SeqRecord(Seq("ATGCGTACGTAGCCGTACGTACGTAGCTAGCTAGCTAGCTAGC"), id="template_1", name="Template 1"),
        SeqRecord(Seq("TGCAGTGCAGTGCTGCAGTCGATGCGTACGTACGTACGTACG"), id="template_2", name="Template 2"),
        SeqRecord(Seq("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"), id="template_3", name="Template 3")
    ]

# Tests

def test_perform_pcr_on_sequences(sample_primer_dataframe, sample_clean_sequences):
    amplicons = perform_pcr_on_sequences(sample_primer_dataframe, sample_clean_sequences)
    
    # Check the number of amplicons generated
    assert len(amplicons) == 3

    # Verify that each amplicon is correctly named and is an instance of Amplicon
    for i, amplicon in enumerate(amplicons):
        assert isinstance(amplicon, Amplicon)
        assert amplicon.name == f"template_{i+1}_amplicon"

    # Verify the lengths of the amplicons
    expected_lengths = [44, 46, 42]  # Example lengths based on sequences
    actual_lengths = [len(amplicon) for amplicon in amplicons]
    assert actual_lengths == expected_lengths

def test_make_amplicons(sample_amplicons):
    amplicons = make_amplicons(sample_amplicons, target_tm=60, limit=10, primer_concentration=0.5, polymerase='phusion')

    # Check that the correct number of amplicons are generated
    assert len(amplicons) == 3

    # Verify that each amplicon is correctly named and is an instance of Amplicon
    for i, amplicon in enumerate(amplicons):
        assert isinstance(amplicon, Amplicon)
        assert amplicon.name == f"Template {i+1}_amplicon"

    # Verify
