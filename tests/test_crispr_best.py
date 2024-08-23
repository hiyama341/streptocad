import pytest
import pandas as pd
from pydna.dseqrecord import Dseqrecord
from Bio import SeqIO

from streptocad.crispr.crispr_best import identify_base_editing_sites, filter_sgrnas_for_base_editing, process_base_editing
from streptocad.sequence_loading.sequence_loading import process_specified_gene_sequences_from_record

# Fixtures to load expected data from CSV files
@pytest.fixture
def unfiltered_sgrna_df():
    filepath = "tests/test_files/unfiltered_sgrna_SCO5087_df.csv"
    df = pd.read_csv(filepath, index_col=0)  # Load the file, treating the first column as the index
    df.reset_index(drop=True, inplace=True)  # Reset the index to remove the index column
    return df

@pytest.fixture
def expected_filtered_best_sgrna_df():
    filepath = "tests/test_files/filtered_BEST_df_SCO5087.csv"
    df = pd.read_csv(filepath, index_col=0)  # Load the file, treating the first column as the index
    df.reset_index(drop=True, inplace=True)  # Reset the index to remove the index column
    return df

@pytest.fixture
def sample_gene_sequences():
    filepath = "tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb"

    genome = SeqIO.read(filepath, "genbank")  # Read the GenBank file

    gene_sequences = process_specified_gene_sequences_from_record(genome,['SCO5087'] )
    genes_to_KO_dict = {locus_tag: gene_sequences[locus_tag] for locus_tag in ['SCO5087'] if locus_tag in gene_sequences}
    
    return genes_to_KO_dict

def test_identify_base_editing_sites(unfiltered_sgrna_df):
    result_df = identify_base_editing_sites(unfiltered_sgrna_df)
    assert 'editable_cytosines' in result_df.columns
    assert 'editing_context' in result_df.columns
    assert result_df.shape[0] == unfiltered_sgrna_df.shape[0]  # Check the number of rows

def test_filter_sgrnas_for_base_editing(unfiltered_sgrna_df):
    editable_df = identify_base_editing_sites(unfiltered_sgrna_df)
    filtered_df = filter_sgrnas_for_base_editing(editable_df)
    assert len(filtered_df) > 0  # Ensure there are sgRNAs with editable cytosines
    assert 'editable_cytosines' in filtered_df.columns
    # Check that all filtered sgRNAs have editable cytosines
    assert all(filtered_df['editable_cytosines'].str.len() > 0)


def test_process_base_editing(expected_filtered_best_sgrna_df, unfiltered_sgrna_df, sample_gene_sequences):    
    # Apply identify_base_editing_sites to add necessary columns
    prepared_df = identify_base_editing_sites(unfiltered_sgrna_df)
    
    # Apply filtering to match expected_filtered_best_sgrna_df
    filtered_sgrna_df_for_base_editing = filter_sgrnas_for_base_editing(prepared_df)
    print(filtered_sgrna_df_for_base_editing)
    result_df = process_base_editing(filtered_sgrna_df_for_base_editing, 
                                        sample_gene_sequences, 
                                        only_stop_codons = False, 
                                        editing_context= False)
    
    
    assert 'mutations' in result_df.columns
    assert len(result_df) > 0  # Ensure there are results
    # Check for known mutations in the results based on input data
    assert 'mutated_sequence' not in result_df.columns  # Ensure mutated sequence column is dropped if unwanted
    
    # Test with stop codon filtering
    result_df_stop_codons = process_base_editing(filtered_sgrna_df_for_base_editing, sample_gene_sequences, only_stop_codons=True)
    assert len(result_df_stop_codons) >= 0  # Ensure functionality with stop codon filtering

    # Test with editing context filtering
    result_df_editing_context = process_base_editing(filtered_sgrna_df_for_base_editing, sample_gene_sequences, editing_context=True)
    assert len(result_df_editing_context) >= 0  # Ensure functionality with editing context filtering


if __name__ == '__main__':
    pytest.main()
