# import pytest
# import pandas as pd
# from pydna.dseqrecord import Dseqrecord

# # Local module imports
# from src.teemi_functions.sgRNA import (
#     make_ssDNA_oligos,
#     get_top_sgRNAs,
#     primers_to_IDT
# )


# # Sample data to be used across tests
# @pytest.fixture
# def sample_dataframe():
#     data = {'Sequence': ['ATGCGT', 'GGATCC', 'TTAACG'],
#             'ID': ['seq1', 'seq2', 'seq3']}
#     return pd.DataFrame(data)

# @pytest.fixture
# def sample_sgRNAs_dataframe():
#     data = {'sgrna': ['NNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNN', 'NNNNNNNNNNNNNNNNNNNN'],
#             'locus_tag': ['gene1', 'gene2', 'gene3'],
#             'off_target_n': [1, 2, 3]}
#     return pd.DataFrame(data)


# # Tests for make_ssDNA_oligos
# def test_make_ssDNA_oligos(sample_sgRNAs_dataframe):
#     # Normal behavior test
#     result = make_ssDNA_oligos(sample_sgRNAs_dataframe)
#     assert len(result) == 3
#     assert all(isinstance(item, Dseqrecord) for item in result)
#     # More tests for edge cases and invalid inputs can be added here

# # Tests for get_top_sgRNAs
# def test_get_top_sgRNAs(sample_sgRNAs_dataframe):
#     # Normal behavior test
#     result = get_top_sgRNAs(sample_sgRNAs_dataframe)
#     assert len(result) <= len(sample_sgRNAs_dataframe)
#     # More tests for edge cases and invalid inputs can be added here

# # Tests for primers_to_IDT
# def test_primers_to_IDT(sample_sgRNAs_dataframe):
#     # Preparing the input list of Dseqrecords
#     sgRNAs_p = make_ssDNA_oligos(sample_sgRNAs_dataframe)
#     # Normal behavior test
#     result = primers_to_IDT(sgRNAs_p)
#     assert len(result) == len(sgRNAs_p)
#     assert 'Name' in result.columns
#     assert 'Sequence' in result.columns
#     assert 'Concentration' in result.columns
#     assert 'Purification' in result.columns
#     # More tests for edge cases and invalid inputs can be added here

# # Here we would add additional tests for each function to cover more scenarios
