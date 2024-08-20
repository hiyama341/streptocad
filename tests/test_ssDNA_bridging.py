from pydna.dseqrecord import Dseqrecord
import pandas as pd
from teemi.design.fetch_sequences import read_genbank_files

# Source functions
from src.teemi_functions.cloning.ssDNA_bridging import assemble_plasmids_by_ssDNA_bridging#, checking_primer_design #extract_sgRNA_midpoint_positions
from src.teemi_functions.sgRNA import make_antismash_df_to_dseq
from src.teemi_functions.sgRNA import make_ssDNA_oligos
from src.teemi_functions.cloning.ssDNA_bridging import assemble_plasmids_by_ssDNA_bridging

## INITIALIZE A Relevant ssDNA assembly
# Input 1:
data_frame_1 = pd.read_csv('data/CRISPy_results/SCO5087_BEST_STOP.csv')

# Input 2: Plasmid pCRISPR–Cas9_plasmid_addgene.gbk works but any that has the relevant sites with the NcoI can be used
input_plasmid = read_genbank_files('data/plasmids/pCRISPR–Cas9_plasmid_addgene.gbk')[0]
input_plasmid = Dseqrecord(input_plasmid, circular=True)

# Input 3: THis should be the standard but the user should be able to change them
overhang_start = 'CGGTTGGTAGGATCGACGGC'
overhang_end = 'GTTTTAGAGCTAGAAATAGC'

# Calculations:
sgRNA_data = make_antismash_df_to_dseq(data_frame_1)

# Make ssDNA oligos
data_frame_1.rename(columns={"Sequence": "sgrna", 'ORF':'locus_tag'}, inplace=True)
ssDNA_list = make_ssDNA_oligos(data_frame_1, upstream_ovh=overhang_start, downstream_ovh=overhang_end)

# Digest backbone plasmid
from Bio.Restriction import NcoI
linearized_input = input_plasmid.cut(NcoI)[0]

# Assemble plasmids
assembled_vectors = assemble_plasmids_by_ssDNA_bridging(ssDNA_list, linearized_input)

# sgRNA midpoint
#sgRNA_midpoints = extract_sgRNA_midpoint_positions(assembled_vectors[0])

#Checking primers
#checking_primer_df = checking_primer_design(assembled_vectors[0].seq ,sgRNA_midpoints, offset= 250,target_tm= 55  )

def test_assemble_plasmids_by_ssDNA_bridging():
    assert len(assembled_vectors) == 5
    assert len(assembled_vectors[0]) == 11279
    assert str(assembled_vectors[0].seq[0:20]) == 'GTTTTAGAGCTAGAAATAGC'

    # check correct features
    for feature in assembled_vectors[0].features:
        if feature.type == 'sgRNA':
            print(feature.type, feature)

            # Check if 'name' is a qualifier and if it has the specific value you're looking for
            if feature.type == 'sgRNA':
                sgRNA_label = feature.qualifiers.get('label', [None])
                assert sgRNA_label is not None, "sgRNA does not have a label qualifier"
                assert sgRNA_label == "sgRNA_SCO5087", f"sgRNA label {sgRNA_label} does not match expected 'sgRNA_SCO5087'"

                # Check location
                assert feature.location.start == 11259, f"sgRNA start {feature.location.start} does not match expected start"
                assert feature.location.end == 11279, f"sgRNA end {feature.location.end} does not match expected start"

#TODO make this function
def test_make_ssDNA_oligos():
    pass

# def test_extract_sgRNA_midpoint_positions():
#     # should be in between the start and end positions of sgRNA
#     assert sgRNA_midpoints == [11259+10]
#     assert sgRNA_midpoints == [11279-10]


# def test_checking_primer_design_output_type():

#     # Check that the output is a DataFrame
#     assert isinstance(checking_primer_df, pd.DataFrame), "The output should be a pandas DataFrame"

#     # Define the expected DataFrame
#     expected_df = pd.DataFrame({
#         'Name': ['forward_checking_primer', 'reverse_checking_primer'],
#         'Sequence': ['GTTAGCTCACTCATTAGG', 'AGAATCGCAGCAACTT'],
#         'Concentration': ['25nm', '25nm'],
#         'Purification': ['STD', 'STD']
#     })

#     # Check that the output DataFrame matches the expected DataFrame
#     pd.testing.assert_frame_equal(checking_primer_df, expected_df, check_dtype=False)