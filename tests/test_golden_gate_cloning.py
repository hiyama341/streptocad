import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from teemi.build.PCR import primer_tm_neb
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
from pydna.dseqrecord import Dseqrecord
import pandas as pd
from teemi.design.fetch_sequences import read_genbank_files
from pydna.amplicon import Amplicon

from streptocad.cloning.golden_gate_cloning import (
    GoldenGateCloning,
    make_golden_gate_overhangs_forward_primers,
    make_golden_gate_overhangs_reverse_primers,
    make_amplicons,
    digest_amplicons_w_BsaI,
    create_overhang_dataframe,
)

from streptocad.utils import dataframe_to_seqrecords
from streptocad.utils import polymerase_dict


# Input
filtered_df = pd.read_csv("tests/test_files/filtered_BEST_df.csv")
pCRISPR_plasmid_MCBE = read_genbank_files(
    "tests/test_files/pCRISPR-MCBE_Csy4_kasopGFP.gb"
)[0]
vector = Dseqrecord(pCRISPR_plasmid_MCBE, circular=True)
vector.name = "pCRISPR-MCBE_Csy4"

# Calculations:
sgRNA_list = dataframe_to_seqrecords(filtered_df)
# This is the 82nt sgRNA handle + 28nt csy4
sgRNA_handle_cys4_sites = [
    Dseqrecord(
        "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTTCACTGCCGTATAGGCAGCTAAGAAA",
        name="sgRNA_handle_cys4",
    )
] * len(sgRNA_list)
cys4 = "gTTCACTGCCGTATAGGCAGCTAAGAAA"
chosen_polymerase = polymerase_dict["Phusion High-Fidelity DNA Polymerase (GC Buffer)"]
primer_concentration = 0.4
primer_number_increment = 0

golden_gate = GoldenGateCloning(
    sgRNA_list,
    sgRNA_handle_cys4_sites,
    target_tm=65,
    cys4=cys4,
    tm_function=primer_tm_neb,
    primer_concentration=primer_concentration,
    polymerase=chosen_polymerase,
)

# make amplicons
list_of_amplicons = golden_gate.simulate_pcrs()

# Making a vector
digest = digest_amplicons_w_BsaI(list_of_amplicons)

# melting temperatures
list_of_melting_temperatures = golden_gate.calculate_melting_temperature()


# Digest vector
from Bio.Restriction import NcoI, NheI

vector_NcoI_NheI = vector.cut(NcoI, NheI)[1]
rec_vec = vector_NcoI_NheI

# assemble vecor
for i in range(len(digest)):
    rec_vec += digest[i]
rec_vec = rec_vec.looped()


#######################
##### TESTING ########
#######################

# TODO make a direct test of Golden gate cloning


def test_amplicon_lengths():
    predicted_length_of_amplicons = [190, 158, 158, 158, 158, 158, 158]
    list_of_amplicons = golden_gate.simulate_pcrs()
    print("lenght of amplicons:   ", list_of_amplicons)
    length_of_amplicons = [len(amplicon) for amplicon in list_of_amplicons]
    assert predicted_length_of_amplicons == length_of_amplicons
    assert all(isinstance(item, Dseqrecord) for item in list_of_amplicons)
    assert len(list_of_amplicons) > 0


def test_melting_temperatures():
    predicted_melting_temperatures = [
        (63, 59),
        (63, 59),
        (63, 59),
        (63, 59),
        (63, 59),
        (63, 59),
        (63, 59),
    ]

    assert predicted_melting_temperatures == list_of_melting_temperatures


def test_pcr_dataframe():
    data = {
        "amplicon_name": [
            "amplicon_0_(SCO5090_140)",
            "amplicon_1_(SCO5089_167)",
            "amplicon_2_(SCO5090_174)",
            "amplicon_3_(SCO5090_208)",
            "amplicon_4_(SCO5087_210)",
            "amplicon_5_(SCO5087_522)",
            "amplicon_6_(SCO5090_640)",
        ],
        "f_primer_id": [
            "primer_0",
            "primer_2",
            "primer_4",
            "primer_6",
            "primer_8",
            "primer_10",
            "primer_12",
        ],
        "r_primer_id": [
            "primer_1",
            "primer_3",
            "primer_5",
            "primer_7",
            "primer_9",
            "primer_11",
            "primer_13",
        ],
        "template_for_amplification": [
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
        ],
        "f_primer_tm": [63, 63, 63, 63, 63, 63, 63],
        "r_primer_tm": [59, 59, 59, 59, 59, 59, 59],
        "ta": [63, 63, 63, 63, 63, 63, 63],
    }

    expected_df = pd.DataFrame(data)
    pd.testing.assert_frame_equal(
        golden_gate.make_pcr_df(), expected_df, check_dtype=False
    )


def test_primer_dataframe():
    data = {
        "Locus Tag": [
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
            "sgRNA_handle_cys4",
        ],
        "f_primer_name": [
            "primer_0",
            "primer_2",
            "primer_4",
            "primer_6",
            "primer_8",
            "primer_10",
            "primer_12",
        ],
        "r_primer_name": [
            "primer_1",
            "primer_3",
            "primer_5",
            "primer_7",
            "primer_9",
            "primer_11",
            "primer_13",
        ],
        "f_primer_sequences(5-3)": [
            "GATCGGGTCTCCCATGGTTCACTGCCGTATAGGCAGCTAAGAAACACGTACAGGTCCTGGAGCGGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC",
            "GATCGGGTCTCCGCGCGACTCGAGAGCCGGTAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC",
            "GATCGGGTCTCCCGCCCAGATCCGGAAGCGTTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC",
            "GATCGGGTCTCCGGGACGTCCAGGTCTTCACCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC",
            "GATCGGGTCTCCGAGCAGTTCCCAGAACTGCCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC",
            "GATCGGGTCTCCGTCGACCTCCCAGTCACGGCGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC",
            "GATCGGGTCTCCGACCGTCCAGGACATGACCAGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC",
        ],
        "r_primer_sequences(5-3)": [
            "GATCAGGTCTCGGCGCTTTCTTAGCTGCCTATACGGC",
            "GATCAGGTCTCGGGCGTTTCTTAGCTGCCTATACGGC",
            "GATCAGGTCTCGTCCCTTTCTTAGCTGCCTATACGGC",
            "GATCAGGTCTCGGCTCTTTCTTAGCTGCCTATACGGC",
            "GATCAGGTCTCGCGACTTTCTTAGCTGCCTATACGGC",
            "GATCAGGTCTCGGGTCTTTCTTAGCTGCCTATACGGC",
            "GATCAGGTCTCGCTAGTTTCTTAGCTGCCTATACGGC",
        ],
        "f_tm": [63, 63, 63, 63, 63, 63, 63],
        "r_tm": [59, 59, 59, 59, 59, 59, 59],
        "ta": [63, 63, 63, 63, 63, 63, 63],
    }

    primer_df = golden_gate.generate_primer_dataframe()
    expected_primer_dataframe = pd.DataFrame(data)
    pd.testing.assert_frame_equal(
        primer_df, expected_primer_dataframe, check_dtype=False
    )


def test_initialization():
    # Test the successful initialization of GoldenGateCloning instance
    assert golden_gate.sgRNAs is not None
    assert golden_gate.list_of_amplicons is not None
    assert isinstance(golden_gate.target_tm, int)
    assert callable(golden_gate.tm_function)


def test_make_amplicons():
    list_of_amplicons = [
        Dseqrecord(
            "TCGGTTGCTACACCCCTGCCGCAACGTTGAAGGTCCCGGATTAGACTGGCTGGATCTATGCCGTGACACCCGTTATACTCCATTACCGTCTGTGGGTCAC"
        ),
        Dseqrecord(
            "AGCTTGTTGTGGACTGGATTGCCATTCTCTCAGTGTATTACGCAGGCCGGCGCACGGGTCCCATATAAACCTGTCATAGCTTACCTGACTCTACTTGGAA"
        ),
        Dseqrecord(
            "ATGTGGCTAGGCCTTTGCCCACGCACCTGATCGGTCCTCGTTTGCTTTTTAGGACCGGATGAACTACAGAGCATTGCAAGAATCTCTACCTGCTTTACAA"
        ),
    ]
    amplicons = make_amplicons(list_of_amplicons, target_tm=55, limit=10)
    assert len(amplicons) == 3
    assert all(isinstance(i, Dseqrecord) for i in list_of_amplicons), (
        "All elements in 'list_of_amplicons' must be Dseqrecord objects"
    )
    assert all(isinstance(i, Amplicon) for i in amplicons), (
        "All elements in 'list_of_amplicons' must be Dseqrecord objects"
    )


def test_make_golden_gate_overhangs_forward_primers():
    # This is the 82nt sgRNA handle + 28nt csy4
    list_of_amplicons = [
        Dseqrecord(
            "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
        )
    ] * 3
    amplicons = make_amplicons(list_of_amplicons, target_tm=55, limit=10)
    assert len(amplicons) == 3

    # protospacers
    sgRNAs = [
        Dseqrecord("aGTGCTGGATCCTATTCCAG"),
        Dseqrecord("CGATGAGAGGAGTATTCGTC"),
        Dseqrecord("CGGGATGTTTTATCTAAACA"),
    ]

    # the restriction enzyme (RE) handle, the BsaI recognition site is in lowercase letters
    restriction_overhang_f = "GATCGggtctcc"
    # 28nt Csy4 recognition site that will be appended to the forward primer only in the first construct
    cys4 = "gTTCACTGCCGTATAGGCAGCTAAGAAA"
    # what fits with the backbone plasmid
    backbone_overhang_f = "cATG"

    amplicons_w_f_overhang = make_golden_gate_overhangs_forward_primers(
        amplicons,
        sgRNAs,
        restriction_overhang_f=restriction_overhang_f,
        backbone_overhang_f=backbone_overhang_f,
        cys4=cys4,
    )
    assert len(amplicons_w_f_overhang) == 3

    # Construct the expected sequence for the first amplicon
    expected_seq_0 = (
        restriction_overhang_f
        + backbone_overhang_f
        + cys4
        + str(sgRNAs[0].seq)
        + "GTTTTAGAGCTAGAAATAGC"
    )
    actual_seq_0 = str(amplicons_w_f_overhang[0].forward_primer.seq)

    assert len(actual_seq_0) == len(expected_seq_0), (
        f"Expected length {len(expected_seq_0)}, got {len(actual_seq_0)}"
    )
    assert actual_seq_0 == expected_seq_0, (
        f"Expected sequence: {expected_seq_0}, got: {actual_seq_0}"
    )

    # Construct the expected sequence for the second amplicon
    expected_seq_1 = (
        restriction_overhang_f + str(sgRNAs[1].seq) + "GTTTTAGAGCTAGAAATAGC"
    )
    actual_seq_1 = str(amplicons_w_f_overhang[1].forward_primer.seq)

    assert len(actual_seq_1) == len(expected_seq_1), (
        f"Expected length {len(expected_seq_1)}, got {len(actual_seq_1)}"
    )
    assert actual_seq_1 == expected_seq_1, (
        f"Expected sequence: {expected_seq_1}, got: {actual_seq_1}"
    )

    # Construct the expected sequence for the third amplicon
    expected_seq_2 = (
        restriction_overhang_f + str(sgRNAs[2].seq) + "GTTTTAGAGCTAGAAATAGC"
    )
    actual_seq_2 = str(amplicons_w_f_overhang[2].forward_primer.seq)

    assert len(actual_seq_2) == len(expected_seq_2), (
        f"Expected length {len(expected_seq_2)}, got {len(actual_seq_2)}"
    )
    assert actual_seq_2 == expected_seq_2, (
        f"Expected sequence: {expected_seq_2}, got: {actual_seq_2}"
    )


def test_make_GG_overhangs_reverse_primers():
    # This is the 82nt sgRNA handle
    list_of_amplicons = [
        Dseqrecord(
            "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT"
        )
    ] * 3
    amplicons = make_amplicons(list_of_amplicons, target_tm=55, limit=10)
    assert len(amplicons) == 3

    # protospacers
    sgRNAs = [
        Dseqrecord("gcggCGAACCAGCCCATCAT"),
        Dseqrecord("cccccgatGAGAGGAGTATTCGTC"),
        Dseqrecord("CGGGATGTTTTATCTAAACA"),
    ]

    # the restriction enzyme (RE) handle, the BsaI recognition site is in lowercase letters
    restriction_overhang_f = "GATCGggtctcc"
    # 28nt Csy4 recognition site that will be appended the the forward primer only in the first construct
    backbone_overhang_f = "cATG"

    amplicons_w_f_overhang = make_golden_gate_overhangs_forward_primers(
        amplicons,
        sgRNAs,
        restriction_overhang_f=restriction_overhang_f,
        backbone_overhang_f=backbone_overhang_f,
    )

    # FORWARD AND REVERSE PRIMERS
    restriction_overhang_r = "GATCAGGTCTCg"
    backbone_overhang_r = "cTAG"

    amplicons_w_f_r_primer_overhang = make_golden_gate_overhangs_reverse_primers(
        amplicons_w_f_overhang, sgRNAs, backbone_overhang_r, restriction_overhang_r
    )

    assert len(amplicons_w_f_r_primer_overhang) == 3

    # Check first amplicon reverse primer
    expected_seq_0 = (
        str(restriction_overhang_r)
        + str(sgRNAs[1][0:4].seq.reverse_complement())
        + str(Dseqrecord("AGTCGGTGCTTTTTT").seq.reverse_complement())
    )
    assert len(amplicons_w_f_r_primer_overhang[0].reverse_primer.seq) == len(
        expected_seq_0
    )
    assert str(amplicons_w_f_r_primer_overhang[0].reverse_primer.seq) == expected_seq_0

    # Check second amplicon reverse primer
    expected_seq_1 = (
        str(restriction_overhang_r)
        + str(sgRNAs[2][0:4].seq.reverse_complement())
        + str(Dseqrecord("AGTCGGTGCTTTTTT").seq.reverse_complement())
    )
    assert len(amplicons_w_f_r_primer_overhang[1].reverse_primer.seq) == len(
        expected_seq_1
    )
    assert str(amplicons_w_f_r_primer_overhang[1].reverse_primer.seq) == expected_seq_1

    # Check third amplicon reverse primer
    expected_seq_2 = (
        str(restriction_overhang_r)
        + str(Dseqrecord(backbone_overhang_r).seq.reverse_complement())
        + str(Dseqrecord("AGTCGGTGCTTTTTT").seq.reverse_complement())
    )
    assert len(amplicons_w_f_r_primer_overhang[2].reverse_primer.seq) == len(
        expected_seq_2
    )
    assert str(amplicons_w_f_r_primer_overhang[2].reverse_primer.seq) == expected_seq_2


def test_digest_amplicons_w_BsaI():
    digested_amplicons_w_BsaI = digest_amplicons_w_BsaI(list_of_amplicons)
    digested_length = [len(digest.seq) for digest in digested_amplicons_w_BsaI]

    assert len(digested_amplicons_w_BsaI) == 7
    assert digested_length == [166, 134, 134, 134, 134, 134, 134]
    # First amplicon
    assert str(digested_amplicons_w_BsaI[0][0:10].seq).upper() == "CATGGTTCAC"
    assert str(digested_amplicons_w_BsaI[0][-10:].seq).upper() == "AAGAAAGCGC"


def test_create_overhang_dataframe():
    overhang_df = create_overhang_dataframe(list_of_amplicons)
    test_data = {
        "Amplicon Name": [
            "Amplicon_SCO5090_140_pcr_products",
            "Amplicon_SCO5089_167_pcr_products",
            "Amplicon_SCO5090_174_pcr_products",
            "Amplicon_SCO5090_208_pcr_products",
            "Amplicon_SCO5087_210_pcr_products",
            "Amplicon_SCO5087_522_pcr_products",
            "Amplicon_SCO5090_640_pcr_products",
        ],
        "5' Overhang": ["CATG", "GCGC", "CGCC", "GGGA", "GAGC", "GTCG", "GACC"],
        "5' Duplicate": ["No", "No", "No", "No", "No", "No", "No"],
        "3' Overhang": ["GCGC", "CGCC", "GGGA", "GAGC", "GTCG", "GACC", "CTAG"],
        "3' Duplicate": ["No", "No", "No", "No", "No", "No", "No"],
    }

    expected_df = pd.DataFrame(test_data)

    pd.testing.assert_frame_equal(overhang_df, expected_df, check_dtype=False)
