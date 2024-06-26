
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

from src.teemi_functions.cloning.golden_gate_cloning import (
    GoldenGateCloning, 
    make_amplicons, 
    make_golden_gate_overhangs_forward_primers, 
    make_golden_gate_overhangs_reverse_primers,
    digest_amplicons_w_BsaI,
    create_overhang_dataframe
)

from src.teemi_functions.sgRNA import (
    make_antismash_df_to_dseq
)


# Input
sco5087 = pd.read_csv('data/CRISPy_results/SCO5087_BEST_STOP.csv')
pCRISPR_plasmid_MCBE = read_genbank_files('data/plasmids/pCRISPR-MCBE_Csy4_kasopGFP.gb')[0]
vector = Dseqrecord(pCRISPR_plasmid_MCBE, circular = True)
vector.name = 'pCRISPR-MCBE_Csy4'

#Calculations: 
sgRNA_list = make_antismash_df_to_dseq(sco5087)
# This is the 82nt sgRNA handle + 28nt csy4
sgRNA_handle_cys4_sites = [Dseqrecord('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTTCACTGCCGTATAGGCAGCTAAGAAA', name = 'sgRNA_handle_cys4')]*len(sgRNA_list)

golden_gate = GoldenGateCloning(sgRNA_list,
                                sgRNA_handle_cys4_sites, 
                                target_tm= 60) 

# make amplicons
list_of_amplicons = golden_gate.simulate_pcrs()

# Making a vector 
digest = digest_amplicons_w_BsaI(list_of_amplicons)

# Digest vector
from Bio.Restriction import NcoI,NheI
vector_NcoI_NheI = vector.cut(NcoI, NheI)[1]
rec_vec = vector_NcoI_NheI

# assemble vecor
for i in range(len(digest)):
    rec_vec += digest[i]
rec_vec = rec_vec.looped()


#######################
##### TESTING ########
#######################

def test_amplicon_lengths():
    predicted_length_of_amplicons = [190, 158, 158, 158, 158]
    list_of_amplicons = golden_gate.simulate_pcrs()
    length_of_amplicons = [len(amplicon) for amplicon in list_of_amplicons]
    assert predicted_length_of_amplicons == length_of_amplicons
    assert all(isinstance(item, Dseqrecord) for item in list_of_amplicons)
    assert len(list_of_amplicons) > 0


def test_melting_temperatures():
    predicted_melting_temperatures = [(59.84033214711462, 60.61636492495495),
                                    (59.84033214711462, 60.61636492495495),
                                    (59.84033214711462, 60.61636492495495),
                                    (59.84033214711462, 60.61636492495495),
                                    (59.84033214711462, 60.61636492495495)]
    list_of_melting_temperatures = golden_gate.calculate_melting_temperature()
    assert predicted_melting_temperatures == list_of_melting_temperatures

def test_pcr_dataframe():
    data = {
        'amplicon_name': [
            'amplicon_CY00000014', 'amplicon_CY00000058',
            'amplicon_CY00000155', 'amplicon_CY00000196',
            'amplicon_CY00000239'
        ],
        'f_primer_id': ['P001', 'P003', 'P005', 'P007', 'P009'],
        'r_primer_id': ['P002', 'P004', 'P006', 'P008', 'P010'],
        'template_for_amplification': ['sgRNA_handle_cys4'] * 5,
        'f_primer tm': [59.840332] * 5,
        'r_primer_tm': [60.616365] * 5,
        'ta': [62] * 5
    }
    expected_df = pd.DataFrame(data)
    pd.testing.assert_frame_equal(golden_gate.make_pcr_df(), expected_df, check_dtype=False)

def test_primer_dataframe():
    expected_primer_dataframe = pd.DataFrame({
        'id': [
            'P001_CY00000014', 'P002_CY00000014',
            'P003_CY00000058', 'P004_CY00000058',
            'P005_CY00000155', 'P006_CY00000155',
            'P007_CY00000196', 'P008_CY00000196',
            'P009_CY00000239', 'P010_CY00000239'
        ],
        'sequence': [
            'GATCGggtctcccATGgTTCACTGCCGTATAGGCAGCTAAGAAAGAGCAGTTCCCAGAACTGCCGTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'GATCAGGTCTCgCGACTTTCTTAGCTGCCTATACGG',
            'GATCGggtctccGTCGACCTCCCAGTCACGGCGTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'GATCAGGTCTCgCGGCTTTCTTAGCTGCCTATACGG',
            'GATCGggtctccGCCGACCGCCCAGGCGACCTGTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'GATCAGGTCTCgACGGTTTCTTAGCTGCCTATACGG',
            'GATCGggtctccCCGTTCACAGGTCGCGGCGGGTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'GATCAGGTCTCgCGGTTTTCTTAGCTGCCTATACGG',
            'GATCGggtctccACCGCCCAGGCGACCTCGGCGTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'GATCAGGTCTCgCTAgTTTCTTAGCTGCCTATACGG'
        ],
        'len_primer': [92, 36, 60, 36, 60, 36, 60, 36, 60, 36],
        'annealing_tm': [59.84033214711462, 
                         60.61636492495495, 
                         59.84033214711462, 
                         60.61636492495495, 
                         59.84033214711462, 
                         60.61636492495495, 
                         59.84033214711462, 
                         60.61636492495495, 
                         59.84033214711462, 
                         60.61636492495495],
        'primer_type': [
            'Forward', 'Reverse', 'Forward', 'Reverse',
            'Forward', 'Reverse', 'Forward', 'Reverse',
            'Forward', 'Reverse'
        ],
        'annealing_sequence': [
            'GTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'TTTCTTAGCTGCCTATACGG',
            'GTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'TTTCTTAGCTGCCTATACGG',
            'GTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'TTTCTTAGCTGCCTATACGG',
            'GTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'TTTCTTAGCTGCCTATACGG',
            'GTTTTAGAGCTAGAAATAGCAAGTTAAA',
            'TTTCTTAGCTGCCTATACGG'
        ],
        'template': ['sgRNA_handle_cys4'] * 10
    })

    pd.testing.assert_frame_equal(golden_gate.make_primer_df(), expected_primer_dataframe)


def test_initialization():
    # Test the successful initialization of GoldenGateCloning instance
    assert golden_gate.sgRNAs is not None
    assert golden_gate.list_of_amplicons is not None
    assert isinstance(golden_gate.target_tm, int)
    assert callable(golden_gate.tm_function)



def test_make_amplicons():

    list_of_amplicons =[Dseqrecord('TCGGTTGCTACACCCCTGCCGCAACGTTGAAGGTCCCGGATTAGACTGGCTGGATCTATGCCGTGACACCCGTTATACTCCATTACCGTCTGTGGGTCAC'), 
                        Dseqrecord('AGCTTGTTGTGGACTGGATTGCCATTCTCTCAGTGTATTACGCAGGCCGGCGCACGGGTCCCATATAAACCTGTCATAGCTTACCTGACTCTACTTGGAA'), 
                        Dseqrecord('ATGTGGCTAGGCCTTTGCCCACGCACCTGATCGGTCCTCGTTTGCTTTTTAGGACCGGATGAACTACAGAGCATTGCAAGAATCTCTACCTGCTTTACAA')] 
    amplicons = make_amplicons(list_of_amplicons, 
                               target_tm = 55,
                                limit= 10)     
    assert len(amplicons) == 3
    assert all(isinstance(i, Dseqrecord) for i in list_of_amplicons), "All elements in 'list_of_amplicons' must be Dseqrecord objects"
    assert all(isinstance(i, Amplicon) for i in amplicons), "All elements in 'list_of_amplicons' must be Dseqrecord objects"



def test_make_golden_gate_overhangs_forward_primers():
    # This is the 82nt sgRNA handle + 28nt csy4
    list_of_amplicons =[Dseqrecord('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT')]*3 
    amplicons = make_amplicons(list_of_amplicons, 
                               target_tm = 55,
                                 limit= 10)
    assert len(amplicons) == 3
    # protospacers
    sgRNAs = [Dseqrecord('aGTGCTGGATCCTATTCCAG'),Dseqrecord('CGATGAGAGGAGTATTCGTC'),Dseqrecord('CGGGATGTTTTATCTAAACA')]
    # the restriction enzyme (RE) handle, the BsaI recognition site is in lowercase letters
    restriction_overhang_f = "GATCGggtctcc"
    # 28nt Csy4 recognition site that will be appended the the forward primer only in the first construct
    cys4= "gTTCACTGCCGTATAGGCAGCTAAGAAA"
    # what fits with the backbone plasmid 
    backbone_overhang_f = "cATG"

    amplicons_w_f_overhang = make_golden_gate_overhangs_forward_primers(amplicons,
                                                           sgRNAs,
            restriction_overhang_f=restriction_overhang_f,
            backbone_overhang_f=backbone_overhang_f,
            cys4=cys4,
        )
    assert len(amplicons_w_f_overhang) == 3
    # First have csy4 sites. The string below is the annealing part of the primer
    assert len(amplicons_w_f_overhang[0].forward_primer.seq) == len(str(restriction_overhang_f)+str(backbone_overhang_f)+str(cys4)+str(sgRNAs[0].seq)+ 'GTTTTAGAGCTAGAAATAGCA')
    assert str(amplicons_w_f_overhang[0].forward_primer.seq) == str(str(restriction_overhang_f)+str(backbone_overhang_f)+str(cys4)+str(sgRNAs[0].seq)+ 'GTTTTAGAGCTAGAAATAGCA')
    # After that no csy4 sites and no backbone overhang bc it comes from the protspacer.
    assert len(amplicons_w_f_overhang[1].forward_primer.seq) == len(str(restriction_overhang_f)+str(sgRNAs[1].seq)+ 'GTTTTAGAGCTAGAAATAGCA')
    assert str(amplicons_w_f_overhang[1].forward_primer.seq) == str(str(restriction_overhang_f)+str(sgRNAs[1].seq)+ 'GTTTTAGAGCTAGAAATAGCA')
    # After that no csy4 sites. The string below is the annealing part of the primer
    assert len(amplicons_w_f_overhang[2].forward_primer.seq) == len(str(restriction_overhang_f)+str(sgRNAs[2].seq)+ 'GTTTTAGAGCTAGAAATAGCA')
    assert str(amplicons_w_f_overhang[2].forward_primer.seq) == str(str(restriction_overhang_f)+str(sgRNAs[2].seq)+ 'GTTTTAGAGCTAGAAATAGCA')


def test_make_GG_overhangs_reverse_primers():
        # This is the 82nt sgRNA handle 
    list_of_amplicons =[Dseqrecord('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT')]*3 
    amplicons = make_amplicons(list_of_amplicons, 
                               target_tm = 55,
                                 limit= 10)
    assert len(amplicons) == 3
    # protospacers
    sgRNAs = [Dseqrecord('gcggCGAACCAGCCCATCAT'),Dseqrecord('cccccgatGAGAGGAGTATTCGTC'),Dseqrecord('CGGGATGTTTTATCTAAACA')]
    # the restriction enzyme (RE) handle, the BsaI recognition site is in lowercase letters
    restriction_overhang_f = "GATCGggtctcc"
    # 28nt Csy4 recognition site that will be appended the the forward primer only in the first construct
    # what fits with the backbone plasmid 
    backbone_overhang_f = "cATG"

    amplicons_w_f_overhang = make_golden_gate_overhangs_forward_primers(amplicons,
                                                           sgRNAs,
            restriction_overhang_f=restriction_overhang_f,
            backbone_overhang_f=backbone_overhang_f,

        )
    # FORWARD AND REVERSE PRIMERS 
    restriction_overhang_r: str = "GATCAGGTCTCg"
    backbone_overhang_r: str = "cTAG"


    amplicons_w_f_r_primer_overhang = make_golden_gate_overhangs_reverse_primers(amplicons_w_f_overhang,
                                               sgRNAs, 
                                               backbone_overhang_r, 
                                               restriction_overhang_r )
    

    assert len(amplicons_w_f_r_primer_overhang) == 3
    # FIRST we have the restriction overhang, then overlapping 4nt sgRNA from the next sgRNA,  and finally the annealing
    assert len(amplicons_w_f_r_primer_overhang[0].reverse_primer.seq) == len(str(restriction_overhang_r)+ str(sgRNAs[1][0:4].seq.reverse_complement())+str(Dseqrecord('CGAGTCGGTGCTTTTTT').seq.reverse_complement()))
    assert str(amplicons_w_f_r_primer_overhang[0].reverse_primer.seq) == str(str(restriction_overhang_r)+ str(sgRNAs[1][0:4].seq.reverse_complement())+str(Dseqrecord('CGAGTCGGTGCTTTTTT').seq.reverse_complement()))
    # FIRST we have the restriction overhang, then overlapping 4nt sgRNA from the next sgRNA,  and finally the annealing
    assert len(amplicons_w_f_r_primer_overhang[1].reverse_primer.seq) == len(str(restriction_overhang_r)+ str(sgRNAs[2][0:4].seq.reverse_complement())+str(Dseqrecord('CGAGTCGGTGCTTTTTT').seq.reverse_complement()))
    assert str(amplicons_w_f_r_primer_overhang[1].reverse_primer.seq) == str(str(restriction_overhang_r)+ str(sgRNAs[2][0:4].seq.reverse_complement())+str(Dseqrecord('CGAGTCGGTGCTTTTTT').seq.reverse_complement()))
    # Finally, we have the restriction overhang, then the last reverse have a fixed backbone_overhang_r, and annealing ofc
    assert len(amplicons_w_f_r_primer_overhang[2].reverse_primer.seq) == len(str(restriction_overhang_r)+ str(Dseqrecord(backbone_overhang_r).seq.reverse_complement())+str(Dseqrecord('CGAGTCGGTGCTTTTTT').seq.reverse_complement()))
    assert str(amplicons_w_f_r_primer_overhang[2].reverse_primer.seq) == str(str(restriction_overhang_r)+ str(Dseqrecord(backbone_overhang_r).seq.reverse_complement())+str(Dseqrecord('CGAGTCGGTGCTTTTTT').seq.reverse_complement()))


def test_digest_amplicons_w_BsaI():
    digested_amplicons_w_BsaI = digest_amplicons_w_BsaI(list_of_amplicons)
    digested_length = [len(digest.seq) for digest in digested_amplicons_w_BsaI]

    assert len(digested_amplicons_w_BsaI) == 5
    assert digested_length == [166, 134, 134, 134, 134]
    # First amplicon
    assert str(digested_amplicons_w_BsaI[0][0:10].seq).upper() == 'CATGGTTCAC'
    assert str(digested_amplicons_w_BsaI[0][-10:].seq).upper() == 'AAGAAAGTCG'


def test_create_overhang_dataframe():
    overhang_df = create_overhang_dataframe(list_of_amplicons)
    test_data = {
        'Amplicon Name': [
            'Amplicon_CY00000014_pcr_products', 'Amplicon_CY00000058_pcr_products',
            'Amplicon_CY00000155_pcr_products', 'Amplicon_CY00000196_pcr_products',
            'Amplicon_CY00000239_pcr_products'    
        ],
        '5\' Overhang': ['cATG', 'GTCG', 'GCCG','CCGT','ACCG' ],
        '5\' Duplicate': ['No']*5,
        '3\' Overhang': ['GTCG', 'GCCG', 'CCGT', 'ACCG', 'cTAG'] ,
        '3\' Duplicate': ['No'] * 5,
    }
    expected_df = pd.DataFrame(test_data)

    pd.testing.assert_frame_equal(overhang_df, expected_df, check_dtype=False)
