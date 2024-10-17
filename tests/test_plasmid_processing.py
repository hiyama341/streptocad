
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from teemi.build.PCR import primer_tm_neb
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO 
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Restriction import StuI
from pydna.dseqrecord import Dseqrecord
import pandas as pd
from teemi.design.fetch_sequences import read_genbank_files

from streptocad.cloning.pcr_simulation import perform_pcr_on_sequences
from streptocad.cloning.plasmid_processing import (
    assemble_and_process_plasmids,
    assemble_plasmids_by_ssDNA_bridging, 
    annotate_plasmid_with_sgrnas, 
    check_plasmid_restriction_sites, 
)


@pytest.fixture
def plasmid_path():
    return 'tests/test_files/pOEX-PkasO.gb'

@pytest.fixture
def clean_plasmid(plasmid_path):
    # Read the plasmid file and return as Dseqrecord
    record = SeqIO.read(plasmid_path, "genbank")
    return Dseqrecord(record.seq, circular=True)


@pytest.fixture
def ssDNA_plasmid_cas9():
    # Read the plasmid file and return as Dseqrecord
    record =read_genbank_files('tests/test_files/pCRISPR-cBEST.gbk')[0]
    return Dseqrecord(record, circular=True)


@pytest.fixture
def goe_regulators_path():
    return 'tests/test_files/GOE_regulators.gb'

@pytest.fixture
def list_of_amplicons(goe_regulators_path, ):
    # Read the GOE regulators GenBank file and return the first record as a Dseqrecord
    clean_seq = read_genbank_files(goe_regulators_path)
    primer_df = pd.read_csv('tests/test_files/primer_df_workflow1.csv')

    list_of_amplicons = perform_pcr_on_sequences(primer_df, 
                                             clean_seq)
    return list_of_amplicons


def test_assemble_and_process_plasmids(clean_plasmid, list_of_amplicons):
    plasmids, assembly_results = assemble_and_process_plasmids(clean_plasmid, 
                                                               list_of_amplicons, 
                                                               enzymes=[StuI], 
                                                               save_plasmids=False)
    print(plasmids)
    # Check if the output is as expected
    assert len(plasmids) == len(list_of_amplicons)
    for plasmid in plasmids:
        assert plasmid.name.startswith("pOEx-KasO_")
        assert isinstance(plasmid, Dseqrecord)
        assert plasmid.circular

def test_assemble_plasmids_by_ssDNA_bridging(ssDNA_plasmid_cas9):

    # cut plasmid
    from Bio.Restriction import NcoI, NheI
    linearized_plasmid = sorted(ssDNA_plasmid_cas9.cut(NcoI,NheI ), key=lambda x: len(x), reverse=True)[0]

    # Convert amplicon to ssDNA primers for the test
    ssDNA_primers = [Dseqrecord('cggttggtaggatcgacggcCCGTTCACAGGTCGCGGCGGtatgtcctccgagaccggcc'), 
                     Dseqrecord('cggttggtaggatcgacggcGTCGACCTCCCAGTCACGGCtatgtcctccgagaccggcc')]
    
    sgrna_plasmids = assemble_plasmids_by_ssDNA_bridging(ssDNA_primers, linearized_plasmid)
    
    #print(sgrna_plasmids[0].seq)
    assert len(sgrna_plasmids) == len(ssDNA_primers)
    for plasmid in sgrna_plasmids:
        assert isinstance(plasmid, Dseqrecord)
        assert plasmid.seq.startswith('TATGTCCTCCGAGACCGGCCCCGTGGCCGTGGACCCGACCCTGCGTCGCCGCATCGAGCCGCACGAGTT') 

@pytest.fixture
def sgrna_annotation_df():
    data = {
        "strain_name": ["NC_003888.3", "NC_003888.3"],
        "locus_tag": ["SCO5087", "SCO5087"],
        "gene_loc": [5529801, 5529801],
        "gene_strand": [1, 1],
        "sgrna_strand": [1, -1],
        "sgrna_loc": [283, 522],
        "gc": [0.75, 0.70],
        "pam": ["AGG", "CGG"],
        "sgrna": ["CCGTTCACAGGTCGCGGCGG", "GTCGACCTCCCAGTCACGGC"],
        "sgrna_seed_sequence": ["CAGGTCGCGGCGG", "TCCCAGTCACGGC"],
        "off_target_count": [0, 0],
        "editing_context": ["0,0", "0,0,0,0,0"],
        "editable_cytosines": ["6,8", "3,6,7,9,10"],
        "mutations": ["S90L, Q91*", "W171*, E172K, V173I, D174N"]
    }

    return pd.DataFrame(data)

def test_annotate_plasmid_with_sgrnas(clean_plasmid, sgrna_annotation_df, ssDNA_plasmid_cas9):
    annotated_plasmid = annotate_plasmid_with_sgrnas(clean_plasmid, sgrna_annotation_df)
    
    # this clean plasmid does not contain an sgrna
    assert len(annotated_plasmid.features) == 0


    # These plasmids shoudl

    # cut plasmid
    from Bio.Restriction import NcoI, NheI
    linearized_plasmid = sorted(ssDNA_plasmid_cas9.cut(NcoI,NheI ), key=lambda x: len(x), reverse=True)[0]

    # Convert amplicon to ssDNA primers for the test
    ssDNA_primers = [Dseqrecord('cggttggtaggatcgacggcCCGTTCACAGGTCGCGGCGGtatgtcctccgagaccggcc'), 
                     Dseqrecord('cggttggtaggatcgacggcGTCGACCTCCCAGTCACGGCtatgtcctccgagaccggcc')]
    
    sgrna_plasmids = assemble_plasmids_by_ssDNA_bridging(ssDNA_primers, linearized_plasmid)
    # annotate plasmids
    for plasmid in sgrna_plasmids: 
        annotate_plasmid_with_sgrnas(plasmid, sgrna_annotation_df)

    for feature in annotated_plasmid.features:
        assert feature.type == "sgRNA"
        assert "label" in feature.qualifiers

def test_check_plasmid_restriction_sites(clean_plasmid):
    enzymes = ["EcoRI", "BamHI"]
    result = check_plasmid_restriction_sites(clean_plasmid, enzymes)
    
    assert isinstance(result, dict)
    for enzyme in enzymes:
        assert enzyme in result
        assert isinstance(result[enzyme], int)