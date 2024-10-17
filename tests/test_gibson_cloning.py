import pytest
from Bio import SeqIO 
from pydna.dseqrecord import Dseqrecord
from unittest.mock import MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import StuI
from streptocad.utils import list_of_objects_in_a_dir
import pydna 
from streptocad.utils import polymerase_dict

# Local module imports
from streptocad.cloning.gibson_cloning import (
    find_up_dw_repair_templates,
    update_primer_names, 
    assemble_single_plasmid_with_repair_templates,
    assemble_multiple_plasmids_with_repair_templates_for_deletion, 
)

@pytest.fixture
def genome():
    genome_path = 'tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb'
    return SeqIO.read(genome_path, "genbank")

@pytest.fixture
def repair_templates():
    return ['SCO5892']

@pytest.fixture
def repair_DNA_templates(genome, repair_templates):
    repair_templates_data = find_up_dw_repair_templates(genome, 
                                                        repair_templates, 
                                                        target_tm=60, 
                                                        
                                                        primer_tm_kwargs={'conc':0.4, 'prodcode':polymerase_dict['Phusion High-Fidelity DNA Polymerase (GC Buffer)']
} )
    return repair_templates_data

@pytest.fixture
def CRISPR_plasmids_wo_repair():

    plasmid = [SeqIO.read(f'tests/test_files/pCas9_SCO5892(27)_0.gb', format='genbank')]
    return plasmid

@pytest.fixture
def digested_plasmids(CRISPR_plasmids_wo_repair):
    return [Dseqrecord(plasmid, circular=True).cut(StuI)[0] for plasmid in CRISPR_plasmids_wo_repair]

@pytest.fixture
def list_of_gene_names():
    return ["SCO5892"]

@pytest.fixture
def assembled_records(list_of_gene_names, digested_plasmids, repair_DNA_templates):

    digested_plasmids[0].name = 'pCas9_SCO5892(27)_0'
    plasmids = assemble_multiple_plasmids_with_repair_templates_for_deletion(list_of_gene_names, 
                                                                         digested_plasmids, 
                                                                         repair_DNA_templates, 
                                                                         overlap=40)

    return plasmids



##### TEST ####

def test_find_up_dw_repair_templates(repair_DNA_templates):
    # Assert the expected output
    assert len(repair_DNA_templates) == 1  # This is a list of dicts
    record1 = repair_DNA_templates[0]
    assert record1["name"] == "SCO5892"
    assert len(record1["up_repair"]) == 1000
    assert str(record1["up_repair"][:20].seq) == 'CGACGAGCTGGACGTCGGCC'
    assert isinstance(record1['up_repair'], Dseqrecord)
    
    # UP Primer Checks
    assert str(record1["up_forwar_p"].seq) == "CGACGAGCTGGAC"
    assert str(record1["up_reverse_p"].seq) == "CTACCGGGCCGTT"
    assert isinstance(record1['up_forwar_p'], pydna.primer.Primer)
    assert isinstance(record1['up_reverse_p'], pydna.primer.Primer)
    
    # DW Primer Checks
    assert len(record1["dw_repair"]) == 1000
    assert str(record1["dw_forwar_p"].seq) == "TCCATGTCCTCACTCAG"
    assert str(record1["dw_reverse_p"].seq) == "CGGTGCGCCGCAT"



def test_assemble_multiple_plasmids_with_repair_templates_for_deletion(assembled_records):    
    # Check the length of the list of records
    print(assembled_records)
    assert len(assembled_records) == 1, "Expected only one record in the output."

    # Check the details of the gene in the first record
    record = assembled_records[0]
    assert record['gene_name'] == 'SCO5892', "The gene name in the record should be 'SCO5892'."
    assert len(record['contig']) == 13279, "The length of the contig should be 13279."

    # UP Primer Checks
    assert record['up_forwar_p_name'] == 'f1000'
    assert record['up_reverse_p_name'] == 'r1000'
    assert str(record['up_forwar_p_anneal']) == 'CGACGAGCTGGAC'
    assert str(record['up_reverse_p_anneal']) == 'CTACCGGGCCGTT'
    
    # DW Primer Checks
    assert record['dw_forwar_p_name'] == 'f1000'
    assert record['dw_reverse_p_name'] == 'r1000'
    assert str(record['dw_forwar_p_anneal']) == 'TCCATGTCCTCACTCAG'
    assert str(record['dw_reverse_p_anneal']) == 'CGGTGCGCCGCAT'

def test_update_primer_names(assembled_records):
    update_primer_names(assembled_records)
    # Check the details of the gene in the first record
    record = assembled_records[0]
    assert record['up_forwar_p_name'] == 'SCO5892_repair_up_forwar_p'
    assert record['up_reverse_p_name'] == 'SCO5892_repair_up_reverse_p'
    assert record['dw_forwar_p_name'] == 'SCO5892_repair_dw_forwar_p'
    assert record['dw_reverse_p_name'] == 'SCO5892_repair_dw_reverse_p'

def test_assemble_single_plasmid_with_repair_templates(CRISPR_plasmids_wo_repair, repair_DNA_templates):
    # to simulate that it will use the correct repair
    vector = Dseqrecord(CRISPR_plasmids_wo_repair[0], circular=True)
    vector.name = 'pCRISPR–1_SCO5892'
    vector_StuI = vector.cut(StuI)[0]
    vector_StuI.name = 'pCRISPR–1_SCO5892'

    repair_names = ['up_SCO5892', 'dw_SCO5892']
    repair_templates_sc000x = [repair_DNA_templates[0]['up_repair'], repair_DNA_templates[0]['dw_repair']]
    repair_templates_sc000x[0].name, repair_templates_sc000x[1].name = repair_names[0], repair_names[1]
    first_vector = assemble_single_plasmid_with_repair_templates(repair_templates_sc000x, vector_StuI, overlap=35)

    # Length assertions
    expected_lengths = [11279, 1053, 1053, 11279]
    for i, fragment in enumerate(first_vector):
        assert len(fragment) == expected_lengths[i]

    # Sequence content assertions
    assert str(first_vector[0].seq).startswith('CCT')
    assert str(first_vector[0].seq).endswith('AGG')

    # Name assertion
    assert first_vector[0].name == 'pCRISPR–1_SCO5892'

def test_find_up_dw_repair_templates_edge_cases():
    # Mock genome as SeqRecord
    mock_genome = SeqRecord(Seq("ATGCGACTACGATCGAGCGATATCGCGATCGATCGAGCGATCGATCGA"), id="MockGenome")

    # Mock feature with locus_tag
    mock_feature = SeqFeature(FeatureLocation(0, 10), type='CDS', qualifiers={'locus_tag': ['MOCK001']})
    mock_genome.features.append(mock_feature)

    # Empty repair_templates list
    empty_repair_templates = []
    results = find_up_dw_repair_templates(mock_genome, empty_repair_templates)
    assert len(results) == 0, "Expected no results for empty repair templates."

    # Non-existent locus tag
    nonexistent_locus_tag = ['MOCK999']
    results = find_up_dw_repair_templates(mock_genome, nonexistent_locus_tag)
    assert len(results) == 0, "Expected no results for non-existent locus tags."

def test_assemble_multiple_plasmids_with_repair_templates_for_deletion_edge_cases():
    # Mock gene names
    mock_gene_names_empty = []
    mock_gene_names_nonexistent = ["MOCK999"]

    # Mock repair templates
    mock_repair_templates = [{'name': 'MOCK001', 'up_repair': MagicMock(), 'dw_repair': MagicMock()}]

    # Mock digested plasmids
    mock_digested_plasmids_empty = []
    mock_digested_plasmids_nonexistent = [MagicMock(name="MockPlasmid")]
    mock_digested_plasmids_nonexistent[0].name = "PlasmidForMOCK999"

    # Test with an empty list of gene names
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        mock_gene_names_empty, mock_digested_plasmids_nonexistent, mock_repair_templates)
    assert len(results) == 0, "Expected no results for empty gene names."

    # Test with gene names that do not have corresponding repair templates
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        mock_gene_names_nonexistent, mock_digested_plasmids_nonexistent, mock_repair_templates)
    assert len(results) == 0, "Expected no results for non-existent gene names."

    # Test with an empty list of digested plasmids
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        mock_gene_names_nonexistent, mock_digested_plasmids_empty, mock_repair_templates)
    assert len(results) == 0, "Expected no results for empty list of digested plasmids."

    # Test with digested plasmids that do not match any gene names
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        mock_gene_names_nonexistent, mock_digested_plasmids_nonexistent, mock_repair_templates)
    assert len(results) == 0, "Expected no results for digested plasmids that do not match any gene names."