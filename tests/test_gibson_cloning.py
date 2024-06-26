from Bio import SeqIO 
from pydna.dseqrecord import Dseqrecord
from src.utils import list_of_objects_in_a_dir
import pydna

# Local module imports
from src.teemi_functions.cloning.gibson_cloning import (
    assemble_multiple_plasmids_with_repair_templates_for_deletion, 
    assemble_single_plasmid_with_repair_templates,
    find_up_dw_repair_templates,
    update_primer_names
)

# read in genome
genome_path = 'data/genomes/Streptomyces_coelicolor_A3_chromosome.gb'
StrepA3 = SeqIO.read(genome_path, "gb")
repair_templates = ['SCO5892']

# Get repair templates
repair_DNA_templates = find_up_dw_repair_templates(StrepA3, repair_templates)

# PLASMIDS
list_of_plasmids= list_of_objects_in_a_dir('data/plasmids/sgRNA_plasmids_pCRISPR–dCas9/SCO05087_SCO005892/')

# lets read in the plasmids
CRISPR_plasmids_wo_repair =[]


# we will test 2 plasmids
for record in list_of_plasmids[0:1]: 
    CRISPR_plasmids_wo_repair.append(SeqIO.read('data/plasmids/sgRNA_plasmids_pCRISPR–dCas9/SCO05087_SCO005892/'+record, format = 'gb'))
CRISPR_plasmids_wo_repair

# digest the plasmids
from Bio.Restriction import StuI
digested_plasmids = [Dseqrecord(digest, circular = True).cut(StuI)[0] for digest in CRISPR_plasmids_wo_repair]

# rename them appropriatly
for i in range(len( digested_plasmids)):
    digested_plasmids[i].name = CRISPR_plasmids_wo_repair[i].name

list_of_gene_names = ["SCO5892", "SCO5087"]
list_of_digested_plasmids = digested_plasmids
list_of_records = assemble_multiple_plasmids_with_repair_templates_for_deletion(list_of_gene_names,list_of_digested_plasmids, repair_DNA_templates, overlap = 40 )





def test_find_up_dw_repair_templates():

    # Assert the expected output
    assert len(repair_DNA_templates) == 1 # This is a list of dicts
    # Assert the content of the first record
    record1 = repair_DNA_templates[0]
    assert record1["name"] == "SCO5892"
    assert len(record1["up_repair"]) == 1000
    assert str(record1["up_repair"][:20].seq) == 'CGACGAGCTGGACGTCGGCC'
    assert isinstance(record1['up_repair'], Dseqrecord)
    # UP
    assert str(record1["up_forwar_p"].seq) == "CGACGAGCTGGACGT"
    assert str(record1["up_reverse_p"].seq) == "CTACCGGGCCGTTCC"
    assert record1["tm_up_forwar_p"] == 65
    assert record1["tm_up_reverse_p"] == 66
    assert record1["location_up_start"] == 6448548
    assert record1["location_up_end"] == 6449548
    # DW
    assert len(record1["dw_repair"]) == 1000
    assert record1["tm_dw_forwar_p"] == 64
    assert record1["tm_dw_reverse_p"] == 67
    assert record1["location_dw_start"] == 6456442
    assert record1["location_dw_end"] == 6457442










def test_assemble_multiple_plasmids_with_repair_templates_for_deletion():    
    # Check the length of the list of records
    assert len(list_of_records) == 1, "Expected only one record in the output."

    # Check the details of the gene in the first record
    record = list_of_records[0]
    assert record['gene_name'] == 'SCO5892', "The gene name in the record should be 'SCO5892'."
    assert len(record['contig']) == 13279, "The length of the contig should be 13279."

    # Check the details of primers and their properties in the record
    # UP
    assert record['up_forwar_p_name'] == 'f1000'
    assert record['up_reverse_p_name'] == 'r1000'
    assert record['up_forwar_primer_str'] == 'AAGGCCGCTTTTGCGGGATCTCGTCGAAGGCACTAGAAGGCGACGAGCTGGACGT'
    assert record['up_reverse_primer_str'] == 'TGGCTGAGTGAGGACATGGACTACCGGGCCGTTCC'
    assert str(record['up_forwar_p_anneal']) == 'CGACGAGCTGGACGT'
    assert str(record['up_reverse_p_anneal']) == 'CTACCGGGCCGTTCC' 
    assert record['tm_up_forwar_p'] == 65
    assert record['tm_up_reverse_p'] == 66

    # DW
    assert record['dw_forwar_p_name'] == 'f1000'
    assert record['dw_reverse_p_name'] == 'r1000'
    assert record['dw_forwar_primer_str'] == 'CGCTCGGAACGGCCCGGTAGTCCATGTCCTCACTCAGC'
    assert record['dw_reverse_primer_str'] == 'TCCCCGTCCGGGACCCGCGCGGTCGATCCCCGCATATAGGCGGTGCGCCGCAT'
    assert str(record['dw_forwar_p_anneal']) == 'TCCATGTCCTCACTCAGC'
    assert str(record['dw_reverse_p_anneal']) == 'CGGTGCGCCGCAT'
    assert record['tm_dw_forwar_p'] == 64
    assert record['tm_dw_reverse_p'] == 67

    # Instances
    assert isinstance(record, dict), "Each record should be a dictionary."
    assert isinstance(record['gene_name'], str), "The gene name should be a string."
    assert isinstance(record['contig'], Dseqrecord), "The contig should be a Dseqrecord object."
    assert isinstance(record['up_forwar_p'], pydna.primer.Primer), "The up forward primer should be a Dseqrecord object."
    assert isinstance(record['up_reverse_p'], pydna.primer.Primer), "The up reverse primer should be a Dseqrecord object."
    assert isinstance(record['tm_up_forwar_p'], (int, float)), "The melting temperature should be a number."
    assert isinstance(record['tm_up_reverse_p'], (int, float)), "The melting temperature should be a number."

    
def test_update_primer_names():
    update_primer_names(list_of_records)
        # Check the length of the list of records
    assert len(list_of_records) == 1

    # Check the details of the gene in the first record
    record = list_of_records[0]
    assert record['gene_name'] == 'SCO5892'
    assert len(record['contig']) == 13279

    # Check the details of primers and their properties in the record
    assert record['up_forwar_p_name'] == 'SCO5892_up_F0'
    assert record['up_reverse_p_name'] == 'SCO5892_up_R0'
    assert record['up_forwar_primer_str'] == 'AAGGCCGCTTTTGCGGGATCTCGTCGAAGGCACTAGAAGGCGACGAGCTGGACGT'
    assert record['up_reverse_primer_str'] == 'TGGCTGAGTGAGGACATGGACTACCGGGCCGTTCC'
    assert str(record['up_forwar_p_anneal']) == 'CGACGAGCTGGACGT'
    assert str(record['up_reverse_p_anneal']) == 'CTACCGGGCCGTTCC' 
    assert record['tm_up_forwar_p'] == 65
    assert record['tm_up_reverse_p'] == 66

    assert record['dw_forwar_p_name'] == 'SCO5892_dw_F0'
    assert record['dw_reverse_p_name'] == 'SCO5892_dw_R0'
    assert record['dw_forwar_primer_str'] == 'CGCTCGGAACGGCCCGGTAGTCCATGTCCTCACTCAGC'
    assert record['dw_reverse_primer_str'] == 'TCCCCGTCCGGGACCCGCGCGGTCGATCCCCGCATATAGGCGGTGCGCCGCAT'
    assert str(record['dw_forwar_p_anneal']) == 'TCCATGTCCTCACTCAGC'
    assert str(record['dw_reverse_p_anneal']) == 'CGGTGCGCCGCAT'
    assert record['tm_dw_forwar_p'] == 64
    assert record['tm_dw_reverse_p'] == 67



def test_assemble_single_plasmid_with_repair_templates():
    vector = Dseqrecord(CRISPR_plasmids_wo_repair[0], circular = True)
    vector.name = 'pCRISPR–1_SCO5892'
    from Bio.Restriction import StuI
    vector_StulI = vector.cut(StuI)[0]
    vector_StulI.name = 'pCRISPR–1_SCO5892'

    repair_names = ['up_SCO5892', 'dw_SCO5892']
    repair_templates_sc000x = [repair_DNA_templates[0]['up_repair'], repair_DNA_templates[0]['dw_repair']]
    repair_templates_sc000x[0].name, repair_templates_sc000x[1].name = repair_names[0], repair_names[1]
    first_vector = assemble_single_plasmid_with_repair_templates(repair_templates_sc000x,vector_StulI , overlap = 35)

    # Length assertions
    expected_lengths = [11279, 1053, 1053, 11279]
    for i, fragment in enumerate(first_vector):
        assert len(fragment) == expected_lengths[i]

    # Sequence content assertions
    assert str(first_vector[0].seq).startswith('CCT')
    assert str(first_vector[0].seq).endswith('AGG')

    # Names 
    assert first_vector[0].name == 'pCRISPR–1_SCO5892'

    # Overlap regions
    assert str(first_vector[0].seq[-35:]) == str(first_vector[1].seq[:35]), "Overlap region mismatch."



from unittest.mock import MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


# Mock genome as SeqRecord
mock_genome = SeqRecord(Seq("ATGCGACTACGATCGAGCGATATCGCGATCGATCGAGCGATCGATCGA"), id="MockGenome")

# Mock feature with locus_tag
mock_feature = SeqFeature(FeatureLocation(0, 10), type='CDS', qualifiers={'locus_tag': ['MOCK001']})
mock_genome.features.append(mock_feature)

# Non-existent locus tag
nonexistent_locus_tag = ['MOCK999']

# Empty repair_templates list
empty_repair_templates = []

def test_find_up_dw_repair_templates_empty():
    # Test with an empty list of repair templates
    results = find_up_dw_repair_templates(mock_genome, empty_repair_templates)
    assert len(results) == 0, "Expected no results for empty repair templates."

def test_find_up_dw_repair_templates_nonexistent_locus_tag():
    # Test with a locus tag that does not exist
    results = find_up_dw_repair_templates(mock_genome, nonexistent_locus_tag)
    assert len(results) == 0, "Expected no results for non-existent locus tags."

def test_assemble_multiple_empty_gene_names():
    # Test with an empty list of gene names
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion([], list_of_digested_plasmids, repair_DNA_templates)
    assert len(results) == 0, "Expected no results with empty gene names."

def test_assemble_multiple_empty_digested_plasmids():
    # Test with an empty list of digested plasmids
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(list_of_gene_names, [], repair_DNA_templates)
    assert len(results) == 0, "Expected no results with empty digested plasmids."

def test_assemble_multiple_empty_repair_templates():
    # Test with an empty list of repair templates
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(list_of_gene_names, list_of_digested_plasmids, [])
    assert len(results) == 0, "Expected no results with empty repair templates."

# Mock gene names
mock_gene_names_empty = []
mock_gene_names_nonexistent = ["MOCK999"]

# Mock repair templates
mock_repair_templates = [{'name': 'MOCK001', 'up_repair': MagicMock(), 'dw_repair': MagicMock()}]

# Mock digested plasmids
mock_digested_plasmids_empty = []
mock_digested_plasmids_nonexistent = [MagicMock(name="MockPlasmid")]

# Assigning a mock name that does not match the repair templates
mock_digested_plasmids_nonexistent[0].name = "PlasmidForMOCK999"

def test_assemble_multiple_empty_gene_names():
    # Test with an empty list of gene names
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        mock_gene_names_empty, mock_digested_plasmids_nonexistent, mock_repair_templates)
    assert len(results) == 0, "Expected no results for empty gene names."

def test_assemble_multiple_nonexistent_gene_names():
    # Test with gene names that do not have corresponding repair templates
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        mock_gene_names_nonexistent, mock_digested_plasmids_nonexistent, mock_repair_templates)
    assert len(results) == 0, "Expected no results for non-existent gene names."

def test_assemble_multiple_empty_digested_plasmids():
    # Test with an empty list of digested plasmids
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        mock_gene_names_nonexistent, mock_digested_plasmids_empty, mock_repair_templates)
    assert len(results) == 0, "Expected no results for empty list of digested plasmids."

def test_assemble_multiple_nonexistent_digested_plasmids():
    # Test with digested plasmids that do not match any gene names
    results = assemble_multiple_plasmids_with_repair_templates_for_deletion(
        mock_gene_names_nonexistent, mock_digested_plasmids_nonexistent, mock_repair_templates)
    assert len(results) == 0, "Expected no results for digested plasmids that do not match any gene names."

