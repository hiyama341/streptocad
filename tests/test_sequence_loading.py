import pytest
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord
from streptocad.sequence_loading.sequence_loading import (
    load_and_process_genome_sequences,
    load_and_process_plasmid,
    load_and_process_gene_sequences,
    process_specified_gene_sequences_from_record,
    validate_range_format,
    convert_to_dict,
    check_and_convert_input,
    annotate_dseqrecord,
)
import os

@pytest.fixture
def genome_path():
    return 'tests/test_files/Streptomyces_coelicolor_A3_chromosome.gb'

@pytest.fixture
def plasmid_path():
    return 'tests/test_files/pCRISPR-Cas9.gbk'

@pytest.fixture
def coelicolor_genbank_record(genome_path):
    record = Dseqrecord(SeqIO.read(genome_path, "genbank"))
    return record

def test_load_and_process_genome_sequences(genome_path):
    sequences = load_and_process_genome_sequences(genome_path)
    assert isinstance(sequences, list)
    assert all(isinstance(seq, Dseqrecord) for seq in sequences)
    assert len(sequences) > 0  # Check that there are sequences loaded

def test_load_and_process_plasmid(plasmid_path):
    plasmid = load_and_process_plasmid(plasmid_path)
    assert isinstance(plasmid, Dseqrecord)
    assert plasmid.circular is True  # Check if the plasmid is marked as circular

def test_load_and_process_gene_sequences(genome_path):
    gene_sequences = load_and_process_gene_sequences(genome_path)
    assert isinstance(gene_sequences, dict)
    assert len(gene_sequences) > 0  # Ensure that gene sequences are extracted

def test_process_specified_gene_sequences_from_record(coelicolor_genbank_record):
    specified_locus_tags = ["SCO5087", "SCO5090"]
    gene_sequences = process_specified_gene_sequences_from_record(coelicolor_genbank_record, specified_locus_tags)
    assert isinstance(gene_sequences, dict)
    assert "SCO5087" in gene_sequences
    assert "SCO5090" in gene_sequences

def test_validate_range_format():
    start, end = validate_range_format("10-50")
    assert start == 10
    assert end == 50

def test_convert_to_dict():
    input_list = ["10-50", "60-100"]
    target_dict, target_keys, is_range = convert_to_dict(input_list)
    assert isinstance(target_dict, list)
    assert target_keys == ["region1", "region2"]
    assert is_range is True

def test_check_and_convert_input_ranges():
    input_list = ["10-50", "60-100"]
    result, keys, is_range = check_and_convert_input(input_list)
    assert isinstance(result, list)
    assert keys == ["region1", "region2"]
    assert is_range is True

def test_check_and_convert_input_genes():
    input_list = ["gene1", "gene2"]
    result, keys, is_range = check_and_convert_input(input_list)
    assert result == input_list
    assert keys == input_list
    assert is_range is False


def test_annotate_dseqrecord(coelicolor_genbank_record):
    input_list = ["10-50", "60-100"]
    target_dict, target_keys, is_range = convert_to_dict(input_list)
    annotated_record = annotate_dseqrecord(coelicolor_genbank_record, target_dict)

    # Check if the specific features 'region1' and 'region2' are added
    feature_labels = [feat.qualifiers.get("locus_tag", [None])[0] for feat in annotated_record.features]

    # Verify that 'region1' and 'region2' are present in the feature labels
    assert "region1" in feature_labels, "Feature 'region1' was not added to the annotated record."
    assert "region2" in feature_labels, "Feature 'region2' was not added to the annotated record."
