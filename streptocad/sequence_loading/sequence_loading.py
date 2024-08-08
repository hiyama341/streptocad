import os
from typing import List, Dict
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord
from teemi.design.fetch_sequences import read_fasta_files, read_genbank_files
from Bio.SeqFeature import SeqFeature, FeatureLocation
import re
import tempfile

def load_and_process_genome_sequences(path_to_file: str) -> List[Dseqrecord]:
    """
    Loads gene sequences from a given file path. The file can be in FASTA or GenBank format.
    Returns a list of Dseqrecord objects representing the sequences.

    Parameters
    ----------
    path_to_file : str
        Path to the file containing gene sequences. Supported formats are FASTA and GenBank (gb or gbk).

    Returns
    -------
    list of Dseqrecord
        A list containing Dseqrecord objects for each gene sequence in the input file.

    Raises
    ------
    ValueError
        If the file format is not supported (not FASTA or GenBank).
    """
    # Check the file extension to determine the format
    file_extension = os.path.splitext(path_to_file)[1].lower()
    if file_extension in ['.fasta', '.fa']:
        # Handle FASTA file
        sequences = read_fasta_files(path_to_file)
    elif file_extension in ['.gb', '.gbk']:
        # Handle GenBank file
        sequences = read_genbank_files(path_to_file
                                       )
    else:
        raise ValueError("Unsupported file format. Please provide a FASTA or GenBank file.")

    # Convert each sequence to a Dseqrecord object
    clean_seq = [Dseqrecord(seq) for seq in sequences]

    return clean_seq

import os
from pydna.dseqrecord import Dseqrecord
from typing import List

def load_and_process_plasmid(path_to_file: str) -> List[Dseqrecord]:

    """
    Loads a plasmid sequence from a given file path. The file can be in FASTA or GenBank format.
    Returns a Dseqrecord object representing the plasmid sequence, marked as circular.

    Parameters
    ----------
    path_to_file : str
        Path to the file containing the plasmid sequence. Supported formats are FASTA and GenBank (gb or gbk).

    Returns
    -------
    Dseqrecord
        A Dseqrecord object representing the plasmid sequence, marked as circular.

    Raises
    ------
    ValueError
        If the file format is not supported (not FASTA or GenBank).
    """
    # Check the file extension to determine the format
    file_extension = os.path.splitext(path_to_file)[1].lower()
    if file_extension in ['.fasta', '.fa']:
        # Handle FASTA file
        sequences = read_fasta_files(path_to_file)
    elif file_extension in ['.gb', '.gbk']:
        # Handle GenBank file
        sequences = read_genbank_files(path_to_file)
    else:
        raise ValueError("Unsupported file format. Please provide a FASTA or GenBank file.")

    # Assuming the file contains only one plasmid sequence, mark it as circular
    if sequences:
        clean_plasmid = Dseqrecord(sequences[0], circular=True)
        clean_plasmid.name = sequences[0].name
        clean_plasmid.id = sequences[0].id
    else:
        raise ValueError("No sequences found in the file.")

    return clean_plasmid



def load_and_process_gene_sequences(path_to_genome: str) -> dict:
    """
    Load and process gene sequences from a GenBank or FASTA file.
    
    Parameters
    ----------
    path_to_genome : str
        Path to the genome file.
        
    Returns
    -------
    dict
        Dictionary mapping locus tags to gene sequences.
    """
    gene_sequences = {}
    
    for record in SeqIO.parse(path_to_genome, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                if locus_tag:
                    gene_seq = feature.extract(record.seq)
                    gene_sequences[locus_tag] = gene_seq
                    
    return gene_sequences



def process_specified_gene_sequences_from_record(seq_record: Dseqrecord, specified_locus_tags: list) -> dict:
    """
    Process gene sequences from a SeqRecord object for specified locus tags.
    
    Parameters
    ----------
    seq_record : SeqRecord
        A SeqRecord object containing genome data.
    specified_locus_tags : list
        A list of locus tags to collect sequences for.
        
    Returns
    -------
    dict
        Dictionary mapping specified locus tags to their gene sequences if found.
    """
    gene_sequences = {}
    
    for feature in seq_record.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
            if locus_tag in specified_locus_tags:
                gene_seq = feature.extract(seq_record.seq)
                gene_sequences[locus_tag] = gene_seq
                
    return gene_sequences


def validate_range_format(item):
    """Validates the range format 'start-end' and ensures 'start' and 'end' are valid."""
    positions = item.split('-')
    
    # Ensure the length of positions is even and contains exactly two numbers
    assert len(positions) == 2, f"Input format error: {item} should contain exactly two numbers separated by a dash."
    
    # Ensure all positions are numbers
    for pos in positions:
        assert pos.isdigit(), f"Input format error: {pos} is not a valid number."
    
    # Convert positions to integers
    start, end = int(positions[0]), int(positions[1])
    
    # Ensure the second number is not smaller than the first number
    assert start <= end, f"Input format error: the second number {end} is smaller than the first number {start} in {item}."
    
    return start, end

def convert_to_dict(input_list):
    """Converts a list of correctly formatted ranges into a list of dictionaries."""
    target_dict = []
    target_keys = []
    for index, item in enumerate(input_list):
        start, end = validate_range_format(item)
        
        # Create the dictionary for this region
        region_key = f'region{index + 1}'#_({start}-{end})'
        region_dict = {region_key: [start, end]}
        
        # Add the dictionary to the list
        target_keys.append(region_key)
        target_dict.append(region_dict)
    
    return target_dict, target_keys, True

def check_and_convert_input(input_list):
    """Checks the input type and converts it appropriately."""
    if all(re.fullmatch(r"\d+-\d+", item) for item in input_list):
        return convert_to_dict(input_list)
    
    elif all(re.fullmatch(r"[A-Za-z0-9_]+", item) for item in input_list):
        return input_list, input_list,  False
    
    else:
        raise AssertionError("Input format error: The input should be either a list of ranges or a list of gene names.")


def annotate_dseqrecord(dseqrecord, target_dict):
    """Annotates a SeqRecord with features specified in target_dict."""
    if not isinstance(dseqrecord, Dseqrecord):
        raise ValueError("dseqrecord must be an instance of SeqRecord.")
    
    print('this is the target dict',target_dict)
    for annotation in target_dict:
        for region, positions in annotation.items():
            print(region, positions)
            start, end = positions
            
            # Ensure start and end are within the bounds of the sequence
            if not (0 <= start < len(dseqrecord.seq)) or not (0 < end <= len(dseqrecord.seq)):
                raise ValueError(f"Positions for region {region} are out of bounds.")
            
            feature_location = FeatureLocation(start, end, strand=1)
            feature = SeqFeature(feature_location, type="CDS", qualifiers={"locus_tag": [region]})
            dseqrecord.features.append(feature)
    
    return dseqrecord