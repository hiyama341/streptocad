import os
from typing import List
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord
from teemi.design.fetch_sequences import read_fasta_files, read_genbank_files


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

    #Assuming the file contains only one plasmid sequence, mark it as circular
    if sequences:
        clean_plasmid = Dseqrecord(sequences[0], circular=True)
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