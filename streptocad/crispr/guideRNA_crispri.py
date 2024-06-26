
#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd
from typing import Counter
import re
import math
from collections import Counter
from Bio.Seq import Seq
from typing import Callable, Dict, Optional, List, Tuple
from Bio.SeqFeature import SeqFeature
from ..crispr.guideRNAcas3_9_12 import find_off_target_hits, revcomp, parse_genbank_file


class SgRNAargs:
    """
    A class used to store the parameters for sgRNA generation and filtering.

    Attributes
    ----------
    genbank : str
        Path to the genbank file.
    locus_tag : list of str
        List of locus tags (i.e., genes) to be analyzed.
    model_param : str, optional
        Path to .csv file containing model parameters (default is 'path/to/model_param.csv').
    step : list of str, optional
        Pipeline components to execute (default is ['find', 'filter', 'mismatch']).
    downstream : int, optional
        Number of downstream nucleotides to consider when filtering sgRNAs (default is 15).
    off_target_seed : int, optional
        How many PAM-adjacent bases to use as seed-sequence for off-target hits (default is 9).
    pam_remove : list of str, optional
        List of PAM sequences to filter (default is None).
    downstream_remove : list of str, optional
        List of downstream sequences to filter (default is None).
    sgrna_remove : list of str, optional
        List of specific sgRNA targets to filter (default is None).
    gc_upper : float, optional
        GC cutoff filter (maximum-value) (default is 1).
    gc_lower : float, optional
        GC cutoff filter (minimum-value) (default is 0).
    offset_upper : int, optional
        Maximum distance from start codon past which to filter (default is 1000).
    off_target_upper : int, optional
        Maximum number of seed nucleotides allowed for off-target (default is 10).
    cas_type : str
        The type of the CRISPR-Cas system used ('cas9', 'cas12a', or 'cas13').
    """
    def __init__(self, genbank: str, locus_tag: List[str], cas_type: str = 'cas9',
                 step: List[str] = ['find', 'filter'], 
                 downstream: int = 15, 
                 off_target_seed: int = 13, 
                 pam_remove: Optional[List[str]] = None, 
                 downstream_remove: Optional[List[str]] = None, 
                 sgrna_remove: Optional[List[str]] = None, 
                 gc_upper: float = 1, gc_lower: float = 0, 
                 off_target_upper: int = 10, 
                 extension_to_promoter_region = 0, 
                 upstream_tss = 100,
                 dwstream_tss = 100, 
                 target_non_template_strand = True ):
        
        self.genbank = genbank
        self.locus_tag = locus_tag if isinstance(locus_tag, list) else [locus_tag]
        self.step = step
        self.downstream = downstream
        self.off_target_seed = off_target_seed
        self.pam_remove = pam_remove
        self.downstream_remove = downstream_remove
        self.sgrna_remove = sgrna_remove
        self.gc_upper = gc_upper
        self.gc_lower = gc_lower
        self.off_target_upper = off_target_upper
        self.cas_type = cas_type
        self.extension_to_promoter_region = extension_to_promoter_region
        self.upstream_tss = upstream_tss
        self.dwstream_tss = dwstream_tss
        self.target_non_template_strand = target_non_template_strand
        self.strain_name = self.extract_strain_name(genbank)
    
    @staticmethod
    def extract_strain_name(genbank_path: str) -> str:
        """
        Extracts the strain name from the genbank file path.
        
        Parameters
        ----------
        genbank_path : str
            The path to the genbank file.
            
        Returns
        -------
        str
            The extracted strain name.
        """
        # Assuming the file is named in a 'Genus_species_strain_chromosome.gb' format
        # and located directly in the directory without subdirectories in the path
        filename = genbank_path.split('/')[-1]  # Get the filename from the path
        strain_name_parts = filename.split('_')[:-1]  # Remove the 'chromosome.gb' part
        strain_name = '_'.join(strain_name_parts)  # Rejoin the remaining parts
        return strain_name

def filter_crispri_guides(args: SgRNAargs, hitframe: pd.DataFrame) -> pd.DataFrame:
    """
    Filter the sgRNAs based on the specified criteria.

    Parameters
    ----------
    args : SgRNAargs
        An object of class SgRNAargs containing the parameters for sgRNA filtering.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the filtered sgRNAs.
    """

    def exclude_rows_based_on_patterns(frame: pd.DataFrame, column: str, patterns: List[str]) -> pd.DataFrame:
        """Exclude rows from frame where frame[column] contains any of the provided patterns"""
        if not patterns:
            return frame
        mask = frame[column].apply(lambda x: not any(pattern in x for pattern in patterns))
        return frame[mask]

    # Filter based on the range relative to TSS
    filtered_frame = hitframe.loc[
        (hitframe['sgrna_loc'] >= -args.upstream_tss) & 
        (hitframe['sgrna_loc'] <= args.dwstream_tss)
    ]

    # Filter rows where gene_strand and sgrna_strand are different
    if args.target_non_template_strand:
        filtered_frame = filtered_frame[filtered_frame['gene_strand'] != filtered_frame['sgrna_strand']]

    # Apply other filters
    if args.pam_remove:
        filtered_frame = exclude_rows_based_on_patterns(filtered_frame, "pam", args.pam_remove)
    
    if args.sgrna_remove:
        filtered_frame = exclude_rows_based_on_patterns(filtered_frame, "sgrna", args.sgrna_remove)

    if args.downstream_remove:
        filtered_frame = exclude_rows_based_on_patterns(filtered_frame, "downstream", args.downstream_remove)

    filtered_frame = filtered_frame.loc[filtered_frame.gc >= args.gc_lower]
    filtered_frame = filtered_frame.loc[filtered_frame.gc <= args.gc_upper]
    filtered_frame = filtered_frame.loc[filtered_frame.off_target_count <= args.off_target_upper]


    # TODO: targeting non-template strand 
    #filtered_frame = filtered_frame[filtered_frame['gene_strand'] != filtered_frame['sgrna_strand']]

    return filtered_frame



def find_sgrna_hits_cas9_crispri(filepath: str, strain_name: str, locus_tags: List[str], off_target_counter: Counter, off_target_seed: int, revcomp: callable, extension_to_promoter_region) -> pd.DataFrame:
    """
    Parse a genbank file to find sgRNA hits.

    Parameters
    ----------
    filepath : str
        The path to the genbank file to parse.
    locus_tags : List[str]
        List of locus tags to find in the genbank file.
    off_target_counter : Counter
        Counter object containing the frequency of each off-target hit.
    off_target_seed : int
        The length of the off-target seed sequence to match.
    downstream : int
        The number of downstream base pairs to include in the downstream sequence.
    revcomp : callable
        Function to get the reverse complement of a sequence.

    Returns
    -------
    sgrna_df : pd.DataFrame
        A DataFrame of sgRNA hits information.

    """
    # Initialize list to store sgRNA hits
    sgrna_hits = list()

    # Cas9-specific parameters
    pam_pattern = r"(?=CC)"
    protospacer_len = 20
    pam_len = 3

    # Parse the genbank file again to find sgRNAs
    for record in SeqIO.parse(filepath, "gb"):
        if record.features:
            for feature in record.features:
                if feature.type == "CDS" and feature.qualifiers.get("locus_tag", [""])[0] in locus_tags:
                    locus_tag = feature.qualifiers["locus_tag"][0]
                    location = feature.location
                    gene_strand = feature.location.strand

                    if gene_strand == 1:  # Positive strand
                        start = max(feature.location.start - extension_to_promoter_region, 0)  # Ensure start doesn't go below 0
                        end = min(feature.location.end, len(record.seq))  # Ensure end doesn't go beyond the sequence length
                        coding_sequence = str(record.seq[start:end])
                        coding_sequence_revcomp = revcomp(coding_sequence)
                        watson = 1
                        crick = -1

                    else:  # Negative strand
                        start = max(feature.location.start, 0)  # Here, 'start' is effectively downstream, but the logic remains the same
                        end = min(feature.location.end + extension_to_promoter_region, len(record.seq))  # And 'end' is upstream in terms of genomic coordinates
                        coding_sequence_revcomp = str(record.seq[start:end]) 
                        coding_sequence = revcomp(coding_sequence_revcomp)
                        watson = -1
                        crick = 1

                    # Find potential sgRNAs in both the coding sequence and its reverse complement
                    for sequence in [(crick, coding_sequence), (watson, coding_sequence_revcomp)]:

                        for match in re.finditer(pam_pattern, sequence[1]):
                            strand_sgrna = sequence[0]
                            sgrna_pam = str(sequence[1][match.start():(match.start() + pam_len + protospacer_len)])
                
                            # Get reverse complement of the sgRNA and PAM sequence
                            sgrna = revcomp(sgrna_pam)[0:protospacer_len]
                            pam = revcomp(sgrna_pam)[protospacer_len:protospacer_len+pam_len]

                            if not sgrna:
                                print(f"No sgRNA found for locus tag {locus_tag}. Skipping to next locus tag.")
                                continue  # This skips the rest of the current iteration and moves to the next feature
                            if len(sgrna) != protospacer_len or len(pam) != pam_len:  # If either sgRNA or PAM length is incorrect
                                print(f"sgRNA or PAM generated were outside the designated border in {locus_tag}. Skipping to next locus tag.")
                                continue  # Skip the rest of the current iteration

                            # Calculate GC content of the sgRNA
                            gc_content = len([base for base in sgrna if base in ["C", "G"]]) / len(sgrna) if len(sgrna) > 0 else 0

                            # Calculate genomic location of the sgRNA depending on the strand
                            genome_location = (int(location.start)) +1

                            if sequence[0] == 1:
                                position_sgrna = len(sequence[1])-extension_to_promoter_region -match.start()-3

                            if sequence[0] == -1:
                                position_sgrna = match.end() + protospacer_len + 3 -extension_to_promoter_region

                            sgrna_seed = sgrna[(protospacer_len - off_target_seed):protospacer_len]

                            # Get number of off-target hits for the seed sequence
                            off_target_count = (
                                off_target_counter[sgrna_seed] - 1 # if it turns out to be minus 1 it has not found the seed and there is a mistake.
                            )

                            # Store sgRNA hit
                            sgrna_hits.append(
                                (strain_name, #
                                    locus_tag, #
                                    genome_location, # 
                                    gene_strand, #
                                    strand_sgrna,#
                                    position_sgrna,#
                                    gc_content,
                                    pam,
                                    sgrna,
                                    sgrna_seed,
                                    off_target_count,
                                )
                            )

    # Convert sgRNA hits into a DataFrame
    sgrna_df = pd.DataFrame(
        sgrna_hits,
        columns=['strain_name',
            "locus_tag",
            "gene_loc",
            "gene_strand",
            "sgrna_strand",
            "sgrna_loc",
            "gc",
            "pam",
            "sgrna",
            "sgrna_seed_sequence",
            "off_target_count",
        ]
    )
    return sgrna_df


def extract_sgRNAs_for_crispri(args: SgRNAargs) -> Tuple[pd.DataFrame, Counter, pd.DataFrame]:
    """
    Execute all three functions together to extract gene information, off-target hits, 
    and sgRNA hits from a given genbank file.

    Parameters
    ----------
    args : SgRNAargs
        An instance of the SgRNAargs class.

    Returns
    -------
    gene_df : pd.DataFrame
        A DataFrame of gene information containing locus tag, gene name, strand, start, end.
    off_target_counter : Counter
        Counter object containing the frequency of each off-target hit.
    sgrna_df : pd.DataFrame
        A DataFrame of sgRNA hits information containing genbank file path, locus tag, 
        gene name, strand, offset, position offset, GC content, sgrna, PAM, 
        downstream sequence, sgrna_pam_downstream, seed, off-target count.
    """
    # Extract gene information
    sequences, gene_df = parse_genbank_file(args.genbank)


    # If 'find' is in the steps, execute off-target and sgRNA finding
    if "find" in args.step:
        # Find all potential off-target hits
        off_target_counter = find_off_target_hits(sequences, args.off_target_seed, cas_type=args.cas_type)

        if 'cas9' in args.cas_type:
            # Find all potential sgRNA hits
            sgrna_df = find_sgrna_hits_cas9_crispri(args.genbank,args.strain_name, 
                                                    args.locus_tag, 
                                                    off_target_counter, 
                                                    args.off_target_seed, 
                                                    revcomp, 
                                                    extension_to_promoter_region=args.extension_to_promoter_region)

        # Sort sgrna_df by 'off-targets' in ascending order
        sgrna_df.sort_values(by='sgrna_loc', ascending=True, inplace=True)
    
    # Filter guides if 'filter' is in the steps
    if "filter" in args.step:
        sgrna_df = filter_crispri_guides(args, sgrna_df)

    return  sgrna_df



