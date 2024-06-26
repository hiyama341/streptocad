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

from pydna.dseqrecord import Dseqrecord
from pydna.assembly import Assembly
import pandas as pd
from pydna.design import primer_design
from teemi.build.PCR import primer_tm_neb
from pydna.tm import tm_default
from Bio.SeqFeature import SeqFeature, FeatureLocation



def assemble_plasmids_by_ssDNA_bridging(ssDNA_primers:list, vector:Dseqrecord)->list:
    ''' Assembles plasmids based on homology. 
        
    Parameters
    ----------
    ssDNA_primers : list
        a list of pydna.Dseqrecords
    vector : Dseqrecord 
  

    Returns
    --------
    sgRNA_vectors : list of pydna.Contigs
        A list of sgRNA_vectors with sgRNA incorporated. 
    '''
    
    sgRNA_vectors = []
    for sgRNA in ssDNA_primers: 
        new_vector = Assembly((vector,sgRNA), limit=20) 
        sgRNA_vectors.append(new_vector.assemble_circular()[0])
        
    return sgRNA_vectors


def checking_primer_design(plasmid: Dseqrecord, sgRNA_midpoints:int, offset:int, target_tm: int = 55, tm_func:callable=tm_default):
    """
    Design checking primers based on specified sgRNA midpoints and offset.

    Parameters
    ----------
    plasmid : Dseqrecord
        The plasmid sequence for which the primers are to be designed.
    
    sgRNA_midpoints : int
        The midpoint position of the sgRNA on the plasmid.
    
    offset : int
        The distance from the midpoint to define the start and end points for primer design.
    
    target_tm : int, optional, default=55
        The desired melting temperature (Tm) for the primers.
    
    tm_func : callable, optional, default=tm_default
        The function to calculate the Tm for the primers.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns 'name', 'Sequence', 'Concentration', and 'Purification'. Contains details for both forward and reverse checking primers.

    Notes
    -----
    - The function assumes that the plasmid sequence is circular.
    - The function returns primers for checking if the sgRNA has been inserted at the expected location in the plasmid.

    Examples
    --------
    >>> plasmid = Dseqrecord("ATCGATCGATCG")
    >>> sgRNA_midpoints = 6
    >>> offset = 3
    >>> checking_primer_design(plasmid, sgRNA_midpoints, offset)
       name                   Sequence Concentration Purification
    0  forward_checking_primer  ATCGATC        25nm          STD
    1  reverse_checking_primer  ATCGATC        25nm          STD

    """
    if sgRNA_midpoints[0]+offset >= len(plasmid): 
        start_point = sgRNA_midpoints[0]-offset
        end_point = sgRNA_midpoints[0]-len(plasmid)+ offset
    else: 
        start_point = sgRNA_midpoints[0]-offset
        end_point = sgRNA_midpoints[0]+ offset

    # make the plasmid circular
    plasmid = Dseqrecord(plasmid, circular = True)
    print(plasmid)


    checking_primer_amplicon = primer_design(Dseqrecord(str(plasmid.seq[start_point:end_point]),
                                        name= f"Checking_primer_amplicon"),
                                        target_tm=target_tm, tm_func=tm_func)
    
    forward_primer = str(checking_primer_amplicon.forward_primer.seq)
    reverse_primer = str(checking_primer_amplicon.reverse_primer.seq)

    df = pd.DataFrame([{'Name':'forward_checking_primer', 'Sequence':forward_primer, 'Concentration':'25nm','Purification':'STD'}, 
          {'Name':'reverse_checking_primer', 'Sequence':reverse_primer, 'Concentration':'25nm','Purification':'STD'}])

    return df


def extract_sgRNA_midpoint_positions(dseq_record):
    """
    Extracts the midpoint positions of 'sgRNA' features from a Dseqrecord.

    Parameters
    ----------
    dseq_record : Dseqrecord
        The input Dseqrecord object.

    Returns
    -------
    list of int
        A list of midpoints for 'sgRNA' features.

    Examples
    --------
    >>> from Bio.Seq import Dseq
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.SeqFeature import SeqFeature, FeatureLocation
    >>> record = SeqRecord(Dseq("ACGTACGT", "TGCA"))
    >>> feature = SeqFeature(FeatureLocation(2, 6), type='sgRNA')
    >>> record.features.append(feature)
    >>> extract_sgRNA_midpoint_positions(record)
    [4]
    """
    sgRNA_midpoints = []
    for feature in dseq_record.features:
        if feature.type == 'sgRNA':
            midpoint = (feature.location.start + feature.location.end) // 2
            sgRNA_midpoints.append(midpoint)
    
    return sgRNA_midpoints


def make_ssDNA_oligos(best_gRNAs:pd.DataFrame,
                      upstream_ovh:str = 'CGGTTGGTAGGATCGACGGC',
                      downstream_ovh:str='GTTTTAGAGCTAGAAATAGC' )-> list:
    
    ''' Makes ssDNA_primers with incorporated sgRNA.
    Incorporates these overhangs:
    CGGTTGGTAGGATCGACGGC **-N20-** GTTTTAGAGCTAGAAATAGC
    For more information please visit the excellent paper for more information: 
    "CRISPR–Cas9, CRISPRi and CRISPR-BEST-mediated genetic manipulation in streptomycetes". 
        
    Parameters
    ----------
    best_gRNAs : pd.DataFrame
        should have column called "sgrna" and one with "locus_tag".
    
    
    upstream_ovh : str 
        optional
    downstream_ovh : str
        optional

    Returns
    --------
    sgRNAs_p : list of pydna.Dseqrecord
        A list of ssDNA primers with sgRNA incorporated. 
    
    '''
    sgRNAs_p = []
    for index, row in best_gRNAs.iterrows():
        record = Dseqrecord(upstream_ovh + row['sgrna']+downstream_ovh, id = f"{row['locus_tag']}_oligo{str(index)}", name= f"{row['locus_tag']}_loc_{row['sgrna_loc']}")
        record.name = f"{row['locus_tag']}_loc_{row['sgrna_loc']}"
        # adding sgRNA feature
        record.features.append(SeqFeature(FeatureLocation(20, 40),
                        type = "sgRNA",
                        qualifiers={"label":f"sgRNA_{row['locus_tag']}_loc_{row['sgrna_loc']}"})                 
                                         )
        sgRNAs_p.append(record)
    
    return sgRNAs_p