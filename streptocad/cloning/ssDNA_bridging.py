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