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
from pydna.primer import Primer
from teemi.build.PCR import primer_ta_neb, primer_tm_neb
from pydna.design import primer_design
from Bio.Seq import Seq
import pandas as pd
from typing import List, Callable
from dataclasses import dataclass, field
from pydna.amplify import pcr
from Bio.Restriction import BsaI
from pydna.tm import tm_default as _tm_default
import collections
from Bio.SeqFeature import SeqFeature, FeatureLocation
import json
import os
import pandas as pd
import requests


@dataclass
class GoldenGateCloning:
    """Data class for making golden gate cloning a brease.

    Parameters:
        sgRNAs: list
            list of sgRNA sequences
        list_of_amplicons: list
            list of amplicon sequences
        target_tm: int, default = 55
            target melting temperature for the primers
        tm_function : function, default = primer_tm_neb
            function to calculate melting temperature
        restriction_overhang_f: str, default = 'GATCGggtctcc'
            overhang to be added to the forward primer
        restriction_overhang_r: str, default = 'GATCAGGTCTCg'
            overhang to be added to the reverse primer
        backbone_overhang_f:str, default = 'cATG'
            overhang to be added to the forward primer
        backbone_overhang_r:str, default = 'cTAG'
            overhang to be added to the reverse primer
        cys4: str, default = 'gTTCACTGCCGTATAGGCAGCTAAGAAA'
            to be apended to the forward primers
    Attributes:
        amplicons : list
            list of amplicon objects
        amplicons_w_f_primer_overhang : list
            list of amplicon objects with forward primer overhangs
        amplicons_w_f_r_primer_overhang : list
            list of amplicon objects with forward and reverse primer overhangs

    Methods:
        simulate_pcrs()
            simulates PCR reactions and returns list of PCR products
        calculate_melting_temperature()
            calculates melting temperature of the primers and returns list of tuples
        make_pcr_df()
            makes a dataframe of PCR amplicon details
        make_primer_df()
            makes a dataframe of primer details
    """

    sgRNAs: list
    list_of_amplicons: list
    target_tm: int = 55
    tm_function: Callable = _tm_default  # u can use your own
    polymerase:str = 'onetaq-3'
    primer_concentration:float = 0.4
    primer_incrementation:int = 1

    restriction_overhang_f: str = "GATCGggtctcc"
    restriction_overhang_r: str = "GATCAGGTCTCg"

    backbone_overhang_f: str = "cATG"
    backbone_overhang_r: str = "cTAG"

    # to be appended to the forward primers
    cys4: str = "gTTCACTGCCGTATAGGCAGCTAAGAAA"

    # for post_init
    amplicons: list = field(init=False)
    amplicons_w_f_primer_overhang: list = field(init=False)
    amplicons_w_f_r_primer_overhang: list = field(init=False)

    def __post_init__(self) -> None:
        # Error handling
        assert isinstance(self.sgRNAs, list), "The input 'sgRNAs' must be a list"
        assert all(isinstance(i, Dseqrecord) for i in self.sgRNAs), "All elements in 'sgRNAs' must be Dseqrecord objects"
        assert isinstance(self.list_of_amplicons, list), "The input 'list_of_amplicons' must be a list"
        assert all(isinstance(i, Dseqrecord) for i in self.list_of_amplicons), "All elements in 'list_of_amplicons' must be Dseqrecord objects"
        assert isinstance(self.target_tm, int), "The input 'target_tm' must be an integer"
        assert isinstance(self.tm_function, Callable), "The input 'tm_function' must be a function"
        assert isinstance(self.restriction_overhang_f, str), "The input 'restriction_overhang_f' must be a string"
        assert isinstance(self.restriction_overhang_r, str), "The input 'restriction_overhang_r' must be a string"
        assert isinstance(self.backbone_overhang_f, str), "The input 'backbone_overhang_f' must be a string"
        assert isinstance(self.backbone_overhang_r, str), "The input 'backbone_overhang_r' must be a string"
        
        # generate the amplicons
        self.amplicons = make_amplicons(
                                        self.list_of_amplicons,
                                        target_tm=self.target_tm,
                                        polymerase=self.polymerase,
                                        primer_concentration=self.primer_concentration,
                                        primer_tm_function=self.tm_function,  # Assuming make_amplicons accepts and uses this
                                        primer_tm_kwargs={'conc': self.primer_concentration, 'prodcode': self.polymerase}  # Pass extra params if needed
                                        )

        # make f primer overhangs
        self.amplicons_w_f_primer_overhang = make_golden_gate_overhangs_forward_primers(
            self.amplicons,
            self.sgRNAs,
            restriction_overhang_f=self.restriction_overhang_f,
            backbone_overhang_f=self.backbone_overhang_f,
            cys4=self.cys4,
        )
        # make f and r primer overhangs
        self.amplicons_w_f_r_primer_overhang = make_golden_gate_overhangs_reverse_primers(
            self.amplicons_w_f_primer_overhang,
            self.sgRNAs,
            backbone_overhang_r=self.backbone_overhang_r,
            restriction_overhang_r=self.restriction_overhang_r,
        )
        # Updating primer names
        primer_counter = self.primer_incrementation # Initialize primer counter
        for amplicon in self.amplicons_w_f_r_primer_overhang:
            amplicon.forward_primer.id  = f"primer_{primer_counter}"  # Assuming a primer_counter is defined somewhere
            amplicon.forward_primer.name = f"{amplicon.forward_primer.id}"
            primer_counter += 1  # Increment the counter for the next primer
            amplicon.reverse_primer.id = f"primer_{primer_counter}"
            amplicon.reverse_primer.name = f"{amplicon.reverse_primer.id}"
            primer_counter += 1

    def simulate_pcrs(self):
        '''Simulates PCR reactions and returns list of PCR products'''
        # Simulation 
        pcr_products = []
        for amplicon in self.amplicons_w_f_r_primer_overhang:
            pcr_product = pcr(
                amplicon.forward_primer.seq.upper(),
                amplicon.reverse_primer.seq.upper(),
                amplicon.template,
            )
            pcr_products.append(pcr_product)

        #Adding features
        for i in range(len(pcr_products)):
            pcr_products[i].name = self.sgRNAs[i].name+'_pcr_products'
            feature = SeqFeature(FeatureLocation(0, len(pcr_products[i])), type="misc_feature", qualifiers={"label": self.sgRNAs[i].name})
            pcr_products[i].features.append(feature)


        return pcr_products
    
    def calculate_melting_temperature(self):
        tm_list = []
        for amplicon in self.amplicons_w_f_r_primer_overhang:
            forward_tm = self.tm_function(str(amplicon.forward_primer.footprint), conc=self.primer_concentration, prodcode=self.polymerase)
            reverse_tm = self.tm_function(str(amplicon.reverse_primer.footprint), conc=self.primer_concentration, prodcode=self.polymerase)
            tm_list.append((forward_tm, reverse_tm))
        return tm_list
    
    def generate_primer_dataframe(self):
        # Initialize an empty list to store dictionaries that represent each row of the DataFrame
        list_of_dicts = []
        
        # Iterate over each amplicon to extract and compute the necessary details
        for i, amplicon in enumerate(self.amplicons_w_f_r_primer_overhang):
            # Compute melting temperatures (Tm) for forward and reverse primers
            f_tm = self.tm_function(str(amplicon.forward_primer.footprint), conc=self.primer_concentration, prodcode=self.polymerase)
            r_tm = self.tm_function(str(amplicon.reverse_primer.footprint), conc=self.primer_concentration, prodcode=self.polymerase)
            
            # Compute annealing temperature (Ta), assuming primer_ta_neb function exists and computes Ta
            ta = primer_ta_neb(
                str(amplicon.forward_primer.footprint),
                str(amplicon.reverse_primer.footprint),
                conc=self.primer_concentration, prodcode=self.polymerase
            )
            
            # Construct the dictionary representing the current amplicon's primer details
            record = {
                "Locus Tag": amplicon.template.name,
                "f_primer_name": amplicon.forward_primer.name,
                "r_primer_name": amplicon.reverse_primer.name,
                "f_primer_sequences(5-3)": str(amplicon.forward_primer.seq).upper(),
                "r_primer_sequences(5-3)": str(amplicon.reverse_primer.seq).upper(),
                "f_tm": f_tm,
                "r_tm": r_tm,
                "ta": ta,
            }
            
            # Append the record to the list
            list_of_dicts.append(record)
        
        # Convert the list of dictionaries to a DataFrame
        df = pd.DataFrame(list_of_dicts, columns=["Locus Tag", "f_primer_name", "r_primer_name", "f_primer_sequences(5-3)", "r_primer_sequences(5-3)", "f_tm", "r_tm", "ta"])
        
        return df



    def make_pcr_df(self):
        list_of_dicts = []
        for i, amplicon in enumerate(self.amplicons_w_f_r_primer_overhang):
            record = {
                "amplicon_name": f'amplicon_{i}_({amplicon.name})',
                "f_primer_id": str(amplicon.forward_primer.name),
                "r_primer_id": str(amplicon.reverse_primer.name),
                "template_for_amplification": amplicon.template.name,
                "f_primer_tm": self.tm_function(str(amplicon.forward_primer.footprint), conc=self.primer_concentration, prodcode=self.polymerase),
                "r_primer_tm": self.tm_function(str(amplicon.reverse_primer.footprint), conc=self.primer_concentration, prodcode=self.polymerase),
                "ta": primer_ta_neb(
                    str(amplicon.forward_primer.footprint ),
                    str(amplicon.reverse_primer.footprint),
                    conc=self.primer_concentration, prodcode=self.polymerase 
                ),
            }
            list_of_dicts.append(record)
        df = pd.DataFrame(list_of_dicts)

        return df

    def make_primer_df(self) -> pd.DataFrame:
        """Makes a dataframe of primer details."""
        list_of_dicts = []
        seen_sequences = set()  # Set to keep track of seen primer sequences

        for amplicon in self.amplicons_w_f_r_primer_overhang:
            # Record for forward primer
            f_sequence = str(amplicon.forward_primer.seq)
            if f_sequence not in seen_sequences:
                record_f = {
                    "id":  f'{str(amplicon.forward_primer.name)}',
                    "sequence": str(amplicon.forward_primer.seq).upper(),
                    "len_primer": len(amplicon.forward_primer.seq),   
                    "annealing_tm": self.tm_function(str(amplicon.forward_primer.footprint), conc=self.primer_concentration, prodcode=self.polymerase),        
                    "primer_type": "Forward",
                    "annealing_sequence": str(amplicon.forward_primer.footprint),
                }
                list_of_dicts.append(record_f)
            
            # Check and add reverse primer if not seen
            r_sequence = str(amplicon.reverse_primer.seq)
            if r_sequence not in seen_sequences:
                record_r = {
                    "id":  f'{str(amplicon.reverse_primer.name)}',
                    "sequence": str(amplicon.reverse_primer.seq).upper(),
                    "len_primer": len(amplicon.reverse_primer.seq),
                    "annealing_tm": self.tm_function(str(amplicon.reverse_primer.footprint), conc=self.primer_concentration, prodcode=self.polymerase),
                    "primer_type": "Reverse",
                    "annealing_sequence": str(amplicon.reverse_primer.footprint),
                }
                list_of_dicts.append(record_r)

        df = pd.DataFrame(list_of_dicts)
        return df


def make_golden_gate_overhangs_forward_primers(
    list_of_amplicons: List[Dseqrecord],
    sgRNA_list: List[Dseqrecord],
    restriction_overhang_f: str = "GATCGggtctcc",
    backbone_overhang_f: str = "cATG",
    cys4: str = "gTTCACTGCCGTATAGGCAGCTAAGAAA",
) -> List[Dseqrecord]:
    """Generates forward primer overhangs to accomodate GoldenGate cloning.

    Parameters
    ----------
    list_of_amplicons: list
        list of pydna.Dseqrecords
    sgRNA_list: list
        list of pydna.Dseqrecords
    - restriction_overhang_f: str
        representing the restriction overhang sequence for the forward primer (default='GATCGggtctcc')
    - backbone_overhang_f: str
        representing the backbone overhang sequence for the forward primer (default='cATG')
    - cys4: str
        representing the Cys4 sequence (default='gTTCACTGCCGTATAGGCAGCTAAGAAA')

    Returns
    -------
    list_of_amplicons: list
        lsit of pydna.amplicon objects with modified forward primer sequences

    """
    # First overhang
    list_of_amplicons[0].forward_primer.seq = (
        restriction_overhang_f
        + backbone_overhang_f
        + cys4
        + sgRNA_list[0].seq.watson
        + list_of_amplicons[0].forward_primer.seq
    )

    # the rest - we start at index 1.
    for i in range(1, len(list_of_amplicons)):
        list_of_amplicons[i].forward_primer.seq = (
            restriction_overhang_f
            + sgRNA_list[i].seq.watson
            + list_of_amplicons[i].forward_primer.seq
        )

    return list_of_amplicons


def make_golden_gate_overhangs_reverse_primers(
    list_of_amplicons: List[Dseqrecord],
    sgRNA_list: List[Dseqrecord],
    backbone_overhang_r: str = "cTAG",
    restriction_overhang_r: str = "gATCAGGTCTCG",
) -> List[Dseqrecord]:
    """Generates reverse primer overhangs to accommodate Golden Gate cloning.

    Parameters
    ----------
    list_of_amplicons: list
        list of pydna.Dseqrecords
    sgRNA_list: list
        list of pydna.Dseqrecords

    backbone_overhang_r: str
        representing the backbone overhang sequence for the reverse primer (default='cTAG')
    restriction_overhang_r: str
        representing the restriction overhang sequence for the reverse primer (default='gATCAGGTCTCG')

    Returns
    -------
    list_of_amplicons: list
        list of pydna.amplicon objects with modified reverse primer sequences
    """
    golden_gate_overhangs = []
    primer_name_mapping = {}  # Dictionary to track primer sequences and their names

    # Find overhangs - last is fixed
    for i in range(1, len(sgRNA_list)):
        golden_gate_overhangs.append(sgRNA_list[i].seq.watson[0:4])
    # last amplicon has fixed overhang
    golden_gate_overhangs.append(backbone_overhang_r)

    # We add restriction site, overhangs, primer seq
    for i, amplicon in enumerate(list_of_amplicons):
        new_seq = (
            restriction_overhang_r
            + Seq(golden_gate_overhangs[i]).reverse_complement()
            + amplicon.reverse_primer.seq
        )   
        if new_seq in primer_name_mapping:
            # If sequence has been seen, use the existing primer name
            amplicon.reverse_primer.name = primer_name_mapping[new_seq]
        else:
            # If sequence hasn't been seen, add it to the dictionary
            primer_name_mapping[new_seq] = amplicon.reverse_primer.name

        amplicon.reverse_primer.seq = new_seq
        amplicon.name = sgRNA_list[i].name 

    return list_of_amplicons



def make_amplicons(
    list_of_amplicons: list, target_tm=58, limit=10,
    primer_concentration: float = 0.4, polymerase: str = 'onetaq-3',
    primer_tm_function=None, primer_tm_kwargs=None):
    """Generates amplicons with designed primers, allowing customization of primer TM calculation,
    and reuses primers when possible, including checking reverse primers."""

    if primer_tm_function is None:
        primer_tm_function = primer_tm_neb
    if primer_tm_kwargs is None:
        primer_tm_kwargs = {'conc': primer_concentration, 'prodcode': polymerase}

    amplicons = []
    previous_forward_primer = None
    previous_reverse_primer = None

    for amplicon_template in list_of_amplicons:
        # Check if both the previous forward and reverse primer sequences can be reused
        can_reuse_forward = previous_forward_primer and str(amplicon_template.seq).startswith(str(previous_forward_primer.seq))
        if previous_reverse_primer:
            previous_reverse_primer_revcomp = str(Seq(previous_reverse_primer.seq).reverse_complement())
            can_reuse_reverse = str(amplicon_template.seq).endswith(str(previous_reverse_primer_revcomp))
        else:
            can_reuse_reverse = False
        if can_reuse_forward and can_reuse_reverse:
            # If reusable, perform PCR with the previous primers and current template
            pcr_product = pcr(
                previous_forward_primer.seq.upper(),
                previous_reverse_primer.seq.upper(),
                amplicon_template,
            )
            # Assuming pcr_product mimics the structure of amplicon_template, adjust as needed
            pcr_product.name = amplicon_template.name + "_PCR"
            pcr_product.id = amplicon_template.id  # Adjust according to your ID logic

            amplicons.append(pcr_product)  # Append PCR product as a new amplicon
        else:
            # Design new primers if not reusable
            amplicon = primer_design(
                amplicon_template,
                target_tm=target_tm,
                limit=limit,
                tm_func=primer_tm_function,
                **primer_tm_kwargs
            )
            amplicon.name = amplicon_template.name + "_amplicon"
            amplicon.id = amplicon_template.id

            # Update the previous primers for the next iteration check
            previous_forward_primer = amplicon.forward_primer
            previous_reverse_primer = amplicon.reverse_primer

            amplicons.append(amplicon)

    return amplicons




def digest_amplicons_w_BsaI(list_of_amplicons: List[Dseqrecord]) -> List[Dseqrecord]:
    """
    Digests a list of amplicons with BsaI and selects the largest fragment from each digestion.

    Parameters
    ----------
    list_of_amplicons : List[Dseqrecord]
        List of pydna.Dseqrecords to be digested.

    Returns
    -------
    List[Dseqrecord]
        List of the largest fragments from each digested amplicon.
    """
    list_of_digested_amplicons = []

    for amplicon in list_of_amplicons:
        # Perform the cut with BsaI
        fragments = amplicon.cut(BsaI)
        
        # Select the largest fragment
        largest_fragment_seq = max(fragments, key=len)
        
        # Convert the largest fragment to a Dseqrecord if it's not already one
        if isinstance(largest_fragment_seq, Dseqrecord):
            largest_fragment = Dseqrecord(largest_fragment_seq)
        else:
            largest_fragment = largest_fragment_seq

        list_of_digested_amplicons.append(largest_fragment)

    return list_of_digested_amplicons




def create_overhang_dataframe(list_of_amplicons):
    # Digest the amplicons 
    digest = digest_amplicons_w_BsaI(list_of_amplicons)
    five_overhangs = []
    three_overhangs = []
    amplicon_names = []

    for i, d in enumerate(digest):
        five_overhangs.append(str(d.seq[0:4]))
        three_overhangs.append(str(d.seq[-4:]))
        amplicon_names.append(f"Amplicon_{list_of_amplicons[i].name}")
    
    # Check for duplicates in 5' and 3' overhangs
    duplicates_5 = [item for item, count in collections.Counter(five_overhangs).items() if count > 1]
    duplicates_3 = [item for item, count in collections.Counter(three_overhangs).items() if count > 1]
    
    # Determine duplicate status for each overhang
    duplicate_status_5 = ['Yes' if overhang in duplicates_5 else 'No' for overhang in five_overhangs]
    duplicate_status_3 = ['Yes' if overhang in duplicates_3 else 'No' for overhang in three_overhangs]
    
    # Create the DataFrame
    df = pd.DataFrame({
        'Amplicon Name': amplicon_names,
        '5\' Overhang': five_overhangs,
        '5\' Duplicate': duplicate_status_5,
        '3\' Overhang': three_overhangs,
        '3\' Duplicate': duplicate_status_3
    })
    
    return df