
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from teemi.build.PCR import primer_tm_neb
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO 
from Bio.SeqUtils import MeltingTemp as mt
from pydna.dseqrecord import Dseqrecord
import pandas as pd
from teemi.design.fetch_sequences import read_genbank_files
from pydna.amplicon import Amplicon

from streptocad.cloning.plasmid_processing import (
    assemble_and_process_plasmids,
    assemble_plasmids_by_ssDNA_bridging, 
    annotate_plasmid_with_sgrnas, 
    check_plasmid_restriction_sites, 
)



def test_assemble_and_process_plasmids(): 
    pass


# TODO move to ssDNA bridging 
def test_assemble_plasmids_by_ssDNA_bridging(): 
    pass

def test_annotate_plasmid_with_sgrnas(): 
    pass

def test_check_plasmid_restriction_sites(): 
    pass