
from typing import List
import pandas as pd
from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.amplify import pcr
from pydna.assembly import Assembly
from Bio.Restriction import  BstBI, NdeI


def generate_cas3_primers(spacer_table: pd.DataFrame) -> pd.DataFrame:
    """
    Generate forward and reverse primers for the given spacer sequences.

    Parameters
    ----------
    spacer_table : pd.DataFrame
        DataFrame containing spacer sequences in the 'sgrna' column.

    Returns
    -------
    pd.DataFrame
        DataFrame with additional columns for forward and reverse primers.
    """
    spacer_table = spacer_table.copy()  # Create a copy to avoid the SettingWithCopyWarning
    spacer_table.loc[:, "Fwd Primer"] = spacer_table["sgrna"] + "GTCGCcCggCaaAaccGg"
    spacer_table.loc[:, "Rev Primer"] = spacer_table["sgrna"].map(lambda x: str(Seq(x).reverse_complement())) + "GTTTCAATCCACGCGCCCGT"
    return spacer_table


def cas3_plasmid_pcrs(plasmid: Dseqrecord, filtered_df: pd.DataFrame, universal_fwd_seq: str = "GAGCTCATAAGTTCCTATTCCGAAG", universal_rev_seq: str = "aagaagtgggtgtcggacgc") -> List[List[Dseqrecord]]:
    """
    Build a plasmid using the provided primers and filtered DataFrame of spacer sequences.

    Parameters
    ----------
    plasmid : Dseqrecord
        The plasmid sequence to be used in the PCR reaction.
    filtered_df : pd.DataFrame
        DataFrame containing the filtered spacer sequences.
    universal_fwd_seq : str, optional
        Sequence of the universal forward primer, by default "GAGCTCATAAGTTCCTATTCCGAAG".
    universal_rev_seq : str, optional
        Sequence of the universal reverse primer, by default "aagaagtgggtgtcggacgc".

    Returns
    -------
    List[List[Dseqrecord]]
        A list of lists of PCR amplicons generated from the reactions.

    Examples
    --------
    >>> plasmid = Dseqrecord("ATGC...")
    >>> filtered_df = pd.DataFrame({
    ...     'sgrna': ['ATCG', 'GCTA'],
    ...     'Fwd Primer': ['ATCGGTCGCcCggCaaAaccGg', 'GCTAGTCGCcCggCaaAaccGg'],
    ...     'Rev Primer': ['CGATGTTTCAATCCACGCGCCCGT', 'TAGCGTTTCAATCCACGCGCCCGT']
    ... })
    >>> amplicons = cas3_plasmid_pcrs(plasmid, filtered_df)
    >>> for amplicon_pair in amplicons:
    >>>     print(amplicon_pair[0].name, amplicon_pair[1].name)
    """
    # Generate primers for the filtered DataFrame
    filtered_df = generate_cas3_primers(filtered_df)
    
    amplicons_list: List[List[Dseqrecord]] = []
    
    # Define the universal forward and reverse primers
    universal_fwd = Primer(universal_fwd_seq)
    universal_rev = Primer(universal_rev_seq)
    
    # Iterate over each row in the filtered DataFrame
    for index, row in filtered_df.iterrows():
        fwd_primer_seq = row["Fwd Primer"]
        rev_primer_seq = row["Rev Primer"]
        
        # Create Primer objects for the forward and reverse primers from the row
        fwd_primer = Primer(fwd_primer_seq)
        rev_primer = Primer(rev_primer_seq)
        
        # Perform PCR reactions on the plasmid
        pcr1 = pcr(universal_fwd, rev_primer, plasmid)
        pcr2 = pcr(fwd_primer, universal_rev, plasmid)
        
        # Add the amplicons to the list
        amplicons_list.append([pcr1, pcr2])
    
    return amplicons_list


def assemble_cas3_plasmids(plasmid: Dseqrecord, amplicons: list) -> list:
    """
    Assemble new plasmids using the provided plasmid and PCR amplicons.

    Parameters
    ----------
    plasmid : Dseqrecord
        The original plasmid sequence to be cut and assembled.
    amplicons : list
        A list of lists of PCR amplicons generated from the reactions.

    Returns
    -------
    list
        A list of assembled sgRNA vectors.
    """
    # Cut the plasmid with restriction enzymes
    linearized_plasmid = sorted(plasmid.cut(BstBI, NdeI), key=lambda x: len(x), reverse=True)[0]
    
    # assemble plamids w. gibson
    sgRNA_vectors = []
    for pcrs in amplicons:
        list_of_parts = [linearized_plasmid] + pcrs
        new_vector = Assembly(tuple(list_of_parts), limit=20)
        sgRNA_vectors.append(new_vector.assemble_circular()[0])
    
    return sgRNA_vectors

