
from typing import List, Dict, Any
import pandas as pd
from typing import List
from pydna.dseqrecord import Dseqrecord
from pydna.amplicon import Amplicon
from pydna.primer import Primer
from pydna.amplify import pcr
from pydna.design import primer_design
from teemi.build.PCR import primer_tm_neb, primer_ta_neb
from Bio.SeqRecord import SeqRecord


def perform_pcr_on_sequences(df: pd.DataFrame, clean_seq: List[Dseqrecord]) -> List[Amplicon]:
    """
    Performs PCR amplification on a list of sequences using primers generated from a DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the primer sequences and templates. Expected to have columns
        'f_primer_sequences(5-3)', 'r_primer_sequences(5-3)', and 'template'.
    clean_seq : List[Dseqrecord]
        List of clean sequences (as Dseqrecord objects) to be amplified.

    Returns
    -------
    List[Amplicon]
        A list containing the Amplicon objects generated from the PCR reactions.

    Notes
    -----
    The function iterates through each sequence in `clean_seq`, performing PCR using the
    corresponding forward and reverse primers from `df`. Each amplicon is named using the ID
    of the input sequence with "_amplicon" appended. The function prints the name and length
    of each amplicon, along with a graphical representation of the amplification. If PCR fails
    (no product), it prints a warning and the sequence ID for which amplification failed.

    Examples
    --------
    >>> df = generate_primer_dataframe(...) # Assume df is already created
    >>> clean_seq = [Dseqrecord(...), Dseqrecord(...)] # List of sequences
    >>> amplicons = perform_pcr_on_sequences(df, clean_seq)
    >>> for amplicon in amplicons:
    >>>     print(amplicon.name)
    """
    f_primers_list: List[Primer] = [Primer(seq, id=template) for seq, template in zip(df['f_primer_sequences(5-3)'], df['template'])]
    r_primers_list: List[Primer] = [Primer(seq, id=template) for seq, template in zip(df['r_primer_sequences(5-3)'], df['template'])]

    list_of_amplicons: List[Amplicon] = []
    no_pcr_product_count: int = 0

    for i, template_seq in enumerate(clean_seq):
        try:
            amplicon: Amplicon = pcr(f_primers_list[i], r_primers_list[i], template_seq)
            amplicon.name = template_seq.id + "_amplicon"
            print(f"{amplicon.name} , Length: {len(amplicon)}")
            print(amplicon.figure() + '\n')
            list_of_amplicons.append(amplicon)
        except ValueError:
            no_pcr_product_count += 1
            print(f"\n######## No PCR product for {template_seq.id}! Count: {no_pcr_product_count} ########")
            print(template_seq)
            print('###############################')

    return list_of_amplicons




def make_amplicons(list_of_amplicons: List[SeqRecord], 
                   target_tm: int = 58, 
                   limit: int = 10, 
                   primer_concentration: float = 0.4, 
                   polymerase: str = 'onetaq-3', 
                   **primer_kwargs: Dict[str, Any]) -> List[Any]:
    """
    Generates pydna.amplicons which contains primers with a target temperature.

    Parameters
    ----------
    list_of_amplicons : list of SeqRecord
        list of pydna.Dseqrecords
    target_tm : int
        representing the target melting temperature for the primers (default=55)
    limit: int
        representing the minimum primer size (default=5)
    primer_concentration : float
        Concentration of primers used in the reaction (default=0.4)
    polymerase : str
        This is the product code that we can use for NEB_primer_calculator (default='onetaq-3')
    **primer_kwargs : dict
        Additional keyword arguments for the primer design function.

    Returns
    -------
    amplicons : list
        list of amplicon objects with designed primer sequences
    """
    # Set default values for primer_kwargs if not provided
    if not primer_kwargs:
        primer_kwargs = {'conc': primer_concentration, 'prodcode': polymerase}

    amplicons = []
    for i in range(len(list_of_amplicons)):
        amplicon = primer_design(
            list_of_amplicons[i],
            target_tm=target_tm,
            limit=limit,
            tm_function=primer_tm_neb,
            **primer_kwargs
        )
        amplicon.name = list_of_amplicons[i].name + "_amplicon"
        amplicon.id = list_of_amplicons[i].id

        amplicons.append(amplicon)

    return amplicons
