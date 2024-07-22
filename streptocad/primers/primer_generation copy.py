import pandas as pd
from typing import List
from Bio.SeqRecord import SeqRecord
from teemi.build.PCR import primer_tm_neb, primer_ta_neb
from pydna.design import primer_design
import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pydna.dseqrecord import Dseqrecord
from Bio import SeqIO


from ..cloning.pcr_simulation import make_amplicons


def generate_primer_dataframe(clean_seq: List[SeqRecord], 
                              melting_temperature: float, 
                              polymerase: str, 
                              primer_concentration: float, 
                              up_homology: str, 
                              dw_homology: str, 
                              primer_number_increment: int) -> pd.DataFrame:
    """
    Generates a DataFrame containing primer sequences, melting temperatures,
    annealing temperature, and incremental primer names, appending specified homology
    sequences to each primer.

    Parameters
    ----------
    clean_seq : list of SeqRecord
        List of clean sequences for which to generate primers.
    melting_temperature : float
        Target melting temperature for the primers.
    polymerase : str
        Polymerase used for primer Tm calculations.
    primer_concentration : float
        Concentration of primers used in the reaction.
    up_homology : str
        Upstream homology sequence to append to forward primers.
    dw_homology : str
        Downstream homology sequence to append to reverse primers.
    primer_number_increment : int
        Starting number for incremental primer naming.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the columns: template, f_primer_anneal(5-3),
        r_primer_anneal(5-3), f_tm, r_tm, ta, f_primer_sequences(5-3),
        r_primer_sequences(5-3), f_primer_name, r_primer_name.

    Examples
    --------
    >>> from Bio.SeqRecord import SeqRecord
    >>> from Bio.Seq import Seq
    >>> clean_seq = [SeqRecord(Seq("ATGC"), id="gene1"), SeqRecord(Seq("GCAT"), id="gene2")]
    >>> df = generate_primer_dataframe(clean_seq, 60, 'Taq', 0.5, 'ATCG', 'CGTA', 1)
    >>> print(df.head())
    """
    amplicons = make_amplicons(clean_seq, target_tm=melting_temperature, limit=15, polymerase=polymerase, primer_concentration=primer_concentration)
    data = []

    for amplicon, gene_id in zip(amplicons, clean_seq):
        forward_primer_seq = str(amplicon.forward_primer.seq)
        reverse_primer_seq = str(amplicon.reverse_primer.seq)
        f_tm = primer_tm_neb(forward_primer_seq, prodcode=polymerase, conc=primer_concentration)
        r_tm = primer_tm_neb(reverse_primer_seq, prodcode=polymerase, conc=primer_concentration)
        ta = primer_ta_neb(forward_primer_seq, reverse_primer_seq, prodcode=polymerase, conc=primer_concentration)
        
        data.append({
            'template': f'{str(gene_id.id)}',
            'f_primer_anneal(5-3)': forward_primer_seq,
            'r_primer_anneal(5-3)': reverse_primer_seq,
            'f_tm': f_tm,
            'r_tm': r_tm,
            'ta': ta,
        })
    
    df = pd.DataFrame(data)

    df['f_primer_sequences(5-3)'] = up_homology + df['f_primer_anneal(5-3)']
    df['r_primer_sequences(5-3)'] = dw_homology + df['r_primer_anneal(5-3)'] 

    f_primer_names = [f'primer_{i}' for i in range(primer_number_increment, primer_number_increment + len(df) * 2, 2)]
    r_primer_names = [f'primer_{i}' for i in range(primer_number_increment + 1, primer_number_increment + 1 + len(df) * 2, 2)]

    df['f_primer_name'] = f_primer_names
    df['r_primer_name'] = r_primer_names

    return df


def create_idt_order_dataframe(df: pd.DataFrame, concentration: str = "25nm", purification: str = "STD") -> pd.DataFrame:
    """
    Creates a DataFrame formatted for IDT (Integrated DNA Technologies) orders,
    containing primer names, sequences, specified concentration, and purification method.
    Allows for customization of the concentration and purification while providing
    default values.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing the primer information. Expected to have columns
        'f_primer_name', 'r_primer_name', 'f_primer_sequences(5-3)', and 'r_primer_sequences(5-3)'.
    concentration : str, optional
        The concentration for all primers, by default "25nm".
    purification : str, optional
        The purification method for all primers, by default "STD".

    Returns
    -------
    pd.DataFrame
        A new DataFrame with columns 'Name', 'Sequences', 'Concentration', and 'Purification',
        formatted for IDT orders, with the specified concentration and purification methods.

    Examples
    --------
    >>> df = pd.DataFrame({
    ...     'f_primer_name': ['f_primer1', 'f_primer2'],
    ...     'r_primer_name': ['r_primer1', 'r_primer2'],
    ...     'f_primer_sequences(5-3)': ['ATCG', 'CGTA'],
    ...     'r_primer_sequences(5-3)': ['GATC', 'TACG']
    ... })
    >>> idt_df = create_idt_order_dataframe(df, concentration="100nm", purification="PAGE")
    >>> print(idt_df)
    """
    primer_names = df['f_primer_name'].tolist() + df['r_primer_name'].tolist()
    primer_sequences = df['f_primer_sequences(5-3)'].tolist() + df['r_primer_sequences(5-3)'].tolist()

    idt_data = {
        "Name": primer_names,
        "Sequence": primer_sequences,
        "Concentration": [concentration] * len(primer_names),
        "Purification": [purification] * len(primer_names)
    }

    idt_df = pd.DataFrame(idt_data)

    return idt_df


# TODO almost exactly like the funciton above
def primers_to_IDT(sgRNAs_p:list, 
                   concentration:str = '25nm',
                   purification:str = 'STD')->pd.DataFrame:
    
    '''"Generates a DataFrame containing primers that 
    are ready for submission to IDT for synthesis,
        
    Parameters
    ----------
    sgRNAs_p : list
        list of pydna.Dseqrecords
    concentration : str
    purification : str

    Returns
    --------
    best_gRNAs : pd.DataFrame
        Top sgRNAs
    
    '''
    
    name_lst = [name.name for name in sgRNAs_p]
    seq_lst = [str(seq.seq) for seq in sgRNAs_p]

    d = {'Name':name_lst,'Sequence':seq_lst, 'Concentration': concentration, 'Purification':purification}
    df = pd.DataFrame(d)
    return df


# TODO rename and merge with the dataframe_to_seqrecords funciton in utils.py
def make_primer_records(filtered_df: pd.DataFrame) -> List[Dseqrecord]:
    """
    Create Dseqrecord objects for the forward and reverse primers with names according to gene and location.

    Parameters
    ----------
    filtered_df : pd.DataFrame
        DataFrame containing the filtered spacer sequences with columns for forward and reverse primers.

    Returns
    -------
    List[Dseqrecord]
        A list of Dseqrecord objects for the forward and reverse primers.

    Examples
    --------
    >>> data = {
    ...     'locus_tag': ['gene1', 'gene2'],
    ...     'sgrna_loc': ['loc1', 'loc2'],
    ...     'Fwd Primer': ['ATCGGTCGCcCggCaaAaccGg', 'GCTAGTCGCcCggCaaAaccGg'],
    ...     'Rev Primer': ['CGATGTTTCAATCCACGCGCCCGT', 'TAGCGTTTCAATCCACGCGCCCGT']
    ... }
    >>> df = pd.DataFrame(data)
    >>> primer_records = make_primer_records(df)
    >>> for record in primer_records:
    >>>     print(record)
    """
    primer_records: List[Dseqrecord] = []
    
    for index, row in filtered_df.iterrows():
        # Create Dseqrecord for forward primer
        fwd_primer_record = Dseqrecord(row['Fwd Primer'], 
                                       id=f"{row['locus_tag']}_fwd_{str(index)}", 
                                       name=f"{row['locus_tag']}_loc_{row['sgrna_loc']}_fwd")
        fwd_primer_record.features.append(SeqFeature(FeatureLocation(0, len(row['Fwd Primer'])),
                                                     type="primer_bind",
                                                     qualifiers={"label": f"Fwd_{row['locus_tag']}_loc_{row['sgrna_loc']}"}))
        
        # Create Dseqrecord for reverse primer
        rev_primer_record = Dseqrecord(row['Rev Primer'], 
                                       id=f"{row['locus_tag']}_rev_{str(index)}", 
                                       name=f"{row['locus_tag']}_loc_{row['sgrna_loc']}_rev")
        rev_primer_record.features.append(SeqFeature(FeatureLocation(0, len(row['Rev Primer'])),
                                                     type="primer_bind",
                                                     qualifiers={"label": f"Rev_{row['locus_tag']}_loc_{row['sgrna_loc']}"}))
        
        primer_records.extend([fwd_primer_record, rev_primer_record])
    
    return primer_records

def checking_primers(genbank_file, locus_tags, flanking_region=500, 
                     target_tm=58, 
                     limit=10, 
                     primer_concentration=0.4, 
                     polymerase='onetaq-3', 
                     **primer_kwargs):
    

    # Parse the GenBank file
    record = SeqIO.read(genbank_file, "genbank")
    

    # Set default values for primer_kwargs if not provided
    if not primer_kwargs:
        primer_kwargs = {'conc': primer_concentration, 'prodcode': polymerase}

    # Prepare a list to store primer information
    primer_info = []

    for locus_tag in locus_tags:
        # Find the feature with the specified locus_tag
        feature = None
        for f in record.features:
            if f.type == "CDS" and "locus_tag" in f.qualifiers and f.qualifiers["locus_tag"][0] == locus_tag:
                feature = f
                break

        if feature is None:
            raise ValueError(f"Locus tag {locus_tag} not found in the GenBank file.")

        # Extract the sequence including flanking regions
        start = max(0, feature.location.start - flanking_region)
        end = min(len(record), feature.location.end + flanking_region)
        target_seq = record.seq[start:end]
        
        # Convert to pydna Dseqrecord
        dseqrecord = Dseqrecord(target_seq)
        
        # Design primers with specified parameters and kwargs
        primers = primer_design(
            dseqrecord, 
            target_tm=target_tm, 
            limit=limit, 
            **primer_kwargs
        )
        
        # Rename the primers
        forward_primer_name = f"{locus_tag}_fwd_checking_primer"
        reverse_primer_name = f"{locus_tag}_rev_checking_primer"
        primers.forward_primer.name = forward_primer_name
        primers.reverse_primer.name = reverse_primer_name
        
        # Collect primer information in the desired format
        primer_info.append({
            "Locus Tag": locus_tag,
            "f_primer_name": forward_primer_name,
            "r_primer_name": reverse_primer_name,
            "f_primer_sequences(5-3)": str(primers.forward_primer.seq),
            "r_primer_sequences(5-3)": str(primers.reverse_primer.seq),
            "f_tm": primer_tm_neb(str(primers.forward_primer.seq), **primer_kwargs),
            "r_tm": primer_tm_neb(str(primers.reverse_primer.seq), **primer_kwargs),
            "ta": primer_ta_neb(str(primers.forward_primer.seq), str(primers.reverse_primer.seq), **primer_kwargs)
        })

    # Create a DataFrame from the collected primer information
    primer_df = pd.DataFrame(primer_info)
    
    return primer_df



def checking_primers1(record, locus_tags, flanking_region=500, 
                     target_tm=65, 
                     limit=10, 
                     primer_concentration=0.4, 
                     polymerase='onetaq-3', 
                     **primer_kwargs):
    # Set default values for primer_kwargs if not provided
    if not primer_kwargs:
        primer_kwargs = {'conc': primer_concentration, 'prodcode': polymerase}

    # Prepare a list to store primer information
    primer_info = []

    for locus_tag in locus_tags:
        # Find the feature with the specified locus_tag
        feature = None
        for f in record.features:
            if f.type == "CDS" and "locus_tag" in f.qualifiers and f.qualifiers["locus_tag"][0] == locus_tag:
                feature = f
                break

        if feature is None:
            raise ValueError(f"Locus tag {locus_tag} not found in the GenBank file.")

        # Extract the sequence including flanking regions
        start = max(0, feature.location.start - flanking_region)
        end = min(len(record), feature.location.end + flanking_region)
        target_seq = record.seq[start:end]
        
        # Convert to pydna Dseqrecord
        dseqrecord = Dseqrecord(target_seq)
        
        # Design primers with specified parameters and kwargs
        primers = primer_design(
            dseqrecord, 
            target_tm=target_tm, 
            limit=limit, 
            **primer_kwargs
        )
        
        # Rename the primers
        forward_primer_name = f"{locus_tag}_fwd_checking_primer"
        reverse_primer_name = f"{locus_tag}_rev_checking_primer"
        primers.forward_primer.name = forward_primer_name
        primers.reverse_primer.name = reverse_primer_name
        
        # Collect primer information in the desired format
        primer_info.append({
            "Locus Tag": locus_tag,
            "f_primer_name": forward_primer_name,
            "r_primer_name": reverse_primer_name,
            "f_primer_sequences(5-3)": str(primers.forward_primer.seq),
            "r_primer_sequences(5-3)": str(primers.reverse_primer.seq),
            "f_tm": primer_tm_neb(str(primers.forward_primer.seq), **primer_kwargs),
            "r_tm": primer_tm_neb(str(primers.reverse_primer.seq), **primer_kwargs),
            "ta": primer_ta_neb(str(primers.forward_primer.seq), str(primers.reverse_primer.seq), **primer_kwargs)
        })

    # Create a DataFrame from the collected primer information
    primer_df = pd.DataFrame(primer_info)
    
    return primer_df


def create_dseqrecords_from_df(primer_df: pd.DataFrame) -> pd.DataFrame:
    """
    Iterates through a DataFrame and converts primer sequences to Dseqrecords with appropriate names.

    Parameters
    ----------
    primer_df : pd.DataFrame
        DataFrame containing primer information with columns:
        'Locus Tag', 'Forward Primer Sequence', 'Forward Primer Tm',
        'Reverse Primer Sequence', 'Reverse Primer Tm', 'Reverse Primer Ta'.

    Returns
    -------
    pd.DataFrame
        DataFrame with Dseqrecord objects for the primer sequences.
    """
    primer_dseqrecords = []
    for _, row in primer_df.iterrows():
        fwd_dseqrecord = Dseqrecord(row['Forward Primer Sequence'])
        fwd_dseqrecord.name = row['Locus Tag'] + "_fwd_checking_primer"
        
        rev_dseqrecord = Dseqrecord(row['Reverse Primer Sequence'])
        rev_dseqrecord.name = row['Locus Tag'] + "_rev_checking_primer"
        
        primer_dseqrecords.append(fwd_dseqrecord)
        primer_dseqrecords.append(rev_dseqrecord)
    
    return primer_dseqrecords