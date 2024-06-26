import os
from typing import List, Optional
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.amplicon import Amplicon
from pydna.assembly import Assembly
from Bio.Restriction import StuI
import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation


def assemble_and_process_plasmids(clean_plasmid: Dseq, list_of_amplicons: List[Amplicon], enzyme: Optional[object] = StuI, save_plasmids: bool = False, save_path: str = "data_for_tf_activation_project/plasmids/") -> List[Dseqrecord]:
    """
    Digests a vector with a specified restriction enzyme, assembles it with a list of amplicons,
    and processes and names the assembled plasmids. Optionally saves the assembled plasmids to disk.

    Parameters
    ----------
    clean_plasmid : Dseq
        The clean plasmid sequence to be digested and used as a vector.
    list_of_amplicons : List[Amplicon]
        A list of amplicons to be assembled into the digested vector.
    enzyme : object, optional
        The restriction enzyme used for digesting the vector, by default StuI.
    save_plasmids : bool, optional
        Whether to save the assembled plasmids to disk, by default False.
    save_path : str, optional
        The path to save the assembled plasmids if save_plasmids is True, by default "data_for_tf_activation_project/plasmids/".

    Returns
    -------
    List[Dseqrecord]
        A list of Dseqrecord objects representing the final assembled plasmids.

    Examples
    --------
    >>> clean_plasmid = Dseq("ATGC...")
    >>> list_of_amplicons = [Amplicon(...), Amplicon(...)]
    >>> plasmids = assemble_and_process_plasmids(clean_plasmid, list_of_amplicons)
    >>> for plasmid in plasmids:
    >>>     print(plasmid.name)
    """
    digested_vector = Dseqrecord(clean_plasmid, circular=True).cut(enzyme)[0]
    list_of_assemblies = [(digested_vector, genetic_part) for genetic_part in list_of_amplicons]
    assembled_plasmids = [Assembly(parts, limit=20) for parts in list_of_assemblies]
    
    # figure of the assembly
    assembly_results = []  

    final_assembled_plasmids: List[Dseqrecord] = []
    for i, assemblyobj in enumerate(assembled_plasmids):
        try:
            assembly_result = assemblyobj.assemble_circular()[0]
            assembly_results.append(assembly_result)

            my_plasmid = Dseqrecord(assembly_result.seq, circular=True)
            my_plasmid.name = f'p{i+1}_{list_of_amplicons[i].name[:14]}'
            my_plasmid.id = f'p{i+1}_{list_of_amplicons[i].name[:14]}'

            # Copy annotations from digested vector and assembly result
            my_plasmid.features = digested_vector.features[:]
            my_plasmid.features.extend(assembly_result.features)


            final_assembled_plasmids.append(my_plasmid)
            
            if save_plasmids:
                os.makedirs(save_path, exist_ok=True)
                file_path = os.path.join(save_path, f"{my_plasmid.name}.gb")
                my_plasmid.write(file_path)
        except IndexError:
            print(f"No circular assembly found for assembly {i+1}.")

    return final_assembled_plasmids, assembly_results



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


def annotate_plasmid_with_sgrnas(plasmid_record: Dseqrecord, df: pd.DataFrame) -> Dseqrecord:
    """
    Annotate a plasmid Dseqrecord with sgRNA sequences based on information from a DataFrame.

    Parameters
    ----------
    plasmid_record : Dseqrecord
        The Dseqrecord object representing the plasmid to annotate.
    df : pd.DataFrame
        DataFrame with columns ['strain_name', 'locus_tag', 'gene_loc', 'gene_strand', 'sgrna_strand',
                                'sgrna_loc', 'gc', 'pam', 'sgrna', 'sgrna_seed_sequence', 'off_target_count',
                                'editable_cytosines', 'mutations']

    Returns
    -------
    Dseqrecord
        The annotated Dseqrecord object.
    """
    plasmid_seq = str(plasmid_record.seq)
    
    for index, row in df.iterrows():
        sgrna_seq = row['sgrna']
        
        start_pos = plasmid_seq.find(sgrna_seq)
        if start_pos != -1:
            end_pos = start_pos + len(sgrna_seq)
            sgrna_feature = SeqFeature(
                FeatureLocation(start_pos, end_pos),
                type="sgRNA",
                qualifiers={"label": f"{row['locus_tag']}_pos_{row['sgrna_loc']}_sgRNA"}
            )
            plasmid_record.features.append(sgrna_feature)
    
    return plasmid_record