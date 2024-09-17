import os
from typing import List, Optional, Dict
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.amplicon import Amplicon
from pydna.assembly import Assembly
from Bio.Restriction import StuI
import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Restriction import RestrictionBatch



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
            my_plasmid.name = f'pOEx-KasO_{i+1}'
            my_plasmid.id = f'pOEx-KasO_{i+1}'

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



def check_plasmid_restriction_sites(plasmid: Dseqrecord, enzymes: List[str]) -> Dict[str, int]:
    """
    Check if the plasmid has specific enzyme cutters and how many times each enzyme cuts.

    Parameters
    ----------
    plasmid : Dseqrecord
        The plasmid sequence as a Dseqrecord object.
    enzymes : List[str]
        A list of enzyme names to check for.

    Returns
    -------
    dict
        A dictionary with enzyme names as keys and the number of cutting sites as values.
    """
    # Create a RestrictionBatch with the specified enzymes
    enzyme_batch = RestrictionBatch(enzymes)
    
    # Analyze the plasmid sequence for restriction sites
    analysis = enzyme_batch.search(plasmid.seq)
    
    # Create a dictionary with enzyme names and the number of cutting sites
    results = {str(enzyme): len(sites) for enzyme, sites in analysis.items()}
    
    return results

def determine_workflow_order_for_plasmids(
    sgRNA_plasmids: List[Dseqrecord], 
    repair_template_plasmids: List[Dseqrecord], 
    sgRNA_enzymes: List[str], 
    repair_template_enzymes: List[str]
) -> pd.DataFrame:
    """
    Determine the workflow for each combination of sgRNA and repair template plasmids.

    Parameters
    ----------
    sgRNA_plasmids : List[Dseqrecord]
        List of plasmids with sgRNA integrated.
    repair_template_plasmids : List[Dseqrecord]
        List of plasmids with repair templates integrated.
    sgRNA_enzymes : List[str]
        List of enzyme names to check for in sgRNA plasmids.
    repair_template_enzymes : List[str]
        List of enzyme names to check for in repair template plasmids.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns indicating the sgRNA plasmid, repair template plasmid, and the workflow decision.
    """
    results = []

    for sgRNA_plasmid in sgRNA_plasmids:
        sgRNA_sites = check_plasmid_restriction_sites(sgRNA_plasmid, sgRNA_enzymes)
        sgRNA_base_name = sgRNA_plasmid.name.split("_")[1]
        
        for repair_template_plasmid in repair_template_plasmids:
            repair_template_base_name = repair_template_plasmid.name.split("_")[1]
            
            if sgRNA_base_name == repair_template_base_name:
                repair_template_sites = check_plasmid_restriction_sites(repair_template_plasmid, repair_template_enzymes)

                # Determine the workflow based on the restriction sites
                # This logic assumes the first enzyme in each list is critical for decision making.
                sgRNA_critical_enzyme = sgRNA_enzymes[0]
                repair_template_critical_enzyme = repair_template_enzymes[0]

                if sgRNA_sites.get(sgRNA_critical_enzyme, 0) > 1 and repair_template_sites.get(repair_template_critical_enzyme, 0) > 0:
                    workflow = "Plasmid cannot be used in this workflow"
                elif sgRNA_sites.get(sgRNA_critical_enzyme, 0) > 1:
                    workflow = "Integrate repair templates first"
                else:
                    workflow = "Proceed with sgRNA integration first"

                result_row = {
                    "sgRNA plasmid": sgRNA_plasmid.name,
                    "plasmid_name": repair_template_plasmid.name,
                    "which workflow to proceed with": workflow,
                }

                # Add sgRNA enzyme sites to the row
                for enzyme in sgRNA_enzymes:
                    result_row[f"sgRNA plasmid #{enzyme} sites"] = sgRNA_sites.get(enzyme, 0)

                # Add repair template enzyme sites to the row
                for enzyme in repair_template_enzymes:
                    result_row[f"repair template plasmid #{enzyme} sites"] = repair_template_sites.get(enzyme, 0)

                results.append(result_row)

    df = pd.DataFrame(results)
    return df