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

import os
from typing import List, Dict, Any
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord
import pandas as pd
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import json
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pydna.dseqrecord import Dseqrecord
from datetime import datetime



def list_of_objects_in_a_dir(dir_path:str): 
    """List all objects in a directory.
    
    Parameters
    ----------
    dir_path : str
        The path to the directory.
        
    Returns
    -------
    list
        A list of the objects in the directory.
    """
    # list to store files
    list_of_files = []

    # Iterate directory
    for path in os.listdir(dir_path):
        # check if current path is a file
        if os.path.isfile(os.path.join(dir_path, path)):
            list_of_files.append(path)
            
    return list_of_files




# TODO fix the naming here. merge with the make_primer_records funciton in primer_generation.py
def dataframe_to_seqrecords(df: pd.DataFrame) -> List[Dseqrecord]:
    """
    Convert sgRNA sequences in a DataFrame to SeqRecord objects.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns ['locus_tag', 'sgrna_loc', 'sgrna']

    Returns
    -------
    List[Dseqrecord]
        List of Dseqrecord objects representing the sgRNA sequences.

    Examples
    --------
    >>> data = {
    ...     'locus_tag': ['tag1', 'tag2'],
    ...     'sgrna_loc': ['loc1', 'loc2'],
    ...     'sgrna': ['ATCG', 'CGTA']
    ... }
    >>> df = pd.DataFrame(data)
    >>> seq_records = dataframe_to_seqrecords(df)
    >>> for record in seq_records:
    >>>     print(record)
    """
    seq_records: List[Dseqrecord] = []

    for index, row in df.iterrows():
        sequence: Seq = Seq(row['sgrna'])
        record_id: str = f"{row['locus_tag']}_{row['sgrna_loc']}"
        record: Dseqrecord = Dseqrecord(SeqRecord(sequence, name=record_id, id=record_id, description=""))
        record.features.append(SeqFeature(FeatureLocation(0, len(sequence)),
                                            type="sgRNA",
                                            qualifiers={"label": f"{row['locus_tag']}_{row['sgrna_loc']}"}))
        seq_records.append(record)

    return seq_records


def create_primer_df_from_dict(records: List[Dict[str, Any]]) -> pd.DataFrame:
    """
    Create a DataFrame from a list of records containing primer and gene information,
    including both upstream and downstream primers, calculating Ta using primer_ta_neb.

    Parameters
    ----------
    records : List[Dict[str, Any]]
        List of dictionaries containing gene, primer, and sequence information.

    Returns
    -------
    pd.DataFrame
        DataFrame with specified columns, representing primer and gene information.

    Examples
    --------
    >>> records = [
    ...     {
    ...         'gene_name': 'gene1',
    ...         'up_forwar_p_anneal': 'ATCG',
    ...         'up_reverse_p_anneal': 'CGTA',
    ...         'tm_up_forwar_p': 60.0,
    ...         'tm_up_reverse_p': 60.0,
    ...         'ta_up': 55.0,
    ...         'up_forwar_primer_str': 'ATCGATCG',
    ...         'up_reverse_primer_str': 'CGTACGTA',
    ...         'up_forwar_p_name': 'up_fwd1',
    ...         'up_reverse_p_name': 'up_rev1',
    ...         'dw_forwar_p_anneal': 'GCTA',
    ...         'dw_reverse_p_anneal': 'TAGC',
    ...         'tm_dw_forwar_p': 58.0,
    ...         'tm_dw_reverse_p': 58.0,
    ...         'ta_dw': 53.0,
    ...         'dw_forwar_primer_str': 'GCTAGCTA',
    ...         'dw_reverse_primer_str': 'TAGCTAGC',
    ...         'dw_forwar_p_name': 'dw_fwd1',
    ...         'dw_reverse_p_name': 'dw_rev1'
    ...     }
    ... ]
    >>> df = create_primer_df_from_dict(records)
    >>> print(df)
    """
    rows = []

    for rec in records:
        template = rec['gene_name']

        f_primer_anneal_up = str(rec['up_forwar_p_anneal'])
        r_primer_anneal_up = str(rec['up_reverse_p_anneal'])
        f_tm_up = rec['tm_up_forwar_p']
        r_tm_up = rec['tm_up_reverse_p']
        ta_up = rec['ta_up']
        f_primer_seq_up = rec['up_forwar_primer_str']
        r_primer_seq_up = rec['up_reverse_primer_str']
        f_primer_name_up = rec['up_forwar_p_name']
        r_primer_name_up = rec['up_reverse_p_name']

        f_primer_anneal_dw = str(rec['dw_forwar_p_anneal'])
        r_primer_anneal_dw = str(rec['dw_reverse_p_anneal'])
        f_tm_dw = rec['tm_dw_forwar_p']
        r_tm_dw = rec['tm_dw_reverse_p']
        ta_dw = rec['ta_dw']
        f_primer_seq_dw = rec['dw_forwar_primer_str']
        r_primer_seq_dw = rec['dw_reverse_primer_str']
        f_primer_name_dw = rec['dw_forwar_p_name']
        r_primer_name_dw = rec['dw_reverse_p_name']

        rows.append([
            template, 'upstream', f_primer_anneal_up, r_primer_anneal_up, f_tm_up, r_tm_up, ta_up,
            f_primer_seq_up, r_primer_seq_up, f_primer_name_up, r_primer_name_up
        ])

        rows.append([
            template, 'downstream', f_primer_anneal_dw, r_primer_anneal_dw, f_tm_dw, r_tm_dw, ta_dw,
            f_primer_seq_dw, r_primer_seq_dw, f_primer_name_dw, r_primer_name_dw
        ])

    df = pd.DataFrame(rows, columns=[
        'template', 'direction', 'f_primer_anneal(5-3)', 'r_primer_anneal(5-3)',
        'f_tm', 'r_tm', 'ta',
        'f_primer_sequences(5-3)', 'r_primer_sequences(5-3)',
        'f_primer_name', 'r_primer_name'
    ])

    return df



# TODO might delete function
def make_antismash_df_to_dseq(df):
    """Convert the sequences in a dataframe to pydnas's Dseqrecord objects. 
    
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing the sequences to be converted. Must contain two columns: 'Sequence' and 'ID' 
        where the sequences and their respective names are located.
    
    Returns
    -------
    sgRNA_Dseqrecords : List[pydna.Dseqrecord]
        List of Dseqrecord objects with the sequences and their respective names as attributes.
    """
    sgRNA_list = df['Sequence'].to_list()
    sgRNA_names = df['ID'].to_list()
    
    sgRNA_Dseqrecords = []
    for i in range(len(sgRNA_list)): 
        sgRNA_Dseqrecords.append(Dseqrecord(sgRNA_list[i], name = sgRNA_names[i]))
    return sgRNA_Dseqrecords



def get_top_sgRNAs(gRNAs:pd.DataFrame, number_of_sgRNAs:int = 4)->pd.DataFrame:
    ''' Retrieves top sgRNAs from a pd.DataFrame.
        
    Parameters
    ----------
    gRNAs : pd.DataFrame
        should have column called "sgrna" and one with "off_target_n".
    
    
    number_of _sgRNAs : int 
        determines how many sgRNAs you want per gene:

    Returns
    --------
    best_gRNAs : pd.DataFrame
        Top sgRNAs
    '''
    # Get the sgRNAs with fewest off targets 
    best_gRNAs = gRNAs.groupby('locus_tag')['off_target_n'].nsmallest(number_of_sgRNAs)
    
    # Make a dataframe with the sgRNAs
    indexes = []
    for i in best_gRNAs.index: 
        indexes.append(i[1])

    best_gRNAs = gRNAs.loc[indexes].reset_index()
    
    return best_gRNAs


polymerase_dict = {
 'Hot Start Taq DNA Polymerase': 'hstaq-0',
 'Hot Start Taq DNA Polymerase 2X Master Mix': 'hstaq-1',
 'Taq DNA Polymerase with Standard Taq Buffer': 'taq-0',
 'Taq DNA Polymerase with Standard (Mg-free) Buffer': 'taq-1',
 'Taq DNA Polymerase with Thermopol Buffer': 'taq-2',
 'Taq DNA Polymerase with Thermopol II (Mg-free) Buffer': 'taq-3',
 'Taq 2X Master Mix': 'taqmaster-0',
 'Taq 5X Master Mix': 'taqmaster-1',
 'QuickLoad Taq 2X Master Mix': 'taqquickload-0',
 'Taq PCR Kit (with Controls)': 'taqpcrkit-0',
 'Crimson Taq DNA Polymerase': 'ctaq-0',
 'Crimson Taq DNA Polymerase with Mg-free Buffer': 'ctaq-1',
 'LongAmp Taq DNA Polymerase': 'lataq-0',
 'LongAmp Taq 2X Master Mix': 'lataq-1',
 'LongAmp Taq PCR Kit': 'lataq-2',
 'Crimson LongAmp Taq DNA Polymerase': 'lataq-3',
 'LongAmp Hot Start Taq DNA Polymerase': 'lahstaq-0',
 'LongAmp Hot Start Taq 2X Master Mix': 'lahstaq-1',
 'EpiMark Hot Start Taq DNA Polymerase': 'epimarkhs-0',
 'Phusion High-Fidelity DNA Polymerase (GC Buffer)': 'phusion-1',
 'Phusion High-Fidelity DNA Polymerase (HF Buffer)': 'phusion-0',
 'Phusion High-Fidelity PCR Kit (GC Buffer)': 'phusion-3',
 'Phusion High-Fidelity PCR Kit (HF Buffer)': 'phusion-2',
 'Phusion High-Fidelity PCR Master Mix (GC Buffer)': 'phusion-5',
 'Phusion High-Fidelity PCR Master Mix (HF Buffer)': 'phusion-4',
 'Phusion Hot Start Flex DNA Polymerase (HF Buffer)': 'phusionflex-0',
 'Phusion Hot Start Flex DNA Polymerase (GC Buffer)': 'phusionflex-1',
 'Phusion Hot Start Flex 2X Master Mix': 'phusionflex-2',
 'Q5 High-Fidelity DNA Polymerase': 'q5-0',
 'Q5 High-Fidelity 2X Master Mix': 'q5-1',
 'Q5 Hot Start High-Fidelity DNA Polymerase': 'q5hs-0',
 'Q5 Hot Start High-Fidelity 2X Master Mix': 'q5hs-1',
 'Q5U Hot Start High-Fidelity DNA Polymerase': 'q5u-0',
 'Q5 Blood Direct 2X Master Mix': 'q5bd-0',
 'default': 'default'
}


def generate_project_directory_structure(
    project_name,
    input_files,
    output_files,
    input_values,
    markdown_file_paths=None,
    create_directories=True
):
    project_dir_structure = {}

    # Define the main project directory
    project_dir = f"./{project_name}"
    
    # Define input and output directories
    inputs_dir = os.path.join(project_dir, "inputs")
    outputs_dir = os.path.join(project_dir, "outputs")
    
    # Create directories if specified
    if create_directories:
        os.makedirs(inputs_dir, exist_ok=True)
        os.makedirs(outputs_dir, exist_ok=True)
    
    # Process and save input files
    for input_file in input_files:
        file_save_path = os.path.join(inputs_dir, input_file['name'])
        if isinstance(input_file['content'], list) and all(isinstance(seq, Dseqrecord) for seq in input_file['content']):
            # Save each Dseqrecord in the list as a separate GenBank file
            for i, seq in enumerate(input_file['content']):
                seq_file_save_path = os.path.join(inputs_dir, f"{input_file['name'].split('.')[0]}_{i}.gb")
                if create_directories:
                    SeqIO.write(seq, seq_file_save_path, "genbank")
                project_dir_structure[seq_file_save_path] = str(seq)
        elif isinstance(input_file['content'], SeqRecord):
            if create_directories:
                SeqIO.write(input_file['content'], file_save_path, "genbank")
            project_dir_structure[file_save_path] = str(input_file['content'])
        else:
            if create_directories:
                with open(file_save_path, 'w') as file:
                    file.write(input_file['content'])
            project_dir_structure[file_save_path] = input_file['content']
    
    # Save the input values as a JSON file
    input_values_path = os.path.join(inputs_dir, "input_values.json")
    if create_directories:
        with open(input_values_path, 'w') as json_file:
            json.dump(input_values, json_file, indent=4)
    project_dir_structure[input_values_path] = input_values
    
    # Process and save output files
    for output_file in output_files:
        file_save_path = os.path.join(outputs_dir, output_file['name'])
        if output_file['name'].endswith(".csv"):
            if isinstance(output_file['content'], pd.DataFrame):
                if create_directories:
                    output_file['content'].to_csv(file_save_path, index=False)
                project_dir_structure[file_save_path] = output_file['content'].to_dict()
            else:
                raise TypeError(f"Expected a DataFrame for {output_file['name']}, but got {type(output_file['content'])}")
        elif output_file['name'].endswith(".gb"):
            if isinstance(output_file['content'], list) and all(isinstance(record, SeqRecord) for record in output_file['content']):
                # Save each SeqRecord in the list separately
                for i, record in enumerate(output_file['content']):
                    record_file_save_path = os.path.join(outputs_dir, f"{record.id}_{i}.gb")
                    if create_directories:
                        SeqIO.write(record, record_file_save_path, "genbank")
                    project_dir_structure[record_file_save_path] = str(record)
            elif isinstance(output_file['content'], SeqRecord):
                if create_directories:
                    SeqIO.write(output_file['content'], file_save_path, "genbank")
                project_dir_structure[file_save_path] = str(output_file['content'])
        else:
            if create_directories:
                with open(file_save_path, 'w') as file:
                    file.write(output_file['content'])
            project_dir_structure[file_save_path] = output_file['content']

    # Process and save markdown files from paths
    if markdown_file_paths:
        for md_file_path in markdown_file_paths:
            md_file_name = os.path.basename(md_file_path)
            md_file_save_path = os.path.join(outputs_dir, md_file_name)
            if create_directories:
                with open(md_file_path, 'r') as file:
                    md_content = file.read()
                with open(md_file_save_path, 'w') as file:
                    file.write(md_content)
            project_dir_structure[md_file_save_path] = md_file_path
    
    if not create_directories:
        print(f"Project structure for '{project_name}':")
        for path, content in project_dir_structure.items():
            print(f"{path}: {str(content)[:100]}")  # Print first 100 characters of the content
    
    return project_dir_structure



import io
import zipfile


class ProjectDirectory:
    def __init__(self, project_name, input_files, output_files, input_values, markdown_file_paths=None):
        self.project_name = project_name
        self.input_files = input_files
        self.output_files = output_files
        self.input_values = input_values
        self.markdown_file_paths = markdown_file_paths
        self.project_dir_structure = {}
        self.zip_buffer = io.BytesIO()

    def create_directory_structure(self, create_directories=True):
        # Define the main project directory
        project_dir = f"{self.project_name}"
        
        # Define input and output directories
        inputs_dir = os.path.join(project_dir, "inputs")
        outputs_dir = os.path.join(project_dir, "outputs")
        
        with zipfile.ZipFile(self.zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            # Create directories if specified
            if create_directories:
                zip_file.writestr(inputs_dir + "/", "")
                zip_file.writestr(outputs_dir + "/", "")
            
            # Process and save input files
            for input_file in self.input_files:
                file_save_path = os.path.join(inputs_dir, input_file['name'])
                if isinstance(input_file['content'], list) and all(isinstance(seq, Dseqrecord) for seq in input_file['content']):
                    for i, seq in enumerate(input_file['content']):
                        seq_file_save_path = os.path.join(inputs_dir, f"{input_file['name'].split('.')[0]}_{i}.gb")
                        with io.StringIO() as seq_buffer:
                            SeqIO.write(seq, seq_buffer, "genbank")
                            zip_file.writestr(seq_file_save_path, seq_buffer.getvalue())
                        self.project_dir_structure[seq_file_save_path] = str(seq)
                elif isinstance(input_file['content'], SeqRecord):
                    with io.StringIO() as seq_buffer:
                        SeqIO.write(input_file['content'], seq_buffer, "genbank")
                        zip_file.writestr(file_save_path, seq_buffer.getvalue())
                    self.project_dir_structure[file_save_path] = str(input_file['content'])
                else:
                    zip_file.writestr(file_save_path, input_file['content'])
                    self.project_dir_structure[file_save_path] = input_file['content']
            
            # Save the input values as a JSON file
            input_values_path = os.path.join(inputs_dir, "input_values.json")
            input_values_content = json.dumps(self.input_values, indent=4)
            zip_file.writestr(input_values_path, input_values_content)
            self.project_dir_structure[input_values_path] = self.input_values
            
            # Process and save output files
            for output_file in self.output_files:
                file_save_path = os.path.join(outputs_dir, output_file['name'])
                if output_file['name'].endswith(".csv"):
                    if isinstance(output_file['content'], pd.DataFrame):
                        csv_content = output_file['content'].to_csv(index=False)
                        zip_file.writestr(file_save_path, csv_content)
                        self.project_dir_structure[file_save_path] = output_file['content'].to_dict()
                    else:
                        raise TypeError(f"Expected a DataFrame for {output_file['name']}, but got {type(output_file['content'])}")
                elif output_file['name'].endswith(".gb"):
                    if isinstance(output_file['content'], list) and all(isinstance(record, SeqRecord) for record in output_file['content']):
                        for i, record in enumerate(output_file['content']):
                            record_file_save_path = os.path.join(outputs_dir, f"{record.id}_{i}.gb")
                            with io.StringIO() as record_buffer:
                                SeqIO.write(record, record_buffer, "genbank")
                                zip_file.writestr(record_file_save_path, record_buffer.getvalue())
                            self.project_dir_structure[record_file_save_path] = str(record)
                    elif isinstance(output_file['content'], SeqRecord):
                        with io.StringIO() as record_buffer:
                            SeqIO.write(output_file['content'], record_buffer, "genbank")
                            zip_file.writestr(file_save_path, record_buffer.getvalue())
                        self.project_dir_structure[file_save_path] = str(output_file['content'])
                else:
                    zip_file.writestr(file_save_path, output_file['content'])
                    self.project_dir_structure[file_save_path] = output_file['content']

            # Process and save markdown files from paths
            if self.markdown_file_paths:
                for md_file_path in self.markdown_file_paths:
                    md_file_name = os.path.basename(md_file_path)
                    md_file_save_path = os.path.join(outputs_dir, md_file_name)
                    with open(md_file_path, 'r') as file:
                        md_content = file.read()
                    zip_file.writestr(md_file_save_path, md_content)
                    self.project_dir_structure[md_file_save_path] = md_content
        
        self.zip_buffer.seek(0)
        return self.zip_buffer.getvalue()

    def get_zip_file(self):
        return self.zip_buffer.getvalue()


def extract_metadata_to_dataframe(seqrecords:Dseqrecord, previous_plasmid:Dseqrecord, integration_list):
    # Initialize lists to store metadata
    names = []
    dates = []
    plasmid_origins = []
    integrations = []
    lengths = []
    
    # Get today's date
    today_date = datetime.today().strftime('%Y-%m-%d')
    
    # Extract metadata from each SeqRecord
    for i, record in enumerate(seqrecords):
        names.append(record.name)
        dates.append(today_date)  # Use today's date
        plasmid_origins.append(previous_plasmid.id)
        integrations.append(integration_list[i] if i < len(integration_list) else 'Unknown')
        lengths.append(len(record.seq))  # Length of the plasmid
    
    # Create a DataFrame from the metadata
    data = {
        'plasmid_name': names,
        'date': dates,
        'original_plasmid': plasmid_origins,
        'integration': integrations,
        'size': lengths
    }
    df = pd.DataFrame(data)
    return df