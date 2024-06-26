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

# Standard library imports
import sys
import os
import io
import zipfile
import base64
import csv

# Third-party imports
import pandas as pd
from Bio import SeqIO
from Bio.Restriction import NcoI, StuI
from pydna.dseqrecord import Dseqrecord
from pydna.tm import tm_default as _tm_default
from teemi.build.PCR import primer_ta_neb, primer_tm_neb
from teemi.design.fetch_sequences import read_genbank_files

import dash
from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
from dash.dash_table.Format import Group
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import quote

# Local module imports
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.cloning.gibson_cloning import (
    assemble_multiple_plasmids_with_repair_templates_for_deletion, 
    find_up_dw_repair_templates,
    update_primer_names
)

from streptocad.primers.primer_generation import primers_to_IDT


# Main function for data callbacks
def register_workflow_5_callbacks(app):

    @app.callback(
        Output('uploaded-integrated-sgrna-filename', 'children'),
        [Input('upload-integrated-sgrna', 'filename')]
    )
    def display_uploaded_sgrna_filename(filename):
        if not filename:
            return "No file selected."

        # Check file extension
        for file in filename:
            if not (file.endswith('.gb') or file.endswith('.gbk')):
                return "Invalid file type. Please upload a .gb/.gbk file."

        return f"Uploaded file: {filename}"

    @app.callback(
        Output('uploaded-annotated-genome-filename', 'children'),
        [Input('upload-annotated-genome', 'filename')]
    )
    def display_uploaded_genome_filename(filename):
        if not filename:
            return "No file selected."

        # Check file extension
        if not (filename.endswith('.gb') or filename.endswith('.gbk')):
            return "Invalid file type. Please upload a GenBank (.gb or .gbk) file."

        return f"Uploaded file: {filename}"

    ############################################################ 
    # CALLBacks for getting primers for repair templates
    ############################################################ 

    def parse_contents(contents):
        parts = contents.split(',')
        if len(parts) != 2:
            raise ValueError(f"Unexpected format for contents: {contents}")
        
        content_type, content_string = parts
        decoded = base64.b64decode(content_string)
        return io.StringIO(decoded.decode('utf-8'))

    @app.callback(
        [
            Output('primer-table-integrated', 'data'),
            Output('pcr-table-integrated', 'data'),
            Output('genbank-files-integrated', 'href'),
            Output('download-primers', 'href'),
            Output('download-pcr-schema', 'href')
        ],
        [
            Input('submit-sgrna-genome', 'n_clicks'),
            Input('upload-integrated-sgrna', 'contents'),
            Input('upload-integrated-sgrna', 'filename'),
            Input('upload-annotated-genome', 'contents'),
            Input('primer-melting-temp-input', 'value'),  
            Input('overlap-length-input', 'value')       
        ]
    )
    def plasmid_construction_for_in_frame_deletion(n, sgrna_data, sgrna_filenames, genome_data, primer_melting_temp, overlap):
        if n is None:
            raise dash.exceptions.PreventUpdate
        
        plasmid_records = []
        # Process uploaded sgrna_data (plasmids)
        if sgrna_data:
            for file_data in sgrna_data:
                sgrna_file = parse_contents(file_data)
                plasmid_records.append(SeqIO.read(sgrna_file, format='gb'))
        
        # Process uploaded genome_data
        if genome_data:
            genome_file_content = parse_contents(genome_data)  # Assuming one genome file is uploaded
            genome_record = SeqIO.read(genome_file_content, "gb")
        else:
            raise dash.exceptions.PreventUpdate("Genome file not uploaded.")

        # Extract unique gene names from the filenames
        genes_of_interest = list(set([file.split('_')[0] for file in sgrna_filenames]))

        # Make repair templates
        repair_templates_data = find_up_dw_repair_templates(genome_record, genes_of_interest, target_tm=primer_melting_temp)

        # Digest the plasmids
        processed_records = [Dseqrecord(record, circular=True).cut(StuI)[0] for record in plasmid_records]

        # Rename them appropriately
        for i in range(len(processed_records)):
            processed_records[i].name = plasmid_records[i].name

        # Assembly 
        assembly_data = assemble_multiple_plasmids_with_repair_templates_for_deletion(genes_of_interest, processed_records, repair_templates_data, overlap=overlap)
        update_primer_names(assembly_data)

        # PCR dataframe
        pcr_df = pd.DataFrame(assembly_data)

        # Primers to IDT  - takes in Dseqrecords
        all_primers = pcr_df[['up_forwar_p', 'up_reverse_p', 'dw_forwar_p', 'dw_reverse_p']].values.flatten().tolist()
        primers_to_order_df = primers_to_IDT(all_primers)
        primers_to_order_df = primers_to_order_df.drop_duplicates(subset=['Sequence'])

        # Contigs for assemblies and gb files
        assembled_contigs = []
        for data in assembly_data: 
            contig_record = data['contig']
            contig_record.id = f"{data['name']}_w_rep"
            assembled_contigs.append(contig_record)

        # Create a zip archive in memory
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
            for contig in assembled_contigs:
                # Convert each contig to GenBank format
                genbank_content = contig.format("genbank")
                zip_file.writestr(f"{contig.id}.gb", genbank_content)

        # Prepare the zip archive for download
        zip_buffer.seek(0)
        zip_data = base64.b64encode(zip_buffer.read()).decode('utf-8')
        genbank_download_link = f"data:application/zip;base64,{zip_data}"

        # Remove contig column for some reason pandas doesnt like the contig
        pcr_df.drop(columns=["up_forwar_p", "up_reverse_p","dw_forwar_p", "dw_reverse_p", 'up_forwar_primer_str','up_reverse_primer_str','dw_forwar_primer_str','dw_reverse_primer_str', 'up_forwar_p_anneal', 'up_reverse_p_anneal','dw_forwar_p_anneal', 'dw_reverse_p_anneal' ], inplace=True)
        pcr_df = pcr_df.astype(str)

        # Encoding primer and PCR tables
        primer_csv_string = primers_to_order_df.to_csv(index=False)
        primer_data_url = f'data:text/csv;base64,{base64.b64encode(primer_csv_string.encode()).decode()}'

        pcr_csv_string = pcr_df.to_csv(index=False)
        pcr_data_url = f'data:text/csv;base64,{base64.b64encode(pcr_csv_string.encode()).decode()}'

        return primers_to_order_df.to_dict('records'), pcr_df.to_dict('records'), genbank_download_link, primer_data_url, pcr_data_url
