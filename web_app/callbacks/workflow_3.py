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

from streptocad.cloning.golden_gate_cloning import (
    GoldenGateCloning, 
    create_overhang_dataframe, 
    digest_amplicons_w_BsaI
)

from streptocad.utils import make_antismash_df_to_dseq

# Main function for data callbacks
def register_workflow_3_callbacks(app):

    @app.callback(
        Output('uploaded-genome-filename', 'children'),
        Input('upload-genome-file', 'filename')
    )
    def display_uploaded_genome_filename(genome_filename):
        if not genome_filename:
            return "No file selected."

        # Check file extension
        if not (genome_filename.endswith('.gb') or genome_filename.endswith('.gbk') or genome_filename.endswith('.fasta') or genome_filename.endswith('.fa')):
            return "Invalid file type. Please upload a GenBank (.gb or .gbk) or FASTA (.fasta or .fa) file."

        return f"Uploaded file: {genome_filename}"

    @app.callback(
        Output('uploaded-workflow3-vector-filename', 'children'),
        Input('upload-single-vector', 'filename')
    )
    def display_uploaded_vector_filename(filename):
        if not filename:
            return "No file selected."

        # Check file extension
        if not (filename.endswith('.gb') or filename.endswith('.gbk')):
            return "Invalid file type. Please upload a GenBank (.gb or .gbk) file."

        return f"Uploaded file: {filename}"

    ############################################################ 
    ## CALLBacks for Multiple sgRNA integration
    ############################################################ 

    # Callback function for Multiple sgRNA integration
    @app.callback(
        [
            Output('genbank-file', 'href'),
            Output('primer-table', 'data'),
            Output('primer-table', 'columns'),
            Output('pcr-table', 'data'),
            Output('pcr-table', 'columns'),
            Output('overhang-table', 'data'),
            Output('overhang-table', 'columns'),
            Output('download-primers-link', 'href'),
            Output('download-pcr-link', 'href')
        ],
        Input('submit-button', 'n_clicks'),
        State('upload-genome-file', 'contents'),
        State('upload-genome-file', 'filename'),
        State('upload-single-vector', 'contents'),
        State('upload-single-vector', 'filename'),
        State('genes-to-KO', 'value'),
        State('sgRNA-handle-input', 'value'),
        State('input-tm', 'value'),
        State('gc-upper', 'value'),
        State('gc-lower', 'value'),
        State('off-target-seed', 'value'),
        State('off-target-upper', 'value'),
        State('cas-type', 'value'),
        State('number-of-sgRNAs-per-group', 'value'),
        State('only-stop-codons-checkbox', 'value'),
        State('show-advanced-settings-checkbox', 'value'),
        State('chosen-polymerase', 'value'),
        State('melting-temperature', 'value'),
        State('primer-concentration', 'value'),
        State('primer-number-increment', 'value'),
        State('flanking-region-number', 'value')
    )
    def process_inputs(n_clicks, genome_file_contents, genome_filename, vector_file_contents, vector_filename, genes_to_KO, sgRNA_handle, tm, gc_upper, gc_lower, off_target_seed, off_target_upper, cas_type, number_of_sgRNAs, only_stop_codons, show_advanced, chosen_polymerase, melting_temp, primer_conc, primer_increment, flanking_region):
        if n_clicks is None:
            raise PreventUpdate
        
        # Process the files and inputs here
        # Initialize the outputs
        primer_df = pd.DataFrame()
        pcr_df = pd.DataFrame()
        genbank_data_url = ""
        primer_data_url = ""
        pcr_data_url = ""
        overhang_data = []
        overhang_df = pd.DataFrame()

        # Process genome file
        if genome_file_contents:
            content_type, content_string = genome_file_contents.split(',')
            decoded = base64.b64decode(content_string)
            try:
                if genome_filename.endswith('.gb') or genome_filename.endswith('.gbk'):
                    genome_record = SeqIO.read(io.StringIO(decoded.decode('utf-8')), 'genbank')
                elif genome_filename.endswith('.fasta') or genome_filename.endswith('.fa'):
                    genome_record = SeqIO.read(io.StringIO(decoded.decode('utf-8')), 'fasta')
                else:
                    return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update
            except Exception as e:
                print(e)
                return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        # Process vector file
        if vector_file_contents:
            content_type, content_string = vector_file_contents.split(',')
            decoded = base64.b64decode(content_string)
            try:
                vector_record = SeqIO.read(io.StringIO(decoded.decode('utf-8')), 'genbank')
                vector = Dseqrecord(vector_record, circular=True)
                vector.name = 'pCRISPR-MCBE_Csy4'
            except Exception as e:
                print(e)
                return dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update, dash.no_update

        # Process genes to knock out
        genes_to_KO_list = genes_to_KO.split(',')

        # Initialize the script
        sgRNA_list = make_antismash_df_to_dseq(pd.DataFrame())  # Replace with actual dataframe processing
        sgRNA_handle_cys4_sites = [Dseqrecord('GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTT', name = 'sgRNA_handle_cys4')]*len(sgRNA_list)
        
        golden_gate = GoldenGateCloning(sgRNA_list, sgRNA_handle_cys4_sites, target_tm=tm, tm_function=primer_tm_neb) 
        list_of_amplicons = golden_gate.simulate_pcrs()

        # Overhangs
        overhang_df = create_overhang_dataframe(list_of_amplicons)

        digest = digest_amplicons_w_BsaI(list_of_amplicons)

        from Bio.Restriction import NcoI, NheI
        vector_NcoI_NheI, smaller_frag = sorted(vector.cut(NcoI, NheI), reverse=True)
        rec_vec = vector_NcoI_NheI

        for i in range(len(digest)):
            rec_vec += digest[i]
        rec_vec = rec_vec.looped()

        # Initialize output
        genbank_string_io = io.StringIO()
        SeqIO.write(rec_vec, genbank_string_io, 'genbank')
        genbank_string = genbank_string_io.getvalue()
        
        pcr_df = golden_gate.make_pcr_df()
        primer_df = golden_gate.make_primer_df()
        overhang_data = overhang_df.to_dict('records') if isinstance(overhang_df, pd.DataFrame) else []

        # Convert primer_df and pcr_df to CSV strings
        primer_csv_string = primer_df.to_csv(index=False) if isinstance(primer_df, pd.DataFrame) else ""
        pcr_csv_string = pcr_df.to_csv(index=False) if isinstance(pcr_df, pd.DataFrame) else ""

        # Encode these CSV strings to Base64
        primer_b64 = base64.b64encode(primer_csv_string.encode())
        primer_data_url = f'data:text/csv;base64,{primer_b64.decode()}'
        pcr_b64 = base64.b64encode(pcr_csv_string.encode())
        pcr_data_url = f'data:text/csv;base64,{pcr_b64.decode()}'

        # Encode Genbank file
        genbank_b64 = base64.b64encode(genbank_string.encode())
        genbank_data_url = f'data:text/plain;base64,{genbank_b64.decode()}'

        return (
            genbank_data_url, 
            primer_df.to_dict('records') if isinstance(primer_df, pd.DataFrame) else [], 
            [{"name": i, "id": i} for i in primer_df.columns] if isinstance(primer_df, pd.DataFrame) else [], 
            pcr_df.to_dict('records') if isinstance(pcr_df, pd.DataFrame) else [], 
            [{"name": i, "id": i} for i in pcr_df.columns] if isinstance(pcr_df, pd.DataFrame) else [], 
            overhang_data,  # use the overhang_data here
            [{"name": i, "id": i} for i in overhang_df.columns] if isinstance(overhang_df, pd.DataFrame) else [],
            primer_data_url,  
            pcr_data_url  
        )

