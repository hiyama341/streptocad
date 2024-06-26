#!/usr/bin/env python
# MIT License

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


from streptocad.cloning.ssDNA_bridging import assemble_plasmids_by_ssDNA_bridging, make_ssDNA_oligos
from streptocad.utils import polymerase_dict
from streptocad.primers.primer_generation import primers_to_IDT



# Main function for data callbacks
def register_workflow_2_callbacks(app):

    # callback for toggle_temp_input and toggle_distance_input
    @app.callback(
        [Output('input-tm-ssDNA', 'style'),
        Output('input-distance-ssDNA', 'style'),
        Output('label-tm-ssDNA', 'style'),
        Output('label-distance-ssDNA', 'style')],
        Input('check-primers-checkbox', 'value')
    )
    def toggle_input_fields(checkbox_value):
        if checkbox_value:
            return {"display": "block"}, {"display": "block"}, {"display": "block"}, {"display": "block"}
        else:
            return {"display": "none"}, {"display": "none"}, {"display": "none"}, {"display": "none"}


    # callback for filename for the sgRNA file
    @app.callback(
        Output('uploaded-workflow2-sgrna-filename', 'children'),
        [Input('upload-single-sgrna', 'filename')]
    )
    def display_uploaded_sgrna_filename(filename):
        if not filename:
            return "No file selected."
        # Check file extension
        if not filename.endswith('.csv'):
            return "Invalid file type. Please upload a .csv file."
        return f"Uploaded file: {filename}"

    # callback for filename for the plasmid
    @app.callback(
        Output('uploaded-workflow2-vector-filename', 'children'),
        [Input('upload-single-vector', 'filename')]
    )

    def display_uploaded_vector_filename(filename):
        if not filename:
            return "No file selected."
        # Check file extension
        if not (filename.endswith('.gb') or filename.endswith('.gbk')):
            return "Invalid file type. Please upload a GenBank (.gb or .gbk) file. For example: pCRISPR–Cas9_plasmid_addgene.gbk"
        return f"Uploaded file: {filename}"


    ############################################################ 
    ## CALLBacks for SINGLE sgRNA integration tab
    ############################################################ 

    # Callback for single-sgRNA-integration
    @app.callback(
        [
            Output('primers-output-table', 'data'),
            Output('csv_download_link', 'href'),
            Output('genbank-file-single', 'href')
        ],
        Input('submit-overhangs-button', 'n_clicks'),
        State('upload-single-sgrna', 'contents'),
        State('upload-single-vector', 'contents'),
        State('forward-overhang-input', 'value'),
        State('reverse-overhang-input', 'value')
    )
    def generate_ssDNA_integration(n_clicks, sgrna_csv_content, genbank_content, overhang_start, overhang_end):
        # Check if button is clicked and contents are uploaded
        if n_clicks and sgrna_csv_content and genbank_content:
            # Decode the uploaded CSV and GenBank files
            sgrna_csv_decoded = base64.b64decode(sgrna_csv_content.split(",")[1]).decode("utf-8")
            genbank_decoded = base64.b64decode(genbank_content.split(",")[1]).decode("utf-8")
            
            # Convert decoded strings to pandas dataframe and genbank object
            data_frame_1 = pd.read_csv(io.StringIO(sgrna_csv_decoded))
            input_plasmid = read_genbank_files(io.StringIO(genbank_decoded))[0]
            input_plasmid = Dseqrecord(input_plasmid, circular=True)

            # Make ssDNA oligos
            data_frame_1.rename(columns={"Sequence": "sgrna", 'ORF':'locus_tag'}, inplace=True)
            ssDNA_list = make_ssDNA_oligos(data_frame_1, upstream_ovh=overhang_start, downstream_ovh=overhang_end)

            # Perform calculations
            output_df = primers_to_IDT(ssDNA_list)

            # For Primer display:
            primer_display_data = output_df.to_dict('records')

            # Digest backbone plasmid and assemble
            linearized_input = input_plasmid.cut(NcoI)[0]
            assembled_vectors = assemble_plasmids_by_ssDNA_bridging(ssDNA_list, linearized_input)

            # Renaming the assembled plasmids
            for i in range(len(assembled_vectors)):
                assembled_vectors[i].name = f"{ssDNA_list[i].name}"
                assembled_vectors[i].id = f"{ssDNA_list[i].name}"

            # make the output into a csv file
            csv_string = output_df.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
            csv_data = quote(csv_string)
            csv_download_link = f"data:text/csv;charset=utf-8,{csv_data}"

            # Provide a link for downloading GenBank files:
            encoded_genbank_files = [base64.b64encode(str(vector).encode()).decode() for vector in assembled_vectors]
            genbank_download_link = f"data:application/genbank;base64,{encoded_genbank_files[0]}"

            # Create a zip archive in memory
            zip_buffer = io.BytesIO()
            zip_data = None  # Initialize zip_data here
            with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                for vector in assembled_vectors:
                    # Convert each vector to GenBank format
                    genbank_content = vector.format("genbank")
                    zip_file.writestr(f"{vector.name}.gb", genbank_content)

            # Prepare the zip archive for download
            zip_buffer.seek(0)
            zip_data = base64.b64encode(zip_buffer.read()).decode('utf-8')
            genbank_download_link = f"data:application/zip;base64,{zip_data}"

            return primer_display_data, csv_download_link, genbank_download_link
        else:
            raise PreventUpdate
