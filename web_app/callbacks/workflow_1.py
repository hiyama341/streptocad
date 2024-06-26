from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import pandas as pd
from Bio.Restriction import StuI
import io
import contextlib
import pkg_resources
import os 
import sys
import base64
from teemi.design.fetch_sequences import read_genbank_files
from pydna.dseqrecord import Dseqrecord
import zipfile
import csv
from urllib.parse import quote



# functions from StreptoCAD
# Local module imports
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.wet_lab.gel_simulation import simulate_gel_electrophoresis
from streptocad.primers.primer_analysis import analyze_primers_and_hairpins
from streptocad.cloning.plasmid_processing import assemble_and_process_plasmids
from streptocad.cloning.pcr_simulation import perform_pcr_on_sequences
from streptocad.sequence_loading.sequence_loading import load_and_process_plasmid, load_and_process_genome_sequences
from streptocad.primers.primer_generation import generate_primer_dataframe, create_idt_order_dataframe
from streptocad.utils import polymerase_dict, format_and_print_values



def register_workflow_1_callbacks(app):
    @app.callback(
        [
            Output('primers-output-table_1', 'data'),
            Output('primers-output-table_1', 'columns'),
            Output('pcr-table_1', 'data'),
            Output('pcr-table_1', 'columns'),
            Output('analyzed-primers-table_1', 'data'),
            Output('analyzed-primers-table_1', 'columns'),
            Output('genbank-file-single_1', 'href'),
            Output('primers_download_link_1', 'href'),
            Output('download-pcr-link_1', 'href'),
            Output('download-analyzed-primers-link_1', 'href'),
            Output('download-data-and-protocols-link_1', 'href')
        ],
        [
            Input('submit-settings-button_1', 'n_clicks')
        ],
        [
            State({'type': 'upload-input', 'index': 'sequences'}, 'contents'),
            State({'type': 'upload-input', 'index': 'plasmid'}, 'contents'),
            State({'type': 'upload-input', 'index': 'sequences'}, 'filename'),
            State({'type': 'upload-input', 'index': 'plasmid'}, 'filename'),
            State('forward-overhang-input_1', 'value'),
            State('reverse-overhang-input_1', 'value'),
            State('chosen-polymerase_1', 'value'),
            State('melting-temperature_1', 'value'),
            State('primer-concentration_1', 'value'),
            State('primer-number-increment_1', 'value')
        ]
    )
    def run_workflow(n_clicks, sequences_content, plasmid_content, sequences_filename, plasmid_filename, up_homology, dw_homology, 
                     chosen_polymerase, melting_temperature, primer_concentration, primer_number_increment):
        if n_clicks is None:
            raise PreventUpdate

        try:
            # Decode the uploaded files
            def decode_contents(contents, filename):
                content_type, content_string = contents.split(',')
                decoded = base64.b64decode(content_string)

                if 'gb' in filename or 'gbk' in filename:
                    return decoded.decode('utf-8')
                elif 'fa' in filename or 'fasta' in filename:
                    return decoded.decode('utf-8')
                else:
                    return None
            
            # Decode the sequences
            sequences_decoded = decode_contents(sequences_content, sequences_filename)
            plasmid_decoded = decode_contents(plasmid_content, plasmid_filename)

            if not sequences_decoded or not plasmid_decoded:
                raise ValueError("Unsupported file format. Please provide a FASTA or GenBank file.")

            # Load the uploaded sequences and plasmid
            input_sequences = read_genbank_files(io.StringIO(sequences_decoded))
            sequences = [Dseqrecord(seq) for seq in input_sequences]
            input_plasmid = read_genbank_files(io.StringIO(plasmid_decoded))[0]
            plasmid = Dseqrecord(input_plasmid, circular=True)

            # Generate primers
            print("Generating primers")
            primer_df = generate_primer_dataframe(sequences, 
                                                  melting_temperature, 
                                                  chosen_polymerase, 
                                                  primer_concentration,
                                                  up_homology, dw_homology, 
                                                  primer_number_increment)
            print("Primer DataFrame:", primer_df)

            # Perform PCR
            print("Performing PCR")
            list_of_amplicons = perform_pcr_on_sequences(primer_df, sequences)
            print("List of amplicons:", list_of_amplicons)

            # Create IDT order DataFrame
            print("Creating IDT order DataFrame")
            idt_df = create_idt_order_dataframe(primer_df, concentration="25nm", purification="STD")
            print("IDT DataFrame:", idt_df)

            # Analyze primers and hairpins
            print("Analyzing primers and hairpins")
            analyzed_primers_df = analyze_primers_and_hairpins(primer_df)
            print("Analyzed primers DataFrame:", analyzed_primers_df)

            # Simulate gel electrophoresis
            print("Simulating gel electrophoresis")
            gel_simulation = simulate_gel_electrophoresis(list_of_amplicons)
            print("Gel simulation completed")

            # Assemble plasmids
            print("Assembling plasmids")
            assembled_plasmids, assembly_results = assemble_and_process_plasmids(plasmid, list_of_amplicons, 
                                                                                 enzyme=StuI, 
                                                                                 save_plasmids=False, 
                                                                                 save_path="../../data/plasmids/pOEX_overexpression_plasmids")
            print("Assembly results:", assembly_results)

            # Prepare outputs for the DataTable
            primers_columns = [{"name": col, "id": col} for col in idt_df.columns]
            primers_data = idt_df.to_dict('records')

            pcr_columns = [{"name": col, "id": col} for col in primer_df.columns]
            pcr_data = primer_df.to_dict('records')

            analyzed_primers_columns = [{"name": col, "id": col} for col in analyzed_primers_df.columns]
            analyzed_primers_data = analyzed_primers_df.to_dict('records')

            # Create download links for CSV data
            primer_df_string = idt_df.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
            primer_data_encoded = quote(primer_df_string)
            primer_download_link = f"data:text/csv;charset=utf-8,{primer_data_encoded}"

            pcr_df_string = primer_df.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
            pcr_data_encoded = quote(pcr_df_string)
            pcr_download_link = f"data:text/csv;charset=utf-8,{pcr_data_encoded}"

            analyzed_primers_string = analyzed_primers_df.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
            analyzed_primers_data_encoded = quote(analyzed_primers_string)
            analyzed_primers_download_link = f"data:text/csv;charset=utf-8,{analyzed_primers_data_encoded}"

            # Create download link for GenBank files
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                for vector in assembled_plasmids:
                    genbank_content = vector.format("genbank")
                    zip_file.writestr(f"{vector.name}.gb", genbank_content)

            zip_buffer.seek(0)
            zip_data = base64.b64encode(zip_buffer.read()).decode('utf-8')
            genbank_download_link = f"data:application/zip;base64,{zip_data}"

            # Create data package download link
            input_dict = {
                '######Inputs####': '',
                "Name of sequences": [seq.name for seq in sequences],
                "Plasmid name": plasmid.name,
                "Up homology": up_homology,
                "DW homology": dw_homology,
                "Chosen polymerase": chosen_polymerase,
                "Melting temperature": melting_temperature,
                "Primer_concentration": primer_concentration,
                "Primer number increment": primer_number_increment,
                '####### OUTPUTS #########' : '',
                "Primer df": primer_df,
                'Primer analysis': analyzed_primers_df, 
                "Amplicons":[(amplicon.figure()) for amplicon in list_of_amplicons],
                'Lenght of amplicons' : [(amplicon.name, len(amplicon)) for amplicon in list_of_amplicons],
                'Plasmid assembly':[plasmid.figure() for plasmid in assembly_results],
                }
            
            description = (
                    "Your very own FAIR data file for your overexpression workflow.\n"
                    "This document captures all the values that were put into the application "
                    "to ensure that the experiment can be replicated accurately.\n"
                    "Below are the details of the inputs/output:\n"
                )

            data_and_protocols_package = format_and_print_values(input_dict, description, spacing = 10)
            data_package_encoded = base64.b64encode(data_and_protocols_package.encode('utf-8')).decode('utf-8')
            data_package_download_link = f"data:text/plain;base64,{data_package_encoded}"  
        
            return primers_data, primers_columns, pcr_data, pcr_columns, analyzed_primers_data, analyzed_primers_columns, genbank_download_link, primer_download_link, pcr_download_link, analyzed_primers_download_link, data_package_download_link

        except Exception as e:
            print("An error occurred:", str(e))
            raise PreventUpdate