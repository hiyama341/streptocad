# callback1.py
from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
import pandas as pd
from Bio.Restriction import StuI
import io
import os
import sys
import base64
from teemi.design.fetch_sequences import read_genbank_files, read_fasta_files
from pydna.dseqrecord import Dseqrecord
import zipfile
import csv
from urllib.parse import quote
from datetime import datetime, timezone
import tempfile
import logging
from Bio.Restriction import *
from Bio import Restriction

import sys
import io

from .workflow_logging import (
    capture_workflow_output,
    configure_logging,
    workflow_log,
)

# Log to stdout, and capture each workflow run's own output for the UI
configure_logging()

# Create a logger
logger = logging.getLogger(__name__)

# functions from StreptoCAD
# Local module imports
module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.primers.primer_analysis import analyze_primers_and_hairpins
from streptocad.cloning.plasmid_processing import assemble_and_process_plasmids
from streptocad.cloning.pcr_simulation import perform_pcr_on_sequences
from streptocad.sequence_loading.sequence_loading import (
    load_and_process_plasmid,
    load_and_process_genome_sequences,
)
from streptocad.primers.primer_generation import (
    generate_primer_dataframe,
    create_idt_order_dataframe,
)
from streptocad.utils import (
    extract_metadata_to_dataframe,
    generate_header,
)
from streptocad.output_packaging import OutputPackage
from streptocad.primers.idt_plates import idt_order_df_to_idt_plates


def register_workflow_1_callbacks(app):
    @app.callback(
        [
            Output("primers-output-table_1", "data"),
            Output("primers-output-table_1", "columns"),
            Output("pcr-table_1", "data"),
            Output("pcr-table_1", "columns"),
            Output("analyzed-primers-table_1", "data"),
            Output("analyzed-primers-table_1", "columns"),
            Output("download-data-and-protocols-link_1", "href"),
            Output("error-dialog_1", "message"),
            Output("error-dialog_1", "displayed"),
            Output("plasmid-metadata-table_1", "data"),
            Output("plasmid-metadata-table_1", "columns"),
        ],
        [Input("submit-settings-button_1", "n_clicks")],
        [
            State({"type": "upload-component", "index": "sequences"}, "contents"),
            State({"type": "upload-component", "index": "plasmid"}, "contents"),
            State({"type": "upload-component", "index": "sequences"}, "filename"),
            State({"type": "upload-component", "index": "plasmid"}, "filename"),
            State("forward-overhang-input_1", "value"),
            State("reverse-overhang-input_1", "value"),
            State("chosen-polymerase_1", "value"),
            State("melting-temperature_1", "value"),
            State("primer-concentration_1", "value"),
            State("restriction-enzymes_1", "value"),
            State("primer-lenght_1", "value"),
        ],
    )
    @capture_workflow_output
    def run_workflow(
        n_clicks,
        sequences_content,
        plasmid_content,
        sequences_filename,
        plasmid_filename,
        up_homology,
        dw_homology,
        chosen_polymerase,
        melting_temperature,
        primer_concentration,
        restriction_enzymes,
        primer_anneal_len,
    ):
        if n_clicks is None:
            raise PreventUpdate

        try:
            logging.info("Workflow 1 started")
            # Create a temporary directory
            with tempfile.TemporaryDirectory() as tempdir:
                # Decode the uploaded files
                def decode_contents(contents, filename):
                    content_type, content_string = contents.split(",")
                    decoded = base64.b64decode(content_string)

                    if "gb" in filename or "gbk" in filename:
                        return decoded.decode("utf-8")
                    elif "fa" in filename or "fasta" in filename:
                        return decoded.decode("utf-8")
                    else:
                        return None

                # Decode the sequences
                sequences_decoded = decode_contents(
                    sequences_content, sequences_filename
                )
                plasmid_decoded = decode_contents(plasmid_content, plasmid_filename)

                if not sequences_decoded or not plasmid_decoded:
                    raise ValueError(
                        "Unsupported file format. Please provide a FASTA or GenBank file."
                    )

                # Load the uploaded sequences and plasmid
                # Get file extension in lowercase
                _, ext = os.path.splitext(sequences_filename)
                ext = ext.lower()
                if ext in (".gb", ".gbk"):
                    input_sequences = read_genbank_files(io.StringIO(sequences_decoded))
                elif ext in (".fa", ".fasta"):
                    input_sequences = read_fasta_files(io.StringIO(sequences_decoded))
                else:
                    input_sequences = None

                sequences = [Dseqrecord(seq) for seq in input_sequences]
                input_plasmid = read_genbank_files(io.StringIO(plasmid_decoded))[0]
                plasmid = Dseqrecord(input_plasmid, circular=True)

                # Generate primers
                logging.info("Generating primers")
                primer_df = generate_primer_dataframe(
                    sequences,
                    melting_temperature,
                    chosen_polymerase,
                    primer_concentration,
                    up_homology,
                    dw_homology,
                    min_primer_len=primer_anneal_len,
                )
                logging.info(f"Primer DataFrame: {primer_df}")

                # Perform PCR
                logging.info("Performing PCR")
                list_of_amplicons = perform_pcr_on_sequences(primer_df, sequences)
                logging.info(f"List of amplicons: {list_of_amplicons}")

                # Create IDT order DataFrame
                logging.info("Creating IDT order DataFrame")
                idt_df = create_idt_order_dataframe(
                    primer_df, concentration="25nm", purification="STD"
                )
                logging.info(f"IDT DataFrame: {idt_df}")

                # Analyze primers and hairpins
                logging.info("Analyzing primers and hairpins")
                analyzed_primers_df = analyze_primers_and_hairpins(primer_df)
                logging.info(f"Analyzed primers DataFrame: {analyzed_primers_df}")

                restriction_enzymes = restriction_enzymes.split(",")
                enzymes_for_repair_template_integration = [
                    getattr(Restriction, str(enzyme)) for enzyme in restriction_enzymes
                ]

                # Assemble plasmids
                logging.info("Assembling plasmids")
                assembled_plasmids, assembly_results = assemble_and_process_plasmids(
                    plasmid,
                    list_of_amplicons,
                    enzymes=enzymes_for_repair_template_integration,
                    save_plasmids=False,
                    save_path="../../data/plasmids/pOEX_overexpression_plasmids",
                )
                logging.info(f"Assembly results: {assembly_results}")

                integration_names = [names.name for names in list_of_amplicons]
                plasmid_metadata_df = extract_metadata_to_dataframe(
                    assembled_plasmids, plasmid, integration_names
                )

                # Prepare outputs for the DataTable
                primers_columns = [{"name": col, "id": col} for col in idt_df.columns]
                primers_data = idt_df.to_dict("records")

                pcr_columns = [{"name": col, "id": col} for col in primer_df.columns]
                pcr_data = primer_df.to_dict("records")

                analyzed_primers_columns = [
                    {"name": col, "id": col} for col in analyzed_primers_df.columns
                ]
                analyzed_primers_data = analyzed_primers_df.to_dict("records")

                # Prepare columns and data for the plasmid metadata DataTable
                plasmid_metadata_columns = [
                    {"name": col, "id": col} for col in plasmid_metadata_df.columns
                ]
                plasmid_metadata_data = plasmid_metadata_df.to_dict("records")

                ### DATA PACKAGE:
                # for the assembly overview.
                header_text = generate_header(
                    assembled_plasmids,
                    sequences,
                    idt_df,
                    workflow_log(),
                )

                input_values = {
                    "polymerase_settings": {
                        "chosen_polymerase": chosen_polymerase,
                        "melting_temperature": melting_temperature,
                        "primer_concentration": primer_concentration,
                        "restriction_enzymes": restriction_enzymes,
                    },
                    "overlapping_sequences": {
                        "up_homology": str(up_homology),
                        "dw_homology": str(dw_homology),
                    },
                }

                # Data and time
                timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%SZ")

                package = OutputPackage(
                    workflow_id="w1",
                    outputs=[
                        {"role": "plasmid", "content": assembled_plasmids},
                        {"role": "plasmid.index", "content": plasmid_metadata_df},
                        {"role": "primer.pcr", "content": primer_df},
                        {"role": "primer.order_idt", "content": idt_df},
                        {
                            "role": "primer.order_idt_plate",
                            "content": idt_order_df_to_idt_plates(idt_df),
                        },
                        {"role": "analysis.primer_qc", "content": analyzed_primers_df},
                        {"role": "analysis.assembly_overview", "content": header_text},
                        {"role": "analysis.run_log", "content": workflow_log()},
                    ],
                    inputs=[
                        {"role": "input.sequences", "content": sequences},
                        {"role": "input.plasmid", "content": plasmid},
                    ],
                    parameters=input_values,
                    protocols=["conjugation", "overexpression", "troubleshooting"],
                )

                zip_content = package.to_zip_bytes(timestamp)

                # Encode the zip content for download
                data_package_encoded = base64.b64encode(zip_content).decode("utf-8")
                data_package_download_link = (
                    f"data:application/zip;base64,{data_package_encoded}"
                )

                logging.info("Workflow 1 completed successfully")

            return (
                primers_data,
                primers_columns,
                pcr_data,
                pcr_columns,
                analyzed_primers_data,
                analyzed_primers_columns,
                data_package_download_link,
                "",
                False,
                plasmid_metadata_data,
                plasmid_metadata_columns,
            )

        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
            error_message = f"An error occurred: {str(e)}\n\nLog:\n{workflow_log()}"
            return [], [], [], [], [], [], "", error_message, True, [], []
