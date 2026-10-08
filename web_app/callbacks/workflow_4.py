import sys
import os
import zipfile
import base64
import csv
import logging
from datetime import datetime, timezone

import pandas as pd
from Bio import SeqIO
from Bio.Restriction import NcoI
from pydna.dseqrecord import Dseqrecord
from teemi.design.fetch_sequences import read_genbank_files
from Bio.Restriction import *
from Bio import Restriction

from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import quote

import tempfile

module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.sequence_loading.sequence_loading import (
    load_and_process_gene_sequences,
    load_and_process_plasmid,
    load_and_process_genome_sequences,
    check_and_convert_input,
    annotate_dseqrecord,
    process_specified_gene_sequences_from_record,
)
from streptocad.utils import extract_metadata_to_dataframe
from streptocad.output_packaging import OutputPackage, RunLogCapture
from streptocad.crispr.guideRNA_crispri import extract_sgRNAs_for_crispri, SgRNAargs
from streptocad.cloning.ssDNA_bridging import (
    assemble_plasmids_by_ssDNA_bridging,
    make_ssDNA_oligos,
)
from streptocad.primers.primer_generation import (
    create_idt_order_dataframe,
    primers_to_IDT,
)
from streptocad.primers.idt_plates import idt_order_df_to_idt_plates


# Setup logging
logging.basicConfig(
    level=logging.INFO,  # Set to INFO to capture INFO, WARNING, ERROR, and CRITICAL messages
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),  # Log to the console (stdout)
    ],
)

# Create a logger
logger = logging.getLogger(__name__)


def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(name, "wb") as fp:
        fp.write(base64.decodebytes(data))


def register_workflow_4_callbacks(app):
    @app.callback(
        [
            Output("primer-table_4", "data"),
            Output("primer-table_4", "columns"),
            Output("download-data-and-protocols-link_4", "href"),
            Output("mutated-sgrna-table_4", "data"),
            Output("mutated-sgrna-table_4", "columns"),
            Output("plasmid-metadata-table_4", "data"),
            Output("plasmid-metadata-table_4", "columns"),
            Output("error-dialog_4", "message"),
            Output("error-dialog_4", "displayed"),
        ],
        [Input("submit-settings-button_4", "n_clicks")],
        [
            State(
                {"type": "upload-component", "index": "genome-file-4"}, "contents"
            ),  # This matches the layout
            State(
                {"type": "upload-component", "index": "single-vector-4"}, "contents"
            ),  # This matches the layout
            State(
                {"type": "upload-component", "index": "genome-file-4"}, "filename"
            ),  # This matches the layout
            State(
                {"type": "upload-component", "index": "single-vector-4"}, "filename"
            ),  # This matches the layout
            State("genes-to-KO_4", "value"),
            State("forward-overhang-input_4", "value"),
            State("reverse-overhang-input_4", "value"),
            State("gc-upper_4", "value"),
            State("gc-lower_4", "value"),
            State("off-target-seed_4", "value"),
            State("off-target-upper_4", "value"),
            State("cas-type_4", "value"),
            State("number-of-sgRNAs-per-group_4", "value"),
            State("extension-to-promoter-region_4", "value"),
            State("restriction-enzymes_4", "value"),
        ],
    )
    def run_workflow(
        n_clicks,
        genome_content,
        vector_content,
        genome_filename,
        vector_filename,
        genes_to_KO,
        up_homology,
        dw_homology,
        gc_upper,
        gc_lower,
        off_target_seed,
        off_target_upper,
        cas_type,
        number_of_sgRNAs_per_group,
        extension_to_promoter_region,
        restriction_enzymes,
    ):
        if n_clicks is None:
            raise PreventUpdate

        log_capture = RunLogCapture().start()
        log_stream = log_capture.stream

        try:
            logging.info("Workflow 4 started")

            with tempfile.TemporaryDirectory() as tempdir:
                genome_path = os.path.join(tempdir, genome_filename)
                vector_path = os.path.join(tempdir, vector_filename)

                logging.info(f"Saving uploaded files to temporary directory: {tempdir}")
                save_file(genome_path, genome_content)
                save_file(vector_path, vector_content)

                logging.info("Reading genome and vector files")
                genome = load_and_process_genome_sequences(genome_path)[0]
                clean_plasmid = load_and_process_plasmid(vector_path)

                logging.info("Processing genes to KO")
                target_dict, genes_to_KO_list, annotation_input = (
                    check_and_convert_input(genes_to_KO)
                )
                if annotation_input:
                    genome = annotate_dseqrecord(genome, target_dict)
                logging.info(f"Genes to knock out: {genes_to_KO_list}")

                # Extract sgRNAs
                args = SgRNAargs(
                    genome,
                    genes_to_KO_list,
                    step=["find", "filter"],
                    gc_upper=gc_upper,
                    gc_lower=gc_lower,
                    off_target_seed=off_target_seed,
                    off_target_upper=off_target_upper,
                    cas_type=cas_type,
                    extension_to_promoter_region=extension_to_promoter_region,
                    target_non_template_strand=True,
                )

                logging.info("Extracting sgRNAs for CRISPRi")
                sgrna_df = extract_sgRNAs_for_crispri(args)
                logging.debug(f"sgRNA DataFrame: {sgrna_df}")

                filtered_df = sgrna_df.groupby("locus_tag").head(
                    number_of_sgRNAs_per_group
                )
                logging.debug(f"Filtered DataFrame: {filtered_df}")

                logging.info("Making ssDNA oligos")
                list_of_ssDNAs = make_ssDNA_oligos(
                    filtered_df,
                    upstream_ovh=Dseqrecord(up_homology),
                    downstream_ovh=Dseqrecord(dw_homology),
                )
                logging.debug(f"List of ssDNAs: {list_of_ssDNAs}")

                logging.info("Cutting plasmid")
                restriction_enzymes = restriction_enzymes.split(",")
                enzymes_for_repair_template_integration = [
                    getattr(Restriction, str(enzyme)) for enzyme in restriction_enzymes
                ]

                linearized_plasmid = sorted(
                    clean_plasmid.cut(enzymes_for_repair_template_integration),
                    key=lambda x: len(x),
                    reverse=True,
                )[0]
                logging.debug(f"Linearized plasmid: {linearized_plasmid}")

                logging.info("Assembling plasmid")
                sgRNA_vectors = assemble_plasmids_by_ssDNA_bridging(
                    list_of_ssDNAs, linearized_plasmid
                )
                logging.debug(f"sgRNA vectors: {sgRNA_vectors}")

                targeting_info = []
                for index, row in filtered_df.iterrows():
                    formatted_str = f"CRISPRi_{row['locus_tag']}_p{row['sgrna_loc']}"
                    targeting_info.append(formatted_str)

                for i in range(len(sgRNA_vectors)):
                    sgRNA_vectors[i].name = f"p{targeting_info[i]}_#{i + 1}"
                    sgRNA_vectors[i].id = sgRNA_vectors[i].name
                    sgRNA_vectors[
                        i
                    ].description = f"Assembled plasmid targeting {', '.join(genes_to_KO_list)} for single gene KNOCK-DOWN, assembled using StreptoCAD."

                logging.info("Generating primers for IDT")
                idt_primers = primers_to_IDT(list_of_ssDNAs)
                logging.debug(f"IDT primers DataFrame: {idt_primers}")

                # Prepare DataTables outputs
                primers_columns = [
                    {"name": col, "id": col} for col in idt_primers.columns
                ]
                primers_data = idt_primers.to_dict("records")

                # Metadata df
                integration_names = filtered_df.apply(
                    lambda row: f"sgRNA_{row['locus_tag']}({row['sgrna_loc']})", axis=1
                ).tolist()
                plasmid_metadata_df = extract_metadata_to_dataframe(
                    sgRNA_vectors, clean_plasmid, integration_names
                )
                # The ssDNA order sheet has Name/Sequence columns, not the
                # forward/reverse layout of a primer dataframe.
                idt_primer_plates = idt_order_df_to_idt_plates(idt_primers)

                input_values = {
                    "genes_to_knockout": genes_to_KO_list,
                    "filtering_metrics": {
                        "gc_upper": gc_upper,
                        "gc_lower": gc_lower,
                        "off_target_seed": off_target_seed,
                        "off_target_upper": off_target_upper,
                        "cas_type": cas_type,
                        "number_of_sgRNAs_per_group": number_of_sgRNAs_per_group,
                        "extension_to_promoter_region": extension_to_promoter_region,
                    },
                    "overlapping_sequences": {
                        "up_homology": str(up_homology),
                        "dw_homology": str(dw_homology),
                    },
                }

                # filtered sgrnas
                filtered_df_columns = [
                    {"name": col, "id": col} for col in filtered_df.columns
                ]
                filtered_df_data = filtered_df.to_dict("records")

                # metadata table
                plasmid_metadata_df_columns = [
                    {"name": col, "id": col} for col in plasmid_metadata_df.columns
                ]
                plasmid_metadata_df_data = plasmid_metadata_df.to_dict("records")

                timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%SZ")

                logging.info("Packaging outputs")
                package = OutputPackage(
                    workflow_id="w4",
                    outputs=[
                        {"role": "plasmid", "content": sgRNA_vectors},
                        {"role": "plasmid.index", "content": plasmid_metadata_df},
                        {"role": "primer.order_idt", "content": idt_primers},
                        {
                            "role": "primer.order_idt_plate",
                            "content": idt_primer_plates,
                        },
                        {"role": "sgrna.all", "content": sgrna_df},
                        {"role": "sgrna.selected", "content": filtered_df},
                        {"role": "analysis.run_log", "content": log_stream.getvalue()},
                    ],
                    inputs=[
                        {"role": "input.genome", "content": genome},
                        {"role": "input.plasmid", "content": clean_plasmid},
                    ],
                    parameters=input_values,
                    protocols=[
                        "conjugation",
                        "crispr_single_target",
                        "troubleshooting",
                    ],
                )

                zip_content = package.to_zip_bytes(timestamp)
                data_package_encoded = base64.b64encode(zip_content).decode("utf-8")
                data_package_download_link = (
                    f"data:application/zip;base64,{data_package_encoded}"
                )

                logging.info("Workflow completed successfully")

                return (
                    primers_data,
                    primers_columns,
                    data_package_download_link,
                    filtered_df_data,
                    filtered_df_columns,
                    plasmid_metadata_df_data,
                    plasmid_metadata_df_columns,
                    "",  # Empty message if no error occurred
                    False,  # Error dialog should not be displayed
                )

        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
            print(f"An error occurred: {str(e)}")  # Fallback print

            error_message = (
                f"An error occurred: {str(e)}\n\nLog:\n{log_stream.getvalue()}"
            )
            display_error = True
            # Must match the 9 Outputs declared on this callback, in order:
            # primer table data/columns, download href, filtered sgRNA data/columns,
            # plasmid metadata data/columns, error message, error displayed.
            return (
                [],
                [],
                "",
                [],
                [],
                [],
                [],
                error_message,
                display_error,
            )
        finally:
            log_capture.stop()
