#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Standard library imports
import os
import sys
import io
import zipfile
import base64
import csv
from datetime import datetime

# Third-party imports
import pandas as pd
from Bio import SeqIO
from Bio.Restriction import NcoI, NheI
from pydna.dseqrecord import Dseqrecord
from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import quote
import logging
import tempfile
from teemi.build.PCR import primer_tm_neb
import logging
import sys
import io
from Bio.Restriction import *
from Bio import Restriction


# Create a StringIO object to capture logs in memory
log_stream = io.StringIO()

# Remove any existing handlers
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

# Setup logging
logging.basicConfig(
    level=logging.INFO,  # Set to INFO to capture INFO, WARNING, ERROR, and CRITICAL messages
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),  # Log to the console (stdout)
        logging.StreamHandler(log_stream),  # Capture logs in StringIO
    ],
)

# Create a logger
logger = logging.getLogger(__name__)


# Local module imports
module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.sequence_loading.sequence_loading import (
    load_and_process_genome_sequences,
    load_and_process_plasmid,
    check_and_convert_input,
    annotate_dseqrecord,
    process_specified_gene_sequences_from_record,
)

from streptocad.utils import (
    polymerase_dict,
    dataframe_to_seqrecords,
    ProjectDirectory,
    extract_metadata_to_dataframe,
)
from streptocad.crispr.guideRNAcas3_9_12 import extract_sgRNAs, SgRNAargs
from streptocad.crispr.crispr_best import (
    identify_base_editing_sites,
    filter_sgrnas_for_base_editing,
    process_base_editing,
)

from streptocad.cloning.golden_gate_cloning import (
    GoldenGateCloning,
    create_overhang_dataframe,
    digest_amplicons_w_BsaI,
)

from streptocad.primers.primer_analysis import analyze_primers_and_hairpins
from streptocad.primers.primer_generation import (
    create_idt_order_dataframe,
    find_best_check_primers_from_genome,
)
from streptocad.cloning.plasmid_processing import annotate_plasmid_with_sgrnas


def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(name, "wb") as fp:
        fp.write(base64.decodebytes(data))


def register_workflow_3_callbacks(app):
    @app.callback(
        [
            Output("primer-table_3", "data"),
            Output("primer-table_3", "columns"),
            Output("pcr-table_3", "data"),
            Output("pcr-table_3", "columns"),
            Output("overhang-table_3", "data"),
            Output("overhang-table_3", "columns"),
            Output("download-data-and-protocols-link_3", "href"),
            Output("mutated-sgrna-table_3", "data"),
            Output("mutated-sgrna-table_3", "columns"),
            Output("plasmid-metadata-table_3", "data"),
            Output("plasmid-metadata-table_3", "columns"),
            Output("error-dialog_3", "message"),
            Output("error-dialog_3", "displayed"),
        ],
        [Input("submit-button_3", "n_clicks")],
        [
            State(
                {"type": "upload-component", "index": "genome-file-3"}, "contents"
            ),  # This matches the layout
            State(
                {"type": "upload-component", "index": "single-vector-3"}, "contents"
            ),  # This matches the layout
            State(
                {"type": "upload-component", "index": "genome-file-3"}, "filename"
            ),  # This matches the layout
            State(
                {"type": "upload-component", "index": "single-vector-3"}, "filename"
            ),  # This matches the layout
            State("genes-to-KO_3", "value"),
            State("sgRNA-handle-input_3", "value"),
            State("input-tm_3", "value"),
            State("gc-upper_3", "value"),
            State("gc-lower_3", "value"),
            State("off-target-seed_3", "value"),
            State("off-target-upper_3", "value"),
            State("cas-type_3", "value"),
            State("number-of-sgRNAs-per-group_3", "value"),
            State("only-stop-codons-checkbox_3", "value"),
            State("chosen-polymerase_3", "value"),
            State("melting-temperature_3", "value"),
            State("primer-concentration_3", "value"),
            State("flanking-region-number_3", "value"),
            State("restriction-overhang-f", "value"),
            State("restriction-overhang-r", "value"),
            State("backbone-overhang-f", "value"),
            State("backbone-overhang-r", "value"),
            State("cys4-sequence", "value"),
            State("editing_context_3", "value"),
            State("restriction-enzymes_3", "value"),
            State("checking-primer-length_3", "value"),
        ],
    )
    def run_workflow(
        n_clicks,
        genome_content,
        vector_content,
        genome_filename,
        vector_filename,
        genes_to_KO,
        sgRNA_handle_input,
        input_tm,
        gc_upper,
        gc_lower,
        off_target_seed,
        off_target_upper,
        cas_type,
        number_of_sgRNAs_per_group,
        only_stop_codons,
        chosen_polymerase,
        melting_temperature,
        primer_concentration,
        flanking_region_number,
        restriction_overhang_f,
        restriction_overhang_r,
        backbone_overhang_f,
        backbone_overhang_r,
        cys4_sequence,
        editing_context,
        restriction_enzymes,
        checking_primer_length,
    ):
        if n_clicks is None:
            raise PreventUpdate

        try:
            logger.info("Workflow 3 started")

            # Create a temporary directory
            with tempfile.TemporaryDirectory() as tempdir:
                genome_path = os.path.join(tempdir, genome_filename)
                vector_path = os.path.join(tempdir, vector_filename)

                # Save uploaded files to the temporary directory
                save_file(genome_path, genome_content)
                save_file(vector_path, vector_content)

                # Read the GenBank files from the saved paths
                logging.info("Reading genome and vector files")
                genome = load_and_process_genome_sequences(genome_path)[0]
                clean_plasmid = load_and_process_plasmid(vector_path)

                # Process genes to KO
                genes_to_KO_list = [gene.strip() for gene in genes_to_KO.split(",")]
                logging.info(f"Genes to knock out: {genes_to_KO_list}")

                # Check and convert input genes for annotation
                target_dict, genes_to_KO_list, annotation_input = (
                    check_and_convert_input(genes_to_KO_list)
                )
                if annotation_input:
                    genome = annotate_dseqrecord(genome, target_dict)
                logging.info("Annotation completed.")

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
                )
                sgrna_df = extract_sgRNAs(args)
                logging.info("sgRNA extraction completed.")

                # Further processing of gene sequences and base editing
                gene_sequences = process_specified_gene_sequences_from_record(
                    genome, genes_to_KO_list
                )
                logger.info(f"Gene sequences processed: {list(gene_sequences.keys())}")

                genes_to_KO_dict = {
                    locus_tag: gene_sequences[locus_tag]
                    for locus_tag in genes_to_KO_list
                    if locus_tag in gene_sequences
                }
                sgrna_df_with_editing = identify_base_editing_sites(sgrna_df)
                logger.info(
                    f"Base editing sites identified: {sgrna_df_with_editing.shape[0]} rows"
                )

                filtered_sgrna_df_for_base_editing = filter_sgrnas_for_base_editing(
                    sgrna_df_with_editing
                )
                logger.info(
                    f"sgRNAs filtered for base editing: {filtered_sgrna_df_for_base_editing.shape[0]} rows"
                )

                mutated_sgrna_df = process_base_editing(
                    filtered_sgrna_df_for_base_editing,
                    genes_to_KO_dict,
                    only_stop_codons=bool(only_stop_codons),
                    editing_context=bool(editing_context),
                )

                logger.info(f"Base editing applied: {mutated_sgrna_df.shape[0]} rows")
                filtered_df = mutated_sgrna_df.groupby("locus_tag").head(
                    number_of_sgRNAs_per_group
                )
                logger.info(f"sgRNAs filtered by group: {filtered_df.shape[0]} rows")

                # Prepare sgRNA list and perform Golden Gate Cloning
                sgRNA_list = dataframe_to_seqrecords(filtered_df)
                logger.info(f"sgRNA list prepared: {len(sgRNA_list)} sequences")

                sgRNA_handle_cys4_sites = [
                    Dseqrecord(sgRNA_handle_input, name="sgRNA_handle_cys4")
                ] * len(sgRNA_list)
                golden_gate = GoldenGateCloning(
                    sgRNA_list,
                    sgRNA_handle_cys4_sites,
                    target_tm=input_tm,
                    restriction_overhang_f=restriction_overhang_f,
                    restriction_overhang_r=restriction_overhang_r,
                    backbone_overhang_f=backbone_overhang_f,
                    backbone_overhang_r=backbone_overhang_r,
                    cys4=cys4_sequence,
                    tm_function=primer_tm_neb,
                    polymerase=chosen_polymerase,
                )
                logger.info("Golden Gate Cloning setup completed.")

                # Primer generation, analysis
                primer_df = golden_gate.generate_primer_dataframe()
                logger.info(f"Primer dataframe generated: {primer_df.shape[0]} rows")

                # TODO maybe remove this redundancy.
                analyze_primers_and_hairpins(primer_df)
                logger.info("Primers analyzed for hairpins.")

                idt_df = create_idt_order_dataframe(
                    primer_df, concentration="25nm", purification="STD"
                )

                # Generate checkig primers
                checking_primers_df = find_best_check_primers_from_genome(
                    genome,
                    genes_to_KO_list,
                    flanking_region=flanking_region_number,
                    target_tm=melting_temperature,
                    primer_concentration=primer_concentration,
                    polymerase=chosen_polymerase,
                    limit=checking_primer_length,
                )
                idt_df2 = create_idt_order_dataframe(checking_primers_df)
                full_idt = pd.concat([idt_df, idt_df2])
                logger.info(f"IDT order dataframe created: {idt_df.shape[0]} rows")

                list_of_amplicons = golden_gate.simulate_pcrs()

                # Digest and assemble plasmid
                overhangs = create_overhang_dataframe(list_of_amplicons)
                logger.info(f"Overhang dataframe created: {overhangs.shape[0]} rows")

                digest_amplicons = digest_amplicons_w_BsaI(list_of_amplicons)
                logger.info(f"Amplicons digested: {len(digest_amplicons)} fragments")
                # Cut plasmid
                restriction_enzymes = restriction_enzymes.split(",")
                enzymes_for_repair_template_integration = [
                    getattr(Restriction, str(enzyme)) for enzyme in restriction_enzymes
                ]

                linear_plasmid = sorted(
                    clean_plasmid.cut(enzymes_for_repair_template_integration),
                    key=lambda x: len(x),
                    reverse=True,
                )[0]
                logger.info("Plasmid linearized.")

                for amplicon in digest_amplicons:
                    linear_plasmid.seq += amplicon.seq
                rec_vec = linear_plasmid.looped()
                logger.info("Plasmid recircularized.")

                # Annotate and generate plasmid metadata
                annotate_plasmid_with_sgrnas(rec_vec, filtered_df)
                logger.info("Plasmid annotated with sgRNAs.")

                integration_names = filtered_df.apply(
                    lambda row: f"sgRNA_{row['locus_tag']}({row['sgrna_loc']})", axis=1
                ).tolist()
                integration_names = [";".join(integration_names)]
                rec_vec.name = f"{clean_plasmid.name}_{integration_names}"

                plasmid_metadata_df = extract_metadata_to_dataframe(
                    [rec_vec], clean_plasmid, integration_names
                )

                logger.info(
                    f"Plasmid metadata extracted: {plasmid_metadata_df.shape[0]} rows"
                )

                # Prepare DataTables outputs
                primers_columns = [{"name": col, "id": col} for col in full_idt.columns]
                primers_data = full_idt.to_dict("records")

                pcr_columns = [{"name": col, "id": col} for col in primer_df.columns]
                pcr_data = primer_df.to_dict("records")

                overhang_columns = [
                    {"name": col, "id": col} for col in overhangs.columns
                ]
                overhang_data = overhangs.to_dict("records")

                filtered_sgrna_columns = [
                    {"name": col, "id": col} for col in filtered_df.columns
                ]
                filtered_sgrna_data = filtered_df.to_dict("records")

                plasmid_metadata_columns = [
                    {"name": col, "id": col} for col in plasmid_metadata_df.columns
                ]
                plasmid_metadata_data = plasmid_metadata_df.to_dict("records")

                # Generate download link for data package
                input_files = [
                    {"name": "input_genome.gb", "content": genome},
                    {"name": "input_plasmid.gb", "content": clean_plasmid},
                ]
                output_files = [
                    {"name": "mcBEST_w_sgRNAs.gb", "content": rec_vec},
                    {"name": "01_pcr_df.csv", "content": primer_df},
                    {"name": "02_full_idt.csv", "content": full_idt},
                    {"name": "03_mutated_sgrna_df.csv", "content": mutated_sgrna_df},
                    {"name": "04_filtered_sgrna_df.csv", "content": filtered_df},
                    {
                        "name": "05_plasmid_metadata_df.csv",
                        "content": plasmid_metadata_df,
                    },
                    {"name": "06_overhang_df.csv", "content": overhangs},
                ]
                input_values = {
                    "genes_to_knockout": genes_to_KO_list,
                    "filtering_metrics": {
                        "gc_upper": gc_upper,
                        "gc_lower": gc_lower,
                        "off_target_seed": off_target_seed,
                        "off_target_upper": off_target_upper,
                        "cas_type": cas_type,
                        "number_of_sgRNAs_per_group": number_of_sgRNAs_per_group,
                    },
                    "polymerase_settings": {
                        "chosen_polymerase": chosen_polymerase,
                        "melting_temperature": melting_temperature,
                        "primer_concentration": primer_concentration,
                        "flanking_region": flanking_region_number,
                    },
                    "overlapping_sequences": {
                        "restriction_overhang_f": restriction_overhang_f,
                        "restriction_overhang_r": restriction_overhang_r,
                        "backbone_overhang_f": backbone_overhang_f,
                        "backbone_overhang_r": backbone_overhang_r,
                        "cys4": cys4_sequence,
                        "sgRNA_handle_cys4_site": str(sgRNA_handle_cys4_sites[0].seq),
                    },
                }
                markdown_file_paths = [
                    "protocols/conjugation_protcol.md",
                    "protocols/multi_target_crispr_plasmid_protcol.md",
                ]

                project_directory = ProjectDirectory(
                    project_name=f"CRISPR_mcBEST_workflow_{datetime.utcnow().isoformat()}",
                    input_files=input_files,
                    output_files=output_files,
                    input_values=input_values,
                    markdown_file_paths=markdown_file_paths,
                )

                zip_content = project_directory.create_directory_structure(
                    create_directories=True
                )
                data_package_encoded = base64.b64encode(zip_content).decode("utf-8")
                data_package_download_link = (
                    f"data:application/zip;base64,{data_package_encoded}"
                )

                logger.info("Workflow 3 completed successfully")

            return (
                primers_data,
                primers_columns,
                pcr_data,
                pcr_columns,
                overhang_data,
                overhang_columns,
                data_package_download_link,
                filtered_sgrna_data,
                filtered_sgrna_columns,
                plasmid_metadata_data,
                plasmid_metadata_columns,
                "",  # Empty message if no error occurred
                False,  # Error dialog should not be displayed
            )

        except Exception as e:
            logger.error(f"An error occurred: {str(e)}")
            print(f"An error occurred: {str(e)}")  # Fallback print

            error_message = (
                f"An error occurred: {str(e)}\n\nLog:\n{log_stream.getvalue()}"
            )
            display_error = True
            return (
                [],
                [],
                [],
                [],
                "",
                [],
                [],
                [],
                [],
                [],
                [],
                error_message,
                display_error,
            )
