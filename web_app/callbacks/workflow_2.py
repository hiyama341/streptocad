#!/usr/bin/env python
# MIT License

# Standard library imports
import sys
import os
import io
import zipfile
import base64
import csv
from datetime import datetime

# Third-party imports
import pandas as pd
from Bio import SeqIO
from Bio.Restriction import NcoI
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord
from teemi.design.fetch_sequences import read_genbank_files

import dash
from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import quote
import logging
import tempfile

# Create a StringIO object to capture logs in memory
log_stream = io.StringIO()

# Remove any existing handlers
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

# Setup logging
logging.basicConfig(
    level=logging.INFO,  # Set to INFO to capture INFO, WARNING, ERROR, and CRITICAL messages
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),  # Log to the console (stdout)
        logging.StreamHandler(log_stream)   # Capture logs in StringIO
    ]
)

# Create a logger
logger = logging.getLogger(__name__)

# Local module imports
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.sequence_loading.sequence_loading import load_and_process_plasmid, load_and_process_genome_sequences, annotate_dseqrecord, check_and_convert_input, process_specified_gene_sequences_from_record
from streptocad.utils import polymerase_dict, ProjectDirectory, extract_metadata_to_dataframe
from streptocad.crispr.guideRNAcas3_9_12 import extract_sgRNAs, SgRNAargs
from streptocad.cloning.ssDNA_bridging import assemble_plasmids_by_ssDNA_bridging, make_ssDNA_oligos
from streptocad.crispr.crispr_best import identify_base_editing_sites, filter_sgrnas_for_base_editing, process_base_editing
from streptocad.cloning.plasmid_processing import annotate_plasmid_with_sgrnas
from streptocad.primers.primer_generation import find_best_check_primers_from_genome, create_idt_order_dataframe, primers_to_IDT

def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(name, "wb") as fp:
        fp.write(base64.decodebytes(data))

def register_workflow_2_callbacks(app):

    @app.callback(
        [
            Output('primers-output-table_2', 'data'),
            Output('primers-output-table_2', 'columns'),
            Output('pcr-table_2', 'data'),
            Output('pcr-table_2', 'columns'),
            Output('download-data-and-protocols-link_2', 'href'),
            Output('filtered-df-table', 'data'),
            Output('filtered-df-table', 'columns'),
            Output('error-dialog_2', 'message'),
            Output('error-dialog_2', 'displayed'),
            Output('plasmid-metadata-table_2', 'data'),  # New Output for plasmid metadata DataTable
            Output('plasmid-metadata-table_2', 'columns') # New Output for plasmid metadata DataTable columns
        ],
        [
            Input('submit-settings-button_2', 'n_clicks')
        ],
        [
            State('upload-genome-file_2', 'contents'),
            State('upload-single-vector_2', 'contents'),
            State('upload-genome-file_2', 'filename'),
            State('upload-single-vector_2', 'filename'),
            State('genes-to-KO_2', 'value'),
            State('forward-overhang-input_2', 'value'),
            State('reverse-overhang-input_2', 'value'),
            State('gc-upper_2', 'value'),
            State('gc-lower_2', 'value'),
            State('off-target-seed_2', 'value'),
            State('off-target-upper_2', 'value'),
            State('cas-type_2', 'value'),
            State('number-of-sgRNAs-per-group_2', 'value'),
            State('only-stop-codons-checkbox_2', 'value'),
            State('chosen-polymerase_2', 'value'),
            State('melting-temperature_2', 'value'),
            State('primer-concentration_2', 'value'),
            State('primer-number-increment_2', 'value'),
            State('flanking-region-number_2', 'value')
        ]
    )
    def run_workflow(n_clicks, genome_content, vector_content, genome_filename, vector_filename, genes_to_KO, 
                     up_homology, dw_homology, gc_upper, gc_lower, off_target_seed, off_target_upper, cas_type, 
                     number_of_sgRNAs_per_group, only_stop_codons, chosen_polymerase, melting_temperature, 
                     primer_concentration, primer_number_increment, flanking_region_number):
        if n_clicks is None:
            raise PreventUpdate

        try:
            logging.info("Workflow 2 started")

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
                genes_to_KO_list = [gene.strip() for gene in genes_to_KO.split(',')]
                logging.info(f"Genes to knock out: {genes_to_KO_list}")

                # Check and convert input genes for annotation
                target_dict, genes_to_KO_list, annotation_input = check_and_convert_input(genes_to_KO_list)
                if annotation_input:
                    genome = annotate_dseqrecord(genome, target_dict)
                logging.info("Annotation completed.")

                # Extract sgRNAs
                args = SgRNAargs(genome, 
                                genes_to_KO_list,
                                step=['find', 'filter'],
                                gc_upper=gc_upper,
                                gc_lower=gc_lower,
                                off_target_seed=off_target_seed,
                                off_target_upper=off_target_upper,
                                cas_type=cas_type)
                sgrna_df = extract_sgRNAs(args)
                logging.info("sgRNA extraction completed.")

                # Load gene sequences
                gene_sequences = process_specified_gene_sequences_from_record(genome, genes_to_KO_list)
                genes_to_KO_dict = {locus_tag: gene_sequences[locus_tag] for locus_tag in genes_to_KO_list if locus_tag in gene_sequences}

                # Identify and annotate base editing sites
                sgrna_df_with_editing = identify_base_editing_sites(sgrna_df)

                # Filter out only sgRNAs that result in base-editing
                filtered_sgrna_df_for_base_editing = filter_sgrnas_for_base_editing(sgrna_df_with_editing)
                logging.info("Filtering sgRNAs for base editing completed.")

                # Process the DataFrame to apply C-to-T mutations
                mutated_sgrna_df = process_base_editing(filtered_sgrna_df_for_base_editing, genes_to_KO_dict, only_stop_codons=bool(only_stop_codons))

                # Filter the DataFrame to retain only up to 5 sgRNA sequences per locus_tag
                filtered_df = mutated_sgrna_df.groupby('locus_tag').head(number_of_sgRNAs_per_group)
                logging.info(f"Filtered DataFrame: {filtered_df}")

                # Make oligos
                list_of_ssDNAs = make_ssDNA_oligos(filtered_df, upstream_ovh=up_homology, downstream_ovh=dw_homology)

                # Cut plasmid
                linearized_plasmid = sorted(clean_plasmid.cut(NcoI), key=lambda x: len(x), reverse=True)[0]

                # Assemble plasmid
                sgRNA_vectors = assemble_plasmids_by_ssDNA_bridging(list_of_ssDNAs, linearized_plasmid)

                # Construct meaningful names, IDs, and descriptions
                targeting_info = [f"pCRISPR-BEST_{row['locus_tag']}_p{row['sgrna_loc']}" for _, row in filtered_df.iterrows()]
                for i, vector in enumerate(sgRNA_vectors):
                    vector.name = f'{targeting_info[i]}_#{i+1}'
                    vector.id = vector.name
                    vector.description = f'Assembled plasmid targeting {", ".join(genes_to_KO_list)} for base-editing, assembled using StreptoCAD.'

                # Annotate plasmids
                for plasmid in sgRNA_vectors: 
                    annotate_plasmid_with_sgrnas(plasmid, filtered_df)

                # Extract metadata for plasmid
                integration_names = filtered_df.apply(lambda row: f"sgRNA_{row['locus_tag']}({row['sgrna_loc']})", axis=1).tolist()
                plasmid_metadata_df = extract_metadata_to_dataframe(sgRNA_vectors, clean_plasmid, integration_names)
                
                # Generate primers
                idt_df1 = primers_to_IDT(list_of_ssDNAs)
                checking_primers_df = find_best_check_primers_from_genome(genome, genes_to_KO_list, flanking_region=flanking_region_number, target_tm=melting_temperature, primer_concentration=primer_concentration, polymerase=chosen_polymerase)
                idt_df2 = create_idt_order_dataframe(checking_primers_df)
                full_idt = pd.concat([idt_df1, idt_df2])

                # Prepare outputs for the DataTable
                primers_columns = [{"name": col, "id": col} for col in full_idt.columns]
                primers_data = full_idt.to_dict('records')

                pcr_columns = [{"name": col, "id": col} for col in checking_primers_df.columns]
                pcr_data = checking_primers_df.to_dict('records')

                # Generate download link for filtered_df
                filtered_df_columns = [{"name": col, "id": col} for col in filtered_df.columns]
                filtered_df_data = filtered_df.to_dict('records')

                # Generate project directory structure
                input_files = [
                    {"name": "input_genome.gb", "content": genome},
                    {"name": "input_plasmid.gb", "content": clean_plasmid}
                ]
                output_files = [
                    {"name": "cBEST_w_sgRNAs.gb", "content": sgRNA_vectors},
                    {"name": "primer_df.csv", "content": checking_primers_df},
                    {"name": "full_idt.csv", "content": full_idt},
                    {"name": "mutated_sgrna_df.csv", "content": mutated_sgrna_df},
                    {"name": "filtered_df.csv", "content": filtered_df},
                    {"name": "plasmid_metadata_df.csv", "content": plasmid_metadata_df},
                ]
                input_values = {
                    "genes_to_knockout": genes_to_KO_list,
                    "polymerase_settings": {
                        "chosen_polymerase": chosen_polymerase,
                        "melting_temperature": melting_temperature,
                        "primer_concentration": primer_concentration,
                        "primer_number_increment": primer_number_increment,
                        "flanking_region_number": flanking_region_number
                    },
                    "filtering_metrics": {
                        "gc_upper": gc_upper,
                        "gc_lower": gc_lower,
                        "off_target_seed": off_target_seed,
                        "off_target_upper": off_target_upper,
                        "cas_type": cas_type,
                        "number_of_sgRNAs_per_group": number_of_sgRNAs_per_group
                    },
                    "overlapping_sequences": {
                        "up_homology": str(up_homology),
                        "dw_homology": str(dw_homology)
                    }
                }
                markdown_file_paths = ["../protocols/conjugation_protcol.md", "../protocols/single_target_crispr_plasmid_protcol.md"]
                timestamp = datetime.utcnow().isoformat()

                project_directory = ProjectDirectory(
                    project_name=f"CRISPR_cBEST_workflow_{timestamp}",
                    input_files=input_files,
                    output_files=output_files,
                    input_values=input_values,
                    markdown_file_paths=markdown_file_paths
                )
                zip_content = project_directory.create_directory_structure(create_directories=True)
                data_package_encoded = base64.b64encode(zip_content).decode('utf-8')
                data_package_download_link = f"data:application/zip;base64,{data_package_encoded}"

                logging.info("Workflow 2 completed successfully")

                # Clear the log stream after successful execution
                log_stream.truncate(0)
                log_stream.seek(0)
                
                # Prepare columns and data for the plasmid metadata DataTable
                plasmid_metadata_columns = [{"name": col, "id": col} for col in plasmid_metadata_df.columns]
                plasmid_metadata_data = plasmid_metadata_df.to_dict('records')

            return (
                primers_data, 
                primers_columns, 
                pcr_data, 
                pcr_columns, 
                data_package_download_link, 
                filtered_df_data, 
                filtered_df_columns, 
                "", 
                False, 
                plasmid_metadata_data,  # Data for the new plasmid metadata DataTable
                plasmid_metadata_columns # Columns for the new plasmid metadata DataTable
            )
        
        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
            error_message = f"An error occurred: {str(e)}\n\nLog:\n{log_stream.getvalue()}"
            display_error = True
            return [], [], [], [], "", [], [], error_message, display_error, [], []
