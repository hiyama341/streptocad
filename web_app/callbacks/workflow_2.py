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
from Bio.Restriction import NcoI, StuI
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord
from teemi.design.fetch_sequences import read_genbank_files

import dash
from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
from dash.dash_table.Format import Group
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import quote
import logging
import tempfile

# Set up logging to use a StringIO buffer
log_stream = io.StringIO()
logging.basicConfig(stream=log_stream, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Local module imports
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.primers.primer_generation import primers_to_IDT
from streptocad.cloning.ssDNA_bridging import assemble_plasmids_by_ssDNA_bridging, make_ssDNA_oligos
from streptocad.utils import ProjectDirectory

from streptocad.sequence_loading.sequence_loading import load_and_process_gene_sequences, load_and_process_plasmid, load_and_process_genome_sequences
from streptocad.utils import generate_project_directory_structure
from streptocad.crispr.guideRNAcas3_9_12 import extract_sgRNAs, SgRNAargs
from streptocad.crispr.crispr_best import identify_base_editing_sites, filter_sgrnas_for_base_editing, process_base_editing
from streptocad.cloning.plasmid_processing import annotate_plasmid_with_sgrnas
from streptocad.primers.primer_generation import checking_primers, create_idt_order_dataframe, primers_to_IDT

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
            Output('genbank-file-single_2', 'href'),
            Output('primers_download_link_2', 'href'),
            Output('download-pcr-link_2', 'href'),
            Output('download-data-and-protocols-link_2', 'href'),
            Output('filtered-df-table', 'data'),   # New output for filtered_df
            Output('filtered-df-table', 'columns'), # New output for filtered_df columns
            Output('download-filtered-df-link_2', 'href'), # New output for filtered_df download link
            Output('error-dialog_2', 'message'),
            Output('error-dialog_2', 'displayed')
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
                input_genome = list(SeqIO.parse(genome_path, "genbank"))
                input_plasmid = list(SeqIO.parse(vector_path, "genbank"))

                if not input_genome or not input_plasmid:
                    raise ValueError("Unsupported file format. Please provide a valid GenBank file.")

                genome = Dseqrecord(input_genome[0].seq, id=input_genome[0].id, description=input_genome[0].description)
                plasmid = Dseqrecord(input_plasmid[0], circular=True)

                genes_to_KO_list = [gene.strip() for gene in genes_to_KO.split(',')]
                logging.info(f"Genes to knock out: {genes_to_KO_list}")

                # Initialize SgRNAargs with desired parameters
                args = SgRNAargs(genome_path, 
                                 genes_to_KO_list,
                                 step=['find', 'filter'],
                                 gc_upper=gc_upper,
                                 gc_lower=gc_lower,
                                 off_target_seed=off_target_seed,
                                 off_target_upper=off_target_upper,
                                 cas_type=cas_type)

                logging.info("Extracting sgRNAs")
                sgrna_df = extract_sgRNAs(args)
                logging.info(f"sgRNA DataFrame: {sgrna_df}")

                # Load gene sequences
                logging.info("Loading gene sequences")
                gene_sequences = load_and_process_gene_sequences(genome_path)
                genes_to_KO_dict = {locus_tag: gene_sequences[locus_tag] for locus_tag in genes_to_KO_list if locus_tag in gene_sequences}
                logging.info(f"Genes to KO dictionary: {genes_to_KO_dict}")

                # Identify and annotate base editing sites
                logging.info("Identifying base editing sites")
                sgrna_df_with_editing = identify_base_editing_sites(sgrna_df)

                # Filter out only sgRNAs that result in base-editing
                logging.info("Filtering sgRNAs for base editing")
                filtered_sgrna_df_for_base_editing = filter_sgrnas_for_base_editing(sgrna_df_with_editing)
                logging.info(f"Filtered sgRNAs for base editing: {filtered_sgrna_df_for_base_editing}")

                # Process the DataFrame to apply C-to-T mutations
                logging.info("Processing base editing")
                mutated_sgrna_df = process_base_editing(filtered_sgrna_df_for_base_editing, 
                                                        genes_to_KO_dict, 
                                                        only_stop_codons=bool(only_stop_codons))
                logging.info(f"Mutated sgRNAs: {mutated_sgrna_df}")

                # Filter the DataFrame to retain only up to 5 sgRNA sequences per locus_tag
                logging.info("Filtering sgRNAs to retain only up to 5 sequences per locus tag")
                filtered_df = mutated_sgrna_df.groupby('locus_tag').head(number_of_sgRNAs_per_group)
                logging.info(f"Filtered DataFrame: {filtered_df}")

                # Make oligos
                logging.info("Making ssDNA oligos")
                list_of_ssDNAs = make_ssDNA_oligos(filtered_df, upstream_ovh=Dseqrecord(up_homology),
                                                   downstream_ovh=Dseqrecord(dw_homology))
                logging.info(f"List of ssDNAs: {list_of_ssDNAs}")

                # Cut plasmid
                logging.info("Cutting plasmid")
                linearized_plasmid = sorted(plasmid.cut(NcoI), key=lambda x: len(x), reverse=True)[0]
                logging.info(f"Linearized plasmid: {linearized_plasmid}")

                # Assemble plasmid
                logging.info("Assembling plasmid")
                sgRNA_vectors = assemble_plasmids_by_ssDNA_bridging(list_of_ssDNAs, linearized_plasmid)
                logging.info(f"sgRNA vectors: {sgRNA_vectors}")

                # Constructing a meaningful name, ID, and description for the assembled plasmid using user input
                targeting_info = []
                for index, row in filtered_df.iterrows():
                    formatted_str = f"pCRISPR-BEST_{row['locus_tag']}_p{row['sgrna_loc']}"
                    targeting_info.append(formatted_str)

                for i in range(len(sgRNA_vectors)):
                    sgRNA_vectors[i].name = f'{targeting_info[i]}_#{i+1}'
                    sgRNA_vectors[i].id = sgRNA_vectors[i].name  # Using the same value for ID as for name for simplicity
                    sgRNA_vectors[i].description = f'Assembled plasmid targeting {", ".join(genes_to_KO_list)} for base-editing, assembled using StreptoCAD.'

                logging.info("Annotating plasmids")
                # Annotate plasmids
                for plasmid in sgRNA_vectors: 
                    annotate_plasmid_with_sgrnas(plasmid, filtered_df)

                # Generate primers
                logging.info("Generating primers for IDT")
                idt_df1 = primers_to_IDT(list_of_ssDNAs)
                logging.info(f"IDT primers DataFrame: {idt_df1}")

                # Getting checking primers
                logging.info("Generating checking primers")
                checking_primers_df = checking_primers(genome_path, genes_to_KO_list, 
                                                       flanking_region=flanking_region_number,
                                                       target_tm=melting_temperature, 
                                                       primer_concentration=primer_concentration, 
                                                       polymerase=chosen_polymerase)
                logging.info(f"Checking primers DataFrame: {checking_primers_df}")

                idt_df2 = create_idt_order_dataframe(checking_primers_df)
                logging.info(f"IDT order DataFrame: {idt_df2}")
                full_idt = pd.concat([idt_df1, idt_df2])
                logging.info(f"Full IDT DataFrame: {full_idt}")

                # Prepare outputs for the DataTable
                logging.info("Preparing data for DataTables")
                primers_columns = [{"name": col, "id": col} for col in full_idt.columns]
                primers_data = full_idt.to_dict('records')

                pcr_columns = [{"name": col, "id": col} for col in checking_primers_df.columns]
                pcr_data = checking_primers_df.to_dict('records')

                # IDT 
                pcr_df_string = checking_primers_df.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
                pcr_data_encoded = quote(pcr_df_string)
                pcr_download_link = f"data:text/csv;charset=utf-8,{pcr_data_encoded}"            
                
                # PCR-primer-df
                primer_df_string = full_idt.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
                primer_data_encoded = quote(primer_df_string)
                primer_download_link = f"data:text/csv;charset=utf-8,{primer_data_encoded}" 

                # Provide a link for downloading GenBank files:
                encoded_genbank_files = [base64.b64encode(str(vector).encode()).decode() for vector in sgRNA_vectors]
                genbank_download_link = f"data:application/genbank;base64,{encoded_genbank_files[0]}"
                # Create a zip archive in memory
                zip_buffer = io.BytesIO()
                zip_data = None  # Initialize zip_data here
                with zipfile.ZipFile(zip_buffer, 'a', zipfile.ZIP_DEFLATED, False) as zip_file:
                    for vector in sgRNA_vectors:
                        # Convert each vector to GenBank format
                        genbank_content = vector.format("genbank")
                        zip_file.writestr(f"{vector.name}.gb", genbank_content)

                # Prepare the zip archive for download
                zip_buffer.seek(0)
                zip_data = base64.b64encode(zip_buffer.read()).decode('utf-8')
                genbank_download_link = f"data:application/zip;base64,{zip_data}"


                # New function call to generate project directory structure
                input_files = [
                    {"name": "input_genome.gb", "content": genome},
                    {"name": "input_plasmid.gb", "content": plasmid}
                ]

                output_files = [
                    {"name": "cBEST_w_sgRNAs.gb", "content": sgRNA_vectors}, # LIST OF Dseqrecords
                    {"name": "primer_df.csv", "content": checking_primers_df},
                    {"name": "full_idt.csv", "content": full_idt},
                    {"name": "mutated_sgrna_df.csv", "content": mutated_sgrna_df},
                    {"name": "filtered_df.csv", "content": filtered_df}
                ]

                input_values = {
                    "genes_to_knockout": genes_to_KO,
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


                # Paths to Markdown files
                markdown_file_paths = [
                    "../protocols/conjugation_protcol.md",
                    "../protocols/single_target_crispr_plasmid_protcol.md"

                ]

                # Data and time
                timestamp = datetime.utcnow().isoformat()

                # Create project directory structure
                project_directory = ProjectDirectory(
                    project_name=f"CRISPR_cBEST_workflow_{timestamp}",
                    input_files=input_files,
                    output_files=output_files,
                    input_values=input_values,
                    markdown_file_paths=markdown_file_paths
                )

                # Generate the project directory structure and get the zip content
                zip_content = project_directory.create_directory_structure(create_directories=True)

                # Encode the zip content for download
                data_package_encoded = base64.b64encode(zip_content).decode('utf-8')
                data_package_download_link = f"data:application/zip;base64,{data_package_encoded}"

                # Prepare data for filtered_df
                filtered_df_columns = [{"name": col, "id": col} for col in filtered_df.columns]
                filtered_df_data = filtered_df.to_dict('records')
                
                # Generate download link for filtered_df
                filtered_df_string = filtered_df.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
                filtered_df_data_encoded = quote(filtered_df_string)
                filtered_df_download_link = f"data:text/csv;charset=utf-8,{filtered_df_data_encoded}"

                logging.info("Workflow 2 completed successfully")

                # Clear the log stream after successful execution
                log_stream.truncate(0)
                log_stream.seek(0)

                return (primers_data, primers_columns, pcr_data, pcr_columns, genbank_download_link, primer_download_link, 
                        pcr_download_link, data_package_download_link, filtered_df_data, filtered_df_columns, filtered_df_download_link, "", False)

        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
            error_message = f"An error occurred: {str(e)}\n\nLog:\n{log_stream.getvalue()}"
            return [], [], [], [], "", "", "", "", [], [], "", error_message, True
