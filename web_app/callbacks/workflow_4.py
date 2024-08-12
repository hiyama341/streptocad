import sys
import os
import io
import zipfile
import base64
import csv
import logging
from datetime import datetime

import pandas as pd
from Bio import SeqIO
from Bio.Restriction import NcoI
from pydna.dseqrecord import Dseqrecord
from teemi.design.fetch_sequences import read_genbank_files

from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import quote

import tempfile

module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.sequence_loading.sequence_loading import (
    load_and_process_gene_sequences, 
    load_and_process_plasmid, 
    load_and_process_genome_sequences,
    check_and_convert_input,
    annotate_dseqrecord,
    process_specified_gene_sequences_from_record
)
from streptocad.utils import ProjectDirectory, extract_metadata_to_dataframe
from streptocad.crispr.guideRNA_crispri import extract_sgRNAs_for_crispri, SgRNAargs
from streptocad.cloning.ssDNA_bridging import assemble_plasmids_by_ssDNA_bridging, make_ssDNA_oligos
from streptocad.primers.primer_generation import create_idt_order_dataframe, primers_to_IDT

# Setup logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)s %(message)s',
                    handlers=[
                        logging.FileHandler("workflow4_debug.log"),
                        logging.StreamHandler()
                    ])

def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(name, "wb") as fp:
        fp.write(base64.decodebytes(data))

def register_workflow_4_callbacks(app):

    @app.callback(
        [
            Output('primer-table_4', 'data'),
            Output('primer-table_4', 'columns'),
            Output('download-data-and-protocols-link_4', 'href'),
            Output('mutated-sgrna-table_4', 'data'),
            Output('mutated-sgrna-table_4', 'columns'),
        ],
        [
            Input('submit-settings-button_4', 'n_clicks')
        ],
        [
            State('upload-genome-file_4', 'contents'),
            State('upload-crispri-vector_4', 'contents'),
            State('upload-genome-file_4', 'filename'),
            State('upload-crispri-vector_4', 'filename'),
            State('genes-to-KO_4', 'value'),
            State('forward-overhang-input_4', 'value'),
            State('reverse-overhang-input_4', 'value'),
            State('gc-upper_4', 'value'),
            State('gc-lower_4', 'value'),
            State('off-target-seed_4', 'value'),
            State('off-target-upper_4', 'value'),
            State('cas-type_4', 'value'),
            State('number-of-sgRNAs-per-group_4', 'value'),
            State('extension-to-promoter-region_4', 'value')
        ]
    )
    def run_workflow(n_clicks, genome_content, vector_content, genome_filename, vector_filename, genes_to_KO, 
                     up_homology, dw_homology, gc_upper, gc_lower, off_target_seed, off_target_upper, cas_type, 
                     number_of_sgRNAs_per_group, extension_to_promoter_region):
        if n_clicks is None:
            raise PreventUpdate

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
                target_dict, genes_to_KO_list, annotation_input = check_and_convert_input(genes_to_KO)
                if annotation_input:
                    genome = annotate_dseqrecord(genome, target_dict)
                logging.info(f"Genes to knock out: {genes_to_KO_list}")

                # Extract sgRNAs
                args = SgRNAargs(
                    genome,
                    genes_to_KO_list,
                    step=['find', 'filter'],
                    gc_upper=gc_upper,
                    gc_lower=gc_lower,
                    off_target_seed=off_target_seed,
                    off_target_upper=off_target_upper,
                    cas_type=cas_type,
                    extension_to_promoter_region=extension_to_promoter_region,
                    target_non_template_strand=True
                )

                logging.info("Extracting sgRNAs for CRISPRi")
                sgrna_df = extract_sgRNAs_for_crispri(args)
                logging.debug(f"sgRNA DataFrame: {sgrna_df}")

                filtered_df = sgrna_df.groupby('locus_tag').head(number_of_sgRNAs_per_group)
                logging.debug(f"Filtered DataFrame: {filtered_df}")

                logging.info("Making ssDNA oligos")
                list_of_ssDNAs = make_ssDNA_oligos(
                    filtered_df, 
                    upstream_ovh=Dseqrecord(up_homology),
                    downstream_ovh=Dseqrecord(dw_homology)
                )
                logging.debug(f"List of ssDNAs: {list_of_ssDNAs}")

                logging.info("Cutting plasmid")
                linearized_plasmid = sorted(clean_plasmid.cut(NcoI), key=lambda x: len(x), reverse=True)[0]
                logging.debug(f"Linearized plasmid: {linearized_plasmid}")

                logging.info("Assembling plasmid")
                sgRNA_vectors = assemble_plasmids_by_ssDNA_bridging(list_of_ssDNAs, linearized_plasmid)
                logging.debug(f"sgRNA vectors: {sgRNA_vectors}")

                targeting_info = []
                for index, row in filtered_df.iterrows():
                    formatted_str = f"CRISPRi_{row['locus_tag']}_p{row['sgrna_loc']}"
                    targeting_info.append(formatted_str)

                for i in range(len(sgRNA_vectors)):
                    sgRNA_vectors[i].name = f'p{targeting_info[i]}_#{i+1}'
                    sgRNA_vectors[i].id = sgRNA_vectors[i].name
                    sgRNA_vectors[i].description = f'Assembled plasmid targeting {", ".join(genes_to_KO_list)} for single gene KNOCK-DOWN, assembled using StreptoCAD.'

                logging.info("Generating primers for IDT")
                idt_primers = primers_to_IDT(list_of_ssDNAs)
                logging.debug(f"IDT primers DataFrame: {idt_primers}")

                # Prepare DataTables outputs
                primers_columns = [{"name": col, "id": col} for col in idt_primers.columns]
                primers_data = idt_primers.to_dict('records')

                # Prepare download link for the data package
                input_files = [
                    {"name": "input_genome.gb", "content": genome},
                    {"name": "input_plasmid.gb", "content": clean_plasmid}
                ]

                output_files = [
                    {"name": "cBEST_w_sgRNAs.gb", "content": sgRNA_vectors},
                    {"name": "full_idt.csv", "content": idt_primers},
                    {"name": "sgrna_df.csv", "content": sgrna_df},
                    {"name": "filtered_df.csv", "content": filtered_df}
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
                        'extension_to_promoter_region': extension_to_promoter_region,
                    },
                    "overlapping_sequences": {
                        "up_homology": str(up_homology),
                        "dw_homology": str(dw_homology)
                    }
                }

                markdown_file_paths = [
                    "../protocols/conjugation_protcol.md",
                    "../protocols/single_target_crispr_plasmid_protcol.md"
                ]

                timestamp = datetime.utcnow().isoformat()

                logging.info("Creating project directory structure")
                project_directory = ProjectDirectory(
                    project_name=f"CRISPRi_workflow_{timestamp}",
                    input_files=input_files,
                    output_files=output_files,
                    input_values=input_values,
                    markdown_file_paths=markdown_file_paths
                )

                zip_content = project_directory.create_directory_structure(create_directories=True)
                data_package_encoded = base64.b64encode(zip_content).decode('utf-8')
                data_package_download_link = f"data:application/zip;base64,{data_package_encoded}"

                filtered_df_columns = [{"name": col, "id": col} for col in filtered_df.columns]
                filtered_df_data = filtered_df.to_dict('records')

                logging.info("Workflow completed successfully")

                return (primers_data, primers_columns, data_package_download_link, filtered_df_data, filtered_df_columns)

        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
            raise PreventUpdate
