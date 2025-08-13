import sys
import os
import io
import zipfile
import base64
import tempfile
import logging
import pandas as pd
from datetime import datetime
from Bio import SeqIO
from Bio.Restriction import StuI
from pydna.dseqrecord import Dseqrecord
from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

from dash.exceptions import PreventUpdate
from Bio.Restriction import *
from Bio import Restriction


# Local module imports
module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.sequence_loading.sequence_loading import (
    load_and_process_gene_sequences,
    load_and_process_genome_sequences,
    load_and_process_plasmid,
    check_and_convert_input,
    annotate_dseqrecord,
    process_specified_gene_sequences_from_record,
)

from streptocad.utils import (
    polymerase_dict,
    create_primer_df_from_dict,
    ProjectDirectory,
    extract_metadata_to_dataframe,
)
from streptocad.primers.primer_generation import create_idt_order_dataframe
from streptocad.cloning.ssDNA_bridging import (
    assemble_plasmids_by_ssDNA_bridging,
    make_ssDNA_oligos,
)
from streptocad.crispr.guideRNAcas3_9_12 import extract_sgRNAs, SgRNAargs
from streptocad.cloning.gibson_cloning import (
    find_up_dw_repair_templates,
    assemble_multiple_plasmids_with_repair_templates_for_deletion,
    update_primer_names,
)

from streptocad.cloning.plasmid_processing import (
    check_plasmid_restriction_sites,
    determine_workflow_order_for_plasmids,
)
from streptocad.primers.primer_generation import (
    checking_primers,
    primers_to_IDT,
    find_best_check_primers_from_genome,
)

# Logging setup similar to Workflow 2
log_stream = io.StringIO()
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout), logging.StreamHandler(log_stream)],
)
logger = logging.getLogger(__name__)


def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(name, "wb") as fp:
        fp.write(base64.decodebytes(data))


def register_workflow_5_callbacks(app):
    @app.callback(
        [
            Output("primers-output-table_5", "data"),
            Output("primers-output-table_5", "columns"),
            Output("pcr-table_5", "data"),
            Output("pcr-table_5", "columns"),
            Output("download-data-and-protocols-link_5", "href"),
            Output("filtered-df-table_5", "data"),
            Output("filtered-df-table_5", "columns"),
            Output("error-dialog_5", "message"),
            Output("error-dialog_5", "displayed"),
            Output("plasmid-metadata-table_5", "data"),
            Output("plasmid-metadata-table_5", "columns"),
        ],
        [Input("submit-settings-button_5", "n_clicks")],
        [
            State({"type": "upload-component", "index": "genome-file-5"}, "contents"),
            State({"type": "upload-component", "index": "single-vector-5"}, "contents"),
            State({"type": "upload-component", "index": "genome-file-5"}, "filename"),
            State({"type": "upload-component", "index": "single-vector-5"}, "filename"),
            State("genes-to-KO_5", "value"),
            State("forward-overhang-input_5", "value"),
            State("reverse-overhang-input_5", "value"),
            State("gc-upper_5", "value"),
            State("gc-lower_5", "value"),
            State("off-target-seed_5", "value"),
            State("off-target-upper_5", "value"),
            State("cas-type_5", "value"),
            State("number-of-sgRNAs-per-group_5", "value"),
            State("show-inframe-deletions-settings-checkbox_5", "value"),
            State("chosen-polymerase_5", "value"),
            State("melting-temperature_5", "value"),
            State("primer-concentration_5", "value"),
            State("flanking-region-number_5", "value"),
            State("repair_templates_length_5", "value"),
            State("overlap_for_gibson_length_5", "value"),
            State("restriction-enzymes_5", "value"),
            State("restriction-enzymes_5-2", "value"),
            State("checking-primer-length_5", "value"),
            State("gibson-primer-length_5", "value"),
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
        in_frame_deletion,
        chosen_polymerase,
        melting_temperature,
        primer_concentration,
        flanking_region_number,
        repair_templates_length,
        overlap_for_gibson_length,
        restriction_enzymes,
        restriction_enzymes_2,
        checking_primer_length,
        gibson_primer_length,
    ):
        if n_clicks is None:
            raise PreventUpdate

        try:
            logging.info("Workflow 5 started")
            print(f"########  TEST OF INFRAME DELETION: {in_frame_deletion}########")

            # Initialize variables to avoid UnboundLocalError
            primer_df = pd.DataFrame()
            idt_df = pd.DataFrame()

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

                # Filter the DataFrame to retain only up to 5 sgRNA sequences per locus_tag
                logging.info("Filtering sgRNAs ")
                filtered_df = sgrna_df.groupby("locus_tag").head(
                    number_of_sgRNAs_per_group
                )

                # MAke oligoes
                list_of_ssDNAs = make_ssDNA_oligos(
                    filtered_df, upstream_ovh=up_homology, downstream_ovh=dw_homology
                )

                # cut plasmid
                restriction_enzymes = restriction_enzymes.split(",")
                enzymes_for_repair_template_integration = [
                    getattr(Restriction, str(enzyme)) for enzyme in restriction_enzymes
                ]
                linearized_plasmid = sorted(
                    clean_plasmid.cut(enzymes_for_repair_template_integration),
                    key=lambda x: len(x),
                    reverse=True,
                )[0]
                # print(linearized_plasmid)

                sgRNA_vectors = assemble_plasmids_by_ssDNA_bridging(
                    list_of_ssDNAs, linearized_plasmid
                )

                # Constructing a meaningful name, ID, and description for the assembled plasmid using user input
                targeting_info = []
                for index, row in filtered_df.iterrows():
                    formatted_str = f"pCas9_{row['locus_tag']}({row['sgrna_loc']})"
                    targeting_info.append(formatted_str)

                for i in range(len(sgRNA_vectors)):
                    sgRNA_vectors[i].name = f"{targeting_info[i]}"
                    sgRNA_vectors[i].id = sgRNA_vectors[i].name
                    sgRNA_vectors[
                        i
                    ].description = f"CRISPR-Cas9 targeting {', '.join(genes_to_KO)} for single gene knockout, assembled using StreptoCAD."

                logging.info("Processing idt primers")
                idt_primers = primers_to_IDT(list_of_ssDNAs)

                # **Check if in-frame deletion is enabled**
                if in_frame_deletion == [1]:
                    # Make repair templates
                    repair_templates_data = find_up_dw_repair_templates(
                        genome,
                        genes_to_KO_list,
                        target_tm=melting_temperature,
                        primer_tm_kwargs={
                            "conc": primer_concentration,
                            "prodcode": chosen_polymerase,
                        },
                        repair_length=repair_templates_length,
                        min_primer_length=gibson_primer_length,
                    )
                    logging.info(f"Repair templates data: {repair_templates_data}")

                    # Digest the plasmids
                    restriction_enzymes_2 = restriction_enzymes_2.split(",")
                    enzymes_for_repair_template_integration_2 = [
                        getattr(Restriction, str(enzyme))
                        for enzyme in restriction_enzymes_2
                    ]
                    processed_records = [
                        Dseqrecord(record, circular=True).cut(
                            enzymes_for_repair_template_integration_2
                        )[0]
                        for record in sgRNA_vectors
                    ]
                    logging.info(f"processed_records: {processed_records}")

                    # Rename them appropriately
                    for i in range(len(processed_records)):
                        processed_records[i].name = sgRNA_vectors[i].name
                        logging.info(
                            f"processed_records names: {processed_records[i].name}"
                        )

                    logging.info(f"processed_records names: {processed_records}")

                    # Assembly
                    logging.info(f"genes_to_KO: {genes_to_KO}")
                    logging.info(f"processed_records: {processed_records}")
                    logging.info(f"repair_templates_data: {repair_templates_data}")
                    logging.info(
                        f"overlap_for_gibson_length: {overlap_for_gibson_length}"
                    )

                    # Assembly
                    assembly_data = (
                        assemble_multiple_plasmids_with_repair_templates_for_deletion(
                            genes_to_KO_list,
                            processed_records,
                            repair_templates_data,
                            overlap=overlap_for_gibson_length,
                        )
                    )
                    logging.info(f"assembly_data: {assembly_data}")

                    # updating the primer names to something systematic.
                    update_primer_names(assembly_data)

                    # Parse through the primer df
                    primer_df = create_primer_df_from_dict(assembly_data)
                    unique_df = primer_df.drop_duplicates(keep="first")
                    logging.info(f"unique_df: {unique_df}")

                    # IDT df
                    idt_df = create_idt_order_dataframe(
                        unique_df, concentration="25nm", purification="STD"
                    )

                    # Contigs
                    assembled_contigs = []
                    for data in assembly_data:
                        contig_record = data["contig"]
                        contig_record.id = f"{data['name']}_w_rep"
                        contig_record.name = f"{data['name']}_w_rep"
                        assembled_contigs.append(contig_record)

                # outside in-frame deletion condition
                logging.info("Plasmid metadata df")
                if 1 in in_frame_deletion:
                    integration_names = filtered_df.apply(
                        lambda row: f"sgRNA_{row['locus_tag']}({row['sgrna_loc']})",
                        axis=1,
                    ).tolist()
                    plasmid_metadata_df = extract_metadata_to_dataframe(
                        assembled_contigs, clean_plasmid, integration_names
                    )

                    workflow_df = determine_workflow_order_for_plasmids(
                        sgRNA_vectors,
                        assembled_contigs,
                        restriction_enzymes_2,
                        restriction_enzymes,
                    )

                    plasmid_metadata_df = pd.merge(
                        plasmid_metadata_df, workflow_df, on="plasmid_name", how="inner"
                    )

                else:
                    integration_names = filtered_df.apply(
                        lambda row: f"sgRNA_{row['locus_tag']}({row['sgrna_loc']})",
                        axis=1,
                    ).tolist()
                    plasmid_metadata_df = extract_metadata_to_dataframe(
                        sgRNA_vectors, clean_plasmid, integration_names
                    )

                logging.info("Checking primers")
                checking_primers_df = find_best_check_primers_from_genome(
                    genome,
                    genes_to_KO_list,
                    flanking_region=flanking_region_number,
                    target_tm=melting_temperature,
                    primer_concentration=primer_concentration,
                    polymerase=chosen_polymerase,
                    limit=checking_primer_length,
                )

                logging.info(
                    "Creating a copy of checking primers DataFrame and renaming columns."
                )
                checking_primers_df_copy = checking_primers_df.copy()
                checking_primers_df_copy = checking_primers_df_copy.rename(
                    columns={"locus tag": "template"}
                )
                checking_primers_df_copy = checking_primers_df_copy.loc[
                    :, ~checking_primers_df_copy.columns.duplicated(keep="first")
                ]

                # Build IDT input from the cleaned/renamed DF (not the original)
                checking_primers_df_idt = create_idt_order_dataframe(
                    checking_primers_df_copy
                )
                checking_primers_df_idt = checking_primers_df_idt.loc[
                    :, ~checking_primers_df_idt.columns.duplicated(keep="first")
                ]

                logging.info("Making final IDT order sheet")
                if 1 in in_frame_deletion:
                    # De-dup columns on all frames going into concat
                    idt_primers = idt_primers.loc[
                        :, ~idt_primers.columns.duplicated(keep="first")
                    ]
                    idt_df = idt_df.loc[:, ~idt_df.columns.duplicated(keep="first")]

                    full_idt = pd.concat(
                        [idt_primers, idt_df, checking_primers_df_idt],
                        ignore_index=True,
                    )

                    # pcr_table
                    unique_df = unique_df.loc[
                        :, ~unique_df.columns.duplicated(keep="first")
                    ]
                    pcr_table = pd.concat(
                        [unique_df, checking_primers_df_copy], ignore_index=True
                    )
                else:
                    # Ensure clean columns before concat
                    idt_primers = idt_primers.loc[
                        :, ~idt_primers.columns.duplicated(keep="first")
                    ]
                    checking_primers_df_idt = checking_primers_df_idt.loc[
                        :, ~checking_primers_df_idt.columns.duplicated(keep="first")
                    ]
                    full_idt = pd.concat(
                        [idt_primers, checking_primers_df_idt], ignore_index=True
                    )

                    # For the else branch you previously returned the raw copy; keep it clean too
                    pcr_table = checking_primers_df_copy

                input_files = [
                    {"name": "input_genome.gb", "content": genome},
                    {"name": "input_plasmid.gb", "content": clean_plasmid},
                ]

                if in_frame_deletion:
                    output_files = [
                        {
                            "name": "Cas9_w_sgRNAs.gb",
                            "content": assembled_contigs,
                        },  # LIST OF Dseqrecords
                        {"name": "01_primer_df.csv", "content": unique_df},
                        {"name": "02_full_idt.csv", "content": full_idt},
                        {"name": "03_sgrna_df.csv", "content": sgrna_df},
                        {"name": "04_filtered_sgrna_df.csv", "content": filtered_df},
                        {
                            "name": "05_plasmid_metadata_df.csv",
                            "content": plasmid_metadata_df,
                        },
                        {"name": "06_workflow_order_df.csv", "content": workflow_df},
                    ]
                else:
                    output_files = [
                        {
                            "name": "Cas9_sgRNAs.gb",
                            "content": sgRNA_vectors,
                        },  # LIST OF Dseqrecords
                        {"name": "01_full_idt.csv", "content": full_idt},
                        {"name": "02_sgrna_df.csv", "content": sgrna_df},
                        {"name": "03_filtered_sgrna_df.csv", "content": filtered_df},
                        {
                            "name": "04_plasmid_metadata_df.csv",
                            "content": plasmid_metadata_df,
                        },
                    ]

                input_values = {
                    "genes_to_knockout": genes_to_KO,
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
                        "up_homology": str(up_homology),
                        "dw_homology": str(dw_homology),
                    },
                }

                # Paths to Markdown files
                markdown_file_paths = [
                    "protocols/conjugation_protcol.md",
                    "protocols/single_target_crispr_plasmid_protcol.md",
                    "protocols/trouble_shooting_tips.md",
                ]

                timestamp = datetime.utcnow().isoformat()

                # Create project directory structure
                project_directory = ProjectDirectory(
                    project_name=f"CRISPR_cas9_plasmid_workflow_{timestamp}",
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

                # Prepare outputs for the DataTable
                primer_columns = [{"name": col, "id": col} for col in full_idt.columns]
                primer_data = full_idt.to_dict("records")

                pcr_columns = [{"name": col, "id": col} for col in pcr_table.columns]
                pcr_data = pcr_table.to_dict("records")

                plasmid_metadata_df_columns = [
                    {"name": col, "id": col} for col in plasmid_metadata_df.columns
                ]
                plasmid_metadata_df_data = plasmid_metadata_df.to_dict("records")

                # Generate download link for filtered_df
                filtered_df_columns = [
                    {"name": col, "id": col} for col in filtered_df.columns
                ]
                filtered_df_data = filtered_df.to_dict("records")

                logging.info("Workflow 5 completed successfully")

                # Clear the log stream after successful execution
                log_stream.truncate(0)
                log_stream.seek(0)

            return (
                primer_data,  # primers-output-table_5.data
                primer_columns,  # primers-output-table_5.columns
                pcr_data,  # pcr-table_5.data
                pcr_columns,  # pcr-table_5.columns
                data_package_download_link,  # download-data-and-protocols-link_5.href
                filtered_df_data,  # filtered-df-table_5.data
                filtered_df_columns,  # filtered-df-table_5.columns
                "",  # error-dialog_5.message (empty if no error)
                False,  # error-dialog_5.displayed (set to False if no error)
                plasmid_metadata_df_data,  # plasmid-metadata-table_5.data
                plasmid_metadata_df_columns,  # plasmid-metadata-table_5.columns
            )

        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
            error_message = (
                f"An error occurred: {str(e)}\n\nLog:\n{log_stream.getvalue()}"
            )
            return (
                [],  # primers-output-table_5.data
                [],  # primers-output-table_5.columns
                [],  # pcr-table_5.data
                [],  # pcr-table_5.columns
                "",  # download-data-and-protocols-link_5.href
                [],  # filtered-df-table_5.data
                [],  # filtered-df-table_5.columns
                error_message,  # error-dialog_5.message
                True,  # error-dialog_5.displayed
                [],  # plasmid-metadata-table_5.data
                [],  # plasmid-metadata-table_5.columns
            )
