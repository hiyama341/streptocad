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
from Bio.Restriction import *  # we import all enzymes
from Bio import Restriction

from pydna.dseqrecord import Dseqrecord
from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc

from dash.exceptions import PreventUpdate


# Local module imports
module_path = os.path.abspath(os.path.join(".."))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.sequence_loading.sequence_loading import (
    load_and_process_plasmid,
    load_and_process_genome_sequences,
    check_and_convert_input,
    annotate_dseqrecord,
)

from streptocad.utils import (
    polymerase_dict,
    create_primer_df_from_dict,
    ProjectDirectory,
    extract_metadata_to_dataframe,
)
from streptocad.primers.primer_generation import (
    create_idt_order_dataframe,
    make_primer_records,
    primers_to_IDT,
    find_best_check_primers_from_genome,
)
from streptocad.crispr.guideRNAcas3_9_12 import extract_sgRNAs, SgRNAargs
from streptocad.cloning.cas3_plasmid_cloning import (
    generate_cas3_protospacer_primers,
    cas3_plasmid_pcrs,
    assemble_cas3_plasmids,
)
from streptocad.cloning.gibson_cloning import (
    find_up_dw_repair_templates,
    assemble_multiple_plasmids_with_repair_templates_for_deletion,
    update_primer_names,
)
from streptocad.cloning.plasmid_processing import determine_workflow_order_for_plasmids


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


def register_workflow_6_callbacks(app):
    @app.callback(
        [
            Output("primers-output-table_6", "data"),
            Output("primers-output-table_6", "columns"),
            Output("pcr-table_6", "data"),
            Output("pcr-table_6", "columns"),
            Output("download-data-and-protocols-link_6", "href"),
            Output("filtered-df-table_6", "data"),
            Output("filtered-df-table_6", "columns"),
            Output("error-dialog_6", "message"),
            Output("error-dialog_6", "displayed"),
            Output("plasmid-metadata-table_6", "data"),
            Output("plasmid-metadata-table_6", "columns"),
        ],
        [Input("submit-settings-button_6", "n_clicks")],
        [
            State({"type": "upload-component", "index": "genome-file-6"}, "contents"),
            State({"type": "upload-component", "index": "single-vector-6"}, "contents"),
            State({"type": "upload-component", "index": "genome-file-6"}, "filename"),
            State({"type": "upload-component", "index": "single-vector-6"}, "filename"),
            State("genes-to-KO_6", "value"),
            State("forward-overhang-input_6", "value"),
            State("reverse-overhang-input_6", "value"),
            State("gc-upper_6", "value"),
            State("gc-lower_6", "value"),
            State("off-target-seed_6", "value"),
            State("off-target-upper_6", "value"),
            State("cas-type_6", "value"),
            State("number-of-sgRNAs-per-group_6", "value"),
            State("show-inframe-deletions-settings-checkbox_6", "value"),
            State("chosen-polymerase_6", "value"),
            State("melting-temperature_6", "value"),
            State("primer-concentration_6", "value"),
            State("flanking-region-number_6", "value"),
            State("restriction_enzyme_for_repair_templates_integration_6", "value"),
            State("repair_templates_length_6", "value"),
            State("overlap_for_gibson_length_6", "value"),
            State("backbone-forward-overhang-input_6", "value"),
            State("backbone-reverse-overhang-input_6", "value"),
            State("checking-primer-length_6", "value"),
            State("gibson-primer-length_6", "value"),
        ],
    )
    def run_workflow(
        n_clicks,
        genome_content,
        vector_content,
        genome_filename,
        vector_filename,
        genes_to_KO,
        forward_protospacer_overhang,
        reverse_protospacer_overhang,
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
        restriction_enzymes_2,
        repair_templates_length,
        overlap_for_gibson_length,
        backbone_fwd_overhang,
        backbone_rev_overhang,
        checking_primer_length,
        gibson_primer_length,
    ):
        if n_clicks is None:
            raise PreventUpdate

        try:
            logging.info("Workflow 6 started.")
            print(f"########  TEST OF INFRAME DELETION: {in_frame_deletion}########")

            # Initialize variables to avoid UnboundLocalError
            logging.info("Initializing variables.")
            primer_df = pd.DataFrame()
            unique_df = pd.DataFrame()

            # Create a temporary directory
            logging.info("Creating a temporary directory.")
            with tempfile.TemporaryDirectory() as tempdir:
                genome_path = os.path.join(tempdir, genome_filename)
                vector_path = os.path.join(tempdir, vector_filename)

                # Save uploaded files to the temporary directory
                logging.info(f"Saving genome file to {genome_path}.")
                save_file(genome_path, genome_content)
                logging.info(f"Saving vector file to {vector_path}.")
                save_file(vector_path, vector_content)

                # Read the GenBank files from the saved paths
                logging.info("Reading genome and vector files.")
                genome = load_and_process_genome_sequences(genome_path)[0]
                clean_plasmid = load_and_process_plasmid(vector_path)

                # Process genes to KO
                logging.info(f"Processing genes to knock out: {genes_to_KO}.")
                genes_to_KO_list = [gene.strip() for gene in genes_to_KO.split(",")]
                logging.info(f"Parsed genes to KO: {genes_to_KO_list}.")

                # Check and convert input genes for annotation
                logging.info("Checking and converting input genes for annotation.")
                target_dict, genes_to_KO_list, annotation_input = (
                    check_and_convert_input(genes_to_KO_list)
                )
                if annotation_input:
                    logging.info("Annotating genome with target genes.")
                    genome = annotate_dseqrecord(genome, target_dict)
                logging.info("Annotation completed.")

                # Extract sgRNAs
                logging.info("Extracting sgRNAs.")
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

                # Filter the DataFrame to retain only up to 6 sgRNA sequences per locus_tag
                logging.info("Filtering sgRNAs.")
                filtered_df = sgrna_df.groupby("locus_tag").head(
                    number_of_sgRNAs_per_group
                )

                logging.info("Generating CAS3 primers.")
                filtered_df_w_primers = generate_cas3_protospacer_primers(
                    filtered_df,
                    fwd_overhang=forward_protospacer_overhang,
                    rev_overhang=reverse_protospacer_overhang,
                )
                logging.info("CAS3 primers generated.")

                logging.info("Generating CAS3 plasmid PCR amplicons.")
                amplicons = cas3_plasmid_pcrs(
                    clean_plasmid,
                    filtered_df,
                    universal_fwd_seq=backbone_fwd_overhang,
                    universal_rev_seq=backbone_rev_overhang,
                )

                logging.info("Assembling CAS3 plasmids.")
                assembled_cas3_plasmids = assemble_cas3_plasmids(
                    clean_plasmid, amplicons
                )

                # Constructing a meaningful name, ID, and description for the assembled plasmid using user input
                logging.info(
                    "Constructing names, IDs, and descriptions for assembled plasmids."
                )
                targeting_info = []
                for index, row in filtered_df.iterrows():
                    formatted_str = f"pCas3_{row['locus_tag']}({row['sgrna_loc']})"
                    targeting_info.append(formatted_str)

                for i in range(len(assembled_cas3_plasmids)):
                    assembled_cas3_plasmids[i].name = f"{targeting_info[i]}"
                    assembled_cas3_plasmids[i].id = assembled_cas3_plasmids[i].name
                    assembled_cas3_plasmids[
                        i
                    ].description = f"Assembled plasmid targeting {', '.join(genes_to_KO_list)} for single gene knockout, assembled using StreptoCAD."

                # ### IDT primers
                logging.info("Making IDT primer records.")
                primer_records = make_primer_records(filtered_df_w_primers)
                idt_primers_cas3 = primers_to_IDT(primer_records)

                # In-frame deletion logic
                logging.info("Processing in-frame deletion if enabled.")
                if in_frame_deletion == [1]:
                    logging.info("Generating repair templates.")
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
                    # Convert the lists to a list of enzyme object
                    restriction_enzymes_2 = restriction_enzymes_2.split(",")
                    enzymes_for_repair_template_integration = [
                        getattr(Restriction, str(enzyme))
                        for enzyme in restriction_enzymes_2
                    ]

                    logging.info(
                        "Digesting plasmids with enzyme_for_repair_template_integration."
                    )
                    processed_records = [
                        sorted(
                            Dseqrecord(record, circular=True).cut(
                                enzymes_for_repair_template_integration
                            ),
                            key=lambda x: len(x),
                            reverse=True,
                        )[0]
                        for record in assembled_cas3_plasmids
                    ]

                    print(processed_records)
                    logging.info("Renaming processed records.")
                    for i in range(len(processed_records)):
                        processed_records[i].name = assembled_cas3_plasmids[i].name

                    logging.info(
                        "Assembling plasmids with repair templates for deletion."
                    )
                    assembly_data = (
                        assemble_multiple_plasmids_with_repair_templates_for_deletion(
                            genes_to_KO_list,
                            processed_records,
                            repair_templates_data,
                            overlap=overlap_for_gibson_length,
                        )
                    )
                    logging.info("Updating primer names.")
                    update_primer_names(assembly_data)

                    logging.info("Creating primer DataFrame from assembly data.")
                    primer_df = create_primer_df_from_dict(assembly_data)

                    logging.info("Finding unique primers.")
                    unique_df = primer_df.drop_duplicates(keep="first")
                    idt_df = create_idt_order_dataframe(
                        unique_df, concentration="25nm", purification="STD"
                    )

                    logging.info("Generating assembled contigs.")
                    assembled_contigs = []
                    for data in assembly_data:
                        contig_record = data["contig"]
                        contig_record.id = f"{data['name']}_w_rep"
                        contig_record.name = f"{data['name']}_w_rep"
                        assembled_contigs.append(contig_record)

                if in_frame_deletion == [1]:
                    logging.info(
                        "Extracting metadata for plasmid assembly with in-frame deletion."
                    )
                    integration_names = filtered_df.apply(
                        lambda row: f"sgRNA_{row['locus_tag']}({row['sgrna_loc']})",
                        axis=1,
                    ).tolist()
                    plasmid_metadata_df = extract_metadata_to_dataframe(
                        assembled_contigs, clean_plasmid, integration_names
                    )

                    workflow_df = determine_workflow_order_for_plasmids(
                        assembled_cas3_plasmids,
                        assembled_contigs,
                        restriction_enzymes_2,
                        ["NcoI", "BstBI"],
                    )

                    plasmid_metadata_df = pd.merge(
                        plasmid_metadata_df, workflow_df, on="plasmid_name", how="inner"
                    )

                else:
                    logging.info(
                        "Extracting metadata for plasmid assembly without in-frame deletion."
                    )
                    integration_names = filtered_df.apply(
                        lambda row: f"sgRNA_{row['locus_tag']}({row['sgrna_loc']})",
                        axis=1,
                    ).tolist()
                    plasmid_metadata_df = extract_metadata_to_dataframe(
                        assembled_cas3_plasmids, clean_plasmid, integration_names
                    )

                logging.info("Finding best checking primers from genome.")

                # Getting checking primers
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

                logging.info("Generating IDT order DataFrame.")
                checking_primers_df_idt = create_idt_order_dataframe(
                    checking_primers_df
                )

                if in_frame_deletion == [1]:
                    logging.info(
                        "Concatenating full IDT DataFrame with in-frame deletion primers."
                    )
                    full_idt = pd.concat(
                        [idt_primers_cas3, idt_df, checking_primers_df_idt],
                        ignore_index=True,
                    )
                    pcr_table = pd.concat([unique_df, checking_primers_df_copy])

                else:
                    logging.info(
                        "Concatenating full IDT DataFrame without in-frame deletion primers."
                    )
                    full_idt = pd.concat(
                        [idt_primers_cas3, checking_primers_df_idt], ignore_index=True
                    )
                    pcr_table = checking_primers_df_copy.copy()

                logging.info("Preparing input files for project directory.")
                input_files = [
                    {"name": "input_genome.gb", "content": genome},
                    {"name": "input_plasmid.gb", "content": clean_plasmid},
                ]

                if in_frame_deletion:
                    logging.info(
                        "Preparing output files for project directory with in-frame deletion."
                    )
                    output_files = [
                        {
                            "name": "Cas3_w_sgRNAs.gb",
                            "content": assembled_contigs,
                        },  # LIST OF Dseqrecords
                        {"name": "01_primer_df.csv", "content": primer_df},
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
                    logging.info(
                        "Preparing output files for project directory without in-frame deletion."
                    )
                    output_files = [
                        {
                            "name": "Cas3_sgRNAs.gb",
                            "content": assembled_cas3_plasmids,
                        },  # LIST OF Dseqrecords
                        {"name": "01_full_idt.csv", "content": full_idt},
                        {"name": "02_sgrna_df.csv", "content": sgrna_df},
                        {"name": "03_filtered_sgrna_df.csv", "content": filtered_df},
                        {
                            "name": "04_plasmid_metadata_df.csv",
                            "content": plasmid_metadata_df,
                        },
                    ]

                logging.info("Preparing input values for project directory.")
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
                        "forward_protospacer_overhang": str(
                            forward_protospacer_overhang
                        ),
                        "reverse_protospacer_overhang": str(
                            reverse_protospacer_overhang
                        ),
                        "backbone_fwd_overhang": str(backbone_fwd_overhang),
                        "backbone_rev_overhang": str(backbone_rev_overhang),
                    },
                }

                logging.info("Preparing Markdown file paths for project directory.")
                markdown_file_paths = [
                    "protocols/conjugation_protcol.md",
                    "protocols/single_target_crispr_plasmid_protcol.md",
                    "protocols/trouble_shooting_tips.md",
                ]

                logging.info("Creating project directory structure.")
                timestamp = datetime.utcnow().isoformat()

                project_directory = ProjectDirectory(
                    project_name=f"CRISPR_cas9_inframe_deletion_workflow_{timestamp}",
                    input_files=input_files,
                    output_files=output_files,
                    input_values=input_values,
                    markdown_file_paths=markdown_file_paths,
                )

                logging.info("Creating ZIP file for project directory.")
                zip_content = project_directory.create_directory_structure(
                    create_directories=True
                )
                data_package_encoded = base64.b64encode(zip_content).decode("utf-8")
                data_package_download_link = (
                    f"data:application/zip;base64,{data_package_encoded}"
                )

                logging.info("Preparing outputs for DataTable.")
                primer_columns = [{"name": col, "id": col} for col in full_idt.columns]
                primer_data = full_idt.to_dict("records")

                pcr_columns = [{"name": col, "id": col} for col in pcr_table.columns]
                pcr_data = pcr_table.to_dict("records")

                plasmid_metadata_df_columns = [
                    {"name": col, "id": col} for col in plasmid_metadata_df.columns
                ]
                plasmid_metadata_df_data = plasmid_metadata_df.to_dict("records")

                logging.info("Generating download link for filtered_df.")
                filtered_df_columns = [
                    {"name": col, "id": col} for col in filtered_df.columns
                ]
                filtered_df_data = filtered_df.to_dict("records")

                logging.info("Workflow 6 completed successfully.")

                # Clear the log stream after successful execution
                log_stream.truncate(0)
                log_stream.seek(0)

            return (
                primer_data,  # primers-output-table_6.data
                primer_columns,  # primers-output-table_6.columns
                pcr_data,  # pcr-table_6.data
                pcr_columns,  # pcr-table_6.columns
                data_package_download_link,  # download-data-and-protocols-link_6.href
                filtered_df_data,  # filtered-df-table_6.data
                filtered_df_columns,  # filtered-df-table_6.columns
                "",  # error-dialog_6.message (empty if no error)
                False,  # error-dialog_6.displayed (set to False if no error)
                plasmid_metadata_df_data,  # plasmid-metadata-table_6.data
                plasmid_metadata_df_columns,  # plasmid-metadata-table_6.columns
            )

        except Exception as e:
            logging.error(f"An error occurred: {str(e)}")
            error_message = (
                f"An error occurred: {str(e)}\n\nLog:\n{log_stream.getvalue()}"
            )
            return (
                [],  # primers-output-table_6.data
                [],  # primers-output-table_6.columns
                [],  # pcr-table_6.data
                [],  # pcr-table_6.columns
                "",  # download-data-and-protocols-link_6.href
                [],  # filtered-df-table_6.data
                [],  # filtered-df-table_6.columns
                error_message,  # error-dialog_6.message
                True,  # error-dialog_6.displayed
                [],  # plasmid-metadata-table_6.data
                [],  # plasmid-metadata-table_6.columns
            )
