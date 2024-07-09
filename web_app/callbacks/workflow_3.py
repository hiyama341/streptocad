#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

import os
import sys
import io
import zipfile
import base64
import csv
import pandas as pd
from Bio import SeqIO
from Bio.Restriction import NcoI, NheI
from Bio.SeqRecord import SeqRecord
from pydna.dseqrecord import Dseqrecord
from teemi.design.fetch_sequences import read_genbank_files
from teemi.build.PCR import primer_tm_neb 
import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
from dash.dash_table.Format import Group
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import quote
import tempfile


# Local module imports
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

from streptocad.primers.primer_generation import primers_to_IDT
from streptocad.cloning.golden_gate_cloning import GoldenGateCloning, create_overhang_dataframe, digest_amplicons_w_BsaI
from streptocad.primers.primer_analysis import analyze_primers_and_hairpins
from streptocad.sequence_loading.sequence_loading import load_and_process_gene_sequences, load_and_process_plasmid
from streptocad.utils import polymerase_dict, dataframe_to_seqrecords, generate_project_directory_structure
from streptocad.crispr.guideRNAcas3_9_12 import extract_sgRNAs, SgRNAargs
from streptocad.crispr.crispr_best import identify_base_editing_sites, filter_sgrnas_for_base_editing, process_base_editing
from streptocad.cloning.plasmid_processing import annotate_plasmid_with_sgrnas
from streptocad.wet_lab.gel_simulation import simulate_gel_electrophoresis
from streptocad.primers.primer_generation import create_dseqrecords_from_df, checking_primers, create_idt_order_dataframe

def save_file(name, content):
    """Decode and store a file uploaded with Plotly Dash."""
    data = content.encode("utf8").split(b";base64,")[1]
    with open(name, "wb") as fp:
        fp.write(base64.decodebytes(data))

def register_workflow_3_callbacks(app):

    @app.callback(
    [
        Output('primer-table_3', 'data'),
        Output('primer-table_3', 'columns'),
        Output('pcr-table_3', 'data'),
        Output('pcr-table_3', 'columns'),
        Output('overhang-table_3', 'data'),
        Output('overhang-table_3', 'columns'),
        Output('genbank-file_3', 'href'),
        Output('download-primers-link_3', 'href'),
        Output('download-pcr-link_3', 'href'),
        Output('download-data-and-protocols-link_3', 'href'),
    ],
    [
        Input('submit-button_3', 'n_clicks')
    ],
    [
        State('upload-genome-file_3', 'contents'),
        State('upload-single-vector_3', 'contents'),
        State('upload-genome-file_3', 'filename'),
        State('upload-single-vector_3', 'filename'),
        State('genes-to-KO_3', 'value'),
        State('sgRNA-handle-input_3', 'value'),
        State('input-tm_3', 'value'),
        State('gc-upper_3', 'value'),
        State('gc-lower_3', 'value'),
        State('off-target-seed_3', 'value'),
        State('off-target-upper_3', 'value'),
        State('cas-type_3', 'value'),
        State('number-of-sgRNAs-per-group_3', 'value'),
        State('only-stop-codons-checkbox_3', 'value'),
        State('chosen-polymerase_3', 'value'),
        State('melting-temperature_3', 'value'),
        State('primer-concentration_3', 'value'),
        State('primer-number-increment_3', 'value'),
        State('flanking-region-number_3', 'value')
    ]
)
    def run_workflow(n_clicks, genome_content, vector_content, genome_filename, vector_filename, genes_to_KO, 
                     sgRNA_handle_input, input_tm, gc_upper, gc_lower, off_target_seed, off_target_upper, cas_type, 
                     number_of_sgRNAs_per_group, only_stop_codons, chosen_polymerase, melting_temperature, 
                     primer_concentration, primer_number_increment, flanking_region_number):
        if n_clicks is None:
            raise PreventUpdate

        try:
            print("Workflow 3 started")

            # Create a temporary directory
            with tempfile.TemporaryDirectory() as tempdir:
                genome_path = os.path.join(tempdir, genome_filename)
                vector_path = os.path.join(tempdir, vector_filename)

                # Save uploaded files to the temporary directory
                save_file(genome_path, genome_content)
                save_file(vector_path, vector_content)

                # Read the GenBank files from the saved paths
                print("Reading genome and vector files")
                genome = load_and_process_plasmid(genome_path)
                plasmid = load_and_process_plasmid(vector_path)

                genes_to_KO_list = [gene.strip() for gene in genes_to_KO.split(',')]
                print(f"Genes to knock out: {genes_to_KO_list}")

                # Initialize SgRNAargs with desired parameters
                args = SgRNAargs(genome_path, 
                                 genes_to_KO_list,
                                 step=['find', 'filter'],
                                 gc_upper=gc_upper,
                                 gc_lower=gc_lower,
                                 off_target_seed=off_target_seed,
                                 off_target_upper=off_target_upper,
                                 cas_type=cas_type)

                print("Extracting sgRNAs")
                sgrna_df = extract_sgRNAs(args)
                print(f"sgRNA DataFrame: {sgrna_df}")

                # Load gene sequences
                print("Loading gene sequences")
                gene_sequences = load_and_process_gene_sequences(genome_path)
                genes_to_KO_dict = {locus_tag: gene_sequences[locus_tag] for locus_tag in genes_to_KO_list if locus_tag in gene_sequences}
                print(f"Genes to KO dictionary: {genes_to_KO_dict}")

                # Identify and annotate base editing sites
                print("Identifying base editing sites")
                sgrna_df_with_editing = identify_base_editing_sites(sgrna_df)

                # Filter out only sgRNAs that result in base-editing
                print("Filtering sgRNAs for base editing")
                filtered_sgrna_df_for_base_editing = filter_sgrnas_for_base_editing(sgrna_df_with_editing)
                print(f"Filtered sgRNAs for base editing: {filtered_sgrna_df_for_base_editing}")

                # Process the DataFrame to apply C-to-T mutations
                print("Processing base editing")
                mutated_sgrna_df = process_base_editing(filtered_sgrna_df_for_base_editing, 
                                                        genes_to_KO_dict, 
                                                        only_stop_codons=bool(only_stop_codons))
                print(f"Mutated sgRNAs: {mutated_sgrna_df}")

                # Filter the DataFrame to retain only up to 5 sgRNA sequences per locus_tag
                print("Filtering sgRNAs to retain only up to 5 sequences per locus tag")
                filtered_df = mutated_sgrna_df.groupby('locus_tag').head(number_of_sgRNAs_per_group)
                print(f"Filtered DataFrame: {filtered_df}")

                # Generate sgRNA list and handle sites
                print("Generating sgRNA list")
                sgRNA_list = dataframe_to_seqrecords(filtered_df)
                sgRNA_handle_cys4_sites = [Dseqrecord(sgRNA_handle_input, name='sgRNA_handle_cys4')] * len(sgRNA_list)

                # Golden Gate Cloning
                print("Performing Golden Gate Cloning")
                golden_gate = GoldenGateCloning(sgRNA_list,
                                                sgRNA_handle_cys4_sites, 
                                                target_tm=input_tm,
                                                restriction_overhang_f="GATCGggtctcc",
                                                restriction_overhang_r="GATCAGGTCTCg",
                                                backbone_overhang_f="cATG",
                                                backbone_overhang_r="cTAG",
                                                cys4="gTTCACTGCCGTATAGGCAGCTAAGAAA",
                                                tm_function=primer_tm_neb, 
                                                primer_incrementation=primer_number_increment, 
                                                polymerase=chosen_polymerase)

                # Generate primers
                print("Generating primer DataFrame")
                primer_df = golden_gate.generate_primer_dataframe()
                print(f"Primer DataFrame: {primer_df}")

                # Primer analysis
                print("Analyzing primers")
                analysis_of_primers = analyze_primers_and_hairpins(primer_df)
                print(f"Analysis of primers: {analysis_of_primers}")

                # Generate IDT order DataFrame
                print("Generating IDT order DataFrame")
                idt_df = create_idt_order_dataframe(primer_df, concentration="25nm", purification="STD")
                print(f"IDT DataFrame: {idt_df}")

                # Simulate PCR
                print("Simulating PCRs")
                list_of_amplicons = golden_gate.simulate_pcrs()

                # Simulate gel electrophoresis
                print("Simulating gel electrophoresis")
                simulate_gel_electrophoresis(list_of_amplicons)

                # Generate overhang DataFrame
                print("Generating overhang DataFrame")
                overhangs = create_overhang_dataframe(list_of_amplicons)
                print(f"Overhang DataFrame: {overhangs}")

                # Digest amplicons
                print("Digesting amplicons")
                digest_amplicons = digest_amplicons_w_BsaI(list_of_amplicons)
                for digest in digest_amplicons: 
                    print(digest.figure())

                # Digest and assemble plasmid
                print("Digesting and assembling plasmid")
                linear_plasmid, _ = sorted(plasmid.cut(NcoI, NheI), key=lambda x: len(x), reverse=True)
                for amplicon in digest_amplicons:
                    linear_plasmid.seq += amplicon.seq
                
                rec_vec = linear_plasmid.looped()

                # Annotate plasmid
                print("Annotating plasmid")
                annotate_plasmid_with_sgrnas(rec_vec, filtered_df)

                # Generate checking primers
                print("Generating checking primers")
                checking_primers_df = checking_primers(genome_path, genes_to_KO_list, 
                                                       flanking_region=flanking_region_number,
                                                       target_tm=melting_temperature, 
                                                       primer_concentration=primer_concentration, 
                                                       polymerase=chosen_polymerase)
                checking_primers_df_idt = create_idt_order_dataframe(checking_primers_df)

                # Combine IDT DataFrames
                print("Combining IDT DataFrames")
                full_idt = pd.concat([idt_df, checking_primers_df_idt])
                print(f"Full IDT DataFrame: {full_idt}")

                # Prepare outputs for the DataTable
                print("Preparing data for DataTables")
                primers_columns = [{"name": col, "id": col} for col in full_idt.columns]
                primers_data = full_idt.to_dict('records')

                pcr_columns = [{"name": col, "id": col} for col in checking_primers_df.columns]
                pcr_data = checking_primers_df.to_dict('records')

                overhang_columns = [{"name": col, "id": col} for col in overhangs.columns]
                overhang_data = overhangs.to_dict('records')

                # IDT download links
                primer_df_string = full_idt.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
                primer_data_encoded = quote(primer_df_string)
                primer_download_link = f"data:text/csv;charset=utf-8,{primer_data_encoded}" 

                pcr_df_string = checking_primers_df.to_csv(index=False, quoting=csv.QUOTE_NONNUMERIC)
                pcr_data_encoded = quote(pcr_df_string)
                pcr_download_link = f"data:text/csv;charset=utf-8,{pcr_data_encoded}" 

                # Provide a link for downloading GenBank files
                encoded_genbank_file = base64.b64encode(str(rec_vec).encode()).decode()
                genbank_download_link = f"data:application/genbank;base64,{encoded_genbank_file}"

                # Data and protocols download link
                input_dict = {
                    '######Inputs####': '',
                    "Path_to_file": genome_filename,
                    "Name of genome": genome.name,
                    "Path to plasmid": vector_filename,
                    "Plasmid name": plasmid.name,
                    'Genes to Base-edit': genes_to_KO_list, 
                    'Flanking region for checking primers': flanking_region_number,
                    'Upper GC content for sgRNAs': gc_upper,
                    'Lower GC content for sgRNAs': gc_lower,
                    'Off target seed': off_target_seed, 
                    'Off target upper': off_target_upper, 
                    'Cas type': cas_type, 
                    'Number of sgRNAs per locus tag': number_of_sgRNAs_per_group,
                    'Stop codons': bool(only_stop_codons), 
                    "Chosen polymerase": chosen_polymerase,
                    "Melting temperature": melting_temperature,
                    "Primer_concentration": primer_concentration,
                    "Primer number increment": primer_number_increment,
                    '####### OUTPUTS #########' : '',
                    'Duplication of overhangs': overhangs,
                    'PCR dataframe': primer_df,
                    'Primer analysis dataframe': analysis_of_primers,
                    "Primer df": full_idt,
                    "All mutations dataframe": mutated_sgrna_df,
                    'Filtered dataframe' : filtered_df,
                    'Plasmid assembly': rec_vec.figure(),
                }

                description = (
                    "This document captures all the values that were put into the application "
                    "to ensure that the experiment can be replicated accurately.\n"
                    "Below are the details of the inputs used:\n"
                )

                data_and_protocols_package = format_and_print_values(input_dict, description, spacing=10)
                data_package_encoded = base64.b64encode(data_and_protocols_package.encode('utf-8')).decode('utf-8')
                data_package_download_link = f"data:text/plain;base64,{data_package_encoded}"  

                print("Workflow 3 completed successfully")

                return (primers_data, primers_columns, pcr_data, pcr_columns, overhang_data, overhang_columns,
                        genbank_download_link, primer_download_link, pcr_download_link, data_package_download_link)

        except Exception as e:
            print("An error occurred:", str(e))
            raise PreventUpdate
