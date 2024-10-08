flowchart TD
%% Workflow 6 - CRISPR-Cas3 Plasmid Generation with Conditional In-frame Deletion
subgraph W6 [CRISPR-Cas3 Plasmid Generation with Conditional In-frame Deletion]
%% Inputs
A1["Genome File"] --> B1
A2["Plasmid File"] --> C1
A3["Genes to Knock Out"] --> D1
A4["GC Upper"] --> F1
A5["GC Lower"] --> F1
A6["Off-Target Seed"] --> F1
A7["Off-Target Upper"] --> F1
A8["Cas Type"] --> F1
A9["Forward Protospacer Overhang"] --> G1
A10["Reverse Protospacer Overhang"] --> G1
A11["Number of sgRNAs per Group"] --> G1

        A13["Repair Templates Length"] --> M1
        A14["Overlap for Gibson Length"] --> M1
        A15["Backbone Forward Overhang"] --> G1
        A16["Backbone Reverse Overhang"] --> G1
        A17["Enzyme for Repair Template Integration"] --> N1

        %% Process Steps
        B1["Load and Process Genome (load_and_process_genome_sequences)"]
        C1["Load and Process Plasmid (load_and_process_plasmid)"]
        D1["Check and Convert Input (check_and_convert_input)"] --> E1
        E1["Annotate Genome (annotate_dseqrecord)"] --> F1
        F1["Extract sgRNAs (extract_sgRNAs)"] --> G1
        G1["Filter sgRNAs (filter_sgRNAs_for_knockout)"] --> H1
        H1["Generate CAS3 Primers (generate_cas3_protospacer_primers)"] --> I1
        I1["Perform CAS3 Plasmid PCR (cas3_plasmid_pcrs)"] --> J1
        J1["Assemble CAS3 Plasmids (assemble_cas3_plasmids)"]

        %% Conditional In-Frame Deletion Path
        J1 -->|If In-Frame Deletion Enabled| M1
        M1["Find Repair Templates (find_up_dw_repair_templates)"] --> N1
        J1 --> N1["Digest Plasmids for Repair Template Integration (cut with specified enzyme)"] --> O1
        O1["Assemble Plasmids with Repair Templates (assemble_multiple_plasmids_with_repair_templates_for_deletion)"] --> P1

        %% Else Path for No In-Frame Deletion
        K1 -->|else| P1["Annotate Plasmid with sgRNAs (annotate_plasmid_with_sgrnas)"]

        %% Plasmid Annotation and Metadata Extraction
        P1 --> Q1["Extract Plasmid Metadata (extract_metadata_to_dataframe)"] --> R1

        %% Additional Steps for Project Directory
        B1 --> S1["Generate Primers for IDT (primers_to_IDT)"]
        S1 --> T1["Create IDT Order Format (create_idt_order_dataframe)"]
        T1 --> R1["Project Directory (ProjectDirectory)"] --> U1["FAIR Datafolder"]
    end

    %% Connections
    E1 -->|Genome File| F1
    C1 -->|Plasmid File| J1
