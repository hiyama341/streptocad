flowchart TD
%% Workflow 5 - CRISPR-Cas9 Plasmid Generation with Conditional In-frame Deletion
subgraph W5 [CRISPR-Cas9 Plasmid Generation with Conditional In-frame Deletion]
%% Inputs
A1["Genome File"] --> B1
A2["Plasmid File"] --> C1
A3["Genes to Knock Out"] --> D1
A4["GC Upper"] --> F1
A5["GC Lower"] --> F1
A6["Off-Target Seed"] --> F1
A7["Off-Target Upper"] --> F1
A8["Cas Type"] --> F1
A9["Upstream Homology"] --> H1
A10["Downstream Homology"] --> H1
A11["Number of sgRNAs per Group"] --> G1
A12["In-Frame Deletion Option"] --> I1
A13["Repair Templates Length"] --> K1
A14["Overlap for Gibson Length"] --> K1

        %% Process Steps
        B1["Load and Process Genome (load_and_process_genome_sequences)"]
        C1["Load and Process Plasmid (load_and_process_plasmid)"]
        D1["Check and Convert Input (check_and_convert_input)"] --> E1
        E1["Annotate Genome (annotate_dseqrecord)"] --> F1
        F1["Extract sgRNAs (extract_sgRNAs)"] --> G1
        G1["Filter sgRNAs (filter_sgRNAs_for_knockout)"] --> H1
        H1["Make ssDNA Oligos (make_ssDNA_oligos)"] --> I1
        I1["Cut Plasmid with NcoI (cut)"] --> J1
        J1["Assemble Plasmid by ssDNA Bridging (assemble_plasmids_by_ssDNA_bridging)"]

        %% Conditional In-Frame Deletion
        I1 -->|If In-Frame Deletion Enabled| K1
        K1["Find Repair Templates (find_up_dw_repair_templates)"] --> L1
        L1["Assemble Plasmids with Repair Templates (assemble_multiple_plasmids_with_repair_templates_for_deletion)"] --> M1

        %% Else Path for No In-Frame Deletion
        I1 -->|else| M1["Annotate Plasmid with sgRNAs (annotate_plasmid_with_sgrnas)"]

        %% Plasmid Annotation and Metadata Extraction
        M1 --> N1["Extract Plasmid Metadata (extract_metadata_to_dataframe)"] --> O1

        %% Additional Steps for Project Directory
        B1 --> P1["Generate Primers for IDT (primers_to_IDT)"]
        P1 --> Q1["Create IDT Order Format (create_idt_order_dataframe)"]
        Q1 --> O1["Project Directory (ProjectDirectory)"] --> R1["FAIR Datafolder"]
    end

    %% Connections
    E1 -->|Genome File| F1
    C1 -->|Plasmid File| J1
