flowchart TD
%% Workflow 4 - CRISPRi Plasmid Generation
subgraph W4 [CRISPRi Plasmid Generation]
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
A12["Extension to Promoter Region"] --> F1

        %% Process Steps
        B1["Load and Process Genome (load_and_process_genome_sequences)"]
        C1["Load and Process Plasmid (load_and_process_plasmid)"]
        D1["Check and Convert Input (check_and_convert_input)"] --> E1
        E1["Annotate Genome (annotate_dseqrecord)"] --> F1
        F1["Extract sgRNAs for CRISPRi (extract_sgRNAs_for_crispri)"] --> G1
        G1["Filter sgRNAs for Promoter Regions (filter_sgrnas_for_promoters)"] --> H1
        H1["Make ssDNA Oligos (make_ssDNA_oligos)"] --> I1
        I1["Cut Plasmid with NcoI (cut)"] --> J1
        J1["Assemble Plasmid by ssDNA Bridging (assemble_plasmids_by_ssDNA_bridging)"] --> K1
        K1["Annotate Plasmid with sgRNAs (annotate_plasmid_with_sgrnas)"] --> L1
        L1["Extract Plasmid Metadata (extract_metadata_to_dataframe)"] --> O1

        %% Additional Steps for Project Directory
        B1 --> M1["Generate Primers for IDT (primers_to_IDT)"]
        M1 --> N1["Create IDT Order Format (create_idt_order_dataframe)"]
        N1 --> O1["Project Directory (ProjectDirectory)"] --> P1["FAIR Datafolder"]
    end

    %% Connections
    E1 -->|Genome File| F1
    C1 -->|Plasmid File| J1
