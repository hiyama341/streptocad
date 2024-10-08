%%{init:{"theme":"neutral"}}%%
flowchart TD
%% Workflow 2 - CRISPR-BEST Plasmid Generation
subgraph W2 ["CRISPR-BEST Plasmid Generation"]
%% Inputs
A1["Input:<br>Genome File"] --> B1
A2["Input:<br>Plasmid File"] --> C1
A3["Input:<br>Genes to Knock Out"] --> D1
A4["Input:<br>GC Upper"] --> F1
A5["Input:<br>GC Lower"] --> F1
A6["Input:<br>Off-Target Seed"] --> F1
A7["Input:<br>Off-Target Upper"] --> F1
A8["Input:<br>Cas Type"] --> F1
A12["Input:<br>Number of sgRNAs per Group"] --> H1
A13["Input:<br>Upstream Homology"] --> J1
A14["Input:<br>Downstream Homology"] --> J1
A9["Input:<br>Melting Temperature"] --> O1
A10["Input:<br>Primer Concentration"] --> O1
A15["Input:<br>Flanking Region Number"] --> O1
A16["Input:<br>Only Stop Codons"] --> I1
A17["Input:<br>Editing Context"] --> I1

        %% Process Steps
        B1["Step 1:<br>Load and Process Genome<br>(load_and_process_genome_sequences)"]
        C1["Step 2:<br>Load and Process Plasmid<br>(load_and_process_plasmid)"]
        D1["Step 3:<br>Check and Convert Input<br>(check_and_convert_input)"] --> E1
        E1["Step 4:<br>Annotate Genome<br>(annotate_dseqrecord)"] --> F1
        F1["Step 5:<br>Extract sgRNAs<br>(extract_sgRNAs)"] --> G1
        G1["Step 6:<br>Identify Base Editing Sites<br>(identify_base_editing_sites)"] --> H1
        H1["Step 7:<br>Filter sgRNAs for Base Editing<br>(filter_sgrnas_for_base_editing)"] --> I1
        I1["Step 8:<br>Process Base Editing<br>(process_base_editing)"] --> J1
        J1["Step 9:<br>Make ssDNA Oligos<br>(make_ssDNA_oligos)"] --> K1
        K1["Step 10:<br>Cut Plasmid with NcoI<br>(cut)"] --> L1
        L1["Step 11:<br>Assemble Plasmid by ssDNA Bridging<br>(assemble_plasmids_by_ssDNA_bridging)"] --> M1
        M1["Step 12:<br>Annotate Plasmid with sgRNAs<br>(annotate_plasmid_with_sgrnas)"] --> N1
        N1["Step 13:<br>Extract Plasmid Metadata<br>(extract_metadata_to_dataframe)"]

        %% Additional Steps for Project Directory
        B1 --> O1["Generate Check Primers<br>(find_best_check_primers_from_genome)"]
        O1 --> P1["Create IDT Order Format<br>(create_idt_order_dataframe)"]
        P1 --> Q1["Project Directory<br>(ProjectDirectory)"] --> R1["Final Output:<br>FAIR Datafolder"]
    end

    %% Styling
    classDef input fill:#f9f9f9,stroke:#333,stroke-width:2px,color:#333,stroke-dasharray: 5 5;
    classDef process fill:#dae8fc,stroke:#6c8ebf,stroke-width:2px,color:#333,border-radius:8px;
    classDef output fill:#d5e8d4,stroke:#82b366,stroke-width:2px,color:#333,border-radius:8px;
    classDef final fill:#fff2cc,stroke:#d6b656,stroke-width:3px,color:#333,border-radius:8px;
    classDef cluster fill:#ffffff,stroke:#999999,stroke-width:3px,color:#999999;

    class A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A12,A13,A14,A15,A16,A17 input;
    class B1,C1,D1,E1,F1,G1,H1,I1,J1,K1,L1,M1,N1,O1,P1,Q1 process;
    class R1 final;
    class W2 cluster;
