%%{init:{"theme":"neutral"}}%%
flowchart TD
%% Workflow 1 - Overexpressing Plasmid Generation
subgraph W1 pOXe-PkasO
%% Inputs
A1["Input:<br>Sequence File"]
A2["Input:<br>Plasmid File"]
A3["Input:<br>Melting Temperature"]
A4["Input:<br>Chosen Polymerase"]
A5["Input:<br>Primer Concentration"]
A6["Input:<br>Upstream Homology"]
A7["Input:<br>Downstream Homology"]
A8["Input:<br>Primer Number Increment"]

        %% Process Steps
        A1 --> B1["Step 1:<br>Generate Primer Dataframe<br>(generate_primer_dataframe)"]
        A3 --> B1
        A4 --> B1
        A5 --> B1
        A6 --> B1
        A7 --> B1
        A8 --> B1
        B1 --> C1["Step 2:<br>Perform PCR on Sequences<br>(perform_pcr_on_sequences)"]
        C1 --> D1["Step 3:<br>Analyze Primers and Hairpins<br>(analyze_primers_and_hairpins)"]
        A2 --> E1["Step 4:<br>Plasmid Assembly<br>(assemble_and_process_plasmids)"]
        D1 --> E1
        E1 --> F1["Step 5:<br>IDT Order Format<br>(create_idt_order_dataframe)"]
        F1 --> G1["Step 6:<br>Project Directory<br>(ProjectDirectory)"]
        G1 --> H1["Final Output:<br>FAIR Datafolder"]
    end

    %% Styling
    classDef input fill:#f9f9f9,stroke:#333,stroke-width:2px,color:#333,stroke-dasharray: 5 5;
    classDef process fill:#dae8fc,stroke:#6c8ebf,stroke-width:2px,color:#333,border-radius:8px;
    classDef output fill:#d5e8d4,stroke:#82b366,stroke-width:2px,color:#333,border-radius:8px;
    classDef final fill:#fff2cc,stroke:#d6b656,stroke-width:3px,color:#333,border-radius:8px;
    classDef cluster fill:#ffffff,stroke:#999999,stroke-width:3px,color:#999999;

    class A1,A2,A3,A4,A5,A6,A7,A8 input;
    class B1,C1,D1,E1,F1,G1 process;
    class H1 final;
    class W1 cluster;
