# tooltips.py
import dash_bootstrap_components as dbc


def create_tooltip(content, target_id, placement="top", style=None):
    return dbc.Tooltip(content, target=target_id, placement=placement, style=style)


# Generalized Tooltips for Common Use Cases
tooltips = [
    # Polymerase
    create_tooltip(
        "The selected polymerase is used to calculate primer melting temperature with NEB's api",
        target_id="polymerase-tooltip",
    ),
    create_tooltip(
        "The selected polymerase is used to calculate primer melting temperature with NEB's api",
        target_id="polymerase-tooltip-1",
    ),
    # Melting Temperature
    create_tooltip(
        "Set the desired melting temperature for PCR reactions. The melting temperature is used to calculate the Tm/Ta of primers (with NEBs api).",
        target_id="melting-temp-tooltip",
    ),
    # Primer Concentration
    create_tooltip(
        "Specify the concentration of primers (in μM) for PCR reactions. A normal range is between 0.1-0.5 μM.",
        target_id="primer-concentration-tooltip",
    ),
    # Flanking Region
    create_tooltip(
        "Define the size (in bp) of flanking regions for checking primers. For example, if you are checking primers for a 1000 bp gene, you and you set flanking region to 100 bp, then the primers will be made 100 bp upstream and downstream of the gene. It the algorithm actually checks for secondary structure, primer dimers etc. so the final flanking region migt be higer than specified.",
        target_id="flanking-region-tooltip",
    ),
    # GC Content Boundaries
    create_tooltip(
        "Set the upper bound for GC content in sgRNA sequences. This means that sequences with GC content above this will be excluded. ",
        target_id="gc-upper-tooltip",
    ),
    create_tooltip(
        "Set the lower bound for GC content in sgRNA sequences. This means that sequences with GC content below this will be excluded.",
        target_id="gc-lower-tooltip",
    ),
    # Off-Target Seed Length and Upper Bound
    create_tooltip(
        "Specify the length of the seed sequence for off-target filtering. The seed sequence is used to filter off-target sgRNAs. It is defined as the seed_length upstream and downstream of the target site upstream of the PAM site (Cas9). For other Cas-systems like Cas3 it is the upstream. If you want have hard criteria on you set this lower and if more laxed criterias are wanted, set it high.",
        target_id="off-target-seed-tooltip",
    ),
    create_tooltip(
        "Specify the maximum allowed off-target score for sgRNA filtering. For example if you dont want any off-targets set this to 0. If you dont care so much u can use the default setting.",
        target_id="off-target-upper-tooltip",
    ),
    # Cas Type
    create_tooltip(
        "Cas9 has a PAM sequence that is NGG and Cas3 has TTC. Cas12a has a PAM sequence that is TTT[ACG].",
        target_id="cas-type-tooltip",
    ),
    # Number of sgRNAs per Region
    create_tooltip(
        "Define the number of sgRNAs for each targeted region or locus. For example, if you have added 10 genes to to be deleted, and you want to test 5 sgRNAs per gene, set this to 5. This will produce 50 sgRNAs.",
        target_id="number-sgrnas-tooltip",
    ),
    # Only Stop Codons
    create_tooltip(
        "Filter to include only stop codons when base editing.",
        target_id="only-stop-codons-tooltip",
    ),
    # Editing Context
    create_tooltip(
        "Enable sequence context editing to increase success rate of base editing. It refines sgRNA selection by excluding sgRNAs with an upstream “G” adjacent to the editable “C,” ensuring higher editing efficiency.",
        target_id="editing-context-tooltip",
    ),
    # Repair Template and Overlap for Gibson Cloning
    create_tooltip(
        "Specify the length of repair templates. If you select 1000 bp, the algorithm will produce 1000 bp repair templates with overlaps to your plasmid backbone.",
        target_id="repair-templates-length-tooltip",
    ),
    create_tooltip(
        "Define the overlap length between templates for Gibson assembly. Typically, 30-40 bp produces good results experimentally.",
        target_id="overlap-gibson-tooltip",
    ),
    # Overhangs (Forward and Reverse)
    create_tooltip(
        "Enter the 5' overhang sequence for ssDNA bridging. This DNA will be added upstream to the sgRNA that is found and made into an oligo.",
        target_id="forward-overhang-tooltip",
    ),
    create_tooltip(
        "Enter the 3' overhang sequence for ssDNA bridging. This DNA will be added downstream to the sgRNA that is found and made into an oligo.",
        target_id="reverse-overhang-tooltip",
    ),
    # sgRNA Handle
    create_tooltip(
        "Specify the sgRNA handle sequence for the CRISPR system. Primers will be made annealing to this sequence.",
        target_id="sgRNA-handle-tooltip",
    ),
    # Advanced Settings
    create_tooltip(
        "Show advanced settings for more functionalities.",
        target_id="show-advanced-settings-tooltip",
    ),
    # Golden Gate
    create_tooltip(
        "Specify the forward overhang for restriction-based cloning. This overhang will be integrated into the fwd-primers that BsaI will recognize cut so it is compatible with Golden Gate cloning.",
        target_id="restriction-overhang-f-tooltip",
    ),
    create_tooltip(
        "Specify the reverse overhang for restriction-based cloning. This overhang will be integrated into the rev-primers that BsaI will recognize cut so it is compatible with Golden Gate cloning.",
        target_id="restriction-overhang-r-tooltip",
    ),
    create_tooltip(
        "Define the forward overhang sequence for the plasmid backbone when it has been digested. Depending on the plasmid this can change (the default works with the plasmid you can download above)",
        target_id="backbone-overhang-f-tooltip",
    ),
    create_tooltip(
        "Define the reverse overhang sequence for the plasmid backbone when it has been digested. Depending on the plasmid this can change (the default works with the plasmid you can download above)",
        target_id="backbone-overhang-r-tooltip",
    ),
    create_tooltip(
        "Enter the Csy4 sequence for multiplexing in CRISPR-BEST. This sequence is important for the processing of all sgRNAs and will be integrated in between the sgRNA array.",
        target_id="cys4-sequence-tooltip",
    ),
    # Melting Temperature
    create_tooltip(
        "Set the desired melting temperature for your checking primers PCR reactions. The algorithms will match the primers within a few degrees. ",
        target_id="checking-melting-temp-tooltip",
    ),
    # Protospacer 5' Anneal
    create_tooltip(
        "Enter the forward 5' annealing sequence for the protospacer, used in Gibson assembly. This is the primer that has to overlap with the plasmid and where half of the sgRNA is integrated.",
        target_id="forward-overhang-tooltip-6",
    ),
    # Protospacer 3' Anneal
    create_tooltip(
        "Enter the reverse 3' annealing sequence for the protospacer, used in Gibson assembly. This is the primer that has to overlap with the plasmid and where half of the sgRNA is integrated.",
        target_id="reverse-overhang-tooltip-6",
    ),
    # Backbone 5' Anneal
    create_tooltip(
        "Standardized primer: Specify the forward 5' annealing sequence for the backbone in the Gibson assembly process.",
        target_id="backbone-forward-overhang-tooltip-6",
    ),
    # Backbone 3' Anneal
    create_tooltip(
        "Standardized primer: Specify the reverse 3' annealing sequence for the backbone in the Gibson assembly process.",
        target_id="backbone-reverse-overhang-tooltip-6",
    ),
    # Restriction Enzymes
    create_tooltip(
        """Choose the restriction enzyme(s) for the plasmid digestion.
            This is just for your convenience. 
            By default, we use NcoI for the digestion of  pCRISPR-cBEST.
            If you would like to use a different enzyme, please enter the name of the enzyme(s) in the following format:
            enzyme1,enzyme2,enzyme3,...
            Tip: Remember to update overlapping sequences for the assembly if you change the enzymes.
            """,
        target_id="restriction-enzymes-tooltip-2",
    ),
    create_tooltip(
        """Choose the restriction enzyme(s) for the plasmid digestion.
            This is just for your convenience. 
            By default, we use NcoI and NheI for the digestion of  pCRISPR-MCBE_Csy4_kasopGFP.
            If you would like to use a different enzyme, please enter the name of the enzyme(s) in the following format:
            enzyme1,enzyme2,enzyme3,...
            Tip: Remember to update overlapping sequences for the assembly if you change the enzymes.
            """,
        target_id="restriction-enzymes-tooltip-3",
    ),
    create_tooltip(
        """Choose the restriction enzyme(s) for the plasmid digestion.
            This is just for your convenience. 
            By default, we use NcoI  the digestion of  pCRISPR-dCas9.
            If you would like to use a different enzyme, please enter the name of the enzyme(s) in the following format:
            enzyme1,enzyme2,enzyme3,...
            Tip: Remember to update overlapping sequences for the assembly if you change the enzymes.
            """,
        target_id="restriction-enzymes-tooltip-4",
    ),
    create_tooltip(
        """Choose the restriction enzyme(s) for the plasmid digestion.
            This is just for your convenience. 
            By default, we use NcoI  the digestion of  CRISPR-Cas9.
            If you would like to use a different enzyme, please enter the name of the enzyme(s) in the following format:
            enzyme1,enzyme2,enzyme3,...
            Tip: Remember to update overlapping sequences for the assembly if you change the enzymes.
            """,
        target_id="restriction-enzymes-tooltip-5",
    ),
    create_tooltip(
        """Choose the restriction enzyme(s) for the plasmid digestion.
            This is just for your convenience. 
            
            By default, we use EcoRI for the digestion of  CRISPR-Cas9 for the integration of repair templates. But you can use any enzyme from the MCS.
            If you would like to use a different enzyme, please enter the name of the enzyme(s) in the following format:
            enzyme1,enzyme2,enzyme3,...
            Tip: Remember to update overlapping sequences for the assembly if you change the enzymes.
            """,
        target_id="restriction-enzymes-tooltip-5-2",
    ),
    create_tooltip(
        """Choose the restriction enzyme(s) for the plasmid digestion.
            This is just for your convenience. 
            By default, we use EcoRI for the digestion of  CRISPR-Cas3. But you can use any enzyme from the MCS e.g. HindIII etc
            If you would like to use a different enzyme, please enter the name of the enzyme(s) in the following format:
            enzyme1,enzyme2,enzyme3,...
            Tip: Remember to update overlapping sequences for the assembly if you change the enzymes.
            """,
        target_id="restriction-enzymes-tooltip-6",
    ),
    create_tooltip(
        """ This setting overrides your desired melting temperature for improved primer specificity.
            Tip: Use GC enhancers if the melting temperature remains too high.
            Additional Tips:
            • Experiment with different minimum primer lengths and melting temperatures to optimize primer design for your specific sequences.
            • Adjust these parameters to find the best balance between specificity and efficiency for your amplification needs. 
            """,
        target_id="primer-anneal-length-tooltip",
    ),
]
