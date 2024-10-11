# tooltips.py
import dash_bootstrap_components as dbc

def create_tooltip(content, target_id, placement="top", style=None):
    return dbc.Tooltip(content, target=target_id, placement=placement, style=style)

# Generalized Tooltips for Common Use Cases
tooltips = [
    # Polymerase
    create_tooltip(
        "The selected polymerase is used to calculate primer melting temperature with NEB's api",
        target_id="polymerase-tooltip"
    ),
    
    # Melting Temperature
    create_tooltip(
        "Set the desired melting temperature for PCR reactions. The melting temperature is used to calculate the Tm of primers.",
        target_id="melting-temp-tooltip"
    ),

    # Primer Concentration
    create_tooltip(
        "Specify the concentration of primers (in μM) for PCR reactions.",
        target_id="primer-concentration-tooltip"
    ),


    # Flanking Region
    create_tooltip(
        "Define the size (in bp) of flanking regions for checking primers. For example, if you are checking primers for a 1000 bp gene, you and you set flanking region to 100 bp, then the primers will be checked for 100 bp upstream and downstream of the gene.",
        target_id="flanking-region-tooltip"
    ),

    # GC Content Boundaries
    create_tooltip(
        "Set the upper bound for GC content in sgRNA sequences.",
        target_id="gc-upper-tooltip"
    ),
    create_tooltip(
        "Set the lower bound for GC content in sgRNA sequences.",
        target_id="gc-lower-tooltip"
    ),

    # Off-Target Seed Length and Upper Bound
    create_tooltip(
        "Specify the length of the seed sequence for off-target filtering. The seed sequence is used to filter off-target sgRNAs. It is defined as the seed_length upstream and downstream of the target site upstream of the PAM site (Cas9).",
        target_id="off-target-seed-tooltip"
    ),
    create_tooltip(
        "Specify the maximum allowed off-target score for sgRNA filtering. For example if you dont want any off-targets set this to 0.",
        target_id="off-target-upper-tooltip"
    ),

    # Cas Type
    create_tooltip(
        "Cas9 has a PAM sequence that is NGG and Cas3 has TTC. Cas12a has a PAM sequence that is TTT[ACG].",
        target_id="cas-type-tooltip"
    ),

    # Number of sgRNAs per Region
    create_tooltip(
        "Define the number of sgRNAs for each targeted region or locus. For example, if you have added 10 genes to to be deleted, and you want to test 5 sgRNAs per gene, set this to 5. This will produce 50 sgRNAs.",
        target_id="number-sgrnas-tooltip"
    ),

    # Only Stop Codons
    create_tooltip(
        "Filter to include only stop codons when base editing.",
        target_id="only-stop-codons-tooltip"
    ),

    # Editing Context
    create_tooltip(
        "Enable sequence context editing to increase success rate of base editing. It refines sgRNA selection by excluding sgRNAs with an upstream “G” adjacent to the editable “C,” ensuring higher editing efficiency.",
        target_id="editing-context-tooltip"
    ),

    # Repair Template and Overlap for Gibson Cloning
    create_tooltip(
        "Specify the length of repair templates. If you select 1000 bp, the algorithm will produce 1000 bp repair templates with overlaps to your plasmid backbone.",
        target_id="repair-templates-length-tooltip"
    ),
    create_tooltip(
        "Define the overlap length between templates for Gibson assembly. Typically, 30-40 bp produces good results experimentally.",
        target_id="overlap-gibson-tooltip"
    ),

    # Overhangs (Forward and Reverse)
    create_tooltip(
        "Enter the 5' overhang sequence for ssDNA bridging. This DNA will be added upstream to the sgRNA that is found and made into an oligo.",
        target_id="forward-overhang-tooltip"
    ),
    create_tooltip(
        "Enter the 3' overhang sequence for ssDNA bridging. This DNA will be added downstream to the sgRNA that is found and made into an oligo.",
        target_id="reverse-overhang-tooltip"
    ),

    # sgRNA Handle
    create_tooltip(
        "Specify the sgRNA handle sequence for the CRISPR system.",
        target_id="sgRNA-handle-tooltip"
    ),

    # Advanced Settings
    create_tooltip(
        "Show advanced settings for checking primers and repair templates.",
        target_id="show-advanced-settings-tooltip"
    ),


        # Advanced Settings
    create_tooltip(
        "Show advanced settings for checking primers and repair templates.",
        target_id="show-advanced-settings-tooltip"
    ),
    
    # Golden Gate 
    create_tooltip(
        "Specify the forward overhang for restriction-based cloning.",
        target_id="restriction-overhang-f-tooltip"
    ),
    create_tooltip(
        "Specify the reverse overhang for restriction-based cloning.",
        target_id="restriction-overhang-r-tooltip"
    ),
    create_tooltip(
        "Define the forward overhang sequence for the plasmid backbone.",
        target_id="backbone-overhang-f-tooltip"
    ),
    create_tooltip(
        "Define the reverse overhang sequence for the plasmid backbone.",
        target_id="backbone-overhang-r-tooltip"
    ),
    create_tooltip(
        "Enter the Csy4 sequence for multiplexing in CRISPR-BEST.",
        target_id="cys4-sequence-tooltip"
    ),
        # Melting Temperature
    create_tooltip(
        "Set the desired melting temperature for your checking primers PCR reactions.",
        target_id="checking-melting-temp-tooltip"
    ),

    # second restriction enzyme
    create_tooltip(
    "Specify the restriction enzyme used for integrating repair templates during cloning.",
    target_id="restriction-enzyme-tooltip"
    ),
    # Protospacer 5' Anneal
    create_tooltip(
        "Enter the forward 5' annealing sequence for the protospacer, used in Gibson assembly.",
        target_id="forward-overhang-tooltip-6"
    ),

    # Protospacer 3' Anneal
    create_tooltip(
        "Enter the reverse 3' annealing sequence for the protospacer, used in Gibson assembly.",
        target_id="reverse-overhang-tooltip-6"
    ),

    # Backbone 5' Anneal
    create_tooltip(
        "Standardized primer: Specify the forward 5' annealing sequence for the backbone in the Gibson assembly process.",
        target_id="backbone-forward-overhang-tooltip-6"
    ),

    # Backbone 3' Anneal
    create_tooltip(
        "Standardized primer: Specify the reverse 3' annealing sequence for the backbone in the Gibson assembly process.",
        target_id="backbone-reverse-overhang-tooltip-6"
    ),

]


