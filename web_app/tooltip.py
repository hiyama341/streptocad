# tooltips.py
import dash_bootstrap_components as dbc

def create_tooltip(content, target_id, placement="top", style=None):
    return dbc.Tooltip(content, target=target_id, placement=placement, style=style)

# Generalized Tooltips for Common Use Cases
tooltips = [
    # Polymerase
    create_tooltip(
        "The selected polymerase is used for high-fidelity DNA synthesis.",
        target_id="polymerase-tooltip"
    ),
    
    # Melting Temperature
    create_tooltip(
        "Set the desired melting temperature for PCR reactions.",
        target_id="melting-temp-tooltip"
    ),

    # Primer Concentration
    create_tooltip(
        "Specify the concentration of primers (in μM) for PCR reactions.",
        target_id="primer-concentration-tooltip"
    ),

    # Primer Number Increment
    create_tooltip(
        "Adjust the numbering increment for primer sets.",
        target_id="primer-increment-tooltip"
    ),

    # Flanking Region
    create_tooltip(
        "Define the size (in bp) of flanking regions for PCR reactions.",
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
        "Specify the length of the seed sequence for off-target filtering.",
        target_id="off-target-seed-tooltip"
    ),
    create_tooltip(
        "Specify the maximum allowed off-target score for sgRNA filtering.",
        target_id="off-target-upper-tooltip"
    ),

    # Cas Type
    create_tooltip(
        "Select the Cas protein type for the genome editing workflow.",
        target_id="cas-type-tooltip"
    ),

    # Number of sgRNAs per Region
    create_tooltip(
        "Define the number of sgRNAs for each targeted region or locus.",
        target_id="number-sgrnas-tooltip"
    ),

    # Only Stop Codons
    create_tooltip(
        "Filter to include only stop codons in the sgRNA analysis.",
        target_id="only-stop-codons-tooltip"
    ),

    # Editing Context
    create_tooltip(
        "Enable sequence context editing for greater precision in targeting.",
        target_id="editing-context-tooltip"
    ),

    # Repair Template and Overlap for Gibson Cloning
    create_tooltip(
        "Specify the length of repair templates for precise deletions.",
        target_id="repair-templates-length-tooltip"
    ),
    create_tooltip(
        "Define the overlap length between templates for Gibson assembly.",
        target_id="overlap-gibson-tooltip"
    ),

    # Overhangs (Forward and Reverse)
    create_tooltip(
        "Enter the 5' overhang sequence for Golden Gate or Gibson assembly.",
        target_id="forward-overhang-tooltip"
    ),
    create_tooltip(
        "Enter the 3' overhang sequence for Golden Gate or Gibson assembly.",
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


