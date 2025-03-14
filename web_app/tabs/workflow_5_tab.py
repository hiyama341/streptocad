from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable

from styling import (
    text_style,
    upload_button_style,
    card_style,
    link_style,
    table_style,
    table_header_style,
    table_row_style,
)
from streptocad.utils import polymerase_dict
from tooltip import create_tooltip, tooltips  # Import tooltips
from components import upload_component_with_display


# Reference content used in multiple tabs
reference_content = html.Div(
    [
        html.P(
            "Note: For more information on CRISPR techniques and details, please visit the Nature protocols article below."
        ),
        html.A(
            "CRISPR–Cas9, CRISPRi and CRISPR-BEST-mediated genetic manipulation in streptomycetes",
            href="https://www.nature.com/articles/s41596-020-0339-z",
            target="_blank",
            style=link_style,
        ),
    ],
    style=text_style,
)

# Dropdown options for polymerases
dropdown_options = [
    {"label": key, "value": value} for key, value in polymerase_dict.items()
]

# Tab content for CRISPR–Cas9 plasmid construction for in-frame deletion
gibson_tab = dcc.Tab(
    label="CRISPR–Cas9 plasmid construction",
    children=[
        dbc.Row(
            [
                dbc.Col(
                    [
                        dcc.Markdown(
                            """
            ## **What is CRISPR-Cas9?**
            - A method that can be used to perform random-sized mutations or full in-frame deletions with repair templates.

            ## **Why use it?**
            - When you aim to delete a specific gene.

            ## **Getting Started**
            - Find your plasmids (or download an example below).
            - Fetch your organism's genome (or download an example below).
            - Figure out what genes you want to target. For example, the actinorhodin cluster (SCO5087).

            ## **Instructions**
            - Upload your genome file and CRISPR plasmid files
            - Click 'Submit' to generate your assembly.
            """,
                            style=text_style,
                        ),
                        reference_content,
                    ],
                    style={"padding": "20px", "backgroundColor": "#2C3E50"},
                ),
            ],
            className="mb-4",
        ),
        html.Img(
            src="/assets/w5_pic.png", style={"width": "60%", "margin": "20px auto"}
        ),
        # Section 1: Upload Genome File with Download Button
        dbc.Row(
            [
                dbc.Col(
                    dbc.Row(
                        [
                            # Header Text
                            dbc.Col(
                                html.H4("1) Upload your genome file", style=text_style),
                                width="auto",
                                className="align-self-center",
                            ),
                            # Download Button
                            dbc.Col(
                                html.A(
                                    dbc.Button(
                                        "Download Example Genome File", color="primary"
                                    ),
                                    href="/assets/Streptomyces_coelicolor_A3_chromosome.gb",
                                    download="Streptomyces_coelicolor_A3_chromosome.gb",
                                ),
                                width="auto",
                                className="align-self-center ml-3",
                            ),
                        ]
                    ),
                    width=6,
                    className="d-flex align-items-center",
                ),
            ],
            className="mb-4",
        ),  # Margin bottom for spacing
        # Upload Component for Genome File
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Card(
                            [
                                dbc.CardBody(
                                    [
                                        html.H5(
                                            "Genome File (GenBank format)",
                                            className="card-title",
                                            style=text_style,
                                        ),
                                        upload_component_with_display(
                                            {
                                                "type": "upload-component",
                                                "index": "genome-file-5",
                                            },
                                            text_style,
                                            link_style,
                                            upload_button_style,
                                        ),
                                    ]
                                )
                            ],
                            style=card_style,
                        ),
                    ],
                    width=6,
                ),
            ],
            className="mb-5",
        ),  # Margin bottom for spacing
        # Section 2: Upload Plasmid with Download Button
        dbc.Row(
            [
                dbc.Col(
                    dbc.Row(
                        [
                            # Header Text
                            dbc.Col(
                                html.H4(
                                    "2) Upload the plasmid of choice", style=text_style
                                ),
                                width="auto",
                                className="align-self-center",
                            ),
                            # Download Button
                            dbc.Col(
                                html.A(
                                    dbc.Button(
                                        "Download Example CRISPR Plasmid File",
                                        color="primary",
                                    ),
                                    href="/assets/pCRISPR-Cas9.gbk",
                                    download="pCRISPR-Cas9.gbk",
                                ),
                                width="auto",
                                className="align-self-center ml-3",
                            ),
                        ]
                    ),
                    width=6,
                    className="d-flex align-items-center",
                ),
            ],
            className="mb-4",
        ),  # Margin bottom for spacing
        # Upload Component for Plasmid
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Card(
                            [
                                dbc.CardBody(
                                    [
                                        html.H5(
                                            "CRISPR Plasmid (GenBank format)",
                                            className="card-title",
                                            style=text_style,
                                        ),
                                        upload_component_with_display(
                                            {
                                                "type": "upload-component",
                                                "index": "single-vector-5",
                                            },
                                            text_style,
                                            link_style,
                                            upload_button_style,
                                        ),
                                    ]
                                )
                            ],
                            style=card_style,
                        ),
                    ],
                    width=6,
                ),
            ],
            className="mb-5",
        ),  # Margin bottom for spacing
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4("3) Choose genes/regions to delete", style=text_style),
                        dbc.Card(
                            [
                                dbc.CardBody(
                                    [
                                        html.H5(
                                            "Example for genes: SCO5087, SCO5087,... (comma-separated)",
                                            className="card-title",
                                            style=text_style,
                                        ),
                                        html.H5(
                                            "Example for regions: 1000-2000, 100000-101000,... (comma-separated)",
                                            className="card-title",
                                            style=text_style,
                                        ),
                                        dbc.Textarea(
                                            id="genes-to-KO_5",
                                            placeholder="Enter genes to knock out, e.g., SCO5087",
                                            value="SCO5087",
                                            style={"width": "100%", "height": "100px"},
                                        ),
                                    ]
                                )
                            ],
                            style=card_style,
                        ),
                    ],
                    width=6,
                ),
            ],
            className="mb-3",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4(
                            "4) Select overhangs for oligo bridging", style=text_style
                        ),
                        html.P(
                            "Please enter the 5' and 3' overhangs below for the oligo nucleotide to be made.",
                            className="lead",
                            style=text_style,
                        ),
                        html.P(
                            "Per default the overhangs work with pCRISPR–Cas9_plasmid_addgene.gbk",
                            className="lead",
                            style=text_style,
                        ),
                    ],
                    width=12,
                ),
            ],
            className="mb-3",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label(
                            [
                                "5 prime Overhang: ",
                                html.Span(
                                    "ⓘ",
                                    id="forward-overhang-tooltip-5",
                                    style=link_style,
                                ),
                            ],
                            style=text_style,
                        ),
                        dbc.Input(
                            id="forward-overhang-input_5",
                            type="text",
                            placeholder="Enter Forward Overhang",
                            value="CGGTTGGTAGGATCGACGGC",
                        ),
                    ],
                    width=6,
                    style={"marginRight": "10px", "marginLeft": "10px"},
                ),
                create_tooltip(
                    "Enter the 5' overhang sequence for oligo nucleotide synthesis.",
                    "forward-overhang-tooltip-5",
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Label(
                            [
                                "3 prime Overhang: ",
                                html.Span(
                                    "ⓘ",
                                    id="reverse-overhang-tooltip-5",
                                    style=link_style,
                                ),
                            ],
                            style=text_style,
                        ),
                        dbc.Input(
                            id="reverse-overhang-input_5",
                            type="text",
                            placeholder="Enter Reverse Overhang",
                            value="GTTTTAGAGCTAGAAATAGC",
                        ),
                    ],
                    width=6,
                    style={"marginRight": "10px", "marginLeft": "10px"},
                ),
                create_tooltip(
                    "Enter the 3' overhang sequence for oligo nucleotide synthesis.",
                    "reverse-overhang-tooltip-5",
                ),
            ],
            className="mb-5",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4("5) Filtering metrics for sgRNAs", style=text_style),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "GC Content Upper Bound ",
                                                html.Span(
                                                    "ⓘ",
                                                    id="gc-upper-tooltip-5",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="gc-upper_5",
                                            type="number",
                                            value=0.99,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                        create_tooltip(
                                            "Set the upper bound for GC content in sgRNA sequences.",
                                            "gc-upper-tooltip-5",
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "GC Content Lower Bound ",
                                                html.Span(
                                                    "ⓘ",
                                                    id="gc-lower-tooltip-5",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="gc-lower_5",
                                            type="number",
                                            value=0.01,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                        create_tooltip(
                                            "Set the lower bound for GC content in sgRNA sequences.",
                                            "gc-lower-tooltip-5",
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "Off-Target Seed Length ",
                                                html.Span(
                                                    "ⓘ",
                                                    id="off-target-seed-tooltip-5",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="off-target-seed_5",
                                            type="number",
                                            value=13,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                        create_tooltip(
                                            "Specify the length of the seed sequence for off-target filtering.",
                                            "off-target-seed-tooltip-5",
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "Off-Target Upper Bound ",
                                                html.Span(
                                                    "ⓘ",
                                                    id="off-target-upper-tooltip-5",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="off-target-upper_5",
                                            type="number",
                                            value=10,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                        create_tooltip(
                                            "Specify the maximum allowed off-target score for sgRNA filtering.",
                                            "off-target-upper-tooltip-5",
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "Cas Type ",
                                                html.Span(
                                                    "ⓘ",
                                                    id="cas-type-tooltip-5",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dcc.Dropdown(
                                            id="cas-type_5",
                                            options=[
                                                {"label": "Cas9", "value": "cas9"},
                                            ],
                                            value="cas9",
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                        create_tooltip(
                                            "Select the Cas protein type for the genome editing workflow.",
                                            "cas-type-tooltip-5",
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "Number of sgRNAs per region/locus tag ",
                                                html.Span(
                                                    "ⓘ",
                                                    id="number-sgrnas-tooltip-5",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="number-of-sgRNAs-per-group_5",
                                            type="number",
                                            value=5,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                        create_tooltip(
                                            "Define the number of sgRNAs for each targeted region or locus.",
                                            "number-sgrnas-tooltip-5",
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                            ]
                        ),
                    ],
                    width=6,
                    className="mb-4",
                ),
            ],
            className="mb-3",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4("6) Show advanced settings", style=text_style),
                        dbc.Checklist(
                            options=[{"label": "", "value": 1}],
                            value=[],
                            id="show-advanced-settings-checkbox",
                            inline=True,
                            switch=True,
                            className="big-switch",
                        ),
                    ],
                    width=6,
                ),
            ],
            className="mb-3",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Div(
                            id="advanced-settings-container",
                            children=[
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Choose Polymerase ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="polymerase-tooltip-5",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dcc.Dropdown(
                                                    id="chosen-polymerase_5",
                                                    options=dropdown_options,
                                                    value=polymerase_dict[
                                                        "Q5 High-Fidelity 2X Master Mix"
                                                    ],
                                                    style={"color": "#000"},
                                                ),
                                                create_tooltip(
                                                    "Select the polymerase type for high-fidelity DNA synthesis.",
                                                    "polymerase-tooltip-5",
                                                ),
                                            ],
                                            width=12,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Target Melting Temperature (°C) ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="melting-temp-tooltip-5",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="melting-temperature_5",
                                                    type="number",
                                                    value=65,
                                                    style={"color": "#000"},
                                                ),
                                                create_tooltip(
                                                    "Set the desired melting temperature for the PCR reactions.",
                                                    "melting-temp-tooltip-5",
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Primer Concentration (μM) ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="primer-concentration-tooltip-5",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="primer-concentration_5",
                                                    type="number",
                                                    value=0.4,
                                                    style={"color": "#000"},
                                                ),
                                                create_tooltip(
                                                    "Specify the concentration of primers used in the PCR reaction.",
                                                    "primer-concentration-tooltip-5",
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        # Primer length
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Minimum Checking Primer Anneal Length",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="primer-anneal-length-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="primer-lenght_5",
                                                    type="number",
                                                    value=18,
                                                    style={
                                                        "color": "#000",
                                                        "width": "100%",
                                                    },
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Flanking Region Number ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="flanking-region-tooltip-5",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="flanking-region-number_5",
                                                    type="number",
                                                    value=500,
                                                    style={"color": "#000"},
                                                ),
                                                create_tooltip(
                                                    "Define the size (in bp) of flanking regions for PCR reactions.",
                                                    "flanking-region-tooltip-5",
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Length of repair templates ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="repair-templates-length-tooltip-5",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="repair_templates_length_5",
                                                    type="number",
                                                    value=1000,
                                                    style={"color": "#000"},
                                                ),
                                                create_tooltip(
                                                    "Specify the length of repair templates for precise deletions.",
                                                    "repair-templates-length-tooltip-5",
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Overlap between templates for Gibson Cloning ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="overlap-gibson-tooltip-5",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="overlap_for_gibson_length_5",
                                                    type="number",
                                                    value=40,
                                                    style={"color": "#000"},
                                                ),
                                                create_tooltip(
                                                    "Define the overlap length between templates for Gibson assembly.",
                                                    "overlap-gibson-tooltip-5",
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "First digestion with restriction enzyme(s)",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="restriction-enzymes-tooltip-5",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="restriction-enzymes_5",
                                                    type="text",
                                                    value="NcoI",
                                                    style={
                                                        "color": "#000",
                                                        "width": "100%",
                                                    },
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Second digestion with restriction enzyme(s) for integration of repair templates",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="restriction-enzymes-tooltip-5-2",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="restriction-enzymes_5-2",
                                                    type="text",
                                                    value="StuI",
                                                    style={
                                                        "color": "#000",
                                                        "width": "100%",
                                                    },
                                                ),
                                            ],
                                            width=6,
                                        ),
                                    ]
                                )
                            ],
                            style={"display": "none"},
                        )
                    ],
                    width=6,
                ),
            ],
            className="mb-3",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4("7) Generate Repair Templates", style=text_style),
                        html.P(
                            "Note: adding repair templates to your plasmids (default 1000 bp)",
                            style={"color": "#ddd", "fontSize": "1rem"},
                        ),
                        dbc.Checklist(
                            options=[{"label": "", "value": 1}],
                            value=[],
                            id="show-inframe-deletions-settings-checkbox_5",
                            inline=True,
                            switch=True,
                            className="big-switch",
                        ),
                    ],
                    width=10,
                ),
            ],
            className="mb-3",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Button(
                            "Submit",
                            id="submit-settings-button_5",
                            color="primary",
                            className="mt-3",
                            size="lg",
                        ),
                    ],
                    width=12,
                ),
            ],
            className="mb-4",
        ),
        # Placeholder for the output
        dbc.Row(
            [
                dbc.Col(
                    [
                        dcc.Loading(
                            id="loading-overlay",
                            type="circle",
                            style={
                                "height": "80px",
                                "width": "80px",
                                "margin": "auto",
                            },  # Adjust size as needed
                            children=[
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(
                                                    "Filtered sgRNAs",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="filtered-df-table_5",
                                                    **table_style,
                                                ),
                                            ]
                                        )
                                    ],
                                    style=card_style,
                                ),
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(
                                                    "Primers Output Table",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="primers-output-table_5",
                                                    **table_style,
                                                ),
                                            ]
                                        )
                                    ],
                                    style=card_style,
                                ),
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(
                                                    "PCR Table",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="pcr-table_5", **table_style
                                                ),
                                            ]
                                        )
                                    ],
                                    style=card_style,
                                ),
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(
                                                    "Overview of plasmids",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="plasmid-metadata-table_5",
                                                    **table_style,
                                                ),
                                            ]
                                        )
                                    ],
                                    style=card_style,
                                ),
                                dbc.Card(
                                    [
                                        dbc.CardBody(
                                            [
                                                html.H5(
                                                    "Download All Data & Protocols",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(id="all_data", **table_style),
                                                html.A(
                                                    "Download All Data & Protocols",
                                                    id="download-data-and-protocols-link_5",
                                                    download="all_data",
                                                    href="",
                                                    target="_blank",
                                                    className="btn btn-primary",
                                                ),
                                            ]
                                        )
                                    ],
                                    style=card_style,
                                ),
                            ],
                        )
                    ],
                    width=10,
                ),
            ]
        ),
        # Include all tooltips
        *tooltips,
    ],
)
