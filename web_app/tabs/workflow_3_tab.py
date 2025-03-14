from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable

from styling import text_style, upload_button_style, card_style, link_style, table_style
from streptocad.utils import polymerase_dict
from components import upload_component_with_display
from tooltip import create_tooltip, tooltips  # Import the tooltips

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

# Tab content for Multiple sgRNA-integration
golden_gate_tab = dcc.Tab(
    label="Multiple sgRNA-integration",
    children=[
        dbc.Row(
            [
                dbc.Col(
                    [
                        dcc.Markdown(
                            """
            ## **What is multiplexing of CRISPR-BEST?**
            - A method that employs multiple single guide RNAs (sgRNAs) to target several genomic locations simultaneously.

            ## **Why use it?**
            - When you're looking to study or modify multiple genes or pathways in one go, this is your go-to approach.

            ## **Getting Started**
            - Find your plasmid of choice. We recommend that you use pCRISPR-MCBE_Csy4_kasopGFP.gb.
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
            src="/assets/w3_pic.png", style={"width": "60%", "margin": "20px auto"}
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
                                                "index": "genome-file-3",
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
                                    href="/assets/pCRISPR-MCBE_Csy4_kasopGFP.gb",
                                    download="pCRISPR-MCBE_Csy4_kasopGFP.gb",
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
                                            "CRISPR-mcBEST plasmid (GenBank format)",
                                            className="card-title",
                                            style=text_style,
                                        ),
                                        upload_component_with_display(
                                            {
                                                "type": "upload-component",
                                                "index": "single-vector-3",
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
                        html.H4("3) Choose genes/regions to edit", style=text_style),
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
                                            id="genes-to-KO_3",
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
                            [
                                "4) Enter sgRNA handle ",
                                html.Span(
                                    "ⓘ", id="sgRNA-handle-tooltip", style=link_style
                                ),  # Tooltip icon
                            ],
                            style=text_style,
                        ),
                        html.P(
                            "As a default, we use the sgRNA handle with Csy4 (see figure above)",
                            style=text_style,
                        ),
                        dbc.Input(
                            id="sgRNA-handle-input_3",
                            type="text",
                            placeholder="Enter sgRNA handle: e.g., 110-nt (pJET1.2–sgRNA handle)",
                            value="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTTCACTGCCGTATAGGCAGCTAAGAAA",
                        ),
                    ],
                    width=6,
                ),
            ],
            className="mb-5",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4(
                            [
                                "5) Enter Target Melting Temperature ",
                                html.Span(
                                    "ⓘ", id="melting-temp-tooltip", style=link_style
                                ),  # Tooltip icon
                            ],
                            className="mb-3",
                            style=text_style,
                        ),
                        dbc.Input(
                            id="input-tm_3",
                            type="number",
                            placeholder="Enter Melting Temperature",
                            value=60,
                        ),
                        html.Br(),
                    ],
                    width=6,
                ),
            ],
            className="mb-5",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4("6) Filtering metrics for sgRNAs", style=text_style),
                        dbc.Row(
                            [
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "GC Content Upper Bound",
                                                html.Span(
                                                    "ⓘ",
                                                    id="gc-upper-tooltip",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="gc-upper_3",
                                            type="number",
                                            value=0.99,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "GC Content Lower Bound",
                                                html.Span(
                                                    "ⓘ",
                                                    id="gc-lower-tooltip",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="gc-lower_3",
                                            type="number",
                                            value=0.01,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "Off-Target Seed Length",
                                                html.Span(
                                                    "ⓘ",
                                                    id="off-target-seed-tooltip",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="off-target-seed_3",
                                            type="number",
                                            value=13,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "Off-Target Upper Bound",
                                                html.Span(
                                                    "ⓘ",
                                                    id="off-target-upper-tooltip",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="off-target-upper_3",
                                            type="number",
                                            value=10,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "Cas Type",
                                                html.Span(
                                                    "ⓘ",
                                                    id="cas-type-tooltip",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dcc.Dropdown(
                                            id="cas-type_3",
                                            options=[
                                                {"label": "Cas9", "value": "cas9"}
                                            ],
                                            value="cas9",
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Label(
                                            [
                                                "Number of sgRNAs per region/locus tag",
                                                html.Span(
                                                    "ⓘ",
                                                    id="number-sgrnas-tooltip",
                                                    style=link_style,
                                                ),
                                            ],
                                            style={"color": "#ddd"},
                                        ),
                                        dbc.Input(
                                            id="number-of-sgRNAs-per-group_3",
                                            type="number",
                                            value=5,
                                            style={"color": "#000", "width": "100%"},
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Checklist(
                                            options=[
                                                {
                                                    "label": "Only Stop Codons",
                                                    "value": 1,
                                                }
                                            ],
                                            value=[1],
                                            id="only-stop-codons-checkbox_3",
                                            inline=True,
                                            switch=True,
                                            className="big-switch",
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                                dbc.Col(
                                    [
                                        dbc.Checklist(
                                            options=[
                                                {
                                                    "label": html.Span(
                                                        [
                                                            "Editing Sequence Context ",
                                                            html.Span(
                                                                "ⓘ",
                                                                id="editing-context-tooltip",
                                                                style=link_style,
                                                            ),
                                                        ]
                                                    ),
                                                    "value": 1,
                                                }
                                            ],
                                            value=[1],
                                            id="editing_context_3",
                                            inline=True,
                                            switch=True,
                                            className="big-switch",
                                        ),
                                    ],
                                    width=12,
                                    className="mb-3",
                                ),
                            ]
                        ),
                    ],
                    width=6,
                ),
            ],
            className="mb-4",
        ),
        # Advanced settings section with one button to toggle visibility
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4("7) Show advanced settings", style=text_style),
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
                                                        "Choose Polymerase",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="polymerase-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dcc.Dropdown(
                                                    id="chosen-polymerase_3",
                                                    options=dropdown_options,
                                                    value=polymerase_dict[
                                                        "Q5 High-Fidelity 2X Master Mix"
                                                    ],
                                                    style={"color": "#000"},
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Target Melting Temperature (°C)",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="checking-melting-temp-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="melting-temperature_3",
                                                    type="number",
                                                    value=65,
                                                    style={"color": "#000"},
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Primer Concentration (μM)",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="primer-concentration-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="primer-concentration_3",
                                                    type="number",
                                                    value=0.4,
                                                    style={"color": "#000"},
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
                                                    id="primer-lenght_3",
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
                                                        "Flanking Region Number",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="flanking-region-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="flanking-region-number_3",
                                                    type="number",
                                                    value=500,
                                                    style={"color": "#000"},
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Restriction Overhang Forward ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="restriction-overhang-f-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="restriction-overhang-f",
                                                    type="text",
                                                    value="GATCGggtctcc",
                                                    style={"color": "#000"},
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Restriction Overhang Reverse ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="restriction-overhang-r-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="restriction-overhang-r",
                                                    type="text",
                                                    value="GATCAGGTCTCg",
                                                    style={"color": "#000"},
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Backbone Overhang Forward ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="backbone-overhang-f-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="backbone-overhang-f",
                                                    type="text",
                                                    value="cATG",
                                                    style={"color": "#000"},
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Backbone Overhang Reverse ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="backbone-overhang-r-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="backbone-overhang-r",
                                                    type="text",
                                                    value="cTAG",
                                                    style={"color": "#000"},
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Csy4 Sequence ",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="cys4-sequence-tooltip",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="cys4-sequence",
                                                    type="text",
                                                    value="gTTCACTGCCGTATAGGCAGCTAAGAAA",
                                                    style={"color": "#000"},
                                                ),
                                            ],
                                            width=6,
                                        ),
                                        dbc.Col(
                                            [
                                                dbc.Label(
                                                    [
                                                        "Restriction enzyme(s)",
                                                        html.Span(
                                                            "ⓘ",
                                                            id="restriction-enzymes-tooltip-3",
                                                            style=link_style,
                                                        ),
                                                    ],
                                                    style={"color": "#ddd"},
                                                ),
                                                dbc.Input(
                                                    id="restriction-enzymes_3",
                                                    type="text",
                                                    value="NcoI,NheI",
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
            className="mb-4",
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dbc.Button(
                            "Submit",
                            id="submit-button_3",
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
                            id="loading-overlay-3",
                            type="circle",
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
                                                    id="mutated-sgrna-table_3",
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
                                                    "Primers",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="primer-table_3", **table_style
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
                                                    "Overhang",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="overhang-table_3", **table_style
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
                                                    "PCR",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="pcr-table_3", **table_style
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
                                                    "Overview of plasmids generated",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="plasmid-metadata-table_3",
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
                                                    "Download folder with all data & protocols",
                                                    className="card-title",
                                                    style=text_style,
                                                ),
                                                DataTable(
                                                    id="all_data_4", **table_style
                                                ),
                                                html.A(
                                                    "Download All Data & protocols",
                                                    id="download-data-and-protocols-link_3",
                                                    download="data_package",
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
        # Include the tooltips at the end of the layout
        *tooltips,
    ],
)
