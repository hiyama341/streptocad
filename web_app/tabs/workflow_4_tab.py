from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable

from styling import text_style, upload_button_style, card_style, link_style, table_style
from streptocad.utils import polymerase_dict
from components import upload_component_with_display
from tooltip import create_tooltip, tooltips  # Import the tooltips

# Reference content used in multiple tabs
reference_content = html.Div([
    html.P("Note: For more information on CRISPR techniques and details, please visit the Nature protocols article below."),
    html.A("CRISPR–Cas9, CRISPRi and CRISPR-BEST-mediated genetic manipulation in streptomycetes", href="https://www.nature.com/articles/s41596-020-0339-z", target="_blank"),
], style=text_style)

# Dropdown options for polymerases
dropdown_options = [{'label': key, 'value': value} for key, value in polymerase_dict.items()]

# Tab content for CRISPRi
crispri_tab = html.Div(children=[
    dbc.Row([
        dbc.Col([
            dcc.Markdown("""
            ## **What is CRISPRi?**
            - CRISPR interference (CRISPRi) is a method that uses a catalytically dead Cas9 (dCas9) to block transcription of target genes.

            ## **Why use it?**
            - It allows for the precise and reversible silencing of genes without altering the DNA sequence.

            ## **Getting Started**
            - Find your plasmid of choice. We recommend that you use a dCas9 vector suitable for CRISPRi.
            - Figure out what genes you want to target. For example, the actinorhodin cluster (SCO5087).

            ## **Instructions**
            Upload your genome file and CRISPRi plasmid files, then click 'Submit' to generate your assembly.
            """, style=text_style),
            reference_content
        ], style={'padding': '20px', 'backgroundColor': '#2C3E50'}),
    ], className="mb-4"),

    html.Img(src='/assets/w4_pic.png', style={'width': '60%', 'margin': '20px auto'}),

    dbc.Row([
        dbc.Col(
            html.A(
                dbc.Button("Download Example Genome File", color="primary"),
                href="/assets/Streptomyces_coelicolor_A3_chromosome.gb",
                download="Streptomyces_coelicolor_A3_chromosome.gb"
            ),
            width={"size": 3, "order": 1}
        ),
        dbc.Col(
            html.A(
                dbc.Button("Download Example CRISPRi Plasmid File", color="primary", className="mb-4"),
                href="/assets/pCRISPR-dCas9.gbk",
                download="pCRISPR-dCas9.gbk"
            ),
            width={"size": 3, "order": 2}
        )
    ], className="mb-4", justify="start"),

    # UPLOAD genome
    dbc.Row([
        dbc.Col([
            html.H4("1) Upload your genome file", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Genome File (GenBank format)", className="card-title", style=text_style),
                    upload_component_with_display(
                        {'type': 'upload-component', 'index': 'genome-file-4'},
                        text_style, 
                        link_style, 
                        upload_button_style
                    )
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.H4("2) Upload the plasmid of choice", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("CRISPRi plasmid (GenBank format)", className="card-title", style=text_style),
                    upload_component_with_display(
                        {'type': 'upload-component', 'index': 'single-vector-4'},
                        text_style, 
                        link_style, 
                        upload_button_style
                    )
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-5"),
        
    dbc.Row([
        dbc.Col([
            html.H4("3) Choose genes/regions to knock down", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Example for genes: SCO5087, SCO5087,... (comma-separated)", className="card-title", style=text_style),
                    html.H5("Example for regions: 1000-2000, 100000-101000,... (comma-separated)", className="card-title", style=text_style),
                    dbc.Textarea(
                        id='genes-to-KO_4',
                        placeholder='Enter genes to knock out, e.g., SCO5087',
                        value='SCO5087',
                        style={'width': '100%', 'height': '100px'}
                    )
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.H4("4) Select overhangs", style=text_style),
            html.P("Please enter the 5' and 3' overhangs below for the oligo nucleotide to be made.", className="lead", style=text_style),
            html.P("Per default the overhangs work with the example pCRISPR-dCas9.gbk plasmid.", className="lead", style=text_style)
        ], width=12),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.Label(["5 prime Overhang: ", html.Span("ⓘ", id="forward-overhang-tooltip-4", style=link_style)], style=text_style),
            dbc.Input(
                id='forward-overhang-input_4',
                type='text',
                placeholder='Enter Forward Overhang',
                value='CGGTTGGTAGGATCGACGGC'
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
        create_tooltip("Enter the 5' overhang sequence for oligo nucleotide synthesis.", "forward-overhang-tooltip-4"),
    ]),

    dbc.Row([
        dbc.Col([
            html.Label(["3 prime Overhang: ", html.Span("ⓘ", id="reverse-overhang-tooltip-4", style=link_style)], style=text_style),
            dbc.Input(
                id='reverse-overhang-input_4',
                type='text',
                placeholder='Enter Reverse Overhang',
                value='GTTTTAGAGCTAGAAATAGC'
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
        create_tooltip("Enter the 3' overhang sequence for oligo nucleotide synthesis.", "reverse-overhang-tooltip-4"),
    ], className="mb-5"),
    
    dbc.Row([
        dbc.Col([
            html.H4("5) Filtering metrics for sgRNAs", style=text_style),
            dbc.Row([
                dbc.Col([
                    dbc.Label(["GC Content Upper Bound ", html.Span("ⓘ", id="gc-upper-tooltip-4", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='gc-upper_4',
                        type='number',
                        value=0.99,
                        style={'color': '#000', 'width': '100%'}
                    ),
                    create_tooltip("Set the upper bound for GC content in sgRNA sequences.", "gc-upper-tooltip-4")
                ], width=12, className="mb-3"),

                dbc.Col([
                    dbc.Label(["GC Content Lower Bound ", html.Span("ⓘ", id="gc-lower-tooltip-4", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='gc-lower_4',
                        type='number',
                        value=0.01,
                        style={'color': '#000', 'width': '100%'}
                    ),
                    create_tooltip("Set the lower bound for GC content in sgRNA sequences.", "gc-lower-tooltip-4")
                ], width=12, className="mb-3"),

                dbc.Col([
                    dbc.Label(["Off-Target Seed Length ", html.Span("ⓘ", id="off-target-seed-tooltip-4", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='off-target-seed_4',
                        type='number',
                        value=13,
                        style={'color': '#000', 'width': '100%'}
                    ),
                    create_tooltip("Specify the length of the seed sequence for off-target filtering.", "off-target-seed-tooltip-4")
                ], width=12, className="mb-3"),

                dbc.Col([
                    dbc.Label(["Off-Target Upper Bound ", html.Span("ⓘ", id="off-target-upper-tooltip-4", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='off-target-upper_4',
                        type='number',
                        value=10,
                        style={'color': '#000', 'width': '100%'}
                    ),
                    create_tooltip("Set the maximum allowed off-target score for sgRNA filtering.", "off-target-upper-tooltip-4")
                ], width=12, className="mb-3"),

                dbc.Col([
                    dbc.Label(["Cas Type ", html.Span("ⓘ", id="cas-type-tooltip-4", style=link_style)], style={'color': '#ddd'}),  
                    dcc.Dropdown(
                        id='cas-type_4',
                        options=[
                            {'label': 'Cas9', 'value': 'cas9'},
                        ],
                        value='cas9',
                        style={'color': '#000', 'width': '100%'}
                    ),
                    create_tooltip("Select the Cas protein type for the genome editing workflow.", "cas-type-tooltip-4")
                ], width=12, className="mb-3"),

                dbc.Col([
                    dbc.Label(["Number of sgRNAs per region/locus tag ", html.Span("ⓘ", id="number-sgrnas-tooltip-4", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='number-of-sgRNAs-per-group_4',
                        type='number',
                        value=5,
                        style={'color': '#000', 'width': '100%'}
                    ),
                    create_tooltip("Define the number of sgRNAs for each targeted region or locus.", "number-sgrnas-tooltip-4")
                ], width=12, className="mb-3"),

                dbc.Col([
                    dbc.Label(["Extension to Promoter Region ", html.Span("ⓘ", id="extension-to-promoter-region-tooltip", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='extension-to-promoter-region_4',
                        type='number',
                        value=100,
                        style={'color': '#000', 'width': '100%'}
                    ),
                    create_tooltip("Set the length for promoter region extension in gene silencing.", "extension-to-promoter-region-tooltip")
                ], width=12, className="mb-3"),
            ]),
        ], width=6, className="mb-4"),
    ], className="mb-3"),
    
    dbc.Col([
            dbc.Label(["Restriction enzyme(s)", html.Span("ⓘ", id="restriction-enzymes-tooltip-4", style=link_style)], style={'color': '#ddd'}),
            dbc.Input(
                id='restriction-enzymes_4',
                type='text',
                value='NcoI',
                style={'color': '#000', 'width': '100%'}
            ),
        ], width=6),

    dbc.Row([
        dbc.Col([
            dbc.Button('Submit', id='submit-settings-button_4', color="primary", className="mt-3"),
        ], width=12),
    ], className="mb-4"),
    
    ### Output
    dbc.Row([
        dbc.Col([
            dcc.Loading(
                id="loading-overlay-4",
                type="circle",
                style={'height': '80px', 'width': '80px', 'margin': 'auto'},  # Adjust size as needed
                children=[
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Filtered sgRNAs", className="card-title", style=text_style),
                            DataTable(id='mutated-sgrna-table_4', **table_style),
                        ])
                    ], style=card_style),

                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Primers", className="card-title", style=text_style),
                            DataTable(id='primer-table_4', **table_style),
                        ])
                    ], style=card_style),
                    
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Overview of plasmids generated", className="card-title", style=text_style),
                            DataTable(id='plasmid-metadata-table_4', **table_style),
                        ])
                    ], style=card_style),

                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Download folder with all data & protocols", className="card-title", style=text_style),
                            DataTable(id='all_data_4', **table_style),
                            html.A(
                                'Download All Data & protocols',
                                id='download-data-and-protocols-link_4',
                                download="data_package",
                                href="",
                                target="_blank",
                                className="btn btn-primary"
                            )
                        ])
                    ], style=card_style),
                ]
            )
        ], width=10),
    ]),

    # Include all tooltips
    *tooltips
])
