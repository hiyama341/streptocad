from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable
from styling import text_style, upload_button_style, card_style, link_style, table_style
from streptocad.utils import polymerase_dict
from tooltip import create_tooltip, tooltips

# Reference content
reference_content = html.Div([
    html.P("Note: For more information on CRISPR techniques and details, please visit the article below (the figure and protocol is from there).", style=text_style),
    html.A("CASCADE-Cas3 Enables Highly Efficient Genome Engineering in Streptomyces Species", href="https://www.biorxiv.org/content/10.1101/2023.05.09.539971v1.full", target="_blank", style=text_style),
], style=text_style)

# Dropdown options for polymerases
dropdown_options = [{'label': key, 'value': value} for key, value in polymerase_dict.items()]

# Tab content for CRISPR–Cas3 plasmid construction
cas3_tab = dcc.Tab(label="CRISPR–Cas3 plasmid construction", children=[
    dbc.Row([
        dbc.Col([
            dcc.Markdown("""
            ## **What is CRISPR-Cas3?**
            - A method that can be used to perform random-sized mutations or full in-frame deletions with repair templates.

            ## **Why use it?**
            - When you aim to delete a specific gene.

            ## **Getting Started**
            - Find your Cas3 plasmid (or download an example genome below).
            - Fetch your organism's genome (or download an example genome below).
            - Figure out what genes you want to target. For example, the actinorhodin cluster (SCO5087).    

            ## **Instructions**
            Upload your files, then click 'Submit' to generate your assembly.
            """, style=text_style),
            reference_content
        ], style={'padding': '20px', 'backgroundColor': '#2C3E50'}),
    ], className="mb-4"),
    
    html.Img(src='/assets/w6_pic.png', style={'width': '60%', 'margin': '20px auto'}),

    # Download example files
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
                dbc.Button("Download Example Plasmid File", color="primary", className="mb-4"),
                href="/assets/pCRISPR_cas3.gbk",
                download="pCRISPR_cas3.gbk"
            ),
            width={"size": 3, "order": 2}
        )
    ], className="mb-4", justify="start"),

    # Genome file upload
    dbc.Row([
        dbc.Col([
            html.H4("1) Upload your genome file", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Genome File (GenBank format)", className="card-title", style=text_style),
                    dcc.Upload(
                        id={'type': 'upload-component', 'index': 'genome-file-6'},
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Genome File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id={'type': 'filename-display', 'index': 'genome-file-6'}, children=[], style=text_style),
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-3"),

    # Plasmid file upload
    dbc.Row([
        dbc.Col([
            html.H4("2) Upload the plasmid of choice", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("CRISPR Plasmid (GenBank format)", className="card-title", style=text_style),
                    dcc.Upload(
                        id={'type': 'upload-component', 'index': 'single-vector-6'},
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select CRISPR Plasmid File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id={'type': 'filename-display', 'index': 'single-vector-6'}, children=[], style=text_style),
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-5"),

    # Genes/Regions selection
    dbc.Row([
        dbc.Col([
            html.H4("3) Choose genes/regions to delete", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Example for genes: SCO5087, SCO5087,... (comma-separated)", className="card-title", style=text_style),
                    html.H5("Example for regions: 1000-2000, 100000-101000,... (comma-separated)", className="card-title", style=text_style),
                    dbc.Textarea(
                        id='genes-to-KO_6',
                        placeholder='Enter genes to knock out, e.g., SCO5087',
                        value='SCO5087',
                        style={'width': '100%', 'height': '100px'}
                    )
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-3"),

    # Overhangs
    dbc.Row([
        dbc.Col([
            html.H4("4) Select overhangs", style=text_style),
            html.P("Please enter the 5' and 3' annealing sequences below for the first Gibson reaction.", className="lead", style=text_style),
            html.P("Per default, the overhangs work with pCRISPR–Cas3.gbk", className="lead", style=text_style)
        ], width=12),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.Label(['Protospacer - 5 prime anneal:', html.Span("ⓘ", id="forward-overhang-tooltip-6", style=link_style)], style=text_style),
            dbc.Input(
                id='forward-overhang-input_6',
                type='text',
                placeholder='Enter Forward Overhang',
                value='GTCGCcCggCaaAaccGg'.upper()
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Label(['Protospacer - 3 prime anneal:', html.Span("ⓘ", id="reverse-overhang-tooltip-6", style=link_style)], style=text_style),
            dbc.Input(
                id='reverse-overhang-input_6',
                type='text',
                placeholder='Enter Reverse Overhang',
                value='GTTTCAATCCACGCGCCCGT'
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
    ], className="mb-5"),

    # Backbone Overhangs
    dbc.Row([
        dbc.Col([
            html.Label(['Backbone - 5 prime anneal:', html.Span("ⓘ", id="backbone-forward-overhang-tooltip-6", style=link_style)], style=text_style),
            dbc.Input(
                id='backbone-forward-overhang-input_6',
                type='text',
                placeholder='Enter Forward Overhang',
                value='GAGCTCATAAGTTCCTATTCCGAAG'
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Label(['Backbone - 3 prime anneal:', html.Span("ⓘ", id="backbone-reverse-overhang-tooltip-6", style=link_style)], style=text_style),
            dbc.Input(
                id='backbone-reverse-overhang-input_6',
                type='text',
                placeholder='Enter Reverse Overhang',
                value='aagaagtgggtgtcggacgc'.upper()
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
    ], className="mb-5"),

    # Filtering metrics for sgRNAs
    dbc.Row([
        dbc.Col([
            html.H4("5) Filtering metrics for sgRNAs", style=text_style),
            dbc.Row([
                dbc.Col([
                    dbc.Label(["GC Content Upper Bound", html.Span("ⓘ", id="gc-upper-tooltip", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='gc-upper_6',
                        type='number',
                        value=0.99,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label(["GC Content Lower Bound", html.Span("ⓘ", id="gc-lower-tooltip", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='gc-lower_6',
                        type='number',
                        value=0.01,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label(["Off-Target Seed Length", html.Span("ⓘ", id="off-target-seed-tooltip", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='off-target-seed_6',
                        type='number',
                        value=13,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label(["Off-Target Upper Bound", html.Span("ⓘ", id="off-target-upper-tooltip", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='off-target-upper_6',
                        type='number',
                        value=10,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label(["Cas Type", html.Span("ⓘ", id="cas-type-tooltip", style=link_style)], style={'color': '#ddd'}),  
                    dcc.Dropdown(
                        id='cas-type_6',
                        options=[{'label': 'cas3', 'value': 'cas3'}],
                        value='cas3',
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label(["Number of sgRNAs per region/locus tag", html.Span("ⓘ", id="number-sgrnas-tooltip", style=link_style)], style={'color': '#ddd'}),  
                    dbc.Input(
                        id='number-of-sgRNAs-per-group_6',
                        type='number',
                        value=5,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
            ]),
        ], width=6, className="mb-4"),
    ], className="mb-3"),

    # Switch for Advanced Settings
    dbc.Row([
        dbc.Col([
            html.H4("6) Show advanced settings for checking primers and repair templates", style=text_style),
            dbc.Checklist(
                options=[{"label": "", "value": 1}],
                value=[],
                id="show-advanced-settings-checkbox",
                inline=True,
                switch=True,
                className="big-switch"
            ),
        ], width=6),
    ], className="mb-3"),

    # Advanced Settings
    dbc.Row([
        dbc.Col([
            html.Div(id='advanced-settings-container', children=[
                dbc.Row([
                    dbc.Col([
                        dbc.Label(["Choose Polymerase", html.Span("ⓘ", id="polymerase-tooltip", style=link_style)], style={'color': '#ddd'}),
                        dcc.Dropdown(
                            id='chosen-polymerase_6',
                            options=dropdown_options,
                            value=polymerase_dict['Q5 High-Fidelity 2X Master Mix'],
                            style={'color': '#000'}
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label(["Target Melting Temperature (°C)", html.Span("ⓘ", id="melting-temp-tooltip", style=link_style)], style={'color': '#ddd'}),
                        dbc.Input(
                            id='melting-temperature_6',
                            type='number',
                            value=65,
                            style={'color': '#000'}
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label(["Primer Concentration (μM)", html.Span("ⓘ", id="primer-concentration-tooltip", style=link_style)], style={'color': '#ddd'}),
                        dbc.Input(
                            id='primer-concentration_6',
                            type='number',
                            value=0.4,
                            style={'color': '#000'}
                        ),
                    ], width=6),

                    dbc.Col([
                        dbc.Label(["Flanking Region Number", html.Span("ⓘ", id="flanking-region-tooltip", style=link_style)], style={'color': '#ddd'}),
                        dbc.Input(
                            id='flanking-region-number_6',
                            type='number',
                            value=500,
                            style={'color': '#000'}
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label(["Length of repair templates", html.Span("ⓘ", id="repair-templates-length-tooltip", style=link_style)], style={'color': '#ddd'}),
                        dbc.Input(
                            id='repair_templates_length_6',
                            type='number',
                            value=1000,
                            style={'color': '#000'}
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label(["Overlap between templates for Gibson cloning", html.Span("ⓘ", id="overlap-gibson-tooltip", style=link_style)], style={'color': '#ddd'}),
                        dbc.Input(
                            id='overlap_for_gibson_length_6',
                            type='number',
                            value=40,
                            style={'color': '#000'}
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label(["Restriction enzyme for integration of repair templates", html.Span("ⓘ", id="restriction-enzymes-tooltip-6", style=link_style)], style={'color': '#ddd'}),
                        dbc.Input(
                            id='restriction_enzyme_for_repair_templates_integration_6',
                            type='text',
                            value='EcoRI',
                            style={'color': '#000'}
                        ),
                    ], width=6),
                ])
            ], style={"display": "none"})  # Hidden by default
        ], width=6),
    ], className="mb-3"),

    # Switch for In-Frame Deletions
    dbc.Row([
        dbc.Col([
            html.H4("7) Generate in-frame deletions", style=text_style),
            html.P("Note: adding repair templates to your plasmids", style={'color': '#ddd', 'fontSize': '1rem'}),
            dbc.Checklist(
                options=[{"label": "", "value": 1}],
                value=[],
                id="show-inframe-deletions-settings-checkbox_6",
                inline=True,
                switch=True,
                className="big-switch"
            ),
        ], width=10),
    ], className="mb-3"),

    # Submit Button
    dbc.Row([
        dbc.Col([
            dbc.Button('Submit', id='submit-settings-button_6', color="primary", className="mt-3",size="lg"),
        ], width=12),
    ], className="mb-4"),

    # Output Tables
    dbc.Row([
        dbc.Col([
            dcc.Loading(
                id="loading-overlay",
                type="circle",
                children=[
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Filtered sgRNAs", className="card-title", style=text_style),
                            DataTable(id='filtered-df-table_6', **table_style),
                        ])
                    ], style=card_style),
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Primers Output Table", className="card-title", style=text_style),
                            DataTable(id='primers-output-table_6', **table_style),
                        ])
                    ], style=card_style),
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("PCR Table", className="card-title", style=text_style),
                            DataTable(id='pcr-table_6', **table_style),
                        ])
                    ], style=card_style),
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Overview of Plasmids", className="card-title", style=text_style),
                            DataTable(id='plasmid-metadata-table_6', **table_style),
                        ])
                    ], style=card_style),
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Download All Data & Protocols", className="card-title", style=text_style),
                            DataTable(id='all_data_6', **table_style),
                            html.A(
                                'Download All Data & Protocols',
                                id='download-data-and-protocols-link_6',
                                download="all_data",
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
    # Attach tooltips to the layout
    *tooltips
])
