from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable
from components import upload_component


from styling import text_style, upload_button_style, card_style, link_style, table_style, table_header_style, table_row_style
from streptocad.utils import polymerase_dict


# Reference content used in multiple tabs
reference_content = html.Div([
    html.P("Note: For more information on CRISPR techniques and details, please visit the Nature protocols article below (the figure and protocol is from there).", style=text_style),
    html.A("CRISPR–Cas9, CRISPRi and CRISPR-BEST-mediated genetic manipulation in streptomycetes", href="https://www.nature.com/articles/s41596-020-0339-z", target="_blank", style=text_style),
], style=text_style)

# Dropdown options for polymerases
dropdown_options = [{'label': key, 'value': value} for key, value in polymerase_dict.items()]

# Tab content for CRISPR-cBEST workflow
crispr_cb_tab = html.Div(children=[
    dbc.Row([
        dbc.Col([
            dcc.Markdown("""
            ## **What is it?**
            - A precise method that uses a single sgRNA to target a specific genomic location for base-editing.

            ## **Why use it?**
            - Perfect for those times when you have a single gene in your crosshairs.

            ## **Getting Started**
            - Find your plasmid of choice. We recommend that you use pCRISPR-cBEST.gbk.
            - Find your genome of choice (or download an example genome below, S. coelicor A3).
            - Figure out what genes you want to target. For example, the actinorhodin cluster (SCO5087).

            ## **Instructions**
            Upload your genome file and CRISPR plasmid files, then click 'Submit' to generate your assembly.
            """, style=text_style),
            reference_content
        ], style={'padding': '20px', 'backgroundColor': '#2C3E50'}),
    ], className="mb-4"),

    html.Img(src='/assets/w2_pic.png', 
             style={'width': '50%', 'margin': '20px auto'}),

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
                dbc.Button("Download Example CRISPR Plasmid File", color="primary", className="mb-4"),
                href="/assets/ pCRISPR-cBEST.gbk",
                download="pCRISPR-cBEST.gbk"
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
                    html.H5("Genome File", className="card-title", style=text_style),
                    dcc.Upload(
                        id={'type': 'upload-component', 'index': 'genome-file-2'},  # Updated to use pattern matching ID
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Genome File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id={'type': 'filename-display', 'index': 'genome-file-2'}, children=[], style=text_style),  # Updated to use pattern matching ID
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.H4("2) Upload the plasmid of choice", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("CRISPR-BEST plasmid", className="card-title", style=text_style),
                    dcc.Upload(
                        id={'type': 'upload-component', 'index': 'single-vector-2'},  # Updated to use pattern matching ID
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select CRISPR-BEST plasmid File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id={'type': 'filename-display', 'index': 'single-vector-2'}, children=[], style=text_style),  # Updated to use pattern matching ID
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-5"),

        
    dbc.Row([
        dbc.Col([
            html.H4("3) Choose genes/regions to knock out", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Example for genes: SCO5087, SCO5087,... (comma-separated)", className="card-title", style=text_style),
                    html.H5("Example for regions: 1000-2000, 100000-101000,... (comma-separated)", className="card-title", style=text_style),

                    dbc.Textarea(
                        id='genes-to-KO_2',
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
            html.P("Per default the overhangs work with pCRISPR-cBEST.gbk", className="lead", style=text_style)
        ], width=12),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.Label('5 prime Overhang:', style=text_style),
            dbc.Input(
                id='forward-overhang-input_2',
                type='text',
                placeholder='Enter Forward Overhang',
                value='CGGTTGGTAGGATCGACGGC'
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
    ]),

    dbc.Row([
        dbc.Col([
            html.Label('3 prime Overhang:', style=text_style),
            dbc.Input(
                id='reverse-overhang-input_2',
                type='text',
                placeholder='Enter Reverse Overhang',
                value='GTTTTAGAGCTAGAAATAGC'
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
    ], className="mb-5"),

    dbc.Col([
        html.H4("5) Filtering metrics for sgRNAs", style=text_style),
        dbc.Row([
            dbc.Col([
                dbc.Label("GC Content Upper Bound", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='gc-upper_2',
                    type='number',
                    value=0.99,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Label("GC Content Lower Bound", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='gc-lower_2',
                    type='number',
                    value=0.01,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Label("Off-Target Seed Length", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='off-target-seed_2',
                    type='number',
                    value=13,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Label("Off-Target Upper Bound", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='off-target-upper_2',
                    type='number',
                    value=10,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Label("Cas Type", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dcc.Dropdown(
                    id='cas-type_2',
                    options=[
                        {'label': 'Cas9', 'value': 'cas9'},
                    ],
                    value='cas9',
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Label("Number of sgRNAs per region/locus tag", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='number-of-sgRNAs-per-group_2',
                    type='number',
                    value=5,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Checklist(
                    options=[
                        {"label": "Only Stop Codons", "value": 1},
                    ],
                    value=[1] if True else [],  # default value
                    id="only-stop-codons-checkbox_2",
                    inline=True,
                    switch=True,
                    className="big-switch"
                ),
            ], width=12, className="mb-3"),
        ]),
        dbc.Row([
            dbc.Col([
                dbc.Checklist(
                    options=[
                        {"label": "Editing Sequence Context", "value": 1},
                    ],
                    value=[1] if True else [],  # default value
                    id="editing_context_2",
                    inline=True,
                    switch=True,
                    className="big-switch"
                ),
            ], width=12, className="mb-3"),
        ]),
    ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),



    dbc.Row([
        dbc.Col([
            html.H4("7) Show advanced settings for checking primers", style=text_style),
            dbc.Checklist(
                options=[
                    {"label": "", "value": 1},
                ],
                value=[],
                id="show-advanced-settings-checkbox",
                inline=True,
                switch=True,
                className="big-switch"
            ),
        ], width=6),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.Div(id='advanced-settings-container', children=[
                dbc.Row([
                    dbc.Col([
                        dbc.Label("Choose Polymerase", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dcc.Dropdown(
                            id='chosen-polymerase_2',
                            options=dropdown_options,
                            value=polymerase_dict['Q5 High-Fidelity 2X Master Mix'],  # Set default value
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Target Melting Temperature (°C)", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dbc.Input(
                            id='melting-temperature_2',
                            type='number',
                            value=65,
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Primer Concentration (μM)", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dbc.Input(
                            id='primer-concentration_2',
                            type='number',
                            value=0.4,
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Primer Number Increment", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dbc.Input(
                            id='primer-number-increment_2',
                            type='number',
                            value=1,
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Flanking Region Number", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dbc.Input(
                            id='flanking-region-number_2',
                            type='number',
                            value=500,
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                ])
            ], style={"display": "none"})  # Hidden by default
        ], width=6),
    ], className="mb-4"),


    dbc.Row([
        dbc.Col([
            dbc.Button('Submit', id='submit-settings-button_2', color="primary", className="mt-3"),
        ], width=12),
    ], className="mb-4"),
    
    # Placeholder for the output
    dbc.Row([
        dbc.Col([
            dcc.Loading(
                id="loading-overlay-2",
                type="circle",
                children=[
                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Primers", className="card-title", style=text_style),
                            DataTable(id='primers-output-table_2', **table_style),
                        ])
                    ], style=card_style),

                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Filtered sgRNA DataFrame", className="card-title", style=text_style),
                            DataTable(id='filtered-df-table', **table_style),
                        ])
                    ], style=card_style),

                    dbc.Card([
                        dbc.CardBody([
                            html.H5("PCR", className="card-title", style=text_style),
                            DataTable(id='pcr-table_2', **table_style),
                        ])
                    ], style=card_style),

                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Overview of plasmids generated", className="card-title", style=text_style),
                            DataTable(id='plasmid-metadata-table_2', **table_style),
                        ])
                    ], style=card_style),

                    dbc.Card([
                        dbc.CardBody([
                            html.H5("Download folder with all data & protocols", className="card-title", style=text_style),
                            DataTable(id='all_data_2', **table_style),
                            html.A(
                                'Download All Data & protocols',
                                id='download-data-and-protocols-link_2',
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

])

