from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable

from styling import text_style, upload_button_style, card_style, link_style, table_style, table_header_style, table_row_style
from streptocad.utils import polymerase_dict
# Reference content used in multiple tabs
reference_content = html.Div([
    html.P("Note: For more information on CRISPR techniques and details, please visit the article  below (the figure and protocol is from there).", style={'fontSize': '1.5rem'}),
    html.A("CASCADE-Cas3 Enables Highly Efficient Genome Engineering in Streptomyces Species", href="https://www.biorxiv.org/content/10.1101/2023.05.09.539971v1.full", target="_blank", style={'fontSize': '1.5rem'}),
], style={'fontSize': '1.5rem'})

# Dropdown options for polymerases
dropdown_options = [{'label': key, 'value': value} for key, value in polymerase_dict.items()]

# Tab content for CRISPR–Cas9 plasmid construction for in-frame deletion
cas3_tab = dcc.Tab(label="CRISPR–Cas3 plasmid construction", children=[
    dbc.Row([
        dbc.Col([
            # Introductory text
            html.P("What is it? ", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([html.Li("A method that retrieves repair regions up/downstream and generates primers and PCRs for the assembly.", style={'color': '#ddd', 'fontSize': '1.5rem'})]),
            html.P("Why use it? ", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([html.Li("When you aim to delete a specific gene.", style={'color': '#ddd', 'fontSize': '1.5rem'})]),
            html.P("Getting Started:", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("Find your Cas3 plasmid (Or download an example genome below).", style={'color': '#ddd', 'fontSize': '1.5rem'}),
                html.Li("Fetch your organisms genome (Or download an example genome below).", style={'color': '#ddd', 'fontSize': '1.5rem'})
            ]),
            html.P("Upload your files, then click 'Submit' to generate your assembly.", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            reference_content
        ], style={'padding': '20px', 'backgroundColor': '#2C3E50'}),
    ], className="mb-4"),
    
    html.Img(src='/assets/workflow_6_pic.webp', 
             style={'width': '60%', 'margin': '20px auto'}),  # Adjust the width as needed
    
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
                dbc.Button("Download Example CRISPR Vector File", color="primary", className="mb-4"),
                href="/assets/pCRISPR_cas3.gbk",
                download="pCRISPR_cas3.gbk"
            ),
            width={"size": 3, "order": 2}
        )
    ], className="mb-4", justify="start"),

    dbc.Row([
        dbc.Col([
            html.H4("1) Upload your genome file", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Genome File", className="card-title", style=text_style),
                    dcc.Upload(
                        id='upload-genome-file',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Genome File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id='uploaded-genome-filename', children=[], style=text_style),
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.H4("2) Upload the plasmid of choice", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("CRISPR Vector", className="card-title", style=text_style),
                    dcc.Upload(
                        id='upload-single-vector',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select CRISPR Vector File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id='uploaded-single-vector-filename', children=[], style=text_style),
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-5"),
        
    dbc.Row([
        dbc.Col([
            html.H4("3) Choose genes to knock out", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Genes to Knock Out (comma-separated)", className="card-title", style=text_style),
                    dbc.Textarea(
                        id='genes-to-KO',
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
            html.P("Per default the overhangs work with pCRISPR–Cas9_plasmid_addgene.gbk", className="lead", style=text_style)
        ], width=12),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.Label('5 prime Overhang:', style=text_style),
            dbc.Input(
                id='forward-overhang-input',
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
                id='reverse-overhang-input',
                type='text',
                placeholder='Enter Reverse Overhang',
                value='GTTTTAGAGCTAGAAATAGC'
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
    ], className="mb-5"),

    dbc.Row([
    dbc.Col([
        html.H4("5) Filtering metrics for sgRNAs", style=text_style),
        dbc.Row([
            dbc.Col([
                dbc.Label("GC Content Upper Bound", style={'color': '#ddd'}),  
                dbc.Input(
                    id='gc-upper',
                    type='number',
                    value=0.99,
                    style={'color': '#000', 'width': '100%'}
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("GC Content Lower Bound", style={'color': '#ddd'}),  
                dbc.Input(
                    id='gc-lower',
                    type='number',
                    value=0.01,
                    style={'color': '#000', 'width': '100%'}
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("Off-Target Seed Length", style={'color': '#ddd'}),  
                dbc.Input(
                    id='off-target-seed',
                    type='number',
                    value=13,
                    style={'color': '#000', 'width': '100%'}
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("Off-Target Upper Bound", style={'color': '#ddd'}),  
                dbc.Input(
                    id='off-target-upper',
                    type='number',
                    value=10,
                    style={'color': '#000', 'width': '100%'}
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("Cas Type", style={'color': '#ddd'}),  
                dcc.Dropdown(
                    id='cas-type',
                    options=[
                        {'label': 'Cas9', 'value': 'cas9'},
                        {'label': 'Cas12a', 'value': 'cas12a'}
                    ],
                    value='cas9',
                    style={'color': '#000', 'width': '100%'}
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("Number of sgRNAs per Group", style={'color': '#ddd'}),  
                dbc.Input(
                    id='number-of-sgRNAs-per-group',
                    type='number',
                    value=5,
                    style={'color': '#000', 'width': '100%'}
                ),
            ], width=12, className="mb-3"),
        ]),
    ], width=6, className="mb-4"),
], className="mb-3"),


    dbc.Row([
        dbc.Col([
            html.H4("6) Show advanced settings for checking primers and repair templates", style=text_style),
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
                            id='chosen-polymerase',
                            options=dropdown_options,
                            value=polymerase_dict['Phusion High-Fidelity DNA Polymerase (GC Buffer)'],  # Set default value
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Target Melting Temperature (°C)", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dbc.Input(
                            id='melting-temperature',
                            type='number',
                            value=65,
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Primer Concentration (μM)", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dbc.Input(
                            id='primer-concentration',
                            type='number',
                            value=0.4,
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Primer Number Increment", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dbc.Input(
                            id='primer-number-increment',
                            type='number',
                            value=1,
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                    dbc.Col([
                        dbc.Label("Flanking Region Number", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                        dbc.Input(
                            id='flanking-region-number',
                            type='number',
                            value=500,
                            style={'color': '#000'}  # Ensure text is black
                        ),
                    ], width=6),
                ])
            ], style={"display": "none"})  # Hidden by default
        ], width=6),
    ], className="mb-3"),


    dbc.Row([
        dbc.Col([
            dbc.Button('Submit', id='submit-settings-button', color="primary", className="mt-3"),
        ], width=12),
    ], className="mb-4"),
    
    # Placeholder for the output
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H5("Primers", className="card-title", style=text_style),
                    DataTable(id='primers-output-table', **table_style),
                    html.A(
                        'Download CSV File',
                        id='csv_download_link',
                        download="ssDNA_bridging_oligos.csv",
                        href="",
                        target="_blank",
                        className="btn btn-primary"
                    )
                ])
            ], style=card_style),

            dbc.Card([
                dbc.CardBody([
                    html.H5("PCR", className="card-title", style=text_style),
                    DataTable(id='pcr-table', **table_style),
                    html.A(
                        'Download PCR CSV',
                        id='download-pcr-link',
                        download="pcr.csv",
                        href="",
                        target="_blank",
                        className="btn btn-primary"
                    )
                ])
            ], style=card_style),

            dbc.Card([
                dbc.CardBody([
                    html.H5("GenBank File", className="card-title", style=text_style),
                    html.A(
                        'Download GenBank File',
                        id='genbank-file-single',
                        download="plasmid-single.zip",
                        href="",
                        target="_blank",
                        className="btn btn-primary"
                    ), 
                ])
            ], style=card_style),


        dbc.Card([
            dbc.CardBody([
                html.H5("Download folder with all data & protocols", className="card-title", style=text_style),
                DataTable(id='all_data', **table_style),
                html.A(
                    'Data & protocols',
                    id='download-data-and-protocols-link',
                    download="all_data.csv",
                    href="",
                    target="_blank",
                    className="btn btn-primary"
                )
            ])
        ], style=card_style),
        
        ], width=6),
    ]),
])
