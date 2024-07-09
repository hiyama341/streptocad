from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable

from styling import text_style, upload_button_style, card_style, link_style, table_style
from streptocad.utils import polymerase_dict

# Reference content used in multiple tabs
reference_content = html.Div([
    html.P("Note: For more information on CRISPR techniques and details, please visit the Nature protocols article below (the figure and protocol is from there)."),
    html.A("CRISPR–Cas9, CRISPRi and CRISPR-BEST-mediated genetic manipulation in streptomycetes", href="https://www.nature.com/articles/s41596-020-0339-z", target="_blank"),
], style={'fontSize': '1.5rem'})

# Dropdown options for polymerases
dropdown_options = [{'label': key, 'value': value} for key, value in polymerase_dict.items()]

# Tab content for Multiple sgRNA-integration
golden_gate_tab = dcc.Tab(label="Multiple sgRNA-integration", children=[
    dbc.Row([
        dbc.Col([
            html.P("What is it? ", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("A method that employs multiple single guide RNAs (sgRNAs) to target several genomic locations simultaneously.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Why use it?", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("When you're looking to study or modify multiple genes or pathways in one go, this is your go-to approach.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Getting Started:", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("Find your plasmid of choice. We recommend that you use pCRISPR-MCBE_Csy4_kasopGFP.gb", style={'color': '#ddd', 'fontSize': '1.5rem'}),
                html.Li("Figure out what genes you want to target. For example the actinorhodin cluster (SCO5087)", style={'color': '#ddd', 'fontSize': '1.5rem'})
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Upload your genome file and CRISPR plasmid files, then click 'Submit' to generate your assembly.", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            reference_content
        ], style={'padding': '20px', 'backgroundColor': '#2C3E50'}),
    ], className="mb-4"),

    html.Img(src='/assets/workflow_3_pic.webp', 
             style={'width': '60%', 'margin': '20px auto'}),

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
                href="/assets/pCRISPR-MCBE_Csy4_kasopGFP.gb",
                download="pCRISPR-MCBE_Csy4_kasopGFP.gb"
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
            html.H4("4) Enter sgRNA handle", className="mb-3", style=text_style),
            html.P('As a default we use the sgRNA handle with cys4 (see figure above)', style=text_style),
            dbc.Input(
                id='sgRNA-handle-input',
                type='text',
                placeholder='Enter sgRNA handle: e.g., 110-nt (pJET1.2–sgRNA handle)',
                value='GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGTTCACTGCCGTATAGGCAGCTAAGAAA'
            ),
        ], width=6),
    ], className="mb-5"),

    dbc.Row([
        dbc.Col([
            html.H4("5) Enter Target Melting Temperature", className="mb-3", style=text_style),
            dbc.Input(id='input-tm', type='number', placeholder='Enter Melting Temperature', value=55),
            html.Br(),
        ], width=6),
    ], className="mb-5"),

    dbc.Row([
    dbc.Col([
        html.H4("6) Filtering metrics for sgRNAs", style=text_style),
        dbc.Row([
            dbc.Col([
                dbc.Label("GC Content Upper Bound", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='gc-upper',
                    type='number',
                    value=0.99,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("GC Content Lower Bound", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='gc-lower',
                    type='number',
                    value=0.01,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("Off-Target Seed Length", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='off-target-seed',
                    type='number',
                    value=13,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("Off-Target Upper Bound", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='off-target-upper',
                    type='number',
                    value=10,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("Cas Type", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dcc.Dropdown(
                    id='cas-type',
                    options=[
                        {'label': 'Cas9', 'value': 'cas9'},
                    ],
                    value='cas9',
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Label("Number of sgRNAs per Group", style={'color': '#ddd'}),  # Explicitly set color to ensure visibility
                dbc.Input(
                    id='number-of-sgRNAs-per-group',
                    type='number',
                    value=5,
                    style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
                ),
            ], width=12, className="mb-3"),
            dbc.Col([
                dbc.Checklist(
                    options=[
                        {"label": "Only Stop Codons", "value": 1},
                    ],
                    value=[1] if True else [],  # default value
                    id="only-stop-codons-checkbox",
                    inline=True,
                    switch=True,
                    className="big-switch"
                ),
            ], width=12, className="mb-3"),
        ]),
    ], width=6),
], className="mb-4"),


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
    ], className="mb-4"),

    dbc.Row([
        dbc.Col([
            dbc.Button('Submit', id='submit-button', color="primary", className="mt-3"),
        ], width=12),
    ], className="mb-5"),

    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardBody([
                    html.H5("Primers", className="card-title", style=text_style),
                    DataTable(id='primer-table', **table_style),
                    html.A(
                        'Download Primers CSV',
                        id='download-primers-link',
                        download="primers.csv",
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
                    html.H5("Overhang", className="card-title", style=text_style),
                    DataTable(id='overhang-table', **table_style),
                ])
            ], style=card_style),

            dbc.Card([
                dbc.CardBody([
                    html.H5("GenBank File", className="card-title", style=text_style),
                    html.A(
                        'Download GenBank File',
                        id='genbank-file',
                        download="plasmid.gb",
                        href="",
                        target="_blank",
                        className="btn btn-primary"
                    )
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

        ], width=12),
    ]),
])

