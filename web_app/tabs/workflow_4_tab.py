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

# Tab content for CRISPRi
crispri_tab = html.Div(children=[
    dbc.Row([
        dbc.Col([
            html.P("What is it? ", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("CRISPR interference (CRISPRi) is a method that uses a catalytically dead Cas9 (dCas9) to block transcription of target genes.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Why use it?", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("It allows for the precise and reversible silencing of genes without altering the DNA sequence.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Getting Started:", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("Find your plasmid of choice. We recommend that you use a dCas9 vector suitable for CRISPRi.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
                html.Li("Figure out what genes you want to target. For example the actinorhodin cluster (SCO5087)", style={'color': '#ddd', 'fontSize': '1.5rem'})
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Upload your genome file and CRISPRi plasmid files, then click 'Submit' to generate your assembly.", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            reference_content
        ], style={'padding': '20px', 'backgroundColor': '#2C3E50'}),
    ], className="mb-4"),

    html.Img(src='/assets/workflow_4_pic.webp', 
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
                dbc.Button("Download Example CRISPRi Vector File", color="primary", className="mb-4"),
                href="/assets/pCRISPR-dCas9.gbk",
                download="pCRISPR-dCas9.gbk"
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
                        id='upload-genome-file_4',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Genome File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id='uploaded-genome-filename_4', children=[], style=text_style),
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.H4("2) Upload the plasmid of choice", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("CRISPRi Vector", className="card-title", style=text_style),
                    dcc.Upload(
                        id='upload-crispri-vector_4',
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select CRISPRi Vector File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id='uploaded-crispri-vector-filename_4', children=[], style=text_style),
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
            html.Label('5 prime Overhang:', style=text_style),
            dbc.Input(
                id='forward-overhang-input_4',
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
                id='reverse-overhang-input_4',
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
                        id='gc-upper_4',
                        type='number',
                        value=0.99,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label("GC Content Lower Bound", style={'color': '#ddd'}),  
                    dbc.Input(
                        id='gc-lower_4',
                        type='number',
                        value=0.01,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label("Off-Target Seed Length", style={'color': '#ddd'}),  
                    dbc.Input(
                        id='off-target-seed_4',
                        type='number',
                        value=13,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label("Off-Target Upper Bound", style={'color': '#ddd'}),  
                    dbc.Input(
                        id='off-target-upper_4',
                        type='number',
                        value=10,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label("Cas Type", style={'color': '#ddd'}),  
                    dcc.Dropdown(
                        id='cas-type_4',
                        options=[
                            {'label': 'Cas9', 'value': 'cas9'},
                        ],
                        value='cas9',
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label("Number of sgRNAs per Group", style={'color': '#ddd'}),  
                    dbc.Input(
                        id='number-of-sgRNAs-per-group_4',
                        type='number',
                        value=5,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
                dbc.Col([
                    dbc.Label("Extension to Promoter Region", style={'color': '#ddd'}),  
                    dbc.Input(
                        id='extension-to-promoter-region_4',
                        type='number',
                        value=100,
                        style={'color': '#000', 'width': '100%'}
                    ),
                ], width=12, className="mb-3"),
            ]),
        ], width=6, className="mb-4"),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            dbc.Button('Submit', id='submit-settings-button_4', color="primary", className="mt-3"),
        ], width=12),
    ], className="mb-4"),
    
    ### Output
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.Card([

                dbc.CardBody([
                    html.H5("Filtered sgRNAs", className="card-title", style=text_style),
                    DataTable(id='mutated-sgrna-table_4', **table_style)
                ])
            ], style=card_style),


                dbc.CardBody([
                    html.H5("Primers", className="card-title", style=text_style),
                    DataTable(id='primer-table_4', **table_style),
                    html.A(
                        'Download Primers CSV',
                        id='download-primers-link_4',
                        download="primers.csv",
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
                        id='genbank-file-single_4',
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
                    DataTable(id='all_data_4', **table_style),
                    html.A(
                        'Data & protocols',
                        id='download-data-and-protocols-link_4',
                        download="all_data.csv",
                        href="",
                        target="_blank",
                        className="btn btn-primary"
                    )
                ])
            ], style=card_style),
            
        ], width=10),
    ]),
])