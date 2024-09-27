# workflow_1_tab.py
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable
from styling import text_style, upload_button_style, card_style, link_style, table_style, table_header_style, table_row_style
from streptocad.utils import polymerase_dict
from components import upload_component

# Dropdown options for polymerases
dropdown_options = [{'label': key, 'value': value} for key, value in polymerase_dict.items()]

# Reference content used in multiple tabs
reference_content = html.Div([
    html.P("Note: For more information on overexpression workflows and details, please visit our paper.", style=text_style),
    html.A("StreptoCAD: An open-source software toolbox supporting genome engineering workflows in Streptomyces", href="https://example.com/overexpression-protocol", target="_blank", style=link_style),
], style=text_style)  

# Workflow 1 Tab content
workflow_1_tab = dcc.Tab(label="Workflow 1: Overexpression library construction", children=[
    dbc.Row([
        dbc.Col([
            dcc.Markdown("""
            ## **What is it?**
            - A method to overexpress specific genes within a host organism using plasmid vectors.

            ## **Why use it?**
            - Useful for studying the function of genes by increasing their expression levels.
            - Enables the production of proteins at higher levels for research and industrial applications.

            ## **Getting Started**
            - Select your plasmid vector that is suitable for the host organism.
            - Prepare the gene sequences you want to overexpress.

            ## **Instructions**
            Upload your gene sequences and plasmid files, then configure the settings to generate your constructs.
            """, style=text_style),
            reference_content
        ], style={'padding': '20px', 'backgroundColor': '#2C3E50'}),
    ], className="mb-4"),
    
    html.Img(src='/assets/w1_pic.png', 
             style={'width': '50%', 'margin': '20px auto'}),

    dbc.Row([
        dbc.Col(
            html.A(
                dbc.Button("Download Example Sequence File", color="primary"),
                href="/assets/GOE_regulators.gb", 
                download="example_sequence_file(LuxR_&_SARPs).csv"
            ),
            width={"size": 3, "order": 1}
        ),
        dbc.Col(
            html.A(
                dbc.Button("Download Example Plasmid File", color="primary", className="mb-4"),
                href="/assets/pOEX-PkasO.gb", 
                download="pOEX-PkasO.gb"
            ),
            width={"size": 3, "order": 2}
        )
    ], className="mb-4", justify="start"),

    dbc.Row([
        dbc.Col([
            html.H4("1) Upload your gene sequences", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Sequences (GenBank file)", className="card-title", style=text_style),
                    dcc.Upload(
                        id={'type': 'upload-component', 'index': 'sequences'},  # Updated to use pattern matching ID
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Sequences File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id={'type': 'filename-display', 'index': 'sequences'}, children=[], style=text_style),  # Updated to use pattern matching ID
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.H4("2) Upload your plasmid", style=text_style),
            dbc.Card([
                dbc.CardBody([
                    html.H5("Plasmid File", className="card-title", style=text_style),
                    dcc.Upload(
                        id={'type': 'upload-component', 'index': 'plasmid'},  # Updated to use pattern matching ID
                        children=html.Div([
                            'Drag and Drop or ',
                            html.A('Select Plasmid File', style=link_style)
                        ], style=text_style),
                        style=upload_button_style,
                        multiple=False
                    ),
                    html.Div(id={'type': 'filename-display', 'index': 'plasmid'}, children=[], style=text_style),  # Updated to use pattern matching ID
                ])
            ], style=card_style),
        ], width=6),
    ], className="mb-5"),
    
    dbc.Row([
        dbc.Col([
            html.H4("3) Choose overlapping sequences", style=text_style),
            html.P("Please enter the 5' and 3' overhangs below for the oligo nucleotide to be made.", className="lead", style=text_style),
            html.P("Per default the overhangs work with pOEX-PkasO", className="lead", style=text_style)
        ], width=5),
    ], className="mb-3"),

    dbc.Row([
        dbc.Col([
            html.Label('5 prime Overhang:', style=text_style),
            dbc.Input(
                id='forward-overhang-input_1',
                type='text',
                placeholder='Enter Forward Overhang',
                value='GGCGAGCAACGGAGGTACGGACAGG',
                style={'color': '#000'}  # Ensure text is black
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
        
        dbc.Col([
            html.Label('3 prime Overhang:', style=text_style),
            dbc.Input(
                id='reverse-overhang-input_1',
                type='text',
                placeholder='Enter Reverse Overhang',
                value='CGCAAGCCGCCACTCGAACGGAAGG',
                style={'color': '#000'}  # Ensure text is black
            ),
        ], width=6, style={"marginRight": "10px", "marginLeft": "10px"}),
    ], className="mb-5"),  # Increase margin-bottom here to add more space

    dbc.Col([
        html.H4("4) Customizable Settings", style=text_style),
        dbc.Label("Choose Polymerase", style={'color': '#ddd'}),  
        dcc.Dropdown(
            id='chosen-polymerase_1',
            options=dropdown_options,
            value=polymerase_dict['Q5 High-Fidelity 2X Master Mix'],  # Set default value
            style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
        ),
        dbc.Label("Target Melting Temperature (°C)", style={'color': '#ddd'}),  
        dbc.Input(
            id='melting-temperature_1',
            type='number',
            value=65,
            style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
        ),
        dbc.Label("Primer Concentration (μM)", style={'color': '#ddd'}),  
        dbc.Input(
            id='primer-concentration_1',
            type='number',
            value=0.4,
            style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
        ),
        dbc.Label("Primer Number Increment", style={'color': '#ddd'}),  
        dbc.Input(
            id='primer-number-increment_1',
            type='number',
            value=1,
            style={'color': '#000', 'width': '100%'}  # Ensure text is black and set width
        ),
    ], width=6, className="mb-3"),

    dbc.Row([
        dbc.Col([
            dbc.Button('Submit', id='submit-settings-button_1', color="primary", className="mt-3"),
        ], width=12),
    ], className="mb-4"),
    
# Placeholder for the output
dbc.Row([
    dbc.Col([
        dcc.Loading(
            id="loading-overlay",
            type="circle",
            children=[
                dbc.Card([
                    dbc.CardBody([
                        html.H5("Primers", className="card-title", style=text_style),
                        DataTable(id='primers-output-table_1', **table_style),
                    ])
                ], style=card_style),
                
                # Add Analyzed Primers DataTable
                dbc.Card([
                    dbc.CardBody([
                        html.H5("Analyzed Primers", className="card-title", style=text_style),
                        DataTable(id='analyzed-primers-table_1', **table_style),
                    ])
                ], style=card_style),

                dbc.Card([
                    dbc.CardBody([
                        html.H5("PCR", className="card-title", style=text_style),
                        DataTable(id='pcr-table_1', **table_style),
                    ])
                ], style=card_style),

                dbc.Card([
                    dbc.CardBody([
                        html.H5("Overview of plasmids generated", className="card-title", style=text_style),
                    ])
                ], style=card_style),

                dbc.Card([
                    dbc.CardBody([
                        html.H5("Download folder with all data & protocols", className="card-title", style=text_style),
                        DataTable(id='all_data_1', **table_style),
                        html.A(
                            'Download All Data & protocols',
                            id='download-data-and-protocols-link_1',
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
