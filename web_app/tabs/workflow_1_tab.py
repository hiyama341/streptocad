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
], style={'fontSize': '1.5rem'})  

# Workflow 1 Tab content
workflow_1_tab = dcc.Tab(label="Workflow 1: Overexpression library construction", children=[
    dbc.Row([
        dbc.Col([
            html.P("What is it? ", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("A method to overexpress specific genes within a host organism using plasmid vectors.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Why use it?", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("Useful for studying the function of genes by increasing their expression levels.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
                html.Li("Enables the production of proteins at higher levels for research and industrial applications.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Getting Started:", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.Ul([
                html.Li("Select your plasmid vector that is suitable for the host organism.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
                html.Li("Prepare the gene sequences you want to overexpress.", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            ], style={'color': '#ddd', 'fontSize': '1.5rem'}),
            html.P("Upload your gene sequences and plasmid files, then configure the settings to generate your constructs.", className="lead", style={'color': '#ddd', 'fontSize': '1.5rem'}),
            reference_content
        ], style={'padding': '20px', 'backgroundColor': '#2C3E50'}),
    ], className="mb-4"),
    
    html.Img(src='/assets/workflow_1_pic.webp', 
             style={'width': '60%', 'margin': '20px auto'}),

    dbc.Row([
        dbc.Col(
            html.A(
                dbc.Button("Download Example Sequence File", color="primary"),
                href="/assets/GOE_regulators.gb", # GOE regulators
                download="example_sequence_file(LuxR_&_SARPs).csv"
            ),
            width={"size": 3, "order": 1}
        ),
        dbc.Col(
            html.A(
                dbc.Button("Download Example Plasmid File", color="primary", className="mb-4"),
                href="/assets/pOEX-PkasO.gb", # PKasO file
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
                    upload_component({'type': 'upload-input', 'index': 'sequences'}, text_style, link_style, upload_button_style),
                    html.Div(id={'type': 'uploaded-filename', 'index': 'sequences'}, children=[], style=text_style),
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
                    upload_component({'type': 'upload-input', 'index': 'plasmid'}, text_style, link_style, upload_button_style),
                    html.Div(id={'type': 'uploaded-filename', 'index': 'plasmid'}, children=[], style=text_style),
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
            value=polymerase_dict['Phusion High-Fidelity DNA Polymerase (GC Buffer)'],  # Set default value
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
            dbc.Card([
                dbc.CardBody([
                    html.H5("Primers", className="card-title", style=text_style),
                    DataTable(id='primers-output-table_1', **table_style),
                    # Removed other download buttons
                ])
            ], style=card_style),
            
            # Add Analyzed Primers DataTable
            dbc.Card([
                dbc.CardBody([
                    html.H5("Analyzed Primers", className="card-title", style=text_style),
                    DataTable(id='analyzed-primers-table_1', **table_style),
                    # Removed other download buttons
                ])
            ], style=card_style),

            dbc.Card([
                dbc.CardBody([
                    html.H5("PCR", className="card-title", style=text_style),
                    DataTable(id='pcr-table_1', **table_style),
                    # Removed other download buttons
                ])
            ], style=card_style),

            dbc.Card([
                dbc.CardBody([
                    html.H5("GenBank File", className="card-title", style=text_style),
                    # Removed other download buttons
                ])
            ], style=card_style),

            dbc.Card([
                dbc.CardBody([
                    html.H5("Download folder with all data & protocols", className="card-title", style=text_style),
                    DataTable(id='all_data_1', **table_style),
                    html.A(
                        'Data & protocols',
                        id='download-data-and-protocols-link_1',
                        download="data_package",
                        href="",
                        target="_blank",
                        className="btn btn-primary"
                    )
                ])
            ], style=card_style),

        ], width=10),
    ]),
])
