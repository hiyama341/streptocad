# workflow_1_tab.py
from dash import dcc, html
import dash_bootstrap_components as dbc
from dash.dash_table import DataTable
from styling import text_style, upload_button_style, card_style, link_style, table_style, table_header_style, table_row_style
from streptocad.utils import polymerase_dict
from components import upload_component
from dash import dcc, html
import dash_bootstrap_components as dbc

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
            ## **What is overexpression library construction?**
            - A method to overexpress specific genes within a host organism using plasmid vectors.

            ## **Why use it?**
            - Useful for studying the function of genes by increasing their expression levels.
            - Enables the production of proteins at higher levels for research and industrial applications.

            ## **Getting Started**
            - Find your plasmid of choice. We recommend that you use pOEx_PkasO.gbk.
            - Prepare the gene sequences you want to overexpress into one GenBank file.

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
                    html.H5("Sequences (GenBank format)", className="card-title", style=text_style),
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
                    html.H5("Plasmid File (GenBank format)", className="card-title", style=text_style),
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
            # Main heading with consolidated tooltip
            html.H4([
                "3) Choose overlapping sequences ",
                html.Span("ⓘ", id="overlapping-sequences-info-icon", style=link_style)
            ], style=text_style),
            
            # Consolidated tooltip with all the details
            dbc.Tooltip(
                """
                - Please enter the 5' and 3' overhangs below for the oligo nucleotide to be made.

                - These overhangs are added as 5' and 3' overhangs to the primers being generated.

                - By default, these overhangs are compatible with the pOEX-PkasO plasmid.
                """,
                target="overlapping-sequences-info-icon",
                placement="top",
            ),
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
        
        # Choose Polymerase with tooltip
        dbc.Label([
            "Choose Polymerase ",
            html.Span("ⓘ", id="polymerase-tooltip-icon", style=link_style)
        ], style={'color': '#ddd'}),
        dbc.Tooltip(
            "Select the polymerase type. The default is optimized for high fidelity.",
            target="polymerase-tooltip-icon",
            placement="top",
        ),
        dcc.Dropdown(
            id='chosen-polymerase_1',
            options=dropdown_options,
            value=polymerase_dict['Q5 High-Fidelity 2X Master Mix'],
            style={'color': '#000', 'width': '100%'}
        ),
        
        # Target Melting Temperature with tooltip
        dbc.Label([
            "Target Melting Temperature (°C) ",
            html.Span("ⓘ", id="melting-temp-tooltip-icon", style=link_style)
        ], style={'color': '#ddd'}),
        dbc.Tooltip(
            "Set the desired melting temperature for the PCR reactions.",
            target="melting-temp-tooltip-icon",
            placement="top",
        ),
        dbc.Input(
            id='melting-temperature_1',
            type='number',
            value=65,
            style={'color': '#000', 'width': '100%'}
        ),
        
        # Primer Concentration with tooltip
        dbc.Label([
            "Primer Concentration (μM) ",
            html.Span("ⓘ", id="primer-concentration-tooltip-icon", style=link_style)
        ], style={'color': '#ddd'}),
        dbc.Tooltip(
            "Specify the concentration of primers used in the PCR reaction.",
            target="primer-concentration-tooltip-icon",
            placement="top",
        ),
        dbc.Input(
            id='primer-concentration_1',
            type='number',
            value=0.4,
            style={'color': '#000', 'width': '100%'}
        ),
        
        # restriction enzymes with tooltip
        dbc.Label([
            "Choose restriction enzymes for the plasmid digestion",
            html.Span("ⓘ", id="restriction-enzymes-tooltip-icon", style=link_style)
        ], style={'color': '#ddd'}),
        dbc.Tooltip(
            """Choose the restriction enzyme(s) for the plasmid digestion.
            This is just for your convenience. 
            By default, we use StuI for the digestion of pOEX-PkasO.
            If you would like to use a different enzyme, please enter the name of the enzyme(s) in the following format:
            enzyme1,enzyme2,enzyme3,...
            Tip: Remember to update overlapping sequences for the assembly if you change the enzymes.
            """,
            target="restriction-enzymes-tooltip-icon",
            placement="top",
        ),
        dbc.Input(
            id='restriction-enzymes_1',
            type='text',
            value='StuI',
            style={'color': '#000', 'width': '100%'}
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
                        DataTable(id='plasmid-metadata-table_1', **table_style),
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
