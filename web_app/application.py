#!/usr/bin/env python
# MIT License

# DASH imports
import dash
from dash import dcc, html, dash_table, exceptions
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
from dash.exceptions import PreventUpdate
from urllib.parse import quote

# Styling module
from styling import tab_style, active_tab_style, vertical_tab_style, text_style, upload_button_style, card_style, link_style, button_style_darkly, table_style, table_header_style, table_row_style

# Import callback registration functions
from callbacks.workflow_1 import register_workflow_1_callbacks  # Import the new callback registration function
from callbacks.workflow_2 import register_workflow_2_callbacks
from callbacks.workflow_3 import register_workflow_3_callbacks
from callbacks.workflow_4 import register_workflow_4_callbacks  # Added callback registration for workflow 4
from callbacks.workflow_5 import register_workflow_5_callbacks
from callbacks_interactivity import register_interactivity_callbacks  # Import the interactivity callbacks

# Import tab content
from tabs.workflow_1_tab import workflow_1_tab
from tabs.workflow_2_tab import crispr_cb_tab
from tabs.workflow_3_tab import golden_gate_tab
from tabs.workflow_4_tab import crispri_tab  # Corrected to crispri_tab
from tabs.workflow_5_tab import gibson_tab
from tabs.workflow_6_tab import cas3_tab


external_stylesheets = [
    dbc.themes.DARKLY,
    'https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css',
    '/assets/style.css'  # Include custom CSS for styling switches
]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

# Header content
header_content = dbc.Row(
    [
        # Column for the Header Text
        dbc.Col(
            html.H2(
                "StreptoCAD: An open-source software toolbox supporting genome engineering workflows in Streptomyces - version 0.1.0", 
                style={'color': '#ddd', 'fontSize': '2.5rem'}  # Adjust font size
            ),
            width=12,  # Use the full width of the row
            className="d-flex align-items-center",  # Left-align the text horizontally and center it vertically
        ),
    ],
    style={'marginBottom': '10px', 'borderBottom': '1px solid #444', 'backgroundColor': '#1A242F', 'padding': '10px 0'}  # Adjust padding and margin
)

# Introduction page content
introduction_page = html.Div(
    [
        html.H1("Welcome to StreptoCAD", style={'textAlign': 'center', 'marginTop': '70px', 'color': '#ddd'}),
        html.P(
            "StreptoCAD is a toolbox supporting genome engineering workflows in Streptomyces.",
            style={'textAlign': 'center', 'fontSize': '30px', 'marginBottom': '10px', 'color': '#ddd'}
        ),

        html.P(
            "If you have used StreptoCAD in a publication, please remember to cite us.",
            style={'textAlign': 'center', 'fontSize': '20px', 'marginBottom': '50px', 'color': '#ddd'}
        ),
        dbc.Button("Let's get started", id='start-button', color="primary", size="lg", className="d-block mx-auto"),
    ],
    style={'height': '100vh', 'display': 'flex', 'flexDirection': 'column', 'justifyContent': 'center', 'alignItems': 'center', 'backgroundColor': '#2C3E50'}
)
# Footer content
footer = html.Footer([
    html.Div([
        html.A("Contact Us", href="mailto:your-email@example.com", style={"marginRight": "20px", "fontSize": "30px"}),  # Email link
        html.A([html.I(className="fab fa-twitter")], href="https://twitter.com/YourTwitterHandle", style={"marginRight": "10px", "fontSize": "30px"}),  # Twitter link
        html.A([html.I(className="fab fa-github")], href="https://github.com/YourGitHubProfile", style={"fontSize": "30px"}),  # GitHub link
    ], style={"textAlign": "center", "padding": "40px 0", "backgroundColor": "#1A242F", "color": "#ddd", "fontSize": "30px"})  # Center-align text, add padding, set background and text colors, and increase font size
])

# Main VERTICAL app layout
main_layout = dcc.Loading(
    id="loading",
    type="cube",
    children=[
        dbc.Container([
            header_content,
            dbc.Row([
                # Main Tabs on the left side
                dbc.Col([
                    html.H2("Workflows", style={'color': '#ddd', 'padding': '20px 0'}),  # Add a heading
                    dcc.Tabs(id="main-tabs", vertical=True, style={'borderRight': '0.1px solid black'}, children=[
                        dcc.Tab(label="Workflow 1: Overexpression library construction", value="workflow_1", style=vertical_tab_style, selected_style=active_tab_style),
                        dcc.Tab(label="Workflow 2: Single CRISPR-BEST plasmid generation", value="workflow_2", style=vertical_tab_style, selected_style=active_tab_style),
                        dcc.Tab(label="Workflow 3: Multiplexed CRISPR-BEST plasmid generation", value="workflow_3", style=vertical_tab_style, selected_style=active_tab_style),
                        dcc.Tab(label="Workflow 4: CRISPRi plasmid generation", value="workflow_4", style=vertical_tab_style, selected_style=active_tab_style),
                        dcc.Tab(label="Workflow 5: In-frame deletion with CRISPR-Cas9", value="workflow_5", style=vertical_tab_style, selected_style=active_tab_style),
                        dcc.Tab(label="Workflow 6: In-frame deletion with CRISPR-Cas3", value="workflow_6", style=vertical_tab_style, selected_style=active_tab_style),
                    ])
                ], width=3),
                dbc.Col([
                    dcc.Loading(  # Add a loading spinner specifically for the tab content
                        id="loading-tab-content",
                        type="circle",
                        children=[
                            html.Div(id="tab-content", style={'padding': '20px'})  # Content will be loaded here based on the selected tab
                        ]
                    )
                ], width=9),
            ]),
            footer,
        ], fluid=True, style={'backgroundColor': '#2C3E50', 'padding': '20px', 'color': '#ddd'})
    ]
)

# Overall layout with introduction page
app.layout = html.Div(id="page-content", children=[introduction_page])

# Callbacks
@app.callback(
    Output("page-content", "children"),
    [Input("start-button", "n_clicks")],
    [State("page-content", "children")]
)
def display_main_page(n_clicks, children):
    if n_clicks:
        return main_layout
    return children

@app.callback(
    Output("tab-content", "children"),
    [Input("main-tabs", "value")]
)
def render_tab_content(tab):
    if tab == "workflow_1":
        return workflow_1_tab
    elif tab == "workflow_2":
        return crispr_cb_tab
    elif tab == "workflow_3":
        return golden_gate_tab
    elif tab == "workflow_4":
        return crispri_tab
    elif tab == "workflow_5":
        return gibson_tab
    elif tab == "workflow_6":
        return cas3_tab
    else:
        welcome_message = """

### Welcome to StreptoCAD! 🧬

StreptoCAD is an open-source software toolbox designed to help you **build biology easier**. 

It supports genome engineering workflows in Streptomyces, providing a user-friendly interface to design your experiments.

👉 To get started, please select a workflow from the tabs on the left. 

We hope you find StreptoCAD useful in your research endeavors. If you have any questions or feedback, please don't hesitate to reach out to our team.

**Happy bioengineering!** 🔬


"""
    return html.Div([
        dcc.Markdown(
            welcome_message,
            style={
                'fontSize': '1.5rem',
                'lineHeight': '1.5',
                'textAlign': 'left',
                'padding': '20px',
                'borderRadius': '10px',
                'boxShadow': '0 4px 6px rgba(0, 0, 0, 0.1)'
            }
        )
    ])
# Callback to show/hide advanced settings
@app.callback(
    Output('advanced-settings-container', 'style'),
    [Input('show-advanced-settings-checkbox', 'value')]
)
def toggle_advanced_settings(checkbox_value):
    if checkbox_value:
        return {'display': 'block'}
    else:
        return {'display': 'none'}

# Register callbacks to the app
register_workflow_1_callbacks(app)
register_workflow_2_callbacks(app)
register_workflow_3_callbacks(app)
register_workflow_4_callbacks(app)
register_workflow_5_callbacks(app)
#register_workflow_6_callbacks(app)
register_interactivity_callbacks(app)  # Register the interactivity callbacks



if __name__ == '__main__':
    app.run_server(debug=True)
