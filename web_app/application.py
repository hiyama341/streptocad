#!/usr/bin/env python
# MIT License

# DASH imports
import dash
from dash import dcc, html, no_update
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import time
from urllib.parse import quote

# Bootstrap components import
import dash_bootstrap_components as dbc

# Styling module
from styling import tab_style, active_tab_style, vertical_tab_style, text_style, upload_button_style, card_style, link_style, button_style_darkly, table_style

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

# Import external components
from header import header_content
from introduction_page import introduction_page
from welcome_message import welcome_message_content
from footer import footer_content

external_stylesheets = [
    dbc.themes.DARKLY,
    'https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css',
    '/assets/style.css',  # Include custom CSS for styling switches
    '/assets/custom_style.css'  # Include custom CSS for the popup
]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions=True)

# Main VERTICAL app layout with the loading spinner for submit button
main_layout = dbc.Container([
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
                dcc.Tab(label="Workflow 5: CRISPR-Cas9 plasmid generation", value="workflow_5", style=vertical_tab_style, selected_style=active_tab_style),
                dcc.Tab(label="Workflow 6: CRISPR-Cas3 plasmid generation", value="workflow_6", style=vertical_tab_style, selected_style=active_tab_style),
                dcc.Tab(label="About StreptoCAD", value="about", style=vertical_tab_style, selected_style=active_tab_style)
            ])
        ], width=3),
        dbc.Col([
            html.Div(id="tab-content", style={'padding': '20px'})
        ], width=9),
    ]),
    footer_content,

    #### ERROR BOX 
    dcc.ConfirmDialog(
        id='error-dialog_1',
        message=""
    ),
    dcc.ConfirmDialog(
        id='error-dialog_2',
        message=""
    ),
    dcc.ConfirmDialog(
        id='error-dialog_3',
        message=""
    ),
], fluid=True, style={'backgroundColor': '#2C3E50', 'padding': '20px', 'color': '#ddd'})

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
    if tab == "about":
        return welcome_message_content
    elif tab == "workflow_1":
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
        return welcome_message_content

@app.callback(
    Output("submit-output", "children"),
    [Input("submit-settings-button_1", "n_clicks")]
)
def render_loading_spinner(n_clicks):
    if n_clicks:
        # Simulate a loading process
        time.sleep(2)
        return f"Process completed {n_clicks} times."
    return no_update

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
