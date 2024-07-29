from dash import dcc, html
import dash_bootstrap_components as dbc

introduction_text = """
# **Welcome to StreptoCAD**

StreptoCAD is a toolbox supporting genome engineering workflows in Streptomyces.

If you have used StreptoCAD in a publication, please remember to cite us.
"""

introduction_page = html.Div(
    [
        dcc.Markdown(
            introduction_text,
            style={
                'textAlign': 'center',
                'fontSize': '1.5rem',
                'lineHeight': '1.5',
                'padding': '20px',
                'borderRadius': '10px',
                'boxShadow': '0 4px 6px rgba(0, 0, 0, 0.1)',
                'color': '#ddd',
                'backgroundColor': '#2C3E50',
                'width': '80%',
                'margin': '0 auto'
            }
        ),
        html.Div(
            [
                html.I(className="fas fa-brain", style={'fontSize': '3rem', 'color': '#ddd', 'marginRight': '10px'}),
                html.Span("+", style={'fontSize': '2rem', 'color': '#ddd', 'marginRight': '10px'}),
                html.I(className="fas fa-laptop-code", style={'fontSize': '3rem', 'color': '#ddd', 'marginRight': '10px'}),
                html.Span("=", style={'fontSize': '2rem', 'color': '#ddd', 'marginRight': '10px'}),
                html.I(className="fas fa-dna", style={'fontSize': '3rem', 'color': '#ddd'}),
            ],
            style={'textAlign': 'center', 'marginTop': '20px'}
        ),
        dbc.Button("Let's get started", id='start-button', color="primary", size="lg", className="d-block mx-auto mt-4"),
    ],
    style={'height': '100vh', 'display': 'flex', 'flexDirection': 'column', 'justifyContent': 'center', 'alignItems': 'center', 'backgroundColor': '#2C3E50'}
)

