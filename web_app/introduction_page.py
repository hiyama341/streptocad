from dash import dcc, html
import dash_bootstrap_components as dbc
from styling import text_style

# Introduction text
introduction_text = """
# **Welcome to StreptoCAD**

StreptoCAD is a toolbox supporting genome engineering workflows in Streptomyces.

If you have used StreptoCAD in a publication, please remember to cite us.
"""

# Introduction page
introduction_page = html.Div(
    [
        # Logo at the top center
        html.Img(
            src="/assets/StreptoCAD_logo_round1.png",
            style={
                'width': '150px',
                'marginBottom': '20px',
                'display': 'block',
                'marginLeft': 'auto',
                'marginRight': 'auto',
            }
        ),
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
                html.I(className="fas fa-brain", style={'fontSize': '5rem', 'color': '#ddd', 'marginRight': '10px'}),
                html.Span("+", style={'fontSize': '4rem', 'color': '#ddd', 'marginRight': '10px'}),
                html.I(className="fas fa-laptop-code", style={'fontSize': '5rem', 'color': '#ddd', 'marginRight': '10px'}),
                html.Span("=", style={'fontSize': '4rem', 'color': '#ddd', 'marginRight': '10px'}),
                html.I(className="fas fa-dna", style={'fontSize': '5rem', 'color': '#ddd'}),
            ],
            style={'textAlign': 'center', 'marginTop': '20px'}
        ),
        dbc.Button("Let's get started", id='start-button', color="primary", size="lg", className="d-block mx-auto mt-4"),
        html.Div(
            [
                html.Span("Citation: ", style={"fontSize": "18px", "color": "#ddd"}),
                html.Span(
                    "Levassor, Whitford, Petersen, Madsen, Blin, Weber, Frandsen; TBA journal 2024 ",
                    style={"fontSize": "18px", "color": "#ddd"}
                ),
                html.A(
                    "doi: TBA",
                    href="https://doi.org/TBA",
                    style={"fontSize": "18px", "color": "#00c39a", "marginLeft": "10px"}
                ),
            ],
            style={'textAlign': 'center', 'marginTop': '40px', 'padding': '20px', 'backgroundColor': '#2C3E50', 'borderRadius': '10px', 'boxShadow': '0 4px 6px rgba(0, 0, 0, 0.1)', 'width': '80%', 'margin': '20px auto'}
        )
    ],
    style={'height': '100vh', 'display': 'flex', 'flexDirection': 'column', 'justifyContent': 'center', 'alignItems': 'center', 'backgroundColor': '#2C3E50'}
)
