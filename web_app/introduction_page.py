from dash import html
import dash_bootstrap_components as dbc

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
