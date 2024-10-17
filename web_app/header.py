
import dash_bootstrap_components as dbc
from dash import html

header_content = dbc.Row(
    [
        dbc.Col(
            html.Div(
                [
                    html.Img(src='/assets/StreptoCAD_logo_round1.png', style={'height': '200px', 'marginRight': '20px'}),  # Adjust the height and margin to create space between the image and the text
                    html.Div(
                        [
                            html.H1(
                                "StreptoCAD",
                                style={'color': '#ddd', 'fontSize': '5rem', 'marginBottom': '0.5rem', 'fontWeight': 'bold'}  # Adjust font size, margin, and make text bold
                            ),
                            html.H2(
                                "An open-source software toolbox automating genome engineering design workflows in streptomycetes - version 0.1.0",
                                style={'color': '#ddd', 'fontSize': '2rem'}  # Adjust font size
                            ),
                        ],
                        style={'display': 'flex', 'flexDirection': 'column'}  # Stack title and text vertically
                    ),
                ],
                className="d-flex align-items-center justify-content-start"  # Use flexbox to align items horizontally and start them to the left
            ),
            width=12,  # Use the full width of the row
        ),
    ],
    style={'marginBottom': '10px', 'borderBottom': '1px solid #444', 'backgroundColor': '#1A242F', 'padding': '10px 0'}  # Adjust padding and margin
)

