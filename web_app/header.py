import dash_bootstrap_components as dbc
from dash import html

header_content = dbc.Row(
    [
        # Single column for the header text and the logo
        dbc.Col(
            html.Div(
                [
                    html.Img(src='/assets/StreptoCAD_logo_round.png', style={'height': '200px', 'marginRight': '20px'}),  # Adjust the height and margin to create space between the image and the text
                    html.H2(
                        "StreptoCAD: An open-source software toolbox automating genome engineering design workflows in streptomycetes - version 0.1.0",
                        style={'color': '#ddd', 'fontSize': '2.5rem'}  # Adjust font size, removed marginRight to avoid pushing text
                    ),
                ],
                className="d-flex align-items-center justify-content-start"  # Use flexbox to align items horizontally and start them to the left
            ),
            width=12,  # Use the full width of the row
        ),
    ],
    style={'marginBottom': '10px', 'borderBottom': '1px solid #444', 'backgroundColor': '#1A242F', 'padding': '10px 0'}  # Adjust padding and margin
)
