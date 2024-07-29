import dash_bootstrap_components as dbc
from dash import html

header_content = dbc.Row(
    [
        # Single column for the header text and the logo
        dbc.Col(
            html.Div(
                [
                    html.H2(
                        "StreptoCAD: An open-source software toolbox supporting genome engineering workflows in Streptomyces - version 0.1.0",
                        style={'color': '#ddd', 'fontSize': '2.5rem', 'marginRight': 'auto'}  # Adjust font size and margin to push the logo to the right
                    ),
                    html.Img(src='/assets/dtu_logo.png', style={'height': '100px', 'marginLeft': 'auto'})  # Adjust the height and margin to push the logo to the right
                ],
                className="d-flex align-items-center justify-content-between"  # Use flexbox to align items horizontally and center them vertically
            ),
            width=12,  # Use the full width of the row
        ),
    ],
    style={'marginBottom': '10px', 'borderBottom': '1px solid #444', 'backgroundColor': '#1A242F', 'padding': '10px 0'}  # Adjust padding and margin
)
