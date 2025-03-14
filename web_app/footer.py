from dash import html
import dash_bootstrap_components as dbc

# Footer content
footer_content = html.Footer(
    [
        html.Div(
            [
                html.Img(
                    src="/assets/dtu_logo.png",
                    style={"height": "100px", "marginRight": "20px"},
                ),  # DTU logo
                html.Span(
                    "Levassor, Whitford, Petersen, Blin, Weber, Frandsen; BioXiv 2024 ",
                    style={"fontSize": "18px", "color": "#ddd"},
                ),
                html.A(
                    "doi: https://doi.org/10.1101/2024.12.19.629370",
                    href="https://www.biorxiv.org/content/10.1101/2024.12.19.629370v1",
                    style={
                        "fontSize": "18px",
                        "color": "#00c39a",
                        "marginLeft": "10px",
                    },
                ),
                html.Img(
                    src="/assets/cfb_logo.png",
                    style={"height": "100px", "marginLeft": "20px"},
                ),  # CRB logo
            ],
            style={
                "textAlign": "center",
                "padding": "20px 0",
                "backgroundColor": "#1A242F",
                "color": "#ddd",
                "fontSize": "30px",
            },
        ),  # Center-align text, add padding, set background and text colors, and increase font size
        html.Div(
            [
                html.A(
                    "Contact Us",
                    href="mailto:luclev@dtu.dk",
                    style={
                        "marginRight": "20px",
                        "fontSize": "22px",
                        "color": "#ddd",
                        "textDecoration": "none",
                    },
                ),
                html.A(
                    [
                        html.I(
                            className="fab fa-github",
                            style={"fontSize": "30px", "color": "#ddd"},
                        )
                    ],
                    href="https://github.com/hiyama341",
                    style={"marginRight": "20px"},
                ),  # GitHub link
            ],
            style={
                "textAlign": "center",
                "padding": "10px 0",
                "backgroundColor": "#1A242F",
                "color": "#ddd",
                "fontSize": "18px",
            },
        ),  # Center-align text, add padding, set background and text colors, and set font size
    ]
)
