from dash import html
import dash_bootstrap_components as dbc

# Footer content
footer_content = html.Footer([
    html.Div([
        html.Img(src="/assets/dtu_logo.png", style={"height": "50px", "marginRight": "20px"}),  # DTU logo
        html.Span(
            "Please cite: Levassor, Whitford, Petersen, Madsen, Blin, Weber, Frandsen; TBA journal 2024 ",
            style={"fontSize": "18px", "color": "#ddd"}
        ),
        html.A(
            "doi: TBA",
            href="https://doi.org/TBA",
            style={"fontSize": "18px", "color": "#00c39a", "marginLeft": "10px"}
        ),
        html.Img(src="/assets/cfb_logo.png", style={"height": "50px", "marginLeft": "20px"})  # CRB logo
    ], style={"textAlign": "center", "padding": "20px 0", "backgroundColor": "#1A242F", "color": "#ddd", "fontSize": "30px"}),  # Center-align text, add padding, set background and text colors, and increase font size
    html.Div([
        html.A("Contact Us", href="mailto:your-email@example.com", style={"marginRight": "20px", "fontSize": "18px", "color": "#ddd"}),  # Email link
        html.A([html.I(className="fab fa-github", style={"fontSize": "30px", "color": "#ddd"})], href="https://github.com/YourGitHubProfile", style={"marginRight": "20px"}),  # GitHub link
    ], style={"textAlign": "center", "padding": "10px 0", "backgroundColor": "#1A242F", "color": "#ddd", "fontSize": "18px"})  # Center-align text, add padding, set background and text colors, and set font size
])