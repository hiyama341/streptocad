from dash import dcc, html
from styling import text_style

intro_message_content = html.Div([
    dcc.Markdown(
        """
        ## **Welcome to StreptoCAD**
        """,
        style={**text_style, "marginBottom": "20px", "lineHeight": "2", "marginTop": "30px"}
    ),
    html.Img(
        src="assets/intro_fig.png",  
        style={"width": "30%", "height": "auto", "display": "block", "marginLeft": "0", "marginRight": "auto"}  
    )
])
