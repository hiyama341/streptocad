# components.py
from dash import dcc, html


def upload_component(component_id, text_style, link_style, upload_button_style):
    return dcc.Upload(
        id=component_id,
        children=html.Div([
            'Drag and Drop or ',
            html.A('Select File', style=link_style)
        ], style=text_style),
        style=upload_button_style,
        multiple=False
    )