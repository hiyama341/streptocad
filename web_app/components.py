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

def display_uploaded_filenames(upload_id, filename):
    """
    Display the uploaded filename.
    
    Parameters:
    upload_id (str): The ID of the upload component.
    filename (str): The filename to display.
    
    Returns:
    html.Div: A div containing the filename.
    """
    if filename:
        return html.Div(f"Uploaded file: {filename}", id=upload_id, style={'color': '#ddd'})
    else:
        return html.Div("No file uploaded", id=upload_id, style={'color': '#ddd'})
