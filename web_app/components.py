# components.py
from dash import dcc, html

# TODO Remove everywhere
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

# TODO Remove everywhere
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


# components.py
def upload_component_with_display(id, text_style, link_style, upload_button_style):
    """
    Returns an upload component with an associated display area for the filename.

    Parameters:
    - component_id (str): The ID for the upload component.
    - text_style (dict): The style for the text in the component.
    - link_style (dict): The style for the link in the component.
    - upload_button_style (dict): The style for the upload button.
    - display_id (str): The ID for the filename display Div.

    Returns:
    - html.Div: A Div containing the upload component and the filename display area.
    """
    return html.Div([
        dcc.Upload(
            id=id,  # Directly use the pattern-matching ID here
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select File', style=link_style)
            ], style=text_style),
            style=upload_button_style,
            multiple=False
        ),
        html.Div(id={'type': 'filename-display', 'index': id['index']}, children=[], style=text_style)  # Ensure the filename display also uses pattern-matching IDs
    ])
