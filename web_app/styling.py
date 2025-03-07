#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# styling.py

##########################
# Define the common styles
##########################
tab_style = {
    'borderBottom': '2px solid #888',
    'padding': '10px 20px',
    'fontWeight': 'bold',
    'borderRadius': '4px',
    'backgroundColor': '#3A475B',
    'marginRight': '5px',
    'fontSize': '18px',  # Increase font size
    'color': '#ddd',
    'boxShadow': '0px 0px 10px rgba(0,0,0,0.3)',
}
active_tab_style = {
    'padding': '10px 20px',
    'backgroundColor': '#1A242F',
    'color': '#33C3F0',
    'borderRadius': '4px 0px 0px 0px',
    'fontSize': '18px',  # Increase font size
    'boxShadow': '0px 0px 15px rgba(0,0,0,0.5)',
    'border': 'none',
    'borderTop': 'none',
    'borderRight': 'none',
    'borderBottom': 'none',
    'borderLeft': 'none'
}
vertical_tab_style = {
    'display': 'block',
    'padding': '10px 20px',
    'fontSize': '18px',  # Increase font size
    'backgroundColor': '#3A475B',
    'color': '#bbb',
    'borderLeft': '4px solid transparent',
    'transition': 'color 0.3s, backgroundColor 0.3s',
    'width': '100%',  # Ensure consistent width
    ':hover': {
        'color': '#ddd',
        'backgroundColor': '#1A242F',
        'borderLeft': '4px solid #F39C12'
    }
}
# Text styles
text_style = {'color': '#ddd'}

# Upload button styles
upload_button_style = {
    'width': '100%',
    'height': '60px',
    'lineHeight': '60px',
    'borderWidth': '1px',
    'borderStyle': 'dashed',
    'borderColor': '#ddd',
    'borderRadius': '5px',
    'textAlign': 'center',
    'margin': '10px 0',
    'backgroundColor': '#2C3E50'
}

# Card styles
card_style = {
    'width': '100%',
    'backgroundColor': '#34495E'
}

# Link styles
link_style = {'color': '#3498DB', 'textDecoration': 'underline'}

# Button styles for DARKLY theme
button_style_darkly = {
    'backgroundColor': '#3498DB',
    'color': '#fff',
    'border': 'none',
    'padding': '10px 20px',
    'textDecoration': 'none',
    'display': 'inline-block',
    'fontSize': '16px',
    'borderRadius': '4px',
    'cursor': 'pointer',
}

# Table styles
table_style = {
    'style_table': {
        'maxHeight': '500px',
        'overflowY': 'auto',
        'border': 'none',  # Remove the border to eliminate the white line
        'backgroundColor': '#2C3E50',
    },
    'style_header': {
        'backgroundColor': '#1A242F',
        'fontWeight': 'bold',
        'color': '#fff',
        'border': 'none',  # Remove the border to eliminate the white line
        'textAlign': 'center'
    },
    'style_cell': {
        'padding': '10px',
        'backgroundColor': '#34495E',
        'color': '#ddd',
        'border': 'none',  # Remove the border to eliminate the white line
        'textAlign': 'left',
    },
    'style_data_conditional': [
        {
            'if': {'row_index': 'odd'},
            'backgroundColor': '#3A475B'
        },
        {
            'if': {'row_index': 'even'},
            'backgroundColor': '#2C3E50'
        },
        {
            'if': {'column_id': 'some_column'},
            'textAlign': 'center'
        },
        {
            'if': {'state': 'active'},  # Active cell styling
            'backgroundColor': '#1A242F',
            'border': '1px solid #3498DB',
            'color': '#fff'
        },
        {
            'if': {'state': 'selected'},  # Selected cell styling
            'backgroundColor': '#1A242F',
            'border': '1px solid #3498DB',
            'color': '#fff'
        }
    ],
    'style_as_list_view': True,
    'style_data': {
        'border': 'none'  # Remove the border to eliminate the white line
    },
}

# Specific header style for tables
table_header_style = {
    'backgroundColor': '#1A242F',
    'fontWeight': 'bold',
    'color': '#FFF'
}

# Specific row style for tables
table_row_style = {
    'backgroundColor': '#2C3E50',
    'color': '#ddd',
    'minWidth': '150px',
    'width': '150px',
    'maxWidth': '150px',
    'whiteSpace': 'normal'
}