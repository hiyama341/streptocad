#!/usr/bin/env python
# MIT License
# Copyright (c) 2024, Technical University of Denmark (DTU)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# callbacks_interactivity.py

from dash.dependencies import Input, Output, State, ALL
from dash.exceptions import PreventUpdate

def register_interactivity_callbacks(app):
    ##### FILE EXTENSION functions for error handling
    def get_file_extension(filename):
        """Get the file extension."""
        return filename.split('.')[-1] if filename else ""

    def validate_file_extension(filename):
        """Validate the file extension and return appropriate message."""
        if not filename:
            return "No file selected."
        extension = get_file_extension(filename)
        if extension not in ['csv', 'gb', 'gbk']:
            return "Invalid file type. Please upload a .csv, .gb, or .gbk file."
        return f"Uploaded file: {filename}"

    @app.callback(
        Output({'type': 'uploaded-filename', 'index': ALL}, 'children'),
        [Input({'type': 'upload-input', 'index': ALL}, 'filename')]
    )
    def display_uploaded_filename(filenames):
        if not filenames:
            raise PreventUpdate
        return [validate_file_extension(filename) for filename in filenames]

    # Add the callback specifically for tab2
    @app.callback(
        [
            Output('uploaded-genome-filename_2', 'children'),
            Output('uploaded-single-vector-filename_2', 'children')
        ],
        [
            Input('upload-genome-file_2', 'filename'),
            Input('upload-single-vector_2', 'filename')
        ]
    )
    def display_uploaded_filenames(genome_filename, vector_filename):
        return genome_filename or "No file selected", vector_filename or "No file selected"
