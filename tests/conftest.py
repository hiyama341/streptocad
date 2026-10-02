"""Shared fixtures, including a harness for driving the Dash workflow callbacks.

The callbacks are registered through ``app.callback(...)``, so calling them via
``app.callback_map`` would need a live Dash request context. ``CallbackRecorder``
stands in for the app and simply keeps the undecorated function, which lets the
tests call a whole workflow the way the browser does -- base64 uploads in,
DataTable payloads out.
"""

import base64
import contextlib
from pathlib import Path

import pytest

TEST_FILES = Path(__file__).parent / "test_files"


class CallbackRecorder:
    """Captures callbacks instead of registering them with a real Dash app."""

    def __init__(self):
        self.callbacks = []
        self.outputs = []

    def callback(self, *args, **kwargs):
        self.outputs.append(_flatten_outputs(args) + _flatten_outputs(kwargs.values()))

        def decorator(func):
            self.callbacks.append(func)
            return func

        return decorator

    @property
    def only_callback(self):
        assert (
            len(self.callbacks) == 1
        ), f"expected 1 callback, got {len(self.callbacks)}"
        return self.callbacks[0]

    @property
    def output_keys(self):
        """``component_id.property`` for each declared Output, in order."""
        assert len(self.outputs) == 1, f"expected 1 callback, got {len(self.outputs)}"
        return [f"{o.component_id}.{o.component_property}" for o in self.outputs[0]]

    def call(self, **kwargs):
        """Invoke the callback and label its results by declared Output."""
        result = self.only_callback(**kwargs)
        keys = self.output_keys
        assert len(result) == len(
            keys
        ), f"callback returned {len(result)} values for {len(keys)} declared outputs"
        return dict(zip(keys, result))


def _flatten_outputs(values):
    from dash.dependencies import Output

    found = []
    for value in values:
        if isinstance(value, Output):
            found.append(value)
        elif isinstance(value, (list, tuple)):
            found.extend(_flatten_outputs(value))
    return found


@pytest.fixture
def recorder():
    return CallbackRecorder()


def as_upload(filename):
    """Encode a fixture file the way dcc.Upload hands contents to a callback."""
    payload = base64.b64encode((TEST_FILES / filename).read_bytes()).decode()
    return f"data:application/octet-stream;base64,{payload}"


@pytest.fixture
def upload():
    return as_upload


@pytest.fixture(scope="session")
def genome_upload():
    # 25 MB of GenBank; encode it once for the whole session.
    return as_upload("Streptomyces_coelicolor_A3_chromosome.gb")


class _StubHTMLExporter:
    """Minimal stand-in for nbconvert's HTMLExporter."""

    def from_notebook_node(self, notebook, **kwargs):
        source = notebook.cells[0].source if notebook.cells else ""
        return f"<html><body>{source}</body></html>", {}


@contextlib.contextmanager
def stub_html_exporter():
    """Isolate the protocol markdown -> HTML step of ProjectDirectory.

    The installed nbconvert (7.16.4) is incompatible with the resolved mistune
    (3.2.0) and raises ``'MathBlockParser' object has no attribute
    'parse_axt_heading'``. Rendering the protocol markdown is incidental to the
    workflow logic under test, so the tests stub the exporter out.
    """
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr("nbconvert.HTMLExporter", _StubHTMLExporter)
        yield


@pytest.fixture
def html_exporter_stub():
    with stub_html_exporter():
        yield
