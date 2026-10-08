"""Tests for the Dash app shell: layout wiring, routing and upload helpers."""

import pytest

pytest.importorskip("dash", reason="the 'app' dependency group is not installed")

from dash.exceptions import PreventUpdate  # noqa: E402

from tests.conftest import CallbackRecorder  # noqa: E402

# --- the app shell -------------------------------------------------------


@pytest.fixture(scope="module")
def app_module():
    import application

    return application


def test_the_app_exposes_a_wsgi_server(app_module):
    assert app_module.application is app_module.app.server
    assert app_module.app.title == "StreptoCAD"


def test_every_workflow_callback_is_registered(app_module):
    registered = app_module.app.callback_map
    for n in range(1, 7):
        assert any(f"error-dialog_{n}.message" in key for key in registered), n


def test_the_health_check_endpoint_reports_ok(app_module):
    client = app_module.app.server.test_client()
    response = client.get("/health-check")
    assert response.status_code == 200
    assert b"Health Check OK" in response.data


def test_the_index_template_keeps_the_dash_placeholders(app_module):
    for placeholder in ("{%app_entry%}", "{%config%}", "{%scripts%}", "{%renderer%}"):
        assert placeholder in app_module.app.index_string
    assert "favicon.ico" in app_module.app.index_string


# --- routing and small callbacks ----------------------------------------
#
# These are defined inline in application.py, so they are reached through the
# app's callback map rather than a recorder.


def callback_named(app, name):
    for entry in app.callback_map.values():
        inner = entry["callback"]
        wrapped = getattr(inner, "__wrapped__", inner)
        if getattr(wrapped, "__name__", None) == name:
            return wrapped
    raise AssertionError(f"no callback named {name}")


def test_the_landing_page_switches_to_the_main_layout(app_module):
    display = callback_named(app_module.app, "display_main_page")

    assert display(None, ["intro"]) == ["intro"]
    assert display(1, ["intro"]) is app_module.main_layout


@pytest.mark.parametrize(
    "tab, expected",
    [
        ("workflow_1", "workflow_1_tab"),
        ("workflow_2", "crispr_cb_tab"),
        ("workflow_3", "golden_gate_tab"),
        ("workflow_4", "crispri_tab"),
        ("workflow_5", "gibson_tab"),
        ("workflow_6", "cas3_tab"),
        ("about", "about_streptocad_and_mission_statement"),
    ],
)
def test_each_tab_renders_its_own_content(app_module, tab, expected):
    render = callback_named(app_module.app, "render_tab_content")
    assert render(tab) is getattr(app_module, expected)


def test_an_unknown_tab_falls_back_to_the_welcome_message(app_module):
    render = callback_named(app_module.app, "render_tab_content")
    assert render("nope") is app_module.welcome_message_content


def test_advanced_settings_toggle(app_module):
    toggle = callback_named(app_module.app, "toggle_advanced_settings")
    assert toggle([1]) == {"display": "block"}
    assert toggle([]) == {"display": "none"}


def test_the_filename_display_reports_the_upload(app_module):
    update = callback_named(app_module.app, "update_filename_display")
    assert update("genome.gb") == "Uploaded file: genome.gb"
    assert update(None) == "No file uploaded"


def test_the_loading_spinner_only_reports_after_a_click(app_module, monkeypatch):
    import time

    monkeypatch.setattr(time, "sleep", lambda _seconds: None)
    spinner = callback_named(app_module.app, "render_loading_spinner")

    assert spinner(2) == "Process completed 2 times."
    # Without a click the output is left untouched.
    assert spinner(None) is not None


# --- the interactivity callbacks ----------------------------------------


@pytest.fixture
def interactivity():
    from callbacks_interactivity import register_interactivity_callbacks

    recorder = CallbackRecorder()
    register_interactivity_callbacks(recorder)
    return {cb.__name__: cb for cb in recorder.callbacks}


def test_three_interactivity_callbacks_are_registered(interactivity):
    assert set(interactivity) == {
        "display_uploaded_filename",
        "display_uploaded_filenames_tab2",
        "display_uploaded_filenames_tab3",
    }


@pytest.mark.parametrize(
    "filename, expected",
    [
        ("genome.gb", "Uploaded file: genome.gb"),
        ("guides.csv", "Uploaded file: guides.csv"),
        ("plasmid.gbk", "Uploaded file: plasmid.gbk"),
        ("notes.txt", "Invalid file type. Please upload a .csv, .gb, or .gbk file."),
        ("", "No file selected."),
        (None, "No file selected."),
    ],
)
def test_upload_filenames_are_validated_by_extension(interactivity, filename, expected):
    display = interactivity["display_uploaded_filename"]
    assert display([filename]) == [expected]


def test_an_empty_upload_list_prevents_an_update(interactivity):
    with pytest.raises(PreventUpdate):
        interactivity["display_uploaded_filename"]([])


@pytest.mark.parametrize("tab", ["tab2", "tab3"])
def test_the_per_tab_filename_displays_fall_back_to_a_placeholder(interactivity, tab):
    display = interactivity[f"display_uploaded_filenames_{tab}"]

    assert display("genome.gb", "vector.gbk") == ("genome.gb", "vector.gbk")
    assert display(None, None) == ("No file selected", "No file selected")


# --- the shared components ----------------------------------------------


def test_upload_component_carries_its_id():
    from components import upload_component

    component = upload_component("upload-1", {}, {}, {})
    assert component.id == "upload-1"
    assert component.multiple is False


def test_upload_component_with_display_pairs_an_output_div():
    from components import upload_component_with_display

    component = upload_component_with_display(
        {"type": "upload", "index": "genome"}, {}, {}, {}
    )
    upload, display = component.children
    assert upload.id == {"type": "upload", "index": "genome"}
    assert display.id == {"type": "filename-display", "index": "genome"}


@pytest.mark.parametrize(
    "filename, expected",
    [("genome.gb", "Uploaded file: genome.gb"), (None, "No file uploaded")],
)
def test_display_uploaded_filenames_helper(filename, expected):
    from components import display_uploaded_filenames

    div = display_uploaded_filenames("display-1", filename)
    assert div.children == expected
    assert div.id == "display-1"
