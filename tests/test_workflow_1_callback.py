"""End-to-end exercise of the workflow 1 Dash callback."""

import base64
import io
import zipfile

import pytest

pytest.importorskip("dash", reason="the 'app' dependency group is not installed")

from dash.exceptions import PreventUpdate  # noqa: E402

from streptocad.utils import polymerase_dict  # noqa: E402
from tests.conftest import CallbackRecorder, as_upload, stub_html_exporter  # noqa: E402
from callbacks.workflow_1 import register_workflow_1_callbacks  # noqa: E402

# The defaults the workflow 1 tab ships with.
UP_HOMOLOGY = "GGCGAGCAACGGAGGTACGGACAGG"
DW_HOMOLOGY = "CGCAAGCCGCCACTCGAACGGAAGG"
POLYMERASE = polymerase_dict["Q5 High-Fidelity 2X Master Mix"]

OUTPUTS = (
    "primers_data",
    "primers_columns",
    "pcr_data",
    "pcr_columns",
    "analyzed_data",
    "analyzed_columns",
    "download_link",
    "error_message",
    "display_error",
    "metadata_data",
    "metadata_columns",
)


def run_once(**overrides):
    recorder = CallbackRecorder()
    register_workflow_1_callbacks(recorder)

    kwargs = {
        "n_clicks": 1,
        "sequences_content": as_upload("GOE_regulators.gb"),
        "plasmid_content": as_upload("pOEX-PkasO.gb"),
        "sequences_filename": "GOE_regulators.gb",
        "plasmid_filename": "pOEX-PkasO.gb",
        "up_homology": UP_HOMOLOGY,
        "dw_homology": DW_HOMOLOGY,
        "chosen_polymerase": POLYMERASE,
        "melting_temperature": 65,
        "primer_concentration": 0.4,
        "restriction_enzymes": "StuI",
        "primer_anneal_len": 18,
    }
    kwargs.update(overrides)

    with stub_html_exporter():
        result = recorder.only_callback(**kwargs)
    return dict(zip(OUTPUTS, result))


@pytest.fixture(scope="module")
def result():
    """A single successful run; most assertions below read from it."""
    return run_once()


def test_registers_exactly_one_callback():
    recorder = CallbackRecorder()
    register_workflow_1_callbacks(recorder)
    assert recorder.only_callback.__name__ == "run_workflow"


def test_without_a_click_the_callback_bails_out():
    with pytest.raises(PreventUpdate):
        run_once(n_clicks=None)


def test_the_run_succeeds(result):
    assert result["display_error"] is False
    assert result["error_message"] == ""


@pytest.mark.parametrize(
    "data_key, columns_key",
    [
        ("primers_data", "primers_columns"),
        ("pcr_data", "pcr_columns"),
        ("analyzed_data", "analyzed_columns"),
        ("metadata_data", "metadata_columns"),
    ],
)
def test_every_table_has_rows_matching_its_columns(result, data_key, columns_key):
    data, columns = result[data_key], result[columns_key]
    assert data, f"{data_key} has no rows"
    assert columns, f"{columns_key} has no column definitions"
    assert {c["id"] for c in columns} <= set(data[0])


def test_primer_tables_cover_every_input_sequence(result):
    # GOE_regulators.gb holds 18 records, each needing a forward and reverse primer.
    templates = {row["template"] for row in result["pcr_data"]}
    assert len(templates) == 18
    assert len(result["primers_data"]) == 2 * len(templates)


def test_download_link_is_a_zip_with_the_expected_members(result):
    assert result["download_link"].startswith("data:application/zip;base64,")

    archive = base64.b64decode(result["download_link"].split(",", 1)[1])
    with zipfile.ZipFile(io.BytesIO(archive)) as zf:
        names = zf.namelist()
        basenames = {n.rsplit("/", 1)[-1] for n in names}

        # Results are shelved by what each file is, so paths carry the meaning.
        shelved = {n.split("/", 1)[1] for n in names if "/" in n}
        assert {
            "1_plasmids/plasmid_index.csv",
            "2_primers/pcr_primers.csv",
            "2_primers/oligo_order_idt.csv",
            "2_primers/oligo_order_idt_plate.xlsx",
            "4_analysis/primer_qc_hairpins.csv",
            "4_analysis/assembly_overview.txt",
            "4_analysis/run.log",
            "4_analysis/environment.json",
            "6_inputs/input_sequences.gb",
            "6_inputs/input_plasmid.gb",
            "6_inputs/run_parameters.json",
            "00_READ_ME_FIRST.md",
            "00_READ_ME_FIRST.html",
        } <= shelved

        # Workflow 1 designs no guide RNAs, so that shelf is left out entirely.
        assert not [n for n in shelved if n.startswith("3_sgrnas/")]

        # One GenBank file per assembled plasmid, numbered from 01 so they sort
        # in design order and join to plasmid_index.csv row n.
        plasmids = sorted(
            n for n in shelved if n.startswith("1_plasmids/") and n.endswith(".gb")
        )
        assert len(plasmids) == 18
        assert plasmids[0].startswith("1_plasmids/01_")
        assert plasmids[-1].startswith("1_plasmids/18_")

        # 18 plasmids plus the two inputs.
        assert len([n for n in basenames if n.endswith(".gb")]) == 18 + 2
        # Three protocols, each as .html and .md, plus the generated README.
        assert len([n for n in shelved if n.startswith("5_protocols/")]) == 6
        assert len([n for n in names if n.endswith(".html")]) == 3 + 1

        overview = zf.read(
            next(n for n in names if n.endswith("4_analysis/assembly_overview.txt"))
        ).decode()
        log = zf.read(
            next(n for n in names if n.endswith("4_analysis/run.log"))
        ).decode()

    # The overview is the header plus this run's prints.
    assert overview.startswith("StreptoCAD generated ")
    # The log is this run's own records, shipped separately from the overview.
    assert "Workflow 1 started" in log
    assert "Generating primers" in log
    assert "Assembling plasmids" in log


def test_the_log_file_only_covers_the_current_run():
    first = run_once()
    second = run_once(restriction_enzymes="StuI")

    def log_of(result):
        archive = base64.b64decode(result["download_link"].split(",", 1)[1])
        with zipfile.ZipFile(io.BytesIO(archive)) as zf:
            name = next(
                n for n in zf.namelist() if n.endswith("4_analysis/run.log")
            )
            return zf.read(name).decode()

    # Before the per-run capture, stdout was restored after the first run and
    # the second run shipped the first run's output.
    assert log_of(first).count("Workflow 1 started") == 1
    assert log_of(second).count("Workflow 1 started") == 1


def test_an_unsupported_upload_is_reported_with_the_run_log():
    result = run_once(
        sequences_filename="sequences.txt", plasmid_filename="plasmid.txt"
    )

    assert result["display_error"] is True
    assert "Unsupported file format" in result["error_message"]
    # The error dialog carries this run's log, which was always empty before.
    assert "Workflow 1 started" in result["error_message"]
    assert result["download_link"] == ""
    assert result["primers_data"] == []


def test_an_unknown_restriction_enzyme_is_reported():
    result = run_once(restriction_enzymes="NotAnEnzyme")

    assert result["display_error"] is True
    assert "An error occurred" in result["error_message"]
    assert "Workflow 1 started" in result["error_message"]
