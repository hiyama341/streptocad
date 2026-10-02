"""End-to-end exercise of the workflow 2-6 Dash callbacks.

Each case feeds the callback the defaults its tab ships with, so a run here is
the same work the browser triggers. Results are keyed by the ``Output`` the
callback declares, which keeps the assertions readable as the tuples grow.
"""

import base64
import io
import zipfile

import pytest

pytest.importorskip("dash", reason="the 'app' dependency group is not installed")

from dash.exceptions import PreventUpdate  # noqa: E402

from streptocad.cloning.cas3_plasmid_cloning import (  # noqa: E402
    CAS3_BACKBONE_FWD_PRIMER,
    CAS3_BACKBONE_REV_PRIMER,
    CAS3_PROTOSPACER_FWD_OVERHANG,
    CAS3_PROTOSPACER_REV_OVERHANG,
)
from streptocad.utils import polymerase_dict  # noqa: E402
from tests.conftest import CallbackRecorder, as_upload, stub_html_exporter  # noqa: E402
from callbacks.workflow_2 import register_workflow_2_callbacks  # noqa: E402
from callbacks.workflow_3 import register_workflow_3_callbacks  # noqa: E402
from callbacks.workflow_4 import register_workflow_4_callbacks  # noqa: E402
from callbacks.workflow_5 import register_workflow_5_callbacks  # noqa: E402
from callbacks.workflow_6 import register_workflow_6_callbacks  # noqa: E402

GENOME = "Streptomyces_coelicolor_A3_chromosome.gb"
POLYMERASE = polymerase_dict["Q5 High-Fidelity 2X Master Mix"]

# Settings every CRISPy workflow shares, at their tab defaults.
SGRNA_DEFAULTS = {
    "genes_to_KO": "SCO5087",
    "gc_upper": 0.8,
    "gc_lower": 0.2,
    "off_target_seed": 13,
    "off_target_upper": 10,
    "number_of_sgRNAs_per_group": 5,
}
PCR_DEFAULTS = {
    "chosen_polymerase": POLYMERASE,
    "melting_temperature": 65,
    "primer_concentration": 0.4,
    "checking_primer_length": 18,
}

CASES = {
    "workflow_2": dict(
        register=register_workflow_2_callbacks,
        vector="pCRISPR-cBEST.gbk",
        kwargs={
            **SGRNA_DEFAULTS,
            **PCR_DEFAULTS,
            "up_homology": "CGGTTGGTAGGATCGACGGC",
            "dw_homology": "GTTTTAGAGCTAGAAATAGC",
            "cas_type": "cas9",
            "only_stop_codons": [],
            "flanking_region_number": 200,
            "editing_context": [],
            "restriction_enzymes": "NcoI",
        },
    ),
    "workflow_3": dict(
        register=register_workflow_3_callbacks,
        vector="pCRISPR-MCBE_Csy4_kasopGFP.gb",
        kwargs={
            **SGRNA_DEFAULTS,
            **PCR_DEFAULTS,
            "sgRNA_handle_input": (
                "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGG"
                "CACCGAGTCGGTGCTTTTTTGTTCACTGCCGTATAGGCAGCTAAGAAA"
            ),
            "input_tm": 60,
            "cas_type": "cas9",
            "only_stop_codons": [],
            "flanking_region_number": 200,
            "restriction_overhang_f": "GATCGggtctcc",
            "restriction_overhang_r": "GATCAGGTCTCg",
            "backbone_overhang_f": "cATG",
            "backbone_overhang_r": "cTAG",
            "cys4_sequence": "gTTCACTGCCGTATAGGCAGCTAAGAAA",
            "editing_context": [],
            "restriction_enzymes": "NcoI,NheI",
        },
    ),
    "workflow_4": dict(
        register=register_workflow_4_callbacks,
        vector="pCRISPR-cBEST.gbk",
        kwargs={
            **SGRNA_DEFAULTS,
            "up_homology": "CGGTTGGTAGGATCGACGGC",
            "dw_homology": "GTTTTAGAGCTAGAAATAGC",
            "cas_type": "cas9",
            "extension_to_promoter_region": 100,
            "restriction_enzymes": "NcoI",
        },
    ),
    "workflow_5": dict(
        register=register_workflow_5_callbacks,
        vector="pCRISPR-Cas9.gbk",
        kwargs={
            **SGRNA_DEFAULTS,
            **PCR_DEFAULTS,
            "up_homology": "CGGTTGGTAGGATCGACGGC",
            "dw_homology": "GTTTTAGAGCTAGAAATAGC",
            "cas_type": "cas9",
            "in_frame_deletion": [],
            "flanking_region_number": 1000,
            "repair_templates_length": 1000,
            "overlap_for_gibson_length": 40,
            "restriction_enzymes": "NcoI",
            "restriction_enzymes_2": "StuI",
            "gibson_primer_length": 18,
        },
    ),
    "workflow_6": dict(
        register=register_workflow_6_callbacks,
        vector="pCRISPR_cas3.gbk",
        kwargs={
            **SGRNA_DEFAULTS,
            **PCR_DEFAULTS,
            "forward_protospacer_overhang": CAS3_PROTOSPACER_FWD_OVERHANG.upper(),
            "reverse_protospacer_overhang": CAS3_PROTOSPACER_REV_OVERHANG,
            "backbone_fwd_overhang": CAS3_BACKBONE_FWD_PRIMER,
            "backbone_rev_overhang": CAS3_BACKBONE_REV_PRIMER,
            "cas_type": "cas3",
            "in_frame_deletion": [],
            "flanking_region_number": 1000,
            "restriction_enzymes_2": "EcoRI",
            "repair_templates_length": 1000,
            "overlap_for_gibson_length": 40,
            "gibson_primer_length": 18,
        },
    ),
}


def run_case(name, **overrides):
    case = CASES[name]
    recorder = CallbackRecorder()
    case["register"](recorder)

    kwargs = {
        "n_clicks": 1,
        "genome_content": as_upload(GENOME),
        "vector_content": as_upload(case["vector"]),
        "genome_filename": GENOME,
        "vector_filename": case["vector"],
        **case["kwargs"],
    }
    kwargs.update(overrides)

    with stub_html_exporter():
        return recorder.call(**kwargs)


@pytest.fixture(scope="module")
def results():
    """Run each workflow once and share the results across the assertions."""
    return {name: run_case(name) for name in CASES}


def number_of(name):
    return name.rsplit("_", 1)[1]


@pytest.mark.parametrize("name", list(CASES))
def test_registers_exactly_one_callback(name):
    recorder = CallbackRecorder()
    CASES[name]["register"](recorder)
    assert recorder.only_callback.__name__ == "run_workflow"
    assert recorder.output_keys, "callback declares no outputs"


@pytest.mark.parametrize("name", list(CASES))
def test_every_return_matches_the_declared_outputs(name):
    """Dash rejects a callback whose return arity differs from its Outputs.

    Workflow 4's error path used to return 11 values against 9 declared
    outputs, so its error dialog could never render.
    """
    import ast
    import inspect

    module = inspect.getmodule(CASES[name]["register"])
    tree = ast.parse(inspect.getsource(module))
    run = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.FunctionDef) and node.name == "run_workflow"
    )
    returns = [
        node
        for node in ast.walk(run)
        if isinstance(node, ast.Return)
        and isinstance(node.value, (ast.Tuple, ast.List))
    ]

    recorder = CallbackRecorder()
    CASES[name]["register"](recorder)
    expected = len(recorder.output_keys)

    assert returns, "no tuple returns found"
    for node in returns:
        assert len(node.value.elts) == expected, (
            f"{name} line {node.lineno} returns {len(node.value.elts)} "
            f"values for {expected} declared outputs"
        )


@pytest.mark.parametrize("name", list(CASES))
def test_protocol_paths_exist(name):
    """A missing comma once concatenated two of workflow 4's protocol paths."""
    import ast
    import inspect
    from pathlib import Path

    module = inspect.getmodule(CASES[name]["register"])
    tree = ast.parse(inspect.getsource(module))

    paths = [
        element.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(getattr(t, "id", "") == "markdown_file_paths" for t in node.targets)
        for element in node.value.elts
        if isinstance(element, ast.Constant)
    ]

    assert paths, "no protocol paths found"
    for path in paths:
        assert Path(path).is_file(), f"{name} references a missing protocol: {path}"


@pytest.mark.parametrize("name", list(CASES))
def test_without_a_click_the_callback_bails_out(name):
    recorder = CallbackRecorder()
    CASES[name]["register"](recorder)
    with pytest.raises(PreventUpdate):
        recorder.only_callback(**{"n_clicks": None, **_placeholders(name)})


def _placeholders(name):
    """Every argument other than n_clicks; none is read before the bail-out."""
    kwargs = dict(CASES[name]["kwargs"])
    kwargs.update(
        genome_content=None,
        vector_content=None,
        genome_filename=None,
        vector_filename=None,
    )
    return kwargs


@pytest.mark.parametrize("name", list(CASES))
def test_the_run_succeeds(results, name):
    result = results[name]
    n = number_of(name)
    assert result[f"error-dialog_{n}.displayed"] in (False, None), result[
        f"error-dialog_{n}.message"
    ]
    assert result[f"error-dialog_{n}.message"] == ""


@pytest.mark.parametrize("name", list(CASES))
def test_every_table_has_rows_matching_its_columns(results, name):
    result = results[name]
    table_keys = sorted({k.rsplit(".", 1)[0] for k in result if k.endswith(".data")})
    assert table_keys, "no tables declared"

    for table in table_keys:
        data, columns = result[f"{table}.data"], result[f"{table}.columns"]
        assert data, f"{table} has no rows"
        assert columns, f"{table} has no column definitions"
        assert {c["id"] for c in columns} <= set(data[0]), table


@pytest.mark.parametrize("name", list(CASES))
def test_download_link_carries_a_zip_with_the_run_log(results, name):
    n = number_of(name)
    link = results[name][f"download-data-and-protocols-link_{n}.href"]
    assert link.startswith("data:application/zip;base64,")

    archive = base64.b64decode(link.split(",", 1)[1])
    with zipfile.ZipFile(io.BytesIO(archive)) as zf:
        names = zf.namelist()
        assert any(name_.endswith(".csv") for name_ in names), names
        assert any(name_.endswith(".gb") for name_ in names), names
        # The protocol markdown is rendered into the package.
        assert any(name_.endswith(".html") for name_ in names), names


@pytest.mark.parametrize("name", list(CASES))
def test_an_unsupported_upload_is_reported_with_the_run_log(results, name):
    n = number_of(name)
    result = run_case(name, genome_filename="genome.txt", vector_filename="vector.txt")

    assert result[f"error-dialog_{n}.displayed"] is True
    message = result[f"error-dialog_{n}.message"]
    assert "error occurred" in message.lower()
    # The dialog carries this run's captured log, which was empty before the
    # per-run capture replaced the shared module-level buffer.
    assert "Log:" in message
    assert message.split("Log:", 1)[1].strip(), "the log section is empty"


@pytest.mark.parametrize("name", list(CASES))
def test_an_unknown_gene_yields_an_error_or_no_rows(results, name):
    """A gene that is not in the genome must not produce fabricated output.

    Most workflows raise and surface the error dialog. Workflow 4 instead
    returns empty tables and a downloadable package with no sgRNAs in it.
    """
    n = number_of(name)
    result = run_case(name, genes_to_KO="NOT_A_GENE")

    if result[f"error-dialog_{n}.displayed"]:
        assert result[f"error-dialog_{n}.message"]
        return

    tables = {k.rsplit(".", 1)[0] for k in result if k.endswith(".data")}
    for table in tables:
        assert (
            result[f"{table}.data"] == []
        ), f"{table} produced rows for a missing gene"


# The tab switches (in-frame deletion, stop-codon-only, editing context) each
# select a different branch through the workflow, so they get their own runs.
TOGGLES = [
    ("workflow_2", {"only_stop_codons": [1]}),
    ("workflow_2", {"editing_context": [1]}),
    ("workflow_3", {"only_stop_codons": [1]}),
    ("workflow_5", {"in_frame_deletion": [1]}),
    ("workflow_6", {"in_frame_deletion": [1]}),
]


@pytest.mark.parametrize(
    "name, toggle", TOGGLES, ids=[f"{n}-{list(t)[0]}" for n, t in TOGGLES]
)
def test_switches_select_a_working_branch(name, toggle):
    n = number_of(name)
    result = run_case(name, **toggle)

    assert result[f"error-dialog_{n}.displayed"] in (False, None), result[
        f"error-dialog_{n}.message"
    ]
    assert result[f"download-data-and-protocols-link_{n}.href"].startswith(
        "data:application/zip;base64,"
    )
