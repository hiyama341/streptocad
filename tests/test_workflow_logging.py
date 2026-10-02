import importlib.util
import io
import logging
import os
import subprocess
import sys
import threading
from pathlib import Path

import pytest

from callbacks.workflow_logging import (
    capture_output,
    capture_workflow_output,
    configure_logging,
    workflow_log,
)

configure_logging()

REPO_ROOT = Path(__file__).resolve().parents[1]


class PreventUpdateStub(Exception):
    """Stands in for dash.exceptions.PreventUpdate, which workflows raise."""


# Importing the web app must never swap the process's stdout for a buffer: that
# swallows the container logs and the buffer then grows for the whole life of
# the server process.
STDOUT_CHECK = (
    "import sys, {module}; "
    "assert isinstance(sys.stdout, __import__('io').TextIOWrapper), sys.stdout"
)


def _import_in_subprocess(module, pythonpath):
    env = dict(os.environ, PYTHONPATH=os.pathsep.join(str(p) for p in pythonpath))
    return subprocess.run(
        [sys.executable, "-c", STDOUT_CHECK.format(module=module)],
        env=env,
        capture_output=True,
        text=True,
    )


def test_importing_the_helper_leaves_stdout_alone():
    result = _import_in_subprocess(
        "callbacks.workflow_logging", [REPO_ROOT, REPO_ROOT / "web_app"]
    )
    assert result.returncode == 0, result.stderr


@pytest.mark.skipif(
    importlib.util.find_spec("dash") is None,
    reason="the 'app' dependency group is not installed",
)
def test_importing_the_dash_app_leaves_stdout_alone():
    result = _import_in_subprocess("application", [REPO_ROOT, REPO_ROOT / "web_app"])
    assert result.returncode == 0, result.stderr


def test_capture_collects_prints_and_log_records():
    with capture_output() as buffer:
        print("printed line")
        logging.getLogger(__name__).info("logged line")
        assert "printed line" in workflow_log()

    assert "printed line" in buffer.getvalue()
    assert "logged line" in buffer.getvalue()


def test_capture_restores_stdout():
    before = sys.stdout
    with capture_output():
        assert sys.stdout is not before
    assert sys.stdout is before


def test_capture_leaves_output_on_the_real_stdout(capsys):
    with capture_output() as buffer:
        print("printed line")

    # The run gets its own copy, and the real stdout still sees it.
    assert "printed line" in buffer.getvalue()
    assert "printed line" in capsys.readouterr().out


def test_each_run_starts_from_an_empty_buffer():
    @capture_workflow_output
    def run(tag):
        print(tag)
        return workflow_log()

    assert "first" in run("first")
    second = run("second")
    assert "second" in second
    assert "first" not in second


def test_the_log_survives_into_the_error_handler():
    # Every workflow builds its error dialog from workflow_log() inside an
    # `except` block, so the buffer has to still be live there.
    @capture_workflow_output
    def run():
        try:
            print("work so far")
            raise ValueError("boom")
        except ValueError as exc:
            return f"{exc}\n\nLog:\n{workflow_log()}"

    message = run()
    assert "boom" in message
    assert "work so far" in message


def test_an_escaping_exception_still_restores_stdout():
    @capture_workflow_output
    def run():
        raise PreventUpdateStub

    before = sys.stdout
    with pytest.raises(PreventUpdateStub):
        run()
    assert sys.stdout is before


def test_nothing_is_captured_outside_a_run():
    print("outside a run")
    logging.getLogger(__name__).info("outside a run")
    assert workflow_log() == ""


def test_concurrent_runs_do_not_leak_into_each_other():
    tags = ["w", "x", "y", "z"]
    barrier = threading.Barrier(len(tags))
    results = {}

    @capture_workflow_output
    def run(tag):
        print(f"print-{tag}")
        barrier.wait(timeout=10)  # force the runs to overlap
        logging.getLogger(__name__).info("log-%s", tag)
        return workflow_log()

    def worker(tag):
        results[tag] = run(tag)

    threads = [threading.Thread(target=worker, args=(tag,)) for tag in tags]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    for tag in tags:
        assert f"print-{tag}" in results[tag]
        assert f"log-{tag}" in results[tag]
        for other in set(tags) - {tag}:
            assert f"print-{other}" not in results[tag]


def test_generate_header_appends_the_captured_output():
    from streptocad.utils import generate_header

    header = generate_header(["plasmid"], ["sequence"], [1, 2], "captured text\n")

    assert header.startswith("StreptoCAD generated 1 plasmids from 1 sequences")
    assert header.endswith("captured text\n")


def test_the_tee_delegates_file_attributes_to_the_real_stdout():
    before = sys.stdout
    with capture_output():
        tee = sys.stdout
        assert tee is not before
        # Libraries probe these on sys.stdout; they must reach the real stream.
        assert tee.encoding == before.encoding
        assert tee.isatty() == before.isatty()
        tee.flush()


def test_writelines_reaches_both_the_buffer_and_stdout(capsys):
    with capture_output() as buffer:
        sys.stdout.writelines(["one\n", "two\n"])

    assert buffer.getvalue() == "one\ntwo\n"
    assert "one\ntwo\n" in capsys.readouterr().out


def test_nested_capture_restores_the_outer_buffer():
    with capture_output() as outer:
        print("outer before")
        with capture_output() as inner:
            print("inner only")
        print("outer after")

    assert "inner only" in inner.getvalue()
    assert "inner only" not in outer.getvalue()
    assert "outer before" in outer.getvalue()
    assert "outer after" in outer.getvalue()


def test_configure_logging_is_idempotent():
    import logging as logging_module

    from callbacks.workflow_logging import configure_logging as configure

    before = list(logging_module.root.handlers)
    configure()
    configure()
    assert logging_module.root.handlers == before


def test_only_one_copy_of_the_capture_module_is_live():
    """The app must be imported under a single package name.

    With ``web_app`` on the path it is importable as both ``callbacks.*`` and
    ``web_app.callbacks.*``. Each copy owns its own buffer registry, and the
    second copy's ``configure_logging()`` removes the first copy's handler --
    which silently stops the run log reaching the UI. Production imports
    ``callbacks.*`` (the Dockerfile puts ``web_app`` on PYTHONPATH), so the
    tests must too.
    """
    import logging as logging_module

    import application  # noqa: F401  -- registers all six workflow callbacks

    handlers = [
        h
        for h in logging_module.root.handlers
        if type(h).__name__ == "_RunCaptureHandler"
    ]
    assert len(handlers) == 1, f"expected one capture handler, found {len(handlers)}"

    duplicates = [
        name
        for name in sys.modules
        if name.endswith(".callbacks.workflow_logging")
        or name == "web_app.callbacks.workflow_logging"
    ]
    assert not duplicates, f"the capture module was imported twice as {duplicates}"


def test_the_live_capture_handler_feeds_the_callbacks_buffer():
    """The handler on the root logger must belong to the copy the callbacks use."""
    import logging as logging_module

    import callbacks.workflow_1  # noqa: F401
    from callbacks.workflow_logging import capture_output as app_capture
    from callbacks.workflow_logging import workflow_log as app_log

    with app_capture() as buffer:
        logging_module.getLogger("streptocad.test").info("reaches the UI")
        assert "reaches the UI" in app_log()

    assert "reaches the UI" in buffer.getvalue()
