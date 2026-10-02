"""Per-run capture of everything a workflow prints and logs.

Every workflow surfaces its own output to the user: workflow 1 ships it inside
the downloaded project as ``03_assembly_overview_and_log_file.log``, and all of
them show it in the error dialog when a run fails.

Capture is scoped to a single run and to the thread running it, so the process
keeps its real ``sys.stdout`` -- the container log -- and no buffer outlives the
run that filled it.
"""

import contextlib
import functools
import io
import logging
import sys
import threading

LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

# Buffer of the workflow run currently active on this thread, if any.
_local = threading.local()

# sys.stdout is only swapped while at least one run is in flight.
_stdout_lock = threading.Lock()
_active_runs = 0
_real_stdout = None

_configured = False


def _current_buffer():
    return getattr(_local, "buffer", None)


class _TeeStdout:
    """Forwards writes to the real stdout and to the active run's buffer."""

    def __init__(self, stdout):
        self._stdout = stdout

    def write(self, text):
        buffer = _current_buffer()
        if buffer is not None:
            buffer.write(text)
        return self._stdout.write(text)

    def writelines(self, lines):
        for line in lines:
            self.write(line)

    def flush(self):
        self._stdout.flush()

    def __getattr__(self, name):
        return getattr(self._stdout, name)


class _RunCaptureHandler(logging.Handler):
    """Copies each record into the buffer of the thread that emitted it."""

    def emit(self, record):
        buffer = _current_buffer()
        if buffer is None:
            return
        buffer.write(self.format(record) + "\n")


def configure_logging():
    """Point the root logger at the real stdout, plus the per-run capture.

    Idempotent: the workflow modules each call it at import time, but only the
    first call installs handlers.
    """
    global _configured
    if _configured:
        return

    # Drop only handlers a previous call installed, so re-configuring stays
    # idempotent without stripping handlers someone else owns (pytest's capture,
    # or a second copy of this module imported under another package name).
    for handler in logging.root.handlers[:]:
        if getattr(handler, "_streptocad_handler", False):
            logging.root.removeHandler(handler)

    formatter = logging.Formatter(LOG_FORMAT)

    # Goes to the container log; bound to the real stdout, never to a buffer.
    console_handler = logging.StreamHandler(sys.__stdout__ or sys.stdout)
    console_handler.setFormatter(formatter)

    capture_handler = _RunCaptureHandler()
    capture_handler.setFormatter(formatter)

    console_handler._streptocad_handler = True
    capture_handler._streptocad_handler = True

    logging.root.setLevel(logging.INFO)
    logging.root.addHandler(console_handler)
    logging.root.addHandler(capture_handler)

    _configured = True


@contextlib.contextmanager
def capture_output():
    """Capture this thread's prints and log records for the duration of a run."""
    global _active_runs, _real_stdout

    buffer = io.StringIO()
    previous = _current_buffer()
    _local.buffer = buffer

    with _stdout_lock:
        if _active_runs == 0:
            _real_stdout = sys.stdout
            sys.stdout = _TeeStdout(_real_stdout)
        _active_runs += 1

    try:
        yield buffer
    finally:
        with _stdout_lock:
            _active_runs -= 1
            if _active_runs == 0:
                sys.stdout = _real_stdout
                _real_stdout = None
        _local.buffer = previous


def capture_workflow_output(func):
    """Wrap a Dash callback so its run gets a fresh, isolated output buffer."""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        with capture_output():
            return func(*args, **kwargs)

    return wrapper


def workflow_log():
    """Everything printed or logged so far by the run active on this thread."""
    buffer = _current_buffer()
    return buffer.getvalue() if buffer is not None else ""
