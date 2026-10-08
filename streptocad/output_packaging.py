"""Package a workflow's results into the zip a user downloads.

Every web-app workflow funnels its results through :class:`OutputPackage`, so
the layout of the download is decided here and nowhere else. Results are
declared by *role* rather than by filename::

    OutputPackage(
        workflow_id="w1",
        slug="overexpression",
        outputs=[
            {"role": "plasmid", "content": assembled_plasmids},
            {"role": "primer.pcr", "content": primer_df},
            {"role": "primer.order_idt", "content": idt_df},
        ],
        parameters={"genes": [...]},
    ).to_zip_bytes(timestamp)

A role fixes the artifact's shelf, its filename and the line describing it in
the generated README, all from :data:`ROLES`. Callbacks therefore cannot drift
from one another on any of the three, which is what happened while each one
spelled its own ``00_``/``01_`` prefixes.

The shelves read in the order you use them at the bench:

====================  =========================================================
``1_plasmids/``       maps to open in SnapGene or Benchling, plus their index
``2_primers/``        the oligos to order, in tube and plate format
``3_sgrnas/``         guide RNA candidates and the ones that were built
``4_analysis/``       QC tables, the assembly overview, the run record
``5_protocols/``      bench instructions, as HTML and Markdown
``6_inputs/``         exactly what went in, so the run can be repeated
====================  =========================================================

A shelf a workflow does not use is omitted rather than shipped empty, and the
README says why it is absent.
"""

import copy
import io
import json
import logging
import os
import platform
import posixpath
import re
import sys
from typing import Any, Dict, Iterable, List, Optional, Sequence

import nbformat
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from nbconvert import HTMLExporter

SHELF_PLASMIDS = "1_plasmids"
SHELF_PRIMERS = "2_primers"
SHELF_SGRNAS = "3_sgrnas"
SHELF_ANALYSIS = "4_analysis"
SHELF_PROTOCOLS = "5_protocols"
SHELF_INPUTS = "6_inputs"

SHELF_ORDER = [
    SHELF_PLASMIDS,
    SHELF_PRIMERS,
    SHELF_SGRNAS,
    SHELF_ANALYSIS,
    SHELF_PROTOCOLS,
    SHELF_INPUTS,
]

SHELF_BLURBS = {
    SHELF_PLASMIDS: "Plasmid maps. Open these in SnapGene or Benchling.",
    SHELF_PRIMERS: "Oligos to order. Upload the IDT file to idtdna.com.",
    SHELF_SGRNAS: "Guide RNA candidates, and the ones that ended up in the plasmids.",
    SHELF_ANALYSIS: "Quality control, the assembly overview, and the record of this run.",
    SHELF_PROTOCOLS: "Bench instructions. Open the .html files; the .md files are the same text.",
    SHELF_INPUTS: "Exactly what you submitted, so this run can be repeated.",
}

# Why a shelf can be legitimately absent. Keyed by (workflow_id, shelf).
SHELF_NOT_USED = {
    ("w1", SHELF_SGRNAS): "W1 designs overexpression plasmids and uses no guide RNAs.",
}

# role -> (shelf, filename, kind, blurb)
#
# ``kind`` selects the writer:
#   table    pandas DataFrame            -> .csv
#   records  SeqRecord or list of them   -> .gb
#   plate    list of DataFrames          -> .xlsx, one worksheet per plate
#   text     str                         -> written as-is
#   json     anything json-serialisable  -> .json
#
# ``filename`` of None means the writer names the files itself, which only the
# exploded plasmid role does.
ROLES: Dict[str, tuple] = {
    "plasmid": (
        SHELF_PLASMIDS,
        None,
        "records",
        "One GenBank file per designed plasmid, numbered to match plasmid_index.csv.",
    ),
    "plasmid.all": (
        SHELF_PLASMIDS,
        "00_all_plasmids.gb",
        "records",
        "Every designed plasmid in one file, to import in a single step.",
    ),
    "plasmid.index": (
        SHELF_PLASMIDS,
        "plasmid_index.csv",
        "table",
        "What each plasmid file is: name, size, origin and the edit it carries.",
    ),
    "primer.order_idt": (
        SHELF_PRIMERS,
        "oligo_order_idt.csv",
        "table",
        "Tube-format oligo order. Paste into IDT's bulk entry form.",
    ),
    "primer.order_idt_plate": (
        SHELF_PRIMERS,
        "oligo_order_idt_plate.xlsx",
        "plate",
        "Plate-format oligo order. Upload directly to IDT; one worksheet per 96-well plate.",
    ),
    "primer.pcr": (
        SHELF_PRIMERS,
        "pcr_primers.csv",
        "table",
        "Amplification primers with melting and annealing temperatures for the thermocycler.",
    ),
    "primer.check": (
        SHELF_PRIMERS,
        "check_primers.csv",
        "table",
        "Verification primers for colony or junction PCR.",
    ),
    "sgrna.all": (
        SHELF_SGRNAS,
        "sgrna_all_candidates.csv",
        "table",
        "Every guide RNA found in the target genes.",
    ),
    "sgrna.selected": (
        SHELF_SGRNAS,
        "sgrna_selected.csv",
        "table",
        "The guide RNAs that are in the plasmids in 1_plasmids/.",
    ),
    "sgrna.base_edit_predictions": (
        SHELF_SGRNAS,
        "sgrna_base_edit_predictions.csv",
        "table",
        "Predicted base-editing outcome for each candidate guide.",
    ),
    "analysis.primer_qc": (
        SHELF_ANALYSIS,
        "primer_qc_hairpins.csv",
        "table",
        "Hairpin and dimer check for the designed primers.",
    ),
    "analysis.golden_gate_overhangs": (
        SHELF_ANALYSIS,
        "golden_gate_overhangs.csv",
        "table",
        "Golden Gate overhangs used to order the multiplexed assembly.",
    ),
    "analysis.assembly_overview": (
        SHELF_ANALYSIS,
        "assembly_overview.txt",
        "text",
        "Counts and per-plasmid notes from the assembly simulation.",
    ),
    "analysis.run_log": (
        SHELF_ANALYSIS,
        "run.log",
        "text",
        "Log of this run, for reporting a problem.",
    ),
    "analysis.environment": (
        SHELF_ANALYSIS,
        "environment.json",
        "json",
        "Package versions this run used, so it can be reproduced exactly.",
    ),
    "bench.order": (
        SHELF_PROTOCOLS,
        "bench_order.csv",
        "table",
        "The order to build the constructs in.",
    ),
    "input.sequences": (
        SHELF_INPUTS,
        "input_sequences.gb",
        "records",
        "The sequences you uploaded, in one file.",
    ),
    "input.genome": (
        SHELF_INPUTS,
        "input_genome.gb",
        "records",
        "The genome you uploaded.",
    ),
    "input.plasmid": (
        SHELF_INPUTS,
        "input_plasmid.gb",
        "records",
        "The plasmid backbone you uploaded, before any edit.",
    ),
    "input.parameters": (
        SHELF_INPUTS,
        "run_parameters.json",
        "json",
        "Every setting you chose, as submitted.",
    ),
}

# Protocol key -> (source filename under the protocols directory, shipped stem).
#
# The shipped name comes from this table and not from the source basename, so
# the long-standing "protcol" misspelling on disk does not reach users.
PROTOCOLS: Dict[str, tuple] = {
    "conjugation": ("conjugation_protcol.md", "conjugation_protocol"),
    "overexpression": ("overexpression_protocol.md", "overexpression_protocol"),
    "crispr_single_target": (
        "single_target_crispr_plasmid_protcol.md",
        "crispr_single_target_cloning_protocol",
    ),
    "crispr_multi_target": (
        "multi_target_crispr_plasmid_protcol.md",
        "crispr_multi_target_cloning_protocol",
    ),
    "cas3_single_target": (
        "cas3_single_target_crispr_plasmid_protocol.md",
        "cas3_single_target_cloning_protocol",
    ),
    "troubleshooting": ("trouble_shooting_tips.md", "troubleshooting_tips"),
}

WORKFLOW_SLUGS = {
    "w1": "overexpression",
    "w2": "crispr-best-single",
    "w3": "crispr-best-multiplex",
    "w4": "crispri",
    "w5": "cas9-deletion",
    "w6": "cas3-deletion",
}

#: Protocol markdown lives at the repository root, beside the ``streptocad``
#: package. Resolved from ``__file__`` rather than the process working
#: directory, which is what let a stale second copy be served depending on
#: where the app was started from.
DEFAULT_PROTOCOLS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "protocols"
)

_UNSAFE_PATH_CHARS = re.compile(r"[^A-Za-z0-9_.+-]+")


def sanitize_path_component(value: str, fallback: str = "unnamed") -> str:
    """Reduce a string to something safe to use as a file or directory name.

    Record ids and plasmid names reach the archive straight from GenBank
    annotations, where they can carry ``/``, ``#``, parentheses, spaces or
    ``..``. Any of those either breaks extraction or, in the case of a path
    separator, silently creates a directory nobody asked for.

    Parameters
    ----------
    value : str
        The raw string, typically ``record.name`` or ``record.id``.
    fallback : str
        Returned when sanitising leaves nothing usable.

    Returns
    -------
    str
        A name made only of letters, digits, ``_``, ``.``, ``+`` and ``-``.
    """
    safe = _UNSAFE_PATH_CHARS.sub("_", str(value))
    safe = re.sub(r"_{2,}", "_", safe).strip("._-")
    return safe or fallback


def run_id(workflow_id: str, timestamp: str, slug: Optional[str] = None) -> str:
    """Build the archive's root directory name, which is also its filename.

    Parameters
    ----------
    workflow_id : str
        ``"w1"`` through ``"w6"``.
    timestamp : str
        A ``YYYYMMDD-HHMMSSZ`` stamp. Passed in rather than read from the clock
        so a caller can tie the archive to its own log lines.
    slug : str, optional
        Overrides the slug from :data:`WORKFLOW_SLUGS`.

    Returns
    -------
    str
        For example ``streptocad_w1-overexpression_20261007-121503Z``. Contains
        no colons, so it extracts on Windows -- unlike the ISO-8601 stamps the
        callbacks used to build.
    """
    slug = slug or WORKFLOW_SLUGS.get(workflow_id, workflow_id)
    return sanitize_path_component(f"streptocad_{workflow_id}-{slug}_{timestamp}")


def environment_record() -> Dict[str, str]:
    """Capture the package versions this run used."""
    versions = {"python": sys.version.split()[0], "platform": platform.platform()}
    for name in ("streptocad", "pandas", "Bio", "pydna", "teemi", "openpyxl"):
        try:
            module = __import__(name)
            versions[name] = str(getattr(module, "__version__", "unknown"))
        except Exception:  # pragma: no cover - a missing optional package
            versions[name] = "not installed"
    return versions


def _records_to_genbank(records: Sequence[SeqRecord]) -> str:
    """Serialise records to GenBank, repairing what ``SeqIO`` would reject.

    The caller's records are left untouched: a record may be declared under
    more than one role, and the LOCUS repair below must not leak back into the
    objects the workflow is still holding.
    """
    prepared = []
    for record in records:
        prepared_record = copy.copy(record)
        # SeqIO raises on whitespace in the LOCUS name, which aborts the whole
        # download. A sanitised, length-capped name is better than no archive.
        prepared_record.name = sanitize_path_component(
            record.name or record.id or "plasmid"
        )[:16]
        prepared_record.annotations = dict(record.annotations)
        prepared_record.annotations.setdefault("molecule_type", "DNA")
        prepared.append(prepared_record)

    with io.StringIO() as buffer:
        SeqIO.write(prepared, buffer, "genbank")
        return buffer.getvalue()


class RunLogCapture:
    """Collect this run's log records, and only this run's.

    Each callback module used to build its own ``io.StringIO`` and install it
    through ``logging.basicConfig`` at import time. Because ``basicConfig`` is a
    no-op once the root logger has handlers -- and because some modules cleared
    the root handlers first -- only one module's buffer was ever attached, and
    which one depended on import order. Every other workflow then shipped an
    empty log.

    Capturing per run instead of per module removes the ordering question, and
    means a run's log holds that run's records rather than everything since the
    server started.

    Use it around the body of a callback::

        capture = RunLogCapture().start()
        try:
            ...
        finally:
            capture.stop()

        ...or as a context manager, which is equivalent::

        with RunLogCapture() as capture:
            ...
    """

    FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    def __init__(self, level: int = logging.INFO):
        self.level = level
        self.stream = io.StringIO()
        self._handler: Optional[logging.Handler] = None
        self._previous_level: Optional[int] = None

    def start(self) -> "RunLogCapture":
        """Attach to the root logger. Returns self, so it can be chained."""
        self._handler = logging.StreamHandler(self.stream)
        self._handler.setFormatter(logging.Formatter(self.FORMAT))
        self._handler.setLevel(self.level)

        root = logging.getLogger()
        self._previous_level = root.level
        # A root logger left at WARNING would drop the INFO lines the
        # workflows emit, so lower it for the duration and put it back after.
        if root.level > self.level or root.level == logging.NOTSET:
            root.setLevel(self.level)
        root.addHandler(self._handler)
        return self

    def stop(self) -> str:
        """Detach from the root logger and return what was captured."""
        root = logging.getLogger()
        if self._handler is not None:
            root.removeHandler(self._handler)
            self._handler.close()
            self._handler = None
        if self._previous_level is not None:
            root.setLevel(self._previous_level)
            self._previous_level = None
        return self.stream.getvalue()

    def getvalue(self) -> str:
        """What has been captured so far, without detaching."""
        return self.stream.getvalue()

    def __enter__(self) -> "RunLogCapture":
        return self.start()

    def __exit__(self, *exc_info) -> bool:
        self.stop()
        return False


def _markdown_to_html(markdown: str) -> str:
    """Render protocol markdown the way the app has always rendered it."""
    notebook = nbformat.v4.new_notebook()
    notebook.cells.append(nbformat.v4.new_markdown_cell(markdown))
    html, _ = HTMLExporter().from_notebook_node(notebook)
    return html


class OutputPackage:
    """Collect a workflow's results and write them out as a zip archive.

    Parameters
    ----------
    workflow_id : str
        ``"w1"`` through ``"w6"``. Selects the slug in the archive name and
        decides which absent shelves are expected.
    outputs : Iterable[dict]
        One ``{"role": ..., "content": ...}`` mapping per artifact. Roles must
        appear in :data:`ROLES`; an unknown one raises, so a typo fails here
        rather than quietly dropping a file in the wrong shelf.
    inputs : Iterable[dict], optional
        The uploaded files, declared the same way with ``input.*`` roles.
    parameters : dict, optional
        Everything the user chose, written to ``6_inputs/run_parameters.json``.
    protocols : Sequence[str], optional
        Keys from :data:`PROTOCOLS`. Each ships as both ``.html`` and ``.md``.
    slug : str, optional
        Overrides the workflow's slug in the archive name.
    protocols_dir : str, optional
        Where to read protocol markdown from. Defaults to the repository's
        ``protocols/`` directory, resolved from this module's location.
    """

    def __init__(
        self,
        workflow_id: str,
        outputs: Iterable[Dict[str, Any]],
        inputs: Iterable[Dict[str, Any]] = (),
        parameters: Optional[Dict[str, Any]] = None,
        protocols: Sequence[str] = (),
        slug: Optional[str] = None,
        protocols_dir: Optional[str] = None,
    ):
        self.workflow_id = workflow_id
        self.outputs = list(outputs)
        self.inputs = list(inputs)
        self.parameters = parameters or {}
        self.protocol_keys = list(protocols)
        self.slug = slug
        self.protocols_dir = protocols_dir or DEFAULT_PROTOCOLS_DIR

        self.warnings: List[str] = []
        #: Shipped path -> role, filled in during :meth:`to_zip_bytes`.
        self.written: Dict[str, str] = {}

        unknown = [
            entry.get("role")
            for entry in self.outputs + self.inputs
            if entry.get("role") not in ROLES
        ]
        if unknown:
            raise KeyError(
                f"unknown output role(s) {unknown}; "
                f"add them to streptocad.output_packaging.ROLES"
            )

        unknown_protocols = [k for k in self.protocol_keys if k not in PROTOCOLS]
        if unknown_protocols:
            raise KeyError(
                f"unknown protocol key(s) {unknown_protocols}; "
                f"known keys are {sorted(PROTOCOLS)}"
            )

    # -- writers ---------------------------------------------------------

    def _write_table(self, zip_file, path: str, content: Any) -> None:
        if not isinstance(content, pd.DataFrame):
            raise TypeError(
                f"{path} is declared as a table, so its content must be a "
                f"DataFrame, not {type(content).__name__}"
            )
        if content.empty:
            self.warnings.append(f"{path} is empty -- no rows were produced.")
        zip_file.writestr(path, content.to_csv(index=False))

    def _write_plate(self, zip_file, path: str, content: Any, role: str) -> None:
        # Imported here so the package's primer modules stay importable
        # without openpyxl present.
        from .primers.idt_plates import idt_plates_to_xlsx_bytes

        plates = list(content) if not isinstance(content, pd.DataFrame) else [content]
        if not plates:
            self.warnings.append(f"{path} was not written -- there were no oligos.")
            return
        zip_file.writestr(path, idt_plates_to_xlsx_bytes(plates))
        self.written[path] = role

    def _write_records(
        self, zip_file, shelf: str, filename: Optional[str], content, role: str
    ):
        records = content if isinstance(content, (list, tuple)) else [content]
        records = [r for r in records if r is not None]

        if not records:
            # An all-failed assembly used to ship a clean-looking archive with
            # no plasmids in it and no error anywhere.
            note = posixpath.join(shelf, "NO_PLASMIDS_ASSEMBLED.txt")
            zip_file.writestr(
                note,
                "No plasmids were assembled for this run.\n"
                "Check 4_analysis/run.log and the parameters in "
                "6_inputs/run_parameters.json.\n",
            )
            self.written[note] = "plasmid"
            self.warnings.append(
                "No plasmids were assembled. See 1_plasmids/NO_PLASMIDS_ASSEMBLED.txt."
            )
            return

        if filename is not None:
            path = posixpath.join(shelf, filename)
            zip_file.writestr(path, _records_to_genbank(records))
            self.written[path] = role
            return

        # The exploded plasmid role. A leading, zero-padded index sorts
        # correctly and joins to plasmid_index.csv row n -- the old trailing
        # "_0".."_19" suffix sorted 0,1,10,11,... and ordered nothing.
        width = max(2, len(str(len(records))))
        for index, record in enumerate(records, start=1):
            stem = sanitize_path_component(
                record.name or record.id or f"plasmid_{index}", f"plasmid_{index}"
            )
            path = posixpath.join(shelf, f"{str(index).zfill(width)}_{stem}.gb")
            if path in self.written:
                path = posixpath.join(
                    shelf, f"{str(index).zfill(width)}_{stem}_{index}.gb"
                )
            zip_file.writestr(path, _records_to_genbank([record]))
            self.written[path] = "plasmid"

    def _write_text(self, zip_file, path: str, content: Any) -> None:
        zip_file.writestr(path, "" if content is None else str(content))

    def _write_json(self, zip_file, path: str, content: Any) -> None:
        zip_file.writestr(path, json.dumps(content, indent=2, default=str))

    def _write_entry(self, zip_file, entry: Dict[str, Any]) -> None:
        role = entry["role"]
        shelf, filename, kind, _ = ROLES[role]
        content = entry.get("content")

        if kind == "records":
            self._write_records(zip_file, shelf, filename, content, role)
            return

        path = posixpath.join(shelf, filename)
        if kind == "table":
            self._write_table(zip_file, path, content)
        elif kind == "plate":
            self._write_plate(zip_file, path, content, role)
            return
        elif kind == "text":
            self._write_text(zip_file, path, content)
        elif kind == "json":
            self._write_json(zip_file, path, content)
        else:  # pragma: no cover - guarded by the ROLES table
            raise ValueError(f"role {role!r} has an unknown kind {kind!r}")

        self.written[path] = role

    def _write_protocols(self, zip_file) -> None:
        for key in self.protocol_keys:
            source_name, stem = PROTOCOLS[key]
            source = os.path.join(self.protocols_dir, source_name)
            try:
                with open(source, "r", encoding="utf-8") as handle:
                    markdown = handle.read()
            except OSError as error:
                # A missing protocol is worth a warning, never a failed
                # download -- the designs are the valuable part.
                self.warnings.append(
                    f"Protocol {key!r} could not be read ({error.strerror}); "
                    f"it is missing from 5_protocols/."
                )
                continue

            md_path = posixpath.join(SHELF_PROTOCOLS, f"{stem}.md")
            html_path = posixpath.join(SHELF_PROTOCOLS, f"{stem}.html")
            zip_file.writestr(md_path, markdown)
            zip_file.writestr(html_path, _markdown_to_html(markdown))
            self.written[md_path] = f"protocol.{key}"
            self.written[html_path] = f"protocol.{key}"

    def _check_plate_covers_the_order(self) -> None:
        """Warn when the plate file holds fewer oligos than the tube sheet.

        The two files are the same order in two formats, so they must list the
        same oligos. A workflow whose tube sheet is a ``pd.concat`` of several
        frames can easily build its plate from just one of them, and the result
        is a plate that looks right and orders an incomplete set -- which is a
        failed experiment, not an inconvenience. Counting them here catches
        that for every workflow at once, including ones not written yet.
        """
        by_role = {entry.get("role"): entry.get("content") for entry in self.outputs}
        tube = by_role.get("primer.order_idt")
        plate = by_role.get("primer.order_idt_plate")
        if tube is None or plate is None:
            return

        tube_count = len(tube) if hasattr(tube, "__len__") else 0
        plates = list(plate) if not isinstance(plate, pd.DataFrame) else [plate]
        plate_count = sum(len(p) for p in plates)

        if plate_count < tube_count:
            self.warnings.append(
                f"The plate order in {SHELF_PRIMERS}/oligo_order_idt_plate.xlsx lists "
                f"{plate_count} oligos but the tube order lists {tube_count}. "
                f"Order from oligo_order_idt.csv until this is fixed."
            )

    # -- readme ----------------------------------------------------------

    def _readme_markdown(self, root: str) -> str:
        """Build the README from the same table that placed the files."""
        lines = [
            "# StreptoCAD results",
            "",
            f"Run `{root}`.",
            "",
        ]

        if self.warnings:
            lines += ["## Warnings", ""]
            lines += [f"- {w}" for w in self.warnings]
            lines.append("")

        lines += [
            "## Do this next",
            "",
            f"1. Upload `{SHELF_PRIMERS}/oligo_order_idt_plate.xlsx` to idtdna.com,"
            " or paste `oligo_order_idt.csv` into their bulk entry form.",
            f"2. Open `{SHELF_PLASMIDS}/` in SnapGene or Benchling.",
            f"3. Follow the protocols in `{SHELF_PROTOCOLS}/`.",
            "",
            "## What is in this folder",
            "",
        ]

        shipped_by_shelf: Dict[str, List[str]] = {}
        for path in sorted(self.written):
            shelf = path.split("/")[0]
            shipped_by_shelf.setdefault(shelf, []).append(path)

        for shelf in SHELF_ORDER:
            if shelf not in shipped_by_shelf:
                reason = SHELF_NOT_USED.get((self.workflow_id, shelf))
                if reason:
                    lines += [f"### {shelf}/ -- not used", "", reason, ""]
                continue

            lines += [f"### {shelf}/", "", SHELF_BLURBS[shelf], ""]
            lines += ["| file | what it is |", "| --- | --- |"]

            plasmid_files = [
                p
                for p in shipped_by_shelf[shelf]
                if self.written[p] == "plasmid" and p.endswith(".gb")
            ]
            for path in shipped_by_shelf[shelf]:
                name = path.split("/", 1)[1]
                role = self.written[path]
                if role == "plasmid" and path.endswith(".gb"):
                    continue
                blurb = (
                    ROLES[role][3]
                    if role in ROLES
                    else "Bench instructions for this workflow."
                )
                lines.append(f"| `{name}` | {blurb} |")
            if plasmid_files:
                lines.append(
                    f"| `01_*.gb` .. `{plasmid_files[-1].split('/', 1)[1]}` "
                    f"| {len(plasmid_files)} designed plasmid(s), "
                    f"numbered to match `plasmid_index.csv`. |"
                )
            lines.append("")

        return "\n".join(lines)

    # -- entry point -----------------------------------------------------

    def to_zip_bytes(self, timestamp: str) -> bytes:
        """Write the archive and return its bytes.

        Parameters
        ----------
        timestamp : str
            A ``YYYYMMDD-HHMMSSZ`` stamp for the archive's root directory.

        Returns
        -------
        bytes
            The zip payload, ready to base64-encode for a download link.
        """
        import zipfile

        self.warnings = []
        self.written = {}

        root = run_id(self.workflow_id, timestamp, self.slug)
        inner = io.BytesIO()

        with zipfile.ZipFile(inner, "w", zipfile.ZIP_DEFLATED) as zip_file:
            for entry in self.outputs:
                self._write_entry(zip_file, entry)
            for entry in self.inputs:
                self._write_entry(zip_file, entry)

            self._check_plate_covers_the_order()

            self._write_entry(
                zip_file, {"role": "input.parameters", "content": self.parameters}
            )
            self._write_entry(
                zip_file,
                {"role": "analysis.environment", "content": environment_record()},
            )
            self._write_protocols(zip_file)

            readme = self._readme_markdown(root)
            zip_file.writestr("00_READ_ME_FIRST.md", readme)
            zip_file.writestr("00_READ_ME_FIRST.html", _markdown_to_html(readme))

        # Everything is written relative to the shelves, then re-rooted under a
        # single directory so the archive expands into one folder.
        outer = io.BytesIO()
        with zipfile.ZipFile(inner, "r") as source, zipfile.ZipFile(
            outer, "w", zipfile.ZIP_DEFLATED
        ) as target:
            for name in source.namelist():
                target.writestr(posixpath.join(root, name), source.read(name))

        return outer.getvalue()

    @property
    def archive_name(self) -> str:
        """The filename to offer the browser, without the ``.zip``."""
        return run_id(self.workflow_id, "", self.slug).rstrip("_")
