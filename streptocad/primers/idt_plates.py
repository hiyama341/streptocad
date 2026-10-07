"""Format oligos for IDT's plate upload.

IDT accepts two shapes of bulk oligo order. The tube format, produced by
``create_idt_order_dataframe``, carries a concentration and purification column
per oligo. The plate format carries neither -- those are chosen once for the
whole plate in IDT's web form -- and instead names the well each oligo goes in:

    Well Position | Name | Sequence

Wells are filled row-major (A1..A12, B1..B12, ... H12) and empty wells are
simply omitted. An order larger than one plate becomes several worksheets in a
single workbook, which is how IDT expects multi-plate uploads.
"""

import io
import re
from typing import Iterable, List, Sequence

import pandas as pd

# IDT's plate upload expects exactly these three headers, in this order.
IDT_PLATE_COLUMNS = ["Well Position", "Name", "Sequence"]

DEFAULT_PLATE_ROWS = 8
DEFAULT_PLATE_COLS = 12

# IDT does not publish a character whitelist for oligo names. Letters, digits,
# underscore, dot, plus and hyphen are what their own sample templates use, so
# anything else is folded to an underscore rather than risking a rejected order.
_UNSAFE_NAME_CHARS = re.compile(r"[^A-Za-z0-9_.+-]+")


def well_positions(
    rows: int = DEFAULT_PLATE_ROWS, cols: int = DEFAULT_PLATE_COLS
) -> List[str]:
    """List one plate's well labels in fill order.

    Row-major, matching how IDT's sample templates are laid out: A1 through A12,
    then B1 through B12, and so on.

    Parameters
    ----------
    rows : int
        Number of plate rows, lettered from A. Up to 26.
    cols : int
        Number of plate columns, numbered from 1.

    Returns
    -------
    List[str]
        ``["A1", "A2", ..., "H12"]`` for a default 96-well plate.
    """
    if not 1 <= rows <= 26:
        raise ValueError(f"rows must be between 1 and 26, got {rows}")
    if cols < 1:
        raise ValueError(f"cols must be at least 1, got {cols}")

    return [
        f"{chr(ord('A') + row)}{col}"
        for row in range(rows)
        for col in range(1, cols + 1)
    ]


def sanitize_idt_oligo_name(name: str) -> str:
    """Fold an oligo name to the characters IDT's templates use.

    StreptoCAD builds names from locus tags and plasmid ids, which can carry
    ``#``, parentheses or spaces picked up from a GenBank annotation. Those are
    replaced with a single underscore each so the upload is not rejected for a
    reason the user cannot see.

    Parameters
    ----------
    name : str
        Oligo name as StreptoCAD generated it.

    Returns
    -------
    str
        The name with unsafe runs collapsed to ``_`` and no leading or
        trailing underscore.
    """
    safe = _UNSAFE_NAME_CHARS.sub("_", str(name))
    # Collapse runs so a separator already next to a replaced character does not
    # double up: "sgRNA_#1" becomes "sgRNA_1", not "sgRNA__1".
    safe = re.sub(r"_{2,}", "_", safe).strip("_")
    return safe or "unnamed_oligo"


def assign_plate_wells(
    names: Sequence[str],
    sequences: Sequence[str],
    rows: int = DEFAULT_PLATE_ROWS,
    cols: int = DEFAULT_PLATE_COLS,
    sanitize_names: bool = True,
) -> List[pd.DataFrame]:
    """Lay oligos out across as many plates as they need.

    Order is preserved: the first oligo lands in A1 of the first plate and
    overflow starts a new plate at A1 again.

    Parameters
    ----------
    names : Sequence[str]
        Oligo names, in the order they should fill the plate.
    sequences : Sequence[str]
        Oligo sequences, parallel to ``names``.
    rows, cols : int
        Plate geometry. Defaults give a 96-well plate.
    sanitize_names : bool
        Fold names to IDT-safe characters. Leave on unless you have already
        sanitised them.

    Returns
    -------
    List[pd.DataFrame]
        One DataFrame per plate, each with the ``IDT_PLATE_COLUMNS`` headers.
        An empty input gives an empty list, not one empty plate.
    """
    if len(names) != len(sequences):
        raise ValueError(
            f"names and sequences must be the same length, "
            f"got {len(names)} and {len(sequences)}"
        )

    wells = well_positions(rows, cols)
    capacity = len(wells)

    plates = []
    for start in range(0, len(names), capacity):
        chunk_names = names[start : start + capacity]
        chunk_seqs = sequences[start : start + capacity]
        plates.append(
            pd.DataFrame(
                {
                    IDT_PLATE_COLUMNS[0]: wells[: len(chunk_names)],
                    IDT_PLATE_COLUMNS[1]: [
                        sanitize_idt_oligo_name(n) if sanitize_names else str(n)
                        for n in chunk_names
                    ],
                    IDT_PLATE_COLUMNS[2]: [str(s).upper() for s in chunk_seqs],
                }
            )
        )

    return plates


def idt_order_df_to_idt_plates(
    idt_df: pd.DataFrame,
    rows: int = DEFAULT_PLATE_ROWS,
    cols: int = DEFAULT_PLATE_COLS,
) -> List[pd.DataFrame]:
    """Re-shape a tube-format IDT order sheet into plate worksheets.

    Row order is taken as given, which is what makes this the only layout
    function needed: ``create_idt_order_dataframe`` already emits each
    template's forward primer followed by its reverse, so the pair arrives in
    neighbouring wells without this function knowing anything about primer
    pairs. Orders built from ssDNA or sgRNA oligos, which have no pairing to
    preserve, go through unchanged.

    Taking the plate from the same frame as the tube sheet also means the two
    files cannot disagree about which oligos are in the order.

    Parameters
    ----------
    idt_df : pd.DataFrame
        A DataFrame carrying ``Name`` and ``Sequence`` columns, as produced by
        ``create_idt_order_dataframe`` or ``primers_to_IDT``.
    rows, cols : int
        Plate geometry. Defaults give a 96-well plate.

    Returns
    -------
    List[pd.DataFrame]
        One DataFrame per plate, each with the ``IDT_PLATE_COLUMNS`` headers.
    """
    missing = [c for c in ("Name", "Sequence") if c not in idt_df.columns]
    if missing:
        raise KeyError(
            f"idt_df is missing the columns {missing}; got {list(idt_df.columns)}"
        )

    return assign_plate_wells(
        idt_df["Name"].tolist(), idt_df["Sequence"].tolist(), rows=rows, cols=cols
    )


def idt_plates_to_xlsx_bytes(
    plates: Iterable[pd.DataFrame], sheet_prefix: str = "Plate"
) -> bytes:
    """Write plate worksheets to an in-memory .xlsx for IDT's plate upload.

    IDT reads one plate per worksheet, so a multi-plate order is a single
    workbook with several tabs rather than several files.

    Parameters
    ----------
    plates : Iterable[pd.DataFrame]
        Plate DataFrames, as returned by the ``*_to_idt_plates`` helpers.
    sheet_prefix : str
        Worksheet names become ``{sheet_prefix}1``, ``{sheet_prefix}2``, ...

    Returns
    -------
    bytes
        The .xlsx payload, ready to write into a zip or hand to a browser.
    """
    plates = list(plates)
    if not plates:
        raise ValueError("no plates to write; an IDT upload needs at least one plate")

    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        for index, plate in enumerate(plates, start=1):
            plate.to_excel(
                writer, sheet_name=f"{sheet_prefix}{index}", index=False
            )

    return buffer.getvalue()
