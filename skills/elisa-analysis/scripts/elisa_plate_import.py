#!/usr/bin/env python3
"""Join a plate reader's OD grid with a plate layout into the table elisa_fit.py reads.

Both files hold a grid labelled A-H across and 1-12 down the side and top; reader
preamble above the grid is skipped. Layout cells use these forms:
  std:1000          standard at nominal 1000
  blank             blank well
  qc:QC1:200        QC named QC1, nominal 200 (optional :dilution)
  P01               sample P01, undiluted
  P01:10            sample P01 at 1:10
Empty layout cells are unused wells. A reference wavelength grid can be subtracted.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from elisa_fit import InputError, read_text, safe_console, to_float

ROW_LABELS = "ABCDEFGH"


def read_grid(path: str | Path) -> dict[str, str]:
    """Return {well: cell} for a labelled plate grid."""

    text = read_text(path)
    lines = text.splitlines()
    delimiter = max(",;\t", key=lambda item: sum(line.count(item) for line in lines))
    rows = [[cell.strip() for cell in row] for row in csv.reader(lines, delimiter=delimiter)]

    columns: dict[int, int] | None = None
    start = None
    for index, row in enumerate(rows):
        positions = {}
        for column, cell in enumerate(row):
            if cell.isdigit() and 1 <= int(cell) <= 12:
                positions[int(cell)] = column
        if len(positions) >= 2 and sorted(positions) == list(range(1, len(positions) + 1)):
            columns = positions
            start = index + 1
            break
    if columns is None or start is None:
        raise InputError(f"{path}: no column header 1, 2, 3 ... was found")

    label_column = min(columns.values()) - 1
    if label_column < 0:
        raise InputError(f"{path}: the grid has no row labels A-H")
    grid: dict[str, str] = {}
    for row in rows[start:]:
        if len(row) <= label_column:
            continue
        label = row[label_column].upper()
        if label not in ROW_LABELS or not label:
            if grid:
                break
            continue
        for number, column in columns.items():
            grid[f"{label}{number}"] = row[column] if column < len(row) else ""
    if not grid:
        raise InputError(f"{path}: no rows labelled A-H were found")
    return grid


def parse_cell(cell: str, well: str) -> dict[str, str] | None:
    """Interpret one layout cell."""

    text = cell.strip()
    if not text:
        return None
    parts = [part.strip() for part in text.split(":")]
    head = parts[0].lower()
    if head == "blank" and len(parts) == 1:
        return {"kind": "blank", "id": "blank", "conc": "", "dilution": ""}
    if head == "std":
        if len(parts) != 2 or not parts[1]:
            raise InputError(f"{well}: {text!r} should be std:CONCENTRATION")
        return {"kind": "standard", "id": f"S{parts[1]}", "conc": parts[1], "dilution": ""}
    if head == "qc":
        if len(parts) not in {3, 4}:
            raise InputError(f"{well}: {text!r} should be qc:NAME:NOMINAL or qc:NAME:NOMINAL:DILUTION")
        return {"kind": "qc", "id": parts[1], "conc": parts[2], "dilution": parts[3] if len(parts) == 4 else "1"}
    if len(parts) > 2:
        raise InputError(f"{well}: {text!r} should be NAME or NAME:DILUTION")
    return {"kind": "sample", "id": parts[0], "conc": "", "dilution": parts[1] if len(parts) == 2 else "1"}


def combine(
    od_grid: dict[str, str],
    layout: dict[str, str],
    reference: dict[str, str] | None = None,
) -> list[dict[str, str]]:
    """Build one row per used well, checking that the two grids agree."""

    wells = []
    for well, cell in layout.items():
        parsed = parse_cell(cell, well)
        if parsed is None:
            continue
        raw = od_grid.get(well, "")
        if not raw:
            raise InputError(f"{well}: layout says {cell!r} but the OD grid is empty there")
        value = to_float(raw, line=0, column=f"OD at {well}")
        if reference is not None:
            corrected = reference.get(well, "")
            if not corrected:
                raise InputError(f"{well}: no value in the reference wavelength grid")
            value -= to_float(corrected, line=0, column=f"reference OD at {well}")
        wells.append({**parsed, "well": well, "od": f"{value:.4f}"})
    if not wells:
        raise InputError("the layout has no used wells")
    return wells


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--od", required=True, help="reader grid, e.g. 450 nm")
    parser.add_argument("--layout", required=True, help="plate layout grid")
    parser.add_argument("--reference-od", help="reference wavelength grid, e.g. 570 or 630 nm")
    parser.add_argument("--output", required=True, help="long CSV for elisa_fit.py")
    args = parser.parse_args(argv)

    try:
        reference = read_grid(args.reference_od) if args.reference_od else None
        wells = combine(read_grid(args.od), read_grid(args.layout), reference)
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["kind", "id", "conc", "dilution", "od", "well"])
        writer.writeheader()
        writer.writerows(wells)
    kinds = {}
    for well in wells:
        kinds[well["kind"]] = kinds.get(well["kind"], 0) + 1
    print(f"[OK] {len(wells)} wells: {kinds}")
    print(f"[OK] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
