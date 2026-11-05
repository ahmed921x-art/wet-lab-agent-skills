#!/usr/bin/env python3
"""Convert an instrument qPCR export into the long table qpcr_relative.py reads.

Finds the results table inside a QuantStudio, Bio-Rad CFX, or similar export
(CSV, TSV, or .xlsx), skipping metadata lines, then keeps sample, target, Ct, and
task. Omitted wells and standards are dropped and counted. Groups come from a
sample sheet (sample, group, and optionally pair) or from a regular expression
applied to the sample name.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Any

from _qpcr import InputError, read_table, read_text, safe_console

SAMPLE_NAMES = ("sample name", "sample", "samplename", "sample id", "name")
TARGET_NAMES = ("target name", "target", "gene", "gene name", "detector", "detector name", "assay", "assay name")
CT_NAMES = ("ct", "cт", "cq", "cp", "crt", "c(t)")
TASK_NAMES = ("task", "content", "well type", "type")
OMIT_NAMES = ("omit", "omitted", "excluded")
SPLIT = re.compile(r"[,;\t]")


def normalize_task(raw: str) -> str:
    """Map instrument task labels to unknown, ntc, nrt, or standard."""

    text = raw.strip().lower()
    if not text or text.startswith(("unkn", "unknown", "sample")):
        return "unknown"
    if text.startswith(("ntc", "neg ctrl", "negative", "no template")):
        return "ntc"
    if text.startswith(("nrt", "-rt", "no rt", "no-rt")):
        return "nrt"
    if text.startswith(("std", "standard")):
        return "standard"
    raise InputError(f"task label {raw!r} is not recognized; map it in a sample sheet")


def load_grid(path: Path, sheet: str | None) -> list[list[str]]:
    """Return every row of the export as a list of cell strings."""

    if path.suffix.lower() == ".xls":
        raise InputError("legacy .xls is not read; re-export or save as .xlsx or CSV")
    if path.suffix.lower() in {".xlsx", ".xlsm"}:
        try:
            from openpyxl import load_workbook
        except ImportError as exc:
            raise InputError("reading .xlsx needs openpyxl; install it or export CSV") from exc
        if not path.is_file():
            raise InputError(f"input file not found: {path}")
        book = load_workbook(path, read_only=True, data_only=True)
        names = [sheet] if sheet else book.sheetnames
        for name in names:
            if name not in book.sheetnames:
                raise InputError(f"sheet {name!r} not found; sheets are {book.sheetnames}")
            rows = [
                ["" if cell is None else str(cell).strip() for cell in row]
                for row in book[name].iter_rows(values_only=True)
            ]
            if find_header(rows) is not None:
                return rows
        raise InputError(f"no results table found in sheets {names}")

    text = read_text(path)
    lines = text.splitlines()
    delimiters = {}
    for line in lines:
        for delimiter in ",;\t":
            delimiters[delimiter] = delimiters.get(delimiter, 0) + line.count(delimiter)
    delimiter = max(delimiters, key=delimiters.get) if delimiters else ","
    return [[cell.strip() for cell in row] for row in csv.reader(lines, delimiter=delimiter)]


def locate(header: list[str], names: tuple[str, ...]) -> int | None:
    """First column whose label matches, in order of preference."""

    lowered = [cell.strip().lower() for cell in header]
    for name in names:
        if name in lowered:
            return lowered.index(name)
    return None


def find_header(rows: list[list[str]]) -> int | None:
    """Index of the first row naming sample, target, and Ct columns."""

    for index, row in enumerate(rows):
        if all(locate(row, names) is not None for names in (SAMPLE_NAMES, TARGET_NAMES, CT_NAMES)):
            return index
    return None


def convert(
    rows: list[list[str]],
    groups: dict[str, dict[str, str]] | None,
    group_pattern: str | None,
) -> dict[str, Any]:
    """Extract wells from an export grid and attach group labels."""

    header_index = find_header(rows)
    if header_index is None:
        raise InputError("no header with sample, target, and Ct columns was found")
    header = rows[header_index]
    sample_col = locate(header, SAMPLE_NAMES)
    target_col = locate(header, TARGET_NAMES)
    ct_col = locate(header, CT_NAMES)
    task_col = locate(header, TASK_NAMES)
    omit_col = locate(header, OMIT_NAMES)
    pattern = re.compile(group_pattern) if group_pattern else None

    wells: list[dict[str, str]] = []
    counts = {"omitted": 0, "standard": 0, "blank_name": 0}
    ungrouped: set[str] = set()
    for row in rows[header_index + 1 :]:
        if not any(row) or row[0].startswith("["):
            if wells:
                break
            continue
        if len(row) <= max(sample_col, target_col, ct_col):
            continue
        if omit_col is not None and len(row) > omit_col and row[omit_col].lower() in {"true", "yes", "1"}:
            counts["omitted"] += 1
            continue
        task = normalize_task(row[task_col] if task_col is not None and len(row) > task_col else "")
        if task == "standard":
            counts["standard"] += 1
            continue
        sample = row[sample_col]
        target = row[target_col]
        if not target or (not sample and task == "unknown"):
            counts["blank_name"] += 1
            continue
        if not sample:
            sample = task.upper()
        if sample.upper() == "NTC":
            task = "ntc"
        group, pair = "", ""
        if task == "unknown":
            if groups is not None:
                entry = groups.get(sample)
                if entry is None:
                    ungrouped.add(sample)
                else:
                    group, pair = entry["group"], entry.get("pair", "")
            elif pattern is not None:
                match = pattern.search(sample)
                if match is None or not match.groups():
                    ungrouped.add(sample)
                else:
                    group = match.group(1)
            else:
                ungrouped.add(sample)
        wells.append(
            {"sample": sample, "group": group, "target": target, "ct": row[ct_col], "task": task, "pair": pair}
        )
    if not wells:
        raise InputError("the results table has no usable wells")
    if ungrouped:
        raise InputError(
            f"no group for {len(ungrouped)} sample(s): {', '.join(sorted(ungrouped)[:8])}; "
            "give --sample-sheet or --group-pattern"
        )
    return {
        "wells": wells,
        "omitted_wells": counts["omitted"],
        "standard_wells": counts["standard"],
        "unnamed_wells": counts["blank_name"],
        "samples": len({w["sample"] for w in wells if w["task"] == "unknown"}),
        "targets": sorted({w["target"] for w in wells}),
    }


def read_sample_sheet(path: str) -> dict[str, dict[str, str]]:
    """Read a sample sheet keyed by sample name."""

    return {row["sample"]: row for row in read_table(path, ("sample", "group"))}


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="instrument export: CSV, TSV, or .xlsx")
    parser.add_argument("--output", required=True, help="long CSV for qpcr_relative.py")
    parser.add_argument("--sample-sheet", help="CSV with sample, group, and optional pair")
    parser.add_argument("--group-pattern", help=r"regex with one capture group, e.g. ^(\w+)-\d+$")
    parser.add_argument("--sheet", help="worksheet name for .xlsx")
    args = parser.parse_args(argv)

    try:
        groups = read_sample_sheet(args.sample_sheet) if args.sample_sheet else None
        result = convert(load_grid(Path(args.input), args.sheet), groups, args.group_pattern)
    except (InputError, re.error) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample", "group", "target", "ct", "task", "pair"])
        writer.writeheader()
        writer.writerows(result["wells"])
    print(f"[OK] {len(result['wells'])} wells, {result['samples']} samples, targets: {', '.join(result['targets'])}")
    if result["omitted_wells"]:
        print(f"[WARN] {result['omitted_wells']} well(s) were marked omitted on the instrument and dropped")
    if result["standard_wells"]:
        print(f"[OK] {result['standard_wells']} standard well(s) left out")
    if result["unnamed_wells"]:
        print(f"[WARN] {result['unnamed_wells']} well(s) had no sample or target name and were dropped")
    print(f"[OK] wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
