"""Shared input handling and multiple-testing helpers for the qPCR scripts."""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

UNDETERMINED_TOKENS = {
    "",
    "-",
    "--",
    "n/a",
    "na",
    "nan",
    "none",
    "null",
    "no ct",
    "no cq",
    "undet",
    "undetermined",
    "not detected",
    "nd",
}

_DECIMAL_COMMA = re.compile(r"^[+-]?\d+,\d+$")


class InputError(ValueError):
    """Raised for input problems the analyst has to resolve."""


def read_text(path: str | Path) -> str:
    """Read UTF-8 or UTF-16 text, and explain any other encoding."""

    source = Path(path)
    if not source.is_file():
        raise InputError(f"input file not found: {source}")
    raw = source.read_bytes()
    try:
        if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
            text = raw.decode("utf-16")
        else:
            text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputError(f"{source} is not UTF-8; re-save it as CSV UTF-8") from exc
    if not text.strip():
        raise InputError(f"input file is empty: {source}")
    return text


def safe_console() -> None:
    """Keep non-Latin sample names from crashing output on legacy Windows code pages."""

    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")


def parse_number(raw: str, *, line: int, column: str) -> float:
    """Parse a required number, accepting a decimal comma."""

    text = str(raw).strip()
    if _DECIMAL_COMMA.match(text):
        text = text.replace(",", ".")
    try:
        value = float(text)
    except ValueError as exc:
        raise InputError(f"line {line}: {column} {raw!r} is not a number") from exc
    if not math.isfinite(value):
        raise InputError(f"line {line}: {column} {raw!r} is not a finite number")
    return value


def parse_ct(raw: str | None, *, line: int) -> float | None:
    """Return a Ct value, or None for an undetermined well.

    A Ct of zero or below is not a measurement. Several instruments write 0 for
    "no amplification", so it is treated as undetermined rather than as the
    highest expression on the plate. Unknown text is an error: guessing here
    would silently change a result.
    """

    text = "" if raw is None else str(raw).strip().lower()
    if text in UNDETERMINED_TOKENS:
        return None
    value = parse_number(text, line=line, column="ct")
    return None if value <= 0 else value


def read_table(path: str | Path, required: tuple[str, ...]) -> list[dict[str, Any]]:
    """Read a delimited table into row dicts with normalized header names."""

    source = Path(path)
    text = read_text(source)
    first_line = text.splitlines()[0]
    delimiter = max(",;\t", key=first_line.count)
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    if reader.fieldnames is None:
        raise InputError(f"input file has no header: {source}")
    reader.fieldnames = [name.strip().lower() for name in reader.fieldnames]
    missing = [name for name in required if name not in reader.fieldnames]
    if missing:
        raise InputError(
            f"missing column(s) {', '.join(missing)}; found {', '.join(reader.fieldnames)}"
        )

    rows = []
    for line, row in enumerate(reader, start=2):
        if None in row:
            raise InputError(f"line {line}: more fields than the header")
        if not any((value or "").strip() for value in row.values()):
            continue
        cleaned = {key: (value or "").strip() for key, value in row.items()}
        cleaned["_line"] = line
        rows.append(cleaned)
    if not rows:
        raise InputError(f"input file has no data rows: {source}")
    return rows


def benjamini_hochberg(p_values: list[float | None]) -> list[float | None]:
    """Benjamini-Hochberg adjusted p-values; None entries stay None."""

    indexed = [(p, i) for i, p in enumerate(p_values) if p is not None]
    adjusted: list[float | None] = [None] * len(p_values)
    count = len(indexed)
    running = 1.0
    for rank, (p, index) in zip(
        range(count, 0, -1), sorted(indexed, key=lambda item: item[0], reverse=True)
    ):
        running = min(running, p * count / rank)
        adjusted[index] = running
    return adjusted


def write_json(path: str | Path, payload: dict[str, Any]) -> None:
    """Write a result as UTF-8 JSON."""

    Path(path).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def rounded(value: float | None, digits: int = 4) -> float | None:
    """Round for reporting, passing None through."""

    return None if value is None else round(value, digits)
