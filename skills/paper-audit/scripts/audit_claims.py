#!/usr/bin/env python3
"""Recompute the numbers reported in a paper and write a claim ledger.

Input is a CSV or JSON claim table, one row per claim. Every row needs id, page,
quote (the exact text from the paper) and check, one of:

  p_value   test (t, F, chi2, z, r), statistic, df (df2 for F; n instead of df for r), p, tails
  grim      mean, n, items (default 1); only for means of integer responses
  percent   percent, n, count (optional)
  sum       parts (separated by ';'), total
  ratio_ci  estimate, ci_low, ci_high, level (default 95), p (optional), ci_method (optional)
  diff_ci   estimate, ci_low, ci_high, level, p, df (optional, t-based), ci_method

Each claim ends as consistent, inconsistent, or not checkable, with the reason and
the recomputed values. Rounding of every reported number is allowed for before a
claim is called inconsistent. An inconsistency is something to examine, not
evidence of misconduct. Standard library only; no network access.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from pathlib import Path
from typing import Any

from _audit import ALPHA_DEFAULT, CHECKS, CONSISTENT, INCONSISTENT, NOT_CHECKABLE, InputError, clean_text, result

COLUMNS = [
    "id", "page", "quote", "check", "test", "statistic", "df", "df2", "n", "p", "tails", "mean", "items",
    "percent", "count", "parts", "total", "estimate", "ci_low", "ci_high", "level", "ci_method", "note",
]
MAX_ROWS = 5000


def safe_console() -> None:
    """Keep non-Latin quotes from crashing output on legacy Windows code pages."""

    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")


def read_text(path: str | Path) -> str:
    source = Path(path)
    if not source.is_file():
        raise InputError(f"no such file: {source}")
    raw = source.read_bytes()
    if raw.startswith((b"\xff\xfe", b"\xfe\xff")):
        return raw.decode("utf-16")
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputError(f"{source} is not UTF-8 or UTF-16 text; re-save it as UTF-8") from exc


def read_claims(path: str | Path) -> list[dict[str, Any]]:
    """Read a CSV or JSON claim table into a list of dictionaries with string values."""

    text = read_text(path)
    if not text.strip():
        raise InputError(f"{path} is empty")
    if str(path).lower().endswith(".json") or text.lstrip().startswith(("[", "{")):
        try:
            data = json.loads(text)
        except json.JSONDecodeError as exc:
            raise InputError(f"{path} is not valid JSON: {exc}") from exc
        if isinstance(data, dict) and isinstance(data.get("claims"), list):
            data = data["claims"]
        if not isinstance(data, list) or not all(isinstance(row, dict) for row in data):
            raise InputError("JSON input must be a list of claim objects, or {\"claims\": [...]}")
        rows = [{str(key).strip().lower(): ("" if value is None else str(value)) for key, value in row.items()}
                for row in data]
    else:
        reader = csv.DictReader(io.StringIO(text))
        if not reader.fieldnames:
            raise InputError(f"{path} has no header row")
        reader.fieldnames = [name.strip().lower() for name in reader.fieldnames]
        rows = []
        for row in reader:
            if None in row:
                raise InputError(f"row {reader.line_num} has more cells than the header; quote any text with commas")
            rows.append({key: (value or "") for key, value in row.items()})
    if not rows:
        raise InputError(f"{path} has no claims")
    if len(rows) > MAX_ROWS:
        raise InputError(f"{len(rows)} claims is more than the limit of {MAX_ROWS}")
    for required in ("id", "page", "quote", "check"):
        if required not in rows[0]:
            raise InputError(f"missing column '{required}'")
    return rows


def audit(rows: list[dict[str, Any]], alpha: float = ALPHA_DEFAULT) -> dict[str, Any]:
    """Run every claim through its check and return the ledger."""

    if not 0 < alpha < 1:
        raise InputError("alpha must lie between 0 and 1")
    unknown = sorted({key for row in rows for key in row} - set(COLUMNS))
    warnings = [f"ignored unknown column '{name}'" for name in unknown]
    seen: set[str] = set()
    ledger = []
    for index, row in enumerate(rows, start=1):
        claim_id = clean_text(row.get("id")) or f"row{index}"
        if claim_id in seen:
            raise InputError(f"claim id '{claim_id}' appears more than once")
        seen.add(claim_id)
        check = clean_text(row.get("check")).lower()
        entry = {
            "id": claim_id,
            "page": clean_text(row.get("page")),
            "quote": clean_text(row.get("quote")),
            "check": check,
            "reported": {key: clean_text(value) for key, value in row.items()
                         if key in COLUMNS and key not in {"id", "page", "quote", "check", "note"} and clean_text(value)},
            "note": clean_text(row.get("note")),
        }
        if not entry["quote"] or not entry["page"]:
            outcome = result(NOT_CHECKABLE, "no quote or page, so the finding could not be traced to the paper; "
                                            "add the exact quote and its page")
        elif check not in CHECKS:
            outcome = result(NOT_CHECKABLE, f"check '{check}' is not one of {', '.join(sorted(CHECKS))}")
        else:
            try:
                outcome = CHECKS[check](row, alpha) if check == "p_value" else CHECKS[check](row)
            except InputError as exc:
                outcome = result(NOT_CHECKABLE, f"input problem: {exc}")
            except (ArithmeticError, ValueError) as exc:
                outcome = result(NOT_CHECKABLE, f"could not compute: {exc}")
        entry.update(outcome)
        ledger.append(entry)
    counts = {status: sum(1 for entry in ledger if entry["status"] == status)
              for status in (CONSISTENT, INCONSISTENT, NOT_CHECKABLE)}
    return {
        "tool": "paper-audit",
        "alpha": alpha,
        "counts": counts,
        "decision_errors": sum(1 for entry in ledger if entry["decision_error"]),
        "claims": ledger,
        "warnings": warnings,
    }


def write_csv(ledger: dict[str, Any], path: str | Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["id", "page", "check", "status", "decision_error", "reason", "recomputed", "notes", "quote"])
        for entry in ledger["claims"]:
            writer.writerow([
                entry["id"], entry["page"], entry["check"], entry["status"], "yes" if entry["decision_error"] else "",
                entry["reason"], json.dumps(entry["recomputed"], ensure_ascii=False), " | ".join(entry["notes"]),
                entry["quote"],
            ])


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="claim table, .csv or .json")
    parser.add_argument("--alpha", type=float, default=ALPHA_DEFAULT, help="significance level used for decision checks")
    parser.add_argument("--output-json", help="write the full ledger here")
    parser.add_argument("--output-csv", help="write a flat ledger here")
    args = parser.parse_args(argv)

    try:
        ledger = audit(read_claims(args.input), alpha=args.alpha)
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    for entry in ledger["claims"]:
        mark = {CONSISTENT: "OK  ", INCONSISTENT: "FLAG", NOT_CHECKABLE: "N/A "}[entry["status"]]
        print(f"[{mark}] {entry['id']} (p. {entry['page'] or '?'}, {entry['check']}): {entry['reason']}")
    for warning in ledger["warnings"]:
        print(f"[WARN] {warning}")
    counts = ledger["counts"]
    print(f"  -> {counts[CONSISTENT]} consistent, {counts[INCONSISTENT]} inconsistent, "
          f"{counts[NOT_CHECKABLE]} not checkable")
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(ledger, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[OK] wrote {args.output_json}")
    if args.output_csv:
        write_csv(ledger, args.output_csv)
        print(f"[OK] wrote {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
