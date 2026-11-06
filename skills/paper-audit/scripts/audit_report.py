#!/usr/bin/env python3
"""Turn an audit_claims.py ledger into a Markdown consistency report.

The report lists inconsistencies first, each with its page and exact quote, then
claims that could not be checked and why, then consistent claims, and closes with
what this audit cannot detect. It states that a failed check is something to
examine, not an accusation of misconduct.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from _audit import CONSISTENT, INCONSISTENT, NOT_CHECKABLE, InputError
from audit_claims import read_text, safe_console

CHECK_NAMES = {
    "p_value": "p-value from test statistic",
    "grim": "GRIM (mean of integer data)",
    "percent": "percentage against n",
    "sum": "parts against total",
    "ratio_ci": "ratio with confidence interval",
    "diff_ci": "difference with confidence interval",
}

FRAMING = (
    "A failed check marks an internal inconsistency to examine. It is **not** evidence of misconduct. "
    "Typing and transcription errors, different rounding, one-tailed or corrected tests (Welch, "
    "Greenhouse-Geisser), statistics from another model than assumed here, exclusions that change n, and "
    "mistakes in building the claim table can all produce a flag. Re-read each quoted passage, then raise "
    "remaining questions with the authors or editor in neutral terms."
)

CANNOT_CHECK = [
    "Whether the right test was used, or its assumptions hold.",
    "Whether the data exist or match the analysis; these checks use only the printed numbers.",
    "Numbers that were never entered in the claim table. Claims left out are not audited.",
    "Figures, unless their values were read off and entered as claims.",
    "Standard deviations (GRIMMER and SPRITE-type checks are not implemented).",
    "Multiple-comparison corrections, unless the reported p is the corrected one and the table says so.",
]


def escape(text: Any) -> str:
    return str(text).replace("|", "\\|").replace("\n", " ")


def describe(entry: dict[str, Any]) -> list[str]:
    lines = [
        f"### {escape(entry['id'])} (page {escape(entry['page'] or '?')}): {CHECK_NAMES.get(entry['check'], entry['check'])}",
        "",
        f"> {escape(entry['quote']) or '(no quote given)'}",
        "",
        f"- Finding: {entry['reason']}.",
    ]
    if entry.get("decision_error"):
        lines.append("- The significance decision differs between the reported and the recomputed values.")
    if entry.get("recomputed"):
        shown = "; ".join(f"{key} = {value}" for key, value in entry["recomputed"].items())
        lines.append(f"- Recomputed: {shown}")
    for note in entry.get("notes", []):
        lines.append(f"- Note: {note}")
    if entry.get("note"):
        lines.append(f"- Analyst note: {escape(entry['note'])}")
    lines.append("")
    return lines


def build_markdown(ledger: dict[str, Any], title: str = "Paper consistency audit (draft)") -> str:
    if not isinstance(ledger, dict) or "claims" not in ledger or "counts" not in ledger:
        raise InputError("this is not an audit_claims.py ledger")
    counts = ledger["counts"]
    claims = ledger["claims"]
    lines = [f"# {title}", "", "## How to read this report", "", FRAMING, "",
             "Every finding below quotes the paper and gives its page, so it can be checked against the source.",
             "", "## Summary", "",
             "| Status | Claims |", "| --- | --- |",
             f"| Inconsistent | {counts.get(INCONSISTENT, 0)} |",
             f"| Not checkable | {counts.get(NOT_CHECKABLE, 0)} |",
             f"| Consistent | {counts.get(CONSISTENT, 0)} |",
             "",
             f"Significance decisions that change: {ledger.get('decision_errors', 0)} "
             f"(alpha = {ledger.get('alpha', 0.05)}).", ""]
    for warning in ledger.get("warnings", []):
        lines.append(f"- Warning: {warning}")
    if ledger.get("warnings"):
        lines.append("")
    sections = [
        (INCONSISTENT, "Inconsistencies to examine", "No claim failed its check."),
        (NOT_CHECKABLE, "Not checkable", "Every claim in the table could be checked."),
    ]
    for status, heading, empty in sections:
        lines += [f"## {heading}", ""]
        chosen = [entry for entry in claims if entry["status"] == status]
        if not chosen:
            lines += [empty, ""]
        for entry in chosen:
            lines += describe(entry)
    lines += ["## Consistent", ""]
    consistent = [entry for entry in claims if entry["status"] == CONSISTENT]
    if consistent:
        lines += ["| Claim | Page | Check | Finding |", "| --- | --- | --- | --- |"]
        lines += [f"| {escape(e['id'])} | {escape(e['page'])} | {e['check']} | {escape(e['reason'])} |" for e in consistent]
    else:
        lines.append("No claim was confirmed consistent.")
    lines += ["", "## What this audit cannot detect", ""]
    lines += [f"- {item}" for item in CANNOT_CHECK]
    lines += ["", "A consistent result means only that the printed numbers fit together. It does not show that "
                  "the analysis or the data are right."]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--ledger", required=True, help="JSON written by audit_claims.py --output-json")
    parser.add_argument("--out", default="audit_report.md", help="Markdown file to write")
    parser.add_argument("--title", default="Paper consistency audit (draft)")
    args = parser.parse_args(argv)
    try:
        try:
            ledger = json.loads(read_text(args.ledger))
        except json.JSONDecodeError as exc:
            raise InputError(f"{args.ledger} is not valid JSON: {exc}") from exc
        text = build_markdown(ledger, args.title)
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2
    Path(args.out).write_text(text, encoding="utf-8")
    print(f"[OK] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
