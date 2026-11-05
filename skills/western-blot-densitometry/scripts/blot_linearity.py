#!/usr/bin/env python3
"""Check whether band density tracks protein load in a dilution series.

Input columns: target, load, signal. Fits log2(signal) on log2(load). A slope near
1 means signal doubles when load doubles. The script finds the longest run of
loads, starting from the lowest, that stays within the slope limits, which gives
the highest load that can be run for that antibody.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

from blot_quant import InputError, read_text, safe_console, to_number

SLOPE_LOW = 0.8
SLOPE_HIGH = 1.2
MIN_R_SQUARED = 0.95
MIN_LEVELS = 3


def fit(points: list[tuple[float, float]]) -> dict[str, float]:
    """Least-squares slope and R^2 of log2 signal on log2 load."""

    xs = [math.log2(load) for load, _ in points]
    ys = [math.log2(signal) for _, signal in points]
    mean_x, mean_y = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - mean_x) ** 2 for x in xs)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    syy = sum((y - mean_y) ** 2 for y in ys)
    slope = sxy / sxx
    return {"slope": slope, "r_squared": sxy**2 / (sxx * syy) if syy > 0 else 0.0}


def linear_range(points: list[tuple[float, float]]) -> float | None:
    """Highest load whose run from the lowest level remains log-log linear."""

    levels = sorted({load for load, _ in points})
    best = None
    for count in range(MIN_LEVELS, len(levels) + 1):
        allowed = set(levels[:count])
        result = fit([point for point in points if point[0] in allowed])
        if SLOPE_LOW <= result["slope"] <= SLOPE_HIGH and result["r_squared"] >= MIN_R_SQUARED:
            best = levels[count - 1]
    return best


def analyze(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Assess each target."""

    by_target: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in rows:
        by_target[row["target"]].append((row["load"], row["signal"]))
    results = []
    for target, points in sorted(by_target.items()):
        if len({load for load, _ in points}) < MIN_LEVELS:
            results.append({"target": target, "error": f"needs at least {MIN_LEVELS} load levels"})
            continue
        overall = fit(points)
        top = linear_range(points)
        highest = max(load for load, _ in points)
        flags = []
        if not SLOPE_LOW <= overall["slope"] <= SLOPE_HIGH:
            flags.append(f"slope {overall['slope']:.2f} over the full series is outside {SLOPE_LOW}-{SLOPE_HIGH}")
        if top is None:
            flags.append("no linear range found: lower the load or change the antibody dilution")
        elif top < highest:
            flags.append(f"signal stops tracking load above {top:g}; do not load more than this")
        results.append(
            {
                "target": target,
                "slope": round(overall["slope"], 3),
                "r_squared": round(overall["r_squared"], 4),
                "linear_up_to": top,
                "flags": flags,
            }
        )
    return {"targets": results}


def read_series(path: str | Path) -> list[dict[str, Any]]:
    """Read target, load, and signal."""

    import csv

    text = read_text(path)
    delimiter = max(",;\t", key=text.splitlines()[0].count)
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    reader.fieldnames = [name.strip().lower() for name in reader.fieldnames or []]
    missing = [name for name in ("target", "load", "signal") if name not in reader.fieldnames]
    if missing:
        raise InputError(f"missing column(s): {', '.join(missing)}")
    rows = []
    for line, raw in enumerate(reader, start=2):
        row = {key: (value or "").strip() for key, value in raw.items() if key}
        if not any(row.values()):
            continue
        load = to_number(row["load"], line, "load")
        signal = to_number(row["signal"], line, "signal")
        if load <= 0 or signal <= 0:
            raise InputError(f"line {line}: load and signal must be positive")
        rows.append({"target": row["target"], "load": load, "signal": signal})
    if not rows:
        raise InputError("no rows found")
    return rows


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="CSV with target, load, signal")
    parser.add_argument("--output-json", help="write full results here")
    args = parser.parse_args(argv)
    try:
        result = analyze(read_series(args.input))
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2
    for row in result["targets"]:
        if "error" in row:
            print(f"[WARN] {row['target']}: {row['error']}")
            continue
        marker = "[WARN]" if row["flags"] else "[OK]"
        print(f"{marker} {row['target']}: slope {row['slope']}, R^2 {row['r_squared']}, linear up to {row['linear_up_to']}")
        for flag in row["flags"]:
            print(f"       -> {flag}")
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[OK] wrote {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
