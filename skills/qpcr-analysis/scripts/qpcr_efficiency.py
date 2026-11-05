#!/usr/bin/env python3
"""Estimate qPCR amplification efficiency from a dilution-series standard curve.

Input columns: target, quantity, ct. Quantity is any relative scale (copies, ng,
or 1, 0.1, 0.01 for a serial dilution). Undetermined wells are dropped and
counted. The JSON output can be passed to qpcr_relative.py --efficiency-json.
"""

from __future__ import annotations

import argparse
import math
import statistics
import sys
from collections import defaultdict
from typing import Any

from scipy import stats

from _qpcr import (
    InputError,
    parse_ct,
    parse_number,
    read_table,
    rounded,
    safe_console,
    write_json,
)

EFFICIENCY_LOW_PERCENT = 90.0
EFFICIENCY_HIGH_PERCENT = 110.0
MIN_R_SQUARED = 0.98
RECOMMENDED_LEVELS = 5
MIN_LEVELS = 3


def fit_curve(points: list[tuple[float, float]]) -> dict[str, float]:
    """Ordinary least squares of Ct on log10(quantity), with a slope interval."""

    xs = [math.log10(quantity) for quantity, _ in points]
    ys = [ct for _, ct in points]
    count = len(xs)
    mean_x = statistics.fmean(xs)
    mean_y = statistics.fmean(ys)
    sxx = sum((x - mean_x) ** 2 for x in xs)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    syy = sum((y - mean_y) ** 2 for y in ys)
    if sxx == 0:
        raise InputError("standard curve needs more than one quantity")
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    residual = sum((y - (intercept + slope * x)) ** 2 for x, y in zip(xs, ys))
    r_squared = 1.0 - residual / syy if syy > 0 else 0.0
    se_slope = math.sqrt(residual / (count - 2) / sxx) if count > 2 else math.nan
    critical = stats.t.ppf(0.975, count - 2) if count > 2 else math.nan
    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
        "slope_ci_low": slope - critical * se_slope,
        "slope_ci_high": slope + critical * se_slope,
    }


def amplification_factor(slope: float) -> float:
    """Convert a Ct-vs-log10 slope to the per-cycle amplification factor."""

    if slope >= 0:
        raise InputError(f"slope {slope:.3f} is not negative, so Ct does not fall with input")
    return 10 ** (-1.0 / slope)


def analyze(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Fit one curve per target and flag curves unfit for quantification."""

    by_target: dict[str, list[tuple[float, float]]] = defaultdict(list)
    dropped: dict[str, int] = defaultdict(int)
    for row in rows:
        target = row["target"]
        if not target:
            raise InputError(f"line {row['_line']}: target is empty")
        quantity = parse_number(row["quantity"], line=row["_line"], column="quantity")
        if quantity <= 0:
            raise InputError(f"line {row['_line']}: quantity must be positive")
        ct = parse_ct(row["ct"], line=row["_line"])
        if ct is None:
            dropped[target] += 1
        else:
            by_target[target].append((quantity, ct))

    results = []
    warnings = []
    for target in sorted(set(by_target) | set(dropped)):
        points = by_target.get(target, [])
        levels = len({quantity for quantity, _ in points})
        if levels < MIN_LEVELS:
            warnings.append(
                f"{target}: only {levels} quantity level(s) with signal; no curve fitted"
            )
            continue
        fit = fit_curve(points)
        factor = amplification_factor(fit["slope"])
        percent = (factor - 1.0) * 100.0
        factor_bounds = sorted(
            amplification_factor(bound)
            for bound in (fit["slope_ci_low"], fit["slope_ci_high"])
            if bound < 0
        )
        flags = []
        if not EFFICIENCY_LOW_PERCENT <= percent <= EFFICIENCY_HIGH_PERCENT:
            flags.append(
                f"efficiency {percent:.1f}% is outside "
                f"{EFFICIENCY_LOW_PERCENT:.0f}-{EFFICIENCY_HIGH_PERCENT:.0f}%"
            )
        if fit["r_squared"] < MIN_R_SQUARED:
            flags.append(f"R^2 {fit['r_squared']:.4f} is below {MIN_R_SQUARED}")
        if levels < RECOMMENDED_LEVELS:
            flags.append(f"{levels} levels; {RECOMMENDED_LEVELS} or more are recommended")
        if dropped.get(target):
            flags.append(f"{dropped[target]} undetermined well(s) dropped")
        results.append(
            {
                "target": target,
                "n_wells": len(points),
                "n_levels": levels,
                "slope": rounded(fit["slope"]),
                "intercept": rounded(fit["intercept"]),
                "r_squared": rounded(fit["r_squared"]),
                "amplification_factor": rounded(factor),
                "efficiency_percent": rounded(percent, 2),
                "amplification_factor_ci": [rounded(value) for value in factor_bounds]
                if len(factor_bounds) == 2
                else None,
                "dynamic_range_log10": rounded(
                    math.log10(max(q for q, _ in points) / min(q for q, _ in points)), 2
                ),
                "flags": flags,
            }
        )
    return {
        "curves": results,
        "warnings": warnings,
        "efficiency": {row["target"]: row["amplification_factor"] for row in results},
    }


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="CSV with target, quantity, ct")
    parser.add_argument("--output-json", help="write full results here")
    args = parser.parse_args(argv)

    try:
        result = analyze(read_table(args.input, ("target", "quantity", "ct")))
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    for curve in result["curves"]:
        marker = "[WARN]" if curve["flags"] else "[OK]"
        print(
            f"{marker} {curve['target']}: slope {curve['slope']}, "
            f"R^2 {curve['r_squared']}, efficiency {curve['efficiency_percent']}% "
            f"(factor {curve['amplification_factor']})"
        )
        for flag in curve["flags"]:
            print(f"       -> {flag}")
    for warning in result["warnings"]:
        print(f"[WARN] {warning}")
    if args.output_json:
        write_json(args.output_json, result)
        print(f"[OK] wrote {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
