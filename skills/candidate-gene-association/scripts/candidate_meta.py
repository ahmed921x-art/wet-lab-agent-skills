#!/usr/bin/env python3
"""Meta-analysis of published odds ratios for one variant and one genetic model.

Input is one row per study: study, then either a, b, c, d (effect-carrying and
other counts in cases, then in controls) or or, ci_low, ci_high. An optional
hwe_p column gives the Hardy-Weinberg p in controls, used for a sensitivity
analysis. Reports inverse-variance fixed-effect and DerSimonian-Laird random-effects
estimates, Q, I^2, tau^2, leave-one-out results, and a sensitivity set.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

from _genetics import InputError, read_text, safe_console

Z_95 = 1.959963985
HWE_FLAG = 0.05


def chi2_upper(statistic: float, df: int) -> float:
    """Chi-square upper tail for any df, by the regularized incomplete gamma."""

    if statistic <= 0:
        return 1.0
    shape, x = df / 2.0, statistic / 2.0
    log_prefix = shape * math.log(x) - x - math.lgamma(shape)
    if x < shape + 1:
        term = total = 1.0 / shape
        count = shape
        for _ in range(10000):
            count += 1
            term *= x / count
            total += term
            if abs(term) < abs(total) * 1e-15:
                break
        return max(0.0, 1.0 - total * math.exp(log_prefix))
    tiny = 1e-300
    b = x + 1 - shape
    c = 1 / tiny
    d = 1 / b
    h = d
    for index in range(1, 10000):
        a = -index * (index - shape)
        b += 2
        d = a * d + b
        d = tiny if abs(d) < tiny else d
        c = b + a / c
        c = tiny if abs(c) < tiny else c
        d = 1 / d
        delta = d * c
        h *= delta
        if abs(delta - 1) < 1e-15:
            break
    return min(1.0, h * math.exp(log_prefix))


def read_studies(path: str | Path) -> list[dict[str, Any]]:
    """Read studies and convert each one to a log OR with its variance."""

    text = read_text(path)
    delimiter = max(",;\t", key=text.splitlines()[0].count)
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    reader.fieldnames = [name.strip().lower() for name in reader.fieldnames or []]
    if "study" not in reader.fieldnames:
        raise InputError("a study column is required")
    studies = []
    for line, raw in enumerate(reader, start=2):
        if None in raw:
            raise InputError(f"line {line}: more fields than the header")
        row = {key: (value or "").strip() for key, value in raw.items()}
        if not any(row.values()):
            continue
        studies.append(to_effect(row, line))
    if not studies:
        raise InputError("no studies found")
    return studies


def number(row: dict[str, str], key: str, line: int) -> float:
    """Parse a required numeric field."""

    try:
        return float(row[key])
    except ValueError as exc:
        raise InputError(f"line {line}: {key} {row[key]!r} is not a number") from exc


def to_effect(row: dict[str, str], line: int) -> dict[str, Any]:
    """Convert counts or a published interval to log OR and variance."""

    study: dict[str, Any] = {"study": row["study"], "corrected": False}
    if row.get("hwe_p"):
        study["hwe_p"] = number(row, "hwe_p", line)
    if all(row.get(key) for key in "abcd"):
        counts = [number(row, key, line) for key in "abcd"]
        if min(counts) < 0:
            raise InputError(f"line {line}: counts cannot be negative")
        if (counts[0] == 0 and counts[2] == 0) or (counts[1] == 0 and counts[3] == 0):
            raise InputError(f"line {line}: a double-zero study has no information about the OR")
        if 0 in counts:
            counts = [value + 0.5 for value in counts]
            study["corrected"] = True
        a, b, c, d = counts
        study["log_or"] = math.log(a * d / (b * c))
        study["variance"] = 1 / a + 1 / b + 1 / c + 1 / d
    elif all(row.get(key) for key in ("or", "ci_low", "ci_high")):
        value, low, high = (number(row, key, line) for key in ("or", "ci_low", "ci_high"))
        if not 0 < low < value < high:
            raise InputError(f"line {line}: expected 0 < ci_low < or < ci_high")
        study["log_or"] = math.log(value)
        study["variance"] = ((math.log(high) - math.log(low)) / (2 * Z_95)) ** 2
    else:
        raise InputError(f"line {line}: give a, b, c, d or or, ci_low, ci_high")
    return study


def pool(studies: list[dict[str, Any]]) -> dict[str, Any]:
    """Fixed-effect and DerSimonian-Laird random-effects estimates."""

    weights = [1 / study["variance"] for study in studies]
    total = sum(weights)
    fixed = sum(w * study["log_or"] for w, study in zip(weights, studies)) / total
    q = sum(w * (study["log_or"] - fixed) ** 2 for w, study in zip(weights, studies))
    df = len(studies) - 1
    scale = total - sum(w * w for w in weights) / total
    tau2 = max(0.0, (q - df) / scale) if df > 0 and scale > 0 else 0.0
    random_weights = [1 / (study["variance"] + tau2) for study in studies]
    random_total = sum(random_weights)
    random = sum(w * study["log_or"] for w, study in zip(random_weights, studies)) / random_total

    def summary(center: float, weight_sum: float) -> dict[str, float]:
        se = math.sqrt(1 / weight_sum)
        z = center / se
        return {
            "or": math.exp(center),
            "ci_low": math.exp(center - Z_95 * se),
            "ci_high": math.exp(center + Z_95 * se),
            "p_value": math.erfc(abs(z) / math.sqrt(2)),
        }

    return {
        "k": len(studies),
        "fixed": summary(fixed, total),
        "random": summary(random, random_total),
        "q": q,
        "q_p": chi2_upper(q, df) if df > 0 else None,
        "i2_percent": 100 * max(0.0, (q - df) / q) if q > 0 else 0.0,
        "tau2": tau2,
    }


def analyze(studies: list[dict[str, Any]]) -> dict[str, Any]:
    """Pool studies and add sensitivity analyses and warnings."""

    if len(studies) < 2:
        raise InputError("a meta-analysis needs at least two studies")
    result = pool(studies)
    warnings = []
    if len(studies) < 5:
        warnings.append(f"only {len(studies)} studies: tau^2 and I^2 are poorly estimated")
    if result["i2_percent"] >= 50:
        warnings.append(
            f"I^2 is {result['i2_percent']:.0f}%: studies disagree, so explore ancestry and "
            "design before quoting one pooled number; prefer the random-effects estimate"
        )
    corrected = [study["study"] for study in studies if study["corrected"]]
    if corrected:
        warnings.append(f"0.5 continuity correction used for: {', '.join(corrected)}")
    if len(studies) < 10:
        warnings.append("under 10 studies: funnel plot tests for publication bias are unreliable")

    leave_one_out = []
    if len(studies) > 2:
        for index, study in enumerate(studies):
            others = studies[:index] + studies[index + 1 :]
            pooled = pool(others)["random"]
            leave_one_out.append({"omitted": study["study"], **pooled})

    sensitivity = None
    deviating = [s["study"] for s in studies if s.get("hwe_p") is not None and s["hwe_p"] < HWE_FLAG]
    if deviating:
        kept = [s for s in studies if s["study"] not in deviating]
        warnings.append(f"controls deviate from Hardy-Weinberg in: {', '.join(deviating)}")
        if len(kept) >= 2:
            sensitivity = {"excluded": deviating, **pool(kept)}
    per_study = [
        {
            "study": study["study"],
            "or": math.exp(study["log_or"]),
            "ci_low": math.exp(study["log_or"] - Z_95 * math.sqrt(study["variance"])),
            "ci_high": math.exp(study["log_or"] + Z_95 * math.sqrt(study["variance"])),
            "weight_random_percent": None,
        }
        for study in studies
    ]
    weights = [1 / (study["variance"] + result["tau2"]) for study in studies]
    for row, weight in zip(per_study, weights):
        row["weight_random_percent"] = round(100 * weight / sum(weights), 1)
    return {
        **result,
        "studies": per_study,
        "leave_one_out": leave_one_out,
        "hwe_sensitivity": sensitivity,
        "warnings": warnings,
    }


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="CSV with one row per study")
    parser.add_argument("--output-json", help="write full results here")
    args = parser.parse_args(argv)

    try:
        result = analyze(read_studies(args.input))
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    for label in ("fixed", "random"):
        row = result[label]
        print(
            f"[OK] {label}: OR {row['or']:.3f} ({row['ci_low']:.3f}-{row['ci_high']:.3f}), p {row['p_value']:.4g}"
        )
    print(f"  -> k={result['k']}, Q={result['q']:.2f}, I^2={result['i2_percent']:.0f}%, tau^2={result['tau2']:.4f}")
    if result["hwe_sensitivity"]:
        row = result["hwe_sensitivity"]["random"]
        print(f"  -> without HWE-deviating studies: OR {row['or']:.3f} ({row['ci_low']:.3f}-{row['ci_high']:.3f})")
    for warning in result["warnings"]:
        print(f"[WARN] {warning}")
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[OK] wrote {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
