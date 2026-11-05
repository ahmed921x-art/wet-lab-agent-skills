#!/usr/bin/env python3
"""Normalize Western blot densitometry and compare groups across blots.

Input is one row per band: blot, sample, group, target, signal. Optional columns
are background (subtracted from signal), saturated (yes or true marks a band with
clipped pixels), and pair. The loading control may be a housekeeping protein or
total-protein stain, given as a target name.

Within each blot: band / loading, then log2, then relative to that blot's control
lanes, because raw density is not comparable between membranes. Repeat blots of a
sample are averaged, so the sample is the unit of analysis.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

from scipy import stats

LOADING_CV_LIMIT = 30.0
LOADING_SHIFT_LIMIT = 0.5


class InputError(ValueError):
    """Raised for input problems the analyst has to resolve."""


def safe_console() -> None:
    """Keep non-Latin sample names from crashing output on legacy Windows code pages."""

    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")


def read_text(path: str | Path) -> str:
    """Read UTF-8 or UTF-16 text, and explain any other encoding."""

    source = Path(path)
    if not source.is_file():
        raise InputError(f"input file not found: {source}")
    raw = source.read_bytes()
    try:
        text = raw.decode("utf-16") if raw[:2] in (b"\xff\xfe", b"\xfe\xff") else raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputError(f"{source} is not UTF-8; re-save it as CSV UTF-8") from exc
    if not text.strip():
        raise InputError(f"input file is empty: {source}")
    return text


def to_number(raw: str, line: int, column: str) -> float:
    """Parse a number, accepting a decimal comma."""

    text = raw.strip()
    if text.count(",") == 1 and "." not in text:
        text = text.replace(",", ".")
    try:
        value = float(text)
    except ValueError as exc:
        raise InputError(f"line {line}: {column} {raw!r} is not a number") from exc
    if not math.isfinite(value):
        raise InputError(f"line {line}: {column} {raw!r} is not finite")
    return value


def read_bands(path: str | Path) -> list[dict[str, Any]]:
    """Read band rows from a delimited file."""

    text = read_text(path)
    delimiter = max(",;\t", key=text.splitlines()[0].count)
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    reader.fieldnames = [name.strip().lower() for name in reader.fieldnames or []]
    required = ("blot", "sample", "group", "target", "signal")
    missing = [name for name in required if name not in reader.fieldnames]
    if missing:
        raise InputError(f"missing column(s): {', '.join(missing)}")
    bands = []
    for line, raw in enumerate(reader, start=2):
        if None in raw:
            raise InputError(f"line {line}: more fields than the header")
        row = {key: (value or "").strip() for key, value in raw.items()}
        if not any(row.values()):
            continue
        if not all(row[name] for name in required):
            raise InputError(f"line {line}: blot, sample, group, target, and signal are all required")
        signal = to_number(row["signal"], line, "signal")
        background = to_number(row["background"], line, "background") if row.get("background") else 0.0
        bands.append(
            {
                "blot": row["blot"],
                "sample": row["sample"],
                "group": row["group"],
                "target": row["target"],
                "signal": signal - background,
                "saturated": row.get("saturated", "").lower() in {"yes", "y", "true", "1"},
                "pair": row.get("pair", ""),
            }
        )
    if not bands:
        raise InputError("no bands found")
    return bands


def benjamini_hochberg(p_values: list[float | None]) -> list[float | None]:
    """Benjamini-Hochberg adjusted p-values; None entries stay None."""

    indexed = sorted(
        ((p, i) for i, p in enumerate(p_values) if p is not None), key=lambda item: item[0], reverse=True
    )
    adjusted: list[float | None] = [None] * len(p_values)
    running = 1.0
    count = len(indexed)
    for rank, (p, index) in zip(range(count, 0, -1), indexed):
        running = min(running, p * count / rank)
        adjusted[index] = running
    return adjusted


def compare(group: list[float], control: list[float], pairs, alpha: float) -> dict[str, Any]:
    """Welch or paired t-test on log2 values, with an interval on the difference."""

    if pairs:
        differences = [a - b for a, b in pairs]
        diff = statistics.fmean(differences)
        count = len(differences)
        if count < 2:
            return {"design": "paired", "diff": diff, "p_value": None, "ci": None}
        se = statistics.stdev(differences) / math.sqrt(count)
        df = float(count - 1)
        design = "paired"
    else:
        design = "welch"
        diff = statistics.fmean(group) - statistics.fmean(control)
        if len(group) < 2 or len(control) < 2:
            return {"design": design, "diff": diff, "p_value": None, "ci": None}
        var_g, var_c = statistics.variance(group) / len(group), statistics.variance(control) / len(control)
        se = math.sqrt(var_g + var_c)
        if se == 0:
            return {"design": design, "diff": diff, "p_value": None, "ci": None}
        df = (var_g + var_c) ** 2 / (var_g**2 / (len(group) - 1) + var_c**2 / (len(control) - 1))
    if se == 0:
        return {"design": design, "diff": diff, "p_value": None, "ci": None}
    critical = float(stats.t.ppf(1 - alpha / 2, df))
    return {
        "design": design,
        "diff": diff,
        "p_value": float(2 * stats.t.sf(abs(diff / se), df)),
        "ci": (diff - critical * se, diff + critical * se),
    }


def analyze(
    bands: list[dict[str, Any]],
    loading: str,
    control: str,
    *,
    keep_saturated: bool = False,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Normalize within blots, then compare groups at the sample level."""

    warnings: list[str] = []
    sample_group: dict[str, str] = {}
    sample_pair: dict[str, str] = {}
    for band in bands:
        if sample_group.setdefault(band["sample"], band["group"]) != band["group"]:
            raise InputError(f"sample {band['sample']} appears in more than one group")
        sample_pair.setdefault(band["sample"], band["pair"])
    if control not in sample_group.values():
        raise InputError(f"control group {control!r} not found; groups are {sorted(set(sample_group.values()))}")
    if loading not in {band["target"] for band in bands}:
        raise InputError(f"loading control {loading!r} not found")

    table: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    excluded = []
    for band in bands:
        key = (band["blot"], band["sample"], band["target"])
        if band["saturated"] and not keep_saturated:
            excluded.append({"blot": band["blot"], "sample": band["sample"], "target": band["target"], "reason": "saturated"})
            continue
        if band["signal"] <= 0:
            excluded.append({"blot": band["blot"], "sample": band["sample"], "target": band["target"], "reason": "no signal above background"})
            continue
        table[key].append(band["signal"])
    if any(item["reason"] == "saturated" for item in excluded):
        warnings.append("saturated bands were excluded: a clipped band understates the true signal")
    if keep_saturated and any(band["saturated"] for band in bands):
        warnings.append("saturated bands were kept by request, so their ratios are unreliable")

    blots = sorted({key[0] for key in table})
    targets = sorted({key[2] for key in table} - {loading})
    if not targets:
        raise InputError("no target left after removing the loading control")

    loading_review = []
    for blot in blots:
        values = [statistics.fmean(v) for k, v in table.items() if k[0] == blot and k[2] == loading]
        if len(values) < 2:
            continue
        cv = 100 * statistics.stdev(values) / statistics.fmean(values)
        loading_review.append({"blot": blot, "n_lanes": len(values), "cv_percent": round(cv, 1)})
        if cv > LOADING_CV_LIMIT:
            warnings.append(f"{blot}: loading control CV is {cv:.0f}%; check for uneven loading or transfer")

    log_loading: dict[str, list[float]] = defaultdict(list)
    per_blot: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for (blot, sample, target), values in table.items():
        if target == loading:
            continue
        load = table.get((blot, sample, loading))
        if not load:
            excluded.append({"blot": blot, "sample": sample, "target": target, "reason": "no usable loading band"})
            continue
        per_blot[(blot, target)].append(
            {"sample": sample, "log_ratio": math.log2(statistics.fmean(values) / statistics.fmean(load))}
        )
    for (blot, sample, target), values in table.items():
        if target == loading:
            log_loading[sample_group[sample]].append(math.log2(statistics.fmean(values)))

    sample_values: dict[tuple[str, str], list[float]] = defaultdict(list)
    dropped_blots = set()
    for (blot, target), entries in per_blot.items():
        controls = [e["log_ratio"] for e in entries if sample_group[e["sample"]] == control]
        if not controls:
            dropped_blots.add((blot, target))
            continue
        baseline = statistics.fmean(controls)
        for entry in entries:
            sample_values[(entry["sample"], target)].append(entry["log_ratio"] - baseline)
    for blot, target in sorted(dropped_blots):
        warnings.append(f"{blot} / {target}: no control lane, so this blot cannot be normalized and was left out")

    if len(log_loading) > 1:
        means = {group: statistics.fmean(values) for group, values in log_loading.items()}
        shift = max(means.values()) - min(means.values())
        if shift > LOADING_SHIFT_LIMIT:
            warnings.append(
                f"loading control differs by {shift:.2f} log2 units between groups, so the treatment "
                "may change it, or loading was unequal; consider total-protein normalization"
            )

    samples_out = []
    collapsed: dict[str, dict[str, float]] = defaultdict(dict)
    for (sample, target), values in sorted(sample_values.items()):
        collapsed[target][sample] = statistics.fmean(values)
        samples_out.append(
            {
                "sample": sample,
                "group": sample_group[sample],
                "target": target,
                "log2_relative": round(collapsed[target][sample], 4),
                "n_blots": len(values),
            }
        )

    contrasts = []
    other_groups = sorted(set(sample_group.values()) - {control})
    for target in sorted(collapsed):
        values = collapsed[target]
        control_samples = [s for s in values if sample_group[s] == control]
        for group in other_groups:
            group_samples = [s for s in values if sample_group[s] == group]
            if not group_samples or not control_samples:
                warnings.append(f"{target}: {group} vs {control} has no usable samples")
                continue
            pairs = None
            names = group_samples + control_samples
            if all(sample_pair.get(s) for s in names):
                by_pair = {sample_pair[s]: s for s in control_samples}
                matched = [(values[s], values[by_pair[sample_pair[s]]]) for s in group_samples if sample_pair[s] in by_pair]
                if matched:
                    pairs = matched
            result = compare([values[s] for s in group_samples], [values[s] for s in control_samples], pairs, alpha)
            n_group = len(pairs) if pairs else len(group_samples)
            n_control = len(pairs) if pairs else len(control_samples)
            if min(n_group, n_control) < 3:
                warnings.append(f"{target}: {group} vs {control} has n={n_group} vs n={n_control}; treat as descriptive")
            ci = result["ci"]
            contrasts.append(
                {
                    "target": target,
                    "group": group,
                    "control": control,
                    "n_group": n_group,
                    "n_control": n_control,
                    "design": result["design"],
                    "log2_fold_change": round(result["diff"], 4),
                    "fold_change": round(2 ** result["diff"], 4),
                    "fold_ci_low": round(2 ** ci[0], 4) if ci else None,
                    "fold_ci_high": round(2 ** ci[1], 4) if ci else None,
                    "p_value": result["p_value"],
                }
            )
    for row, adjusted in zip(contrasts, benjamini_hochberg([row["p_value"] for row in contrasts])):
        row["p_adjusted_bh"] = adjusted

    single_blot = sum(1 for row in samples_out if row["n_blots"] == 1)
    return {
        "loading": loading,
        "control": control,
        "ci_level": 1 - alpha,
        "contrasts": contrasts,
        "samples": samples_out,
        "n_samples_single_blot": single_blot,
        "loading_review": loading_review,
        "excluded": excluded,
        "warnings": warnings,
    }


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="CSV with one row per band")
    parser.add_argument("--loading", required=True, help="loading control target name")
    parser.add_argument("--control", required=True, help="control group name")
    parser.add_argument("--keep-saturated", action="store_true", help="keep saturated bands; flagged in the output")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--output-json", help="write full results here")
    args = parser.parse_args(argv)

    try:
        result = analyze(
            read_bands(args.input), args.loading, args.control, keep_saturated=args.keep_saturated, alpha=args.alpha
        )
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    for row in result["contrasts"]:
        interval = f"CI {row['fold_ci_low']}-{row['fold_ci_high']}" if row["fold_ci_low"] is not None else "no CI"
        p_text = "p n/a" if row["p_value"] is None else f"p {row['p_value']:.4g}, BH {row['p_adjusted_bh']:.4g}"
        print(
            f"  -> {row['target']} {row['group']} vs {row['control']}: fold {row['fold_change']} "
            f"({interval}), {p_text}, n={row['n_group']}/{row['n_control']}"
        )
    for item in result["excluded"]:
        print(f"[WARN] left out {item['blot']} / {item['sample']} / {item['target']}: {item['reason']}")
    for warning in result["warnings"]:
        print(f"[WARN] {warning}")
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[OK] wrote {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
