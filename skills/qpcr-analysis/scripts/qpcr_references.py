#!/usr/bin/env python3
"""Rank candidate reference genes with the geNorm stepwise procedure.

Input is the same long table as qpcr_relative.py. Give two or more candidates.
The least stable gene is removed repeatedly to rank them, and pairwise variation
V(n/n+1) shows how many references are enough. Run this on samples that cover the
conditions under study. Ct values are converted with the given efficiencies.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
from collections import defaultdict
from typing import Any

from _qpcr import InputError, read_table, rounded, safe_console, write_json
from qpcr_relative import GENORM_M_LIMIT, load_rows, parse_efficiency_arguments

V_LIMIT = 0.15


def relative_quantities(
    raw_rows: list[dict[str, Any]],
    candidates: list[str],
    efficiency: dict[str, float],
) -> tuple[dict[str, dict[str, float]], list[str]]:
    """log2 quantity per sample and candidate; samples missing a candidate are dropped."""

    rows = load_rows(raw_rows)
    wells: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        if row["task"] == "unknown" and row["ct"] is not None and row["target"] in candidates:
            wells[(row["sample"], row["target"])].append(row["ct"])
    samples = sorted({sample for sample, _ in wells})
    complete: dict[str, dict[str, float]] = {}
    dropped = []
    for sample in samples:
        if all((sample, candidate) in wells for candidate in candidates):
            complete[sample] = {
                candidate: -math.log2(efficiency.get(candidate, 2.0))
                * statistics.fmean(wells[(sample, candidate)])
                for candidate in candidates
            }
        else:
            dropped.append(sample)
    return complete, dropped


def stability(values: dict[str, dict[str, float]], genes: list[str]) -> dict[str, float]:
    """geNorm M of each gene against the others still in the set."""

    result = {}
    for gene in genes:
        spreads = [
            statistics.stdev(sample[gene] - sample[other] for sample in values.values())
            for other in genes
            if other != gene
        ]
        result[gene] = statistics.fmean(spreads)
    return result


def normalization_factor(values: dict[str, float], genes: list[str]) -> float:
    """log2 of the geometric mean quantity across genes."""

    return statistics.fmean(values[gene] for gene in genes)


def rank(values: dict[str, dict[str, float]], candidates: list[str]) -> dict[str, Any]:
    """Stepwise exclusion, then pairwise variation of sequential factors."""

    remaining = list(candidates)
    steps = []
    while len(remaining) > 2:
        scores = stability(values, remaining)
        worst = max(scores, key=scores.get)
        steps.append(
            {
                "genes": sorted(remaining),
                "m": {gene: rounded(score) for gene, score in sorted(scores.items())},
                "removed": worst,
            }
        )
        remaining.remove(worst)
    last = stability(values, remaining)
    steps.append(
        {
            "genes": sorted(remaining),
            "m": {gene: rounded(score) for gene, score in sorted(last.items())},
            "removed": None,
        }
    )
    order = sorted(remaining) + [step["removed"] for step in reversed(steps) if step["removed"]]

    variation = []
    for count in range(2, len(order)):
        first, second = order[:count], order[: count + 1]
        differences = [
            normalization_factor(sample, first) - normalization_factor(sample, second)
            for sample in values.values()
        ]
        variation.append({"comparison": f"V{count}/{count + 1}", "v": rounded(statistics.stdev(differences))})
    recommended = len(order)
    for count, entry in enumerate(variation, start=2):
        if entry["v"] < V_LIMIT:
            recommended = count
            break
    if not variation:
        recommended = 2
    return {
        "ranking_most_stable_first": order,
        "steps": steps,
        "pairwise_variation": variation,
        "recommended_count": recommended,
        "recommended_references": order[:recommended],
    }


def treatment_shifted(
    values: dict[str, dict[str, float]],
    groups: dict[str, str],
    candidates: list[str],
    limit: float = 0.5,
) -> dict[str, float]:
    """Find candidates that move with group, removing the worst one at a time.

    geNorm scores agreement between genes, so a gene that shifts uniformly in one
    condition can still look stable. Each gene is compared with the mean of the
    others within each group. Iterating stops one strong mover from implicating
    the stable genes it is compared against.
    """

    remaining = list(candidates)
    moved: dict[str, float] = {}
    labels = sorted({groups[sample] for sample in values})
    while len(labels) > 1 and len(remaining) > 1:
        shifts = {}
        for gene in remaining:
            others = [other for other in remaining if other != gene]
            means = []
            for label in labels:
                members = [
                    sample[gene] - statistics.fmean(sample[other] for other in others)
                    for name, sample in values.items()
                    if groups[name] == label
                ]
                if members:
                    means.append(statistics.fmean(members))
            shifts[gene] = max(means) - min(means)
        worst = max(shifts, key=shifts.get)
        if shifts[worst] <= limit:
            break
        moved[worst] = shifts[worst]
        remaining.remove(worst)
    return moved


def analyze(
    raw_rows: list[dict[str, Any]],
    candidates: list[str],
    efficiency: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Rank candidates and report the warnings an analyst needs."""

    if len(set(candidates)) < 2:
        raise InputError("give at least two candidate references")
    values, dropped = relative_quantities(raw_rows, candidates, efficiency or {})
    if len(values) < 3:
        raise InputError(f"{len(values)} complete sample(s); stability needs at least 3")
    result = rank(values, candidates)
    warnings = []
    groups = {row["sample"]: row["group"] for row in load_rows(raw_rows) if row["task"] == "unknown"}
    moved = treatment_shifted(values, groups, candidates)
    result["shift_between_groups"] = {gene: rounded(shift) for gene, shift in moved.items()}
    if moved:
        warnings.append(
            "treatment appears to move: "
            + ", ".join(f"{gene} ({shift:.2f} cycles)" for gene, shift in moved.items())
            + "; geNorm alone cannot see this, so these are left out of the recommendation"
        )
        if len(candidates) > 2:
            result["recommended_references"] = [
                gene for gene in result["recommended_references"] if gene not in moved
            ]
            result["recommended_count"] = len(result["recommended_references"])
    if result["pairwise_variation"] and min(row["v"] for row in result["pairwise_variation"]) >= V_LIMIT:
        warnings.append(f"V never fell below {V_LIMIT}; that value is a guide, so judge the trend")
    if not result["recommended_references"]:
        warnings.append("no candidate can be recommended; test other genes")
    if len(values) < 8:
        warnings.append(f"only {len(values)} samples; rankings from few samples are unstable")
    if dropped:
        warnings.append(f"dropped samples missing a candidate: {', '.join(dropped)}")
    best_m = max(result["steps"][-1]["m"].values())
    if best_m >= GENORM_M_LIMIT:
        warnings.append(f"even the best pair has M {best_m}; test other candidates")
    if len(candidates) == 2:
        warnings.append("two candidates can be scored but not ranked: their M is identical")
    return {**result, "n_samples": len(values), "warnings": warnings}


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="CSV of wells")
    parser.add_argument("--candidate", required=True, action="append", help="repeat for each candidate")
    parser.add_argument("--efficiency", action="append", default=[], help="TARGET=FACTOR")
    parser.add_argument("--efficiency-json", help="output of qpcr_efficiency.py")
    parser.add_argument("--output-json", help="write full results here")
    args = parser.parse_args(argv)

    try:
        efficiency = parse_efficiency_arguments(args.efficiency, args.efficiency_json)
        result = analyze(read_table(args.input, ("sample", "group", "target", "ct")), args.candidate, efficiency)
    except (InputError, OSError, json.JSONDecodeError) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    print(f"[OK] {result['n_samples']} samples; most stable first: {', '.join(result['ranking_most_stable_first'])}")
    for step in result["steps"]:
        print(f"  -> M {step['m']}" + (f", removing {step['removed']}" if step["removed"] else ""))
    for entry in result["pairwise_variation"]:
        print(f"  -> {entry['comparison']} = {entry['v']}")
    print(f"[OK] use {result['recommended_count']} references: {', '.join(result['recommended_references'])}")
    for warning in result["warnings"]:
        print(f"[WARN] {warning}")
    if args.output_json:
        write_json(args.output_json, result)
        print(f"[OK] wrote {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
