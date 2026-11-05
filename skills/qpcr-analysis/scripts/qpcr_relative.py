#!/usr/bin/env python3
"""Relative qPCR quantification with plate QC, from a long table of wells.

Input columns: sample, group, target, ct. Optional: task (unknown, ntc, nrt) and
pair (a shared id for matched samples, which switches to a paired test).

Each sample gets log2 relative quantity
    L = mean over references(log2(E_ref) * Ct_ref) - log2(E_target) * Ct_target
which equals -dCt when every E is 2 (Livak) and reproduces the Pfaffl ratio
otherwise. Statistics are run on L per biological sample, never on fold changes
or on technical wells.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from typing import Any

from scipy import stats

from _qpcr import (
    InputError,
    benjamini_hochberg,
    parse_ct,
    read_table,
    rounded,
    safe_console,
    write_json,
)

TASK_ALIASES = {
    "": "unknown",
    "unknown": "unknown",
    "unk": "unknown",
    "sample": "unknown",
    "ntc": "ntc",
    "nrt": "nrt",
    "no-rt": "nrt",
    "-rt": "nrt",
}
GENORM_M_LIMIT = 1.5


def load_rows(raw_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Normalize raw table rows and reject inconsistent labels."""

    rows = []
    sample_group: dict[str, str] = {}
    sample_pair: dict[str, str] = {}
    for raw in raw_rows:
        line = raw["_line"]
        task_text = raw.get("task", "").lower()
        if task_text not in TASK_ALIASES:
            raise InputError(f"line {line}: unknown task {raw.get('task')!r}")
        task = TASK_ALIASES[task_text]
        sample, group, target = raw["sample"], raw["group"], raw["target"]
        if not target or not sample:
            raise InputError(f"line {line}: sample and target are required")
        pair = raw.get("pair", "")
        if task == "unknown":
            if not group:
                raise InputError(f"line {line}: group is required for sample {sample}")
            if sample_group.setdefault(sample, group) != group:
                raise InputError(
                    f"line {line}: sample {sample} appears in groups "
                    f"{sample_group[sample]} and {group}"
                )
            if sample_pair.setdefault(sample, pair) != pair:
                raise InputError(f"line {line}: sample {sample} has more than one pair id")
        rows.append(
            {
                "sample": sample,
                "group": group,
                "target": target,
                "task": task,
                "pair": pair,
                "ct": parse_ct(raw["ct"], line=line),
            }
        )
    return rows


def check_controls(rows, unknown_means, margin, warnings) -> list[dict[str, Any]]:
    """Report signal in no-template and no-RT wells."""

    findings = []
    seen: dict[tuple[str, str], list[float | None]] = defaultdict(list)
    for row in rows:
        if row["task"] != "unknown":
            seen[(row["task"], row["target"])].append(row["ct"])
    for (task, target), cts in sorted(seen.items()):
        detected = [ct for ct in cts if ct is not None]
        entry = {
            "task": task,
            "target": target,
            "n_wells": len(cts),
            "n_detected": len(detected),
            "lowest_ct": rounded(min(detected), 2) if detected else None,
        }
        findings.append(entry)
        if not detected:
            continue
        reference_ct = unknown_means.get(target)
        gap = None if reference_ct is None else min(detected) - reference_ct
        entry["cycles_above_sample_mean"] = rounded(gap, 2)
        if gap is not None and gap < margin:
            warnings.append(
                {
                    "code": "control_contamination",
                    "message": (
                        f"{task.upper()} for {target} amplified at Ct {min(detected):.2f}, only "
                        f"{gap:.2f} cycles after the sample mean: treat this target as "
                        "unreliable until the source is found"
                    ),
                }
            )
        else:
            warnings.append(
                {
                    "code": "control_signal",
                    "message": (
                        f"{task.upper()} for {target} amplified at Ct {min(detected):.2f}; "
                        "inspect the melt curve"
                    ),
                }
            )
    return findings


def summarize_wells(rows, min_wells, max_well_sd, ct_ceiling, warnings) -> dict[tuple[str, str], Any]:
    """Collapse technical wells to one value per sample and target."""

    wells: dict[tuple[str, str], list[float | None]] = defaultdict(list)
    for row in rows:
        if row["task"] == "unknown":
            wells[(row["sample"], row["target"])].append(row["ct"])

    summary = {}
    for (sample, target), cts in sorted(wells.items()):
        detected = [ct for ct in cts if ct is not None]
        entry: dict[str, Any] = {
            "sample": sample,
            "target": target,
            "n_wells": len(cts),
            "n_detected": len(detected),
            "mean_ct": None,
            "sd_ct": None,
            "status": "ok",
            "flags": [],
        }
        if len(detected) < min_wells:
            entry["status"] = "not_detected" if not detected else "too_few_wells"
        else:
            entry["mean_ct"] = statistics.fmean(detected)
            if len(detected) > 1:
                entry["sd_ct"] = statistics.stdev(detected)
                if entry["sd_ct"] > max_well_sd:
                    entry["flags"].append(f"well SD {entry['sd_ct']:.2f} > {max_well_sd}")
            if entry["mean_ct"] > ct_ceiling:
                entry["flags"].append(f"mean Ct {entry['mean_ct']:.2f} > {ct_ceiling}")
            if len(detected) < len(cts):
                entry["flags"].append(
                    f"{len(cts) - len(detected)} of {len(cts)} wells undetermined"
                )
        summary[(sample, target)] = entry
        for flag in entry["flags"]:
            warnings.append(
                {"code": "well_qc", "message": f"{sample} / {target}: {flag}"}
            )
    return summary


def genorm_m(terms: dict[str, dict[str, float]], references: list[str]) -> dict[str, float]:
    """geNorm stability M: mean pairwise SD of log2 ratios to the other references."""

    values = {}
    for reference in references:
        sds = []
        for other in references:
            if other == reference:
                continue
            ratios = [
                sample_terms[reference] - sample_terms[other]
                for sample_terms in terms.values()
            ]
            if len(ratios) > 1:
                sds.append(statistics.stdev(ratios))
        if sds:
            values[reference] = statistics.fmean(sds)
    return values


def compare(values_group, values_control, pairs, alpha) -> dict[str, Any]:
    """Compare log2 relative quantities with a Welch or paired t-test."""

    out: dict[str, Any] = {
        "design": "paired" if pairs else "welch",
        "t": None,
        "df": None,
        "p_value": None,
        "diff_ci": None,
    }
    if pairs:
        differences = [group - control for group, control in pairs]
        diff = statistics.fmean(differences)
        count = len(differences)
        if count < 2:
            return {**out, "diff": diff}
        standard_error = statistics.stdev(differences) / math.sqrt(count)
        df = float(count - 1)
    else:
        diff = statistics.fmean(values_group) - statistics.fmean(values_control)
        n_group, n_control = len(values_group), len(values_control)
        if n_group < 2 or n_control < 2:
            return {**out, "diff": diff}
        var_group = statistics.variance(values_group) / n_group
        var_control = statistics.variance(values_control) / n_control
        standard_error = math.sqrt(var_group + var_control)
        if standard_error == 0:
            return {**out, "diff": diff}
        df = (var_group + var_control) ** 2 / (
            var_group**2 / (n_group - 1) + var_control**2 / (n_control - 1)
        )
    if standard_error == 0:
        return {**out, "diff": diff}
    t_value = diff / standard_error
    critical = float(stats.t.ppf(1 - alpha / 2, df))
    return {
        **out,
        "diff": diff,
        "t": t_value,
        "df": df,
        "p_value": float(2 * stats.t.sf(abs(t_value), df)),
        "diff_ci": (diff - critical * standard_error, diff + critical * standard_error),
    }


def analyze(
    raw_rows: list[dict[str, Any]],
    references: list[str],
    control: str,
    efficiency: dict[str, float] | None = None,
    *,
    min_wells: int = 2,
    max_well_sd: float = 0.5,
    ct_ceiling: float = 35.0,
    ntc_margin: float = 5.0,
    ref_shift_limit: float = 0.5,
    alpha: float = 0.05,
) -> dict[str, Any]:
    """Run plate QC and relative quantification."""

    if min_wells < 1:
        raise InputError("--min-wells must be at least 1")
    rows = load_rows(raw_rows)
    warnings: list[dict[str, str]] = []
    if min_wells == 1:
        warnings.append(
            {
                "code": "single_well_allowed",
                "message": "--min-wells 1 accepts unreplicated wells, so well precision is unchecked",
            }
        )
    unknown = [row for row in rows if row["task"] == "unknown"]
    groups = {row["sample"]: row["group"] for row in unknown}
    pair_ids = {row["sample"]: row["pair"] for row in unknown}
    targets_present = {row["target"] for row in unknown}
    if control not in groups.values():
        raise InputError(
            f"control group {control!r} not found; groups are {sorted(set(groups.values()))}"
        )
    absent = [reference for reference in references if reference not in targets_present]
    if absent:
        raise InputError(f"reference target(s) not on the plate: {', '.join(absent)}")
    targets = sorted(targets_present - set(references))
    if not targets:
        raise InputError("no target left after removing the references")

    efficiency = dict(efficiency or {})
    for target, factor in efficiency.items():
        if not 1.0 < factor <= 3.0:
            raise InputError(
                f"efficiency for {target} is {factor}; give the amplification factor "
                "(1.95 for 95%), not a percentage"
            )
    missing_efficiency = sorted((targets_present) - set(efficiency))
    if efficiency and missing_efficiency:
        warnings.append(
            {
                "code": "efficiency_assumed",
                "message": f"E assumed to be 2.0 for: {', '.join(missing_efficiency)}",
            }
        )
    weight = {name: math.log2(efficiency.get(name, 2.0)) for name in targets_present}
    method = "livak" if all(value == 1.0 for value in weight.values()) else "pfaffl"
    if not efficiency:
        warnings.append(
            {
                "code": "efficiency_unverified",
                "message": (
                    "All efficiencies taken as 100%. This needs a validation curve for "
                    "each primer pair; pass measured values with --efficiency-json."
                ),
            }
        )

    wells = summarize_wells(rows, min_wells, max_well_sd, ct_ceiling, warnings)
    mean_cts: dict[str, list[float]] = defaultdict(list)
    for (_, target), entry in wells.items():
        if entry["mean_ct"] is not None:
            mean_cts[target].append(entry["mean_ct"])
    controls = check_controls(
        rows,
        {target: statistics.fmean(values) for target, values in mean_cts.items()},
        ntc_margin,
        warnings,
    )

    ref_terms: dict[str, dict[str, float]] = {}
    excluded: dict[str, str] = {}
    for sample in sorted(groups):
        terms = {}
        for reference in references:
            entry = wells.get((sample, reference))
            if entry is None or entry["mean_ct"] is None:
                excluded[sample] = f"reference {reference} unusable"
                break
            terms[reference] = weight[reference] * entry["mean_ct"]
        else:
            ref_terms[sample] = terms
    for sample, reason in excluded.items():
        warnings.append(
            {"code": "sample_excluded", "message": f"sample {sample} excluded: {reason}"}
        )

    reference_term = {
        sample: statistics.fmean(terms.values()) for sample, terms in ref_terms.items()
    }
    by_group: dict[str, list[float]] = defaultdict(list)
    for sample, value in reference_term.items():
        by_group[groups[sample]].append(value)
    group_means = {group: statistics.fmean(values) for group, values in by_group.items()}
    shift = max(group_means.values()) - min(group_means.values()) if group_means else None
    if shift is not None and shift > ref_shift_limit:
        warnings.append(
            {
                "code": "reference_shift",
                "message": (
                    f"Reference signal differs by {shift:.2f} cycles between groups, so "
                    "the treatment may regulate it and any fold change inherits that shift"
                ),
            }
        )
    m_values = genorm_m(ref_terms, references)
    for reference, value in m_values.items():
        if value >= GENORM_M_LIMIT:
            warnings.append(
                {
                    "code": "reference_unstable",
                    "message": f"geNorm M for {reference} is {value:.2f} (limit {GENORM_M_LIMIT})",
                }
            )

    contrasts = []
    detection = []
    sample_values = []
    other_groups = sorted(set(groups.values()) - {control})
    if not other_groups:
        raise InputError("only the control group is present, so there is nothing to compare")
    for target in targets:
        log2_values: dict[str, float] = {}
        for sample, ref_value in reference_term.items():
            entry = wells.get((sample, target))
            if entry is None:
                continue
            if entry["mean_ct"] is not None:
                log2_values[sample] = ref_value - weight[target] * entry["mean_ct"]
                sample_values.append(
                    {
                        "sample": sample,
                        "group": groups[sample],
                        "target": target,
                        "log2_relative": rounded(log2_values[sample]),
                    }
                )

        for group in sorted(set(groups.values())):
            members = [sample for sample in reference_term if groups[sample] == group]
            usable = [sample for sample in members if sample in log2_values]
            detection.append(
                {"target": target, "group": group, "n_samples": len(members), "n_usable": len(usable)}
            )
            if len(usable) < len(members):
                warnings.append(
                    {
                        "code": "target_not_quantified",
                        "message": (
                            f"{target} usable in {len(usable)} of {len(members)} {group} "
                            "samples; the rest are left out, not imputed, so report the "
                            "detection rate next to any fold change"
                        ),
                    }
                )

        control_samples = [s for s in log2_values if groups[s] == control]
        for group in other_groups:
            group_samples = [s for s in log2_values if groups[s] == group]
            if not group_samples or not control_samples:
                warnings.append(
                    {
                        "code": "contrast_skipped",
                        "message": f"{target}: {group} vs {control} has no usable samples",
                    }
                )
                continue
            pairs = None
            contrast_samples = group_samples + control_samples
            labelled = [s for s in contrast_samples if pair_ids.get(s)]
            if labelled and len(labelled) < len(contrast_samples):
                warnings.append(
                    {
                        "code": "partial_pair_labels",
                        "message": (
                            f"{target}: only some {group} vs {control} samples have a pair "
                            "id, so the unpaired test was used"
                        ),
                    }
                )
            if len(labelled) == len(contrast_samples):
                control_by_pair = {pair_ids[s]: s for s in control_samples}
                matched = [
                    (log2_values[s], log2_values[control_by_pair[pair_ids[s]]])
                    for s in group_samples
                    if pair_ids[s] in control_by_pair
                ]
                if len(matched) < max(len(group_samples), len(control_samples)):
                    warnings.append(
                        {
                            "code": "unmatched_pairs",
                            "message": (
                                f"{target}: {group} vs {control} kept {len(matched)} "
                                "complete pairs and dropped the unmatched samples"
                            ),
                        }
                    )
                if not matched:
                    continue
                pairs = matched
            result = compare(
                [log2_values[s] for s in group_samples],
                [log2_values[s] for s in control_samples],
                pairs,
                alpha,
            )
            n_group, n_control = len(group_samples), len(control_samples)
            if pairs:
                n_group = n_control = len(pairs)
            if min(n_group, n_control) < 3:
                warnings.append(
                    {
                        "code": "low_n",
                        "message": (
                            f"{target}: {group} vs {control} has n={n_group} vs n={n_control} "
                            "biological samples; the test has very little power"
                        ),
                    }
                )
            diff = result["diff"]
            ci = result["diff_ci"]
            contrasts.append(
                {
                    "target": target,
                    "group": group,
                    "control": control,
                    "n_group": n_group,
                    "n_control": n_control,
                    "design": result["design"],
                    "delta_delta_ct": rounded(-diff),
                    "log2_fold_change": rounded(diff),
                    "fold_change": rounded(2**diff),
                    "fold_ci_low": rounded(2 ** ci[0]) if ci else None,
                    "fold_ci_high": rounded(2 ** ci[1]) if ci else None,
                    "t": rounded(result["t"]),
                    "df": rounded(result["df"], 2),
                    "p_value": result["p_value"],
                }
            )

    adjusted = benjamini_hochberg([row["p_value"] for row in contrasts])
    for row, value in zip(contrasts, adjusted):
        row["p_adjusted_bh"] = value

    return {
        "method": method,
        "references": references,
        "control": control,
        "ci_level": 1 - alpha,
        "n_samples": len(groups),
        "n_samples_used": len(reference_term),
        "contrasts": contrasts,
        "samples": sample_values,
        "detection": detection,
        "reference_stability": {
            "group_mean_reference_term": {k: rounded(v) for k, v in sorted(group_means.items())},
            "between_group_shift_cycles": rounded(shift),
            "genorm_m": {k: rounded(v) for k, v in m_values.items()},
        },
        "controls": controls,
        "wells": [
            {
                **entry,
                "mean_ct": rounded(entry["mean_ct"], 3),
                "sd_ct": rounded(entry["sd_ct"], 3),
            }
            for entry in wells.values()
        ],
        "warnings": warnings,
    }


def parse_efficiency_arguments(pairs: list[str], json_path: str | None) -> dict[str, float]:
    """Merge --efficiency-json with --efficiency TARGET=FACTOR overrides."""

    efficiency: dict[str, float] = {}
    if json_path:
        with open(json_path, encoding="utf-8") as handle:
            payload = json.load(handle)
        efficiency.update({k: float(v) for k, v in payload.get("efficiency", {}).items()})
    for item in pairs:
        target, separator, value = item.partition("=")
        if not separator:
            raise InputError(f"--efficiency needs TARGET=FACTOR, got {item!r}")
        try:
            efficiency[target.strip()] = float(value)
        except ValueError as exc:
            raise InputError(f"--efficiency factor is not a number: {item!r}") from exc
    return efficiency


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="CSV of wells")
    parser.add_argument("--reference", required=True, action="append", help="repeat for each reference")
    parser.add_argument("--control", required=True, help="control group name")
    parser.add_argument("--efficiency", action="append", default=[], help="TARGET=FACTOR, e.g. ACTB=1.94")
    parser.add_argument("--efficiency-json", help="output of qpcr_efficiency.py")
    parser.add_argument("--min-wells", type=int, default=2, help="detected wells needed per sample")
    parser.add_argument("--max-well-sd", type=float, default=0.5, help="flag above this technical SD")
    parser.add_argument("--ct-ceiling", type=float, default=35.0, help="flag mean Ct above this")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--output-json", help="write full results here")
    parser.add_argument("--output-csv", help="write the contrast table here")
    args = parser.parse_args(argv)

    try:
        efficiency = parse_efficiency_arguments(args.efficiency, args.efficiency_json)
        result = analyze(
            read_table(args.input, ("sample", "group", "target", "ct")),
            args.reference,
            args.control,
            efficiency,
            min_wells=args.min_wells,
            max_well_sd=args.max_well_sd,
            ct_ceiling=args.ct_ceiling,
            alpha=args.alpha,
        )
    except (InputError, OSError, json.JSONDecodeError) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    print(
        f"[OK] {result['method']} model, {result['n_samples_used']} of "
        f"{result['n_samples']} samples used"
    )
    for row in result["contrasts"]:
        interval = (
            f"CI {row['fold_ci_low']}-{row['fold_ci_high']}"
            if row["fold_ci_low"] is not None
            else "no CI"
        )
        p_text = "p n/a" if row["p_value"] is None else f"p {row['p_value']:.4g}, BH {row['p_adjusted_bh']:.4g}"
        print(
            f"  -> {row['target']} {row['group']} vs {row['control']}: "
            f"fold {row['fold_change']} ({interval}), {p_text}, "
            f"n={row['n_group']}/{row['n_control']}"
        )
    for warning in result["warnings"]:
        print(f"[WARN] {warning['code']}: {warning['message']}")

    if args.output_json:
        write_json(args.output_json, result)
        print(f"[OK] wrote {args.output_json}")
    if args.output_csv:
        columns = list(result["contrasts"][0]) if result["contrasts"] else []
        with open(args.output_csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(result["contrasts"])
        print(f"[OK] wrote {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
