#!/usr/bin/env python3
"""Turn qpcr_relative.py JSON into a figure and a draft Methods and Results report.

Writes fold_change.png, fold_change.pdf, and report.md to --out-dir. The figure
shows each biological sample with the group mean and confidence interval on a
log2 axis. Open warnings are listed first in the report, because a draft should
not read as finished while they stand. Needs matplotlib.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from _qpcr import InputError, safe_console

COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")


def format_p(value: float | None) -> str:
    """p-value in reporting style."""

    if value is None:
        return "not tested"
    return "p < 0.001" if value < 0.001 else f"p = {value:.3f}"


def build_markdown(result: dict[str, Any], efficiency: dict[str, float] | None, figure: str | None) -> str:
    """Draft report text from the analysis JSON."""

    model = (
        "the efficiency-corrected Pfaffl model"
        if result["method"] == "pfaffl"
        else "the comparative Ct (2^-ddCt) method, assuming 100% efficiency"
    )
    references = " and ".join(result["references"])
    level = round(result["ci_level"] * 100)
    designs = {row["design"] for row in result["contrasts"]}
    tests = " and ".join(
        {"welch": "Welch's t-test", "paired": "a paired t-test"}[design] for design in sorted(designs)
    )
    lines = ["# qPCR results (draft)", ""]
    lines += ["## Check before using this draft", ""]
    if result["warnings"]:
        lines += [f"- **{item['code']}**: {item['message']}" for item in result["warnings"]]
    else:
        lines.append("- No automated warnings. Melt curves still need checking on the instrument software.")
    lines += ["", "## Methods (draft)", ""]
    method_text = (
        f"Relative expression was calculated with {model}, normalized to {references}. "
        "Technical replicates were averaged per biological sample. "
        f"Statistical tests used {tests or 'no test'} on log2 relative quantities, with "
        f"{level}% confidence intervals back-transformed to fold change and "
        "Benjamini-Hochberg correction across all comparisons."
    )
    if efficiency:
        values = ", ".join(f"{gene} {(factor - 1) * 100:.1f}%" for gene, factor in sorted(efficiency.items()))
        method_text += f" Measured amplification efficiencies were {values}."
    lines += [method_text, "", "## Results (draft)", ""]
    for row in result["contrasts"]:
        interval = (
            f"{level}% CI {row['fold_ci_low']:.2f}-{row['fold_ci_high']:.2f}"
            if row["fold_ci_low"] is not None
            else "no interval available"
        )
        adjusted = "" if row["p_adjusted_bh"] is None else f", adjusted p = {row['p_adjusted_bh']:.3g}"
        lines.append(
            f"- {row['target']} in {row['group']} vs {row['control']}: {row['fold_change']:.2f}-fold "
            f"({interval}; {format_p(row['p_value'])}{adjusted}; "
            f"n = {row['n_group']} vs {row['n_control']})."
        )
    lines += ["", "## Reference stability", ""]
    stability = result["reference_stability"]
    lines.append(f"- Largest between-group shift: {stability['between_group_shift_cycles']} cycles.")
    if stability["genorm_m"]:
        lines.append(f"- geNorm M: {stability['genorm_m']}.")
    if figure:
        lines += ["", "## Figure", "", f"![Fold change by group]({figure})", ""]
        lines.append(
            f"Points are biological samples relative to the mean of {result['control']}. "
            f"Bars show each group's mean and {level}% confidence interval for the difference from control."
        )
    return "\n".join(lines) + "\n"


def draw_figure(result: dict[str, Any], out_dir: Path) -> str:
    """Draw one panel per target."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise InputError("figures need matplotlib; install it or use --no-figure") from exc

    samples = result.get("samples", [])
    if not samples:
        raise InputError("result JSON has no per-sample values; re-run qpcr_relative.py")
    targets = sorted({row["target"] for row in samples})
    control = result["control"]
    groups = [control] + sorted({row["group"] for row in samples} - {control})
    figure, axes = plt.subplots(
        1, len(targets), figsize=(max(3.2, 1.1 * len(groups) + 1.2) * len(targets), 3.6), squeeze=False, sharey=True
    )
    for axis, target in zip(axes[0], targets):
        members = [row for row in samples if row["target"] == target]
        baseline_values = [row["log2_relative"] for row in members if row["group"] == control]
        if not baseline_values:
            axis.set_title(f"{target}\n(no control data)", fontsize=10)
            continue
        baseline = sum(baseline_values) / len(baseline_values)
        for index, group in enumerate(groups):
            values = [row["log2_relative"] - baseline for row in members if row["group"] == group]
            color = COLORS[index % len(COLORS)]
            count = len(values)
            jitter = [(position - (count - 1) / 2) * min(0.32 / max(count - 1, 1), 0.08) for position in range(count)]
            axis.scatter([index + shift for shift in jitter], values, s=26, color=color, alpha=0.8, linewidth=0, zorder=3)
            contrast = next(
                (row for row in result["contrasts"] if row["target"] == target and row["group"] == group), None
            )
            if contrast is not None and contrast["fold_ci_low"] is not None:
                import math

                center = contrast["log2_fold_change"]
                low, high = math.log2(contrast["fold_ci_low"]), math.log2(contrast["fold_ci_high"])
                axis.errorbar(index + 0.26, center, yerr=[[center - low], [high - center]], color="#222222", capsize=3, linewidth=1.2, marker="o", markersize=4, zorder=4)
                label = "n/a" if contrast["p_adjusted_bh"] is None else f"q={contrast['p_adjusted_bh']:.2g}"
                axis.annotate(label, (index, 1.0), xycoords=("data", "axes fraction"), ha="center", va="top", fontsize=8)
        axis.axhline(0, color="#888888", linewidth=0.8, linestyle="--", zorder=1)
        axis.set_xticks(range(len(groups)), groups, rotation=20 if len(groups) > 3 else 0, fontsize=9)
        axis.set_xlim(-0.5, len(groups) - 0.3)
        axis.set_title(target, fontsize=10, pad=12)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0][0].set_ylabel(f"log2 fold change vs {control}", fontsize=9)
    figure.tight_layout()
    figure.savefig(out_dir / "fold_change.png", dpi=300)
    figure.savefig(out_dir / "fold_change.pdf")
    plt.close(figure)
    return "fold_change.png"


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--result", required=True, help="JSON from qpcr_relative.py")
    parser.add_argument("--efficiency-json", help="JSON from qpcr_efficiency.py")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--no-figure", action="store_true", help="write the text report only")
    args = parser.parse_args(argv)

    try:
        with open(args.result, encoding="utf-8") as handle:
            result = json.load(handle)
        efficiency = None
        if args.efficiency_json:
            with open(args.efficiency_json, encoding="utf-8") as handle:
                efficiency = json.load(handle).get("efficiency")
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        figure = None if args.no_figure else draw_figure(result, out_dir)
        (out_dir / "report.md").write_text(build_markdown(result, efficiency, figure), encoding="utf-8")
    except (InputError, OSError, json.JSONDecodeError, KeyError) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2
    print(f"[OK] wrote {out_dir / 'report.md'}")
    if figure:
        print(f"[OK] wrote {out_dir / figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
