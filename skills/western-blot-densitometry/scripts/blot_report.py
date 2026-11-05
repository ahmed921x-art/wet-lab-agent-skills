#!/usr/bin/env python3
"""Turn blot_quant.py JSON into a figure and a draft Methods and Results report.

Writes blot_fold_change.png, blot_fold_change.pdf, and report.md to --out-dir.
Needs matplotlib for the figure.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

from blot_quant import InputError, safe_console

COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")


def format_p(value: float | None) -> str:
    """p-value in reporting style."""

    if value is None:
        return "not tested"
    return "p < 0.001" if value < 0.001 else f"p = {value:.3f}"


def build_markdown(result: dict[str, Any], figure: str | None) -> str:
    """Draft report text from the analysis JSON."""

    level = round(result["ci_level"] * 100)
    lines = ["# Western blot results (draft)", "", "## Check before using this draft", ""]
    lines += [f"- {warning}" for warning in result["warnings"]] or ["- No automated warnings."]
    lines.append("- Confirm that all bands were within the linear range and that uncropped blots are available.")
    lines += ["", "## Methods (draft)", ""]
    lines.append(
        f"Band densities were background-subtracted and normalized to {result['loading']} within each "
        f"membrane, then expressed relative to the mean of {result['control']} lanes on that membrane. "
        "Repeat blots of a sample were averaged. Groups were compared on log2 values, with "
        f"{level}% confidence intervals back-transformed to fold change and Benjamini-Hochberg correction."
    )
    lines += ["", "## Results (draft)", ""]
    for row in result["contrasts"]:
        interval = (
            f"{level}% CI {row['fold_ci_low']:.2f}-{row['fold_ci_high']:.2f}"
            if row["fold_ci_low"] is not None
            else "no interval available"
        )
        lines.append(
            f"- {row['target']} in {row['group']} vs {row['control']}: {row['fold_change']:.2f}-fold "
            f"({interval}; {format_p(row['p_value'])}; n = {row['n_group']} vs {row['n_control']})."
        )
    if figure:
        lines += ["", "## Figure", "", f"![Normalized band density]({figure})"]
    return "\n".join(lines) + "\n"


def draw_figure(result: dict[str, Any], out_dir: Path) -> str:
    """Plot each sample and each contrast's interval, one panel per target."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise InputError("figures need matplotlib; install it or use --no-figure") from exc

    samples = result["samples"]
    if not samples:
        raise InputError("no sample values to plot")
    control = result["control"]
    targets = sorted({row["target"] for row in samples})
    groups = [control] + sorted({row["group"] for row in samples} - {control})
    figure, axes = plt.subplots(
        1, len(targets), figsize=(max(3.0, 1.1 * len(groups) + 1.2) * len(targets), 3.5), squeeze=False, sharey=True
    )
    for axis, target in zip(axes[0], targets):
        members = [row for row in samples if row["target"] == target]
        control_values = [row["log2_relative"] for row in members if row["group"] == control]
        baseline = sum(control_values) / len(control_values) if control_values else 0.0
        for index, group in enumerate(groups):
            values = [row["log2_relative"] - baseline for row in members if row["group"] == group]
            count = len(values)
            jitter = [(position - (count - 1) / 2) * min(0.32 / max(count - 1, 1), 0.08) for position in range(count)]
            axis.scatter([index + shift for shift in jitter], values, s=26, color=COLORS[index % len(COLORS)], alpha=0.85, linewidth=0, zorder=3)
            contrast = next((row for row in result["contrasts"] if row["target"] == target and row["group"] == group), None)
            if contrast is not None and contrast["fold_ci_low"] is not None:
                center = contrast["log2_fold_change"]
                low, high = math.log2(contrast["fold_ci_low"]), math.log2(contrast["fold_ci_high"])
                axis.errorbar(index + 0.26, center, yerr=[[center - low], [high - center]], color="#222222", capsize=3, linewidth=1.2, marker="o", markersize=4, zorder=4)
        axis.axhline(0, color="#888888", linewidth=0.8, linestyle="--", zorder=1)
        axis.set_xticks(range(len(groups)), groups, fontsize=9)
        axis.set_xlim(-0.5, len(groups) - 0.3)
        axis.set_title(target, fontsize=10)
        axis.spines[["top", "right"]].set_visible(False)
    axes[0][0].set_ylabel(f"log2 fold change vs {control}\n(normalized to {result['loading']})", fontsize=9)
    figure.tight_layout()
    figure.savefig(out_dir / "blot_fold_change.png", dpi=300)
    figure.savefig(out_dir / "blot_fold_change.pdf")
    plt.close(figure)
    return "blot_fold_change.png"


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--result", required=True, help="JSON from blot_quant.py")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--no-figure", action="store_true", help="write the text report only")
    args = parser.parse_args(argv)
    try:
        with open(args.result, encoding="utf-8") as handle:
            result = json.load(handle)
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        figure = None if args.no_figure else draw_figure(result, out_dir)
        (out_dir / "report.md").write_text(build_markdown(result, figure), encoding="utf-8")
    except (InputError, OSError, json.JSONDecodeError, KeyError) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2
    print(f"[OK] wrote {out_dir / 'report.md'}")
    if figure:
        print(f"[OK] wrote {out_dir / figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
