#!/usr/bin/env python3
"""Turn association or meta-analysis JSON into a forest plot and a draft report.

--result takes candidate_association.py JSON, and --meta takes candidate_meta.py
JSON. Writes forest.png, forest.pdf, and report.md to --out-dir. Open warnings go
first, because a draft should not read as finished while they stand. Figures need
matplotlib; the text report does not.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

from _genetics import InputError, safe_console

BLUE, ORANGE = "#0072B2", "#D55E00"


def format_p(value: float | None) -> str:
    """p-value in reporting style."""

    if value is None:
        return "not available"
    return "p < 0.001" if value < 0.001 else f"p = {value:.3f}"


def format_or(row: dict[str, Any]) -> str:
    """Odds ratio with its interval."""

    if row.get("or") is None:
        return "not estimable"
    return f"OR {row['or']:.2f} (95% CI {row['ci_low']:.2f}-{row['ci_high']:.2f})"


def association_markdown(result: dict[str, Any]) -> list[str]:
    """Report sections for a case-control run."""

    lines = ["# Genetic association results (draft)", "", "## Check before using this draft", ""]
    warnings = [warning for variant in result["variants"] for warning in variant["warnings"]]
    warnings += [f"{row['variant']} skipped: {row['reason']}" for row in result["skipped"]]
    lines += [f"- {warning}" for warning in warnings] or ["- No automated warnings."]
    models = sorted({variant["primary_model"] for variant in result["variants"]})
    lines += ["", "## Methods (draft)", ""]
    lines.append(
        "Hardy-Weinberg equilibrium was tested in controls with an exact test. Association was tested "
        f"under a pre-specified {', '.join(models)} model. Odds ratios are unadjusted for covariates. "
        f"P-values were Bonferroni-corrected for {len(result['variants'])} variant(s)."
    )
    lines += ["", "## Results (draft)", ""]
    for variant in result["variants"]:
        primary = variant["primary"]
        frequency = variant["effect_allele_frequency"]
        lines.append(
            f"- {variant['variant']} ({variant['effect_allele']} allele, {variant['primary_model']}): "
            f"frequency {frequency['cases']} in cases and {frequency['controls']} in controls; "
            f"{format_or(primary)}, {format_p(primary['p_value'])}, "
            f"corrected {format_p(primary['p_bonferroni'])}. "
            f"Controls Hardy-Weinberg exact p = {variant['hwe_controls']['exact_p']:.3g}, "
            f"F = {variant['hwe_controls']['inbreeding_f']}."
        )
    lines += ["", "## Genotype counts", "", "| Variant | Effect allele doses 0/1/2, cases | Controls |", "| --- | --- | --- |"]
    for variant in result["variants"]:
        lines.append(
            f"| {variant['variant']} | {variant['cases']['counts_by_dose']} | {variant['controls']['counts_by_dose']} |"
        )
    return lines


def meta_markdown(result: dict[str, Any]) -> list[str]:
    """Report sections for a meta-analysis."""

    lines = ["# Meta-analysis results (draft)", "", "## Check before using this draft", ""]
    lines += [f"- {warning}" for warning in result["warnings"]] or ["- No automated warnings."]
    lines += ["", "## Methods (draft)", ""]
    lines.append(
        "Log odds ratios were pooled by inverse variance, with DerSimonian-Laird random effects. "
        "Heterogeneity was assessed with Cochran's Q and I^2."
    )
    lines += ["", "## Results (draft)", ""]
    lines.append(
        f"- {result['k']} studies. Random effects: {format_or(result['random'])}, "
        f"{format_p(result['random']['p_value'])}. Fixed effect: {format_or(result['fixed'])}."
    )
    q_text = "not computed" if result["q_p"] is None else f"{result['q_p']:.3g}"
    lines.append(f"- Heterogeneity: Q = {result['q']:.2f} (p = {q_text}), I^2 = {result['i2_percent']:.0f}%.")
    if result["leave_one_out"]:
        lowest = min(row["or"] for row in result["leave_one_out"])
        highest = max(row["or"] for row in result["leave_one_out"])
        lines.append(f"- Leave-one-out pooled OR ranged from {lowest:.2f} to {highest:.2f}.")
    if result["hwe_sensitivity"]:
        lines.append(
            f"- Without {', '.join(result['hwe_sensitivity']['excluded'])}: "
            f"{format_or(result['hwe_sensitivity']['random'])}."
        )
    return lines


def forest_rows(association: dict[str, Any] | None, meta: dict[str, Any] | None) -> list[dict[str, Any]]:
    """Rows to plot, top to bottom."""

    rows = []
    if association and not meta:
        for variant in association["variants"]:
            primary = variant["primary"]
            if primary.get("or") is not None:
                rows.append(
                    {
                        "label": f"{variant['variant']} {variant['effect_allele']} ({variant['primary_model']})",
                        "or": primary["or"],
                        "low": primary["ci_low"],
                        "high": primary["ci_high"],
                        "summary": False,
                    }
                )
    if meta:
        for study in meta["studies"]:
            rows.append(
                {"label": f"{study['study']} ({study['weight_random_percent']}%)", "or": study["or"], "low": study["ci_low"], "high": study["ci_high"], "summary": False}
            )
        random = meta["random"]
        rows.append({"label": "Random effects", "or": random["or"], "low": random["ci_low"], "high": random["ci_high"], "summary": True})
    return rows


def draw_figure(rows: list[dict[str, Any]], out_dir: Path) -> str:
    """Forest plot of odds ratios on a log axis."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise InputError("figures need matplotlib; install it or use --no-figure") from exc

    figure, axis = plt.subplots(figsize=(6.2, max(2.0, 0.42 * len(rows) + 1.1)))
    positions = list(range(len(rows), 0, -1))
    for position, row in zip(positions, rows):
        color = ORANGE if row["summary"] else BLUE
        axis.plot([row["low"], row["high"]], [position, position], color=color, linewidth=1.6)
        axis.scatter([row["or"]], [position], marker="D" if row["summary"] else "s", s=46 if row["summary"] else 30, color=color, zorder=3)
    axis.axvline(1.0, color="#888888", linewidth=0.9, linestyle="--")
    axis.set_xscale("log")
    axis.set_yticks(positions, [row["label"] for row in rows], fontsize=8.5)
    lows = min(row["low"] for row in rows)
    highs = max(row["high"] for row in rows)
    axis.set_xlim(min(lows, 1.0) / 1.15, max(highs, 1.0) * 1.15)
    ticks = [value for value in (0.1, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 8, 16) if axis.get_xlim()[0] <= value <= axis.get_xlim()[1]]
    axis.set_xticks(ticks, [f"{value:g}" for value in ticks])
    axis.minorticks_off()
    axis.set_xlabel("Odds ratio (95% CI)", fontsize=9)
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(out_dir / "forest.png", dpi=300)
    figure.savefig(out_dir / "forest.pdf")
    plt.close(figure)
    return "forest.png"


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--result", help="JSON from candidate_association.py")
    parser.add_argument("--meta", help="JSON from candidate_meta.py")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--no-figure", action="store_true", help="write the text report only")
    args = parser.parse_args(argv)

    try:
        if not args.result and not args.meta:
            raise InputError("give --result, --meta, or both")
        association = meta = None
        if args.result:
            with open(args.result, encoding="utf-8") as handle:
                association = json.load(handle)
        if args.meta:
            with open(args.meta, encoding="utf-8") as handle:
                meta = json.load(handle)
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        lines: list[str] = []
        if association:
            lines += association_markdown(association)
        if meta:
            lines += ([""] if lines else []) + meta_markdown(meta)
        rows = forest_rows(association, meta)
        figure = None if args.no_figure or not rows else draw_figure(rows, out_dir)
        if figure:
            lines += ["", "## Figure", "", f"![Forest plot]({figure})"]
        (out_dir / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    except (InputError, OSError, json.JSONDecodeError, KeyError) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2
    print(f"[OK] wrote {out_dir / 'report.md'}")
    if figure:
        print(f"[OK] wrote {out_dir / figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
