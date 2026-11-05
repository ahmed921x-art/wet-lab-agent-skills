#!/usr/bin/env python3
"""Turn elisa_fit.py JSON into a standard-curve figure and a draft report.

Writes standard_curve.png, standard_curve.pdf, and report.md to --out-dir. The
figure marks the working range, so anything outside it reads as unquantified.
Needs matplotlib, and numpy to draw the curve.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from elisa_fit import InputError, logistic, safe_console

BLUE, ORANGE, GREEN = "#0072B2", "#D55E00", "#009E73"


def build_markdown(result: dict[str, Any], figure: str | None) -> str:
    """Draft report text from the analysis JSON."""

    curve = result["curve"]
    name = curve["model"].upper()
    lines = ["# ELISA results (draft)", "", "## Check before using this draft", ""]
    if result["warnings"]:
        lines += [f"- {warning}" for warning in result["warnings"]]
    else:
        lines.append("- No automated warnings.")
    lines += ["", "## Methods (draft)", ""]
    blank = (
        f"Mean blank OD ({result['blank_od']}) was subtracted."
        if result["blank_od"]
        else "No blank was subtracted."
    )
    lines.append(
        f"{blank} Concentrations were back-calculated from a {name} standard curve fitted to individual "
        f"wells with {curve['weight']} weighting. The working range was the longest run of standards "
        "meeting the recovery and precision criteria. Values outside it were reported as censored "
        "and were not extrapolated."
    )
    lines += ["", "## Curve", ""]
    lines.append(
        f"- Midpoint {curve['c']:.4g}, slope {curve['b']:.3f}, asymptotes {curve['a']:.3f} and {curve['d']:.3f} "
        f"({curve['direction']})."
    )
    lines.append(f"- Working range: {result['lloq']} to {result['uloq']} (well concentration).")
    lines += ["", "| Standard | Recovery % | OD CV % | In range |", "| --- | --- | --- | --- |"]
    for level in result["standards"]:
        lines.append(
            f"| {level['nominal']} | {level['recovery_percent']} | {level['od_cv_percent']} | "
            f"{'yes' if level.get('in_working_range') else 'no'} |"
        )
    lines += ["", "## Samples", "", "| ID | Dilution | Concentration | 95% CI | CV % | Status |", "| --- | --- | --- | --- | --- | --- |"]
    for row in result["unknowns"]:
        value = row["concentration"] if row["concentration"] is not None else row["censored"]
        interval = row.get("concentration_ci")
        ci_text = f"{interval[0]} to {interval[1]}" if interval else "-"
        lines.append(
            f"| {row['id']} | 1:{row['dilution']:g} | {value} | {ci_text} | {row['od_cv_percent']} | {row['status']} |"
        )
    censored = sum(1 for row in result["unknowns"] if row["censored"])
    if censored:
        lines += ["", f"{censored} result(s) are censored. State how they were handled in any group comparison."]
    if figure:
        lines += ["", "## Figure", "", f"![Standard curve]({figure})"]
    return "\n".join(lines) + "\n"


def draw_figure(result: dict[str, Any], out_dir: Path) -> str:
    """Plot standards, the fitted curve, samples, and the working range."""

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise InputError("figures need matplotlib; install it or use --no-figure") from exc

    curve = result["curve"]
    levels = [level for level in result["standards"] if not level["excluded"]]
    nominal = [level["nominal"] for level in levels]
    low, high = min(nominal) / 2.5, max(nominal) * 2.5
    x = np.geomspace(low, high, 300)
    y = logistic(x, curve["a"], curve["d"], curve["c"], curve["b"], curve["g"])

    figure, axis = plt.subplots(figsize=(5.4, 3.9))
    if result["lloq"] is not None:
        axis.axvspan(result["lloq"], result["uloq"], color=GREEN, alpha=0.10, linewidth=0, label="Working range")
    axis.plot(x, y, color=BLUE, linewidth=1.6, label=f"{curve['model'].upper()} fit")
    axis.scatter(nominal, [level["mean_net_od"] for level in levels], s=30, color=BLUE, zorder=3, label="Standards")
    excluded = [level for level in result["standards"] if level["excluded"]]
    if excluded:
        axis.scatter(
            [level["nominal"] for level in excluded],
            [level["mean_net_od"] for level in excluded],
            s=34,
            facecolors="none",
            edgecolors="#555555",
            zorder=3,
            label="Excluded standard",
        )
    placed = [row for row in result["unknowns"] if row["concentration"] is not None]
    if placed:
        axis.scatter(
            [row["well_concentration"] for row in placed],
            [row["mean_net_od"] for row in placed],
            s=30,
            marker="D",
            color=ORANGE,
            zorder=4,
            label="Samples in range",
        )
    axis.set_xscale("log")
    axis.set_xlabel("Concentration in well")
    axis.set_ylabel("Net OD")
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, fontsize=8)
    figure.tight_layout()
    figure.savefig(out_dir / "standard_curve.png", dpi=300)
    figure.savefig(out_dir / "standard_curve.pdf")
    plt.close(figure)
    return "standard_curve.png"


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--result", required=True, help="JSON from elisa_fit.py")
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
