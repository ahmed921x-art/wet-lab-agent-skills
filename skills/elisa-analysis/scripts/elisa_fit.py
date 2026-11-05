#!/usr/bin/env python3
"""Fit a 4PL or 5PL ELISA standard curve and back-calculate samples from one plate.

Input columns: kind, id, conc, dilution, od. kind is standard, blank, sample, or
qc. conc is nominal for standards and qc. dilution applies to samples and qc and
defaults to 1. Every row is one well.

The fit uses wells rather than means, reports a working range from standard
recovery, and refuses to extrapolate: values outside the range are returned as
censored bounds.
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

import numpy as np
from scipy.optimize import least_squares
from scipy.stats import t as student_t

MIN_LEVELS = {"4pl": 5, "5pl": 7}
HARD_MIN_LEVELS = {"4pl": 4, "5pl": 6}
WEIGHTS = ("1/y2", "1/y", "none")


class InputError(ValueError):
    """Raised for input problems the analyst has to resolve."""


def read_text(path: str | Path) -> str:
    """Read UTF-8 or UTF-16 text, and explain any other encoding."""

    source = Path(path)
    if not source.is_file():
        raise InputError(f"input file not found: {source}")
    raw = source.read_bytes()
    try:
        if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
            text = raw.decode("utf-16")
        else:
            text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputError(f"{source} is not UTF-8; re-save it as CSV UTF-8") from exc
    if not text.strip():
        raise InputError(f"input file is empty: {source}")
    return text


def safe_console() -> None:
    """Keep non-Latin sample names from crashing output on legacy Windows code pages."""

    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")


def sig(value: float | None, digits: int = 6) -> float | None:
    """Round to significant digits, so pg/mL and ng/mL scales both survive."""

    return None if value is None else float(f"{value:.{digits}g}")


def to_float(raw: str, *, line: int, column: str) -> float:
    """Parse a required number with a useful error."""

    text = str(raw).strip()
    try:
        value = float(text.replace(",", ".") if text.count(",") == 1 and "." not in text else text)
    except ValueError as exc:
        raise InputError(f"line {line}: {column} {raw!r} is not a number") from exc
    if not math.isfinite(value):
        raise InputError(f"line {line}: {column} {raw!r} is not finite")
    return value


def read_wells(path: str | Path) -> list[dict[str, Any]]:
    """Read one plate of wells from a delimited file."""

    text = read_text(path)
    delimiter = max(",;\t", key=text.splitlines()[0].count)
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    reader.fieldnames = [name.strip().lower() for name in reader.fieldnames or []]
    missing = [name for name in ("kind", "id", "od") if name not in reader.fieldnames]
    if missing:
        raise InputError(f"missing column(s): {', '.join(missing)}")

    wells = []
    plates = set()
    for line, raw in enumerate(reader, start=2):
        if None in raw:
            raise InputError(f"line {line}: more fields than the header")
        row = {key: (value or "").strip() for key, value in raw.items()}
        if not any(row.values()):
            continue
        kind = row["kind"].lower()
        if kind not in {"standard", "blank", "sample", "qc"}:
            raise InputError(f"line {line}: unknown kind {row['kind']!r}")
        if row.get("plate"):
            plates.add(row["plate"])
        conc = None
        if kind in {"standard", "qc"}:
            if not row.get("conc"):
                raise InputError(f"line {line}: {kind} needs a nominal conc")
            conc = to_float(row["conc"], line=line, column="conc")
            if conc < 0:
                raise InputError(f"line {line}: conc is negative")
        dilution = 1.0
        if kind in {"sample", "qc"} and row.get("dilution"):
            dilution = to_float(row["dilution"], line=line, column="dilution")
            if dilution < 1:
                raise InputError(
                    f"line {line}: dilution {dilution} is below 1; give the dilution "
                    "factor (10 for 1:10), not the fraction"
                )
        if not row["od"]:
            raise InputError(f"line {line}: od is empty")
        wells.append(
            {
                "kind": kind,
                "id": row["id"],
                "conc": conc,
                "dilution": dilution,
                "od": to_float(row["od"], line=line, column="od"),
            }
        )
    if len(plates) > 1:
        raise InputError(
            f"file holds plates {sorted(plates)}; a standard curve belongs to one "
            "plate, so run each plate separately"
        )
    return wells


def logistic(x: np.ndarray, a: float, d: float, c: float, b: float, g: float = 1.0) -> np.ndarray:
    """Response at concentration x. a is the x->0 asymptote, d the x->inf asymptote."""

    return d + (a - d) / (1.0 + (x / c) ** b) ** g


def invert(y: float, a: float, d: float, c: float, b: float, g: float = 1.0) -> float | None:
    """Concentration for response y, or None when y is outside the asymptotes."""

    low, high = min(a, d), max(a, d)
    if not low < y < high:
        return None
    inner = ((a - d) / (y - d)) ** (1.0 / g) - 1.0
    return float(c * inner ** (1.0 / b)) if inner > 0 else None


def sigma_for(y: np.ndarray, weight: str) -> np.ndarray:
    """Per-well sigma for a weighting scheme, floored so low wells stay finite."""

    if weight == "none":
        return np.ones_like(y)
    floor = max(float(np.max(np.abs(y))) * 0.01, 1e-6)
    scale = np.maximum(np.abs(y), floor)
    return scale if weight == "1/y2" else np.sqrt(scale)


def fit_curve(x: np.ndarray, y: np.ndarray, model: str, weight: str) -> dict[str, Any]:
    """Weighted least squares with several starts; the best converged fit wins."""

    sigma = sigma_for(y, weight)
    five = model == "5pl"
    span = float(np.max(y) - np.min(y))
    if span <= 0:
        raise InputError("standards show no change in signal")
    order = np.argsort(x)
    rising = np.mean(y[order][-3:]) > np.mean(y[order][:3])
    a0, d0 = (float(np.min(y)), float(np.max(y))) if rising else (float(np.max(y)), float(np.min(y)))
    centers = np.geomspace(float(np.min(x)), float(np.max(x)), 5)

    def residuals(theta: np.ndarray) -> np.ndarray:
        a, d, log_c, log_b = theta[:4]
        g = math.exp(theta[4]) if five else 1.0
        return (logistic(x, a, d, math.exp(log_c), math.exp(log_b), g) - y) / sigma

    log_c_bounds = (math.log(float(np.min(x))) - 7.0, math.log(float(np.max(x))) + 7.0)
    lower = [-np.inf, -np.inf, log_c_bounds[0], math.log(0.05)] + ([-4.0] if five else [])
    upper = [np.inf, np.inf, log_c_bounds[1], math.log(20.0)] + ([4.0] if five else [])

    best = None
    for center in centers:
        for slope in (0.7, 1.0, 1.5, 2.5):
            start = [a0, d0, math.log(center), math.log(slope)] + ([0.0] if five else [])
            try:
                solved = least_squares(
                    residuals, start, bounds=(lower, upper), x_scale="jac", max_nfev=4000
                )
            except ValueError:
                continue
            if solved.success and (best is None or solved.cost < best.cost):
                best = solved
    if best is None:
        raise InputError("curve fit did not converge")

    a, d = float(best.x[0]), float(best.x[1])
    c = math.exp(best.x[2])
    b = math.exp(best.x[3])
    g = math.exp(best.x[4]) if five else 1.0
    fitted = logistic(x, a, d, c, b, g)
    total = float(np.sum((y - np.mean(y)) ** 2))
    df = len(x) - len(best.x)
    covariance = None
    scale = None
    if df > 0:
        scale = float(2.0 * best.cost / df)
        try:
            covariance = (np.linalg.inv(best.jac.T @ best.jac) * scale).tolist()
        except np.linalg.LinAlgError:
            covariance = None
    return {
        "model": model,
        "weight": weight,
        "direction": "sandwich" if d > a else "competitive",
        "a": a,
        "d": d,
        "c": c,
        "b": b,
        "g": g,
        "r_squared_unweighted": 1.0 - float(np.sum((y - fitted) ** 2)) / total,
        "n_wells": int(len(x)),
        "n_parameters": 5 if five else 4,
        "df": df,
        "residual_variance": scale,
        "signal_floor": max(float(np.max(np.abs(y))) * 0.01, 1e-6),
        "theta": [float(value) for value in best.x],
        "covariance": covariance,
    }


def theta_to_parameters(theta: list[float], five: bool) -> dict[str, float]:
    """Map the fitted vector to model parameters."""

    return {
        "a": theta[0],
        "d": theta[1],
        "c": math.exp(theta[2]),
        "b": math.exp(theta[3]),
        "g": math.exp(theta[4]) if five else 1.0,
    }


def concentration_interval(
    y_mean: float, n_wells: int, curve: dict[str, Any], level: float = 0.95
) -> tuple[float, float] | None:
    """Delta-method interval for a back-calculated concentration, on the log scale.

    Combines curve-parameter uncertainty with the read-out noise the weighted model
    implies for a mean of n_wells. Returns None when it cannot be evaluated.
    """

    covariance = curve.get("covariance")
    if covariance is None:
        return None
    theta = np.array(curve["theta"], dtype=float)
    five = curve["model"] == "5pl"

    def log_x(vector: np.ndarray, signal: float) -> float | None:
        value = invert(signal, **theta_to_parameters(list(vector), five))
        return None if value is None or value <= 0 else math.log(value)

    center = log_x(theta, y_mean)
    if center is None:
        return None
    gradient = np.zeros_like(theta)
    for index in range(len(theta)):
        step = 1e-5 * max(abs(theta[index]), 1e-3)
        up, down = theta.copy(), theta.copy()
        up[index] += step
        down[index] -= step
        upper, lower = log_x(up, y_mean), log_x(down, y_mean)
        if upper is None or lower is None:
            return None
        gradient[index] = (upper - lower) / (2 * step)
    signal_step = 1e-5 * max(abs(y_mean), 1e-3)
    upper, lower = log_x(theta, y_mean + signal_step), log_x(theta, y_mean - signal_step)
    if upper is None or lower is None:
        return None
    slope = (upper - lower) / (2 * signal_step)
    signal = max(abs(y_mean), curve["signal_floor"])
    unit = {"1/y2": signal, "1/y": math.sqrt(signal), "none": 1.0}[curve["weight"]]
    read_variance = curve["residual_variance"] * unit**2 / n_wells
    variance = float(gradient @ np.array(covariance) @ gradient) + slope**2 * read_variance
    if variance < 0:
        return None
    margin = float(student_t.ppf(0.5 + level / 2, curve["df"])) * math.sqrt(variance)
    return math.exp(center - margin), math.exp(center + margin)


def coefficient_of_variation(values: list[float]) -> float | None:
    """Percent CV of replicate readings."""

    if len(values) < 2:
        return None
    mean = statistics.fmean(values)
    return None if mean == 0 else 100.0 * statistics.stdev(values) / abs(mean)


def choose_blank(wells, no_subtract, warnings) -> float:
    """Pick the blank. Zero standards count only on a rising curve."""

    if no_subtract:
        return 0.0
    blanks = [well["od"] for well in wells if well["kind"] == "blank"]
    if blanks:
        return statistics.fmean(blanks)
    zeros = [w["od"] for w in wells if w["kind"] == "standard" and w["conc"] == 0]
    if not zeros:
        warnings.append("no blank wells and no zero standard: OD used without blank subtraction")
        return 0.0
    positive = [w for w in wells if w["kind"] == "standard" and w["conc"] > 0]
    if not positive:
        raise InputError("no standards above zero concentration")
    top_level = max(w["conc"] for w in positive)
    top_mean = statistics.fmean(w["od"] for w in positive if w["conc"] == top_level)
    if statistics.fmean(zeros) < top_mean:
        return statistics.fmean(zeros)
    warnings.append(
        "zero standard reads above the top standard, as in a competitive assay, "
        "so it was not used as a blank"
    )
    return 0.0


def working_range(levels: list[dict[str, Any]], max_recovery_error: float, max_cv: float):
    """Longest contiguous run of standard levels that pass recovery and CV."""

    best: list[dict[str, Any]] = []
    run: list[dict[str, Any]] = []
    for level in levels:
        recovery = level["recovery_percent"]
        cv = level["od_cv_percent"]
        passes = (
            recovery is not None
            and abs(recovery - 100.0) <= max_recovery_error
            and (cv is None or cv <= max_cv)
        )
        level["in_working_range"] = False
        if passes:
            run.append(level)
            if len(run) > len(best):
                best = list(run)
        else:
            run = []
    for level in best:
        level["in_working_range"] = True
    return (best[0]["nominal"], best[-1]["nominal"]) if best else (None, None)


def assess_unknowns(rows, curve, lloq, uloq, max_cv) -> list[dict[str, Any]]:
    """Back-calculate samples or QCs per (id, dilution) group."""

    parameters = {key: curve[key] for key in ("a", "d", "c", "b", "g")}
    grouped: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(row["id"], row["dilution"])].append(row)

    results = []
    for (sample_id, dilution), members in sorted(grouped.items()):
        net = [member["net_od"] for member in members]
        neat = invert(statistics.fmean(net), **parameters)
        cv = coefficient_of_variation([member["od"] for member in members])
        entry: dict[str, Any] = {
            "kind": members[0]["kind"],
            "id": sample_id,
            "dilution": dilution,
            "n_wells": len(net),
            "mean_net_od": round(statistics.fmean(net), 4),
            "od_cv_percent": None if cv is None else round(cv, 2),
            "well_concentration": sig(neat),
            "concentration": None,
            "censored": None,
            "status": "ok",
            "flags": [],
        }
        if lloq is None:
            entry["status"] = "no_working_range"
        elif neat is None:
            low_side = (statistics.fmean(net) <= min(curve["a"], curve["d"])) == (
                curve["direction"] == "sandwich"
            )
            entry["status"] = "below_range" if low_side else "above_range"
        elif neat < lloq:
            entry["status"] = "below_range"
        elif neat > uloq:
            entry["status"] = "above_range"
        if entry["status"] == "below_range":
            entry["censored"] = f"<{sig(lloq * dilution)}"
        elif entry["status"] == "above_range":
            entry["censored"] = f">{sig(uloq * dilution)}"
            entry["flags"].append("re-assay at a higher dilution")
        elif entry["status"] == "ok":
            entry["concentration"] = sig(neat * dilution)
            interval = concentration_interval(statistics.fmean(net), len(net), curve)
            if interval is not None:
                entry["concentration_ci"] = [sig(bound * dilution, 4) for bound in interval]
        if cv is not None and cv > max_cv:
            entry["flags"].append(f"OD CV {cv:.1f}% > {max_cv}%")
        if len(net) < 2:
            entry["flags"].append("single well, so precision is unknown")
        nominal = members[0]["conc"]
        if members[0]["kind"] == "qc" and nominal:
            entry["nominal"] = nominal
            value = entry["concentration"]
            entry["recovery_percent"] = None if value is None else round(100.0 * value / nominal, 2)
        results.append(entry)
    return results


def check_dilution_linearity(results, max_cv, warnings) -> list[dict[str, Any]]:
    """Dilution-corrected agreement for samples read at two or more dilutions."""

    by_id: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for entry in results:
        if entry["kind"] == "sample" and entry["concentration"] is not None:
            by_id[entry["id"]].append(entry)
    findings = []
    for sample_id, entries in sorted(by_id.items()):
        if len(entries) < 2:
            continue
        values = [entry["concentration"] for entry in entries]
        cv = coefficient_of_variation(values)
        findings.append(
            {
                "id": sample_id,
                "dilutions": [entry["dilution"] for entry in entries],
                "concentrations": values,
                "cv_percent": round(cv, 2),
            }
        )
        if cv > max_cv:
            warnings.append(
                f"{sample_id}: dilution-corrected values differ by CV {cv:.1f}%, which "
                "suggests a matrix effect or non-parallelism"
            )
    return findings


def analyze(
    wells: list[dict[str, Any]],
    *,
    model: str = "4pl",
    weight: str = "1/y2",
    exclude_levels: tuple[float, ...] = (),
    no_blank_subtract: bool = False,
    max_recovery_error: float = 20.0,
    max_cv: float = 20.0,
    max_od: float = 3.5,
) -> dict[str, Any]:
    """Fit the curve, derive the working range, and assess unknowns."""

    if model not in MIN_LEVELS:
        raise InputError(f"unknown model {model!r}")
    if weight not in WEIGHTS:
        raise InputError(f"unknown weight {weight!r}")
    warnings: list[str] = []
    blank = choose_blank(wells, no_blank_subtract, warnings)
    rows = [{**well, "net_od": well["od"] - blank} for well in wells]

    standards = [row for row in rows if row["kind"] == "standard" and row["conc"] > 0]
    excluded = [level for level in exclude_levels if any(row["conc"] == level for row in standards)]
    unknown_levels = sorted(set(exclude_levels) - set(excluded))
    if unknown_levels:
        raise InputError(f"--exclude-level {unknown_levels} matches no standard")
    if excluded:
        warnings.append(
            f"standard level(s) {sorted(excluded)} excluded by the analyst; "
            "state the reason in the report"
        )
    fit_rows = [row for row in standards if row["conc"] not in excluded]
    level_count = len({row["conc"] for row in fit_rows})
    if level_count < HARD_MIN_LEVELS[model]:
        raise InputError(
            f"{level_count} standard levels cannot support {model}; "
            f"at least {HARD_MIN_LEVELS[model]} are needed"
        )
    if level_count < MIN_LEVELS[model]:
        warnings.append(
            f"{level_count} standard levels leave few degrees of freedom for {model}; "
            f"{MIN_LEVELS[model]} or more are preferable"
        )

    x = np.array([row["conc"] for row in fit_rows], dtype=float)
    y = np.array([row["net_od"] for row in fit_rows], dtype=float)
    curve = fit_curve(x, y, model, weight)
    parameters = {key: curve[key] for key in ("a", "d", "c", "b", "g")}

    grouped: dict[float, list[dict[str, Any]]] = defaultdict(list)
    for row in standards:
        grouped[row["conc"]].append(row)
    levels = []
    for nominal in sorted(grouped):
        members = grouped[nominal]
        estimate = invert(statistics.fmean(m["net_od"] for m in members), **parameters)
        cv = coefficient_of_variation([m["od"] for m in members])
        levels.append(
            {
                "nominal": nominal,
                "n_wells": len(members),
                "mean_net_od": round(statistics.fmean(m["net_od"] for m in members), 4),
                "od_cv_percent": None if cv is None else round(cv, 2),
                "back_calculated": sig(estimate),
                "recovery_percent": None if estimate is None else round(100 * estimate / nominal, 2),
                "excluded": nominal in excluded,
            }
        )
    lloq, uloq = working_range(
        [level for level in levels if not level["excluded"]], max_recovery_error, max_cv
    )
    if lloq is None:
        warnings.append("no standard level passed recovery and CV, so no working range exists")

    if curve["direction"] == "sandwich" and len(levels) >= 2:
        if levels[-1]["mean_net_od"] < levels[-2]["mean_net_od"]:
            warnings.append(
                "top standard reads below the one under it: possible hook effect or saturation"
            )
    hot = sorted({row["id"] for row in rows if row["od"] > max_od})
    if hot:
        warnings.append(f"OD above {max_od} in {', '.join(hot)}: check the reader's linear range")

    unknowns = assess_unknowns(
        [row for row in rows if row["kind"] in {"sample", "qc"}], curve, lloq, uloq, max_cv
    )
    linearity = check_dilution_linearity(unknowns, max_cv, warnings)
    for entry in unknowns:
        if entry["kind"] == "qc":
            recovery = entry.get("recovery_percent")
            if recovery is None or abs(recovery - 100.0) > max_recovery_error:
                warnings.append(f"QC {entry['id']} recovery is {recovery}; the run may be invalid")

    return {
        "blank_od": round(blank, 4),
        "curve": {
            key: (sig(value, 8) if isinstance(value, float) else value)
            for key, value in curve.items()
            if key not in {"theta", "covariance", "signal_floor", "residual_variance"}
        },
        "lloq": lloq,
        "uloq": uloq,
        "standards": levels,
        "unknowns": unknowns,
        "dilution_linearity": linearity,
        "warnings": warnings,
    }


def compare_fits(wells: list[dict[str, Any]], **options: Any) -> list[dict[str, Any]]:
    """Summarize standard recovery for each model and weight; samples are not shown."""

    rows = []
    for model in ("4pl", "5pl"):
        for weight in WEIGHTS:
            try:
                result = analyze(wells, model=model, weight=weight, **options)
            except InputError as exc:
                rows.append({"model": model, "weight": weight, "error": str(exc)})
                continue
            errors = [
                abs(level["recovery_percent"] - 100.0)
                for level in result["standards"]
                if level["recovery_percent"] is not None and not level["excluded"]
            ]
            rows.append(
                {
                    "model": model,
                    "weight": weight,
                    "mean_recovery_error_percent": round(statistics.fmean(errors), 2) if errors else None,
                    "worst_recovery_error_percent": round(max(errors), 2) if errors else None,
                    "lloq": result["lloq"],
                    "uloq": result["uloq"],
                }
            )
    return rows


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="CSV of wells from one plate")
    parser.add_argument("--model", choices=sorted(MIN_LEVELS), default="4pl")
    parser.add_argument("--weight", choices=WEIGHTS, default="1/y2")
    parser.add_argument("--exclude-level", type=float, action="append", default=[], help="standard level to leave out; justify it")
    parser.add_argument("--no-blank-subtract", action="store_true", help="fit raw OD")
    parser.add_argument("--max-recovery-error", type=float, default=20.0, help="allowed percent error at a standard level")
    parser.add_argument("--max-cv", type=float, default=20.0, help="allowed replicate percent CV")
    parser.add_argument("--max-od", type=float, default=3.5, help="flag OD above this")
    parser.add_argument("--compare", action="store_true", help="compare models and weights on standards only, then stop")
    parser.add_argument("--output-json", help="write full results here")
    parser.add_argument("--output-csv", help="write the sample table here")
    args = parser.parse_args(argv)

    if args.compare:
        try:
            wells = [well for well in read_wells(args.input) if well["kind"] in {"standard", "blank"}]
            rows = compare_fits(
                wells,
                exclude_levels=tuple(args.exclude_level),
                no_blank_subtract=args.no_blank_subtract,
                max_recovery_error=args.max_recovery_error,
                max_cv=args.max_cv,
                max_od=args.max_od,
            )
        except (InputError, OSError) as exc:
            print(f"[FAIL] {exc}", file=sys.stderr)
            return 2
        for row in rows:
            if "error" in row:
                print(f"  -> {row['model']} {row['weight']}: not possible ({row['error']})")
            else:
                print(
                    f"  -> {row['model']} {row['weight']}: mean error {row['mean_recovery_error_percent']}%, "
                    f"worst {row['worst_recovery_error_percent']}%, range {row['lloq']} to {row['uloq']}"
                )
        print("[OK] choose a model and weight here, before looking at sample results")
        return 0

    try:
        result = analyze(
            read_wells(args.input),
            model=args.model,
            weight=args.weight,
            exclude_levels=tuple(args.exclude_level),
            no_blank_subtract=args.no_blank_subtract,
            max_recovery_error=args.max_recovery_error,
            max_cv=args.max_cv,
            max_od=args.max_od,
        )
    except (InputError, OSError) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    curve = result["curve"]
    print(
        f"[OK] {curve['model']} ({curve['weight']}, {curve['direction']}): "
        f"midpoint {curve['c']:.4g}, slope {curve['b']:.3f}, blank {result['blank_od']}"
    )
    print(f"[OK] working range {result['lloq']} to {result['uloq']} in well concentration")
    for level in result["standards"]:
        print(
            f"  -> standard {level['nominal']}: recovery {level['recovery_percent']}%, "
            f"CV {level['od_cv_percent']}%"
        )
    for entry in result["unknowns"]:
        value = entry["concentration"] if entry["concentration"] is not None else entry["censored"]
        print(f"  -> {entry['id']} (1:{entry['dilution']:g}): {value} [{entry['status']}]")
    for warning in result["warnings"]:
        print(f"[WARN] {warning}")

    if args.output_json:
        Path(args.output_json).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(f"[OK] wrote {args.output_json}")
    if args.output_csv:
        columns = [
            "kind", "id", "dilution", "n_wells", "mean_net_od", "od_cv_percent",
            "concentration", "censored", "status",
        ]
        with open(args.output_csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(result["unknowns"])
        print(f"[OK] wrote {args.output_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
