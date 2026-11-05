#!/usr/bin/env python3
"""Case-control association for biallelic candidate variants, with genotype QC.

Two inputs:
  --input table.csv   one row per person, a status column, one column per variant
  --counts            --cases HOM,HET,HOM --controls HOM,HET,HOM for one variant

Genetic models are coded on the effect allele. One model is primary and carries
the multiple-testing correction; the others are returned as exploratory, because
choosing the best of several models after seeing them inflates false positives.
Hardy-Weinberg is tested in controls and reported, and never removes a variant.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

from _genetics import (
    InputError,
    chi2_sf,
    fisher_exact,
    hardy_weinberg,
    odds_ratio,
    parse_genotype,
    pearson_2x2,
    pearson_2x3,
    per_allele_logistic,
    read_text,
    safe_console,
    trend_chi_square,
)

MODELS = ("additive", "allelic", "dominant", "recessive", "codominant")
HWE_FLAG = 0.05
SMALL_EXPECTED = 5.0
MISSINGNESS_GAP = 0.05


def two_by_two(label, a, b, c, d) -> dict[str, Any]:
    """Test and effect size for one 2x2 table."""

    statistic = pearson_2x2(a, b, c, d)
    total = a + b + c + d
    expected = (
        min((a + b), (c + d)) * min((a + c), (b + d)) / total if total else 0.0
    )
    use_exact = expected < SMALL_EXPECTED
    fisher = fisher_exact(a, b, c, d) if total else None
    effect = odds_ratio(a, b, c, d)
    p_value = fisher if use_exact or statistic is None else chi2_sf(statistic, 1)
    return {
        "model": label,
        "table": [[a, b], [c, d]],
        "test": "fisher_exact" if use_exact or statistic is None else "chi_square",
        "statistic": statistic,
        "p_value": p_value,
        "fisher_p": fisher,
        **effect,
    }


def association_tests(cases: tuple[int, int, int], controls: tuple[int, int, int]) -> dict[str, Any]:
    """All models for counts ordered by effect-allele dose 0, 1, 2."""

    r0, r1, r2 = cases
    s0, s1, s2 = controls
    allelic = two_by_two("allelic", 2 * r2 + r1, 2 * r0 + r1, 2 * s2 + s1, 2 * s0 + s1)
    allelic["unit"] = "per allele count (assumes Hardy-Weinberg in the pooled sample)"
    dominant = two_by_two("dominant", r1 + r2, r0, s1 + s2, s0)
    recessive = two_by_two("recessive", r2, r0 + r1, s2, s0 + s1)

    trend = trend_chi_square(cases, controls)
    fit = per_allele_logistic(cases, controls)
    additive = {
        "model": "additive",
        "test": "cochran_armitage_trend",
        "statistic": trend,
        "p_value": None if trend is None else chi2_sf(trend, 1),
        "or": None if fit is None else fit["or"],
        "ci_low": None if fit is None else fit["ci_low"],
        "ci_high": None if fit is None else fit["ci_high"],
        "unit": "per effect allele",
        "haldane_corrected": False,
    }

    overall = pearson_2x3(cases, controls)
    codominant = {
        "model": "codominant",
        "test": "pearson_2df",
        "statistic": overall,
        "p_value": None if overall is None else chi2_sf(overall, 2),
        "het_vs_ref": odds_ratio(r1, r0, s1, s0),
        "hom_vs_ref": odds_ratio(r2, r0, s2, s0),
        "or": None,
        "ci_low": None,
        "ci_high": None,
    }
    return {
        "additive": additive,
        "allelic": allelic,
        "dominant": dominant,
        "recessive": recessive,
        "codominant": codominant,
    }


def analyze_counts(
    variant: str,
    effect_allele: str,
    other_allele: str,
    cases: tuple[int, int, int],
    controls: tuple[int, int, int],
    primary_model: str,
    n_tests: int,
    extra_warnings: list[str] | None = None,
) -> dict[str, Any]:
    """Assemble QC and tests for one variant, given counts by effect-allele dose."""

    warnings = list(extra_warnings or [])
    if sum(cases) == 0 or sum(controls) == 0:
        raise InputError(f"{variant}: needs both cases and controls")

    hwe_controls = hardy_weinberg(*controls)
    hwe_cases = hardy_weinberg(*cases)
    if hwe_controls["exact_p"] < HWE_FLAG:
        reason = (
            "a heterozygote deficit suggests allele dropout, mistyping, or inbreeding"
            if hwe_controls["direction"] == "heterozygote_deficit"
            else "a heterozygote excess suggests paralogous sequence or merged clusters"
        )
        warnings.append(
            f"{variant}: controls deviate from Hardy-Weinberg "
            f"(exact p={hwe_controls['exact_p']:.3g}, F={hwe_controls['inbreeding_f']}); "
            f"{reason}. Recheck the genotyping before trusting any result."
        )
    if hwe_controls["min_expected_count"] < SMALL_EXPECTED:
        warnings.append(
            f"{variant}: an expected genotype count is below {SMALL_EXPECTED:g}, so read the "
            "exact Hardy-Weinberg p and not the chi-square p"
        )

    tests = association_tests(cases, controls)
    primary = tests[primary_model]
    p_value = primary["p_value"]
    corrected = None if p_value is None else min(1.0, p_value * n_tests)
    if primary_model == "allelic" and hwe_controls["exact_p"] < HWE_FLAG:
        warnings.append(
            f"{variant}: the allelic test is unreliable without Hardy-Weinberg; "
            "use the additive trend test"
        )
    if primary_model in {"additive", "codominant"}:
        small = min(cases[i] + controls[i] for i in range(3))
        if small < SMALL_EXPECTED:
            warnings.append(
                f"{variant}: a genotype class has {small} people, so the asymptotic "
                f"{primary_model} p-value is approximate"
            )
    if any(row.get("haldane_corrected") for row in tests.values()):
        warnings.append(f"{variant}: a zero cell needed a 0.5 correction, so its OR is unstable")

    case_alleles = 2 * sum(cases)
    control_alleles = 2 * sum(controls)
    return {
        "variant": variant,
        "effect_allele": effect_allele,
        "other_allele": other_allele,
        "cases": {"counts_by_dose": list(cases), "n": sum(cases)},
        "controls": {"counts_by_dose": list(controls), "n": sum(controls)},
        "effect_allele_frequency": {
            "cases": round((2 * cases[2] + cases[1]) / case_alleles, 4),
            "controls": round((2 * controls[2] + controls[1]) / control_alleles, 4),
        },
        "hwe_controls": hwe_controls,
        "hwe_cases": hwe_cases,
        "primary_model": primary_model,
        "primary": {
            **primary,
            "p_bonferroni": corrected,
            "n_tests_corrected": n_tests,
        },
        "exploratory": {name: row for name, row in tests.items() if name != primary_model},
        "warnings": warnings,
    }


def read_table(path: str | Path, status_column: str) -> tuple[list[str], list[dict[str, Any]]]:
    """Read a one-row-per-person table."""

    source = Path(path)
    text = read_text(source)
    delimiter = max(",;\t", key=text.splitlines()[0].count)
    reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
    columns = [name.strip() for name in reader.fieldnames or []]
    reader.fieldnames = columns
    if status_column not in columns:
        raise InputError(f"status column {status_column!r} not found in {columns}")
    rows = []
    for line, row in enumerate(reader, start=2):
        if None in row:
            raise InputError(f"line {line}: more fields than the header")
        cleaned = {key: (value or "").strip() for key, value in row.items()}
        if any(cleaned.values()):
            cleaned["_line"] = line
            rows.append(cleaned)
    if not rows:
        raise InputError(f"input file has no data rows: {source}")
    return columns, rows


def choose_effect_allele(alleles, control_counts, requested, variant) -> tuple[str, str]:
    """Use the requested allele, otherwise the minor allele in controls."""

    if requested:
        if requested not in alleles:
            raise InputError(f"{variant}: effect allele {requested!r} not in {sorted(alleles)}")
        other = next(allele for allele in sorted(alleles) if allele != requested)
        return requested, other
    first, second = sorted(alleles)
    if control_counts[second] < control_counts[first]:
        return second, first
    return first, second


def analyze_table(
    columns: list[str],
    rows: list[dict[str, Any]],
    *,
    status_column: str,
    case_label: str,
    control_label: str,
    variants: list[str] | None,
    id_column: str | None,
    effect_alleles: dict[str, str],
    primary_model: str,
) -> dict[str, Any]:
    """Count genotypes per variant and run the association tests."""

    labels = {}
    for row in rows:
        status = row[status_column]
        if status == case_label:
            labels[row["_line"]] = "case"
        elif status == control_label:
            labels[row["_line"]] = "control"
        else:
            raise InputError(
                f"line {row['_line']}: status {status!r} is neither "
                f"{case_label!r} nor {control_label!r}"
            )

    skip = {status_column, "_line"} | ({id_column} if id_column else set())
    chosen = variants or [name for name in columns if name not in skip]
    absent = [name for name in chosen if name not in columns]
    if absent:
        raise InputError(f"variant column(s) not found: {', '.join(absent)}")
    if not chosen:
        raise InputError("no variant columns to test")

    results = []
    skipped = []
    dosages: dict[str, dict[int, int | None]] = {}
    for variant in chosen:
        genotypes = {}
        for row in rows:
            try:
                genotypes[row["_line"]] = parse_genotype(row[variant])
            except InputError as exc:
                raise InputError(f"line {row['_line']}, {variant}: {exc}") from exc
        called = {line: call for line, call in genotypes.items() if call is not None}
        alleles = {allele for call in called.values() for allele in call}
        if len(alleles) != 2:
            skipped.append(
                {"variant": variant, "reason": f"{len(alleles)} allele(s) seen; exactly 2 are supported"}
            )
            continue
        control_allele_counts = Counter(
            allele
            for line, call in called.items()
            if labels[line] == "control"
            for allele in call
        )
        effect, other = choose_effect_allele(
            alleles, control_allele_counts, effect_alleles.get(variant), variant
        )
        counts = {"case": [0, 0, 0], "control": [0, 0, 0]}
        dosage: dict[int, int | None] = {line: None for line in genotypes}
        for line, call in called.items():
            dose = call.count(effect)
            counts[labels[line]][dose] += 1
            dosage[line] = dose
        dosages[variant] = dosage

        total_by_status = Counter(labels.values())
        missing = {
            status: 1.0 - sum(counts[status]) / total_by_status[status]
            for status in ("case", "control")
        }
        extra = []
        if abs(missing["case"] - missing["control"]) > MISSINGNESS_GAP:
            extra.append(
                f"{variant}: missing calls differ by status "
                f"({missing['case']:.1%} of cases, {missing['control']:.1%} of controls), "
                "which can create a false association"
            )
        result = {
            "variant": variant,
            "missing_rate": {key: round(value, 4) for key, value in missing.items()},
            "counts": counts,
            "effect": effect,
            "other": other,
            "extra": extra,
        }
        results.append(result)

    n_tests = len(results)
    final = []
    for item in results:
        analysis = analyze_counts(
            item["variant"],
            item["effect"],
            item["other"],
            tuple(item["counts"]["case"]),
            tuple(item["counts"]["control"]),
            primary_model,
            n_tests,
            item["extra"],
        )
        analysis["missing_rate"] = item["missing_rate"]
        final.append(analysis)
    return {"variants": final, "skipped": skipped, "dosages": dosages}


def write_dosages(path, rows, dosages, status_column, id_column) -> None:
    """Write effect-allele doses so covariates can be modelled elsewhere."""

    names = list(dosages)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        lead = [id_column] if id_column else []
        writer.writerow(lead + [status_column] + names)
        for row in rows:
            values = ["" if dosages[name][row["_line"]] is None else dosages[name][row["_line"]] for name in names]
            writer.writerow(([row[id_column]] if id_column else []) + [row[status_column]] + values)


def parse_counts(text: str, flag: str) -> tuple[int, int, int]:
    """Parse three genotype counts."""

    try:
        counts = tuple(int(part) for part in text.split(","))
    except ValueError as exc:
        raise InputError(f"{flag} must be integers, got {text!r}") from exc
    if len(counts) != 3 or min(counts) < 0:
        raise InputError(f"{flag} needs three non-negative counts, got {text!r}")
    return counts  # type: ignore[return-value]


def report(result: dict[str, Any]) -> None:
    """Print a short ASCII summary."""

    for variant in result["variants"]:
        primary = variant["primary"]
        p_text = "n/a" if primary["p_value"] is None else f"{primary['p_value']:.4g}"
        corrected = "n/a" if primary["p_bonferroni"] is None else f"{primary['p_bonferroni']:.4g}"
        effect = "n/a" if primary["or"] is None else (
            f"OR {primary['or']:.3f} ({primary['ci_low']:.3f}-{primary['ci_high']:.3f})"
        )
        marker = "[WARN]" if variant["warnings"] else "[OK]"
        print(
            f"{marker} {variant['variant']} effect allele {variant['effect_allele']} "
            f"[{variant['primary_model']}]: {effect}, p {p_text}, Bonferroni {corrected}, "
            f"HWE controls p {variant['hwe_controls']['exact_p']:.3g}"
        )
        for warning in variant["warnings"]:
            print(f"       -> {warning}")
    for entry in result["skipped"]:
        print(f"[WARN] skipped {entry['variant']}: {entry['reason']}")


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", help="CSV with one row per person")
    parser.add_argument("--status-column", default="status")
    parser.add_argument("--case-label", default="case")
    parser.add_argument("--control-label", default="control")
    parser.add_argument("--id-column", help="column to leave out of the variants")
    parser.add_argument("--variants", help="comma-separated variant columns; default is all others")
    parser.add_argument("--effect-allele", action="append", default=[], help="VARIANT=ALLELE, repeatable")
    parser.add_argument("--model", choices=MODELS, default="additive", help="pre-specified primary model")
    parser.add_argument("--dosage-csv", help="write effect-allele doses for covariate models")
    parser.add_argument("--cases", help="dose 0,1,2 counts in cases, for one variant")
    parser.add_argument("--controls", help="dose 0,1,2 counts in controls, for one variant")
    parser.add_argument("--variant-name", default="variant")
    parser.add_argument("--effect-name", default="effect", help="label for the dose-counted allele")
    parser.add_argument("--other-name", default="other")
    parser.add_argument("--n-tests", type=int, default=1, help="variants tested, for --counts mode")
    parser.add_argument("--output-json", help="write full results here")
    args = parser.parse_args(argv)

    try:
        if args.input and (args.cases or args.controls):
            raise InputError("give --input or --cases/--controls, not both")
        if args.input:
            effect_alleles = {}
            for item in args.effect_allele:
                name, separator, allele = item.partition("=")
                if not separator:
                    raise InputError(f"--effect-allele needs VARIANT=ALLELE, got {item!r}")
                effect_alleles[name] = allele
            columns, rows = read_table(args.input, args.status_column)
            result = analyze_table(
                columns,
                rows,
                status_column=args.status_column,
                case_label=args.case_label,
                control_label=args.control_label,
                variants=args.variants.split(",") if args.variants else None,
                id_column=args.id_column,
                effect_alleles=effect_alleles,
                primary_model=args.model,
            )
            dosages = result.pop("dosages")
            if args.dosage_csv:
                write_dosages(args.dosage_csv, rows, dosages, args.status_column, args.id_column)
        elif args.cases and args.controls:
            analysis = analyze_counts(
                args.variant_name,
                args.effect_name,
                args.other_name,
                parse_counts(args.cases, "--cases"),
                parse_counts(args.controls, "--controls"),
                args.model,
                args.n_tests,
            )
            result = {"variants": [analysis], "skipped": []}
        else:
            raise InputError("give --input, or both --cases and --controls")
    except (InputError, OSError) as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    report(result)
    if args.dosage_csv:
        print(f"[OK] wrote {args.dosage_csv}")
    if args.output_json:
        Path(args.output_json).write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(f"[OK] wrote {args.output_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
