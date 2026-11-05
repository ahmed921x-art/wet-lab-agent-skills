#!/usr/bin/env python3
"""Tests for candidate-gene-association against exact and identity-based oracles."""

from __future__ import annotations

import csv
import json
import math
import os
import random
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest
import skill_contract

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "candidate-gene-association"
SCRIPTS = SKILL_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import _genetics  # noqa: E402
import candidate_association  # noqa: E402
import candidate_meta  # noqa: E402
import candidate_power  # noqa: E402
import candidate_report  # noqa: E402
import hwe_test  # noqa: E402

CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)


def exact_oracle(hom_a: int, het: int, hom_b: int) -> Fraction:
    """Direct enumeration in exact arithmetic, independent of the log-gamma code."""

    n = hom_a + het + hom_b
    n_a = 2 * hom_a + het
    n_b = 2 * hom_b + het

    def probability(hets: int) -> Fraction:
        aa = (n_a - hets) // 2
        bb = (n_b - hets) // 2
        return Fraction(
            2**hets
            * math.factorial(n)
            * math.factorial(n_a)
            * math.factorial(n_b),
            math.factorial(aa)
            * math.factorial(hets)
            * math.factorial(bb)
            * math.factorial(2 * n),
        )

    support = range(n_a % 2, min(n_a, n_b) + 1, 2)
    observed = probability(het)
    assert sum(probability(h) for h in support) == 1
    return sum((probability(h) for h in support if probability(h) <= observed), Fraction(0))


@pytest.mark.parametrize(
    "counts",
    [(12, 40, 8), (20, 10, 20), (1, 0, 9), (0, 3, 0), (5, 5, 5), (21, 13, 16), (30, 2, 1)],
)
def test_exact_hwe_matches_exact_arithmetic(counts):
    assert _genetics.hwe_exact_p(*counts) == pytest.approx(float(exact_oracle(*counts)), rel=1e-9)


def test_exact_hwe_of_a_monomorphic_variant_is_one():
    assert _genetics.hwe_exact_p(50, 0, 0) == pytest.approx(1.0)


def test_chi_square_by_hand():
    result = _genetics.hardy_weinberg(10, 20, 30)
    p = 40 / 120
    expected = [60 * p * p, 120 * p * (1 - p), 60 * (1 - p) ** 2]
    hand = sum((o - e) ** 2 / e for o, e in zip((10, 20, 30), expected))
    assert result["chi_square"] == pytest.approx(hand, abs=1e-4)
    assert result["inbreeding_f"] == pytest.approx(1 - (20 / 60) / (2 * p * (1 - p)), abs=1e-4)
    assert result["direction"] == "heterozygote_deficit"


def test_chi2_tail_matches_table_values():
    assert _genetics.chi2_sf(3.841459, 1) == pytest.approx(0.05, abs=1e-6)
    assert _genetics.chi2_sf(5.991465, 2) == pytest.approx(0.05, abs=1e-6)


def test_fisher_exact_on_the_tea_tasting_table():
    assert _genetics.fisher_exact(3, 1, 1, 3) == pytest.approx(34 / 70, rel=1e-12)


def test_odds_ratio_and_woolf_interval_by_hand():
    result = _genetics.odds_ratio(20, 10, 10, 20)
    se = math.sqrt(1 / 20 + 1 / 10 + 1 / 10 + 1 / 20)
    assert result["or"] == pytest.approx(4.0)
    assert result["ci_low"] == pytest.approx(math.exp(math.log(4) - 1.959963985 * se))
    assert result["haldane_corrected"] is False


def test_zero_cell_uses_the_haldane_correction():
    result = _genetics.odds_ratio(10, 0, 5, 5)
    assert result["or"] == pytest.approx((10.5 * 5.5) / (0.5 * 5.5))
    assert result["haldane_corrected"] is True


def test_trend_statistic_equals_n_r_squared():
    cases, controls = (12, 34, 32), (22, 40, 17)
    dose, status = [], []
    for value, group in ((1, cases), (0, controls)):
        for score, count in enumerate(group):
            dose += [score] * count
            status += [value] * count
    n = len(dose)
    mean_x, mean_y = sum(dose) / n, sum(status) / n
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(dose, status))
    sxx = sum((x - mean_x) ** 2 for x in dose)
    syy = sum((y - mean_y) ** 2 for y in status)
    assert _genetics.trend_chi_square(cases, controls) == pytest.approx(n * sxy**2 / (sxx * syy))


def test_logistic_fit_solves_the_score_equations():
    cases, controls = (12, 34, 32), (22, 40, 17)
    fit = _genetics.per_allele_logistic(cases, controls)
    score_a = score_b = 0.0
    for dose, (r, s) in enumerate(zip(cases, controls)):
        pi = 1 / (1 + math.exp(-(fit["alpha"] + fit["beta"] * dose)))
        score_a += r - (r + s) * pi
        score_b += dose * (r - (r + s) * pi)
    assert score_a == pytest.approx(0.0, abs=1e-6)
    assert score_b == pytest.approx(0.0, abs=1e-6)
    assert fit["ci_low"] < fit["or"] < fit["ci_high"]


def test_separated_logistic_returns_none():
    assert _genetics.per_allele_logistic((0, 0, 20), (20, 0, 0)) is None


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("A/G", ("A", "G")),
        ("G|A", ("A", "G")),
        ("AG", ("A", "G")),
        ("GA", ("A", "G")),
        ("II", ("I", "I")),
        ("D/I", ("D", "I")),
        ("X+/X-", ("X+", "X-")),
        ("X-X+", ("X+", "X-")),
        ("NA", None),
        ("0/0", None),
        ("./.", None),
        ("NN", None),
        ("", None),
    ],
)
def test_parse_genotype(raw, expected):
    assert _genetics.parse_genotype(raw) == expected


@pytest.mark.parametrize("raw", ["A", "A/G/T", "/G"])
def test_parse_genotype_refuses_what_it_cannot_read(raw):
    with pytest.raises(_genetics.InputError):
        _genetics.parse_genotype(raw)


def test_models_count_the_effect_allele_correctly():
    tests = candidate_association.association_tests((10, 20, 30), (30, 20, 10))
    assert tests["dominant"]["table"] == [[50, 10], [30, 30]]
    assert tests["recessive"]["table"] == [[30, 30], [10, 50]]
    assert tests["allelic"]["table"] == [[80, 40], [40, 80]]
    assert tests["allelic"]["or"] == pytest.approx(4.0)
    assert tests["additive"]["or"] > 1


def test_small_tables_use_fisher():
    row = candidate_association.two_by_two("dominant", 2, 8, 1, 9)
    assert row["test"] == "fisher_exact"
    assert row["p_value"] == row["fisher_p"]


def counts_analysis(cases, controls, model="additive", n_tests=1):
    return candidate_association.analyze_counts("v", "T", "C", cases, controls, model, n_tests)


def test_bonferroni_counts_variants_and_caps_at_one():
    result = counts_analysis((12, 34, 32), (22, 40, 17), n_tests=2)
    primary = result["primary"]
    assert primary["p_bonferroni"] == pytest.approx(min(1.0, primary["p_value"] * 2))
    capped = counts_analysis((22, 40, 17), (22, 40, 17), n_tests=5)
    assert capped["primary"]["p_bonferroni"] == 1.0


def test_only_the_primary_model_is_corrected_and_the_rest_are_exploratory():
    result = counts_analysis((12, 34, 32), (22, 40, 17), model="dominant")
    assert result["primary"]["model"] == "dominant"
    assert set(result["exploratory"]) == {"additive", "allelic", "recessive", "codominant"}
    assert all("p_bonferroni" not in row for row in result["exploratory"].values())


def test_heterozygote_deficit_in_controls_warns_about_mistyping():
    result = counts_analysis((30, 20, 30), (40, 10, 40))
    assert result["hwe_controls"]["direction"] == "heterozygote_deficit"
    assert any("deviate from Hardy-Weinberg" in warning for warning in result["warnings"])
    assert any("mistyping" in warning for warning in result["warnings"])


def test_allelic_primary_warns_when_hwe_fails():
    result = counts_analysis((30, 20, 30), (40, 10, 40), model="allelic")
    assert any("allelic test is unreliable" in warning for warning in result["warnings"])


def test_deviation_in_cases_alone_does_not_warn():
    result = counts_analysis((40, 10, 40), (25, 50, 25))
    assert result["hwe_cases"]["exact_p"] < 0.05
    assert not any("deviate" in warning for warning in result["warnings"])


def test_missing_group_is_refused():
    with pytest.raises(_genetics.InputError, match="both cases and controls"):
        counts_analysis((0, 0, 0), (1, 2, 1))


def run_table(tmp_path, rows, *extra):
    path = tmp_path / "people.csv"
    with open(path, "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerows(rows)
    columns, table = candidate_association.read_table(path, "status")
    options = dict(
        status_column="status",
        case_label="case",
        control_label="control",
        variants=None,
        id_column="id",
        effect_alleles={},
        primary_model="additive",
    )
    options.update(*extra)
    return candidate_association.analyze_table(columns, table, **options)


def people(case_genotypes, control_genotypes, column="rs1"):
    rows = [["id", "status", column]]
    for status, genotypes in (("case", case_genotypes), ("control", control_genotypes)):
        for index, genotype in enumerate(genotypes):
            rows.append([f"{status}{index}", status, genotype])
    return rows


def test_default_effect_allele_is_minor_in_controls(tmp_path):
    rows = people(["AA", "AG", "GG"] * 4, ["AA"] * 8 + ["AG"] * 3 + ["GG"])
    variant = run_table(tmp_path, rows)["variants"][0]
    assert variant["effect_allele"] == "G"
    assert variant["controls"]["counts_by_dose"] == [8, 3, 1]


def test_requested_effect_allele_reverses_the_dose_counts(tmp_path):
    rows = people(["II", "ID", "DD"] * 4, ["II"] * 8 + ["ID"] * 3 + ["DD"])
    variant = run_table(tmp_path, rows, {"effect_alleles": {"rs1": "I"}})["variants"][0]
    assert variant["effect_allele"] == "I"
    assert variant["controls"]["counts_by_dose"] == [1, 3, 8]


def test_unknown_effect_allele_is_refused(tmp_path):
    rows = people(["AA", "AG"], ["AG", "GG"])
    with pytest.raises(_genetics.InputError, match="effect allele"):
        run_table(tmp_path, rows, {"effect_alleles": {"rs1": "T"}})


def test_triallelic_variant_is_skipped_with_a_reason(tmp_path):
    rows = people(["AA", "AG", "GT"], ["AA", "AG", "GG"])
    result = run_table(tmp_path, rows)
    assert result["variants"] == []
    assert "3 allele(s)" in result["skipped"][0]["reason"]


def test_differential_missingness_is_flagged(tmp_path):
    rows = people(["NA"] * 5 + ["AA", "AG", "GG"] * 5, ["AA", "AG", "GG"] * 7)
    variant = run_table(tmp_path, rows)["variants"][0]
    assert variant["missing_rate"]["case"] == pytest.approx(0.25)
    assert any("missing calls differ" in warning for warning in variant["warnings"])


def test_unknown_status_is_refused(tmp_path):
    rows = people(["AA", "AG"], ["AG", "GG"])
    rows[1][1] = "Case"
    with pytest.raises(_genetics.InputError, match="neither"):
        run_table(tmp_path, rows)


def test_hwe_script_parses_counts():
    assert hwe_test.parse_counts("1,2,3") == (1, 2, 3)
    with pytest.raises(_genetics.InputError):
        hwe_test.parse_counts("1,2")


def test_cli_table_run_with_dosage_export(tmp_path):
    output = tmp_path / "out.json"
    dosage = tmp_path / "dose.csv"
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "candidate_association.py"),
            "--input",
            str(SKILL_ROOT / "assets" / "example_genotypes.csv"),
            "--id-column",
            "id",
            "--effect-allele",
            "ACE_ID=D",
            "--output-json",
            str(output),
            "--dosage-csv",
            str(dosage),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "ACE_ID effect allele D" in done.stdout
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["variants"][0]["primary"]["n_tests_corrected"] == 2
    header = dosage.read_text(encoding="utf-8").splitlines()
    assert header[0] == "id,status,ACE_ID,SNP_2"
    assert len(header) == 161


def test_cli_counts_mode_and_hwe_json():
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "candidate_association.py"),
            "--cases",
            "12,34,32",
            "--controls",
            "22,40,17",
            "--variant-name",
            "ACE",
            "--effect-name",
            "D",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "ACE effect allele D" in done.stdout
    hwe = subprocess.run(
        [sys.executable, str(SCRIPTS / "hwe_test.py"), "--counts", "22,40,17", "--json"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert json.loads(hwe.stdout)["n"] == 79


def test_cli_refuses_mixed_inputs(tmp_path):
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "candidate_association.py"),
            "--input",
            str(tmp_path / "x.csv"),
            "--cases",
            "1,2,3",
        ],
        capture_output=True,
        text=True,
    )
    assert done.returncode == 2
    assert "[FAIL]" in done.stderr


def test_arabic_labels_survive_a_legacy_console(tmp_path):
    case, control = "\u0645\u0631\u064a\u0636", "\u0633\u0644\u064a\u0645"
    rows = [f"p{i},{case},{g}" for i, g in enumerate(["II", "ID", "DD"] * 10)]
    rows += [f"c{i},{control},{g}" for i, g in enumerate(["II", "II", "ID", "DD"] * 10)]
    path = tmp_path / "people.csv"
    path.write_text("id,status,\u062c\u064a\u0646\n" + "\n".join(rows) + "\n", encoding="utf-8")
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "candidate_association.py"),
            "--input",
            str(path),
            "--id-column",
            "id",
            "--case-label",
            case,
            "--control-label",
            control,
        ],
        capture_output=True,
        env={**os.environ, "PYTHONIOENCODING": "cp1252", "PYTHONUTF8": "0"},
    )
    assert done.returncode == 0, done.stderr
    assert b"effect allele" in done.stdout


def test_other_encodings_get_a_clear_message(tmp_path):
    path = tmp_path / "people.csv"
    path.write_bytes("id,status\n\u0639,case\n".encode("cp1256"))
    with pytest.raises(_genetics.InputError, match="CSV UTF-8"):
        candidate_association.read_table(path, "status")


# ------------------------------------------------------------------- power


def test_case_frequency_reproduces_the_odds_ratio():
    p0 = 0.3
    p1 = candidate_power.case_frequency(p0, 1.8)
    assert (p1 / (1 - p1)) / (p0 / (1 - p0)) == pytest.approx(1.8)


def test_no_effect_has_power_near_alpha():
    assert candidate_power.power(200, 200, 0.3, 1.0) == pytest.approx(0.025, abs=0.003)


def test_power_agrees_with_simulation():
    rng = random.Random(7)
    cases, controls, p0, ratio = 150, 150, 0.3, 1.5
    p1 = candidate_power.case_frequency(p0, ratio)
    hits = 0
    trials = 4000
    for _ in range(trials):
        a = sum(rng.random() < p1 for _ in range(2 * cases))
        c = sum(rng.random() < p0 for _ in range(2 * controls))
        statistic = _genetics.pearson_2x2(a, 2 * cases - a, c, 2 * controls - c)
        hits += statistic is not None and statistic > 3.841459
    assert candidate_power.power(cases, controls, p0, ratio) == pytest.approx(hits / trials, abs=0.03)


def test_required_size_is_the_smallest_that_reaches_the_target():
    count = candidate_power.required_per_group(0.3, 1.5, 0.8, 0.05, 0.0)
    assert candidate_power.power(count, count, 0.3, 1.5) >= 0.8
    assert candidate_power.power(count - 1, count - 1, 0.3, 1.5) < 0.8


def test_inbreeding_raises_the_required_size():
    plain = candidate_power.required_per_group(0.3, 1.5, 0.8, 0.05, 0.0)
    inbred = candidate_power.required_per_group(0.3, 1.5, 0.8, 0.05, 0.05)
    assert inbred == pytest.approx(plain * 1.05, abs=2)


def test_minimum_detectable_or_sits_on_the_target_both_ways():
    risk = candidate_power.minimum_detectable(150, 150, 0.3, 0.8, 0.05, 0.0, 1.0)
    protective = candidate_power.minimum_detectable(150, 150, 0.3, 0.8, 0.05, 0.0, -1.0)
    assert candidate_power.power(150, 150, 0.3, risk) == pytest.approx(0.8, abs=1e-6)
    assert protective < 1 < risk


@pytest.mark.parametrize(
    "arguments",
    [
        ["--frequency", "1.2", "--odds-ratio", "1.5"],
        ["--frequency", "0.3"],
        ["--frequency", "0.3", "--cases", "10"],
        ["--frequency", "0.3", "--odds-ratio", "1"],
    ],
)
def test_power_cli_refuses_bad_requests(arguments):
    assert candidate_power.main(arguments) == 2


# -------------------------------------------------------------------- meta


@pytest.mark.parametrize(
    ("statistic", "df", "expected"),
    [(3.841459, 1, 0.05), (5.991465, 2, 0.05), (11.0705, 5, 0.05), (18.30704, 10, 0.05), (9.34182, 10, 0.5), (0.0, 3, 1.0)],
)
def test_chi2_upper_matches_table_values(statistic, df, expected):
    assert candidate_meta.chi2_upper(statistic, df) == pytest.approx(expected, abs=1e-4)


def effects(*pairs):
    return [
        {"study": f"s{i}", "log_or": effect, "variance": variance, "corrected": False}
        for i, (effect, variance) in enumerate(pairs)
    ]


def test_fixed_effect_is_the_inverse_variance_mean():
    result = candidate_meta.pool(effects((0.2, 0.04), (0.6, 0.01)))
    expected = (0.2 / 0.04 + 0.6 / 0.01) / (1 / 0.04 + 1 / 0.01)
    assert math.log(result["fixed"]["or"]) == pytest.approx(expected)
    se = math.sqrt(1 / (1 / 0.04 + 1 / 0.01))
    assert math.log(result["fixed"]["ci_high"]) == pytest.approx(expected + 1.959963985 * se)


def test_dersimonian_laird_by_hand():
    studies = effects((0.1, 0.02), (0.9, 0.02), (0.5, 0.02))
    result = candidate_meta.pool(studies)
    w = 50.0
    mean = 0.5
    q = w * (0.16 + 0.16 + 0.0)
    tau2 = (q - 2) / (150 - (3 * w * w) / 150)
    assert result["q"] == pytest.approx(q)
    assert result["tau2"] == pytest.approx(tau2)
    assert result["i2_percent"] == pytest.approx(100 * (q - 2) / q)
    assert math.log(result["random"]["or"]) == pytest.approx(mean)
    se = math.sqrt((0.02 + tau2) / 3)
    assert math.log(result["random"]["ci_high"]) == pytest.approx(mean + 1.959963985 * se)


def test_homogeneous_studies_have_zero_tau_and_i2():
    result = candidate_meta.pool(effects((0.5, 0.02), (0.5, 0.03), (0.5, 0.01)))
    assert result["tau2"] == 0.0
    assert result["i2_percent"] == 0.0
    assert result["random"]["or"] == pytest.approx(result["fixed"]["or"])


def test_published_interval_converts_back_to_variance():
    study = candidate_meta.to_effect({"study": "x", "or": "2.0", "ci_low": "1.0", "ci_high": "4.0"}, 2)
    assert math.exp(study["log_or"]) == pytest.approx(2.0)
    assert study["variance"] == pytest.approx((math.log(4) / (2 * 1.959963985)) ** 2)


def test_counts_use_woolf_variance_and_correct_zero_cells():
    plain = candidate_meta.to_effect({"study": "x", "a": "20", "b": "10", "c": "10", "d": "20"}, 2)
    assert plain["log_or"] == pytest.approx(math.log(4))
    assert plain["variance"] == pytest.approx(1 / 20 + 1 / 10 + 1 / 10 + 1 / 20)
    corrected = candidate_meta.to_effect({"study": "x", "a": "0", "b": "10", "c": "4", "d": "20"}, 2)
    assert corrected["corrected"] is True


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ({"study": "x", "a": "0", "b": "5", "c": "0", "d": "9"}, "double-zero"),
        ({"study": "x", "or": "2", "ci_low": "3", "ci_high": "4"}, "ci_low"),
        ({"study": "x"}, "give a, b, c, d"),
    ],
)
def test_bad_studies_are_refused(row, message):
    with pytest.raises(_genetics.InputError, match=message):
        candidate_meta.to_effect(row, 2)


def test_hwe_sensitivity_and_leave_one_out():
    studies = candidate_meta.read_studies(SKILL_ROOT / "assets" / "example_studies.csv")
    result = candidate_meta.analyze(studies)
    assert result["hwe_sensitivity"]["excluded"] == ["Study D"]
    assert result["hwe_sensitivity"]["k"] == 4
    assert len(result["leave_one_out"]) == 5
    assert sum(row["weight_random_percent"] for row in result["studies"]) == pytest.approx(100, abs=0.5)


def test_high_heterogeneity_warns():
    result = candidate_meta.analyze(effects((0.0, 0.01), (1.5, 0.01), (0.7, 0.01), (2.0, 0.01), (0.1, 0.01)))
    assert any("studies disagree" in warning for warning in result["warnings"])


def test_one_study_is_not_a_meta_analysis():
    with pytest.raises(_genetics.InputError, match="at least two"):
        candidate_meta.analyze(effects((0.1, 0.02)))


# ------------------------------------------------------------------ report


def test_report_text_needs_no_plotting_library(tmp_path):
    variant = counts_analysis((12, 34, 32), (22, 40, 17))
    result = {"variants": [variant], "skipped": [{"variant": "x", "reason": "3 allele(s)"}]}
    path = tmp_path / "a.json"
    path.write_text(json.dumps(result), encoding="utf-8")
    out = tmp_path / "report"
    assert candidate_report.main(["--result", str(path), "--out-dir", str(out), "--no-figure"]) == 0
    text = (out / "report.md").read_text(encoding="utf-8")
    assert text.index("x skipped") < text.index("## Methods")
    assert "pre-specified additive model" in text
    assert "unadjusted" in text


def test_meta_report_and_forest_plot(tmp_path):
    pytest.importorskip("matplotlib")
    result = candidate_meta.analyze(candidate_meta.read_studies(SKILL_ROOT / "assets" / "example_studies.csv"))
    path = tmp_path / "m.json"
    path.write_text(json.dumps(result), encoding="utf-8")
    out = tmp_path / "report"
    assert candidate_report.main(["--meta", str(path), "--out-dir", str(out)]) == 0
    assert (out / "forest.png").stat().st_size > 5000
    assert "Leave-one-out" in (out / "report.md").read_text(encoding="utf-8")


def test_report_needs_an_input(tmp_path):
    assert candidate_report.main(["--out-dir", str(tmp_path)]) == 2
