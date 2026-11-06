#!/usr/bin/env python3
"""Tests for paper-audit against independent values.

Distribution values come from scipy; GRIM and percentage checks are compared with
brute-force enumeration in exact fractions; the remaining cases are hand
calculations written out in the comments.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

stats = pytest.importorskip("scipy.stats")
special = pytest.importorskip("scipy.special")

import skill_contract  # noqa: E402

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "paper-audit"
SCRIPTS = SKILL_ROOT / "scripts"
EXAMPLE = SKILL_ROOT / "assets" / "example_claims.csv"
sys.path.insert(0, str(SCRIPTS))

import _audit  # noqa: E402
import audit_claims  # noqa: E402
import audit_report  # noqa: E402

CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)


def claim(check, **fields):
    row = {"id": "x", "page": "1", "quote": "q", "check": check}
    row.update({key: str(value) for key, value in fields.items()})
    return row


def status_of(row):
    return audit_claims.audit([row])["claims"][0]


# ---------------------------------------------------------------------------
# Special functions and p-values against scipy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("a", [0.5, 1.0, 2.5, 14.0, 150.0])
@pytest.mark.parametrize("b", [0.5, 1.0, 3.0, 40.0])
@pytest.mark.parametrize("x", [1e-6, 0.01, 0.3, 0.5, 0.77, 0.999])
def test_incomplete_beta_matches_scipy(a, b, x):
    assert _audit.betainc(a, b, x) == pytest.approx(special.betainc(a, b, x), rel=1e-10, abs=1e-14)


@pytest.mark.parametrize("a", [0.5, 1.0, 1.5, 4.0, 30.0])
@pytest.mark.parametrize("x", [1e-4, 0.4, 1.0, 3.7, 25.0, 80.0])
def test_upper_incomplete_gamma_matches_scipy(a, x):
    assert _audit.gammaincc(a, x) == pytest.approx(special.gammaincc(a, x), rel=1e-9, abs=1e-14)


@pytest.mark.parametrize("t,df", [(0.0, 5), (0.4, 3), (2.2, 28), (-2.2, 28), (1.96, 1000), (3.1, 7.6), (12.0, 4)])
def test_t_p_values_match_scipy(t, df):
    assert _audit.p_from_t(t, df) == pytest.approx(2 * stats.t.sf(abs(t), df), rel=1e-9)
    assert _audit.p_from_t(t, df, tails=1) == pytest.approx(stats.t.sf(abs(t), df), rel=1e-9)


@pytest.mark.parametrize("f,df1,df2", [(0.5, 1, 10), (3.41, 2, 57), (4.0, 3, 120), (25.0, 4, 9), (1.0, 6, 6.5)])
def test_f_p_values_match_scipy(f, df1, df2):
    assert _audit.p_from_f(f, df1, df2) == pytest.approx(stats.f.sf(f, df1, df2), rel=1e-9)


@pytest.mark.parametrize("x,df", [(0.2, 1), (3.84, 1), (6.8, 3), (40.0, 10), (150.0, 2), (0.5, 7)])
def test_chi2_p_values_match_scipy(x, df):
    assert _audit.p_from_chi2(x, df) == pytest.approx(stats.chi2.sf(x, df), rel=1e-9)


@pytest.mark.parametrize("z", [0.0, 0.5, 1.75, 1.959964, 4.5, -2.3])
def test_z_p_values_match_scipy(z):
    assert _audit.p_from_z(z) == pytest.approx(2 * stats.norm.sf(abs(z)), rel=1e-12)


@pytest.mark.parametrize("r,n", [(0.31, 45), (-0.5, 12), (0.05, 400), (0.9, 6)])
def test_r_p_values_match_scipy_pearson(r, n):
    # scipy's exact null distribution of Pearson r is the reference.
    expected = 2 * stats.beta(n / 2 - 1, n / 2 - 1, loc=-1, scale=2).sf(abs(r))
    assert _audit.p_from_r(r, n - 2) == pytest.approx(expected, rel=1e-8)


@pytest.mark.parametrize("q,df", [(0.975, 5), (0.975, 28), (0.95, 2.5), (0.995, 100)])
def test_t_quantile_matches_scipy(q, df):
    assert _audit.t_quantile(q, df) == pytest.approx(stats.t.ppf(q, df), rel=1e-9)


# ---------------------------------------------------------------------------
# p-value check with rounding
# ---------------------------------------------------------------------------


def test_reported_t_test_consistent_within_rounding():
    # scipy: t(28)=2.195 -> .03662, 2.205 -> .03584; reported .036 covers [.0355, .0365].
    entry = status_of(claim("p_value", test="t", statistic="2.20", df=28, p="= .036"))
    assert entry["status"] == "consistent"
    low, high = entry["recomputed"]["p_range_given_rounding"]
    assert low == pytest.approx(stats.t.sf(2.205, 28) * 2, rel=1e-3)
    assert high == pytest.approx(stats.t.sf(2.195, 28) * 2, rel=1e-3)


def test_rounding_of_statistic_rescues_a_borderline_p():
    # Exact t = 2.05 at df 20 gives p = .0537, but the printed 2.05 may be 2.0549, giving p = .0532.
    # Reported p = .053 ([.0525, .0535]) only fits because of the statistic's rounding interval.
    assert 2 * stats.t.sf(2.05, 20) > 0.0535
    assert 2 * stats.t.sf(2.055, 20) < 0.0535
    assert status_of(claim("p_value", test="t", statistic="2.05", df=20, p=".053"))["status"] == "consistent"


def test_wrong_p_is_flagged_and_decision_error_marked():
    # chi2(3) = 6.80 gives p = .0786 (scipy), so "p < .05" is a change of decision.
    entry = status_of(claim("p_value", test="chi2", statistic="6.80", df=3, p="< .05"))
    assert entry["status"] == "inconsistent"
    assert entry["decision_error"] is True


def test_one_tailed_hint_without_changing_status():
    # z = 1.75: two-tailed p = .0801, one-tailed p = .0401 (scipy); reported .04.
    entry = status_of(claim("p_value", test="z", statistic="1.75", p=".04"))
    assert entry["status"] == "inconsistent"
    assert any("one-tailed" in note for note in entry["notes"])
    assert status_of(claim("p_value", test="z", statistic="1.75", p=".04", tails=1))["status"] == "consistent"


def test_non_decision_mismatch_is_not_a_decision_error():
    # t(28) = 2.20 gives p = .036, reported .030: wrong, but both are below .05.
    entry = status_of(claim("p_value", test="t", statistic="2.20", df=28, p=".030"))
    assert entry["status"] == "inconsistent"
    assert entry["decision_error"] is False


def test_inequalities_and_printed_zero():
    assert status_of(claim("p_value", test="t", statistic="5.0", df=30, p="< .001"))["status"] == "consistent"
    assert status_of(claim("p_value", test="t", statistic="5.0", df=30, p="p = .000"))["status"] == "consistent"
    assert status_of(claim("p_value", test="t", statistic="1.0", df=30, p="> .05"))["status"] == "consistent"
    assert status_of(claim("p_value", test="t", statistic="1.0", df=30, p="ns"))["status"] == "consistent"
    assert status_of(claim("p_value", test="t", statistic="1.0", df=30, p="< .05"))["status"] == "inconsistent"


def test_f_and_r_checks():
    # scipy: F(2,57)=3.41 -> p = .0399; r = .31 with n = 45 -> p = .0382.
    assert status_of(claim("p_value", test="F", statistic="3.41", df=2, df2=57, p=".04"))["status"] == "consistent"
    assert status_of(claim("p_value", test="r", statistic=".31", n=45, p=".038"))["status"] == "consistent"
    assert status_of(claim("p_value", test="r", statistic=".31", df=43, p=".038"))["status"] == "consistent"


def test_decimal_df_from_welch_is_handled():
    # Welch df printed as 17.3 covers [17.25, 17.35]; scipy gives p = .0486 at t = 2.13, df 17.3.
    p = 2 * stats.t.sf(2.13, 17.3)
    entry = status_of(claim("p_value", test="t", statistic="2.13", df="17.3", p=f"{p:.3f}"))
    assert entry["status"] == "consistent"


# ---------------------------------------------------------------------------
# GRIM against brute-force enumeration
# ---------------------------------------------------------------------------


def possible_by_enumeration(text: str, cells: int) -> bool:
    """Brute force: does some integer total give a mean that rounds to text (either way at a half)?"""

    decimals = len(text.split(".")[1])
    target = Fraction(text)
    half = Fraction(1, 2 * 10**decimals)
    return any(abs(Fraction(total, cells) - target) <= half for total in range(0, 10 * cells + 1))


@pytest.mark.parametrize("cells", [3, 7, 12, 19, 28, 33, 47, 60, 99])
def test_grim_matches_enumeration(cells):
    for hundredths in range(100, 700, 7):
        text = f"{hundredths / 100:.2f}"
        mean = _audit.parse_reported(text, "mean")
        assert _audit.grim_possible(mean, cells)[0] == possible_by_enumeration(text, cells), (text, cells)


def test_grim_hand_examples():
    # 5.19 x 28 = 145.32; 145/28 = 5.1786 and 146/28 = 5.2143, neither rounds to 5.19.
    assert status_of(claim("grim", mean="5.19", n=28))["status"] == "inconsistent"
    # 87/25 = 3.48 exactly.
    assert status_of(claim("grim", mean="3.48", n=25))["status"] == "consistent"
    # Three items, 20 people: 60 cells; 208/60 = 3.4667 -> 3.47, while 3.49 needs 209.4 (209/60 = 3.483, 210/60 = 3.5).
    assert status_of(claim("grim", mean="3.47", n=20, items=3))["status"] == "consistent"
    assert status_of(claim("grim", mean="3.49", n=20, items=3))["status"] == "inconsistent"


def test_grim_exact_half_accepts_either_rounding():
    # 1/8 = 0.125 prints as 0.12 or 0.13 depending on convention; both must pass.
    assert status_of(claim("grim", mean="0.12", n=8))["status"] == "consistent"
    assert status_of(claim("grim", mean="0.13", n=8))["status"] == "consistent"


def test_grim_not_checkable_when_n_is_large():
    # With 2 decimals every mean is reachable once n >= 100.
    assert status_of(claim("grim", mean="2.47", n=100))["status"] == "not checkable"
    # Just below 100 a few means are still impossible: 148/99 = 1.4949 and 149/99 = 1.5051, so 1.50 cannot occur.
    assert status_of(claim("grim", mean="1.50", n=99))["status"] == "inconsistent"
    assert status_of(claim("grim", mean="2.47", n=99))["status"] == "consistent"  # 245/99 = 2.4747


# ---------------------------------------------------------------------------
# Percentages and sums
# ---------------------------------------------------------------------------


def test_percent_with_count():
    assert status_of(claim("percent", percent="41.7", count=25, n=60))["status"] == "consistent"  # 41.667
    assert status_of(claim("percent", percent="35.0", count=21, n=62))["status"] == "inconsistent"  # 33.871
    assert status_of(claim("percent", percent="50", count=8, n=7))["status"] == "inconsistent"


def test_percent_without_count_matches_enumeration():
    for n in (7, 13, 30, 61):
        for tenths in range(0, 1001, 3):
            text = f"{tenths / 10:.1f}"
            expected = any(abs(Fraction(100 * k, n) - Fraction(text)) <= Fraction(1, 20) for k in range(n + 1))
            status = status_of(claim("percent", percent=text, n=n))["status"]
            assert status == ("consistent" if expected else "inconsistent"), (text, n)


def test_percent_not_checkable_when_fine_grained():
    assert status_of(claim("percent", percent="37", n=500))["status"] == "not checkable"


def test_sums():
    assert status_of(claim("sum", parts="30;31;28", total=90))["status"] == "inconsistent"
    assert status_of(claim("sum", parts="30; 31; 29", total=90))["status"] == "consistent"
    # 33.3 + 33.3 + 33.3 = 99.9; three rounded parts allow 0.15 slack against an exact 100.
    assert status_of(claim("sum", parts="33.3;33.3;33.3", total=100))["status"] == "consistent"
    assert status_of(claim("sum", parts="48.3;50.7", total=100))["status"] == "inconsistent"


# ---------------------------------------------------------------------------
# Confidence intervals
# ---------------------------------------------------------------------------


def test_ratio_ci_symmetric_on_log_scale_and_p_agrees():
    # sqrt(1.21 x 2.95) = 1.889; SE = ln(2.95/1.21)/3.92 = 0.2274; z = ln(1.889)/0.2274 = 2.80; p = .0051.
    entry = status_of(claim("ratio_ci", estimate="1.89", ci_low="1.21", ci_high="2.95", p=".005"))
    assert entry["status"] == "consistent"
    se = (math.log(2.95) - math.log(1.21)) / (2 * stats.norm.ppf(0.975))
    expected = 2 * stats.norm.sf(math.log(1.89) / se)
    assert entry["recomputed"]["p_from_interval"] == pytest.approx(expected, rel=1e-3)


def test_ratio_ci_asymmetry_flagged_unless_method_says_so():
    row = claim("ratio_ci", estimate="2.40", ci_low="1.10", ci_high="3.20")
    assert status_of(row)["status"] == "inconsistent"  # sqrt(1.1 x 3.2) = 1.876
    row["ci_method"] = "exact"
    assert status_of(row)["status"] == "consistent"


def test_ratio_ci_and_p_disagree_on_significance():
    entry = status_of(claim("ratio_ci", estimate="1.45", ci_low="0.92", ci_high="2.28", p=".03"))
    assert entry["status"] == "inconsistent"
    assert entry["decision_error"] is True


def test_interval_bound_at_null_is_not_judged():
    entry = status_of(claim("ratio_ci", estimate="1.50", ci_low="1.00", ci_high="2.25", p=".05"))
    assert entry["decision_error"] is False
    assert any("equals the null" in note for note in entry["notes"])


def test_estimate_outside_interval():
    entry = status_of(claim("diff_ci", estimate="5.0", ci_low="1.0", ci_high="4.0"))
    assert entry["status"] == "inconsistent"
    assert "outside" in entry["reason"]


def test_difference_ci_uses_t_when_df_given():
    # Midpoint 4.2; with df 20, SE = 3.1 / t(.975, 20) and p = 2 * t.sf(4.2 / SE, 20).
    entry = status_of(claim("diff_ci", estimate="4.2", ci_low="1.1", ci_high="7.3", df=20, p=".01"))
    se = 3.1 / stats.t.ppf(0.975, 20)
    assert entry["recomputed"]["p_from_interval"] == pytest.approx(2 * stats.t.sf(4.2 / se, 20), rel=1e-3)
    assert entry["status"] == "consistent"


def test_ninety_percent_interval_uses_alpha_ten_percent():
    # A 90% CI that excludes 0 corresponds to p < .10, so p = .07 is consistent with it.
    entry = status_of(claim("diff_ci", estimate="1.0", ci_low="0.1", ci_high="1.9", level=90, p=".07"))
    assert entry["decision_error"] is False


# ---------------------------------------------------------------------------
# Bad input and traceability
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "row,fragment",
    [
        (claim("p_value", test="t", statistic="2,20", df=28, p=".03"), "decimal comma"),
        (claim("p_value", test="t", statistic="abc", df=28, p=".03"), "not a plain number"),
        (claim("p_value", test="t", statistic="2.2", p=".03"), "missing df"),
        (claim("p_value", test="q", statistic="2.2", df=3, p=".03"), "not one of"),
        (claim("p_value", test="t", statistic="2.2", df=3, p="1.4"), "outside 0 to 1"),
        (claim("p_value", test="r", statistic="1.2", n=20, p=".03"), "between -1 and 1"),
        (claim("p_value", test="F", statistic="-1", df=2, df2=10, p=".03"), "negative"),
        (claim("grim", mean="3.4", n="0"), "at least 1"),
        (claim("grim", mean="3.4", n="12.5"), "whole number"),
        (claim("ratio_ci", estimate="-1", ci_low="0.5", ci_high="2"), "positive"),
        (claim("diff_ci", estimate="1", ci_low="3", ci_high="2"), "larger than"),
        (claim("sum", parts="12", total="12"), "at least two"),
        (claim("nonsense"), "is not one of"),
    ],
)
def test_bad_rows_become_not_checkable_with_a_reason(row, fragment):
    entry = status_of(row)
    assert entry["status"] == "not checkable"
    assert fragment in entry["reason"]


def test_unicode_minus_and_leading_dot_are_read():
    assert status_of(claim("p_value", test="t", statistic="−2.20", df=28, p=".036"))["status"] == "consistent"


def test_claim_without_quote_or_page_is_not_audited():
    row = claim("grim", mean="5.19", n=28)
    row["quote"] = ""
    entry = status_of(row)
    assert entry["status"] == "not checkable"
    assert "traced" in entry["reason"]


def test_duplicate_ids_are_rejected():
    with pytest.raises(_audit.InputError):
        audit_claims.audit([claim("grim", mean="3.48", n=25), claim("grim", mean="3.48", n=25)])


def test_readers(tmp_path):
    rows = [claim("grim", mean="3.48", n=25)]
    json_file = tmp_path / "claims.json"
    json_file.write_text(json.dumps({"claims": rows}), encoding="utf-8")
    assert audit_claims.read_claims(json_file)[0]["mean"] == "3.48"
    utf16 = tmp_path / "claims16.csv"
    utf16.write_text("id,page,quote,check,mean,n\nA,2,\"Mean was 3.48, n 25\",grim,3.48,25\n", encoding="utf-16")
    assert audit_claims.audit(audit_claims.read_claims(utf16))["counts"]["consistent"] == 1
    for content, fragment in [("", "empty"), ("id,page\n1,2\n", "missing column"),
                              ("id,page,quote,check\n1,2,a,b,extra\n", "more cells"), ("[1, 2]", "list of claim")]:
        bad = tmp_path / "bad.csv"
        bad.write_text(content, encoding="utf-8")
        with pytest.raises(_audit.InputError, match=fragment):
            audit_claims.read_claims(bad)
    with pytest.raises(_audit.InputError, match="no such file"):
        audit_claims.read_claims(tmp_path / "missing.csv")


# ---------------------------------------------------------------------------
# Example, CLI, and report
# ---------------------------------------------------------------------------

EXPECTED_EXAMPLE = {
    "C01": "consistent", "C02": "consistent", "C03": "inconsistent", "C04": "consistent", "C05": "inconsistent",
    "C06": "inconsistent", "C07": "consistent", "C08": "not checkable", "C09": "consistent", "C10": "inconsistent",
    "C11": "inconsistent", "C12": "consistent", "C13": "consistent", "C14": "inconsistent", "C15": "inconsistent",
    "C16": "consistent", "C17": "not checkable", "C18": "consistent",
}


def test_example_is_labelled_synthetic_and_gives_expected_ledger(tmp_path):
    text = EXAMPLE.read_text(encoding="utf-8")
    rows = audit_claims.read_claims(EXAMPLE)
    assert all("SYNTHETIC" in row["note"] for row in rows)
    ledger_path = tmp_path / "ledger.json"
    csv_path = tmp_path / "ledger.csv"
    done = subprocess.run(
        [sys.executable, str(SCRIPTS / "audit_claims.py"), "--input", str(EXAMPLE), "--output-json", str(ledger_path),
         "--output-csv", str(csv_path)],
        capture_output=True, text=True, encoding="utf-8", timeout=120,
    )
    assert done.returncode == 0, done.stderr
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    assert {entry["id"]: entry["status"] for entry in ledger["claims"]} == EXPECTED_EXAMPLE
    assert csv_path.read_text(encoding="utf-8").count("\n") == len(EXPECTED_EXAMPLE) + 1
    assert text  # the asset is non-empty

    report_path = tmp_path / "report.md"
    done = subprocess.run(
        [sys.executable, str(SCRIPTS / "audit_report.py"), "--ledger", str(ledger_path), "--out", str(report_path)],
        capture_output=True, text=True, encoding="utf-8", timeout=120,
    )
    assert done.returncode == 0, done.stderr
    report = report_path.read_text(encoding="utf-8")
    assert "**not** evidence of misconduct" in report
    assert "What this audit cannot detect" in report
    assert report.index("Inconsistencies to examine") < report.index("## Consistent")
    for entry in ledger["claims"]:
        if entry["status"] == "inconsistent":
            assert f"(page {entry['page']})" in report
            assert entry["quote"] in report


def test_cli_reports_input_errors(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("id,page\n", encoding="utf-8")
    done = subprocess.run([sys.executable, str(SCRIPTS / "audit_claims.py"), "--input", str(bad)],
                          capture_output=True, text=True, timeout=120)
    assert done.returncode == 2 and "[FAIL]" in done.stderr
    notledger = tmp_path / "x.json"
    notledger.write_text("{}", encoding="utf-8")
    done = subprocess.run([sys.executable, str(SCRIPTS / "audit_report.py"), "--ledger", str(notledger),
                           "--out", str(tmp_path / "r.md")], capture_output=True, text=True, timeout=120)
    assert done.returncode == 2


def test_markdown_escapes_table_pipes():
    ledger = audit_claims.audit([dict(claim("grim", mean="3.48", n=25), quote="a | b")])
    ledger["claims"][0]["reason"] = "x | y"
    text = audit_report.build_markdown(ledger)
    assert "x \\| y" in text and "x | y" not in text
