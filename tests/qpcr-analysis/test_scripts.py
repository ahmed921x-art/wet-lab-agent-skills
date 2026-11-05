#!/usr/bin/env python3
"""Tests for qpcr-analysis with hand-derived and scipy-derived expected values."""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

stats = pytest.importorskip("scipy.stats")

import skill_contract  # noqa: E402

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "qpcr-analysis"
SCRIPTS = SKILL_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import _qpcr  # noqa: E402
import qpcr_efficiency  # noqa: E402
import qpcr_import  # noqa: E402
import qpcr_references  # noqa: E402
import qpcr_relative  # noqa: E402
import qpcr_report  # noqa: E402

CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)


def plate(cts: dict[str, tuple[str, dict[str, list]]]) -> list[dict]:
    """Build raw rows from {sample: (group, {target: [ct, ...]})}."""

    rows = []
    for sample, (group, targets) in cts.items():
        for target, wells in targets.items():
            for ct in wells:
                rows.append(
                    {
                        "sample": sample,
                        "group": group,
                        "target": target,
                        "ct": "Undetermined" if ct is None else str(ct),
                        "_line": len(rows) + 2,
                    }
                )
    return rows


def simple_plate() -> list[dict]:
    """dCt control 7.0, 7.2, 6.8 and treated 5.0, 5.2, 4.8, so ddCt is -2."""

    cts = {}
    for name, group, delta in (
        ("C1", "Ctrl", 7.0),
        ("C2", "Ctrl", 7.2),
        ("C3", "Ctrl", 6.8),
        ("T1", "Drug", 5.0),
        ("T2", "Drug", 5.2),
        ("T3", "Drug", 4.8),
    ):
        cts[name] = (group, {"REF": [18.0, 18.0], "GOI": [18.0 + delta] * 2})
    return plate(cts)


def contrast(result: dict) -> dict:
    assert len(result["contrasts"]) == 1
    return result["contrasts"][0]


def codes(result: dict) -> set[str]:
    return {warning["code"] for warning in result["warnings"]}


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("23.45", 23.45),
        ("23,45", 23.45),
        ("Undetermined", None),
        ("N/A", None),
        ("", None),
        ("0", None),
        ("-1", None),
    ],
)
def test_parse_ct(raw, expected):
    assert _qpcr.parse_ct(raw, line=2) == expected


def test_parse_ct_rejects_text_it_cannot_interpret():
    with pytest.raises(_qpcr.InputError, match="line 7"):
        _qpcr.parse_ct("low", line=7)


def test_benjamini_hochberg_by_hand():
    assert _qpcr.benjamini_hochberg([0.01, 0.04, 0.03, None]) == pytest.approx(
        [0.03, 0.04, 0.04, None]
    )


def test_livak_fold_change_and_welch_match_independent_values():
    row = contrast(qpcr_relative.analyze(simple_plate(), ["REF"], "Ctrl"))
    assert row["delta_delta_ct"] == pytest.approx(-2.0)
    assert row["fold_change"] == pytest.approx(4.0)

    drug, ctrl = [-5.0, -5.2, -4.8], [-7.0, -7.2, -6.8]
    oracle = stats.ttest_ind(drug, ctrl, equal_var=False)
    assert row["p_value"] == pytest.approx(oracle.pvalue)
    assert row["df"] == pytest.approx(oracle.df, abs=0.01)

    standard_error = math.sqrt(0.04 / 3 + 0.04 / 3)
    critical = stats.t.ppf(0.975, 4)
    assert row["fold_ci_low"] == pytest.approx(2 ** (2 - critical * standard_error), abs=1e-3)
    assert row["fold_ci_high"] == pytest.approx(2 ** (2 + critical * standard_error), abs=1e-3)


def test_technical_wells_do_not_inflate_n():
    rows = simple_plate()
    tripled = rows + [dict(row, _line=row["_line"] + 1000) for row in rows]
    single = contrast(qpcr_relative.analyze(rows, ["REF"], "Ctrl"))
    many = contrast(qpcr_relative.analyze(tripled, ["REF"], "Ctrl"))
    assert many["n_group"] == 3
    assert many["p_value"] == pytest.approx(single["p_value"])


def test_efficiency_correction_reproduces_the_pfaffl_ratio():
    e_target, e_reference = 1.9, 2.1
    cts = {
        "C1": ("Ctrl", {"REF": [18.0, 18.0], "GOI": [26.0, 26.0]}),
        "C2": ("Ctrl", {"REF": [18.4, 18.4], "GOI": [26.6, 26.6]}),
        "T1": ("Drug", {"REF": [18.1, 18.1], "GOI": [24.0, 24.0]}),
        "T2": ("Drug", {"REF": [18.3, 18.3], "GOI": [24.2, 24.2]}),
    }
    result = qpcr_relative.analyze(
        plate(cts), ["REF"], "Ctrl", {"GOI": e_target, "REF": e_reference}
    )
    expected = e_target ** (26.3 - 24.1) / e_reference ** (18.2 - 18.2)
    assert result["method"] == "pfaffl"
    assert contrast(result)["fold_change"] == pytest.approx(expected, abs=1e-3)
    assert "efficiency_unverified" not in codes(result)


def test_two_references_use_the_geometric_mean():
    cts = {
        "C1": ("Ctrl", {"R1": [18, 18], "R2": [20, 20], "GOI": [26, 26]}),
        "C2": ("Ctrl", {"R1": [18, 18], "R2": [20, 20], "GOI": [26, 26]}),
        "T1": ("Drug", {"R1": [18, 18], "R2": [21, 21], "GOI": [25, 25]}),
        "T2": ("Drug", {"R1": [18, 18], "R2": [21, 21], "GOI": [25, 25]}),
    }
    result = qpcr_relative.analyze(plate(cts), ["R1", "R2"], "Ctrl")
    assert contrast(result)["log2_fold_change"] == pytest.approx(1.5)


def test_genorm_m_is_the_sd_of_the_log_ratio_for_two_references():
    terms = {"a": {"R1": 18.0, "R2": 20.0}, "b": {"R1": 18.0, "R2": 21.0}, "c": {"R1": 18.0, "R2": 22.0}}
    values = qpcr_relative.genorm_m(terms, ["R1", "R2"])
    assert values["R1"] == pytest.approx(1.0)
    assert values["R2"] == pytest.approx(1.0)


def test_undetermined_is_left_out_not_imputed():
    rows = simple_plate()
    for row in rows:
        if row["sample"] == "T3" and row["target"] == "GOI":
            row["ct"] = "Undetermined"
    result = qpcr_relative.analyze(rows, ["REF"], "Ctrl")
    row = contrast(result)
    assert row["n_group"] == 2
    assert "target_not_quantified" in codes(result)
    assert {"target": "GOI", "group": "Drug", "n_samples": 3, "n_usable": 2} in result["detection"]


def test_lost_reference_excludes_the_sample():
    rows = simple_plate()
    for row in rows:
        if row["sample"] == "T1" and row["target"] == "REF":
            row["ct"] = ""
    result = qpcr_relative.analyze(rows, ["REF"], "Ctrl")
    assert result["n_samples_used"] == 5
    assert "sample_excluded" in codes(result)


def test_one_good_well_is_too_few_by_default():
    cts = {
        "C1": ("Ctrl", {"REF": [18, 18], "GOI": [25, None]}),
        "C2": ("Ctrl", {"REF": [18, 18], "GOI": [25, 25]}),
        "T1": ("Drug", {"REF": [18, 18], "GOI": [24, 24]}),
        "T2": ("Drug", {"REF": [18, 18], "GOI": [24, 24]}),
    }
    result = qpcr_relative.analyze(plate(cts), ["REF"], "Ctrl")
    status = {(w["sample"], w["target"]): w["status"] for w in result["wells"]}
    assert status[("C1", "GOI")] == "too_few_wells"


def test_noisy_technical_wells_are_flagged():
    cts = {
        "C1": ("Ctrl", {"REF": [18, 18], "GOI": [25, 26.5]}),
        "C2": ("Ctrl", {"REF": [18, 18], "GOI": [25, 25]}),
        "T1": ("Drug", {"REF": [18, 18], "GOI": [24, 24]}),
    }
    result = qpcr_relative.analyze(plate(cts), ["REF"], "Ctrl")
    assert any("well SD" in w["message"] for w in result["warnings"])
    assert "low_n" in codes(result)


def test_reference_that_moves_with_treatment_is_flagged():
    cts = {
        "C1": ("Ctrl", {"REF": [18, 18], "GOI": [25, 25]}),
        "C2": ("Ctrl", {"REF": [18, 18], "GOI": [25, 25]}),
        "T1": ("Drug", {"REF": [19.5, 19.5], "GOI": [25, 25]}),
        "T2": ("Drug", {"REF": [19.5, 19.5], "GOI": [25, 25]}),
    }
    result = qpcr_relative.analyze(plate(cts), ["REF"], "Ctrl")
    assert "reference_shift" in codes(result)


def test_ntc_amplification_near_samples_is_contamination():
    rows = simple_plate()
    rows.append({"sample": "NTC", "group": "", "target": "GOI", "ct": "27", "task": "NTC", "_line": 99})
    result = qpcr_relative.analyze(rows, ["REF"], "Ctrl")
    assert "control_contamination" in codes(result)
    rows[-1]["ct"] = "Undetermined"
    assert "control_contamination" not in codes(qpcr_relative.analyze(rows, ["REF"], "Ctrl"))


def test_paired_design_matches_ttest_rel():
    rows = simple_plate()
    pair_of = {"C1": "a", "T1": "a", "C2": "b", "T2": "b", "C3": "c", "T3": "c"}
    for row in rows:
        row["pair"] = pair_of[row["sample"]]
    row = contrast(qpcr_relative.analyze(rows, ["REF"], "Ctrl"))
    assert row["design"] == "paired"
    assert row["n_group"] == 3
    assert row["log2_fold_change"] == pytest.approx(2.0)
    assert row["t"] is None  # identical differences have no variance


def test_paired_p_value_against_scipy():
    rows = []
    for name, group, delta, pair in (
        ("C1", "Ctrl", 7.0, "a"),
        ("C2", "Ctrl", 7.5, "b"),
        ("C3", "Ctrl", 6.5, "c"),
        ("T1", "Drug", 5.1, "a"),
        ("T2", "Drug", 5.2, "b"),
        ("T3", "Drug", 4.9, "c"),
    ):
        for target, ct in (("REF", 18.0), ("GOI", 18.0 + delta)):
            for _ in range(2):
                rows.append(
                    {"sample": name, "group": group, "target": target, "ct": str(ct), "pair": pair, "_line": 2}
                )
    row = contrast(qpcr_relative.analyze(rows, ["REF"], "Ctrl"))
    oracle = stats.ttest_rel([-5.1, -5.2, -4.9], [-7.0, -7.5, -6.5])
    assert row["p_value"] == pytest.approx(oracle.pvalue)


def test_multiple_contrasts_receive_bh_adjustment():
    rows = simple_plate()
    for row in list(rows):
        if row["sample"].startswith("T") and row["target"] == "GOI":
            rows.append({**row, "sample": row["sample"] + "x", "group": "Other"})
        if row["sample"].startswith("T") and row["target"] == "REF":
            rows.append({**row, "sample": row["sample"] + "x", "group": "Other"})
    result = qpcr_relative.analyze(rows, ["REF"], "Ctrl")
    assert len(result["contrasts"]) == 2
    assert all(c["p_adjusted_bh"] >= c["p_value"] for c in result["contrasts"])


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"control": "Missing"}, "control group"),
        ({"references": ["NOPE"]}, "reference target"),
        ({"efficiency": {"GOI": 95.0}}, "amplification factor"),
    ],
)
def test_bad_arguments_are_refused(change, message):
    arguments = {"references": ["REF"], "control": "Ctrl", **change}
    with pytest.raises(_qpcr.InputError, match=message):
        qpcr_relative.analyze(simple_plate(), **arguments)


def test_sample_in_two_groups_is_refused():
    rows = simple_plate()
    rows[0]["group"] = "Drug"
    with pytest.raises(_qpcr.InputError, match="appears in groups"):
        qpcr_relative.analyze(rows, ["REF"], "Ctrl")


def perfect_curve(slope: float, levels: int = 5) -> list[dict]:
    return [
        {
            "target": "GOI",
            "quantity": str(10.0**-index),
            "ct": str(20.0 + slope * -index),
            "_line": index + 2,
        }
        for index in range(levels)
    ]


def test_perfect_doubling_is_100_percent():
    curve = qpcr_efficiency.analyze(perfect_curve(-math.log2(10)))["curves"][0]
    assert curve["efficiency_percent"] == pytest.approx(100.0, abs=0.01)
    assert curve["r_squared"] == pytest.approx(1.0)
    assert curve["flags"] == []


def test_low_efficiency_is_flagged_and_exported_as_a_factor():
    result = qpcr_efficiency.analyze(perfect_curve(-3.9))
    curve = result["curves"][0]
    assert curve["amplification_factor"] == pytest.approx(10 ** (1 / 3.9), abs=1e-3)
    assert any("outside" in flag for flag in curve["flags"])
    assert result["efficiency"]["GOI"] == curve["amplification_factor"]


def test_flat_or_rising_curve_is_refused():
    with pytest.raises(_qpcr.InputError, match="not negative"):
        qpcr_efficiency.amplification_factor(0.2)


def test_too_few_levels_gives_no_curve():
    result = qpcr_efficiency.analyze(perfect_curve(-3.32, levels=2))
    assert result["curves"] == []
    assert result["warnings"]


def test_cli_runs_the_documented_workflow(tmp_path):
    efficiency = tmp_path / "eff.json"
    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "qpcr_efficiency.py"),
            "--input",
            str(SKILL_ROOT / "assets" / "example_dilution_series.csv"),
            "--output-json",
            str(efficiency),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    output = tmp_path / "out.json"
    table = tmp_path / "out.csv"
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "qpcr_relative.py"),
            "--input",
            str(SKILL_ROOT / "assets" / "example_plate.csv"),
            "--reference",
            "GAPDH",
            "--reference",
            "ACTB",
            "--control",
            "Vehicle",
            "--efficiency-json",
            str(efficiency),
            "--output-json",
            str(output),
            "--output-csv",
            str(table),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "IL6 Treated vs Vehicle" in done.stdout
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["method"] == "pfaffl"
    assert payload["contrasts"][0]["p_value"] < 0.05
    assert table.read_text(encoding="utf-8").startswith("target,group,control")


def test_cli_reports_a_missing_file_cleanly(tmp_path):
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "qpcr_relative.py"),
            "--input",
            str(tmp_path / "absent.csv"),
            "--reference",
            "R",
            "--control",
            "C",
        ],
        capture_output=True,
        text=True,
    )
    assert done.returncode == 2
    assert "[FAIL]" in done.stderr


def test_regional_export_with_arabic_names_survives_a_legacy_console(tmp_path):
    lines = ["\ufeffSample;Group;Target;Ct;Task"]
    for sample, group, ct in (
        ("\u0639\u06461", "\u0636\u0627\u0628\u0637", "25,0"),
        ("\u0639\u06462", "\u0636\u0627\u0628\u0637", "25,2"),
        ("\u0639\u06463", "\u0645\u0639\u0627\u0644\u062c", "23,0"),
        ("\u0639\u06464", "\u0645\u0639\u0627\u0644\u062c", "23,2"),
    ):
        for _ in range(2):
            lines.append(f"{sample};{group};REF;18,0;Unknown")
            lines.append(f"{sample};{group};\u062c\u064a\u0646;{ct};Unknown")
    path = tmp_path / "plate.csv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "qpcr_relative.py"),
            "--input",
            str(path),
            "--reference",
            "REF",
            "--control",
            "\u0636\u0627\u0628\u0637",
        ],
        capture_output=True,
        env={**os.environ, "PYTHONIOENCODING": "cp1252", "PYTHONUTF8": "0"},
    )
    assert done.returncode == 0, done.stderr
    assert b"fold 4.0" in done.stdout


def test_utf16_export_is_read(tmp_path):
    path = tmp_path / "plate.csv"
    path.write_text("sample,group,target,ct\ns1,A,REF,18\n", encoding="utf-16")
    assert _qpcr.read_table(path, ("sample", "group", "target", "ct"))[0]["ct"] == "18"


def test_other_encodings_get_a_clear_message(tmp_path):
    path = tmp_path / "plate.csv"
    path.write_bytes("sample,group,target,ct\n\u0639,A,REF,18\n".encode("cp1256"))
    with pytest.raises(_qpcr.InputError, match="CSV UTF-8"):
        _qpcr.read_table(path, ("sample",))


# ------------------------------------------------------------------ import


def sheet() -> dict:
    return {
        name: {"sample": name, "group": group}
        for name, group in (
            ("Veh-1", "Vehicle"), ("Veh-2", "Vehicle"), ("Veh-3", "Vehicle"),
            ("Drug-1", "Treated"), ("Drug-2", "Treated"), ("Drug-3", "Treated"),
        )
    }


def example_grid() -> list[list[str]]:
    return qpcr_import.load_grid(SKILL_ROOT / "assets" / "example_quantstudio_export.csv", None)


def test_quantstudio_export_is_extracted_and_stops_at_the_next_section():
    result = qpcr_import.convert(example_grid(), sheet(), None)
    assert result["omitted_wells"] == 1
    assert result["standard_wells"] == 1
    assert result["samples"] == 6
    assert result["targets"] == ["ACTB", "GAPDH", "HK3", "IL6"]
    assert {w["group"] for w in result["wells"] if w["task"] == "unknown"} == {"Vehicle", "Treated"}
    ntc = [w for w in result["wells"] if w["task"] == "ntc"]
    assert ntc[0]["ct"] == "Undetermined"


def test_replicate_ct_column_is_used_not_the_mean_column():
    rows = [["Sample Name", "Target Name", "Ct Mean", "Ct"], ["s1", "T", "20.5", "20.1"], ["", "", "", ""]]
    result = qpcr_import.convert(rows, None, r"^(s)\d$")
    assert result["wells"][0]["ct"] == "20.1"


def test_cfx_layout_and_semicolons_are_read(tmp_path):
    path = tmp_path / "cfx.csv"
    path.write_text(
        "Well;Fluor;Target;Content;Sample;Biological Set Name;Cq;Cq Mean\n"
        "A01;SYBR;GOI;Unkn-01;A-1;;24,5;24,5\n"
        "A02;SYBR;GOI;Neg Ctrl;NTC;;N/A;\n"
        "A03;SYBR;GOI;Std;S1;;18,0;\n",
        encoding="utf-8",
    )
    result = qpcr_import.convert(qpcr_import.load_grid(path, None), None, r"^(\w)-\d$")
    assert [w["task"] for w in result["wells"]] == ["unknown", "ntc"]
    assert result["wells"][0]["group"] == "A"
    assert result["standard_wells"] == 1


def test_xlsx_export_is_read(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    book = openpyxl.Workbook()
    book.active.title = "Notes"
    book.active.append(["nothing here"])
    results = book.create_sheet("Results")
    results.append(["Instrument", "x"])
    results.append(["Sample Name", "Target Name", "Cq", "Task"])
    results.append(["s1", "GOI", 22.5, "Unknown"])
    path = tmp_path / "run.xlsx"
    book.save(path)
    result = qpcr_import.convert(qpcr_import.load_grid(path, None), None, r"^(s)\d$")
    assert result["wells"][0]["ct"] == "22.5"


def test_legacy_xls_gets_a_clear_message(tmp_path):
    with pytest.raises(_qpcr.InputError, match=".xls"):
        qpcr_import.load_grid(tmp_path / "run.xls", None)


def test_samples_without_a_group_are_refused():
    with pytest.raises(_qpcr.InputError, match="no group for 6"):
        qpcr_import.convert(example_grid(), None, None)
    partial = sheet()
    del partial["Veh-1"]
    with pytest.raises(_qpcr.InputError, match="Veh-1"):
        qpcr_import.convert(example_grid(), partial, None)


def test_unknown_task_and_missing_header_are_refused():
    with pytest.raises(_qpcr.InputError, match="not recognized"):
        qpcr_import.normalize_task("Mystery")
    with pytest.raises(_qpcr.InputError, match="no header"):
        qpcr_import.convert([["a", "b"], ["1", "2"]], None, None)


def test_imported_table_runs_through_quantification(tmp_path):
    output = tmp_path / "wells.csv"
    assert (
        qpcr_import.main(
            [
                "--input", str(SKILL_ROOT / "assets" / "example_quantstudio_export.csv"),
                "--sample-sheet", str(SKILL_ROOT / "assets" / "example_sample_sheet.csv"),
                "--output", str(output),
            ]
        )
        == 0
    )
    rows = _qpcr.read_table(output, ("sample", "group", "target", "ct"))
    result = qpcr_relative.analyze(rows, ["GAPDH", "ACTB"], "Vehicle")
    assert {row["target"] for row in result["contrasts"]} == {"HK3", "IL6"}
    assert result["samples"]


# -------------------------------------------------------------- references


def reference_plate(shift_hk3: float = 0.0) -> list[dict]:
    cts = {}
    wobble = [0.0, 0.3, -0.2, 0.1, -0.3, 0.2, 0.25, -0.15]
    for index, delta in enumerate(wobble):
        group = "A" if index < 4 else "B"
        cts[f"s{index}"] = (
            group,
            {
                "G1": [18 + delta] * 2,
                "G2": [20 + delta + 0.05 * (-1) ** index] * 2,
                "G3": [22 + delta + (index % 3) * 0.6] * 2,
                "HK3": [24 + delta + (shift_hk3 if group == "B" else 0.0) + 0.02 * index] * 2,
            },
        )
    return plate(cts)


def test_stepwise_ranking_puts_the_noisy_gene_last():
    result = qpcr_references.analyze(reference_plate(), ["G1", "G2", "G3"])
    assert result["ranking_most_stable_first"][-1] == "G3"
    assert result["steps"][0]["removed"] == "G3"


def test_stability_matches_the_single_step_genorm_helper():
    values, _ = qpcr_references.relative_quantities(reference_plate(), ["G1", "G2"], {})
    assert qpcr_references.stability(values, ["G1", "G2"])["G1"] == pytest.approx(
        qpcr_relative.genorm_m(values, ["G1", "G2"])["G1"]
    )


def test_pairwise_variation_by_hand():
    values = {
        "a": {"G1": 0.0, "G2": 0.0, "G3": 3.0},
        "b": {"G1": 0.0, "G2": 0.0, "G3": 0.0},
        "c": {"G1": 0.0, "G2": 0.0, "G3": -3.0},
    }
    result = qpcr_references.rank(values, ["G1", "G2", "G3"])
    assert result["pairwise_variation"][0]["v"] == pytest.approx(1.0)
    assert result["recommended_count"] == 3


def test_gene_moved_by_treatment_is_not_recommended_even_if_genorm_likes_it():
    result = qpcr_references.analyze(reference_plate(shift_hk3=1.2), ["G1", "G2", "HK3"])
    assert "HK3" in result["shift_between_groups"]
    assert "HK3" not in result["recommended_references"]
    assert any("treatment appears to move" in warning for warning in result["warnings"])


def test_stable_genes_are_not_implicated_by_a_strong_mover():
    result = qpcr_references.analyze(reference_plate(shift_hk3=3.0), ["G1", "G2", "HK3"])
    assert list(result["shift_between_groups"]) == ["HK3"]


def test_two_candidates_cannot_be_ranked():
    result = qpcr_references.analyze(reference_plate(), ["G1", "G2"])
    assert any("two candidates" in warning for warning in result["warnings"])
    assert result["pairwise_variation"] == []


def test_too_few_samples_are_refused():
    with pytest.raises(_qpcr.InputError, match="at least 3"):
        qpcr_references.analyze(simple_plate()[:4], ["REF", "GOI"])


# ------------------------------------------------------------------ report


def test_report_puts_warnings_first_and_states_the_method():
    result = qpcr_relative.analyze(simple_plate(), ["REF"], "Ctrl")
    text = qpcr_report.build_markdown(result, None, None)
    assert text.index("efficiency_unverified") < text.index("## Methods")
    assert "assuming 100% efficiency" in text
    assert "4.00-fold" in text
    assert "n = 3 vs 3" in text


def test_report_lists_measured_efficiencies():
    result = qpcr_relative.analyze(simple_plate(), ["REF"], "Ctrl", {"GOI": 1.95, "REF": 1.9})
    text = qpcr_report.build_markdown(result, {"GOI": 1.95, "REF": 1.9}, None)
    assert "Pfaffl" in text
    assert "GOI 95.0%" in text


def test_figure_and_report_are_written(tmp_path):
    pytest.importorskip("matplotlib")
    result_path = tmp_path / "result.json"
    _qpcr.write_json(result_path, qpcr_relative.analyze(simple_plate(), ["REF"], "Ctrl"))
    out = tmp_path / "report"
    assert qpcr_report.main(["--result", str(result_path), "--out-dir", str(out)]) == 0
    assert (out / "fold_change.png").stat().st_size > 5000
    assert "fold_change.png" in (out / "report.md").read_text(encoding="utf-8")


def test_old_result_without_samples_gives_a_clear_error(tmp_path):
    pytest.importorskip("matplotlib")
    path = tmp_path / "old.json"
    path.write_text(json.dumps({"control": "C"}), encoding="utf-8")
    assert qpcr_report.main(["--result", str(path), "--out-dir", str(tmp_path / "x")]) == 2


def test_importer_stops_at_the_end_of_the_results_table():
    rows = [
        ["Sample Name", "Target Name", "Ct"],
        ["s1", "GOI", "20.0"],
        [],
        ["[Replicate Group Result]"],
        ["Sample Name", "Target Name", "Ct"],
        ["s1", "GOI", "99.0"],
    ]
    result = qpcr_import.convert(rows, None, r"^(s)\d$")
    assert [well["ct"] for well in result["wells"]] == ["20.0"]
