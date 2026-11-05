#!/usr/bin/env python3
"""Tests for western-blot-densitometry with hand-derived expected values."""

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

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "western-blot-densitometry"
SCRIPTS = SKILL_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import blot_linearity  # noqa: E402
import blot_quant  # noqa: E402
import blot_report  # noqa: E402

CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)


def band(blot, sample, group, target, signal, saturated=False, pair=""):
    return {
        "blot": blot,
        "sample": sample,
        "group": group,
        "target": target,
        "signal": float(signal),
        "saturated": saturated,
        "pair": pair,
    }


def membrane(name="B1", exposure=1.0, ratios=None):
    """Three control and three treated lanes with chosen target:loading ratios."""

    ratios = ratios or {"C1": 0.5, "C2": 0.6, "C3": 0.4, "T1": 1.0, "T2": 1.2, "T3": 0.8}
    bands = []
    for sample, ratio in ratios.items():
        group = "Ctrl" if sample.startswith("C") else "Drug"
        bands.append(band(name, sample, group, "ACTB", 1000 * exposure))
        bands.append(band(name, sample, group, "GOI", 1000 * exposure * ratio))
    return bands


def contrast(result):
    assert len(result["contrasts"]) == 1
    return result["contrasts"][0]


def test_within_blot_normalization_and_welch_match_independent_values():
    row = contrast(blot_quant.analyze(membrane(), "ACTB", "Ctrl"))
    control = [math.log2(x) for x in (0.5, 0.6, 0.4)]
    drug = [math.log2(x) for x in (1.0, 1.2, 0.8)]
    assert row["log2_fold_change"] == pytest.approx(sum(drug) / 3 - sum(control) / 3, abs=1e-3)
    oracle = stats.ttest_ind(drug, control, equal_var=False)
    assert row["p_value"] == pytest.approx(oracle.pvalue)


def test_exposure_differences_between_blots_cancel():
    first = contrast(blot_quant.analyze(membrane(), "ACTB", "Ctrl"))
    merged = membrane("B1") + membrane("B2", exposure=7.0)
    second = contrast(blot_quant.analyze(merged, "ACTB", "Ctrl"))
    assert second["log2_fold_change"] == pytest.approx(first["log2_fold_change"], abs=1e-3)


def test_repeat_blots_do_not_inflate_n():
    result = blot_quant.analyze(membrane("B1") + membrane("B2"), "ACTB", "Ctrl")
    assert contrast(result)["n_group"] == 3
    assert all(row["n_blots"] == 2 for row in result["samples"])


def test_saturated_band_is_excluded_by_default():
    bands = membrane()
    bands[1]["saturated"] = True
    result = blot_quant.analyze(bands, "ACTB", "Ctrl")
    assert {"blot": "B1", "sample": "C1", "target": "GOI", "reason": "saturated"} in result["excluded"]
    assert contrast(result)["n_control"] == 2
    kept = blot_quant.analyze(bands, "ACTB", "Ctrl", keep_saturated=True)
    assert contrast(kept)["n_control"] == 3
    assert any("kept by request" in warning for warning in kept["warnings"])


def test_band_without_signal_is_excluded_not_logged():
    bands = membrane()
    bands[1]["signal"] = -5.0
    result = blot_quant.analyze(bands, "ACTB", "Ctrl")
    assert result["excluded"][0]["reason"] == "no signal above background"


def test_lane_without_loading_control_is_excluded():
    bands = [row for row in membrane() if not (row["sample"] == "T1" and row["target"] == "ACTB")]
    result = blot_quant.analyze(bands, "ACTB", "Ctrl")
    assert result["excluded"][0]["reason"] == "no usable loading band"
    assert contrast(result)["n_group"] == 2


def test_blot_without_control_lanes_is_left_out():
    extra = [row for row in membrane("B2") if row["group"] == "Drug"]
    result = blot_quant.analyze(membrane("B1") + extra, "ACTB", "Ctrl")
    assert any("no control lane" in warning for warning in result["warnings"])
    assert all(row["n_blots"] == 1 for row in result["samples"])


def test_uneven_loading_is_flagged():
    bands = membrane()
    for row in bands:
        if row["target"] == "ACTB" and row["sample"] == "T1":
            row["signal"] = 3000.0
    assert any("loading control CV" in w for w in blot_quant.analyze(bands, "ACTB", "Ctrl")["warnings"])


def test_loading_control_that_moves_with_group_is_flagged():
    bands = membrane()
    for row in bands:
        if row["target"] == "ACTB" and row["group"] == "Drug":
            row["signal"] = 2000.0
    warnings = blot_quant.analyze(bands, "ACTB", "Ctrl")["warnings"]
    assert any("total-protein" in warning for warning in warnings)


def test_paired_design_matches_ttest_rel():
    bands = []
    for pair, control_ratio, drug_ratio in (("a", 0.5, 1.1), ("b", 0.7, 1.2), ("c", 0.4, 0.9)):
        for sample, group, ratio in ((f"C{pair}", "Ctrl", control_ratio), (f"T{pair}", "Drug", drug_ratio)):
            bands.append(band("B1", sample, group, "ACTB", 1000, pair=pair))
            bands.append(band("B1", sample, group, "GOI", 1000 * ratio, pair=pair))
    row = contrast(blot_quant.analyze(bands, "ACTB", "Ctrl"))
    assert row["design"] == "paired"
    oracle = stats.ttest_rel(
        [math.log2(x) for x in (1.1, 1.2, 0.9)], [math.log2(x) for x in (0.5, 0.7, 0.4)]
    )
    assert row["p_value"] == pytest.approx(oracle.pvalue)


def test_small_groups_are_called_descriptive():
    ratios = {"C1": 0.5, "C2": 0.6, "T1": 1.0, "T2": 1.2}
    assert any("descriptive" in w for w in blot_quant.analyze(membrane(ratios=ratios), "ACTB", "Ctrl")["warnings"])


@pytest.mark.parametrize(
    ("loading", "control", "message"),
    [("NOPE", "Ctrl", "loading control"), ("ACTB", "Nope", "control group")],
)
def test_bad_labels_are_refused(loading, control, message):
    with pytest.raises(blot_quant.InputError, match=message):
        blot_quant.analyze(membrane(), loading, control)


def test_sample_in_two_groups_is_refused():
    bands = membrane()
    bands[0]["group"] = "Drug"
    with pytest.raises(blot_quant.InputError, match="more than one group"):
        blot_quant.analyze(bands, "ACTB", "Ctrl")


def test_reader_subtracts_background_and_reads_flags(tmp_path):
    path = tmp_path / "bands.csv"
    path.write_text(
        "Blot,Sample,Group,Target,Signal,Background,Saturated\nB,s,g,T,1000,100,yes\nB,s2,g,T,\"900,5\",,no\n",
        encoding="utf-8",
    )
    rows = blot_quant.read_bands(path)
    assert rows[0]["signal"] == 900.0
    assert rows[0]["saturated"] is True
    assert rows[1]["signal"] == 900.5


def test_reader_refuses_incomplete_rows(tmp_path):
    path = tmp_path / "bands.csv"
    path.write_text("blot,sample,group,target,signal\nB,s,,T,1\n", encoding="utf-8")
    with pytest.raises(blot_quant.InputError, match="required"):
        blot_quant.read_bands(path)


def test_reader_explains_other_encodings(tmp_path):
    path = tmp_path / "bands.csv"
    path.write_bytes("blot,sample,group,target,signal\nB,ع,g,T,1\n".encode("cp1256"))
    with pytest.raises(blot_quant.InputError, match="CSV UTF-8"):
        blot_quant.read_bands(path)


def test_proportional_series_has_unit_slope():
    points = [(load, 100.0 * load) for load in (5, 10, 20, 40)]
    result = blot_linearity.fit(points)
    assert result["slope"] == pytest.approx(1.0)
    assert blot_linearity.linear_range(points) == 40


def test_saturating_series_reports_where_linearity_ends():
    points = [(load, 1000.0 * load / (load + 30.0)) for load in (1, 2, 4, 8, 16, 32, 64)]
    top = blot_linearity.linear_range(points)
    assert top is not None and top < 64
    result = blot_linearity.analyze([{"target": "T", "load": l, "signal": s} for l, s in points])
    assert any("stops tracking" in flag or "outside" in flag for flag in result["targets"][0]["flags"])


def test_too_few_levels_are_reported():
    rows = [{"target": "T", "load": load, "signal": load} for load in (1, 2)]
    assert "error" in blot_linearity.analyze(rows)["targets"][0]


def test_report_states_normalization_and_warnings():
    bands = membrane()
    bands[1]["saturated"] = True
    result = blot_quant.analyze(bands, "ACTB", "Ctrl")
    text = blot_report.build_markdown(json.loads(json.dumps(result)), None)
    assert text.index("saturated") < text.index("## Methods")
    assert "normalized to ACTB within each membrane" in text


def test_figure_and_report_are_written(tmp_path):
    pytest.importorskip("matplotlib")
    path = tmp_path / "result.json"
    path.write_text(json.dumps(blot_quant.analyze(membrane(), "ACTB", "Ctrl")), encoding="utf-8")
    out = tmp_path / "report"
    assert blot_report.main(["--result", str(path), "--out-dir", str(out)]) == 0
    assert (out / "blot_fold_change.png").stat().st_size > 5000


def test_cli_runs_the_documented_workflow(tmp_path):
    output = tmp_path / "out.json"
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "blot_quant.py"),
            "--input", str(SKILL_ROOT / "assets" / "example_bands.csv"),
            "--loading", "ACTB",
            "--control", "Control",
            "--output-json", str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "TP53 Treated vs Control" in done.stdout
    assert json.loads(output.read_text(encoding="utf-8"))["excluded"][0]["reason"] == "saturated"
    linearity = subprocess.run(
        [sys.executable, str(SCRIPTS / "blot_linearity.py"), "--input", str(SKILL_ROOT / "assets" / "example_linearity.csv")],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "ACTB" in linearity.stdout


def test_arabic_labels_survive_a_legacy_console(tmp_path):
    lines = ["blot,sample,group,target,signal"]
    for sample, group, ratio in (("c1", "ضابط", 1.0), ("c2", "ضابط", 1.1), ("t1", "دواء", 2.0), ("t2", "دواء", 2.2)):
        lines.append(f"B,{sample},{group},ACTB,1000")
        lines.append(f"B,{sample},{group},بروتين,{1000 * ratio}")
    path = tmp_path / "bands.csv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    done = subprocess.run(
        [sys.executable, str(SCRIPTS / "blot_quant.py"), "--input", str(path), "--loading", "ACTB", "--control", "ضابط"],
        capture_output=True,
        env={**os.environ, "PYTHONIOENCODING": "cp1252", "PYTHONUTF8": "0"},
    )
    assert done.returncode == 0, done.stderr
    assert b"fold" in done.stdout


def test_cli_reports_a_missing_file_cleanly(tmp_path):
    done = subprocess.run(
        [sys.executable, str(SCRIPTS / "blot_quant.py"), "--input", str(tmp_path / "x.csv"), "--loading", "A", "--control", "C"],
        capture_output=True,
        text=True,
    )
    assert done.returncode == 2
    assert "[FAIL]" in done.stderr


def test_each_membrane_is_centred_on_its_own_control_lanes():
    """A stronger antibody lot on one membrane must not leak into the comparison."""

    first = [
        band("B1", "C1", "Ctrl", "ACTB", 1000), band("B1", "C1", "Ctrl", "GOI", 500),
        band("B1", "T1", "Drug", "ACTB", 1000), band("B1", "T1", "Drug", "GOI", 1000),
    ]
    second = [
        band("B2", "C2", "Ctrl", "ACTB", 1000), band("B2", "C2", "Ctrl", "GOI", 4000),
        band("B2", "T2", "Drug", "ACTB", 1000), band("B2", "T2", "Drug", "GOI", 8000),
        band("B2", "C3", "Ctrl", "ACTB", 1000), band("B2", "C3", "Ctrl", "GOI", 4000),
    ]
    row = contrast(blot_quant.analyze(first + second, "ACTB", "Ctrl"))
    assert row["fold_change"] == pytest.approx(2.0, abs=1e-3)
