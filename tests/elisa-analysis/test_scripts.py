#!/usr/bin/env python3
"""Tests for elisa-analysis on noise-free curves with known parameters."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("scipy")

import numpy as np  # noqa: E402
import skill_contract  # noqa: E402

SKILL_ROOT = Path(__file__).resolve().parents[2] / "skills" / "elisa-analysis"
SCRIPTS = SKILL_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import elisa_fit  # noqa: E402
import elisa_plate_import  # noqa: E402
import elisa_report  # noqa: E402

CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)

LEVELS = (15.6, 31.25, 62.5, 125.0, 250.0, 500.0, 1000.0)


def response(x, a=0.02, d=2.5, c=150.0, b=1.2, g=1.0):
    return float(elisa_fit.logistic(np.array([float(x)]), a, d, c, b, g)[0])


def standards(levels=LEVELS, **parameters):
    wells = []
    for level in levels:
        for _ in range(2):
            wells.append(well("standard", f"S{level}", conc=level, od=response(level, **parameters)))
    return wells


def well(kind, name, *, conc=None, dilution=1.0, od=0.0):
    return {"kind": kind, "id": name, "conc": conc, "dilution": dilution, "od": od}


def unknown(result, name):
    return next(entry for entry in result["unknowns"] if entry["id"] == name)


def write_plate(path: Path, rows: list[str]) -> Path:
    path.write_text("kind,id,conc,dilution,od\n" + "\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_invert_is_the_inverse_of_the_model():
    parameters = {"a": 0.02, "d": 2.5, "c": 150.0, "b": 1.2, "g": 0.8}
    for x in (5.0, 80.0, 900.0):
        assert elisa_fit.invert(response(x, **parameters), **parameters) == pytest.approx(x)


@pytest.mark.parametrize("y", [0.01, 2.6])
def test_invert_returns_none_outside_the_asymptotes(y):
    assert elisa_fit.invert(y, 0.02, 2.5, 150.0, 1.2) is None


def test_fit_recovers_known_4pl_parameters():
    result = elisa_fit.analyze(standards(), no_blank_subtract=True)
    curve = result["curve"]
    assert curve["direction"] == "sandwich"
    assert curve["a"] == pytest.approx(0.02, abs=1e-4)
    assert curve["d"] == pytest.approx(2.5, rel=1e-4)
    assert curve["c"] == pytest.approx(150.0, rel=1e-4)
    assert curve["b"] == pytest.approx(1.2, rel=1e-4)
    assert all(level["recovery_percent"] == pytest.approx(100.0, abs=0.05) for level in result["standards"])
    assert (result["lloq"], result["uloq"]) == (15.6, 1000.0)


@pytest.mark.parametrize("weight", ["1/y2", "1/y", "none"])
def test_every_weighting_recovers_a_noise_free_curve(weight):
    curve = elisa_fit.analyze(standards(), weight=weight, no_blank_subtract=True)["curve"]
    assert curve["c"] == pytest.approx(150.0, rel=1e-3)


def test_competitive_curve_falls_with_concentration():
    parameters = {"a": 2.4, "d": 0.1, "c": 90.0, "b": 1.1}
    wells = standards(**parameters) + [
        well("sample", "S", od=response(100.0, **parameters)),
        well("sample", "S", od=response(100.0, **parameters)),
    ]
    result = elisa_fit.analyze(wells, no_blank_subtract=True)
    assert result["curve"]["direction"] == "competitive"
    assert unknown(result, "S")["concentration"] == pytest.approx(100.0, rel=1e-3)


def test_5pl_recovers_asymmetry():
    levels = (5.0, 15.0, 45.0, 135.0, 405.0, 1215.0, 3645.0, 10935.0)
    result = elisa_fit.analyze(standards(levels=levels, g=0.6), model="5pl", no_blank_subtract=True)
    assert result["curve"]["g"] == pytest.approx(0.6, rel=0.02)


def test_5pl_needs_more_levels():
    with pytest.raises(elisa_fit.InputError, match="cannot support 5pl"):
        elisa_fit.analyze(standards(levels=(10.0, 40.0, 160.0, 640.0, 2560.0)), model="5pl")


def test_three_levels_cannot_fit_4pl():
    with pytest.raises(elisa_fit.InputError, match="cannot support 4pl"):
        elisa_fit.analyze(standards(levels=(10.0, 100.0, 1000.0)))


def test_blank_is_subtracted_before_fitting():
    shift = 0.06
    wells = [dict(row, od=row["od"] + shift) for row in standards()]
    wells += [well("blank", "B", od=shift), well("blank", "B", od=shift)]
    result = elisa_fit.analyze(wells)
    assert result["blank_od"] == pytest.approx(shift)
    assert result["curve"]["a"] == pytest.approx(0.02, abs=1e-3)


def test_zero_standard_acts_as_blank_on_a_rising_curve():
    wells = standards() + [well("standard", "S0", conc=0.0, od=0.05)]
    assert elisa_fit.analyze(wells)["blank_od"] == pytest.approx(0.05)


def test_zero_standard_is_not_a_blank_on_a_competitive_curve():
    parameters = {"a": 2.4, "d": 0.1, "c": 90.0, "b": 1.1}
    wells = standards(**parameters) + [well("standard", "S0", conc=0.0, od=2.4)]
    result = elisa_fit.analyze(wells)
    assert result["blank_od"] == 0.0
    assert any("competitive" in warning for warning in result["warnings"])


def test_sample_concentration_is_multiplied_by_dilution():
    wells = standards() + [well("sample", "S", dilution=10.0, od=response(200.0))] * 2
    entry = unknown(elisa_fit.analyze(wells, no_blank_subtract=True), "S")
    assert entry["well_concentration"] == pytest.approx(200.0, rel=1e-3)
    assert entry["concentration"] == pytest.approx(2000.0, rel=1e-3)


def test_values_past_the_top_standard_are_censored_not_extrapolated():
    wells = standards() + [well("sample", "HI", dilution=5.0, od=response(5000.0))] * 2
    entry = unknown(elisa_fit.analyze(wells, no_blank_subtract=True), "HI")
    assert entry["status"] == "above_range"
    assert entry["concentration"] is None
    assert entry["censored"] == ">5000.0"
    assert "re-assay at a higher dilution" in entry["flags"]


def test_signal_beyond_the_asymptote_is_above_range():
    wells = standards() + [well("sample", "SAT", od=9.0)] * 2
    entry = unknown(elisa_fit.analyze(wells, no_blank_subtract=True), "SAT")
    assert entry["status"] == "above_range"
    assert entry["well_concentration"] is None


def test_signal_under_the_bottom_asymptote_is_below_range():
    wells = standards() + [well("sample", "LOW", od=0.001)] * 2
    assert unknown(elisa_fit.analyze(wells, no_blank_subtract=True), "LOW")["status"] == "below_range"


def test_failing_low_standard_raises_the_lloq():
    wells = standards()
    for row in wells:
        if row["conc"] == 15.6:
            row["od"] *= 2.5
    result = elisa_fit.analyze(wells, no_blank_subtract=True)
    assert result["lloq"] > 15.6


def test_imprecise_level_leaves_the_working_range():
    wells = standards()
    top = [row for row in wells if row["conc"] == 1000.0]
    top[0]["od"] *= 0.6
    top[1]["od"] *= 1.4
    assert elisa_fit.analyze(wells, no_blank_subtract=True)["uloq"] < 1000.0


def test_exclusion_is_recorded_and_unknown_levels_are_refused():
    result = elisa_fit.analyze(standards(), exclude_levels=(1000.0,), no_blank_subtract=True)
    assert any("excluded by the analyst" in warning for warning in result["warnings"])
    assert result["uloq"] == 500.0
    with pytest.raises(elisa_fit.InputError, match="matches no standard"):
        elisa_fit.analyze(standards(), exclude_levels=(7.0,))


def test_qc_recovery_is_reported_and_failures_warn():
    wells = standards()
    wells += [well("qc", "Q1", conc=200.0, od=response(200.0))] * 2
    wells += [well("qc", "Q2", conc=200.0, od=response(400.0))] * 2
    result = elisa_fit.analyze(wells, no_blank_subtract=True)
    assert unknown(result, "Q1")["recovery_percent"] == pytest.approx(100.0, abs=0.2)
    assert any("Q2" in warning for warning in result["warnings"])


def test_dilution_linearity_flags_disagreement():
    wells = standards()
    wells += [well("sample", "OK", dilution=10.0, od=response(200.0))] * 2
    wells += [well("sample", "OK", dilution=20.0, od=response(100.0))] * 2
    wells += [well("sample", "BAD", dilution=10.0, od=response(200.0))] * 2
    wells += [well("sample", "BAD", dilution=20.0, od=response(200.0))] * 2
    result = elisa_fit.analyze(wells, no_blank_subtract=True)
    by_id = {row["id"]: row for row in result["dilution_linearity"]}
    assert by_id["OK"]["cv_percent"] == pytest.approx(0.0, abs=0.5)
    assert by_id["BAD"]["cv_percent"] > 20
    assert any("BAD" in warning for warning in result["warnings"])


def test_replicate_disagreement_is_flagged():
    wells = standards() + [
        well("sample", "S", od=response(100.0) * 0.6),
        well("sample", "S", od=response(100.0) * 1.4),
    ]
    assert any("CV" in flag for flag in unknown(elisa_fit.analyze(wells, no_blank_subtract=True), "S")["flags"])


def test_small_concentration_scales_keep_their_digits():
    parameters = {"a": 0.02, "d": 2.5, "c": 0.015, "b": 1.2}
    levels = tuple(level / 10000 for level in LEVELS)
    wells = standards(levels=levels, **parameters)
    wells += [well("sample", "S", od=response(0.012, **parameters))] * 2
    result = elisa_fit.analyze(wells, no_blank_subtract=True)
    assert result["curve"]["c"] == pytest.approx(0.015, rel=1e-3)
    assert unknown(result, "S")["concentration"] == pytest.approx(0.012, rel=1e-3)


def test_hook_like_top_standard_warns():
    wells = standards()
    for row in wells:
        if row["conc"] == 1000.0:
            row["od"] = response(250.0)
    result = elisa_fit.analyze(wells, no_blank_subtract=True)
    assert any("hook" in warning for warning in result["warnings"])


def test_two_plates_in_one_file_are_refused(tmp_path):
    path = tmp_path / "plates.csv"
    path.write_text(
        "kind,id,conc,dilution,od,plate\nstandard,S1,10,,0.1,A\nstandard,S1,10,,0.1,B\n",
        encoding="utf-8",
    )
    with pytest.raises(elisa_fit.InputError, match="run each plate separately"):
        elisa_fit.read_wells(path)


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("mystery,X,1,,0.1", "unknown kind"),
        ("standard,X,,,0.1", "nominal conc"),
        ("sample,X,,0.1,0.1", "dilution factor"),
        ("sample,X,,1,high", "not a number"),
        ("sample,X,,1,0.1,extra", "more fields"),
    ],
)
def test_bad_rows_name_the_problem(tmp_path, row, message):
    with pytest.raises(elisa_fit.InputError, match=message):
        elisa_fit.read_wells(write_plate(tmp_path / "plate.csv", [row]))


def test_cli_runs_the_documented_workflow(tmp_path):
    output = tmp_path / "out.json"
    table = tmp_path / "out.csv"
    done = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "elisa_fit.py"),
            "--input",
            str(SKILL_ROOT / "assets" / "example_plate.csv"),
            "--output-json",
            str(output),
            "--output-csv",
            str(table),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "working range" in done.stdout
    payload = json.loads(output.read_text(encoding="utf-8"))
    status = {(row["id"], row["dilution"]): row["status"] for row in payload["unknowns"]}
    assert status[("P03", 10.0)] == "above_range"
    assert status[("P04", 2.0)] == "below_range"
    assert status[("P01", 10.0)] == "ok"
    assert table.read_text(encoding="utf-8").startswith("kind,id,dilution")


def test_cli_reports_a_missing_file_cleanly(tmp_path):
    done = subprocess.run(
        [sys.executable, str(SCRIPTS / "elisa_fit.py"), "--input", str(tmp_path / "absent.csv")],
        capture_output=True,
        text=True,
    )
    assert done.returncode == 2
    assert "[FAIL]" in done.stderr


def test_regional_export_with_arabic_names_survives_a_legacy_console(tmp_path):
    def comma(value: float) -> str:
        return f"{value:.4f}".replace(".", ",")

    lines = ["kind\tid\tconc\tdilution\tod"]
    for level in LEVELS:
        for _ in range(2):
            lines.append(f"standard\tS\t{comma(level)}\t\t{comma(response(level))}")
    lines += [f"sample\t\u0645\u0631\u064a\u06361\t\t10\t{comma(response(100.0))}"] * 2
    path = tmp_path / "plate.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    done = subprocess.run(
        [sys.executable, str(SCRIPTS / "elisa_fit.py"), "--input", str(path), "--no-blank-subtract"],
        capture_output=True,
        env={**os.environ, "PYTHONIOENCODING": "cp1252", "PYTHONUTF8": "0"},
    )
    assert done.returncode == 0, done.stderr
    assert b"[ok]" in done.stdout


def test_other_encodings_get_a_clear_message(tmp_path):
    path = tmp_path / "plate.csv"
    path.write_bytes("kind,id,od\nsample,\u0639,1\n".encode("cp1256"))
    with pytest.raises(elisa_fit.InputError, match="CSV UTF-8"):
        elisa_fit.read_wells(path)


# ------------------------------------------------------------------ import


def grid_file(path, cells: dict[str, str], preamble: str = "") -> Path:
    rows = [preamble] if preamble else []
    rows.append("," + ",".join(str(i) for i in range(1, 13)))
    for row in "ABCDEFGH":
        rows.append(row + "," + ",".join(cells.get(f"{row}{i}", "") for i in range(1, 13)))
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_grid_reader_skips_reader_preamble(tmp_path):
    grid = elisa_plate_import.read_grid(grid_file(tmp_path / "g.csv", {"A1": "0.5", "H12": "2.1"}, "Reader 450 nm,,,"))
    assert grid["A1"] == "0.5"
    assert grid["H12"] == "2.1"
    assert len(grid) == 96


@pytest.mark.parametrize(
    ("cell", "expected"),
    [
        ("std:1000", {"kind": "standard", "conc": "1000"}),
        ("blank", {"kind": "blank"}),
        ("qc:Q1:200", {"kind": "qc", "id": "Q1", "conc": "200", "dilution": "1"}),
        ("qc:Q1:200:5", {"kind": "qc", "dilution": "5"}),
        ("P01", {"kind": "sample", "id": "P01", "dilution": "1"}),
        ("P01:10", {"kind": "sample", "id": "P01", "dilution": "10"}),
    ],
)
def test_layout_cells(cell, expected):
    parsed = elisa_plate_import.parse_cell(cell, "A1")
    assert expected.items() <= parsed.items()


def test_empty_layout_cell_is_unused():
    assert elisa_plate_import.parse_cell("  ", "A1") is None


@pytest.mark.parametrize("cell", ["std:", "std:1:2", "qc:Q1", "P:1:2"])
def test_bad_layout_cells_are_refused(cell):
    with pytest.raises(elisa_fit.InputError):
        elisa_plate_import.parse_cell(cell, "A1")


def test_reference_wavelength_is_subtracted():
    wells = elisa_plate_import.combine({"A1": "1.500"}, {"A1": "P1"}, {"A1": "0.040"})
    assert wells[0]["od"] == "1.4600"


def test_layout_without_a_reading_is_refused():
    with pytest.raises(elisa_fit.InputError, match="OD grid is empty"):
        elisa_plate_import.combine({"A1": ""}, {"A1": "P1"})


def test_grid_without_column_numbers_is_refused(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("a,b\nc,d\n", encoding="utf-8")
    with pytest.raises(elisa_fit.InputError, match="no column header"):
        elisa_plate_import.read_grid(path)


def test_example_plate_imports_and_fits(tmp_path):
    output = tmp_path / "wells.csv"
    assets = SKILL_ROOT / "assets"
    assert (
        elisa_plate_import.main(
            [
                "--od", str(assets / "example_od450.csv"),
                "--reference-od", str(assets / "example_od630.csv"),
                "--layout", str(assets / "example_layout.csv"),
                "--output", str(output),
            ]
        )
        == 0
    )
    result = elisa_fit.analyze(elisa_fit.read_wells(output))
    assert result["lloq"] == 15.6
    assert unknown(result, "P03")["status"] == "above_range"
    assert unknown(result, "QC1")["recovery_percent"] == pytest.approx(100, abs=10)


# ----------------------------------------------------------------- intervals


def noisy_plate(rng, cv=0.04):
    wells = []
    for level in LEVELS:
        for _ in range(2):
            wells.append(well("standard", f"S{level}", conc=level, od=response(level) * (1 + rng.normal(0, cv))))
    for _ in range(2):
        wells.append(well("sample", "S", od=response(120.0) * (1 + rng.normal(0, cv))))
    return wells


def test_interval_brackets_the_estimate_and_scales_with_dilution():
    rng = np.random.default_rng(1)
    wells = noisy_plate(rng)
    wells += [dict(row, id="D", dilution=10.0) for row in wells if row["id"] == "S"]
    result = elisa_fit.analyze(wells, no_blank_subtract=True)
    plain, diluted = unknown(result, "S"), unknown(result, "D")
    low, high = plain["concentration_ci"]
    assert low < plain["concentration"] < high
    assert diluted["concentration_ci"][0] == pytest.approx(low * 10, rel=1e-2)


def test_interval_is_wider_with_noisier_data():
    def width(cv):
        result = elisa_fit.analyze(noisy_plate(np.random.default_rng(3), cv=cv), no_blank_subtract=True)
        low, high = unknown(result, "S")["concentration_ci"]
        return high / low

    assert width(0.08) > width(0.02)


def test_interval_coverage_is_close_to_nominal():
    rng = np.random.default_rng(20261106)
    trials, hits = 0, 0
    for _ in range(120):
        result = elisa_fit.analyze(noisy_plate(rng), no_blank_subtract=True)
        entry = unknown(result, "S")
        if entry.get("concentration_ci"):
            trials += 1
            low, high = entry["concentration_ci"]
            hits += low <= 120.0 <= high
    assert trials >= 110
    assert 0.87 <= hits / trials <= 1.0


def test_internal_fit_details_stay_out_of_the_output():
    curve = elisa_fit.analyze(standards(), no_blank_subtract=True)["curve"]
    assert "covariance" not in curve and "theta" not in curve


# ---------------------------------------------------------------- comparison


def test_compare_covers_every_model_and_weight():
    rows = elisa_fit.compare_fits(standards(levels=(5.0, 15.0, 45.0, 135.0, 405.0, 1215.0, 3645.0)), no_blank_subtract=True)
    assert len(rows) == 6
    assert all(row["worst_recovery_error_percent"] < 0.5 for row in rows)


def test_compare_explains_a_model_it_cannot_fit():
    rows = elisa_fit.compare_fits(standards(levels=(10.0, 40.0, 160.0, 640.0, 2560.0)), no_blank_subtract=True)
    assert all("cannot support 5pl" in row["error"] for row in rows if row["model"] == "5pl")


def test_cli_compare_ignores_samples(tmp_path):
    done = subprocess.run(
        [sys.executable, str(SCRIPTS / "elisa_fit.py"), "--input", str(SKILL_ROOT / "assets" / "example_plate.csv"), "--compare"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "4pl 1/y2" in done.stdout
    assert "P01" not in done.stdout


# ------------------------------------------------------------------ report


def test_report_marks_censored_values_and_warnings():
    wells = standards() + [well("sample", "HI", dilution=5.0, od=response(5000.0))] * 2
    result = elisa_fit.analyze(wells, no_blank_subtract=True)
    text = elisa_report.build_markdown(json.loads(json.dumps(result)), None)
    assert ">5000.0" in text
    assert "censored" in text
    assert "No blank was subtracted" in text


def test_figure_is_written(tmp_path):
    pytest.importorskip("matplotlib")
    wells = standards() + [well("sample", "S", od=response(100.0))] * 2
    path = tmp_path / "result.json"
    path.write_text(json.dumps(elisa_fit.analyze(wells, no_blank_subtract=True)), encoding="utf-8")
    out = tmp_path / "report"
    assert elisa_report.main(["--result", str(path), "--out-dir", str(out)]) == 0
    assert (out / "standard_curve.png").stat().st_size > 5000
