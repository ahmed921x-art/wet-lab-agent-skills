---
name: elisa-analysis
description: Analyzes ELISA and other plate immunoassay readouts. Fits 4PL or 5PL standard curves with weighting, back-calculates sample concentrations with dilution correction, derives a working range from standard recovery and replicate CV, and censors values outside it instead of extrapolating. Handles sandwich and competitive formats, QC recovery, dilution linearity, and hook-effect flags. Use for ELISA, OD or absorbance plate readings, standard curves, four-parameter logistic fits, cytokine or biomarker concentrations, kit results, or converting optical density to concentration.
license: MIT
compatibility: Requires Python 3.11+ with numpy and scipy; figures need matplotlib. Tested on Python 3.13 with scipy 1.18.1 and matplotlib 3.11.2. No network access.
allowed-tools: Read Write Edit Bash
metadata:
  version: "1.1"
  skill-author: ahmed921x-art
  last-reviewed: "2026-11-06"
---

# ELISA Analysis

## When to use

Use this skill to turn one plate of optical densities into concentrations, or to work out why a standard
curve or a sample looks wrong. Run it once per plate: a curve belongs to the plate it was read on, so
pooling plates is refused. Formal assay validation belongs to `analytical-method-validation`. Comparing
concentrations between groups belongs to `statistical-analysis`. Multiplex bead assays can use the same
model, but this skill was not tested on them.

## The three rules

**1. Never report a number from outside the working range.** A curve flattens at both ends, where a small
change in OD maps to a huge change in concentration. Those values are returned as `<LLOQ` or `>ULOQ`
bounds with the dilution applied. Re-assay high samples at a higher dilution.

**2. Judge a curve by back-calculated recovery, not by R^2.** A logistic curve has R^2 near 0.99 even when
the low end is wrong by 50%. The test that matters is whether each standard, read as an unknown, comes back
near its nominal value.

**3. Do not delete standards to improve a fit.** Removing a point is a decision that needs a recorded
reason, such as a pipetting error noted at the bench. `--exclude-level` does this openly and the report
records it. Nothing is ever dropped automatically.

## Workflow

### 1. Shape the data

One row per well, one plate per file:

| Column | Required | Meaning |
| --- | --- | --- |
| `kind` | yes | `standard`, `blank`, `sample`, or `qc` |
| `id` | yes | Well label; replicate wells share an id |
| `conc` | standards, qc | Nominal concentration |
| `dilution` | no | Dilution factor, for example 10 for 1:10; default 1 |
| `od` | yes | Measured OD, after wavelength correction if used |

Keep each well as its own row, because the fit uses wells and replicate CV needs them. A QC's `conc` is its
nominal neat value.

Plate readers export an 8x12 grid. `scripts/elisa_plate_import.py` joins that grid with a layout grid and
skips reader preamble. It can also subtract a reference wavelength:

```bash
python scripts/elisa_plate_import.py --od assets/example_od450.csv --reference-od assets/example_od630.csv \
  --layout assets/example_layout.csv --output wells.csv
```

Layout cells are `std:1000`, `blank`, `qc:QC1:200`, `P01`, or `P01:10` for 1:10. Empty cells are unused.

### 2. Fit

```bash
python scripts/elisa_fit.py --input assets/example_plate.csv --output-json result.json --output-csv samples.csv
```

The curve is `y = d + (a - d) / (1 + (x / c)^b)^g`, where `a` is the response at zero, `d` the response at
saturation, `c` the midpoint, and `b` the slope. `g` is 1 in 4PL. The direction is detected, so competitive
assays need no flag.

- **Model:** start with 4PL. Use `--model 5pl` only for a visibly asymmetric curve with at least seven
  levels. With fewer levels the extra parameter follows the noise.
- **Weighting:** the default `1/y2` is used because OD variance grows with signal, and an unweighted fit
  sacrifices the low end, usually where samples sit.
- **Comparing fits:** `--compare` prints standard recovery for each model and weight, and hides samples on
  purpose, so the choice cannot be steered by the answer. Small differences are noise; keep the default
  unless another option clearly fixes a failing end of the curve.
- **Blank:** `blank` wells are averaged and subtracted. With none, a zero standard is used on a rising
  curve. On a competitive curve the zero standard is the maximum signal, so it is not a blank.
  `--no-blank-subtract` fits raw OD, when the kit's instructions say so.

### 3. Read the output

The working range is the longest contiguous run of standards with recovery within 20% and replicate CV up
to 20%. Both limits can be changed, before analysis, with `--max-recovery-error` and `--max-cv`. ICH M10
allows 25% at the range limits for ligand binding assays, so these defaults are slightly stricter there.
A regulated study follows its own validated criteria.

| Sample status | Meaning | Action |
| --- | --- | --- |
| `ok` | Inside the range | Report `concentration`, already corrected for dilution |
| `above_range` | Over the ULOQ or the upper asymptote | Dilute and re-assay; report the `>` bound until then |
| `below_range` | Under the LLOQ or the lower asymptote | Report the `<` bound, or re-assay less diluted |
| `no_working_range` | No standard passed | The plate needs repeating |

Check these before using any concentration:

- Replicate CV above 20% on a sample flags a pipetting or washing issue.
- Failing QC recovery can invalidate the plate, whatever the standards did.
- Dilution linearity is computed when one sample is read at two or more dilutions. Disagreement after
  correction suggests a matrix effect or non-parallelism, and the concentration should not be trusted.
- A top standard below the one under it suggests hook effect or saturation.
- An OD above `--max-od` may be past the reader's linear range. The default of 3.5 is a working value;
  check the instrument's specification.

`concentration_ci` is a 95% interval that combines curve uncertainty with read-out noise. It is wider near
the ends of the curve, so a sample with a wide interval may deserve a re-assay at another dilution.

### 4. Handle censored values downstream

Do not replace `<LLOQ` with zero or with LLOQ/2 without saying so, since that choice can change a group
comparison. State the number of censored samples per group. If many are censored, use a method built for
censored data, and involve `statistical-analysis` or `statsmodels`.

### 5. Report

```bash
python scripts/elisa_report.py --result result.json --out-dir report
```

This writes a 300 dpi standard-curve figure with the working range marked, plus `report.md`, with warnings
first and draft Methods and result tables. Then state the kit and lot, model and weighting, standard range, working range, blank handling, replicate
scheme, QC recovery, and how censored values were treated.

## Limits

- Concentrations are comparable only within what the standard calibrates. Different kits for one analyte
  often disagree in absolute terms.
- The interval assumes the fitted weighting describes the noise, and one plate. Between-plate variation is
  not in it, so compare results across plates only with a shared QC.
- Relative potency needs a parallel-curves analysis, which this skill does not perform.

## Sources

Method details and citations are in `references/curve-fitting-and-sources.md`.
