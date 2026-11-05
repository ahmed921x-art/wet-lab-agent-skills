---
name: qpcr-analysis
description: Analyzes RT-qPCR plate exports for relative gene expression. Runs well and plate QC, ΔΔCt (Livak) and efficiency-corrected (Pfaffl) quantification with one or more reference genes, standard-curve efficiency, confidence intervals, and multiple-testing correction on biological replicates. Use for qPCR, RT-PCR, Ct or Cq tables, fold change, 2^-ΔΔCt, housekeeping or reference gene normalization, primer efficiency, melt or NTC checks, or QuantStudio, CFX, and LightCycler exports, including when a user only says their PCR results need analysing.
license: MIT
compatibility: Requires Python 3.11+ with scipy; figures need matplotlib and .xlsx import needs openpyxl. Tested on Python 3.13 with scipy 1.18.1 and matplotlib 3.11.2. No network access.
allowed-tools: Read Write Edit Bash
metadata:
  version: "1.1"
  skill-author: ahmed921x-art
  last-reviewed: "2026-11-06"
---

# qPCR Analysis

## When to use

Use this skill to turn exported Ct (Cq) values into defensible relative expression results, to check a
primer pair's efficiency from a dilution series, or to find out why a plate looks wrong. For primer
design use `primer-design`. For RNA-seq counts use `bulk-rnaseq`. Absolute copy-number assays and digital
PCR are out of scope.

## The three rules

**1. The unit of analysis is the biological sample.** Technical wells are averaged first and never become
n. Three wells from one mouse is n = 1, and treating them as three inflates significance.

**2. Test on the log scale, then back-transform.** Ct differences are roughly normal; fold changes are
not. Means, tests, and intervals are computed on log2 relative quantity, then converted, so the interval
around a fold change is asymmetric. An error bar built from the SD of 2^-ΔΔCt values is wrong.

**3. An undetermined well is not a number.** Replacing it with 40 or 0 invents a result. The scripts leave
such samples out, count them, and report the detection rate. If a target is off in one group, report that
directly instead of forcing a fold change.

## Workflow

### 1. Shape the data

Convert an instrument export directly. `scripts/qpcr_import.py` finds the results table in QuantStudio,
CFX, and similar exports (CSV, TSV, or .xlsx), skips metadata, and drops omitted wells and standards with
a count. Groups come from a sample sheet, or from `--group-pattern` on the sample name.

```bash
python scripts/qpcr_import.py --input assets/example_quantstudio_export.csv \
  --sample-sheet assets/example_sample_sheet.csv --output wells.csv
```

An omitted well was removed by someone on the instrument, so the report asks whether that was justified.
The importer was built from documented export layouts, so check the first run against the original file.

`scripts/qpcr_relative.py` needs one row per well:

| Column | Required | Meaning |
| --- | --- | --- |
| `sample` | yes | Biological sample id |
| `group` | yes for samples | Condition; empty for NTC and no-RT wells |
| `target` | yes | Gene or assay |
| `ct` | yes | Number, or `Undetermined`, `N/A`, or blank |
| `task` | no | `unknown` (default), `ntc`, or `nrt` |
| `pair` | no | Shared id for matched samples; switches to a paired test |

For any other export, `references/methods-and-sources.md` maps common columns. A Ct of 0 is read as undetermined, because some instruments export it that way. Text the
parser does not recognize stops the run instead of being guessed at.

### 2. Measure efficiency

Run this for each primer pair, once per assay conditions, from a serial dilution with at least five
levels. The JSON is the input to the next step.

```bash
python scripts/qpcr_efficiency.py --input assets/example_dilution_series.csv --output-json efficiency.json
```

Columns are `target`, `quantity`, `ct`. A curve is flagged outside 90-110% efficiency, under R^2 0.98, or
below five levels. Above 110% usually means inhibition at the concentrated end, or primer dimers, not a
better reaction.

### 3. Choose references

Rank three or more candidates on samples that cover every condition:

```bash
python scripts/qpcr_references.py --input wells.csv --candidate GAPDH --candidate ACTB --candidate HK3
```

This runs the geNorm stepwise ranking and pairwise variation. geNorm measures agreement between genes, so
it cannot see a gene that shifts in one condition. The script therefore also tests each candidate against
the others between groups, and leaves movers out of the recommendation. Choose references here, not by
seeing which one gives the desired fold change.

### 4. Quantify

```bash
python scripts/qpcr_relative.py \
  --input assets/example_plate.csv \
  --reference GAPDH --reference ACTB \
  --control Vehicle \
  --efficiency-json efficiency.json \
  --output-json result.json --output-csv contrasts.csv
```

Without `--efficiency-json` every efficiency is assumed to be 100% (Livak), and the output says so. That
assumption has to be checked, not just stated.

Each contrast reports fold change, its interval, Welch or paired t-test, and a Benjamini-Hochberg p-value
across all contrasts in the run.

### 5. Read the warnings before the numbers

| Code | Meaning | What to do |
| --- | --- | --- |
| `control_contamination` | NTC or no-RT signal within 5 cycles of samples | Treat the target as unreliable; find the source |
| `control_signal` | NTC or no-RT signal, far from samples | Check the melt curve; often primer dimers |
| `reference_shift` | Reference differs between groups by over 0.5 cycles | The treatment may regulate it; choose or add references |
| `reference_unstable` | geNorm M at or above 1.5 | Replace the reference |
| `well_qc` | Technical SD over 0.5, Ct over 35, or lost wells | Inspect the wells; do not delete them to tighten a result |
| `target_not_quantified` | Samples left out for a target | Report the detection rate beside the fold change |
| `low_n` | Fewer than three biological samples in a group | Treat as descriptive |

Thresholds are working defaults, not universal limits, and `--max-well-sd`, `--ct-ceiling`, and
`--min-wells` change them. Set them before looking at the outcome. Provenance is listed in
`references/methods-and-sources.md`.

### 6. Report

```bash
python scripts/qpcr_report.py --result result.json --efficiency-json efficiency.json --out-dir report
```

This writes a 300 dpi figure of every biological sample with its interval, plus `report.md`, with open
warnings first and draft Methods and Results text. Edit it into your own words. Then state the model (Livak or Pfaffl), measured efficiencies, references and their stability, n of biological
samples per group, the test, and the interval. Give fold change with its interval, not the mean with SD.
A fold change near 1.5 sits close to the noise of a typical assay even when the p-value is small, so judge
effect size as well as significance. `references/methods-and-sources.md` lists the MIQE items to include.

## Limits

- Fold changes are relative to the control group mean, so the control has its own spread. The contrast's
  interval already includes it; do not show the control as exactly 1 with no error.
- The script compares every group with one control. For factorial designs, take log2 values from
  `result.json` and fit a linear model with the `statsmodels` skill.
- Results do not establish reference stability from a single experiment. Validate references for the
  tissue and treatment.
- Melt curves and amplification plots are not read here. Check them on the instrument software.

## Sources

Methods, thresholds, export mapping, and citations are in `references/methods-and-sources.md`.
