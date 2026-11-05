---
name: western-blot-densitometry
description: Analyzes Western blot densitometry tables. Normalizes bands to a housekeeping protein or total-protein stain within each membrane, combines repeat blots, compares groups on the log scale with confidence intervals and multiple-testing correction, and flags saturated bands, uneven loading, and loading controls that move with treatment. Also checks linear range from a dilution series. Use for Western blot or immunoblot quantification, ImageJ or Fiji gel analysis results, band density, beta-actin or GAPDH normalization, phospho to total ratios, or fold change from blots.
license: MIT
compatibility: Requires Python 3.11+ with scipy; figures need matplotlib. Tested on Python 3.13 with scipy 1.18.1 and matplotlib 3.11.2. No network access.
allowed-tools: Read Write Edit Bash
metadata:
  version: "1.0"
  skill-author: ahmed921x-art
  last-reviewed: "2026-11-06"
---

# Western Blot Densitometry

## When to use

Use this skill after bands have been measured, in ImageJ/Fiji, Image Lab, or similar, to turn the
measurements into normalized results. It does not read blot images or measure bands. For protein
quantification by mass spectrometry use `pyopenms`.

## The three rules

**1. Raw density only means something within one membrane.** Exposure, transfer, and antibody batch all
change between blots. So each band is divided by its loading control, then compared with control lanes on
the same membrane. A blot with no control lane cannot be included.

**2. A saturated band is not a measurement.** Clipped pixels understate the signal, and strong bands often
saturate while the loading control still looks fine. Saturated bands are left out by default. Abundant
housekeeping proteins are frequently saturated at normal loads, which makes them poor normalizers. That is
why the linearity check exists.

**3. The sample is the unit, not the blot.** Running one lysate on three membranes is technical repetition.
Those are averaged, so n is the number of independent biological samples.

## Workflow

### 1. Check linearity once per antibody

Run a dilution series of a pooled lysate and measure it. Columns are `target`, `load`, `signal`.

```bash
python scripts/blot_linearity.py --input assets/example_linearity.csv
```

A slope near 1 on log-log axes means signal doubles when load doubles. `linear up to` gives the highest load
that stays linear. Run target and loading control within this range, or the ratio will mislead.

### 2. Shape the data

One row per band:

| Column | Required | Meaning |
| --- | --- | --- |
| `blot` | yes | Membrane id |
| `sample` | yes | Biological sample id, the same across membranes |
| `group` | yes | Condition |
| `target` | yes | Protein, including the loading control |
| `signal` | yes | Measured density; ImageJ's gel analyzer gives this as `Area` |
| `background` | no | Subtracted from signal; omit if already subtracted |
| `saturated` | no | `yes` for a band with clipped pixels |
| `pair` | no | Shared id for matched samples; switches to a paired test |

For phospho/total ratios, treat each pair as a target and its total, and run once per pair with
`--loading` set to the total protein.

### 3. Quantify

```bash
python scripts/blot_quant.py --input assets/example_bands.csv --loading ACTB --control Control --output-json result.json
```

### 4. Read the warnings

| Warning | Meaning | What to do |
| --- | --- | --- |
| saturated bands excluded | Clipped bands were dropped | Re-image at shorter exposure |
| loading control CV | Lanes differ by over 30% | Inspect the blot; re-run if transfer was uneven |
| loading control differs between groups | Treatment may change it, or loading differed | Use total-protein normalization |
| no control lane | That membrane cannot be normalized | Include control lanes on every membrane |
| descriptive | Under three samples in a group | Report values without a test |

### 5. Draft the report

```bash
python scripts/blot_report.py --result result.json --out-dir report
```

This writes a figure and `report.md`, with warnings first and draft Methods and Results text. Keep uncropped
blots for every panel, since journals ask for them. Do not adjust contrast on images used for measurement.

## Limits

- Densitometry is semi-quantitative. Treat a fold change under about 1.5 with caution, whatever its p-value.
- Comparison is with a single control group. For factorial designs, take log2 values from `result.json` and
  fit a model with `statsmodels`.
- The script trusts the `saturated` column. Check your imager's saturation display, because a flag that was
  never set cannot be caught here.

## Sources

Rationale, thresholds, and citations are in `references/methods-and-sources.md`.
