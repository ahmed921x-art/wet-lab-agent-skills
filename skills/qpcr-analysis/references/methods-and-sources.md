# Methods, thresholds, and sources

## Model

For sample `s`, with amplification factor `E` (2.0 is 100% efficiency):

```text
L(s, target) = mean over references[ log2(E_ref) * Ct_ref(s) ] - log2(E_target) * Ct_target(s)
log2 fold change(group vs control) = mean L(group) - mean L(control)
```

- With every `E = 2` and one reference, `L = -dCt` and the contrast is `-ddCt` (Livak and Schmittgen 2001).
- With measured `E` and one reference, the contrast equals the log2 of the Pfaffl (2001) ratio.
- With several references, the mean above is the log2 geometric mean of reference quantities
  (Vandesompele et al. 2002), in the efficiency-corrected form used by qBase (Hellemans et al. 2007).

Statistics are run on `L` per biological sample (Schmittgen and Livak 2008). Efficiency is
`E = 10^(-1/slope)` from Ct against log10 input, and percent efficiency is `(E - 1) * 100`.

Welch's t-test is the default because groups often differ in variance. A `pair` column switches to a
paired t-test, which is the right choice for the same donor, litter, or experimental day across
conditions. Only complete pairs are used.

## Default thresholds

| Check | Default | Source |
| --- | --- | --- |
| Efficiency | 90-110% | Common laboratory practice, not a MIQE limit |
| Curve R^2 | at least 0.98 | Common laboratory practice |
| Dilution levels | at least 5 | Common laboratory practice |
| Technical SD | at most 0.5 cycles | Common laboratory practice; many labs use 0.3 |
| Ct ceiling | 35 | Common laboratory practice; depends on the assay |
| NTC margin | 5 cycles | Working default of this skill |
| Between-group reference shift | 0.5 cycles | Working default of this skill |
| geNorm M | below 1.5 | Vandesompele et al. 2002 |
| Pairwise variation V | below 0.15 | Vandesompele et al. 2002, as a guide |
| Candidate shift between groups | 0.5 cycles | Working default of this skill |

Set a stricter value when the assay's own validation supports it, and record it before analysis.

## Export mapping

Export the results table, then rename:

| This skill | QuantStudio / Design and Analysis | Bio-Rad CFX Maestro | LightCycler |
| --- | --- | --- | --- |
| `sample` | Sample Name | Sample | Sample Name |
| `target` | Target Name | Target | Target name, or one file per assay |
| `ct` | CT or Cq | Cq | Cp or Cq |
| `task` | Task (`UNKNOWN`, `NTC`) | Content (`Unkn`, `NTC`) | Map by hand |

Column names change between software versions, so read the actual header. The `group` column usually
comes from the experimenter's sample sheet, and must be joined on sample name. `Unkn` needs mapping to
`unknown`. Skip the instrument metadata rows above the header.

## MIQE items to report

Report at least: RNA quality and quantity, reverse transcription conditions, primer sequences or assay
ids, amplicon length, efficiency with its curve, NTC and no-RT results, reference genes and how they were
validated, number of biological and technical replicates, the software, and the statistical method
(Bustin et al. 2009).

## Sources

- Livak KJ, Schmittgen TD. Analysis of relative gene expression data using real-time quantitative PCR and
  the 2^-ddCT method. Methods. 2001;25:402-408. https://doi.org/10.1006/meth.2001.1262
- Pfaffl MW. A new mathematical model for relative quantification in real-time RT-PCR. Nucleic Acids Res.
  2001;29:e45. https://doi.org/10.1093/nar/29.9.e45
- Schmittgen TD, Livak KJ. Analyzing real-time PCR data by the comparative CT method. Nat Protoc.
  2008;3:1101-1108. https://doi.org/10.1038/nprot.2008.73
- Bustin SA, et al. The MIQE guidelines. Clin Chem. 2009;55:611-622.
  https://doi.org/10.1373/clinchem.2008.112797
- Vandesompele J, et al. Accurate normalization of real-time quantitative RT-PCR data by geometric
  averaging of multiple internal control genes. Genome Biol. 2002;3:research0034.
  https://doi.org/10.1186/gb-2002-3-7-research0034
- Hellemans J, et al. qBase relative quantification framework and software for management and automated
  analysis of real-time quantitative PCR data. Genome Biol. 2007;8:R19.
  https://doi.org/10.1186/gb-2007-8-2-r19
