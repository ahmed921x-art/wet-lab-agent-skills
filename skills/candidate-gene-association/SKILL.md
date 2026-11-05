---
name: candidate-gene-association
description: Analyzes small case-control genotype studies of biallelic candidate variants such as SNPs and insertion-deletion polymorphisms. Runs Hardy-Weinberg tests, allele and genotype frequencies, Cochran-Armitage trend, allelic, dominant, recessive, and codominant models, odds ratios, Fisher exact tests, and multiple-testing correction, with QC warnings for genotyping error, inbreeding, and differential missingness. Also computes study power and sample size, and meta-analyzes published odds ratios with forest plots. Use for PCR-RFLP or TaqMan genotype tables, polymorphism and disease association, Hardy-Weinberg equilibrium, or checking genotype counts from a paper. Not for genome-wide data.
license: MIT
compatibility: Requires Python 3.11+. Analysis uses the standard library; figures need matplotlib. Tested on Python 3.13 with matplotlib 3.11.2. No network access.
allowed-tools: Read Write Edit Bash
metadata:
  version: "1.1"
  skill-author: ahmed921x-art
  last-reviewed: "2026-11-06"
---

# Candidate Gene Association

## When to use

Use this skill for a handful of pre-chosen biallelic variants genotyped in cases and controls, whether the
input is a table of people or three genotype counts from a paper. Genome-wide data needs PLINK or a
dedicated pipeline, because it needs ancestry control this skill cannot supply. Variant lookups belong to
`database-lookup`.

Variants with more than two alleles, such as APOE e2/e3/e4, are skipped with a reason. Recode them into
biallelic contrasts first, for a stated hypothesis.

## The three rules

**1. Choose the genetic model before seeing the result.** Testing several models and reporting the best one
is a multiple-testing problem in disguise. One `--model` is primary and carries the correction. The others
are returned as exploratory and should be reported that way. With no biological reason to prefer another,
use the default additive trend test, which stays valid when Hardy-Weinberg does not hold.

**2. Hardy-Weinberg is a genotyping check, read in controls.** A deviation in cases can be true association,
since risk genotypes are enriched there. A deviation in controls means look at the assay before believing
any result. The scripts warn and do not remove anything, because an analyst has to find the cause.

**3. One small study is a hypothesis, not a finding.** Candidate-gene literature has a poor replication
record, largely from small samples and uncontrolled ancestry. Report the interval, not only the p-value,
and treat a nominal result as something to replicate.

## Workflow

### 1. Plan the sample size

Do this before genotyping, or use it to say what a finished study could have detected:

```bash
python scripts/candidate_power.py --frequency 0.3 --cases 150 --controls 150 --odds-ratio 1.5
python scripts/candidate_power.py --frequency 0.3 --cases 150 --controls 150
python scripts/candidate_power.py --frequency 0.3 --odds-ratio 1.5 --inbreeding-f 0.03
```

These give power, the smallest detectable OR, and the people needed per group. Set `--alpha` to the
corrected threshold when testing several variants. Use a realistic OR from the literature, remembering
that first reports tend to overstate effects. `--inbreeding-f` allows for consanguinity, which reduces
effective allele counts. A null result from an underpowered study is not evidence of no effect.

### 2. Check Hardy-Weinberg

```bash
python scripts/hwe_test.py --counts 22,40,17
```

Give the two homozygote counts around the heterozygotes, in controls. Read the exact p whenever an expected
count is below 5. `F` and the direction tell you what kind of deviation you have:

| Pattern in controls | Likely causes | What to do |
| --- | --- | --- |
| Heterozygote deficit, one variant | Allele dropout, preferential amplification, mistyping | Re-genotype a subset by another method |
| Heterozygote deficit, every variant | Consanguinity, population substructure | Report F; consider ancestry-matched controls |
| Heterozygote excess | Paralogous sequence, merged clusters | Check primer specificity and cluster plots |

The ACE I/D polymorphism is a known example of the first row: the shorter D allele can amplify
preferentially and make ID heterozygotes look like DD (Shanmugam et al. 1993). A confirmation PCR with an
insertion-specific primer is the usual fix.

In populations with high consanguinity, a modest heterozygote deficit across all variants is expected. That
is not genotyping error, but it still reduces the number of independent chromosomes in the sample, so an
allelic test overstates its evidence. The default trend test handles this better.

### 3. Run the association

One row per person, one column per variant:

```bash
python scripts/candidate_association.py \
  --input assets/example_genotypes.csv \
  --id-column id \
  --effect-allele ACE_ID=D \
  --model additive \
  --output-json result.json --dosage-csv dosage.csv
```

- Accepted genotype forms are `A/G`, `A|G`, `AG`, `I/D`, `ID`, and `X+/X-`. Missing values include `NA`,
  blank, `0/0`, and `./.`.
- Without `--effect-allele`, the minor allele in controls is used. Name the risk allele from the hypothesis
  so odds ratios point the way the report does.
- `--status-column`, `--case-label`, and `--control-label` handle other coding. Any other status stops the
  run instead of silently dropping people.

To re-check counts from a paper, give them by dose of the effect allele:

```bash
python scripts/candidate_association.py --cases 12,34,32 --controls 22,40,17 \
  --variant-name ACE --effect-name D --other-name I
```

### 4. Read the result

| Model | Odds ratio compares | Notes |
| --- | --- | --- |
| `additive` | Per effect allele | Trend test p; OR from a logistic fit. Robust to HWE deviation |
| `allelic` | Effect vs other allele | Counts chromosomes, so it needs HWE in the pooled sample |
| `dominant` | One or two copies vs none | Fisher exact when an expected cell is below 5 |
| `recessive` | Two copies vs fewer | Often sparse |
| `codominant` | Each genotype vs reference, 2 df | Gives two ORs |

`p_bonferroni` multiplies the primary p by the number of variants tested in the run. If more variants were
tested elsewhere, pass the honest total. In `--counts` mode use `--n-tests`.

Warnings to act on: Hardy-Weinberg deviation in controls, missing calls that differ between cases and
controls (technical artefact can produce false association), and zero cells that forced a 0.5 correction
(the OR is unstable).

### 5. Adjust for covariates

These tests are crude. Age, sex, and ancestry can confound them. Use the `--dosage-csv` output with
logistic regression, here using the `statsmodels` skill. Merge covariates by `id` first. This snippet is
illustrative:

```python
import pandas as pd
import statsmodels.formula.api as smf

data = pd.read_csv("dosage.csv").merge(pd.read_csv("clinical.csv"), on="id")
data["case"] = (data["status"] == "case").astype(int)
fit = smf.logit("case ~ ACE_ID + age + C(sex)", data=data).fit()
print(fit.summary())
```

Keep about ten events per parameter, or the model becomes unstable in a small study.

### 6. Meta-analyze the literature

To place a result among published studies, or to pool them, build one row per study with `a, b, c, d`
counts or `or, ci_low, ci_high`. Add `hwe_p` when known. Keep one variant and one genetic model per file.

```bash
python scripts/candidate_meta.py --input assets/example_studies.csv --output-json meta.json
```

This reports fixed-effect and DerSimonian-Laird random-effects estimates, Q, I^2, tau^2, leave-one-out
results, and a sensitivity analysis without studies whose controls deviate from Hardy-Weinberg. High I^2
means asking why studies differ, often by ancestry, before quoting one number. With few studies, treat
heterogeneity estimates as rough. For a systematic review, use the `literature-review` skill to find
studies and report by PRISMA.

### 7. Report

```bash
python scripts/candidate_report.py --result result.json --meta meta.json --out-dir report
```

This writes a forest plot and `report.md`, with warnings first. Either input can be given alone. Then
follow STREGA (Little et al. 2009): how cases and controls were selected and matched for ancestry, the
genotyping method and its quality checks, call rate, Hardy-Weinberg results in controls, the pre-specified
model, ORs with intervals, and every variant tested, including the null ones.

## Limits

- No correction for population structure is possible from a few variants. Ancestry-matched controls are the
  defence, and the report should describe them.
- Bonferroni is conservative for variants in linkage disequilibrium, and this skill does not model LD or
  haplotypes.
- An odds ratio from a case-control study is not a risk, and does not by itself imply causation.

## Sources

Citations and further notes are in `references/interpretation-and-sources.md`.
