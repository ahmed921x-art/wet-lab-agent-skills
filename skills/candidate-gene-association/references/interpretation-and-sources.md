# Statistical notes and sources

## Formulas

Counts `r0, r1, r2` in cases and `s0, s1, s2` in controls are ordered by dose of the effect allele, with
`n_i = r_i + s_i`, `R` cases, and `N` people.

**Trend test** (Sasieni 1997), 1 df:

```text
chi2 = N * [N(r1 + 2r2) - R(n1 + 2n2)]^2 / ( R(N - R) * [N(n1 + 4n2) - (n1 + 2n2)^2] )
```

This equals `N * r^2`, where `r` is the correlation between dose and status. The per-allele OR comes from a
logistic regression of status on dose, fitted by Newton steps. When genotypes separate perfectly, the fit
has no finite solution and the OR is returned as null.

**Allelic test.** A 2x2 table of allele counts. It treats 2N chromosomes as independent, which holds only
under Hardy-Weinberg in the pooled sample. Sasieni (1997) shows that the variance is wrong otherwise.

**Exact Hardy-Weinberg test** (Wigginton et al. 2005). Conditional on allele counts, the probability of `h`
heterozygotes among `n` people, with allele counts `n_A` and `n_B`, is:

```text
P(h) = 2^h * n! * n_A! * n_B! / ( ((n_A - h)/2)! * h! * ((n_B - h)/2)! * (2n)! )
```

The two-sided p sums every `P(h)` no larger than the observed one.

**Inbreeding coefficient.** `F = 1 - H_observed / H_expected`. It is positive for a heterozygote deficit.

**Odds ratio intervals.** Woolf's logit method. Any zero cell adds 0.5 to every cell of that table, and the
output marks it.

## Choosing the test

- A 2x2 table with a smallest expected cell under 5 uses Fisher's exact test. Otherwise Pearson chi-square,
  with no continuity correction. `fisher_p` is reported either way.
- A genotype class with under five people makes the trend and 2 df p-values approximate, and a warning says
  so.

## Hardy-Weinberg threshold

The default warning line is p < 0.05 in controls, which suits a study of a few variants, where each deserves
a look. Genome-wide studies use far lower thresholds because of how many tests they run. That convention
does not transfer here. The test has little power in small samples, so a non-significant result does not
prove the genotyping was correct, and F and direction are worth reading even when p is large.

## Power

`candidate_power.py` compares effect-allele frequency between groups, with `p1` derived from the control
frequency and per-allele OR. Null and alternative standard errors differ, and only the upper rejection
tail is counted, which is negligible away from OR 1. Sham and Purcell (2014) review these ideas. Allele
count variance is multiplied by `1 + F` when `--inbreeding-f` is given.

## Meta-analysis

Study weights are `1 / variance` of the log OR for fixed effect. DerSimonian and Laird (1986) add `tau^2`,
estimated from Q, to each variance. `I^2 = 100 * (Q - df) / Q`, floored at 0 (Higgins et al. 2003).
A published interval is converted with `var = (ln(high) - ln(low))^2 / (2 * 1.96)^2`. Double-zero studies
are refused because they hold no information on the OR.

## Sources

- Wigginton JE, Cutler DJ, Abecasis GR. A note on exact tests of Hardy-Weinberg equilibrium. Am J Hum Genet.
  2005;76:887-893. https://doi.org/10.1086/429864
- Sasieni PD. From genotypes to genes: doubling the sample size. Biometrics. 1997;53:1253-1261.
  https://doi.org/10.2307/2533494
- Clarke GM, et al. Basic statistical analysis in genetic case-control studies. Nat Protoc. 2011;6:121-133.
  https://doi.org/10.1038/nprot.2010.182
- Little J, et al. STrengthening the REporting of Genetic Association Studies (STREGA). PLoS Med.
  2009;6:e1000022. https://doi.org/10.1371/journal.pmed.1000022
- DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7:177-188.
  https://doi.org/10.1016/0197-2456(86)90046-2
- Higgins JPT, et al. Measuring inconsistency in meta-analyses. BMJ. 2003;327:557-560.
  https://doi.org/10.1136/bmj.327.7414.557
- Sham PC, Purcell SM. Statistical power and significance testing in large-scale genetic studies. Nat Rev
  Genet. 2014;15:335-346. https://doi.org/10.1038/nrg3706
- Rigat B, et al. An insertion/deletion polymorphism in the angiotensin I-converting enzyme gene accounting
  for half the variance of serum enzyme levels. J Clin Invest. 1990;86:1343-1346.
  https://doi.org/10.1172/JCI114844
- Shanmugam V, Sell KW, Saha BK. Mistyping ACE heterozygotes. PCR Methods Appl. 1993;3:120-121.
  https://doi.org/10.1101/gr.3.2.120
