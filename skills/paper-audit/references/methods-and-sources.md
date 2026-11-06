# Methods, rounding rules, and sources

## Rounding intervals

A number printed with `d` decimals stands for every value within half a unit of its last digit. For a printed
value `v` the interval is `[v - 0.5 x 10^-d, v + 0.5 x 10^-d]`. Both ends are included, so a value exactly
halfway (for example 0.125 printed as 0.12 or 0.13) is accepted under either rounding convention. Numbers are
parsed with Python's `Decimal`, so 0.005 is exact and not a binary approximation.

Whole numbers are treated as follows:

- degrees of freedom, `n`, `count`, and `items` are exact;
- a test statistic or mean printed without decimals still has a rounding interval of plus or minus 0.5;
- in a `sum`, whole numbers are exact counts (or 100%), and decimals are rounded values.

Truncation, where 0.0349 would be printed as .03, is not allowed for. If a journal truncates, say so in the
report and treat borderline flags with care.

## p-values from test statistics

The p-value is recomputed at the corners of the rounding intervals of the statistic and, when printed with
decimals, of the df. Because p moves monotonically with the statistic and the df, the corners bound every
value in between. The reported p is consistent when its own rounding interval (or inequality) overlaps the
recomputed range. This follows the logic of statcheck (Nuijten et al. 2016) and Bakker and Wicherts (2011),
reimplemented here without their code.

| Test | p-value (two-tailed unless noted) |
| --- | --- |
| t with df | `I_x(df/2, 1/2)`, with `x = df / (df + t^2)` |
| z | `erfc(|z| / sqrt 2)` |
| r with df = n - 2 | t test of `t = r sqrt(df) / sqrt(1 - r^2)` |
| F with df1, df2 | upper tail `I_x(df2/2, df1/2)`, with `x = df2 / (df2 + df1 F)` |
| chi-square with df | upper tail `Q(df/2, chi2/2)` |

`I` is the regularized incomplete beta function and `Q` the regularized upper incomplete gamma function. Both
are evaluated with continued fractions using the modified Lentz method (Lentz 1976; Thompson and Barnett
1986). For the beta function, the symmetry `I_x(a, b) = 1 - I_(1-x)(b, a)` keeps the fraction in its fast
region. DiDonato and Morris (1992) describe the accuracy issues. The test suite compares these functions with
scipy across a grid that includes small and fractional parameters, and requires a relative error below 1e-9.

A one-tailed p is half the two-tailed value for t, z, and r. F and chi-square p-values are upper-tail by
construction. When a two-tailed check fails but a one-tailed check would pass, the ledger adds a note and
keeps the flag, since only the paper can say which was intended.

A **decision error** is a flag where the printed p and the recomputed range fall on opposite sides of alpha.
When the printed p equals alpha exactly, no decision is assigned.

## GRIM

For a mean of integer responses from `n` people answering `k` items, every possible mean is a multiple of
`1 / (n k)`. The printed mean is possible when some integer total `S` gives `S / (n k)` inside its rounding
interval (Brown and Heathers 2017). The check finds the smallest and largest integers in
`[(m - h) n k, (m + h) n k]` with exact arithmetic. When `1 / (n k)` is no larger than the interval width, every
mean is possible and the claim is not checkable. For two decimals that happens once `n k` reaches 100.

The tests compare this with a brute-force search over all totals, written as an independent function using
`fractions.Fraction`.

## Percentages

With a count, `100 count / n` must fall in the printed percentage's interval. Without a count, the check asks
whether any integer count from 0 to `n` would do, which is GRIM applied to a proportion. When `100 / n` is no
larger than the interval width, every percentage is possible and the claim is not checkable.

## Parts and totals

The parts must sum to the total within the summed half-units of the rounded values. Integer group sizes must
match exactly. Percentages of 33.3, 33.3, and 33.3 against 100 are allowed a slack of 0.15.

## Estimates with confidence intervals

For ratios (odds, risk, and hazard ratios), a Wald interval is symmetric on the log scale, so the estimate
should equal the geometric mean of the bounds, `sqrt(L U)`. For differences it should equal the midpoint
`(L + U) / 2`. The centre is recomputed across the rounding intervals of both bounds and compared with the
printed estimate's interval. Georgescu and Wren (2018) applied this approach to published ratios at scale.

Agreement between the interval and p is checked as a logical rule. A `level`% interval that excludes the null
value (1 for ratios, 0 for differences) corresponds to p below `1 - level/100`. A disagreement is flagged. When
a bound equals the null within rounding, no judgement is made.

The ledger also reports an approximate p from the interval (Altman and Bland 2011a):
`SE = (T(U) - T(L)) / (2 z)`, with `T` the log for ratios and the identity for differences, and `z` the normal
quantile, or a t quantile when `df` is given. The p then follows from `T(estimate) / SE`. Because the printed p
may come from a different test, a mismatch here is noted rather than flagged. Altman and Bland (2011b) give the
reverse calculation.

Intervals from exact, profile-likelihood, score, or bootstrap methods are asymmetric by design. Set `ci_method`
and only containment and the significance rule are checked.

## How the tests were built

- Distribution functions: compared with `scipy.special` and `scipy.stats`, including the exact null
  distribution of Pearson's r.
- GRIM and percentages: compared with independent brute-force enumeration in exact fractions.
- Rounding, intervals, and sums: hand calculations, written out in the test comments.
- Bug injection: 18 deliberate faults, injected one at a time (for example swapped F degrees of freedom, a
  rounding interval ten times too narrow, an inverted significance rule, a dropped misconduct disclaimer). The
  suite failed for each one. One proposed fault, removing a special case for `p = .000`, changed nothing,
  because the rounding interval already handles it. That special case was then removed as redundant.

## Sources

All DOIs were checked against Crossref on 2026-11-06.

- Brown NJL, Heathers JAJ. The GRIM test: a simple technique detects numerous anomalies in the reporting of
  results in psychology. Soc Psychol Personal Sci. 2017;8(4):363-369. https://doi.org/10.1177/1948550616673876
- Nuijten MB, Hartgerink CHJ, van Assen MALM, Epskamp S, Wicherts JM. The prevalence of statistical reporting
  errors in psychology (1985-2013). Behav Res Methods. 2016;48(4):1205-1226.
  https://doi.org/10.3758/s13428-015-0664-2
- Bakker M, Wicherts JM. The (mis)reporting of statistical results in psychology journals. Behav Res Methods.
  2011;43(3):666-678. https://doi.org/10.3758/s13428-011-0089-5
- Georgescu C, Wren JD. Algorithmic identification of discrepancies between published ratios and their
  reported confidence intervals and P-values. Bioinformatics. 2018;34(10):1758-1766.
  https://doi.org/10.1093/bioinformatics/btx811
- Altman DG, Bland JM. How to obtain the P value from a confidence interval. BMJ. 2011;343:d2304.
  https://doi.org/10.1136/bmj.d2304 (2011a)
- Altman DG, Bland JM. How to obtain the confidence interval from a P value. BMJ. 2011;343:d2090.
  https://doi.org/10.1136/bmj.d2090 (2011b)
- DiDonato AR, Morris AH. Algorithm 708: significant digit computation of the incomplete beta function
  ratios. ACM Trans Math Softw. 1992;18(3):360-373. https://doi.org/10.1145/131766.131776
- Lentz WJ. Generating Bessel functions in Mie scattering calculations using continued fractions. Appl Opt.
  1976;15(3):668. https://doi.org/10.1364/AO.15.000668
- Thompson IJ, Barnett AR. Coulomb and Bessel functions of complex arguments and order. J Comput Phys.
  1986;64(2):490-509. https://doi.org/10.1016/0021-9991(86)90046-X
