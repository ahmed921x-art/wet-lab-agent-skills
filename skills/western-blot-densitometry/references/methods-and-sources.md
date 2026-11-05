# Methods, thresholds, and sources

## Calculation

For band `b` of sample `s` on membrane `m`:

```text
net = signal - background
R(s, m) = log2(net_target / net_loading)
L(s, m) = R(s, m) - mean of R over control lanes on membrane m
L(s)    = mean of L(s, m) across membranes
```

Groups are compared on `L(s)` with Welch's t-test, or a paired t-test when `pair` is given. Fold changes and
intervals are back-transformed from log2 and adjusted with Benjamini-Hochberg across all contrasts.

Scaling by exposure cancels within a membrane, so a membrane imaged longer gives the same `L` if it stayed in
the linear range.

## Default thresholds

| Check | Default | Basis |
| --- | --- | --- |
| Log-log slope | 0.8 to 1.2 | Working default of this skill |
| Log-log R^2 | at least 0.95 | Working default of this skill |
| Loading control CV within a membrane | up to 30% | Working default of this skill |
| Between-group loading shift | up to 0.5 log2 units | Working default of this skill |

These are screening values rather than published limits. Set them from your own validation where it exists.

## Housekeeping proteins or total protein

Abundant housekeeping proteins can saturate at loads that suit lower-abundance targets, and some change with
treatment or tissue. Total-protein stains avoid both problems. Pillai-Kastoori et al. (2020) and Taylor and
Posch (2014) discuss the choice. Whichever is used, it needs the same linearity check as the target.

## Sources

- Gassmann M, et al. Quantifying Western blots: pitfalls of densitometry. Electrophoresis. 2009;30:1845-1855.
  https://doi.org/10.1002/elps.200800720
- Taylor SC, Posch A. The design of a quantitative Western blot experiment. Biomed Res Int. 2014;2014:361590.
  https://doi.org/10.1155/2014/361590
- Pillai-Kastoori L, et al. A systematic approach to quantitative Western blot analysis. Anal Biochem.
  2020;593:113608. https://doi.org/10.1016/j.ab.2020.113608
- Janes KA. An analysis of critical factors for quantitative immunoblotting. Sci Signal. 2015;8:rs2.
  https://doi.org/10.1126/scisignal.2005966
