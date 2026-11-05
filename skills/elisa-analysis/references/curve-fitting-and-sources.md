# Curve fitting notes and sources

## Fit details

- Parameters `c`, `b`, and `g` are fitted on the log scale, so they stay positive. `b > 0` means `a` is
  always the zero-dose asymptote, so `a > d` identifies a competitive assay.
- Fitting uses `scipy.optimize.least_squares` from a grid of starting midpoints and slopes, and keeps the
  best converged fit. This avoids the local minima a single start can land in.
- Weighted residuals are `(fit - y) / sigma`, with `sigma = y` for `1/y2` and `sigma = sqrt(y)` for `1/y`.
  `sigma` is floored at 1% of the largest signal, so a blank-corrected well near zero cannot dominate.
- Back-calculation inverts the model:
  `x = c * (((a - d) / (y - d))^(1/g) - 1)^(1/b)`. It has no real solution outside the asymptotes, which is
  why a saturated sample gets a status and no number.
- The mean net OD of a sample's wells is back-calculated, then multiplied by the dilution factor.

## Concentration interval

The interval is built with the delta method on log concentration. Parameter covariance comes from the
weighted fit as `s^2 (J'J)^-1`. Read-out variance comes from the same residual scale divided by the number
of wells, so a plate's precision is pooled across all its wells rather than taken from a noisy two-well SD.
The critical value is Student's t on the fit's residual degrees of freedom.

## Level requirements

| Model | Refused below | Warned below |
| --- | --- | --- |
| 4PL | 4 levels | 5 levels |
| 5PL | 6 levels | 7 levels |

Fitting four parameters to four levels leaves nothing to check the fit against, so a run at the minimum
has little evidence behind it. Kits commonly supply seven or eight levels.

## Acceptance values

| Check | Default | Basis |
| --- | --- | --- |
| Standard recovery | within 20% | ICH M10 ligand binding calibrators: 20%, and 25% at LLOQ and ULOQ |
| Replicate CV | up to 20% | ICH M10 ligand binding precision: 20%, and 25% at the limits |
| QC recovery | within 20% | ICH M10 ligand binding routine run QCs |
| Dilution linearity CV | up to 20% | Working default of this skill |

These numbers are used as sensible research defaults. They do not make a run compliant with M10, which
also sets run-level rules and full validation requirements. See `analytical-method-validation`.

## Sources

- Findlay JWA, Dillard RF. Appropriate calibration curve fitting in ligand binding assays. AAPS J.
  2007;9:E260-E267. https://doi.org/10.1208/aapsj0902029
- Gottschalk PG, Dunn JR. The five-parameter logistic: a characterization and comparison with the
  four-parameter logistic. Anal Biochem. 2005;343:54-65. https://doi.org/10.1016/j.ab.2005.04.035
- International Council for Harmonisation. M10 Bioanalytical Method Validation and Study Sample Analysis.
  2022. https://database.ich.org/sites/default/files/M10_Guideline_Step4_2022_0524.pdf
