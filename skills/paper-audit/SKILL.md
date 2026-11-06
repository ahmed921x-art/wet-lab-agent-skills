---
name: paper-audit
description: Recomputes the numbers reported in a scientific paper and reports which claims are internally consistent. Checks p-values against their t, F, chi-square, z, or r statistics and degrees of freedom; GRIM for means of integer data; percentages against counts and n; subgroups against totals; and odds ratios, risk ratios, or mean differences against their confidence intervals and p-values. Allows for rounding before flagging, and ties every finding to an exact quote and page. Use when checking a manuscript, preprint, or published paper for statistical reporting errors, when reviewing or refereeing, or when asked whether a paper's numbers add up.
license: MIT
compatibility: Requires Python 3.11+ standard library only. Tests use scipy. Tested on Python 3.13. No network access.
allowed-tools: Read Write Edit Bash
metadata:
  version: "1.0"
  skill-author: ahmed921x-art
  last-reviewed: "2026-11-06"
---

# Paper Audit

## When to use

Use this skill to check whether the numbers printed in a paper agree with each other: a test statistic with
its p-value, a mean with its sample size, a percentage with its count, an odds ratio with its interval. It
recomputes what can be recomputed from the printed values alone. It does not read the raw data, judge
whether the right test was used, or assess the science. For those, read the paper as a reviewer would.

## The three rules

**1. An inconsistency is a question, not an accusation.** Typing errors, copy-paste slips, different rounding,
one-tailed or corrected tests, and mistakes in building the claim table all produce flags. Never describe a
flag as misconduct, fraud, or fabrication. Write that the numbers do not agree and say what could explain it.
The report says this at the top, and any summary you write must say it too.

**2. Rounding comes first.** "t = 2.20" means a value between 2.195 and 2.205, and "p = .04" means between
.035 and .045. Every check recomputes across those intervals and flags only when no value inside them fits.
So numbers must be entered exactly as printed, keeping trailing zeros ("2.20", not "2.2").

**3. Every finding traces to the paper.** Each claim carries the exact quote and the page. A row without
both is marked not checkable. Say plainly which claims could not be checked, and that claims not entered
were not audited.

## Workflow

### 1. Extract the claims

Read the paper and copy each checkable number into a claim table, one row per claim. Quote the sentence or
table cell exactly. Use the page number of the PDF or the journal page, and say which in the report.
Enter numbers as printed: a leading dot (`.31`), Unicode minus (`−2.1`), and `p < .001` are all accepted.
Decimal commas are rejected, so convert `2,20` to `2.20` by hand.

| Column | Used by | Meaning |
| --- | --- | --- |
| `id`, `page`, `quote`, `check` | all | Claim id, page, exact text, and one of the checks below |
| `test`, `statistic`, `df`, `df2`, `n`, `tails` | `p_value` | `t`, `F`, `chi2`, `z`, or `r`; F uses `df` and `df2`; r takes `n` or `df`; `tails` 1 or 2 (default 2) |
| `p` | `p_value`, CI checks | As printed: `= .036`, `< .001`, `> .05`, `ns` |
| `mean`, `n`, `items` | `grim` | Mean, participants, and items averaged per participant (default 1) |
| `percent`, `n`, `count` | `percent` | Percentage, its denominator, and the count if given |
| `parts`, `total` | `sum` | Values separated by `;`, and the total they should reach |
| `estimate`, `ci_low`, `ci_high`, `level`, `ci_method`, `df` | `ratio_ci`, `diff_ci` | Estimate, bounds, level (default 95), interval method, and df for t-based intervals |
| `note` | all | Free text carried into the report |

A CSV or a JSON list of objects with the same keys both work. See `assets/example_claims.csv`, which is a
synthetic example built for testing and does not come from any real paper.

What to enter, and what to leave out:

- `p_value`: only when the statistic, its df, and the p-value are all printed. Use the df as printed, including
  decimals from Welch or Greenhouse-Geisser corrections.
- `grim`: only for means of integer responses, such as Likert items or counts, and only when n is the number of
  people behind that mean. It does not apply to continuous measures like reaction time or weight.
- `percent`: use the denominator the paper used for that percentage, which may differ from the total n.
- `sum`: group sizes that should add up to the total, or percentages that should add up to 100.
- `ratio_ci` for odds, risk, or hazard ratios; `diff_ci` for mean or risk differences. Set `ci_method` to
  `exact`, `profile`, `score`, or `bootstrap` when the paper used one, since those intervals are legitimately
  asymmetric.

### 2. Run the checks

```bash
python scripts/audit_claims.py --input assets/example_claims.csv --output-json ledger.json --output-csv ledger.csv
```

Each claim ends with one status:

| Status | Meaning |
| --- | --- |
| consistent | Some value inside the rounding intervals fits every printed number |
| inconsistent | No value inside the rounding intervals fits; the reason says what disagrees |
| not checkable | Information is missing, the input is unreadable, or the check cannot discriminate (for example GRIM with n of 100 or more at two decimals) |

`decision_error` marks a flag where the significance decision itself changes at `--alpha` (default 0.05).

### 3. Write the report

```bash
python scripts/audit_report.py --ledger ledger.json --out audit_report.md
```

The report opens with how to read it, lists inconsistencies first with their page and quote, then claims that
could not be checked and why, then consistent claims, and ends with what the audit cannot detect.

### 4. Review every flag before passing it on

For each inconsistent claim, re-read the quoted passage and check the row against it. Most flags in practice
come from the claim table, not the paper. Then look for an innocent explanation:

| Flag | Explanations to rule out first |
| --- | --- |
| p-value | One-tailed test (the ledger notes when that would fit), corrected df, a different test than entered, a p adjusted for multiple comparisons |
| GRIM | Missing responses for some participants, a multi-item scale entered with `items` = 1, non-integer data |
| percent | A different denominator, such as valid responses rather than all participants |
| sum | Exclusions, dropouts, or categories that overlap |
| interval | An exact, profile, or bootstrap interval; a p from a different test than the interval |

Only then describe what remains, in neutral terms: "the reported p = .04 does not match t(28) = 1.20, which
gives p = .24".

## Limits

- Only printed numbers are checked. A paper can be internally consistent and still wrong.
- The p from an interval (Altman and Bland method) is approximate, so a mismatch with the printed p is a
  note and not a flag. Only a disagreement about significance is flagged.
- Ratios are assumed symmetric on the log scale and differences on the linear scale unless `ci_method`
  says otherwise.
- Whole numbers in a `sum` are treated as exact counts and decimals as rounded values.
- Standard deviations are not checked (no GRIMMER or SPRITE).
- Claim extraction is manual or done by the agent. Check the table against the paper before trusting a run.

## Sources

Methods, rounding rules, numerical accuracy, and citations are in `references/methods-and-sources.md`.
