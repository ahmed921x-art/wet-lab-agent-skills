"""Exact and asymptotic statistics for biallelic genotype counts (standard library)."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Any

MISSING_TOKENS = {"", "na", "n/a", "nan", "none", "null", "-", "--", "0", "00", "0/0", "./.", ".", "nc", "nn", "fail", "failed"}


class InputError(ValueError):
    """Raised for input problems the analyst has to resolve."""


def read_text(path: str | Path) -> str:
    """Read UTF-8 or UTF-16 text, and explain any other encoding."""

    source = Path(path)
    if not source.is_file():
        raise InputError(f"input file not found: {source}")
    raw = source.read_bytes()
    try:
        if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
            text = raw.decode("utf-16")
        else:
            text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise InputError(f"{source} is not UTF-8; re-save it as CSV UTF-8") from exc
    if not text.strip():
        raise InputError(f"input file is empty: {source}")
    return text


def safe_console() -> None:
    """Keep non-Latin sample names from crashing output on legacy Windows code pages."""

    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")


def chi2_sf(statistic: float, df: int) -> float:
    """Upper tail of the chi-square distribution for 1 or 2 degrees of freedom."""

    if df == 1:
        return math.erfc(math.sqrt(statistic / 2.0))
    if df == 2:
        return math.exp(-statistic / 2.0)
    raise ValueError("only df 1 and 2 are supported")


def parse_genotype(raw: str) -> tuple[str, str] | None:
    """Split a genotype into two alleles, or None when it is missing.

    Accepts A/G, A|G, I/D, X+/X-, and unseparated even-length calls (AG, II, ID).
    Allele order is normalized so that AG and GA count as one genotype.
    """

    text = str(raw).strip()
    if text.lower() in MISSING_TOKENS:
        return None
    for separator in ("/", "|"):
        if separator in text:
            parts = [part.strip() for part in text.split(separator)]
            if len(parts) != 2 or not all(parts):
                raise InputError(f"cannot parse genotype {raw!r}")
            return tuple(sorted(parts))  # type: ignore[return-value]
    if len(text) % 2 == 0:
        half = len(text) // 2
        return tuple(sorted((text[:half], text[half:])))  # type: ignore[return-value]
    raise InputError(f"cannot parse genotype {raw!r}")


def hwe_exact_p(hom_a: int, het: int, hom_b: int) -> float:
    """Two-sided exact Hardy-Weinberg p-value (Wigginton et al. 2005 definition).

    Sums the probability of every heterozygote count, given the allele counts,
    that is no more likely than the one observed.
    """

    total = hom_a + het + hom_b
    if total == 0:
        raise InputError("no genotypes to test")
    minor = min(2 * hom_a + het, 2 * hom_b + het)
    major = 2 * total - minor

    def log_probability(hets: int) -> float:
        rare_homs = (minor - hets) // 2
        common_homs = total - hets - rare_homs
        return (
            hets * math.log(2.0)
            + math.lgamma(total + 1)
            - math.lgamma(rare_homs + 1)
            - math.lgamma(hets + 1)
            - math.lgamma(common_homs + 1)
            + math.lgamma(minor + 1)
            + math.lgamma(major + 1)
            - math.lgamma(2 * total + 1)
        )

    supports = range(minor % 2, minor + 1, 2)
    probabilities = {hets: math.exp(log_probability(hets)) for hets in supports}
    normalizer = sum(probabilities.values())
    observed = probabilities[het]
    tail = sum(p for p in probabilities.values() if p <= observed * (1 + 1e-7))
    return min(1.0, tail / normalizer)


def hardy_weinberg(hom_a: int, het: int, hom_b: int) -> dict[str, Any]:
    """Allele frequency, expected counts, chi-square, exact p, and inbreeding F."""

    total = hom_a + het + hom_b
    if total == 0:
        raise InputError("no genotypes to test")
    p = (2 * hom_a + het) / (2 * total)
    q = 1.0 - p
    expected = (p * p * total, 2 * p * q * total, q * q * total)
    monomorphic = p in (0.0, 1.0)
    chi_square = (
        None
        if monomorphic
        else sum((o - e) ** 2 / e for o, e in zip((hom_a, het, hom_b), expected))
    )
    expected_het_fraction = 2 * p * q
    inbreeding = None if monomorphic else 1.0 - (het / total) / expected_het_fraction
    direction = "none"
    if not monomorphic:
        if het < expected[1]:
            direction = "heterozygote_deficit"
        elif het > expected[1]:
            direction = "heterozygote_excess"
    return {
        "n": total,
        "observed": [hom_a, het, hom_b],
        "expected": [round(value, 2) for value in expected],
        "allele_frequency_first": round(p, 4),
        "minor_allele_frequency": round(min(p, q), 4),
        "chi_square": None if chi_square is None else round(chi_square, 4),
        "chi_square_p": None if chi_square is None else chi2_sf(chi_square, 1),
        "exact_p": hwe_exact_p(hom_a, het, hom_b),
        "inbreeding_f": None if inbreeding is None else round(inbreeding, 4),
        "direction": direction,
        "min_expected_count": round(min(expected), 2),
    }


def fisher_exact(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p for [[a, b], [c, d]], by point probability."""

    row1, row2, col1 = a + b, c + d, a + c
    total = row1 + row2

    def log_choose(n: int, k: int) -> float:
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    def probability(cell: int) -> float:
        return math.exp(
            log_choose(row1, cell) + log_choose(row2, col1 - cell) - log_choose(total, col1)
        )

    low, high = max(0, col1 - row2), min(col1, row1)
    observed = probability(a)
    tail = sum(
        value
        for value in (probability(cell) for cell in range(low, high + 1))
        if value <= observed * (1 + 1e-7)
    )
    return min(1.0, tail)


def odds_ratio(a: int, b: int, c: int, d: int) -> dict[str, Any]:
    """Odds ratio with a Woolf interval; adds 0.5 to all cells when one is zero."""

    corrected = 0 in (a, b, c, d)
    pad = 0.5 if corrected else 0.0
    w, x, y, z = a + pad, b + pad, c + pad, d + pad
    if min(w, x, y, z) == 0:
        return {"or": None, "ci_low": None, "ci_high": None, "haldane_corrected": corrected}
    estimate = (w * z) / (x * y)
    se_log = math.sqrt(1 / w + 1 / x + 1 / y + 1 / z)
    return {
        "or": estimate,
        "ci_low": math.exp(math.log(estimate) - 1.959963985 * se_log),
        "ci_high": math.exp(math.log(estimate) + 1.959963985 * se_log),
        "haldane_corrected": corrected,
    }


def pearson_2x2(a: int, b: int, c: int, d: int) -> float | None:
    """Pearson chi-square without continuity correction."""

    total = a + b + c + d
    margins = (a + b) * (c + d) * (a + c) * (b + d)
    return None if margins == 0 else total * (a * d - b * c) ** 2 / margins


def trend_chi_square(cases: tuple[int, int, int], controls: tuple[int, int, int]) -> float | None:
    """Cochran-Armitage trend chi-square with scores 0, 1, 2 (Sasieni 1997 form)."""

    r0, r1, r2 = cases
    s0, s1, s2 = controls
    n0, n1, n2 = r0 + s0, r1 + s1, r2 + s2
    big_r = r0 + r1 + r2
    total = n0 + n1 + n2
    denominator = big_r * (total - big_r) * (total * (n1 + 4 * n2) - (n1 + 2 * n2) ** 2)
    if denominator == 0:
        return None
    numerator = total * (total * (r1 + 2 * r2) - big_r * (n1 + 2 * n2)) ** 2
    return numerator / denominator


def pearson_2x3(cases: tuple[int, int, int], controls: tuple[int, int, int]) -> float | None:
    """Pearson chi-square for genotype by status; None when a column is empty."""

    big_r, big_s = sum(cases), sum(controls)
    total = big_r + big_s
    statistic = 0.0
    for r, s in zip(cases, controls):
        column = r + s
        if column == 0:
            return None
        expected_r = big_r * column / total
        expected_s = big_s * column / total
        statistic += (r - expected_r) ** 2 / expected_r + (s - expected_s) ** 2 / expected_s
    return statistic


def per_allele_logistic(
    cases: tuple[int, int, int], controls: tuple[int, int, int], iterations: int = 60
) -> dict[str, Any] | None:
    """Maximum-likelihood per-allele OR, by Newton steps on grouped data."""

    alpha = beta = 0.0
    for _ in range(iterations):
        score_a = score_b = info_aa = info_ab = info_bb = 0.0
        for dose, (r, s) in enumerate(zip(cases, controls)):
            n = r + s
            if n == 0:
                continue
            eta = alpha + beta * dose
            if abs(eta) > 30:
                return None
            pi = 1.0 / (1.0 + math.exp(-eta))
            residual = r - n * pi
            variance = n * pi * (1.0 - pi)
            score_a += residual
            score_b += dose * residual
            info_aa += variance
            info_ab += dose * variance
            info_bb += dose * dose * variance
        determinant = info_aa * info_bb - info_ab * info_ab
        if determinant <= 1e-12:
            return None
        step_a = (info_bb * score_a - info_ab * score_b) / determinant
        step_b = (info_aa * score_b - info_ab * score_a) / determinant
        alpha += step_a
        beta += step_b
        if abs(step_a) + abs(step_b) < 1e-10:
            se_beta = math.sqrt(info_aa / determinant)
            return {
                "or": math.exp(beta),
                "ci_low": math.exp(beta - 1.959963985 * se_beta),
                "ci_high": math.exp(beta + 1.959963985 * se_beta),
                "alpha": alpha,
                "beta": beta,
            }
    return None
