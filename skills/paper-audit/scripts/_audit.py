"""Shared numerics and checks for paper-audit. Standard library only.

Every reported number is read as text so that its precision is kept: "2.20" means
the true value lies in [2.195, 2.205]. Each check works with those rounding
intervals and only reports an inconsistency when no value inside them fits.
"""

from __future__ import annotations

import math
import re
from decimal import Decimal, InvalidOperation
from statistics import NormalDist
from typing import Any

STANDARD_NORMAL = NormalDist()
ALPHA_DEFAULT = 0.05
CONSISTENT = "consistent"
INCONSISTENT = "inconsistent"
NOT_CHECKABLE = "not checkable"


class InputError(ValueError):
    """Raised for an input problem the analyst has to resolve."""


# ---------------------------------------------------------------------------
# Special functions
# ---------------------------------------------------------------------------

_TINY = 1e-300
_EPS = 1e-15
_MAX_ITER = 10000


def _beta_fraction(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta ratio, evaluated by modified Lentz."""

    c_term = 1.0
    d_term = 1.0 - (a + b) * x / (a + 1.0)
    d_term = 1.0 / (d_term if abs(d_term) > _TINY else _TINY)
    total = d_term
    for m in range(1, _MAX_ITER + 1):
        # even step
        numerator = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        d_term = 1.0 + numerator * d_term
        c_term = 1.0 + numerator / c_term
        d_term = 1.0 / (d_term if abs(d_term) > _TINY else _TINY)
        c_term = c_term if abs(c_term) > _TINY else _TINY
        total *= d_term * c_term
        # odd step
        numerator = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d_term = 1.0 + numerator * d_term
        c_term = 1.0 + numerator / c_term
        d_term = 1.0 / (d_term if abs(d_term) > _TINY else _TINY)
        c_term = c_term if abs(c_term) > _TINY else _TINY
        step = d_term * c_term
        total *= step
        if abs(step - 1.0) < _EPS:
            return total
    raise ArithmeticError("incomplete beta did not converge")


def betainc(a: float, b: float, x: float) -> float:
    """Regularized incomplete beta I_x(a, b) for a, b > 0 and 0 <= x <= 1."""

    if a <= 0 or b <= 0:
        raise ValueError("betainc needs a > 0 and b > 0")
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    log_front = a * math.log(x) + b * math.log1p(-x) - (math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(log_front) * _beta_fraction(a, b, x) / a
    return 1.0 - math.exp(log_front) * _beta_fraction(b, a, 1.0 - x) / b


def gammaincc(a: float, x: float) -> float:
    """Regularized upper incomplete gamma Q(a, x) for a > 0 and x >= 0."""

    if a <= 0:
        raise ValueError("gammaincc needs a > 0")
    if x <= 0:
        return 1.0
    log_front = a * math.log(x) - x - math.lgamma(a)
    if x < a + 1.0:
        term = 1.0 / a
        total = term
        for n in range(1, _MAX_ITER + 1):
            term *= x / (a + n)
            total += term
            if abs(term) < abs(total) * _EPS:
                return max(0.0, 1.0 - total * math.exp(log_front))
        raise ArithmeticError("incomplete gamma series did not converge")
    b_term = x + 1.0 - a
    c_term = 1.0 / _TINY
    d_term = 1.0 / b_term
    total = d_term
    for n in range(1, _MAX_ITER + 1):
        numerator = -n * (n - a)
        b_term += 2.0
        d_term = numerator * d_term + b_term
        d_term = 1.0 / (d_term if abs(d_term) > _TINY else _TINY)
        c_term = b_term + numerator / c_term
        c_term = c_term if abs(c_term) > _TINY else _TINY
        step = d_term * c_term
        total *= step
        if abs(step - 1.0) < _EPS:
            return math.exp(log_front) * total
    raise ArithmeticError("incomplete gamma fraction did not converge")


def p_from_t(t: float, df: float, tails: int = 2) -> float:
    """p-value for Student's t with df degrees of freedom."""

    two_sided = betainc(df / 2.0, 0.5, df / (df + t * t))
    return two_sided if tails == 2 else two_sided / 2.0


def p_from_z(z: float, tails: int = 2) -> float:
    two_sided = math.erfc(abs(z) / math.sqrt(2.0))
    return two_sided if tails == 2 else two_sided / 2.0


def p_from_f(f: float, df1: float, df2: float) -> float:
    """Upper-tail p-value for F(df1, df2)."""

    if f <= 0:
        return 1.0
    return betainc(df2 / 2.0, df1 / 2.0, df2 / (df2 + df1 * f))


def p_from_chi2(chi2: float, df: float) -> float:
    """Upper-tail p-value for chi-square with df degrees of freedom."""

    return gammaincc(df / 2.0, chi2 / 2.0)


def p_from_r(r: float, df: float, tails: int = 2) -> float:
    """p-value for a Pearson correlation via its t transform with df = n - 2."""

    if abs(r) >= 1:
        return 0.0
    return p_from_t(r * math.sqrt(df) / math.sqrt(1.0 - r * r), df, tails)


def t_quantile(probability: float, df: float) -> float:
    """Inverse of the t distribution for 0.5 < probability < 1, by bisection."""

    target = 2.0 * (1.0 - probability)
    low, high = 0.0, 1.0
    while p_from_t(high, df) > target:
        high *= 2.0
    for _ in range(200):
        mid = (low + high) / 2.0
        if p_from_t(mid, df) > target:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


# ---------------------------------------------------------------------------
# Reading reported numbers with their precision
# ---------------------------------------------------------------------------

_MINUS = {"−": "-", "–": "-", "—": "-", "‒": "-"}


class Reported:
    """A reported number and the interval it stands for after rounding."""

    def __init__(self, text: str, value: Decimal, half_unit: Decimal):
        self.text = text
        self.value = value
        self.half_unit = half_unit

    @property
    def low(self) -> float:
        return float(self.value - self.half_unit)

    @property
    def high(self) -> float:
        return float(self.value + self.half_unit)

    @property
    def decimals(self) -> int:
        return max(0, -self.value.as_tuple().exponent)

    def __float__(self) -> float:
        return float(self.value)


def clean_text(raw: Any) -> str:
    text = "" if raw is None else str(raw)
    for bad, good in _MINUS.items():
        text = text.replace(bad, good)
    return text.strip()


def parse_reported(raw: Any, field: str, exact: bool = False) -> Reported:
    """Parse a number as printed in a paper; exact=True treats it as a count."""

    text = clean_text(raw)
    if not text:
        raise InputError(f"{field} is empty")
    if re.fullmatch(r"-?\d+,\d+", text):
        raise InputError(f"{field} '{text}' uses a decimal comma; write it with a dot")
    if not re.fullmatch(r"[-+]?(\d+\.?\d*|\.\d+)", text):
        raise InputError(f"{field} '{text}' is not a plain number")
    try:
        value = Decimal(text if not text.startswith(".") else "0" + text)
    except InvalidOperation as exc:  # pragma: no cover - regex already guards this
        raise InputError(f"{field} '{text}' is not a number") from exc
    exponent = value.as_tuple().exponent
    if exact:
        if exponent < 0 and value != value.to_integral_value():
            raise InputError(f"{field} '{text}' must be a whole number")
        return Reported(text, value, Decimal(0))
    half_unit = Decimal(5) * Decimal(10) ** (exponent - 1)
    return Reported(text, value, half_unit)


def parse_count(raw: Any, field: str, minimum: int = 0) -> int:
    value = int(parse_reported(raw, field, exact=True).value)
    if value < minimum:
        raise InputError(f"{field} must be at least {minimum}, got {value}")
    return value


_P_PATTERN = re.compile(r"^p?\s*(<=|>=|=|<|>|≤|≥)?\s*(.+)$", re.I)


def parse_p(raw: Any) -> tuple[str, Reported]:
    """Parse a reported p-value such as '= .035', '<.001', 'p > 0.05' or 'ns'.

    Returns the relation ('=', '<', '>') and the number.
    """

    text = clean_text(raw).replace(" ", "")
    if not text:
        raise InputError("p is empty")
    if text.lower() in {"ns", "n.s.", "p>.05"}:
        return ">", parse_reported("0.05", "p")
    match = _P_PATTERN.match(text)
    if not match:
        raise InputError(f"p '{raw}' is not readable")
    relation = {"<=": "<", "≤": "<", ">=": ">", "≥": ">", None: "="}.get(match.group(1), match.group(1))
    number = parse_reported(match.group(2), "p")
    if not 0 <= number.value <= 1:
        raise InputError(f"p '{raw}' is outside 0 to 1")
    # A printed "p = .000" needs no special case: its rounding interval is [0, .0005].
    return relation, number


def p_range_fits(relation: str, reported: Reported, low: float, high: float) -> bool:
    """Do the recomputed p-values [low, high] fit the reported p after rounding?"""

    if relation == "<":
        return low < reported.high if reported.half_unit == 0 else low < float(reported.value)
    if relation == ">":
        return high > float(reported.value)
    return high >= reported.low and low <= reported.high


def reported_significant(relation: str, reported: Reported, alpha: float) -> bool | None:
    """Whether the reported p claims significance at alpha; None if it cannot be told."""

    value = float(reported.value)
    if relation == "<":
        return True if value <= alpha else None
    if relation == ">":
        return False if value >= alpha else None
    return value <= alpha if abs(value - alpha) > 1e-12 else None


# ---------------------------------------------------------------------------
# Result helpers
# ---------------------------------------------------------------------------


def result(status: str, reason: str, recomputed: dict[str, Any] | None = None, notes: list[str] | None = None,
           decision_error: bool = False) -> dict[str, Any]:
    return {
        "status": status,
        "reason": reason,
        "recomputed": recomputed or {},
        "notes": notes or [],
        "decision_error": decision_error,
    }


def need(claim: dict[str, Any], *fields: str) -> None:
    missing = [field for field in fields if not clean_text(claim.get(field))]
    if missing:
        raise InputError("missing " + ", ".join(missing))


def sig(value: float) -> float:
    """Round for display without hiding small p-values."""

    return float(f"{value:.4g}")


# ---------------------------------------------------------------------------
# Check 1: p-value from a test statistic
# ---------------------------------------------------------------------------

TESTS = {"t", "f", "chi2", "z", "r"}


def _df_range(claim: dict[str, Any], field: str) -> tuple[float, float, str]:
    reported = parse_reported(claim.get(field), field)
    if reported.value <= 0:
        raise InputError(f"{field} must be positive")
    if reported.decimals == 0:
        return float(reported.value), float(reported.value), reported.text
    return max(reported.low, 1e-9), reported.high, reported.text


def check_p_value(claim: dict[str, Any], alpha: float = ALPHA_DEFAULT) -> dict[str, Any]:
    need(claim, "test", "statistic", "p")
    test = clean_text(claim["test"]).lower().replace("χ", "chi").replace("chisq", "chi2").replace("chi-square", "chi2")
    if test == "chi":
        test = "chi2"
    if test not in TESTS:
        raise InputError(f"test '{claim['test']}' is not one of t, F, chi2, z, r")
    tails_text = clean_text(claim.get("tails")) or "2"
    if tails_text not in {"1", "2"}:
        raise InputError("tails must be 1 or 2")
    tails = int(tails_text)
    notes: list[str] = []
    if tails == 1 and test in {"f", "chi2"}:
        notes.append(f"tails = 1 ignored: {test} tests are upper-tail by construction")
        tails = 2
    statistic = parse_reported(claim["statistic"], "statistic")
    relation, p_reported = parse_p(claim["p"])

    stat_low, stat_high = statistic.low, statistic.high
    if test in {"t", "z", "r"}:
        if stat_low < 0 < stat_high:
            magnitudes = (0.0, max(abs(stat_low), abs(stat_high)))
        else:
            magnitudes = tuple(sorted((abs(stat_low), abs(stat_high))))
    else:
        if statistic.value < 0:
            raise InputError(f"{test} statistic cannot be negative")
        magnitudes = (max(stat_low, 0.0), stat_high)
    if test == "r" and magnitudes[1] >= 1:
        if statistic.value.copy_abs() >= 1:
            raise InputError("r must lie between -1 and 1")
        magnitudes = (magnitudes[0], min(magnitudes[1], 1 - 1e-15))

    def compute(stat: float, df_a: float | None, df_b: float | None, n_tails: int) -> float:
        if test == "t":
            return p_from_t(stat, df_a, n_tails)
        if test == "z":
            return p_from_z(stat, n_tails)
        if test == "r":
            return p_from_r(stat, df_a, n_tails)
        if test == "chi2":
            return p_from_chi2(stat, df_a)
        return p_from_f(stat, df_a, df_b)

    df_options: list[tuple[float | None, float | None]] = [(None, None)]
    df_text = ""
    if test in {"t", "chi2"}:
        need(claim, "df")
        low, high, df_text = _df_range(claim, "df")
        df_options = [(low, None), (high, None)]
    elif test == "r":
        if clean_text(claim.get("df")):
            low, high, df_text = _df_range(claim, "df")
        else:
            need(claim, "n")
            n = parse_count(claim["n"], "n", minimum=3)
            low = high = float(n - 2)
            df_text = f"{n - 2} (n - 2)"
        df_options = [(low, None), (high, None)]
    elif test == "f":
        need(claim, "df", "df2")
        low1, high1, text1 = _df_range(claim, "df")
        low2, high2, text2 = _df_range(claim, "df2")
        df_options = [(a, b) for a in (low1, high1) for b in (low2, high2)]
        df_text = f"{text1}, {text2}"

    def p_bounds(n_tails: int) -> tuple[float, float]:
        values = [compute(stat, a, b, n_tails) for stat in magnitudes for a, b in df_options]
        return min(values), max(values)

    p_low, p_high = p_bounds(tails)
    p_exact = compute(abs(float(statistic.value)), *df_options[0], tails) if test != "f" else compute(
        float(statistic.value), *df_options[0], tails)
    recomputed = {
        "test": test,
        "df": df_text,
        "tails": tails,
        "p_recomputed": sig(p_exact),
        "p_range_given_rounding": [sig(p_low), sig(p_high)],
    }
    p_label = f"p {relation} {p_reported.text}"
    if p_range_fits(relation, p_reported, p_low, p_high):
        return result(CONSISTENT, f"recomputed p {sig(p_low)} to {sig(p_high)} fits reported {p_label}", recomputed,
                      notes)

    claimed = reported_significant(relation, p_reported, alpha)
    computed_sig = True if p_high < alpha else False if p_low >= alpha else None
    decision_error = claimed is not None and computed_sig is not None and claimed != computed_sig
    if tails == 2 and test in {"t", "z", "r"}:
        one_low, one_high = p_bounds(1)
        if p_range_fits(relation, p_reported, one_low, one_high):
            notes.append("would be consistent with a one-tailed test; check whether the paper says one-tailed")
    if decision_error:
        notes.append(f"the significance decision at alpha = {alpha} changes")
    return result(INCONSISTENT, f"recomputed p {sig(p_low)} to {sig(p_high)} does not fit reported {p_label}",
                  recomputed, notes, decision_error)


# ---------------------------------------------------------------------------
# Check 2: GRIM (granularity of means of integer data)
# ---------------------------------------------------------------------------


def grim_possible(mean: Reported, cells: int) -> tuple[bool, list[int]]:
    """Is there an integer total S with S / cells inside the mean's rounding interval?

    Both interval ends are included, so either rounding convention at an exact
    half is accepted. Arithmetic is exact (Decimal) so 0.005 stays 0.005.
    """

    low = (mean.value - mean.half_unit) * cells
    high = (mean.value + mean.half_unit) * cells
    first = int(low.to_integral_value(rounding="ROUND_CEILING"))
    last = int(high.to_integral_value(rounding="ROUND_FLOOR"))
    return first <= last, list(range(first, last + 1))[:5]


def check_grim(claim: dict[str, Any]) -> dict[str, Any]:
    need(claim, "mean", "n")
    mean = parse_reported(claim["mean"], "mean")
    n = parse_count(claim["n"], "n", minimum=1)
    items = parse_count(claim.get("items") or "1", "items", minimum=1)
    cells = n * items
    if mean.decimals == 0:
        return result(NOT_CHECKABLE, "mean is reported without decimals, so every integer total fits")
    step = Decimal(1) / cells
    if step <= mean.half_unit * 2:
        return result(
            NOT_CHECKABLE,
            f"n x items = {cells} is too large for {mean.decimals} decimals: every value is a possible mean",
            {"granularity": float(step)},
        )
    possible, totals = grim_possible(mean, cells)
    recomputed = {"n_x_items": cells, "granularity": sig(float(step))}
    if possible:
        recomputed["example_total"] = totals[0]
        recomputed["example_mean"] = sig(totals[0] / cells)
        return result(CONSISTENT, f"a total of {totals[0]} over {cells} gives {totals[0] / cells:.6g}", recomputed)
    below = int(((mean.value) * cells).to_integral_value(rounding="ROUND_FLOOR"))
    recomputed["nearest_possible_means"] = [round(below / cells, mean.decimals + 2),
                                            round((below + 1) / cells, mean.decimals + 2)]
    notes = ["GRIM assumes every response is an integer; it does not apply to continuous measures",
             "check that n is the number of people contributing to this mean (missing data change it)"]
    return result(
        INCONSISTENT,
        f"no integer total over {cells} rounds to {mean.text}; nearest possible means are "
        f"{below / cells:.{mean.decimals + 2}f} and {(below + 1) / cells:.{mean.decimals + 2}f}",
        recomputed,
        notes,
    )


# ---------------------------------------------------------------------------
# Check 3: percentage against count and n
# ---------------------------------------------------------------------------


def check_percent(claim: dict[str, Any]) -> dict[str, Any]:
    need(claim, "percent", "n")
    percent = parse_reported(claim["percent"], "percent")
    n = parse_count(claim["n"], "n", minimum=1)
    if not 0 <= percent.value <= 100:
        raise InputError("percent must lie between 0 and 100")
    low = (percent.value - percent.half_unit) * n / 100
    high = (percent.value + percent.half_unit) * n / 100
    if clean_text(claim.get("count")):
        count = parse_count(claim["count"], "count")
        if count > n:
            return result(INCONSISTENT, f"count {count} is larger than n {n}", {"count": count, "n": n})
        exact = 100 * count / n
        recomputed = {"percent_recomputed": round(exact, percent.decimals + 2)}
        if low <= count <= high:
            return result(CONSISTENT, f"{count}/{n} = {exact:.{percent.decimals + 2}f}% rounds to {percent.text}%",
                          recomputed)
        return result(INCONSISTENT, f"{count}/{n} = {exact:.{percent.decimals + 2}f}%, not {percent.text}%",
                      recomputed, ["check that n is the denominator used for this percentage"])
    step = Decimal(100) / n
    if step <= percent.half_unit * 2:
        return result(NOT_CHECKABLE, f"no count given and n = {n} is fine-grained enough that any "
                                     f"{percent.decimals}-decimal percentage is possible")
    first = int(low.to_integral_value(rounding="ROUND_CEILING"))
    last = int(high.to_integral_value(rounding="ROUND_FLOOR"))
    if first <= last:
        return result(CONSISTENT, f"a count of {first} out of {n} gives {100 * first / n:.4g}%",
                      {"possible_counts": list(range(first, last + 1))[:5]})
    return result(INCONSISTENT, f"no whole count out of {n} rounds to {percent.text}%",
                  {"nearest_percents": [round(100 * math.floor(low) / n, 3), round(100 * math.ceil(high) / n, 3)]},
                  ["check that n is the denominator used for this percentage"])


# ---------------------------------------------------------------------------
# Check 4: parts that should add up to a total
# ---------------------------------------------------------------------------


def check_sum(claim: dict[str, Any]) -> dict[str, Any]:
    need(claim, "parts", "total")
    pieces = [piece for piece in re.split(r"[;|]", clean_text(claim["parts"])) if piece.strip()]
    if len(pieces) < 2:
        raise InputError("parts needs at least two values separated by ';'")
    parts = [parse_reported(piece.strip(), "parts") for piece in pieces]
    total = parse_reported(claim["total"], "total")
    # A whole number in a sum is a count (or 100%), which is exact; decimals are rounded values.
    parts = [parse_reported(item.text, "parts", exact=True) if item.decimals == 0 else item for item in parts]
    if total.decimals == 0:
        total = parse_reported(total.text, "total", exact=True)
    exact_sum = sum((part.value for part in parts), Decimal(0))
    slack = sum((part.half_unit for part in parts), Decimal(0)) + total.half_unit
    recomputed = {"sum": float(exact_sum), "allowed_rounding_slack": float(slack)}
    if abs(exact_sum - total.value) <= slack:
        why = "exactly" if slack == 0 else f"within rounding slack {slack}"
        return result(CONSISTENT, f"parts sum to {exact_sum}, matching {total.text} {why}", recomputed)
    return result(INCONSISTENT, f"parts sum to {exact_sum}, not {total.text} (allowed rounding slack {slack})",
                  recomputed, ["exclusions, missing data, or overlapping categories can explain a gap; "
                               "look for them in the paper before concluding"])


# ---------------------------------------------------------------------------
# Checks 5 and 6: estimate with a confidence interval
# ---------------------------------------------------------------------------

SYMMETRIC_METHODS = {"", "wald", "normal", "t", "regression", "model"}


def _ci_check(claim: dict[str, Any], ratio: bool, alpha_default: float) -> dict[str, Any]:
    need(claim, "estimate", "ci_low", "ci_high")
    estimate = parse_reported(claim["estimate"], "estimate")
    low = parse_reported(claim["ci_low"], "ci_low")
    high = parse_reported(claim["ci_high"], "ci_high")
    level_text = clean_text(claim.get("level")) or "95"
    level = float(parse_reported(level_text, "level").value)
    if not 50 <= level < 100:
        raise InputError("level is a percentage such as 95")
    alpha = 1 - level / 100
    method = clean_text(claim.get("ci_method")).lower()
    if low.value > high.value:
        raise InputError("ci_low is larger than ci_high")
    if ratio and (estimate.value <= 0 or low.value <= 0):
        if low.value == 0 and low.half_unit > 0:
            pass
        else:
            raise InputError("ratio estimates and interval bounds must be positive")
    null = 1.0 if ratio else 0.0
    notes: list[str] = []
    recomputed: dict[str, Any] = {"scale": "log" if ratio else "linear", "level": level}
    problems: list[str] = []

    if estimate.high < low.low or estimate.low > high.high:
        problems.append(f"estimate {estimate.text} lies outside its own interval {low.text} to {high.text}")

    def transform(value: float) -> float:
        return math.log(max(value, 1e-300)) if ratio else value

    def back(value: float) -> float:
        return math.exp(value) if ratio else value

    low_ends = (max(low.low, 1e-12) if ratio else low.low, low.high)
    high_ends = (high.low, high.high)
    centres = [back((transform(a) + transform(b)) / 2) for a in low_ends for b in high_ends]
    centre_low, centre_high = min(centres), max(centres)
    recomputed["interval_centre"] = sig(back((transform(float(low.value) if float(low.value) > 0 else low_ends[0])
                                              + transform(float(high.value))) / 2))
    recomputed["centre_range_given_rounding"] = [sig(centre_low), sig(centre_high)]
    if method in SYMMETRIC_METHODS:
        if not problems and (estimate.high < centre_low or estimate.low > centre_high):
            scale = "log scale (geometric mean of the bounds)" if ratio else "linear scale (midpoint of the bounds)"
            problems.append(
                f"estimate {estimate.text} is not at the centre of its interval on the {scale}: "
                f"centre {sig(centre_low)} to {sig(centre_high)}")
            notes.append("exact, profile-likelihood, score, or bootstrap intervals are legitimately asymmetric; "
                         "set ci_method if the paper used one")
    else:
        notes.append(f"ci_method '{method}' is not symmetric, so only containment and the p agreement were checked")

    excludes_null = low.low > null or high.high < null
    touches_null = low.low <= null <= low.high or high.low <= null <= high.high
    recomputed["interval_excludes_null"] = None if touches_null else excludes_null

    p_text = clean_text(claim.get("p"))
    decision_error = False
    if p_text:
        relation, p_reported = parse_p(p_text)
        z_crit = STANDARD_NORMAL.inv_cdf(1 - alpha / 2)
        df_text = clean_text(claim.get("df"))
        if df_text and not ratio:
            z_crit = t_quantile(1 - alpha / 2, float(parse_reported(df_text, "df").value))
        width = transform(float(high.value)) - transform(max(float(low.value), 1e-12))
        if width > 0:
            se = width / (2 * z_crit)
            stat = transform(float(estimate.value)) / se
            if df_text and not ratio:
                p_ci = p_from_t(stat, float(parse_reported(df_text, "df").value))
            else:
                p_ci = p_from_z(stat)
            recomputed["p_from_interval"] = sig(p_ci)
            if not p_range_fits(relation, p_reported, p_ci * 0.5, p_ci * 2.0):
                notes.append(f"p from the interval is about {sig(p_ci)} versus reported p {relation} {p_reported.text}; "
                             "this is an approximation (Altman and Bland), so it is noted, not flagged")
        claimed = reported_significant(relation, p_reported, 1 - level / 100)
        if claimed is not None and not touches_null and claimed != excludes_null:
            decision_error = True
            side = "excludes" if excludes_null else "includes"
            problems.append(f"the {level:g}% interval {side} the null value {null:g} but p {relation} {p_reported.text}")
        elif touches_null:
            notes.append("an interval bound equals the null value within rounding, so the p agreement cannot be judged")
        elif claimed is None:
            notes.append(f"reported p cannot be compared with a {level:g}% interval")

    if problems:
        return result(INCONSISTENT, "; ".join(problems), recomputed, notes, decision_error)
    return result(CONSISTENT, "estimate, interval" + (" and p agree" if p_text else " agree"), recomputed, notes)


def check_ratio_ci(claim: dict[str, Any]) -> dict[str, Any]:
    return _ci_check(claim, ratio=True, alpha_default=ALPHA_DEFAULT)


def check_diff_ci(claim: dict[str, Any]) -> dict[str, Any]:
    return _ci_check(claim, ratio=False, alpha_default=ALPHA_DEFAULT)


CHECKS = {
    "p_value": check_p_value,
    "grim": check_grim,
    "percent": check_percent,
    "sum": check_sum,
    "ratio_ci": check_ratio_ci,
    "diff_ci": check_diff_ci,
}
