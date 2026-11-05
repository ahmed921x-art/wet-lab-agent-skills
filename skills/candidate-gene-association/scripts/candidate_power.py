#!/usr/bin/env python3
"""Power and sample size for a case-control test of one biallelic variant.

Uses the normal approximation for comparing effect-allele frequency between cases
and controls under an additive per-allele odds ratio. Three questions:
  power for a planned study    --cases, --controls, and --odds-ratio
  smallest detectable OR       --cases and --controls
  people needed per group      --odds-ratio only (equal groups)
--inbreeding-f inflates allele-count variance by (1 + F), as in consanguineous samples.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from statistics import NormalDist
from typing import Any

from _genetics import InputError, safe_console

NORMAL = NormalDist()


def case_frequency(control_frequency: float, odds_ratio: float) -> float:
    """Effect-allele frequency in cases implied by a per-allele odds ratio."""

    odds = odds_ratio * control_frequency / (1.0 - control_frequency)
    return odds / (1.0 + odds)


def power(
    cases: int,
    controls: int,
    frequency: float,
    odds_ratio: float,
    alpha: float = 0.05,
    inbreeding_f: float = 0.0,
) -> float:
    """Two-sided power of the allele-frequency comparison."""

    p_control = frequency
    p_case = case_frequency(p_control, odds_ratio)
    pooled = (cases * p_case + controls * p_control) / (cases + controls)
    inflation = 1.0 + inbreeding_f
    null_se = math.sqrt(pooled * (1 - pooled) * inflation * (1 / (2 * cases) + 1 / (2 * controls)))
    alt_se = math.sqrt(
        inflation * (p_case * (1 - p_case) / (2 * cases) + p_control * (1 - p_control) / (2 * controls))
    )
    critical = NORMAL.inv_cdf(1 - alpha / 2)
    gap = abs(p_case - p_control)
    return NORMAL.cdf((gap - critical * null_se) / alt_se)


def minimum_detectable(
    cases: int,
    controls: int,
    frequency: float,
    target: float,
    alpha: float,
    inbreeding_f: float,
    direction: float,
) -> float | None:
    """Smallest OR away from 1 reaching target power, by bisection on log OR."""

    low, high = 0.0, math.log(1000.0)
    if power(cases, controls, frequency, math.exp(direction * high), alpha, inbreeding_f) < target:
        return None
    for _ in range(80):
        middle = (low + high) / 2
        if power(cases, controls, frequency, math.exp(direction * middle), alpha, inbreeding_f) < target:
            low = middle
        else:
            high = middle
    return math.exp(direction * high)


def required_per_group(
    frequency: float, odds_ratio: float, target: float, alpha: float, inbreeding_f: float
) -> int:
    """Smallest equal group size reaching target power."""

    if odds_ratio == 1.0:
        raise InputError("an odds ratio of 1 has no detectable effect")
    low, high = 2, 2
    while power(high, high, frequency, odds_ratio, alpha, inbreeding_f) < target:
        high *= 2
        if high > 10**9:
            raise InputError("the required sample size is beyond 10^9 per group")
    while low < high:
        middle = (low + high) // 2
        if power(middle, middle, frequency, odds_ratio, alpha, inbreeding_f) < target:
            low = middle + 1
        else:
            high = middle
    return low


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--frequency", type=float, required=True, help="effect-allele frequency in controls")
    parser.add_argument("--cases", type=int)
    parser.add_argument("--controls", type=int)
    parser.add_argument("--odds-ratio", type=float, help="per-allele odds ratio")
    parser.add_argument("--alpha", type=float, default=0.05, help="after any correction, e.g. 0.05 / variants")
    parser.add_argument("--power", type=float, default=0.8, dest="target", help="target power")
    parser.add_argument("--inbreeding-f", type=float, default=0.0)
    args = parser.parse_args(argv)

    try:
        if not 0 < args.frequency < 1:
            raise InputError("--frequency must be between 0 and 1")
        if not 0 < args.alpha < 1 or not 0 < args.target < 1:
            raise InputError("--alpha and --power must be between 0 and 1")
        if not 0 <= args.inbreeding_f < 1:
            raise InputError("--inbreeding-f must be between 0 and 1")
        if (args.cases is None) != (args.controls is None):
            raise InputError("give both --cases and --controls, or neither")
        if args.cases is not None and min(args.cases, args.controls) < 1:
            raise InputError("group sizes must be positive")
        if args.odds_ratio is not None and args.odds_ratio <= 0:
            raise InputError("--odds-ratio must be positive")
        options = (args.alpha, args.inbreeding_f)
        result: dict[str, Any] = {"frequency": args.frequency, "alpha": args.alpha}
        if args.cases is not None and args.odds_ratio is not None:
            result["power"] = round(power(args.cases, args.controls, args.frequency, args.odds_ratio, *options), 4)
            print(f"[OK] power {result['power']:.1%} for OR {args.odds_ratio}")
        elif args.cases is not None:
            up = minimum_detectable(args.cases, args.controls, args.frequency, args.target, *options, 1.0)
            down = minimum_detectable(args.cases, args.controls, args.frequency, args.target, *options, -1.0)
            result["min_detectable_or"] = {
                "risk": None if up is None else round(up, 3),
                "protective": None if down is None else round(down, 3),
            }
            print(f"[OK] at {args.target:.0%} power, detectable OR: {result['min_detectable_or']}")
        elif args.odds_ratio is not None:
            count = required_per_group(args.frequency, args.odds_ratio, args.target, *options)
            result["required_per_group"] = count
            print(f"[OK] {count} cases and {count} controls for {args.target:.0%} power")
        else:
            raise InputError("give --odds-ratio, or --cases and --controls")
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
