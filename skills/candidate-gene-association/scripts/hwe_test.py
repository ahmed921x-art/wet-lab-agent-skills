#!/usr/bin/env python3
"""Hardy-Weinberg test for one biallelic variant from three genotype counts.

Reports the exact test, the chi-square test, the minor allele frequency, the
inbreeding coefficient F, and whether any deviation is a heterozygote deficit or
excess. Give counts from controls: a deviation in cases can be real association.
"""

from __future__ import annotations

import argparse
import json
import sys

from _genetics import InputError, hardy_weinberg, safe_console


def parse_counts(text: str) -> tuple[int, int, int]:
    """Parse HOM,HET,HOM."""

    parts = text.split(",")
    try:
        counts = tuple(int(part) for part in parts)
    except ValueError as exc:
        raise InputError(f"counts must be integers, got {text!r}") from exc
    if len(counts) != 3 or min(counts) < 0:
        raise InputError(f"give three non-negative counts as HOM,HET,HOM, got {text!r}")
    return counts  # type: ignore[return-value]


def main(argv: list[str] | None = None) -> int:
    safe_console()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--counts", required=True, help="HOM,HET,HOM, heterozygotes in the middle")
    parser.add_argument("--alpha", type=float, default=0.05, help="flag exact p below this")
    parser.add_argument("--json", action="store_true", help="print JSON only")
    args = parser.parse_args(argv)

    try:
        result = hardy_weinberg(*parse_counts(args.counts))
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2

    if args.json:
        print(json.dumps(result, indent=2))
        return 0
    marker = "[WARN]" if result["exact_p"] < args.alpha else "[OK]"
    print(f"{marker} n={result['n']} MAF={result['minor_allele_frequency']} exact p={result['exact_p']:.4g}")
    print(f"  -> observed {result['observed']}, expected {result['expected']}")
    if result["chi_square"] is not None:
        print(f"  -> chi-square {result['chi_square']} (p={result['chi_square_p']:.4g})")
    print(f"  -> F={result['inbreeding_f']} ({result['direction']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
