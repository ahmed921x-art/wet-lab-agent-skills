#!/usr/bin/env python3
"""Run the structure checks, then each skill's suite in its own process."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    targets = ["tests/test_structure.py"] + [
        f"tests/{path.name}" for path in sorted((ROOT / "tests").iterdir()) if (path / "test_scripts.py").is_file()
    ]
    failed = []
    for target in targets:
        print(f"=== {target}", flush=True)
        done = subprocess.run(
            [sys.executable, "-m", "pytest", target, "-q", "--import-mode=importlib", "-p", "no:cacheprovider"],
            cwd=ROOT,
        )
        if done.returncode != 0:
            failed.append(target)
    print(f"{len(targets) - len(failed)} passed, {len(failed)} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
