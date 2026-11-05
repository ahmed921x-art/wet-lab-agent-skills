"""Registers the small shared contract that every skill's test suite imports."""

from __future__ import annotations

import subprocess
import sys
import types
import unittest
from pathlib import Path


def help_test_case(skill_root: Path) -> type[unittest.TestCase]:
    """Build a test case checking that each public script answers --help."""

    class HelpTests(unittest.TestCase):
        def test_every_cli_answers_help(self) -> None:
            scripts = sorted(path for path in (skill_root / "scripts").glob("*.py") if not path.name.startswith("_"))
            self.assertTrue(scripts)
            for script in scripts:
                with self.subTest(script=script.name):
                    done = subprocess.run(
                        [sys.executable, str(script), "--help"], capture_output=True, text=True, timeout=120
                    )
                    self.assertEqual(done.returncode, 0, done.stderr)

    return HelpTests


contract = types.ModuleType("skill_contract")
contract.cli = types.SimpleNamespace(help_test_case=help_test_case)
sys.modules["skill_contract"] = contract
