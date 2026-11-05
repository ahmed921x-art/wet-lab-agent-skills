"""Repository-level checks that need no scientific packages."""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SKILLS = sorted(path for path in (ROOT / "skills").iterdir() if path.is_dir())
ALLOWED = {"name", "description", "license", "compatibility", "allowed-tools", "metadata"}


def frontmatter(text: str) -> dict[str, str]:
    block = text.split("---\n", 2)[1]
    return dict(re.findall(r"^([a-z-]+):\s*(.*)$", block, re.M))


def test_manifest_is_valid_json():
    assert json.loads((ROOT / "plugin.json").read_text(encoding="utf-8"))["license"] == "MIT"


@pytest.mark.parametrize("skill", SKILLS, ids=lambda path: path.name)
def test_skill_frontmatter_and_size(skill):
    text = (skill / "SKILL.md").read_text(encoding="utf-8")
    fields = frontmatter(text)
    assert fields["name"] == skill.name
    assert 0 < len(fields["description"]) <= 1024
    assert set(fields) <= ALLOWED
    assert re.search(r'^  version: "\d+\.\d+"$', text, re.M)
    assert text.count("\n") < 500


@pytest.mark.parametrize("skill", SKILLS, ids=lambda path: path.name)
def test_documented_files_exist(skill):
    text = (skill / "SKILL.md").read_text(encoding="utf-8")
    for match in re.findall(r"(?:assets|references|scripts)/[A-Za-z0-9_./-]+\.[a-z]+", text):
        assert (skill / match).is_file(), match


@pytest.mark.parametrize("skill", SKILLS, ids=lambda path: path.name)
def test_scripts_parse_and_ship_no_tests(skill):
    for script in (skill / "scripts").glob("*.py"):
        ast.parse(script.read_text(encoding="utf-8"))
        assert not script.name.startswith("test_")
