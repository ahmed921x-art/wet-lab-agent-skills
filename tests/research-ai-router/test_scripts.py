#!/usr/bin/env python3
"""Tests for research-ai-router: catalog integrity, routing, exclusions, and bad input."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

import skill_contract  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
SKILL_ROOT = ROOT / "skills" / "research-ai-router"
SCRIPTS = SKILL_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import pick_tool  # noqa: E402

CliHelpTests = skill_contract.cli.help_test_case(SKILL_ROOT)

# The 14 repositories the owner asked about, plus the two reference collections.
REVIEWED = {
    "Future-House/paper-qa", "SakanaAI/AI-Scientist-v2", "SakanaAI/AI-Scientist", "HKUDS/AI-Researcher",
    "zLanqing/codex-claude-academic-skills", "going-doer/Paper2Code", "deepchem/deepchem",
    "hussius/deeplearning-biology", "ai4s-research/awesome-ai-for-science", "aipoch/medical-research-skills",
    "duerrsimon/bioicons", "icebird1998/drawio-scientific-illustrator", "showlab/Paper2Video", "pyqtgraph/pyqtgraph",
    "K-Dense-AI/scientific-agent-skills", "synthetic-sciences/openscience",
}
EXCLUDED = {"SakanaAI/AI-Scientist", "SakanaAI/AI-Scientist-v2", "hussius/deeplearning-biology"}


@pytest.fixture(scope="module")
def tools():
    return pick_tool.load_catalog()


def top(task, tools):
    result = pick_tool.pick(task, tools, limit=1)
    return result["suggestions"][0]["id"] if result["suggestions"] else None


def test_catalog_covers_every_reviewed_repository_once(tools):
    repos = [tool["repo"] for tool in tools if "/" in tool["repo"] and " " not in tool["repo"]]
    assert set(repos) == REVIEWED
    assert len(repos) == len(set(repos))


def test_catalog_matches_comparison_report(tools):
    report = (ROOT / "docs" / "repo-comparison.md").read_text(encoding="utf-8")
    for repo in REVIEWED:
        assert repo in report, repo


def test_excluded_repositories_are_marked_and_have_alternatives(tools):
    by_repo = {tool["repo"]: tool for tool in tools}
    for repo in EXCLUDED:
        assert by_repo[repo]["status"] == "excluded"
        assert by_repo[repo]["reuse"] == "no"
        alternative = by_repo[repo]["alternative"].split(" ")[0]
        assert alternative in {tool["id"] for tool in tools if tool["status"] != "excluded"}
    assert {tool["repo"] for tool in tools if tool["status"] == "excluded"} == EXCLUDED


def test_every_entry_is_complete(tools):
    for tool in tools:
        assert tool["keywords"] and tool["install"] and tool["license"] and tool["what"], tool["id"]
        assert isinstance(tool["needs"], list) and isinstance(tool["cautions"], list)


def test_no_install_command_pipes_an_unread_script(tools):
    for tool in tools:
        assert not re.search(r"\|\s*(ba)?sh", tool["install"]), tool["id"]


@pytest.mark.parametrize(
    "task,expected",
    [
        ("ask questions about my folder of PDFs", "paper-qa"),
        ("ابحث في أوراقي عن مقاومة المضادات", "paper-qa"),
        ("predict toxicity of a molecule from SMILES", "deepchem"),
        ("صمم جزيء دواء جديد", "deepchem"),
        ("an icon of a cell for my graphical abstract", "bioicons"),
        ("أبي أيقونة خلية", "bioicons"),
        ("draw a schematic diagram of the pathway", "drawio-scientific-illustrator"),
        ("plot streaming live data in real-time", "pyqtgraph"),
        ("check whether the p-values in this paper add up", "paper-audit"),
        ("تدقيق أرقام الورقة", "paper-audit"),
        ("single cell rna-seq analysis with scanpy", "k-dense-scientific-agent-skills"),
        ("critical appraisal for a systematic review of a clinical trial", "aipoch-medical-research-skills"),
        ("implement paper as code to reproduce paper results", "paper2code"),
        ("make a presentation video from my paper", "paper2video"),
        ("find a tool in a list of tools", "awesome-ai-for-science"),
    ],
)
def test_routing(task, expected, tools):
    assert top(task, tools) == expected


@pytest.mark.parametrize("task", [
    "write a paper automatically like the AI Scientist",
    "use ai-scientist-v2 agentic tree search",
    "deep learning biology papers list",
])
def test_excluded_tools_are_never_suggested(task, tools):
    result = pick_tool.pick(task, tools, limit=20)
    assert result["excluded_matches"]
    suggested = {item["id"] for item in result["suggestions"]}
    excluded_ids = {tool["id"] for tool in tools if tool["status"] == "excluded"}
    assert not suggested & excluded_ids
    assert all(item["alternative"] for item in result["excluded_matches"])


def test_unrelated_task_gets_no_suggestion(tools):
    result = pick_tool.pick("bake a chocolate cake", tools)
    assert result["suggestions"] == [] and result["excluded_matches"] == []
    assert "No reviewed tool" in pick_tool.render(result)


def test_whole_word_matching_avoids_false_hits(tools):
    # The pyqtgraph keyword "qt" must not match inside "qtl".
    assert pick_tool.matches("map the qtl peaks", "qt") is False
    assert pick_tool.matches("a qt window", "qt") is True
    assert top("map the qtl peaks", tools) is None


def test_arabic_letter_variants_are_folded():
    assert pick_tool.normalize("أَيقونة") == pick_tool.normalize("ايقونه")


def test_limit_and_ranking(tools):
    task = "molecule toxicity smiles, an icon, a diagram, and a presentation video"
    assert len(pick_tool.pick(task, tools, limit=10)["suggestions"]) >= 4
    result = pick_tool.pick(task, tools, limit=2)
    assert len(result["suggestions"]) == 2
    assert result["suggestions"][0]["id"] == "deepchem"  # three keyword hits beat the others
    scores = [item["score"] for item in result["suggestions"]]
    assert scores == sorted(scores, reverse=True)


@pytest.mark.parametrize("task,limit", [("", 3), ("   ", 3), ("x" * 3000, 3), ("molecule", 0)])
def test_bad_input_is_rejected(task, limit, tools):
    with pytest.raises(pick_tool.InputError):
        pick_tool.pick(task, tools, limit)


def test_bad_catalogs_are_rejected(tmp_path):
    cases = {
        "not json": "{",
        "no tools": json.dumps({"tools": []}),
        "missing field": json.dumps({"tools": [{"id": "a"}]}),
    }
    good = json.loads((SKILL_ROOT / "assets" / "tools.json").read_text(encoding="utf-8"))["tools"][0]
    cases["bad status"] = json.dumps({"tools": [dict(good, status="maybe")]})
    cases["duplicate"] = json.dumps({"tools": [good, good]})
    for name, content in cases.items():
        path = tmp_path / f"{name}.json"
        path.write_text(content, encoding="utf-8")
        with pytest.raises(pick_tool.InputError):
            pick_tool.load_catalog(path)
    with pytest.raises(pick_tool.InputError):
        pick_tool.load_catalog(tmp_path / "missing.json")


def test_cli_text_json_and_list():
    def run(*args):
        return subprocess.run([sys.executable, str(SCRIPTS / "pick_tool.py"), *args], capture_output=True,
                              text=True, encoding="utf-8", timeout=120)

    done = run("--task", "predict toxicity of a molecule")
    assert done.returncode == 0 and "pip install deepchem" in done.stdout and "license: MIT" in done.stdout
    done = run("--task", "ابحث في أوراقي", "--json")
    assert done.returncode == 0 and json.loads(done.stdout)["suggestions"][0]["id"] == "paper-qa"
    done = run("--list")
    assert done.returncode == 0 and len(done.stdout.strip().splitlines()) == 17
    done = run("--task", "")
    assert done.returncode == 2 and "[FAIL]" in done.stderr
