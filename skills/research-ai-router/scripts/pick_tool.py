#!/usr/bin/env python3
"""Suggest which of the reviewed research-AI tools fits a task, with install, license, and cautions.

Reads assets/tools.json, scores each tool by how many of its keywords appear in the
task (English or Arabic), and prints the best matches. Tools excluded for license
reasons are never recommended: if one matches, it is named with the reason and,
where the catalog gives one, an allowed alternative. Standard library only; no network access.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import unicodedata
from pathlib import Path
from typing import Any

CATALOG = Path(__file__).resolve().parents[1] / "assets" / "tools.json"
STATUSES = {"recommended", "caution", "excluded"}
REQUIRED = {"id", "repo", "kind", "status", "what", "keywords", "install", "license", "reuse", "needs", "cautions"}
MAX_TASK = 2000


class InputError(ValueError):
    """Raised for input problems the user has to resolve."""


def normalize(text: str) -> str:
    """Lower-case, fold Arabic letter variants, drop diacritics, and collapse spaces."""

    text = unicodedata.normalize("NFKC", text).lower()
    text = re.sub(r"[ً-ٰٟـ]", "", text)  # harakat, dagger alef, tatweel
    text = text.translate(str.maketrans({"أ": "ا", "إ": "ا", "آ": "ا", "ة": "ه", "ى": "ي"}))
    text = re.sub(r"[_\-/]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def load_catalog(path: str | Path = CATALOG) -> list[dict[str, Any]]:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise InputError(f"no catalog at {path}") from exc
    except json.JSONDecodeError as exc:
        raise InputError(f"catalog {path} is not valid JSON: {exc}") from exc
    tools = data.get("tools") if isinstance(data, dict) else None
    if not isinstance(tools, list) or not tools:
        raise InputError("catalog has no tools list")
    seen = set()
    for tool in tools:
        missing = REQUIRED - set(tool)
        if missing:
            raise InputError(f"tool {tool.get('id', '?')} is missing {', '.join(sorted(missing))}")
        if tool["status"] not in STATUSES:
            raise InputError(f"tool {tool['id']} has unknown status {tool['status']}")
        if tool["id"] in seen:
            raise InputError(f"tool id {tool['id']} appears twice")
        seen.add(tool["id"])
    return tools


def matches(task: str, keyword: str) -> bool:
    """Whole-word match for Latin keywords; substring match for Arabic, which attaches prefixes."""

    key = normalize(keyword)
    if re.search(r"[؀-ۿ]", key):
        return key in task
    return re.search(rf"(?<![a-z0-9]){re.escape(key)}(?:s|es)?(?![a-z0-9])", task) is not None


def score(task: str, tool: dict[str, Any]) -> tuple[int, list[str]]:
    hits = [keyword for keyword in tool["keywords"] if matches(task, keyword)]
    # Longer phrases are more specific, so they count for more.
    return sum(1 + normalize(hit).count(" ") for hit in hits), hits


def pick(task: str, tools: list[dict[str, Any]], limit: int = 3) -> dict[str, Any]:
    text = normalize(task or "")
    if not text:
        raise InputError("describe the task in a few words")
    if len(text) > MAX_TASK:
        raise InputError(f"task description is longer than {MAX_TASK} characters")
    if limit < 1:
        raise InputError("limit must be at least 1")
    scored = []
    for order, tool in enumerate(tools):
        points, hits = score(text, tool)
        if points:
            scored.append((-points, order, tool, hits))
    scored.sort(key=lambda item: (item[0], item[1]))
    allowed = [(tool, hits, -points) for points, _, tool, hits in scored if tool["status"] != "excluded"]
    excluded = [(tool, hits) for _, _, tool, hits in scored if tool["status"] == "excluded"]
    return {
        "task": task,
        "suggestions": [
            {"id": tool["id"], "score": points, "matched": hits, **{key: tool[key] for key in sorted(REQUIRED - {"keywords"})}}
            for tool, hits, points in allowed[:limit]
        ],
        "excluded_matches": [
            {"id": tool["id"], "repo": tool["repo"], "matched": hits, "reason": "; ".join(tool["cautions"]),
             "license": tool["license"], "alternative": tool.get("alternative", "")}
            for tool, hits in excluded
        ],
    }


def render(result: dict[str, Any]) -> str:
    lines = []
    if not result["suggestions"]:
        lines.append("[NONE] No reviewed tool matches this task. Describe it differently, or use --list.")
    for rank, tool in enumerate(result["suggestions"], start=1):
        flag = "" if tool["status"] == "recommended" else " (use with caution)"
        lines += [
            f"{rank}. {tool['id']}{flag}: {tool['what']}",
            f"   repo:    https://github.com/{tool['repo'].split(' ')[0]}",
            f"   install: {tool['install']}",
            f"   license: {tool['license']}; reuse: {tool['reuse']}",
        ]
        if tool["needs"]:
            lines.append(f"   needs:   {'; '.join(tool['needs'])}")
        for caution in tool["cautions"]:
            lines.append(f"   caution: {caution}")
        lines.append(f"   matched: {', '.join(tool['matched'])}")
    for tool in result["excluded_matches"]:
        lines.append(f"[EXCLUDED] {tool['id']} ({tool['repo']}) matches but is not recommended: {tool['reason']}.")
        if tool["alternative"]:
            lines.append(f"           allowed alternative: {tool['alternative']}")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="replace")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--task", help="what you want to do, in English or Arabic")
    group.add_argument("--list", action="store_true", help="list every tool in the catalog")
    parser.add_argument("--limit", type=int, default=3, help="how many suggestions to show")
    parser.add_argument("--json", action="store_true", help="print JSON instead of text")
    parser.add_argument("--catalog", default=str(CATALOG), help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    try:
        tools = load_catalog(args.catalog)
        if args.list:
            for tool in tools:
                print(f"{tool['status']:<12} {tool['id']:<34} {tool['kind']:<17} {tool['what']}")
            return 0
        result = pick(args.task, tools, args.limit)
    except InputError as exc:
        print(f"[FAIL] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, ensure_ascii=False) if args.json else render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
