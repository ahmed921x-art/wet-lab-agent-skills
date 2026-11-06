---
name: research-ai-router
description: Picks the right tool for a research task from 16 reviewed research-AI repositories (paper-qa, deepchem, pyqtgraph, bioicons, draw.io illustrator, K-Dense skills, OpenScience, AIPOCH medical skills, Paper2Code, Paper2Video, AI-Researcher, awesome-ai-for-science, and others) plus this repository's paper-audit skill, and gives the official install command, license, reuse terms, requirements, and security cautions for each. Never recommends tools excluded for license reasons and names an allowed alternative instead. Use when asked which AI tool or skill to use for literature questions, molecules, interactive plots, scientific icons or diagrams, paper-to-code, paper-to-video, automated research, or checking a paper's numbers, in English or Arabic.
license: MIT
compatibility: Requires Python 3.11+ standard library only. Tested on Python 3.13. No network access; the catalog is a snapshot reviewed on 2026-11-06.
allowed-tools: Read Bash
metadata:
  version: "1.0"
  skill-author: ahmed921x-art
  last-reviewed: "2026-11-06"
---

# Research AI Router

## When to use

Use this skill when someone asks which research-AI tool fits a task, or how to install one safely. It points to
the original repositories. **It contains no code from any of them.** Each tool keeps its own license, and the
user installs it from its source.

## The three rules

**1. Link, never copy.** Recommend the official repository and its own install command. Do not paste or vendor
code from these repositories into the user's project unless that tool's license allows it, and then keep the
license and notices.

**2. Excluded means excluded.** SakanaAI/AI-Scientist, SakanaAI/AI-Scientist-v2, and hussius/deeplearning-biology
are never recommended. The first two have use-restricting license terms, and the last has no license. If a
task matches one, say so, give the reason, and offer the allowed alternative from the catalog.

**3. Say the risks out loud.** Every suggestion carries its cautions, such as API keys, code execution, pickled
files, or per-icon licenses. Repeat them to the user instead of hiding them.

## Workflow

### 1. Ask the catalog

```bash
python scripts/pick_tool.py --task "predict toxicity of a molecule from SMILES"
python scripts/pick_tool.py --task "ابحث في أوراقي عن مقاومة المضادات"
python scripts/pick_tool.py --list
```

Add `--json` for structured output and `--limit` to show more or fewer suggestions (default 3). The script
matches English keywords as whole words (plural `s` allowed) and Arabic keywords as substrings, after folding
letter variants such as أ, إ, and آ.

### 2. Read the suggestion before passing it on

| Status | Meaning | What to tell the user |
| --- | --- | --- |
| recommended | Permissive license, fits the task, low risk | Install command, license, needs |
| caution | Usable, but untested, heavy, or risky | The cautions, and to run it in a sandbox or check its output |
| excluded | License forbids reuse here | Do not install; give the alternative |

If nothing matches, the script says so. Do not stretch a tool to fit. Say no reviewed tool covers the task.

### 3. Match the task to the tool family

| Task | First choice | Notes |
| --- | --- | --- |
| Check whether a paper's numbers agree | `paper-audit` (this repository) | Inconsistencies are questions, not accusations |
| Ask questions over your own PDFs | paper-qa | Needs model keys, or a local embedding model |
| Molecules, drugs, proteins, materials | deepchem | Choose a back end: torch, tensorflow, or jax |
| A specific analysis skill (single cell, statistics, databases) | K-Dense scientific-agent-skills | Prefer skills that have tests |
| Medical evidence appraisal or protocols | AIPOCH medical-research-skills | Check any numbers against statsmodels |
| Icons for figures | bioicons | Cite every icon; watch share-alike licenses |
| Diagrams and schematics | drawio-scientific-illustrator | Needs the draw.io desktop app |
| Live or very large interactive plots | pyqtgraph | Desktop only; use matplotlib for papers |
| Turn a paper into code or a video | Paper2Code, Paper2Video | No tests; sandbox the generated code |
| Find other tools | awesome-ai-for-science | Links are not checked |

The wet-lab skills in this repository (qPCR, ELISA, Western blot, candidate gene) are not in the catalog,
because they are chosen by the kind of data rather than by tool.

## Limits

- The catalog is a snapshot from 2026-11-06. Licenses, install commands, and maintenance can change, so check the
  repository before relying on an old entry.
- Keyword matching is simple. It can miss a task phrased in unexpected words, so use `--list` to browse.
- Scores and cautions come from sampled files, not full audits.

## Sources

The evidence behind every status and caution is in `references/catalog-notes.md`, which points to the full
comparison report.
