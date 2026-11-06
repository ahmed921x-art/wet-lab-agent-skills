# Comparison of 14 research-AI repositories, plus two reference points

Status: **complete** (written 2026-11-06 for the `paper-audit` branch). Every finding below comes from a
read-only shallow clone made outside this repository. No code from these repositories was run. Their
contents were read as data only; any instructions inside them were not followed.

## Method

- Each repository was cloned with `--depth 1` at the commit listed below, outside this working tree.
- Evidence came from reading files and from searches: file counts, `LICENSE` files, test directories, CI
  workflows, and patterns such as `shell=True`, `exec(`, `eval(`, `pickle.load`, and API-key-shaped strings.
- Spot checks: for collections of hundreds of skills, one or two scripts were read closely. That is a sample,
  not a full audit, and scores for large collections carry that uncertainty.
- Scores run from 0 (absent or harmful) to 5 (strong evidence of good practice). A dash means the dimension
  does not apply, for example scientific correctness for an icon library.

| Repository | Commit | Date |
| --- | --- | --- |
| Future-House/paper-qa | c692d0d | 2026-10-25 |
| SakanaAI/AI-Scientist-v2 | 96bd516 | 2025-12-19 |
| SakanaAI/AI-Scientist | 1de1dbc | 2025-12-19 |
| HKUDS/AI-Researcher | f9a6f84 | 2025-10-16 |
| zLanqing/codex-claude-academic-skills | 7ed6377 | 2026-05-14 |
| going-doer/Paper2Code | ba91699 | 2026-03-25 |
| deepchem/deepchem | 455d07f | 2026-08-20 |
| hussius/deeplearning-biology | 9da4fae | 2026-09-12 |
| ai4s-research/awesome-ai-for-science | d8d6460 | 2026-11-05 |
| aipoch/medical-research-skills | 686e09d | 2026-09-17 |
| duerrsimon/bioicons | d29e766 | 2026-06-10 |
| icebird1998/drawio-scientific-illustrator | 9bbeca9 | 2026-07-30 |
| showlab/Paper2Video | 47beb50 | 2026-03-05 |
| pyqtgraph/pyqtgraph | ae5c86c | 2026-11-05 |
| K-Dense-AI/scientific-agent-skills (reference) | a25cf9e | 2026-11-02 |
| synthetic-sciences/openscience (reference) | a08e0ad | 2026-10-30 |

## What each repository is

| Kind | Repositories |
| --- | --- |
| Agent skills (`SKILL.md`) | K-Dense (177 skills), openscience (372 bundled), aipoch (605 `SKILL.md` files), zLanqing (3 skills plus 21 bundled copies), drawio-scientific-illustrator (1 skill plus an MCP server) |
| Autonomous research systems | AI-Scientist, AI-Scientist-v2, AI-Researcher, Paper2Code, Paper2Video |
| Libraries | paper-qa (RAG over papers), deepchem (molecular ML), pyqtgraph (Qt plotting) |
| Lists | awesome-ai-for-science, deeplearning-biology |
| Assets | bioicons (scientific icons) |

## Scores

Columns: Sci = scientific correctness, Ver = verification and tests, Bad = handling of bad input,
Skill = usability as an agent skill, Lic = license, Sec = security.

| Repository | Sci | Ver | Bad | Skill | Lic | Sec |
| --- | --- | --- | --- | --- | --- | --- |
| K-Dense scientific-agent-skills | 4 | 4 | 4 | 4 | 5 | 4 |
| synthetic-sciences/openscience | 3 | 3 | 3 | 3 | 5 | 4 |
| aipoch/medical-research-skills | 2 | 1 | 2 | 3 | 2 | 2 |
| zLanqing/codex-claude-academic-skills | 2 | 0 | — | 3 | 1 | 3 |
| icebird1998/drawio-scientific-illustrator | — | 2 | — | 3 | 5 | 3 |
| Future-House/paper-qa | 4 | 4 | 3 | 2 | 5 | 4 |
| deepchem/deepchem | 4 | 5 | 3 | 1 | 5 | 3 |
| pyqtgraph/pyqtgraph | — | 4 | — | 0 | 5 | 4 |
| going-doer/Paper2Code | 2 | 0 | 1 | 1 | 5 | 3 |
| showlab/Paper2Video | — | 0 | 1 | 1 | 5 | 2 |
| HKUDS/AI-Researcher | 2 | 0 | 1 | 1 | 2 | 1 |
| SakanaAI/AI-Scientist | — | 0 | — | 1 | 1 | 2 |
| SakanaAI/AI-Scientist-v2 | — | 0 | — | 1 | 1 | 2 |
| ai4s-research/awesome-ai-for-science | — | 1 | — | 0 | 5 | — |
| hussius/deeplearning-biology | — | 0 | — | 0 | 0 | — |
| duerrsimon/bioicons | — | 1 | — | 1 | 4 | 4 |

## Evidence by repository

### K-Dense-AI/scientific-agent-skills (reference)

- **Strongest:** tests use independent oracles. `tests/statistical-power/test_scripts.py` has 51 tests, compares
  results with `statsmodels` (`FTestPowerF2`) and closed-form normal values, and states solver tolerances. There
  are 224 test files and a `skill-tests.yml` workflow with a repository-wide contract and coverage guard.
- Bad input: `skills/peer-review/scripts/_common.py` checks input paths, rejects duplicate JSON keys, and
  validates enums, dates, DOIs, and URLs. The `statistical-analysis` skill documents that it rejects
  infinities and degenerate samples.
- Security: hand-written `SECURITY.md` with private vulnerability reporting, plus `security-scan.yml` and
  `pr-skill-scan.yml`, which scan skills with an LLM-based scanner (`scan_skills.py`). The 21 files with
  `exec(` are all under `tests/`, where they run documented examples. The `sk-...` strings are obvious test
  placeholders.
- License: MIT.
- **Weakest:** breadth over depth. 158 per-skill test directories cover 177 skills, so some skills have none. The
  `peer-review` skill's `audit_statistics_reproducibility.py` checks that a reviewer's checklist is complete
  and well formed. It does not recompute any reported statistic.

### synthetic-sciences/openscience (reference)

- **Strongest:** licensing hygiene. `backend/cli/skills/ATTRIBUTION.md` records the source, license, and degree
  of change for each imported skill (166 from K-Dense), and each `SKILL.md` repeats it under `metadata.upstream*`.
  The repository is Apache-2.0 with a `NOTICE` file. CI includes CodeQL, gitleaks, and OpenSSF Scorecard.
- Its `core/peer-review` skill asks the reviewer to show any numbers they recomputed and to list unchecked
  items as unchecked. That is good guidance, but it has no script that does the recomputation.
- **Weakest:** few skill-level tests. Only 10 test files sit under `skills/` against 372 skills, so most skill
  scripts are untested. Frontmatter names runtime-specific tools (`literature`, `task`, `webfetch`), which
  limits portability outside OpenScience.

### aipoch/medical-research-skills

- **Weakest, scientific correctness:** in
  `scientific-skills/Data Analysis/statistical-analysis-advisor/scripts/main.py`, `calculate_power()` accepts
  an `alpha` argument but hard-codes `z_alpha = 1.96`. For any power other than 0.8, 0.9, or 0.95 it silently
  falls back to 0.84 (`.get(power, 0.84)`). It also truncates the sample size with `int()` rather than
  rounding up. For d = 0.5 at 80% power that gives 62 per group, where the formula gives 62.7, so 63.
- Verification: one regression script (`tests/publication_regressions.py`) covers three named skills. The 279
  `eval_report_*.json` files are rubric scores from the authors' own `skill-auditor`. Under "code usability"
  they record checks such as "Scripts parse without errors", which is not a test of results.
- License: MIT. The README says all skills are "originally designed and developed in-house", yet 19
  `SKILL.md` files mention K-Dense. Whether those are derived could not be settled from the files, so
  attribution is unclear.
- Security: `scripts/openclaw-install.sh`, two files with `shell=True`, and two with `os.system`.
- **Strongest:** breadth of medical-research topics, and some skills, such as `paper-to-claim-verifier`, state
  explicitly that references must never be fabricated.

### zLanqing/codex-claude-academic-skills

- Three skills: Chinese-first research writing, office documents, and a scientific toolkit. The toolkit
  bundles 21 skill folders under `references/scientific-skills/`, most marked `skill-author: K-Dense Inc.`
- **Weakest, license:** the bundled `pdf/` and `xlsx/` folders carry `LICENSE.txt` reading
  "© 2025 Anthropic, PBC. All rights reserved", while the repository is MIT. Proprietary files are being
  redistributed inside an MIT repository. No tests, no CI.
- **Strongest:** clear routing between its three skills in each `SKILL.md`.

### icebird1998/drawio-scientific-illustrator

- One 45-line skill plus a Node MCP server that drives the draw.io desktop app. `npm test` runs syntax checks
  and a smoke test in CI.
- Security: `PRIVACY.md` states that the server binds only to 127.0.0.1 and has no telemetry. However,
  `live-server.mjs` launches draw.io with `--remote-debugging-port`. While it runs, any local process can
  attach to and control that app instance. It is local-only, but worth knowing.
- **Strongest:** honest privacy statement, MIT. **Weakest:** no tests of the drawings produced.

### Future-House/paper-qa

- **Strongest:** a serious, tested library. It has 14 test modules with 47 recorded HTTP cassettes, so tests
  replay network calls offline, and CI in `tests.yml`. Answers are grounded in cited passages.
- Security: `agents/search.py` uses `pickle.loads` on its own index files (marked `noqa: S301`). That is safe
  for indexes you built yourself, but do not load an index from an untrusted source.
- **Weakest as a skill:** it is a Python package and CLI that needs LLM and embedding API keys. It retrieves and
  summarizes text. It does not recompute numbers.

### deepchem/deepchem

- **Strongest:** 243 in-package test files and 12 CI workflows across TensorFlow, PyTorch, and JAX back ends.
  MIT.
- Security: `pickle.load` in 8 files (dataset caching, normal for ML libraries), and an AWS-access-key-shaped
  string in a tutorial notebook's stored output (`examples/tutorials/Multisequence_Alignments.ipynb`). It is
  probably a signed URL fragment. It was not tested and is not reproduced here.
- **Weakest as a skill:** a heavy library, not a skill. K-Dense already ships a `deepchem` skill that wraps it.

### pyqtgraph/pyqtgraph

- Interactive Qt plotting library with 72 test files and CI. MIT (UNC copyright).
- **Weakest:** not relevant to agent skills, since it needs a GUI event loop. **Strongest:** a mature,
  well-tested library.

### going-doer/Paper2Code

- Multi-agent LLM pipeline from paper to code repository (ICLR 2026). Apache-2.0.
- **Weakest:** no tests and no CI. Ten `except Exception` blocks in `codes/`. Generated code is checked only by
  later LLM stages. **Strongest:** clear staged pipeline and permissive license.

### showlab/Paper2Video

- Generates presentation videos from papers. MIT. Needs Gemini or OpenAI keys (README).
- Security: `src/evaluation/PresentQuiz/utils/wei_utils.py` calls `exec(code, ...)`, and `apis.py` calls `eval`
  on model output lines. **Weakest:** no tests. **Strongest:** permissive license and a documented pipeline.

### HKUDS/AI-Researcher

- Autonomous research agent run in Docker.
- **Weakest, security:** API-key-shaped strings are committed in `README.md`,
  `benchmark/process/dataset_candidate/reasoning/math_reasoning/.env`, and `.../reasoning/metaprompt.py`. They
  were not tested and are not reproduced here. There are also four `shell=True` calls, one of which runs
  received commands in `research_agent/inno/environment/tcp_server.py`.
- License: MIT is declared only in `setup.cfg`. There is no `LICENSE` file. No tests, and 190 broad `except`
  blocks. **Strongest:** Docker isolation is the documented default.

### SakanaAI/AI-Scientist and AI-Scientist-v2

- License: "The AI Scientist Source Code License" (based on the Responsible AI License). It restricts use
  (surveillance, health care, and others) and requires a prominent disclosure on any manuscript produced
  (clause 3.2.e). That is not an open-source license. **The owner excluded both repositories, and nothing from
  them was used.**
- Security: v2 executes LLM-written code with `exec(compile(...))` in `ai_scientist/treesearch/interpreter.py`.
  The README warns about this and asks for a sandbox. No tests in either repository.

### ai4s-research/awesome-ai-for-science

- A curated list of about 785 linked entries. MIT. No CI or link checking, so link rot is unchecked. Useful
  for discovery, not as a skill.

### hussius/deeplearning-biology

- A single 72 KB README list. It has **no license file**, so all rights are reserved by default. The owner
  excluded it, and nothing from it was used.

### duerrsimon/bioicons

- The web app is MIT. Icons sit under `static/icons/<license>/` (`cc-0`, `cc-by-3.0`, `cc-by-4.0`,
  `cc-by-sa-3.0`, `cc-by-sa-4.0`, `bsd`, and others), and the README asks users to cite each icon. CI only lints
  the web app.
- **Strongest:** per-icon licensing is explicit. **Weakest:** a skill using it would have to carry each icon's
  attribution and share-alike terms into figures.

## Is the paper-audit gap real?

Yes, as far as these repositories go. Searches across all 16 clones for `GRIM`, `statcheck`, `GRIMMER`,
`SPRITE`, "recompute", "recalculate", and "reported p-value" found no skill or script that recomputes the
numbers printed in a paper. The `SPRITE` hits were sprite graphics and an OpenGL constant. The closest
neighbours:

| Skill | What it does | Recomputes printed numbers? |
| --- | --- | --- |
| K-Dense `peer-review` (`audit_statistics_reproducibility.py`) | Checks a reviewer's 22-item statistics checklist for completeness | No |
| openscience `core/peer-review` | Prose guidance; asks the reviewer to show recomputed numbers | No script |
| openscience `core/acceptance-checks` | Turns a task's stated criteria into tests of the agent's own outputs | No; it is about the agent's work, not a paper |
| aipoch `consistency-checker-across-manuscript` | LLM prompt and reference rules on text and numeric drift between sections | No script |
| aipoch `paper-to-claim-verifier` | Checks whether cited papers support a claim | No; it works on text |
| K-Dense `scholar-evaluation` | Scores a paper against a rubric | No |

So `skills/paper-audit` adds something these lack: deterministic recomputation of p-values, GRIM, percentages,
totals, and interval agreement, with rounding intervals, a traceable claim ledger, and tests against
independent values. Outside these repositories, the R package statcheck and the published GRIM procedure do
similar things. They are cited as methods in the skill's references, and no code was taken from either.

## Lessons taken into this repository

- From K-Dense: test against independent oracles (scipy, exact fractions), not against the code's own output.
- From openscience: record provenance. paper-audit is original and says so; its sources are method papers.
- From aipoch's power function: never accept a parameter and ignore it. paper-audit rejects inputs it cannot
  use and reports them as not checkable instead.
