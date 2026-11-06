# Catalog notes

The catalog, `assets/tools.json`, summarizes `docs/repo-comparison.md` at the root of this repository. That
report holds the evidence: commit hashes, file paths, and scores for each repository. This file explains how
the catalog was built and how to keep it current.

## How each status was decided

- **recommended**: a permissive license (MIT or Apache-2.0, or per-item licenses that are stated), plus tests or
  a clear, low-risk use. Examples are paper-qa, deepchem, pyqtgraph, K-Dense, bioicons, the draw.io illustrator,
  the awesome list, and paper-audit.
- **caution**: a usable license but no tests, executes model-generated code, unclear attribution, or a whole
  agent rather than a drop-in tool. Examples are OpenScience, AIPOCH, Paper2Code, Paper2Video, AI-Researcher,
  and zLanqing.
- **excluded**: license terms the repository owner does not accept. These are the two SakanaAI repositories
  (the AI Scientist Source Code License restricts use and requires disclosure) and hussius/deeplearning-biology
  (no license file, so all rights are reserved).

## Where the install commands came from

Each `install` field was taken from that repository's own README at the reviewed commit. Examples are
`pip install paper-qa`, `pip install deepchem`, `pip install pyqtgraph`, and
`npx skills add K-Dense-AI/scientific-agent-skills`. For tools without a package, the field says to clone the
repository and follow its README. Where a README pipes a remote script into a shell, the catalog asks the user
to read the script first.

## Updating the catalog

1. Re-check the repository's license file, README install section, and tests.
2. Edit its entry in `assets/tools.json`. Keep every field, and keep an `alternative` for excluded tools.
3. Update `docs/repo-comparison.md` with the new evidence.
4. Run `python tests/run_all.py`. The tests check that every entry is complete, that excluded tools are never
   suggested, and that the catalog and the comparison report name the same repositories.
