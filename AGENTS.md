# AGENTS

## Project Purpose
Pulmonary artery enlargement and mortality in an emergency department CTPA population

## Public and Data-Safety Rules
- Treat this repository as public. Do not add PHI, restricted datasets, credentials, private drafts, or publisher-formatted article text.
- Clinical imaging/EHR data likely restricted; verify no PHI
- Manuscript status: Local manuscript candidates exist in PA-size project; publisher-policy check needed before Markdown

## How to Orient Quickly
- Start with `README.md` for project scope, workflow, data notes, citation, and license information.
- Use `CITATION.cff` for structured citation metadata when present.
- Inspect scripts/notebooks before running them; do not assume generated outputs are current.

## Workflow
From the repository root, use this as the initial run guidance:

```bash
stata-mp -b do "Age-Metric Analysis.do"
```

If the command is a placeholder, refine it after reading the local scripts and existing README.

## Verification Before Publishing Changes
- Run `git diff --check`.
- Validate `CITATION.cff` as YAML after citation edits.
- Do not commit generated outputs, logs, caches, virtual environments, `.DS_Store`, or checkpoint files unless intentionally released.
- For clinical or collaborator data, confirm that no row-level restricted data or identifiers are included.
