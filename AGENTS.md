# AGENTS

## Project Purpose

This is a public repository for an unpublished CTPA dose-response analysis presented as a CHEST 2023 abstract. The repository contains legacy Stata workflows, exploratory notebooks, a labeled unpublished author draft, and machine-readable documentation.

Primary scholarly citation: CHEST 2023 abstract DOI `10.1016/j.chest.2023.07.3871`.

Related background paper, not the primary citation target: Pulmonary Circulation DOI `10.1002/pul2.12225`.

## Data-Safety Rules

- Do not add PHI, patient-level datasets, row-level derived `.dta` files, local registry files, credentials, or collaborator-only files.
- Treat `final_noempi.dta`, `cleaned_noempi.dta`, and `ACCT EMPI BMI Weight Height.csv` as restricted clinical data even though direct EMPI identifiers were removed.
- Do not commit publisher PDFs or third-party reference PDFs. Link DOI or publisher pages instead.
- Generated notebook HTML/PDF exports belong in releases, not in the source tree.
- Keep the manuscript Markdown labeled as unpublished author draft text. Do not describe it as peer reviewed, accepted, or journal-published.

## Workflow

Run Stata commands from the repository root. Default input and output paths are `data/private` and `outputs/stata`.

```bash
stata-mp -b do "Data_preprocess.do" "data/private" "outputs/stata"
stata-mp -b do "Age-Metric Analysis.do" "data/private" "outputs/stata"
stata-mp -b do "PA Body-size Analysis.do" "data/private" "outputs/stata"
```

The older `BL CombofilewithoutEMPI.do` script is retained for provenance and can be run with the same two optional arguments. Do not refactor scientific modeling logic unless explicitly requested.

## Documentation Checks

- Keep `README.md`, `llms.txt`, `CITATION.cff`, and the data dictionaries aligned on study status, DOI, run order, and data restrictions.
- Update both `data_dictionary.md` and `data_dictionary.csv` when variable definitions or expected inputs change.
- Validate `CITATION.cff` after edits.
- Run `git diff --check` before publishing.

## Artifact Hygiene

- Expected generated output path: `outputs/stata/<date>/`.
- Expected restricted input path: `data/private/`.
- Do not commit `.DS_Store`, Stata logs, `.gph`, `.dta`, generated `.xlsx`, generated HTML/PDF notebook exports, caches, or virtual environments.
- If preserving large historical generated artifacts, use a GitHub release and label them as historical/provenance artifacts.
