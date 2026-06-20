# CTPA Dose-Response Analysis

[![DOI: CHEST abstract](https://img.shields.io/badge/DOI-10.1016%2Fj.chest.2023.07.3871-blue)](https://doi.org/10.1016/j.chest.2023.07.3871)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![CITATION.cff](https://img.shields.io/badge/citation-CFF%201.2-informational)](CITATION.cff)
[![LLM index](https://img.shields.io/badge/LLM-index-llms.txt-informational)](llms.txt)

Public repository for the unpublished full manuscript and legacy analysis code for **Pulmonary Artery Size on CT and Mortality Risk: A Dose-Response Analysis**, presented as a CHEST 2023 abstract.

The work evaluates whether pulmonary artery diameter (PAd) and pulmonary artery to ascending aorta ratio (PA:AA) have dose-response associations with mortality after emergency-department computed tomography pulmonary angiography (CTPA), and whether age/sex normalization changes prognostic interpretation.

## Project Status and Links

This repository is associated with a conference abstract and an unpublished author manuscript draft. There is no final peer-reviewed journal article for this dose-response analysis.

- CHEST 2023 abstract: [Pulmonary Artery Size on CT and Mortality Risk: A Dose-Response Analysis](https://doi.org/10.1016/j.chest.2023.07.3871), CHEST 2023 annual meeting abstracts, volume 164, issue 4 supplement, pages A6004-A6005.
- Public CV/abstract page: [Brian W. Locke abstracts and presentations](https://reblocke.github.io/talks/).
- Unpublished author draft: [manuscript/ctpa-dose-response-unpublished-author-draft.md](manuscript/ctpa-dose-response-unpublished-author-draft.md).
- Source DOCX for provenance: [manuscript/source/2023-10-10-dose-response-letter.docx](manuscript/source/2023-10-10-dose-response-letter.docx).
- Legacy rendered notebook exports: [legacy-rendered-artifacts-2026-06-04](https://github.com/reblocke/CTPA-Dose-Response/releases/tag/legacy-rendered-artifacts-2026-06-04).

Related background paper, not the primary citation target for this repository:

- Scarpato BM, Locke BW, Bledsoe J, et al. The association between pulmonary artery enlargement and mortality in an Emergency Department population undergoing computed tomography pulmonary angiography. *Pulmonary Circulation*. 2023;13(2):e12225. DOI: [10.1002/pul2.12225](https://doi.org/10.1002/pul2.12225).

## Authors, Funding, and Disclosures

Authors listed on the CHEST abstract and manuscript draft: Brian W. Locke, Brittany M. Scarpato, Joseph Bledsoe, Daniel B. Knox, Karen Conner, Gregory J. Stoddard, Meghan M. Cirulis, C. Gregory Elliott, and Mark W. Dodson.

Maintainer: Brian W. Locke, ORCID [0000-0002-3588-5238](https://orcid.org/0000-0002-3588-5238).

The draft states that Brian W. Locke was supported by NIH Ruth L. Kirschstein National Research Service Award `5T32HL105321`. It also notes unrelated funding through an American Thoracic Society program supported by ResMed, Philips Respironics, and Fisher & Paykel Healthcare.

## Data Access and Privacy

The analytic datasets are restricted clinical imaging/EHR-derived files and are not public. Do not commit patient-level data, row-level derived datasets, identifiers, local registry files, or collaborator-only workbooks.

Expected local-only inputs:

- `data/private/final_noempi.dta`
- `data/private/cleaned_noempi.dta`
- `data/private/ACCT EMPI BMI Weight Height.csv`

The `noempi` names indicate prior removal of direct EMPI identifiers, but these files still represent patient-level clinical data and must remain restricted.

## Repository Layout

| Path | Purpose |
| --- | --- |
| `Data_preprocess.do` | Cleans restricted source data and merges supplemental age/BMI/body-size fields. |
| `Age-Metric Analysis.do` | Main age/sex normalization and survival-analysis workflow for the dose-response manuscript. |
| `PA Body-size Analysis.do` | Exploratory analysis of PA, AA, body size, ratios, and mortality prediction. |
| `BL CombofilewithoutEMPI.do` | Older combined preprocessing/analysis script retained for legacy provenance. |
| `Unit Tests.do`, `calc_ibs.py`, `define_calc_ibs.do` | Experimental integrated Brier score test/prototype materials. |
| `age-sex-norm-letter.ipynb`, `ratio_val_analysis.ipynb` | Historical exploratory notebooks with cleared outputs. |
| `manuscript/` | Unpublished author manuscript draft and source DOCX. |
| `data_dictionary.md`, `data_dictionary.csv` | Machine-readable and human-readable variable documentation. |
| `llms.txt`, `AGENTS.md`, `CITATION.cff` | Machine-readable repository orientation and citation metadata. |

Third-party reference PDFs and generated HTML/PDF notebook exports are intentionally not tracked on the branch tip. Use DOI links and the legacy rendered-artifacts release instead.

## Workflow

Run Stata from the repository root. Full execution requires local restricted inputs.

```bash
stata-mp -b do "Data_preprocess.do" "data/private" "outputs/stata"
stata-mp -b do "Age-Metric Analysis.do" "data/private" "outputs/stata"
stata-mp -b do "PA Body-size Analysis.do" "data/private" "outputs/stata"
```

Default arguments are `data/private` for inputs and `outputs/stata` for generated outputs. The older combined workflow can be run similarly:

```bash
stata-mp -b do "BL CombofilewithoutEMPI.do" "data/private" "outputs/stata"
```

Generated logs, copied do-files, tables, graphs, and derived `.dta` files are written under `outputs/stata/<date>/` and ignored by git.

## Dependencies

| Component | Requirement |
| --- | --- |
| Stata | Tested/developed with Stata 17-18 syntax. |
| Stata packages | `missings`, `mdesc`, `table1_mc`, `coefplot`, `outreg2`, `mkspline2`, `xblc`, `stcoxkm`, `cleanplots`, `schemepack`; some exploratory commands may also need `palettes`, `white_tableau`, or related schemes. |
| Python notebooks | Install `requirements.txt` for notebook support and the integrated Brier score prototype. |
| Restricted data | Not included; supply compatible local files under `data/private/`. |

## Outputs

Expected generated outputs include:

- cleaned or merged local `.dta` files;
- Stata logs and copied do-files;
- table exports such as age-stratified baseline characteristics;
- survival curves, spline plots, quantile-regression plots, and PA/AA/PA:AA figures.

Historical rendered notebook exports are archived in the GitHub release named `legacy-rendered-artifacts-2026-06-04`. They are retained for provenance and reading convenience, not as canonical reproducible outputs.

## Citation

If citing the scholarly work, cite the CHEST abstract:

> Locke BW, Scarpato B, Bledsoe J, Knox D, Conner K, Stoddard GJ, Cirulis MM, Elliott CG, Dodson MW. Pulmonary Artery Size on CT and Mortality Risk: A Dose-Response Analysis. *CHEST*. 2023;164(4 Suppl):A6004-A6005. doi:10.1016/j.chest.2023.07.3871.

If citing repository code or draft materials, cite the repository URL, commit SHA, and `CITATION.cff`.

## License

Repository code and original documentation are released under the [MIT License](LICENSE). Restricted clinical data, third-party reference papers, publisher-formatted PDFs, and external datasets are excluded. The unpublished author draft is included for scholarly indexing and review; verify authorship, permissions, and journal/preprint policy before redistributing outside the repository context.

## Contact

Use GitHub issues or pull requests for public repository questions. Do not post PHI, patient-level data, local registry files, or collaborator-only materials in issues or pull requests.
