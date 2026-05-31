Code for PA size to mortality analyses 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/reblocke/CTPA-Dose-Response/HEAD?urlpath=%2Fdoc%2Ftree%2Fage-sex-norm-letter.ipynb)

## LLM and Repository Readiness Notes

### Description
Pulmonary artery enlargement and mortality in an emergency department CTPA population

### Instructions
Start with this README, then inspect the files listed under Repository Layout. For computational workflows, run commands from the repository root and avoid committing generated outputs unless a release explicitly calls for them.

### Authors, Funding, and Acknowledgments
Maintainer: Brian W. Locke (`@reblocke`, ORCID 0000-0002-3588-5238). Preserve any project-specific author, funding, and acknowledgment details already listed elsewhere in the repository or accompanying publication.

### Repository Layout
- `.DS_Store`
- `2023-10-10 Dose Response Letter.docx`
- `Age-Metric Analysis.do`
- `BL CombofilewithoutEMPI.do`
- `Berger et al PA and AA by BSA nomograms.pdf`
- `CHEST 2025 - Prog Value of MPAD PAH.pdf`
- `Data_preprocess.do`
- `PA Body-size Analysis.do`
- `Pulm  circ - 2023 - Scarpato - The association between pulmonary artery enlargement and mortality in an Emergency (1).pdf`
- `README.md`
- `Unit Tests.do`
- `age-sex-norm-letter.html`
- `age-sex-norm-letter.ipynb`
- `age-sex-norm-letter.pdf`

### Data and Codebook
Clinical imaging/EHR data likely restricted; verify no PHI

### Workflow / Script Order
stata-mp -b do "Age-Metric Analysis.do"

### Dependencies / Environment
Repo README and scripts

### Citation
Preferred scholarly citation: https://doi.org/10.1002/pul2.12225. Cite this repository with the GitHub URL and the commit or release used.

### License
MIT License for repository code; see `LICENSE`. Third-party data, publisher text, and restricted clinical data are excluded.

### Manuscript Status
Local manuscript candidates exist in PA-size project; publisher-policy check needed before Markdown Publisher text not copied; license missing in first audit

### Contact
Maintainer: Brian W. Locke (`@reblocke`). Use GitHub issues or pull requests for repository-specific questions when the repository is public.
