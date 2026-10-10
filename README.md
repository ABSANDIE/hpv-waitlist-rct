# HPV community follow-up trial: analysis code

Analysis code for the randomized controlled trial **"Improving HPV Immunization Coverage Among
Adolescent Girls in Mali and Côte d'Ivoire Through an Innovative Community-Based Follow-Up"**
(Muso Health; principal investigator Arsène Brunelle Sandie).

- Registration: [ClinicalTrials.gov NCT07807202](https://clinicaltrials.gov/study/NCT07807202)
- Statistical analysis plan (registered): [`sap/SAP_v2.1_2026-10-07.pdf`](sap/SAP_v2.1_2026-10-07.pdf), version 2.1, 7 October 2026, uploaded to the ClinicalTrials.gov record before the outcome data were accessed
- Later revision (not registered): [`sap/SAP_v2.5_2026-10-09.pdf`](sap/SAP_v2.5_2026-10-09.pdf), version 2.5, 9 October 2026, which maps every analysis variable to the CHT fields

This repository contains **code and synthetic data only**. It contains no trial data.

## The trial in brief

| | |
|---|---|
| Design | Pragmatic, individually randomized trial with a waitlist crossover |
| Participants | Adolescent girls eligible for HPV vaccination: age 10 in Mali, 9–14 in Côte d'Ivoire |
| Intervention | Community health worker (CHW) follow-up visits supported by a digital package |
| Control | Usual care, then the intervention from 2 July 2026 (waitlist) |
| Launch | 8 April 2026 (Mali), 9 April 2026 (Côte d'Ivoire) |
| Primary outcome | HPV vaccination status on 1 July 2026 (12 weeks), verified up to 1 October 2026 |
| Secondary outcome | Time to HPV vaccination within 12 weeks (difference in restricted mean survival time) |
| Phase 2 | Vaccination status on 1 October 2026: early vs delayed start of the intervention |
| Primary model | Mixed-effects logistic regression: arm + country + age (spline) + random intercept for CHW supervision zone |

Girls expected to age out of eligibility within the first three months were allocated to the
intervention arm without randomization. The primary analysis includes them, and every analysis is
repeated in randomized girls only (SAP, "Girls Allocated by the Age-Out Rule").

## What the synthetic data is, and what it is not

`R/00_simulate_test_data.R` generates fake tables with the **same structure and field names** as the
four Community Health Toolkit (CHT) extracts used in the trial: patients, baseline assessment,
follow-up visits and verification round.

- **It is** a test bed: it reproduces the data problems the code must handle. These include duplicate forms,
  doses reported without a date, deactivated records, girls who age out, refusals and missing visits.
- **It is not** trial data. The sample is about 15% of the real size, and vaccination rates and effect sizes are
  invented, with a built-in intervention effect. **No result produced from it says anything about the trial.**
- No names, phone numbers or real identifiers appear anywhere in this repository.

## How to reproduce the outputs

Requirements: R ≥ 4.3 and these packages:

```r
install.packages(c("dplyr", "tidyr", "readr", "stringr", "purrr", "lubridate", "knitr",
                   "rmarkdown", "ggplot2", "lme4", "survival", "geepack", "mice",
                   "marginaleffects", "MASS", "scales", "readxl"))
```

The `survRM2` package is optional; if installed, it adds the Uno et al. RMST cross-check.

Then, from the repository root:

```sh
Rscript run_all.R quick   # about 5 minutes: 5 imputations, 50 bootstrap replicates
Rscript run_all.R         # full specification: 50 imputations, 500 bootstrap replicates
```

`run_all.R` generates the synthetic data (fixed seed) if it is missing, then renders the two reports.
Everything lands in `outputs/`:

| File | Content |
|---|---|
| `01_prepare_analysis_dataset.html` | Derivation of the girl-level dataset and data-quality checks (no outcomes by arm) |
| `02_prespecified_analysis.html` | All pre-specified analyses (SAP Tables 1–8, Figures 3–6) |
| `table*.csv`, `sens_*.csv`, `fig*.csv` | The same tables as CSV |
| `fig5a_*.png`, `fig5b_*.png` | Subgroup forest plots (coverage and uptake at 12 weeks) |

## Repository layout

```
R/hpv_functions.R            derivations from CHT fields (SAP v2.5, "Derivation from CHT Fields")
R/hpv_analysis_functions.R   models: glmer, marginal standardization, KM, RMST, Cox, MI, bounds
R/00_simulate_test_data.R    synthetic data generator
analysis/01_prepare_analysis_dataset.Rmd
analysis/02_prespecified_analysis.Rmd
sap/                         statistical analysis plan
run_all.R                    one command to reproduce everything
config_real_template.R       template for running on the real extracts (the copy is git-ignored)
```

## Pre-specification and deviations

The **registered SAP is version 2.1** (7 October 2026). It is the version uploaded to ClinicalTrials.gov
before the outcome data were accessed. The git tag `sap-v2.1-registered` marks the commit that adds
this file; the tag message records its SHA-256 checksum, so anyone can check it is the same file as
the registered document. The tag was created when this repository was set up, after the analysis
had been run: the registry upload, not the tag, is the evidence of pre-specification.

Version 2.5 (9 October 2026) was written after version 2.1 was registered. It adds a mapping of the
analysis variables to the CHT fields and two rule changes. The code in this repository was first
written against version 2.5, so it differs from the registered plan in two places:

| Topic | Registered SAP v2.1 | Code as first committed (SAP v2.5) |
|---|---|---|
| Dose reported without a documented date | Unknown status, cause 2 | Vaccinated if reported at a visit on or before 1 July; unknown (cause 2) if reported later |
| Parents refused to give information | No separate category; counted as cause 2 (reached, status not verified) | Separate cause 4, placed before 3a in the hierarchy |

Two further points follow from what the CHT tables contain, and apply in every version of the code:

- **Eligibility review (cause 3b).** Version 2.1 bases it on the CHT audit trail of age and sex
  corrections. CHT keeps only current values, so the code applies the review rule to the current
  date of birth and sex (as described in version 2.5).
- **Age-out allocation.** The allocation list does not record which girls were allocated by the
  age-out rule. The code recovers them from the date of birth and checks the counts against the
  known totals (1,945 in Mali, 165 in Côte d'Ivoire).

Any change to the analysis after the first code commit is visible in the history of this
repository.

### Additions not in either SAP version

- Overall Kaplan-Meier coverage and overall RMST difference also standardized to the country mix,
  because the age-out allocation adds Mali girls to the intervention arm only.
- Phase 2 status recomputed from visits on or after 2 July only (sensitivity).
- Verification-lag model fitted among girls with a verification visit in the 2 July – 1 October
  window: girls known only from a dose dated before 1 July have no verification visit.
- Subgroups plotted only if each arm has at least 5 vaccinated girls (all estimates kept in the CSV).
- Figure 5b, subgroup ORs for uptake at 12 weeks: **post hoc**, added after the primary results
  were seen.

## Licence and contact

Code: MIT licence (see `LICENSE`). Arsène Brunelle Sandie, Muso Health.
