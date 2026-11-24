# discrete-time-competing-risks

Code for the manuscript **"A partial-likelihood approach for the analysis of discrete time competing risks data"** (Seohyun Lee, Minjung Lee, Sangbum Choi).

The manuscript proposes a Poisson-approximation-based partial likelihood for cause-specific hazards with discrete-time competing risks data. Regression coefficients are estimated by a Cox-type partial likelihood with cluster-robust sandwich variance estimation. The baseline hazards are estimated by a Breslow-type estimator, and a trend-filtering penalty smooths them using ADMM. The methods are compared with Lee, Feuer and Fine (2018, *Biometrics*) and Meir and Gorfine (2025, *Biometrics*).

## Repository structure

```
.
├── manuscript/          # PDF of the manuscript
└── code/
    ├── R/               # Core functions (partial likelihood, ADMM trend filtering)
    ├── sim_4.1_PL/      # Section 4.1: simulation for regression coefficients and CIF prediction
    ├── sim_4.2_TF/      # Section 4.2: simulation for trend-filtered baseline hazards
    └── seer/            # Section 5: SEER colorectal cancer application
```

| Folder | Content |
|---|---|
| `code/R` | `admm_tf_long.R`: person-period data utilities and the original ADMM implementation. `admm_bcd_exact.R`: final ADMM solver (exact block coordinate descent primal step) and the regularization path function `tf_path_bcd()`. |
| `code/sim_4.1_PL` | Comparison of the proposed estimator, Lee et al. (2018) and Meir and Gorfine (2025) under two data-generating processes. Reports bias, average and empirical standard errors, and coverage. Also includes the cumulative incidence prediction simulation. |
| `code/sim_4.2_TF` | Baseline hazard estimation with trend filtering under smooth, step and abrupt shapes, for the log and logit links. |
| `code/seer` | SEER application: model fits, baseline hazard and cumulative incidence figures, and bootstrap confidence intervals. |

## Requirements

- R 4.6.0
- Packages: `survival`, `dplyr`, `ggplot2`, `CVXR`

## How to run

All scripts are run from the `code/` directory, since file paths are relative to it.

```bash
cd code

# Section 4.2 simulation (four independent R processes)
sh sim_4.2_TF/run_sim_tf3.sh

# Section 5: SEER application
sh seer/run_seer_tf3.sh
sh seer/run_seer_cif.sh
```

Simulation settings (sample size, number of replications, data-generating process) are set at the top of each script or through environment variables described in the script header. Setting a small number of replications (for example `B=2`) is useful to check that everything runs before a full run.

## SEER data

The SEER data are not included in this repository because of the SEER data use agreement. Access can be requested from the National Cancer Institute (https://seer.cancer.gov/data/), and case listings are extracted with SEER*Stat from the SEER Research Data (12 registries, November 2025 submission).

The following restrictions are applied in SEER*Stat:

- Site recode ICD-O-3/WHO 2008: colon excluding rectum
- AJCC stage III (3rd edition), year of diagnosis 1994–2003
- First primary tumor (sequence number 00 or 01)
- Cancer-directed surgery performed
- Reporting source other than autopsy or death certificate only

The case listing must contain the following variables: Patient ID, Survival months, Survival months flag, SEER cause-specific death classification, SEER other cause of death classification, Sex, Marital status at diagnosis, Site recode ICD-O-3/WHO 2008, Year of diagnosis, Age recode with single ages and 90+, Grade Recode (thru 2017), and Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic).

Export the listing as `export.csv` and run `Rscript seer_prep_base.R`. The script applies the remaining restrictions:

- Age 66 or older
- Survival time computed from complete dates
- Classifiable cause of death
- Appendiceal tumors excluded
