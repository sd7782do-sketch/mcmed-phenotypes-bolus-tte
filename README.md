# Hemodynamic phenotypes and crystalloid boluses in infection-associated hypotension (MC-MED)

Analysis code for the study *Hemodynamic phenotypes and crystalloid boluses in emergency department patients with infection-associated hypotension: a sequential target trial emulation*.

Registration: Open Science Framework, <https://doi.org/10.17605/OSF.IO/EWQDK> (with updates documenting post-registration changes).

## Study in brief

- **Data**: Multimodal Clinical Monitoring in the Emergency Department (MC-MED), PhysioNet, version 1.0.1.
- **Design**: sequential target trial emulation. A new trial starts at the first hypotensive cuff reading and every 15 minutes thereafter, up to 3 hours, in emergency department visits with infection. Eligibility and confounders use only information available at trial start.
- **Exposure**: crystalloid bolus started within 15 minutes of trial start (grace period) versus not.
- **Effect modifier**: hemodynamic phenotype on the trial-start reading:
  - low pulse pressure: diastolic > 50 and pulse pressure < 40 mm Hg;
  - vasodilatory: diastolic ≤ 50 and pulse pressure ≥ 40 mm Hg;
  - mixed: both low.

  Diastolic pressure, pulse pressure and the diastolic shock index are also analyzed as continuous modifiers.
- **Outcomes**:
  - primary: mean arterial pressure at 60 minutes;
  - secondary: perfusion index at 60 minutes, and norepinephrine within 6 hours.
- **Analysis**:
  - stabilized inverse probability weights, weighted outcome models with a bolus-by-modifier interaction, and g-computation;
  - cluster bootstrap resampling patients;
  - Bayesian reanalysis under skeptical, weakly informative and enthusiastic priors;
  - 11 sensitivity analyses, including clone-censor weighting of the grace period and phenotype-stratified propensity models.

## Data access

MC-MED is a credentialed-access PhysioNet database (<https://doi.org/10.13026/wvyw-g663>). **No data are included in this repository.** To run the code, you must obtain access yourself and download the database. Under the PhysioNet data use agreement, record-level data must not be redistributed, and that includes the intermediate file `07_analytic_dataset.csv`.

The scripts read `visits.csv`, `orders.csv`, `labs.csv` and `numerics.csv`.

## Repository structure

| File | Purpose |
|---|---|
| `R/01_analysis_seqtrials_mcmed.R` | Full analysis: cohort, emulated trials, weights, outcome models, bootstrap, Bayesian and sensitivity analyses (1–10), phenotype transitions, trajectories. Writes all result tables (`*.csv`) and most figures. |
| `R/02_figure2_effect_distributions.R` | Figure 2 (distribution of bootstrap estimates by phenotype); reads outputs of script 01. |
| `R/03_figureE1_flowchart.R` | Figure E1 (flow of visits and emulated trials); reads outputs of script 01. |
| `R/04_race_ethnicity_sex.R` | Sex, race and ethnicity for Tables 1 and 2; reads outputs of script 01 and `visits.csv`. |
| `R/legacy/analysis_as_registered.R` | Analysis exactly as originally registered, kept unchanged for transparency; used only for sensitivity analysis 11. |

## How to run

The scripts are designed for R 4.2 or later.

```r
install.packages(c("data.table", "arrow", "dplyr", "lubridate", "rlang",
                   "ggplot2", "ggalluvial", "ragg", "systemfonts"))
# splines and grid ship with R
```

Paths are set by environment variables. Without them, the scripts fall back to `~/Desktop/mimiciv/...`.

```r
Sys.setenv(MCMED_DATA_DIR = "/path/to/mc-med/1.0.1",        # folder with the MC-MED csv files
           MCMED_OUT_DIR  = "/path/to/output")               # results folder
source("R/01_analysis_seqtrials_mcmed.R")   # long run: 300 + 200 bootstrap replicates
source("R/02_figure2_effect_distributions.R")
source("R/03_figureE1_flowchart.R")
source("R/04_race_ethnicity_sex.R")

# optional, sensitivity analysis 11 (writes to MCMED_OUT_DIR_REGISTERED)
source("R/legacy/analysis_as_registered.R")
```

Further notes:

- **Seed**: the seed is fixed in the scripts (`set.seed(20261008)`).
- **Font**: Figure 2 uses the Lato font when it is installed, with a fallback to Arial, Helvetica or Liberation Sans. Run script 02 in a UTF-8 locale so that minus signs render.
- **Figure format**: figures are saved as PNG, PDF and TIFF (LZW compression, 600 dpi) and contain no titles; the legends are in the manuscript.

## Correspondence between outputs and the manuscript

| Manuscript | Output |
|---|---|
| Table 1 / Table 2 | `02_table1_by_exposure.csv`, `03_table_by_phenotype.csv`, `20_race_ethnicity_sex.csv` |
| Figure 1 | `fig4_alluvial_phenotype_transitions.*` (counts in `10_`, `11_`) |
| Figure 2 | `fig2_effect_distributions_map.*` (script 02; `05_effects.csv`, `08_bootstrap_draws_primary.csv`) |
| Figure 3 | `fig5_trajectories_map_perf.*` (`12_`, `13_`) |
| Figure 4 | `fig2_effect_by_dap_pp_dsi.*` (`06_effect_curves.csv`) |
| Table E2 | `15_bolus_volume_by_phenotype.csv` |
| Table E3 | `16_weight_diagnostics.csv` |
| Table E4 | `04_balance.csv` |
| Table E5 | `17_balance_by_phenotype.csv`, `17b_balance_by_phenotype_alt_weights.csv` |
| Table E6 | `18_observed_n_by_outcome.csv`, `19_norepinephrine_timing.csv`, `14_outcome_ascertainment_by_arm.csv` |
| Table E7 | `09_bayesian_posteriors.csv`, `09b_bayesian_large_boluses_map.csv` |
| Table E8 | `05_effects.csv` (sensitivity analyses 1–10) and `05_effects.csv` from the legacy script (analysis 11) |
| Figure E1 | `figS1_flowchart.*` (script 03; `01_flow.csv`, `18_observed_n_by_outcome.csv`) |
| Figure E2 | `fig1_effect_by_phenotype.*` |
| Figure E3 | `figA_map_response_by_baseline.*` |
| Figure E4 | `fig6_effect_by_prior_boluses.*` |
| Figure E5 | `fig3_bayesian_posteriors_map.*` |
| Figure E6 | `figB_probability_of_benefit.*` |

Table E1 (the target trial specification) is text only and has no corresponding output. The Word tables were formatted from these CSV files.

## License

Code released under the MIT License (see `LICENSE`).
