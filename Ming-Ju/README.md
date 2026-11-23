# Temperature forecasting: model selection by cross-validation

This folder contains the code, results and report for forecasting the
daily maximum, mean and minimum temperature 1 to 5 days ahead at 15
target airports. The model for each horizon is a least-squares
regression on recent temperatures, the origin day's weather, airport
effects and annual harmonics. The history length `h` (0 to 15 days)
and the seasonal structure (4 options) are chosen jointly by 10-fold
cross-validation, one fold per calendar year from 2005 to 2014.

## Reports

There are two versions of the report. They differ only in the rule
used to flag temperature outliers.

| Report | Outlier rule | Results folder it reads |
|---|---|---|
| `capstone.Rmd` / `capstone.pdf` | 1.5 × IQR | `weather_selection_ols_h15_iqr1.5_v1/` |
| `capstone_iqr3.Rmd` / `capstone_iqr3.pdf` | 3 × IQR | `weather_selection_ols_h15_iqr3_v1/` |

Both runs select the same models: airport-specific first and second
harmonics, with `h` = 14 (maximum), 13 (mean) and 13 (minimum).

| Outcome | CV-MSE, 1.5 × IQR | CV-MSE, 3 × IQR |
|---|---|---|
| Daily maximum | 56.05 | 56.51 |
| Daily mean | 36.94 | 36.96 |
| Daily minimum | 35.46 | 35.76 |

## R files that produce the reported results

Three R files produce every number in the Results section of the
reports.

| File | Role | Output used by the report |
|---|---|---|
| `model_selection.R` | The whole procedure: data preparation, outlier rule, training tables for the five horizons, least-squares fit, sequential forecasts, selection of `h` and the seasonal structure, and the refit of the selected models on all ten years. | `best_settings.csv`, `horizon_performance.csv`, `final_models.rds` |
| `greatlakes/cv_task.R` | Runs one cross-validation task on the cluster: one outcome, one held-out year, one seasonal structure and one `h`. It loads its functions from sections 1–10 of `model_selection.R`, so the method exists in one place only. | One checkpoint file per task (1,920 per run) |
| `final_predictions.R` | Uses the final models to forecast every scored origin at the 15 target airports. | `final_predictions.rds` (the two forecast figures) |

The order is:

1. `greatlakes/cv_task.R` runs 1,920 times in parallel on the cluster
   (3 outcomes × 10 years × 4 seasonal structures × 16 values of `h`).
2. `model_selection.R` runs once. It finds the checkpoints, combines
   them, selects the model for each outcome and refits it on all ten
   years.
3. `final_predictions.R` runs once to produce the forecasts of the
   final models.

The settings at the top of `model_selection.R` decide which run the
three scripts work on. They are currently set to the 3 × IQR run:

```r
iqr_multiple <- 3
output_directory <- "weather_selection_ols_h15_iqr3_v1"
```

For the 1.5 × IQR run, set `iqr_multiple <- 1.5` and
`output_directory <- "weather_selection_ols_h15_iqr1.5_v1"`.

## What the Rmd files do themselves

- The exploratory analysis and the two harmonic figures are computed
  from `weather.rda` when the report is knitted.
- The code in the "Model Selection" section is displayed but not run.
  It is a shorter rewrite of `model_selection.R` using `lm()`, meant
  to be read. For one fold and three candidate models it gives the
  same squared errors as the cluster run, under both outlier rules.
- The "Results" section only reads the saved files listed above and
  formats them. Knitting does not refit anything.

## Other files

| File or folder | Contents |
|---|---|
| `weather.rda` | The data. |
| `greatlakes/` | Scripts for the Great Lakes cluster: `launch.sh`, `submit.sh`, `watch.sh` and `fetch.sh` copy the project, submit the jobs, show progress and copy results back; `cv_array.sbatch` and `aggregate.sbatch` are the two Slurm jobs; `setup_packages.R` installs the R packages. |
| `weather_selection_ols_h15_iqr1.5_v1/`, `weather_selection_ols_h15_iqr3_v1/` | Results of the two runs: the checkpoints, `cv_summary.csv` (CV-MSE of all 64 candidates per outcome), `cv_curves.pdf`, and the files listed in the table above. |
| `capstone.R` | Early exploratory script. Its code is now part of the Rmd files. |
| `Test/` | Experiments that are not part of the reports (see below). |

## Test folder

- `Test/interaction_cv.R` and `Test/summarize.R` test whether
  airport-by-predictor interactions improve the forecasts, using the
  same folds and scored days as the 3 × IQR run. Interacting the
  weather predictors with airport lowers the CV-MSE by about 1%
  overall (about 7% for the daily mean at SFO), mostly at one and two
  days ahead. Results are in `Test/results/`.
- The `Test/weather_selection_*` folders hold earlier runs: group
  lasso with different maximum history lengths, and month effects in
  place of harmonics. Neither improved on the reported models.

## Reproducing

All scripts are run from this folder.

```bash
# Cross-validation and selection on the cluster
bash greatlakes/launch.sh <uniqname>
bash greatlakes/watch.sh <uniqname>
bash greatlakes/fetch.sh <uniqname>

# Forecasts of the final models, then the report
Rscript final_predictions.R
Rscript -e 'rmarkdown::render("capstone_iqr3.Rmd")'
```

`model_selection.R` can also be run on its own (`Rscript
model_selection.R`); it skips every checkpoint that already exists.
