# Forecasting airport temperatures

Forecast daily minimum, mean and maximum temperature at 15 U.S. airports from the five most recent days of U.S. airport temperature and pressure. The featured model fits a separate ordinary Lasso for each forecast airport, temperature type and lead. Its 1,669 candidate predictors include all three temperatures at each of 97 outside airports. All predictors receive the same Lasso penalty.

Read the [three-page model note](airport_temperature_model_note.pdf) or open the [complete model score table](main_results/experiments/model_scores.csv). Results use training through 2012, expanding-year cross-validation, and separate 2013 and 2014 evaluations.

| Directory | Contents |
| --- | --- |
| data | The 98 stations' three daily temperatures and mean sea-level pressure, with dates in the row names |
| common_design | R preprocessing, model fitting, scoring and plotting scripts, the C Lasso solver, and geographic reference data |
| main_results/baseline | Featured 1,669-input model, signed coefficients, forecasts, MSE tables and selection maps |
| main_results/mean_only | Matched 699-input model using only external daily mean temperatures |
| main_results/experiments | Fits and tables for every additional model in the note, stored under the IDs in model_key.csv |

## Reproduce the results

Run these commands from the project root with R and its C build toolchain installed.

```sh
Rscript common_design/airport_lasso.R
Rscript -e 'source("common_design/airport_lasso.R"); run_experiment("mean")'
Rscript common_design/run_comparisons.R
Rscript common_design/summarize_results.R
Rscript common_design/plot_selected_airports.R
```

Run `Rscript common_design/preprocess_weather.R` to rebuild the full standardized matrix for all 225 featured fits. See the [R usage guide](common_design/README.txt) for loading one matrix or issuing forecasts from five daily rows. Geographic inputs use OpenFlights coordinates and 2023 U.S. Census state outlines.
