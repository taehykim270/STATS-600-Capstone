Airport temperature forecasting in R

Run from the extracted project folder containing data/weather.rda. Rscript can also take the full path to a script. Base R and the R C toolchain are required. The included C solver builds with R CMD SHLIB on first use. On macOS or Linux, AIRPORT_WORKERS controls parallel airport fits (default 4). Windows uses one worker.

Rscript common_design/airport_lasso.R
  Fit the featured ordinary Lasso using minimum, mean and maximum temperature at the other 97 airports. Save 225 final fits, signed nonzero coefficients, rolling-origin forecasts and scored tables under main_results/baseline/.
Rscript -e 'source("common_design/airport_lasso.R"); run_experiment("mean")'
  Repeat the matched comparison with mean temperature alone at the other airports. Results go to main_results/mean_only/.
Rscript common_design/preprocess_weather.R
  Export the featured model's training and evaluation design matrices to baseline/artifacts/weather_preprocessed.rds. The export uses a lossless column dictionary. Extract any matrix and its dated responses with get_model_data().
Rscript common_design/run_comparisons.R
  Refit every additional model in the note. Expanding-year CV chooses a shared penalty separately for Lasso and ridge. Write the development losses and final comparison fits under main_results/experiments/.
Rscript common_design/summarize_results.R
  Recompute the model comparisons, day-15 subset, feature selection summaries and lead/target results from the fitted runs.
Rscript common_design/plot_selected_airports.R
  Make the three featured mean-temperature selection maps for Detroit, Atlanta and San Francisco.
Rscript common_design/score_mse.R main_results/baseline/artifacts/forecast_errors_input.rds
  Recompute the featured MSE tables without refitting. Use the mean_only path to score that design.

Forecast from five consecutive daily rows (the row names must be ISO dates):
  source("common_design/airport_lasso.R")
  load("data/weather.rda")
  model <- readRDS("main_results/baseline/artifacts/airport_lasso.rds")
  dates <- seq(as.Date("2013-07-11"), as.Date("2013-07-15"), by="day")
  rows <- match(dates, as.Date(rownames(weather)))
  forecast_airport(model, weather[rows,,drop=FALSE], "DTW")

Inspect a design matrix after running the preprocessing export:
  source("common_design/preprocess_weather.R")
  d <- readRDS("main_results/baseline/artifacts/weather_preprocessed.rds")
  x <- get_model_data(d, "DTW", "Mean", 1, "train")
  head(x$observations)
  colnames(x$x)

For the other fitted models, source run_comparisons.R, then call forecast_comparison(five_rows, "DTW", "pc8"). See the model IDs in experiments/model_key.csv. The day-15 row uses the featured model. Positive integer lags look backwards from the last observed day. Five direct predictions always use the exact same five observed days.
