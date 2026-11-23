# Summarize the interaction test. Run from the project folder after
# Test/interaction_cv.R:
#   Rscript Test/summarize.R
#
# Writes to Test/results:
#   summary_overall.csv     CV-MSE by outcome and model (all airports,
#                           all horizons), and the change from pooled
#   summary_by_airport.csv  CV-MSE by outcome, airport and model
#   summary_by_horizon.csv  CV-MSE by outcome, horizon and model

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

results_directory <- file.path("Test", "results")

errors <- read.csv(file.path(results_directory, "interaction_cv_errors.csv"))

model_order <- c("pooled", "lags", "weather", "lags_weather", "separate")

cv_mse <- function(data, ...) {
  data %>%
    group_by(..., model) %>%
    summarise(mse = sum(sse) / sum(n), .groups = "drop") %>%
    pivot_wider(names_from = model, values_from = mse) %>%
    select(..., all_of(model_order))
}

# Percentage change in CV-MSE relative to the pooled model.
change <- function(table) {
  table %>%
    mutate(across(
      all_of(model_order[-1]),
      ~ 100 * (.x / pooled - 1),
      .names = "{.col}_pct"
    ))
}

overall <- change(cv_mse(errors, outcome))
by_airport <- change(cv_mse(errors, outcome, airport))
by_horizon <- change(cv_mse(errors, outcome, horizon))

# Number of coefficients per horizon model (largest over folds).
coefficients <- errors %>%
  group_by(outcome, model) %>%
  summarise(n_coefficients = max(n_coefficients), .groups = "drop") %>%
  pivot_wider(names_from = model, values_from = n_coefficients) %>%
  select(outcome, all_of(model_order))

# In how many of the ten folds does each model beat the pooled model?
folds_better <- errors %>%
  group_by(outcome, year, model) %>%
  summarise(mse = sum(sse) / sum(n), .groups = "drop") %>%
  group_by(outcome, year) %>%
  mutate(better = mse < mse[model == "pooled"]) %>%
  group_by(outcome, model) %>%
  summarise(folds_better = sum(better), .groups = "drop") %>%
  filter(model != "pooled") %>%
  pivot_wider(names_from = model, values_from = folds_better) %>%
  select(outcome, all_of(model_order[-1]))

write.csv(overall, file.path(results_directory, "summary_overall.csv"),
          row.names = FALSE)
write.csv(by_airport, file.path(results_directory, "summary_by_airport.csv"),
          row.names = FALSE)
write.csv(by_horizon, file.path(results_directory, "summary_by_horizon.csv"),
          row.names = FALSE)

options(width = 200)

cat("\nCV-MSE by outcome (pct = change from pooled, %)\n")
print(as.data.frame(overall), digits = 4)

cat("\nCoefficients per horizon model\n")
print(as.data.frame(coefficients))

cat("\nFolds (of 10) in which the model beats pooled\n")
print(as.data.frame(folds_better))

cat("\nCV-MSE by airport, daily mean temperature\n")
print(
  as.data.frame(by_airport %>% filter(outcome == "Mean.TemperatureF")),
  digits = 3
)

cat("\nSFO, all outcomes\n")
print(as.data.frame(by_airport %>% filter(airport == "SFO")), digits = 3)

cat("\nCV-MSE by horizon, daily mean temperature\n")
print(
  as.data.frame(by_horizon %>% filter(outcome == "Mean.TemperatureF")),
  digits = 4
)
