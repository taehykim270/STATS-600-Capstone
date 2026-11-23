# Forecasts of the final (all-years) models at the target airports.
#
# Run from the folder containing weather.rda, after model_selection.R
# has written final_models.rds:
#   Rscript final_predictions.R
#
# For every complete forecast origin it produces the five sequential
# forecasts, exactly as in the cross-validation, and saves them next
# to the other results as final_predictions.rds. capstone.Rmd reads
# that file for its figures, so nothing is refitted when knitting.
#
# These are in-sample forecasts: the models were fitted to the same
# ten years.

# Settings, data and helper functions: sections 1-10 of model_selection.R.
script <- readLines("model_selection.R", warn = FALSE)
cut <- grep("^# 11\\. Run outcomes in parallel", script)
stopifnot(length(cut) == 1)

suppressPackageStartupMessages(
  eval(parse(text = script[seq_len(cut - 1)]), envir = globalenv())
)

final_models <- readRDS(file.path(output_directory, "final_models.rds"))

# Outlier fences from all ten years, as in the final refit.
cleaned <- clean_temperature_outliers(weather_data, NA)

predict_outcome <- function(outcome) {

  final <- final_models[[outcome]]
  weather_columns <- setdiff(numeric_weather, outcome)

  tables <- make_horizon_tables(cleaned, outcome)

  # The five tables share their rows; keep the target airports.
  origins <- lapply(
    tables,
    function(x) x[x$airport %in% target_airports, ]
  )

  # Complete origins: weather, 15 observed temperatures, five targets.
  keep <- complete_rows(
    origins[[1]],
    c(
      "airport", "Events", weather_columns,
      paste0("observed_", seq_len(max_history))
    )
  )

  for (horizon in 1:5) {
    keep <- keep & is.finite(origins[[horizon]]$response)
  }

  origins <- lapply(origins, function(x) x[keep, ])

  prepared <- lapply(origins, function(valid) {

    # Start the histories from the observed values.
    for (j in seq_len(max_history)) {
      valid[[paste0("history_", j)]] <- valid[[paste0("observed_", j)]]
    }

    list(
      x_validation = make_matrix(
        valid,
        c(history_names(max_history), weather_columns),
        final$airport_levels,
        final$event_levels,
        final$seasonal %in% c("airport1", "airport2")
      ),
      y_validation = valid$response,
      cell_validation = paste(valid$airport, valid$month),
      observed_history = as.matrix(
        valid[, paste0("observed_", seq_len(max_history))]
      )
    )
  })

  predictions <- predict_sequence(final$horizons, prepared, final$h)

  dplyr::bind_rows(lapply(1:5, function(horizon) {
    data.frame(
      outcome = outcome,
      airport = origins[[horizon]]$airport,
      origin_date = origins[[horizon]]$date,
      horizon = horizon,
      actual = origins[[horizon]]$response,
      predicted = predictions[, horizon]
    )
  }))
}

final_predictions <- dplyr::bind_rows(lapply(outcomes, predict_outcome))

saveRDS(
  final_predictions,
  file.path(output_directory, "final_predictions.rds")
)

print(
  final_predictions %>%
    dplyr::group_by(outcome, horizon) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mse = mean((actual - predicted)^2),
      .groups = "drop"
    ) %>%
    as.data.frame()
)
