# Test: do airport-by-predictor interactions improve the forecasts?
#
# Run from the project folder (the one containing weather.rda):
#   Rscript Test/interaction_cv.R            # all folds, three outcomes
#   Rscript Test/interaction_cv.R check      # one fold, compare with cluster
#
# Everything except the model is taken from model_selection.R: the
# data, the outlier rule, the complete-case rules, the ten folds by
# year, the scored origins and the sequential forecasts. For each
# outcome the history length h is the one selected by the main run.
#
# Models compared (all fitted by least squares):
#   pooled       the selected model: airport intercepts and
#                airport-specific harmonics; all other coefficients
#                are common to all airports
#   lags         pooled + airport x temperature lags
#   weather      pooled + airport x weather predictors
#   lags_weather pooled + airport x (temperature lags and weather)
#   separate     one regression per airport: every coefficient,
#                including the event effects, is airport-specific
#
# Output: Test/results/fold_<outcome>_<year>.rds (one per fold) and
# Test/results/interaction_cv_errors.csv, with the squared error by
# model, outcome, year, airport and horizon.

suppressPackageStartupMessages(library(Matrix))

mode <- commandArgs(trailingOnly = TRUE)
mode <- if (length(mode) == 0) "full" else mode[1]

# Settings, data and helper functions: sections 1-10 of model_selection.R.
script <- readLines("model_selection.R", warn = FALSE)
cut <- grep("^# 11\\. Run outcomes in parallel", script)
stopifnot(length(cut) == 1)

suppressWarnings(suppressPackageStartupMessages(
  eval(parse(text = script[seq_len(cut - 1)]), envir = globalenv())
))

results_directory <- file.path("Test", "results")
dir.create(results_directory, showWarnings = FALSE, recursive = TRUE)

best_settings <- read.csv(file.path(output_directory, "best_settings.csv"))
selected_h <- setNames(best_settings$h, best_settings$outcome)
stopifnot(all(best_settings$seasonal == "airport2"))

model_names <- c("pooled", "lags", "weather", "lags_weather", "separate")

harmonics <- paste0("season_", c("sin1", "cos1", "sin2", "cos2"))

# --------------------------------------------------
# Design matrix
# --------------------------------------------------

# x:       dense matrix from prepare_fold() (raw values)
# airport: airport of each row
# info:    list(h, weather, events, airport_levels, center, scale)
# Returns a sparse matrix. Numeric predictors are centred and scaled
# with the training means and standard deviations; this does not
# change the fitted values and keeps the normal equations stable.
build_design <- function(x, airport, info, model) {

  numeric_columns <- c(history_names(info$h), info$weather)

  z <- sweep(x[, numeric_columns, drop = FALSE], 2, info$center, "-")
  z <- sweep(z, 2, info$scale, "/")

  # Indicator of every airport except the reference one.
  levels_used <- info$airport_levels[-1]
  index <- match(airport, levels_used)
  rows <- which(!is.na(index))

  indicator <- sparseMatrix(
    i = rows,
    j = index[rows],
    x = 1,
    dims = c(nrow(x), length(levels_used)),
    dimnames = list(NULL, paste0("airport_", levels_used))
  )

  # Airport indicator times one column.
  interact <- function(values, name) {
    block <- Diagonal(x = values) %*% indicator
    colnames(block) <- paste0(name, ":", levels_used)
    block
  }

  blocks <- list(
    Matrix(cbind("(Intercept)" = 1, z), sparse = TRUE),
    Matrix(x[, c(harmonics, info$events), drop = FALSE], sparse = TRUE),
    indicator
  )

  # Airport-specific harmonics (the selected seasonal structure).
  for (name in harmonics) {
    blocks <- c(blocks, interact(x[, name], name))
  }

  interacted <- c(
    if (model %in% c("lags", "lags_weather")) history_names(info$h),
    if (model %in% c("weather", "lags_weather")) info$weather
  )

  for (name in interacted) {
    blocks <- c(blocks, interact(z[, name], name))
  }

  do.call(cbind, blocks)
}

# --------------------------------------------------
# Least squares through the normal equations
# --------------------------------------------------

# Columns that are exact linear combinations of earlier ones get a
# zero coefficient, as in lm().
solve_least_squares <- function(design, y) {

  xtx <- as.matrix(crossprod(design))
  xty <- as.numeric(crossprod(design, y))

  size <- sqrt(diag(xtx))
  used <- which(size > 0)

  # Scale to unit diagonal so that the rank tolerance is meaningful.
  scaled <- xtx[used, used] / tcrossprod(size[used])

  factor <- suppressWarnings(chol(scaled, pivot = TRUE, tol = 1e-9))
  rank <- attr(factor, "rank")
  pivot <- attr(factor, "pivot")[seq_len(rank)]

  solution <- backsolve(
    factor[seq_len(rank), seq_len(rank), drop = FALSE],
    forwardsolve(
      t(factor[seq_len(rank), seq_len(rank), drop = FALSE]),
      (xty[used] / size[used])[pivot]
    )
  )

  beta <- numeric(ncol(design))
  beta[used[pivot]] <- solution / size[used[pivot]]

  list(beta = beta, n_aliased = ncol(design) - rank)
}

# --------------------------------------------------
# Fit and forecast
# --------------------------------------------------

# History inputs of the horizon model at validation: observed values
# up to the origin day, forecasts after it (as in predict_sequence()).
fill_history <- function(x, prepared_horizon, predictions, horizon, h) {

  for (j in seq_len(h)) {

    relative_day <- horizon - j

    x[, paste0("history_", j)] <- if (relative_day > 0) {
      predictions[, relative_day]
    } else {
      prepared_horizon$observed_history[, 1 - relative_day]
    }
  }

  x
}

# Pooled model, with or without airport interactions.
forecast_pooled <- function(prepared, h, model, weather, events,
                            airport_levels) {

  predictions <- matrix(NA_real_, length(prepared[[1]]$y_validation), 5)
  n_columns <- n_aliased <- integer(5)

  for (horizon in 1:5) {

    current <- prepared[[horizon]]
    numeric_columns <- c(history_names(h), weather)

    keep <- rowSums(
      !is.finite(current$x_training[, numeric_columns, drop = FALSE])
    ) == 0

    x <- current$x_training[keep, , drop = FALSE]

    info <- list(
      h = h,
      weather = weather,
      events = events,
      airport_levels = airport_levels,
      center = colMeans(x[, numeric_columns, drop = FALSE]),
      scale = apply(x[, numeric_columns, drop = FALSE], 2, sd)
    )

    design <- build_design(
      x, sub(" .*$", "", current$cell_training[keep]), info, model
    )

    fit <- solve_least_squares(design, current$y_training[keep])
    rm(design, x)

    x_new <- fill_history(
      current$x_validation, current, predictions, horizon, h
    )

    design_new <- build_design(
      x_new, sub(" .*$", "", current$cell_validation), info, model
    )

    predictions[, horizon] <- as.numeric(design_new %*% fit$beta)
    n_columns[horizon] <- length(fit$beta)
    n_aliased[horizon] <- fit$n_aliased
  }

  list(predictions = predictions, n_coefficients = n_columns - n_aliased)
}

# One regression per target airport.
forecast_separate <- function(prepared, h, weather, events) {

  predictions <- matrix(NA_real_, length(prepared[[1]]$y_validation), 5)
  columns <- c(history_names(h), weather, harmonics, events)
  n_coefficients <- integer(5)

  for (horizon in 1:5) {

    current <- prepared[[horizon]]
    airport_training <- sub(" .*$", "", current$cell_training)
    airport_validation <- sub(" .*$", "", current$cell_validation)

    x_new <- fill_history(
      current$x_validation, current, predictions, horizon, h
    )

    for (airport in target_airports) {

      rows <- which(airport_training == airport)
      x <- current$x_training[rows, columns, drop = FALSE]
      keep <- rowSums(!is.finite(x)) == 0

      fit <- lm.fit(
        cbind(1, x[keep, , drop = FALSE]),
        current$y_training[rows][keep]
      )

      beta <- fit$coefficients
      n_coefficients[horizon] <- n_coefficients[horizon] + sum(!is.na(beta))
      beta[is.na(beta)] <- 0

      new_rows <- which(airport_validation == airport)

      predictions[new_rows, horizon] <- as.numeric(
        cbind(1, x_new[new_rows, columns, drop = FALSE]) %*% beta
      )
    }
  }

  list(predictions = predictions, n_coefficients = n_coefficients)
}

# --------------------------------------------------
# One fold
# --------------------------------------------------

run_fold <- function(outcome, held_out_year) {

  file <- file.path(
    results_directory,
    paste0("fold_", outcome, "_", held_out_year, ".rds")
  )

  if (file.exists(file)) {
    return(readRDS(file))
  }

  h <- selected_h[[outcome]]
  weather <- setdiff(numeric_weather, c(outcome, harmonics))

  cleaned <- clean_temperature_outliers(weather_data, held_out_year)
  tables <- make_horizon_tables(cleaned, outcome)
  prepared <- prepare_fold(tables, outcome, held_out_year)$horizons
  rm(cleaned, tables)

  all_columns <- colnames(prepared[[1]]$x_training)
  events <- grep("^event_", all_columns, value = TRUE)

  airport_levels <- c(
    setdiff(
      unique(sub(" .*$", "", prepared[[1]]$cell_training)),
      sub("^airport_", "", grep("^airport_", all_columns, value = TRUE))
    ),
    sub("^airport_", "", grep("^airport_", all_columns, value = TRUE))
  )

  airport <- sub(" .*$", "", prepared[[1]]$cell_validation)

  results <- lapply(model_names, function(model) {

    started <- Sys.time()

    forecast <- if (model == "separate") {
      forecast_separate(prepared, h, weather, events)
    } else {
      forecast_pooled(prepared, h, model, weather, events, airport_levels)
    }

    stopifnot(all(is.finite(forecast$predictions)))

    cat(sprintf(
      "%s | %d | %-12s | %5.1f min\n",
      outcome, held_out_year, model,
      as.numeric(difftime(Sys.time(), started, units = "mins"))
    ))

    do.call(rbind, lapply(1:5, function(horizon) {

      squared_error <- (
        prepared[[horizon]]$y_validation - forecast$predictions[, horizon]
      )^2

      data.frame(
        model = model,
        outcome = outcome,
        year = held_out_year,
        h = h,
        horizon = horizon,
        airport = names(tapply(squared_error, airport, sum)),
        sse = as.numeric(tapply(squared_error, airport, sum)),
        n = as.numeric(tapply(squared_error, airport, length)),
        n_coefficients = forecast$n_coefficients[horizon]
      )
    }))
  })

  results <- do.call(rbind, results)
  saveRDS(results, file)
  results
}

# --------------------------------------------------
# Run
# --------------------------------------------------

if (mode == "check") {

  # The pooled model must reproduce the cluster's checkpoint.
  outcome <- "Mean.TemperatureF"
  result <- run_fold(outcome, 2010)

  cluster <- readRDS(file.path(
    output_directory,
    sprintf("%s_year_2010_season_airport2_h_%d.rds",
            outcome, selected_h[[outcome]])
  ))

  pooled <- result[result$model == "pooled", ]

  print(data.frame(
    horizon = 1:5,
    test = as.numeric(tapply(pooled$sse, pooled$horizon, sum)),
    cluster = cluster$sse,
    n_test = as.numeric(tapply(pooled$n, pooled$horizon, sum)),
    n_cluster = cluster$n
  ))

  print(aggregate(sse ~ model, result, sum))

} else {

  tasks <- expand.grid(year = years, outcome = outcomes,
                       stringsAsFactors = FALSE)

  results <- parallel::mclapply(
    seq_len(nrow(tasks)),
    function(i) run_fold(tasks$outcome[i], tasks$year[i]),
    mc.cores = 3,
    mc.preschedule = FALSE
  )

  failed <- vapply(results, inherits, logical(1), what = "try-error")

  if (any(failed)) {
    print(results[failed])
    stop(sum(failed), " folds failed; rerun to retry them.")
  }

  write.csv(
    do.call(rbind, results),
    file.path(results_directory, "interaction_cv_errors.csv"),
    row.names = FALSE
  )

  cat("\nFinished. Results are in:", results_directory, "\n")
}
