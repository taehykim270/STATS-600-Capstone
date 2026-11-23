# Weather model selection
# Run this script from the folder containing weather.rda.

# Install once if needed:
# install.packages(c(
#   "dplyr", "tidyr", "tibble", "ggplot2",
#   "future", "future.apply"
# ))

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(future)
library(future.apply)

# --------------------------------------------------
# 1. Settings
# --------------------------------------------------

data_file <- "weather.rda"

outcomes <- c(
  "Max.TemperatureF",
  "Mean.TemperatureF",
  "Min.TemperatureF"
)

target_airports <- c(
  "ABQ", "ATL", "BOS", "DEN", "DFW",
  "DLH", "DTW", "EYW", "IAD", "LAN",
  "LAX", "SEA", "SFO", "SLC", "TVC"
)

years <- 2005:2014
history_grid <- 0:15
max_history <- max(history_grid)

# Seasonal structures compared in this run, simplest to richest.
#   none:     no seasonal terms
#   shared:   two annual harmonics common to all airports
#   airport1: shared, plus an airport-specific first harmonic
#   airport2: shared, plus airport-specific first and second harmonics
#
# Also supported (used in weather_selection_month_h15_v1, with
# seasonal_grid <- c("month", "month_airport")):
#   month:         month effects common to all airports
#   month_airport: month effects plus airport-specific month effects
seasonal_grid <- c("none", "shared", "airport1", "airport2")

# Each candidate (seasonal structure, h) is fitted by ordinary least
# squares. Earlier group-lasso runs always chose the weakest penalty,
# so no shrinkage is applied.

# Temperatures further than this multiple of the IQR from the
# quartiles are treated as outliers. weather_selection_ols_h15_iqr1.5_v1
# used 1.5.
iqr_multiple <- 3

# "rolling" or "monthly_15".
validation_origins <- "rolling"

# Use 1 if memory is limited.
number_of_workers <- min(
  3L,
  as.integer(future::availableCores())
)

# Change this folder if data or modeling rules change.
output_directory <- "weather_selection_ols_h15_iqr3_v1"

dir.create(
  output_directory,
  showWarnings = FALSE,
  recursive = TRUE
)

# --------------------------------------------------
# 2. Load and reshape data
# --------------------------------------------------

load(data_file)

weather_long <- weather %>%
  tibble::rownames_to_column(var = "date") %>%
  dplyr::mutate(
    dplyr::across(
      matches("^[A-Z]{3}_PrecipitationIn$"),
      ~ {
        x <- trimws(as.character(.x))
        
        # Preserve the precipitation convention from the EDA.
        x[x %in% c("T", "")] <- NA_character_
        
        values <- suppressWarnings(as.numeric(x))
        
        unknown <- unique(x[
          !is.na(x) & is.na(values)
        ])
        
        if (length(unknown) > 0) {
          stop(
            "Unknown precipitation codes: ",
            paste(unknown, collapse = ", ")
          )
        }
        
        values
      }
    )
  ) %>%
  tidyr::pivot_longer(
    cols = matches("^[A-Z]{3}_"),
    names_to = c("airport", ".value"),
    names_pattern = "^([^_]+)_(.*)$"
  ) %>%
  dplyr::mutate(date = as.Date(date))

weather_data <- weather_long %>%
  dplyr::arrange(airport, date) %>%
  dplyr::mutate(
    year = as.integer(format(date, "%Y")),
    day_of_year = as.integer(format(date, "%j")),
    day_of_month = as.integer(format(date, "%d")),
    month = as.integer(format(date, "%m")),
    wind_sin = sin(WindDirDegrees * pi / 180),
    wind_cos = cos(WindDirDegrees * pi / 180),
    season_sin1 = sin(2 * pi * day_of_year / 365.25),
    season_cos1 = cos(2 * pi * day_of_year / 365.25),
    season_sin2 = sin(4 * pi * day_of_year / 365.25),
    season_cos2 = cos(4 * pi * day_of_year / 365.25),
    Events = trimws(as.character(Events))
  )

# Blank events mean no event. NA remains missing.
weather_data$Events[
  !is.na(weather_data$Events) &
    weather_data$Events == ""
] <- "None"

stopifnot(!anyNA(weather_data$date))
stopifnot(!anyNA(weather_data$airport))
stopifnot(all(outcomes %in% names(weather_data)))
stopifnot(all(target_airports %in% weather_data$airport))

duplicates <- weather_data %>%
  dplyr::count(airport, date) %>%
  dplyr::filter(n != 1)

if (nrow(duplicates) > 0) {
  stop("Duplicate airport-date observations.")
}

date_checks <- weather_data %>%
  dplyr::group_by(airport) %>%
  dplyr::summarise(
    has_gap = any(diff(as.integer(date)) != 1),
    .groups = "drop"
  )

if (any(date_checks$has_gap)) {
  stop("Missing calendar dates. Complete the date grid first.")
}

numeric_weather <- names(weather_data)[
  vapply(weather_data, is.numeric, logical(1))
]

numeric_weather <- setdiff(
  numeric_weather,
  c(
    "year",
    "day_of_year",
    "day_of_month",
    "month",
    "WindDirDegrees"
  )
)

# Avoid keeping unnecessary copies during parallel execution.
rm(weather, weather_long, duplicates, date_checks)
invisible(gc())

# --------------------------------------------------
# 3. Checkpoint settings
# --------------------------------------------------

run_settings <- list(
  version = "ols_season_h15_complete_case_v1",
  outcomes = outcomes,
  years = years,
  history_grid = history_grid,
  seasonal_grid = seasonal_grid,
  validation_origins = validation_origins,
  target_airports = target_airports,
  numeric_weather = numeric_weather
)

# Recorded only when it differs from the original 1.5, so that the
# settings saved by earlier runs still match.
if (iqr_multiple != 1.5) {
  run_settings$iqr_multiple <- iqr_multiple
}

settings_file <- file.path(
  output_directory,
  "run_settings.rds"
)

if (file.exists(settings_file)) {
  
  saved_settings <- readRDS(settings_file)
  
  # all.equal tolerates last-digit floating point differences
  # between machines (e.g. laptop vs. cluster).
  if (!isTRUE(all.equal(saved_settings, run_settings))) {
    stop(
      "Saved settings differ. Use a new output directory."
    )
  }
  
} else {
  saveRDS(run_settings, settings_file)
}

# --------------------------------------------------
# 4. Cleaning helpers
# --------------------------------------------------

complete_rows <- function(data, columns) {
  
  keep <- rep(TRUE, nrow(data))
  
  for (column_name in columns) {
    
    values <- data[[column_name]]
    
    if (is.numeric(values)) {
      keep <- keep & is.finite(values)
    } else {
      keep <- keep & !is.na(values)
    }
  }
  
  keep
}

# held_out_year = NA uses every year to set the cutoffs (final refit).
clean_temperature_outliers <- function(data, held_out_year) {
  
  cleaned <- data
  
  for (airport_name in unique(data$airport)) {
    
    airport_rows <- data$airport == airport_name
    
    training_rows <- airport_rows &
      (is.na(held_out_year) | data$year != held_out_year)
    
    for (outcome in outcomes) {
      
      training_values <- data[[outcome]][training_rows]
      training_values <- training_values[
        is.finite(training_values)
      ]
      
      if (length(training_values) == 0) {
        cleaned[[outcome]][airport_rows] <- NA_real_
        next
      }
      
      quartiles <- quantile(
        training_values,
        probs = c(0.25, 0.75),
        names = FALSE
      )
      
      iqr <- quartiles[2] - quartiles[1]
      lower <- quartiles[1] - iqr_multiple * iqr
      upper <- quartiles[2] + iqr_multiple * iqr
      
      values <- cleaned[[outcome]]
      
      bad <- airport_rows &
        is.finite(values) &
        (values < lower | values > upper)
      
      values[bad] <- NA_real_
      values[!is.finite(values)] <- NA_real_
      
      cleaned[[outcome]] <- values
    }
  }
  
  cleaned
}

history_names <- function(h) {
  
  if (h == 0) {
    return(character(0))
  }
  
  paste0("history_", seq_len(h))
}

shift_temperature <- function(x, offset) {
  
  if (offset > 0) {
    return(dplyr::lead(x, offset))
  }
  
  if (offset < 0) {
    return(dplyr::lag(x, -offset))
  }
  
  x
}

# --------------------------------------------------
# 5. Construct all five horizon tables
# --------------------------------------------------

make_horizon_tables <- function(data, outcome) {
  
  weather_columns <- setdiff(numeric_weather, outcome)
  airport_data <- split(data, data$airport)
  
  horizon_tables <- vector("list", 5)
  
  for (horizon in 1:5) {
    
    tables <- vector("list", length(airport_data))
    
    for (i in seq_along(airport_data)) {
      
      current <- airport_data[[i]]
      y <- current[[outcome]]
      
      result <- current[, c(
        "airport", "date", "year",
        "day_of_year", "day_of_month", "month", "Events",
        weather_columns
      )]
      
      result$response <- dplyr::lead(y, horizon)
      
      for (j in seq_len(max_history)) {
        
        result[[paste0("history_", j)]] <- shift_temperature(
          y,
          horizon - j
        )
        
        result[[paste0("observed_", j)]] <- dplyr::lag(
          y,
          j - 1
        )
      }
      
      # Keep histories and targets within the same year.
      keep <- result$year %in% years &
        result$day_of_year >= max_history &
        as.integer(format(result$date + 5, "%Y")) ==
        result$year
      
      tables[[i]] <- result[keep, ]
    }
    
    horizon_tables[[horizon]] <- dplyr::bind_rows(tables)
  }
  
  horizon_tables
}

# --------------------------------------------------
# 6. Predictor matrices
# --------------------------------------------------

make_matrix <- function(data, numeric_columns, airport_levels,
                        event_levels, airport_season = FALSE) {
  
  x <- as.matrix(
    data[, numeric_columns, drop = FALSE]
  )
  
  storage.mode(x) <- "double"
  
  # Missing categories remain missing until rows are excluded.
  airport_values <- data$airport
  event_values <- data$Events
  
  for (level in airport_levels[-1]) {
    
    x <- cbind(
      x,
      as.numeric(airport_values == level)
    )
    
    colnames(x)[ncol(x)] <- paste0("airport_", level)
  }
  
  # Airport-specific seasonal harmonics: airport indicator times
  # each shared harmonic. Named airportseason<k>_<sin|cos>_<airport>.
  if (airport_season) {
    
    airport_columns <- paste0("airport_", airport_levels[-1])
    
    blocks <- lapply(
      c("sin1", "cos1", "sin2", "cos2"),
      function(harmonic) {
        
        block <- x[, airport_columns, drop = FALSE] *
          x[, paste0("season_", harmonic)]
        
        colnames(block) <- paste0(
          "airportseason", substr(harmonic, 4, 4), "_",
          substr(harmonic, 1, 3), "_",
          airport_levels[-1]
        )
        
        block
      }
    )
    
    x <- cbind(x, do.call(cbind, blocks))
  }
  
  # Month of the forecast origin; January is the reference.
  month_block <- vapply(
    2:12,
    function(m) as.numeric(data$month == m),
    numeric(nrow(data))
  )
  
  dim(month_block) <- c(nrow(data), 11)
  colnames(month_block) <- paste0("month_", 2:12)
  
  x <- cbind(x, month_block)
  
  for (level in event_levels[-1]) {
    
    x <- cbind(
      x,
      as.numeric(event_values == level)
    )
    
    colnames(x)[ncol(x)] <- paste0("event_", level)
  }
  
  x
}

# Not used in fitting; kept to describe predictor blocks.
make_group_labels <- function(column_names) {
  
  labels <- column_names
  
  labels[grepl("^history_", column_names)] <- "history"
  labels[grepl("^airport_", column_names)] <- "airport"
  labels[grepl("^season_", column_names)] <- "season"
  labels[grepl("^month_", column_names)] <- "month"
  labels[
    grepl("^airportseason1_", column_names)
  ] <- "season_airport1"
  labels[
    grepl("^airportseason2_", column_names)
  ] <- "season_airport2"
  labels[grepl("^event_", column_names)] <- "events"
  
  labels[
    column_names %in% c("wind_sin", "wind_cos")
  ] <- "wind_direction"
  
  labels
}

# --------------------------------------------------
# 7. Prepare one fold
# --------------------------------------------------

prepare_fold <- function(tables, outcome, held_out_year,
                         airport_season = FALSE) {
  
  weather_columns <- setdiff(numeric_weather, outcome)
  
  training_tables <- lapply(
    tables,
    function(x) x[x$year != held_out_year, ]
  )
  
  validation_tables <- lapply(
    tables,
    function(x) {
      
      keep <- x$year == held_out_year &
        x$airport %in% target_airports
      
      if (validation_origins == "monthly_15") {
        keep <- keep & x$day_of_month == 15
      }
      
      x[keep, ]
    }
  )
  
  # All candidates use the same complete validation origins.
  required_origin <- c(
    "airport",
    "Events",
    weather_columns,
    paste0("observed_", seq_len(max_history))
  )
  
  keep_validation <- complete_rows(
    validation_tables[[1]],
    required_origin
  )
  
  for (horizon in 1:5) {
    
    stopifnot(identical(
      validation_tables[[1]]$date,
      validation_tables[[horizon]]$date
    ))
    
    stopifnot(identical(
      validation_tables[[1]]$airport,
      validation_tables[[horizon]]$airport
    ))
    
    keep_validation <- keep_validation &
      is.finite(validation_tables[[horizon]]$response)
  }
  
  n_before <- length(keep_validation)
  n_after <- sum(keep_validation)
  
  if (n_after == 0) {
    stop(
      "No complete validation sequences for ",
      outcome, " in ", held_out_year,
      ". Check missingness, especially precipitation and gusts."
    )
  }
  
  validation_tables <- lapply(
    validation_tables,
    function(x) x[keep_validation, ]
  )
  
  # Learn categorical levels from training data only.
  airport_levels <- sort(unique(
    training_tables[[1]]$airport
  ))
  
  event_levels <- sort(unique(
    training_tables[[1]]$Events[
      !is.na(training_tables[[1]]$Events)
    ]
  ))
  
  numeric_columns <- c(
    history_names(max_history),
    weather_columns
  )
  
  prepared <- vector("list", 5)
  
  for (horizon in 1:5) {
    
    train <- training_tables[[horizon]]
    valid <- validation_tables[[horizon]]
    
    # Remove rows missing fixed weather inputs or targets.
    keep_training <- complete_rows(
      train,
      c("response", "airport", "Events", weather_columns)
    )
    
    train <- train[keep_training, ]
    
    x_training <- make_matrix(
      train,
      numeric_columns,
      airport_levels,
      event_levels,
      airport_season
    )
    
    # Initialize histories with available observed values.
    for (j in seq_len(max_history)) {
      valid[[paste0("history_", j)]] <-
        valid[[paste0("observed_", j)]]
    }
    
    x_validation <- make_matrix(
      valid,
      numeric_columns,
      airport_levels,
      event_levels,
      airport_season
    )
    
    observed_history <- as.matrix(
      valid[, paste0("observed_", seq_len(max_history))]
    )
    
    prepared[[horizon]] <- list(
      x_training = x_training,
      y_training = train$response,
      cell_training = paste(train$airport, train$month),
      x_validation = x_validation,
      y_validation = valid$response,
      cell_validation = paste(valid$airport, valid$month),
      observed_history = observed_history
    )
  }
  
  list(
    horizons = prepared,
    origins_before = n_before,
    origins_after = n_after
  )
}

# --------------------------------------------------
# 8. Fit one horizon model by ordinary least squares
# --------------------------------------------------

# details = TRUE also returns standard errors, the residual standard
# error and R-squared (used for the final refit only).
fit_ols <- function(prepared, h, seasonal, details = FALSE) {
  
  all_columns <- colnames(prepared$x_training)
  
  is_shared <- grepl("^season_", all_columns)
  is_airport1 <- grepl("^airportseason1_", all_columns)
  is_airport2 <- grepl("^airportseason2_", all_columns)
  is_month <- grepl("^month_", all_columns)
  is_airport <- grepl("^airport_", all_columns)
  
  none <- rep(FALSE, length(all_columns))
  
  use_seasonal <- switch(
    seasonal,
    none = none,
    shared = is_shared,
    airport1 = is_shared | is_airport1,
    airport2 = is_shared | is_airport1 | is_airport2,
    month = is_month,
    month_airport = none,
    stop("Unknown seasonal structure: ", seasonal)
  )
  
  if (seasonal %in% c("airport1", "airport2") &&
      !any(is_airport1)) {
    stop("Airport-specific seasonal columns were not built.")
  }
  
  is_other <- !grepl("^history_", all_columns) &
    !is_shared & !is_airport1 & !is_airport2 & !is_month
  
  # Airport, month and airport-by-month effects together are one
  # intercept per airport-month cell. They are fitted exactly by
  # centering within cells instead of building ~1,100 indicators.
  by_cell <- seasonal == "month_airport"
  
  if (by_cell) {
    is_other <- is_other & !is_airport
  }
  
  columns <- c(
    history_names(h),
    all_columns[is_other | use_seasonal]
  )
  
  x <- prepared$x_training[, columns, drop = FALSE]
  y <- prepared$y_training
  
  # Exclude missing history inputs for this h.
  keep_rows <- rowSums(!is.finite(x)) == 0
  x <- x[keep_rows, , drop = FALSE]
  y <- y[keep_rows]
  
  if (length(y) < 2) {
    stop("Too few complete training observations.")
  }
  
  keep_columns <- apply(
    x,
    2,
    function(values) length(unique(values)) > 1
  )
  
  x <- x[, keep_columns, drop = FALSE]
  
  if (ncol(x) == 0) {
    stop("No nonconstant predictors remain.")
  }
  
  if (by_cell) {
    
    cell <- prepared$cell_training[keep_rows]
    x_means <- rowsum(x, cell)
    y_means <- rowsum(y, cell)[, 1]
    
    # Cell sizes in the row order that rowsum() used.
    cell_size <- as.numeric(table(cell)[rownames(x_means)])
    
    x_means <- x_means / cell_size
    y_means <- y_means / cell_size
    
    fit <- lm.fit(
      x - x_means[cell, , drop = FALSE],
      y - y_means[cell]
    )
    
    slopes <- fit$coefficients
    n_aliased <- sum(is.na(slopes))
    slopes[is.na(slopes)] <- 0
    
    cell_intercepts <- y_means - as.numeric(x_means %*% slopes)
    names(cell_intercepts) <- rownames(x_means)
    
    beta <- c("(Intercept)" = 0, slopes)
    
  } else {
    
    fit <- lm.fit(cbind("(Intercept)" = 1, x), y)
    
    beta <- fit$coefficients
    
    # A predictor that is an exact linear combination of earlier
    # ones is dropped, as in lm(); its coefficient is zero.
    n_aliased <- sum(is.na(beta))
    beta[is.na(beta)] <- 0
    
    cell_intercepts <- NULL
  }
  
  fit_details <- NULL
  
  if (details && !by_cell) {
    
    rank <- fit$rank
    estimated <- fit$qr$pivot[seq_len(rank)]
    
    rss <- sum(fit$residuals^2)
    df_residual <- length(y) - rank
    
    # (X'X)^{-1} from the triangular factor of the QR decomposition.
    unscaled <- chol2inv(
      qr.R(fit$qr)[seq_len(rank), seq_len(rank), drop = FALSE]
    )
    
    std_errors <- rep(NA_real_, length(beta))
    names(std_errors) <- names(beta)
    std_errors[estimated] <- sqrt(diag(unscaled) * rss / df_residual)
    
    fit_details <- list(
      std_errors = std_errors,
      sigma = sqrt(rss / df_residual),
      df_residual = df_residual,
      r_squared = 1 - rss / sum((y - mean(y))^2)
    )
  }
  
  if (any(!is.finite(beta)) || any(!is.finite(cell_intercepts))) {
    stop("Nonfinite coefficients.")
  }
  
  list(
    coefficients = beta,
    cell_intercepts = cell_intercepts,
    columns = colnames(x),
    n_training = length(y),
    n_aliased = n_aliased,
    details = fit_details
  )
}

# --------------------------------------------------
# 9. Sequential validation
# --------------------------------------------------

predict_sequence <- function(models, prepared, h) {
  
  n_origins <- length(prepared[[1]]$y_validation)
  
  predictions <- matrix(NA_real_, nrow = n_origins, ncol = 5)
  
  for (horizon in 1:5) {
    
    x <- prepared[[horizon]]$x_validation[
      ,
      models[[horizon]]$columns,
      drop = FALSE
    ]
    
    for (j in seq_len(h)) {
      
      column_name <- paste0("history_", j)
      
      if (!column_name %in% colnames(x)) {
        next
      }
      
      relative_day <- horizon - j
      
      if (relative_day > 0) {
        
        # Not yet observed at the origin: use the forecast.
        values <- predictions[, relative_day]
        
      } else {
        
        values <- prepared[[horizon]]$observed_history[
          ,
          1 - relative_day
        ]
      }
      
      x[, column_name] <- values
    }
    
    beta <- models[[horizon]]$coefficients
    cell_intercepts <- models[[horizon]]$cell_intercepts
    
    if (is.null(cell_intercepts)) {
      
      intercept <- beta[1]
      
    } else {
      
      intercept <- cell_intercepts[
        prepared[[horizon]]$cell_validation
      ]
      
      if (anyNA(intercept)) {
        stop("A validation airport-month has no training data.")
      }
    }
    
    predictions[, horizon] <- as.numeric(
      intercept + x %*% beta[-1]
    )
  }
  
  if (any(!is.finite(predictions))) {
    stop("Nonfinite validation predictions.")
  }
  
  predictions
}

# --------------------------------------------------
# 10. Cross-validation for one outcome
# --------------------------------------------------

# fold_years, history_values and seasonal_values restrict the work to
# part of the grid, so cluster array tasks can each fit one checkpoint.
checkpoint_path <- function(outcome, held_out_year, seasonal, h) {
  
  file.path(
    output_directory,
    paste0(
      outcome,
      "_year_", held_out_year,
      "_season_", seasonal,
      "_h_", h,
      ".rds"
    )
  )
}

run_outcome_cv <- function(outcome,
                           fold_years = years,
                           history_values = history_grid,
                           seasonal_values = seasonal_grid) {
  
  all_results <- list()
  result_index <- 1
  
  candidates <- expand.grid(
    h = history_values,
    seasonal = seasonal_values,
    stringsAsFactors = FALSE
  )
  
  for (held_out_year in fold_years) {
    
    checkpoint_files <- checkpoint_path(
      outcome,
      held_out_year,
      candidates$seasonal,
      candidates$h
    )
    
    # Avoid preparing an already completed fold.
    if (all(file.exists(checkpoint_files))) {
      
      for (file in checkpoint_files) {
        all_results[[result_index]] <- readRDS(file)
        result_index <- result_index + 1
      }
      
      next
    }
    
    cat(
      "\nPreparing:", outcome,
      "| Year:", held_out_year, "\n"
    )
    
    cleaned <- clean_temperature_outliers(
      weather_data,
      held_out_year
    )
    
    tables <- make_horizon_tables(cleaned, outcome)
    
    # The airport-specific columns are large; build them only
    # when a remaining candidate needs them.
    remaining <- candidates$seasonal[
      !file.exists(checkpoint_files)
    ]
    
    fold <- prepare_fold(
      tables,
      outcome,
      held_out_year,
      airport_season = any(
        remaining %in% c("airport1", "airport2")
      )
    )
    
    rm(cleaned, tables)
    
    cat(
      "Retained validation origins:",
      fold$origins_after,
      "of",
      fold$origins_before,
      "\n"
    )
    
    for (candidate in seq_len(nrow(candidates))) {
      
      h <- candidates$h[candidate]
      seasonal <- candidates$seasonal[candidate]
      checkpoint_file <- checkpoint_files[candidate]
      
      if (file.exists(checkpoint_file)) {
        
        all_results[[result_index]] <- readRDS(
          checkpoint_file
        )
        
        result_index <- result_index + 1
        next
      }
      
      cat("Fitting:", outcome, "| Year:", held_out_year,
          "| Season:", seasonal, "| h:", h, "\n")
      
      models <- lapply(
        fold$horizons,
        function(prepared) fit_ols(prepared, h, seasonal)
      )
      
      predictions <- predict_sequence(
        models,
        fold$horizons,
        h
      )
      
      rows <- lapply(1:5, function(horizon) {
        
        actual <- fold$horizons[[horizon]]$y_validation
        
        data.frame(
          outcome = outcome,
          year = held_out_year,
          seasonal = seasonal,
          h = h,
          horizon = horizon,
          sse = sum((actual - predictions[, horizon])^2),
          n = length(actual),
          n_training = models[[horizon]]$n_training,
          n_aliased = models[[horizon]]$n_aliased,
          origins_before = fold$origins_before,
          origins_after = fold$origins_after
        )
      })
      
      result <- dplyr::bind_rows(rows)
      
      temporary_file <- paste0(checkpoint_file, ".tmp")
      saveRDS(result, temporary_file)
      
      if (!file.rename(temporary_file, checkpoint_file)) {
        stop("Could not save checkpoint.")
      }
      
      all_results[[result_index]] <- result
      result_index <- result_index + 1
      
      rm(models, predictions)
    }
    
    rm(fold)
    invisible(gc())
  }
  
  result <- dplyr::bind_rows(all_results)
  
  # Only a full run produces the per-outcome summary file.
  if (identical(fold_years, years) &&
      identical(history_values, history_grid) &&
      identical(seasonal_values, seasonal_grid)) {
    
    saveRDS(
      result,
      file.path(
        output_directory,
        paste0("results_", outcome, ".rds")
      )
    )
  }
  
  result
}

# --------------------------------------------------
# 11. Run outcomes in parallel
# --------------------------------------------------

cat("Workers:", number_of_workers, "\n")
cat("Evaluation:", validation_origins, "\n")
cat("Checkpoint folder:", output_directory, "\n")

if (number_of_workers > 1) {
  future::plan(
    future::multisession,
    workers = number_of_workers
  )
} else {
  future::plan(future::sequential)
}

start_time <- Sys.time()

cv_results_list <- tryCatch(
  future.apply::future_lapply(
    outcomes,
    run_outcome_cv,
    future.packages = "dplyr",
    future.seed = TRUE
  ),
  finally = {
    future::plan(future::sequential)
  }
)

cv_results <- dplyr::bind_rows(cv_results_list)

saveRDS(
  cv_results,
  file.path(output_directory, "combined_cv_results.rds")
)

cat("Elapsed time:\n")
print(Sys.time() - start_time)

# --------------------------------------------------
# 12. Select seasonal structure and history length
# --------------------------------------------------

cv_summary <- cv_results %>%
  dplyr::mutate(
    seasonal = factor(seasonal, levels = seasonal_grid)
  ) %>%
  dplyr::group_by(outcome, seasonal, h) %>%
  dplyr::summarise(
    sse = sum(sse),
    n = sum(n),
    .groups = "drop"
  ) %>%
  dplyr::mutate(cv_mse = sse / n)

# Ties go to the simpler seasonal structure, then the shorter history.
best_settings <- cv_summary %>%
  dplyr::arrange(outcome, cv_mse, seasonal, h) %>%
  dplyr::group_by(outcome) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

print(best_settings)

write.csv(
  best_settings,
  file.path(output_directory, "best_settings.csv"),
  row.names = FALSE
)

saveRDS(
  best_settings,
  file.path(output_directory, "best_settings.rds")
)

write.csv(
  cv_summary,
  file.path(output_directory, "cv_summary.csv"),
  row.names = FALSE
)

# --------------------------------------------------
# 13. Plot CV curves
# --------------------------------------------------

seasonal_colors <- c(
  none = "#2a78d6",
  shared = "#eb6834",
  airport1 = "#1baf7a",
  airport2 = "#eda100",
  month = "#2a78d6",
  month_airport = "#eb6834"
)

cv_plot <- ggplot(
  cv_summary,
  aes(x = h, y = cv_mse, color = seasonal, shape = seasonal)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  geom_point(
    data = best_settings,
    color = "black",
    shape = 4,
    size = 3.5,
    stroke = 1.2
  ) +
  facet_wrap(~ outcome, scales = "free_y") +
  scale_x_continuous(breaks = history_grid) +
  scale_color_manual(values = seasonal_colors) +
  labs(
    title = "Cross-Validated Error by History Length and Seasonal Structure",
    subtitle = paste(
      "Complete, non-outlier validation sequences;",
      "cross marks the selected model"
    ),
    x = "History length h (days)",
    y = "Cross-validated MSE",
    color = "Seasonal structure",
    shape = "Seasonal structure"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(output_directory, "cv_curves.pdf"),
  plot = cv_plot,
  width = 11,
  height = 4.5
)

# The same curves for the two airport-specific structures only, so
# their small differences are visible.
airport_structures <- c("airport1", "airport2")

if (all(airport_structures %in% seasonal_grid)) {

  cv_plot_airport <- cv_plot %+%
    dplyr::filter(cv_summary, seasonal %in% airport_structures) +
    labs(
      title = paste(
        "Cross-Validated Error by History Length:",
        "Airport-Specific Seasonal Structures"
      )
    )

  # Show the selected model only if it is one of the two.
  cv_plot_airport$layers[[3]]$data <- dplyr::filter(
    best_settings,
    seasonal %in% airport_structures
  )

  ggsave(
    filename = file.path(output_directory, "cv_curves_airport.pdf"),
    plot = cv_plot_airport,
    width = 11,
    height = 4.5
  )
}

# --------------------------------------------------
# 14. Horizon performance and data retention
# --------------------------------------------------

selected_results <- cv_results %>%
  dplyr::inner_join(
    best_settings %>%
      dplyr::mutate(seasonal = as.character(seasonal)) %>%
      dplyr::select(outcome, seasonal, h),
    by = c("outcome", "seasonal", "h")
  )

horizon_performance <- selected_results %>%
  dplyr::group_by(outcome, horizon) %>%
  dplyr::summarise(
    sse = sum(sse),
    n = sum(n),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mse = sse / n,
    rmse = sqrt(mse)
  )

retention <- cv_results %>%
  dplyr::distinct(
    outcome, year, origins_before, origins_after
  ) %>%
  dplyr::mutate(
    fraction_retained = origins_after / origins_before
  )

write.csv(
  horizon_performance,
  file.path(output_directory, "horizon_performance.csv"),
  row.names = FALSE
)

write.csv(
  retention,
  file.path(output_directory, "validation_retention.csv"),
  row.names = FALSE
)

print(horizon_performance)
print(retention)

# --------------------------------------------------
# 15. Refit the selected models on all ten years
# --------------------------------------------------

fit_final_models <- function(outcome, seasonal, h) {
  
  weather_columns <- setdiff(numeric_weather, outcome)
  
  # Outlier cutoffs from all ten years.
  cleaned <- clean_temperature_outliers(weather_data, NA)
  tables <- make_horizon_tables(cleaned, outcome)
  rm(cleaned)
  
  airport_levels <- sort(unique(tables[[1]]$airport))
  
  event_levels <- sort(unique(
    tables[[1]]$Events[!is.na(tables[[1]]$Events)]
  ))
  
  numeric_columns <- c(
    history_names(max_history),
    weather_columns
  )
  
  models <- vector("list", 5)
  
  for (horizon in 1:5) {
    
    train <- tables[[horizon]]
    
    train <- train[complete_rows(
      train,
      c("response", "airport", "Events", weather_columns)
    ), ]
    
    x_training <- make_matrix(
      train,
      numeric_columns,
      airport_levels,
      event_levels,
      seasonal %in% c("airport1", "airport2")
    )
    
    models[[horizon]] <- fit_ols(
      list(
        x_training = x_training,
        y_training = train$response,
        cell_training = paste(train$airport, train$month)
      ),
      h,
      seasonal,
      details = TRUE
    )
    
    rm(train, x_training)
    invisible(gc())
  }
  
  list(
    outcome = outcome,
    seasonal = seasonal,
    h = h,
    airport_levels = airport_levels,
    event_levels = event_levels,
    horizons = models
  )
}

final_models <- lapply(
  seq_len(nrow(best_settings)),
  function(i) {
    
    cat(
      "Refitting on all years:", best_settings$outcome[i],
      "| Season:", as.character(best_settings$seasonal[i]),
      "| h:", best_settings$h[i], "\n"
    )
    
    fit_final_models(
      best_settings$outcome[i],
      as.character(best_settings$seasonal[i]),
      best_settings$h[i]
    )
  }
)

names(final_models) <- best_settings$outcome

saveRDS(
  final_models,
  file.path(output_directory, "final_models.rds")
)

cat("\nFinished. Results are in:", output_directory, "\n")
