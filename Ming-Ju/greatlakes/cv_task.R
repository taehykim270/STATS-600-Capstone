# One cluster array task = one checkpoint
# (outcome, held-out year, seasonal structure, h).
#
# Usage, from the folder containing weather.rda and model_selection.R:
#   Rscript greatlakes/cv_task.R --count   # print the number of tasks
#   Rscript greatlakes/cv_task.R --init    # write run_settings.rds only
#   Rscript greatlakes/cv_task.R 17        # run task 17 (1-based)
#
# With no argument the task number is read from SLURM_ARRAY_TASK_ID.
#
# The settings, data preparation and model code all come from
# model_selection.R (sections 1-10), so there is one source of truth.

script <- readLines("model_selection.R", warn = FALSE)
cut <- grep("^# 11\\. Run outcomes in parallel", script)

if (length(cut) != 1) {
  stop("Could not find the section 11 header in model_selection.R.")
}

argument <- commandArgs(trailingOnly = TRUE)

if (length(argument) == 0) {
  argument <- Sys.getenv("SLURM_ARRAY_TASK_ID")
}

if (length(argument) != 1 || !nzchar(argument)) {
  stop("Give a task number, --count, or --init.")
}

definitions <- parse(
  text = script[seq_len(cut - 1)],
  keep.source = FALSE
)

if (argument == "--count") {

  # Only the settings are needed; skip loading the data.
  is_setting <- vapply(
    definitions,
    function(expression) {
      is.call(expression) &&
        identical(expression[[1]], as.name("<-")) &&
        as.character(expression[[2]])[1] %in%
        c("outcomes", "years", "history_grid", "seasonal_grid")
    },
    logical(1)
  )

  for (expression in definitions[is_setting]) {
    eval(expression, envir = globalenv())
  }

  cat(
    length(outcomes) * length(years) *
      length(history_grid) * length(seasonal_grid),
    "\n"
  )
  quit(save = "no")
}

suppressPackageStartupMessages(
  for (expression in definitions) {
    eval(expression, envir = globalenv())
  }
)

if (argument == "--init") {
  cat("Settings file ready:", settings_file, "\n")
  quit(save = "no")
}

# h varies fastest, so neighbouring tasks share a fold.
task_grid <- expand.grid(
  h = history_grid,
  seasonal = seasonal_grid,
  year = years,
  outcome = outcomes,
  stringsAsFactors = FALSE
)

task_id <- suppressWarnings(as.integer(argument))

if (is.na(task_id) || task_id < 1 || task_id > nrow(task_grid)) {
  stop("Task number must be between 1 and ", nrow(task_grid), ".")
}

task <- task_grid[task_id, ]

cat(
  "Task", task_id, "of", nrow(task_grid),
  "|", task$outcome,
  "| Year:", task$year,
  "| Season:", task$seasonal,
  "| h:", task$h, "\n"
)

start_time <- Sys.time()

result <- run_outcome_cv(
  task$outcome,
  fold_years = task$year,
  history_values = task$h,
  seasonal_values = task$seasonal
)

stopifnot(nrow(result) == 5)

cat("Task finished. Elapsed:\n")
print(Sys.time() - start_time)
