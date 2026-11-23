# Install any missing packages into the personal R library.
# Run on a Great Lakes login node (it has internet access).

packages <- c(
  "dplyr", "tidyr", "tibble", "ggplot2",
  "future", "future.apply"
)

library_path <- Sys.getenv("R_LIBS_USER")

dir.create(library_path, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(library_path, .libPaths()))

missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {

  cat("Installing:", paste(missing, collapse = ", "), "\n")

  install.packages(
    missing,
    lib = library_path,
    repos = "https://cloud.r-project.org",
    Ncpus = 4
  )
}

still_missing <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(still_missing) > 0) {
  stop(
    "Could not install: ",
    paste(still_missing, collapse = ", ")
  )
}

cat(R.version.string, "\n")

for (package in packages) {
  cat(package, as.character(packageVersion(package)), "\n")
}
