# Run the complete analysis from the repository root.

# 1. Define the reproducible script order ----------------------------------

scripts <- c(
  "01_prepare_data.R",
  "02_descriptive_statistics.R",
  "03_genoox_permutation.R",
  "04_genoox_regression.R",
  "05_aggregated_prediction_permutation.R",
  "06_aggregated_prediction_regression.R",
  "07_sensitivity_analyses.R",
  "08_figures.R"
)

# 2. Run every script in an isolated environment ---------------------------

for (script in scripts) {
  message("Running R/", script)
  source(file.path("R", script), local = new.env(parent = globalenv()))
}

# 3. Record the software versions used for this completed run ---------------
dir.create("output", showWarnings = FALSE)
writeLines(trimws(capture.output(sessionInfo()), which = "right"),
           file.path("output", "sessionInfo.txt"))


