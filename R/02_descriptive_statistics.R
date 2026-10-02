# Descriptive overview of the processed OLIDA–Franklin data.

# 1. Setup -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
source(file.path("R", "standard_statistical_tests.R"))

# 2. Read the processed analysis datasets ---------------------------------
# See standard_statistical_tests.R, Function 1: read_analysis_data().

data <- read_analysis_data()
combinations <- data$combinations

# 3. Describe all OLIDA combinations ---------------------------------------

# Keep Unknown for description; inferential scripts select the three defined effects.
combination_counts <- combinations |>
  count(Oligogenic_Effect, name = "n_combinations") |>
  arrange(desc(n_combinations))

# 4. Construct the two primary analytical samples --------------------------
# See standard_statistical_tests.R, Function 2: prepare_endpoint_sample().

genoox_sample <- prepare_endpoint_sample(combinations, "genoox")
ap_sample <- prepare_endpoint_sample(combinations, "aggregated_prediction")

# 5. Summarize combinations, alleles and binary outcomes -------------------

analysis_samples <- bind_rows(
  genoox_sample |>
    group_by(Oligogenic_Effect) |>
    summarise(
      n_combinations = n(), n_alleles = sum(n_alleles),
      positive = sum(n_success), negative = sum(n_failure), .groups = "drop"
    ) |>
    mutate(endpoint = "Genoox P/LP", .before = 1),
  ap_sample |>
    group_by(Oligogenic_Effect) |>
    summarise(
      n_combinations = n(), n_alleles = sum(n_alleles),
      positive = sum(n_success), negative = sum(n_failure), .groups = "drop"
    ) |>
    mutate(endpoint = "Aggregated Prediction: Deleterious", .before = 1)
)

# 6. Verify the expected sample sizes --------------------------------------

stopifnot(
  nrow(genoox_sample) == 147L,
  nrow(ap_sample) == 113L,
  sum(analysis_samples$n_combinations[analysis_samples$endpoint == "Genoox P/LP"]) == 147L,
  sum(analysis_samples$n_combinations[analysis_samples$endpoint == "Aggregated Prediction: Deleterious"]) == 113L
)

# 7. Export the descriptive tables -----------------------------------------
# See standard_statistical_tests.R, Function 13: write_result().

write_result(combination_counts, "02_all_combination_counts.csv")
write_result(analysis_samples, "02_primary_sample_counts.csv")

# 8. Report completion ------------------------------------------------------

cat("Descriptive tables written. Genoox n = 147; Aggregated Prediction n = 113.\n")

