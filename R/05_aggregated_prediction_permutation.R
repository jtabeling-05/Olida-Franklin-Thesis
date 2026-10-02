# Global comparison of Deleterious, Uncertain and Benign profiles.

# 1. Setup -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(vegan)
  library(permute)
})
source(file.path("R", "standard_statistical_tests.R"))

# 2. Construct the complete-prediction AP sample ---------------------------
# See standard_statistical_tests.R:
# Function 1: read_analysis_data().
# Function 2: prepare_endpoint_sample().

data <- read_analysis_data()
sample <- prepare_endpoint_sample(data$combinations, "aggregated_prediction")

# 3. Run the global profile and dispersion tests ---------------------------
# See standard_statistical_tests.R:
# Function 5: run_profile_analysis().
# Internally uses Function 3: profile_permanova()
# and Function 4: profile_dispersion().

result <- run_profile_analysis(
  sample,
  c("n_deleterious", "n_uncertain", "n_benign"),
  sample_name = "Primary",
  global_seed = 20260920L,
  dispersion_seed = 20260923L
)

# 4. Verify the expected primary sample size ---------------------------

stopifnot(nrow(sample) == 113L)

# 5. Export the profile results --------------------------------------------
# See standard_statistical_tests.R, Function 13: write_result().

write_result(result$test, "05_ap_profile_global_test.csv")
write_result(result$summary, "05_ap_profile_group_summary.csv")

# 6. Display the global test ------------------------------------------------

print(result$test)

