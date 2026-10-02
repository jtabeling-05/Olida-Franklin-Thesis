# Grouped binomial regression of Deleterious counts per OLIDA combination.

# 1. Setup -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glmmTMB)
  library(emmeans)
})
source(file.path("R", "standard_statistical_tests.R"))

# 2. Construct the complete-prediction AP sample ---------------------------
# See standard_statistical_tests.R:
# Function 1: read_analysis_data().
# Function 2: prepare_endpoint_sample().

data <- read_analysis_data()
sample <- prepare_endpoint_sample(data$combinations, "aggregated_prediction")

# 3. Fit the grouped binomial regression -----------------------------------
# See standard_statistical_tests.R, Function 8: summarize_binomial_fit().
# Internally uses Function 6: fit_grouped_binomial()
# and Function 7: effect_pairwise_holm().

fit <- summarize_binomial_fit(sample, "Primary")

# 4. Assess residual deviance by parametric bootstrap ----------------------
# See standard_statistical_tests.R, Function 9: deviance_bootstrap().

diagnostic <- deviance_bootstrap(sample, simulations = 4999L, seed = 20260913L)

# 5. Compare binomial and beta-binomial model families ---------------------
# See standard_statistical_tests.R, Function 10: compare_model_families().

families <- compare_model_families(sample)

# 6. Verify the expected primary sample size ---------------------------

stopifnot(nrow(sample) == 113L)

# 7. Export regression and diagnostic results ------------------------------
# See standard_statistical_tests.R, Function 13: write_result().

write_result(fit$global, "06_ap_regression_global_test.csv")
write_result(fit$pairwise, "06_ap_regression_pairwise_holm.csv")
write_result(diagnostic, "06_ap_regression_fit.csv")
write_result(families, "06_ap_regression_model_comparison.csv")

# 8. Display global and pairwise results -----------------------------------

print(fit$global)
print(fit$pairwise)

