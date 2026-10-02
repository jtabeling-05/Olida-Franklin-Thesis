# Two structural restrictions for both outputs and an exploratory Genoox
# regression restricted to the complete AP sample (evaluated separately).

# 1. Setup -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(vegan)
  library(permute)
  library(glmmTMB)
  library(emmeans)
})
source(file.path("R", "standard_statistical_tests.R"))

# 2. Construct the three independently evaluated sample restrictions -------
# See standard_statistical_tests.R:
# Function 1: read_analysis_data().
# Function 2: prepare_endpoint_sample().
# Function 12: recurrent_allele_free().

data <- read_analysis_data()
genoox_primary <- prepare_endpoint_sample(data$combinations, "genoox")
ap_primary <- prepare_endpoint_sample(data$combinations, "aggregated_prediction")

make_samples <- function(primary) {
  list(
    `Exactly two alleles` = filter(primary, n_alleles == 2L),
    `Recurrent-allele-free` = recurrent_allele_free(primary, data$memberships)
  )
}

genoox <- make_samples(genoox_primary)
ap <- make_samples(ap_primary)

# Keep the Genoox outcome and denominators; use AP only to select OLIDA IDs.
common_sample <- semi_join(genoox_primary, ap_primary, by = "OLIDA_ID")
stopifnot(
  !anyDuplicated(common_sample$OLIDA_ID),
  setequal(common_sample$OLIDA_ID, ap_primary$OLIDA_ID),
  nrow(common_sample) == 113L,
  all(common_sample$n_success == common_sample$n_PLP),
  all(common_sample$n_success + common_sample$n_failure == common_sample$n_alleles)
)

# 3. Test profiles in the two-allele and recurrent-allele-free samples -------
# The AP-complete Genoox restriction addresses regression results only.
# See standard_statistical_tests.R:
# Function 5: run_profile_analysis().
# Internally uses Function 3: profile_permanova()
# and Function 4: profile_dispersion().

profile_specs <- list(
  list("Genoox", genoox[[1]], c("n_PLP", "n_possibly_pathogenic", "n_VUS", "n_possibly_benign", "n_BLB"), 20260828L, 20260918L),
  list("Genoox", genoox[[2]], c("n_PLP", "n_possibly_pathogenic", "n_VUS", "n_possibly_benign", "n_BLB"), 20260920L, 20260921L),
  list("Aggregated Prediction", ap[[1]], c("n_deleterious", "n_uncertain", "n_benign"), 20260921L, 20260924L),
  list("Aggregated Prediction", ap[[2]], c("n_deleterious", "n_uncertain", "n_benign"), 20260922L, 20260925L)
)
sample_names <- rep(c("Exactly two alleles", "Recurrent-allele-free"), 2)

profiles <- lapply(seq_along(profile_specs), function(i) {
  z <- profile_specs[[i]]
  result <- run_profile_analysis(z[[2]], z[[3]], sample_names[i], z[[4]], z[[5]])
  result$test |>
    mutate(endpoint = z[[1]], .before = 1)
}) |>
  bind_rows()

# 4. Fit regressions and all pairwise contrasts for every sensitivity sample -
# See standard_statistical_tests.R, Function 8: summarize_binomial_fit().
# Internally uses Function 6: fit_grouped_binomial()
# and Function 7: effect_pairwise_holm().
# Each call below fits the full and intercept-only models, verifies convergence
# and positive-definite Hessians, then returns the LR test and Holm contrasts.
# The primary model specification is retained for every restricted sample.

regression_specs <- list(
  list("Genoox", genoox[[1]], "Exactly two alleles"),
  list("Genoox", genoox[[2]], "Recurrent-allele-free"),
  list("Aggregated Prediction", ap[[1]], "Exactly two alleles"),
  list("Aggregated Prediction", ap[[2]], "Recurrent-allele-free"),
  list("Genoox", common_sample, "AP-complete sample")
)
regressions <- lapply(regression_specs, function(z) {
  result <- summarize_binomial_fit(z[[2]], z[[3]])
  list(
    global = mutate(result$global, endpoint = z[[1]], .before = 1),
    pairwise = mutate(result$pairwise, endpoint = z[[1]], .before = 1)
  )
})

global <- bind_rows(lapply(regressions, `[[`, "global"))
pairwise <- bind_rows(lapply(regressions, `[[`, "pairwise"))

# 5. Combine all sensitivity samples in one concise overview ---------------
# Start from regression results so the AP-complete sample is retained.
# NA profile p-values mean that no profile tests were performed for that row.

sensitivity_summary <- tibble::as_tibble(global) |>
  select(endpoint, sample, n_combinations, regression_p = p) |>
  left_join(
    profiles |> select(endpoint, sample, p_PERMANOVA, p_PERMDISP),
    by = c("endpoint", "sample")
  ) |>
  left_join(
    pairwise |>
      filter(contrast == "Monogenic+Modifier vs True Digenic") |>
      select(endpoint, sample, MM_vs_TD_OR = OR, MM_vs_TD_p_Holm = p_Holm),
    by = c("endpoint", "sample")
  ) |>
  left_join(
    pairwise |>
      filter(contrast == "Monogenic+Modifier vs Dual Molecular Diagnosis") |>
      select(endpoint, sample, MM_vs_DMD_OR = OR, MM_vs_DMD_p_Holm = p_Holm),
    by = c("endpoint", "sample")
  ) |>
  relocate(p_PERMANOVA, p_PERMDISP, .before = regression_p) |>
  arrange(match(endpoint, c("Genoox", "Aggregated Prediction")),
          match(sample, c("Exactly two alleles", "Recurrent-allele-free", "AP-complete sample")))

# 6. Export profile tests, global regression tests and pairwise contrasts ---
# See standard_statistical_tests.R, Function 13: write_result().

write_result(profiles, "07_sensitivity_profile_tests.csv")
write_result(global, "07_sensitivity_regression_tests.csv")
write_result(pairwise, "07_sensitivity_pairwise_holm.csv")

# 7. Display the complete sensitivity overview -----------------------------

print(sensitivity_summary, width = Inf)

