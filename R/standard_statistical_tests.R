# Shared functions used by the four primary analyses and their sensitivities.

effect_levels <- c(
  "True Digenic", "Monogenic+Modifier", "Dual Molecular Diagnosis"
)

# 1. Read the processed analysis datasets ----------------------------------

read_analysis_data <- function() {
  list(
    combinations = readr::read_csv(
      file.path("data", "processed", "combination_analysis.csv"),
      show_col_types = FALSE
    ),
    memberships = readr::read_csv(
      file.path("data", "processed", "membership_long.csv"),
      show_col_types = FALSE
    )
  )
}

# 2. Construct an endpoint-specific primary sample -------------------------

prepare_endpoint_sample <- function(combinations, endpoint) {
  x <- combinations |>
    dplyr::filter(Oligogenic_Effect %in% effect_levels, n_alleles >= 2L)

  if (endpoint == "genoox") {
    x <- x |>
      dplyr::filter(coverage == "Complete coverage") |>
      dplyr::mutate(
        n_success = n_PLP,
        n_failure = n_alleles - n_PLP
      )
  } else if (endpoint == "aggregated_prediction") {
    x <- x |>
      dplyr::mutate(
        n_predictions = n_deleterious + n_uncertain + n_benign
      ) |>
      dplyr::filter(n_predictions == n_alleles) |>
      dplyr::mutate(
        n_success = n_deleterious,
        n_failure = n_uncertain + n_benign
      )
  } else {
    stop("Unknown endpoint: ", endpoint)
  }

  x |>
    dplyr::mutate(
      Oligogenic_Effect = factor(Oligogenic_Effect, levels = effect_levels)
    )
}

# 3. Run the global permutation-based profile comparison -------------------

profile_permanova <- function(y, group, permutations = 9999L) {
  group <- droplevels(factor(group))
  indices <- t(replicate(permutations, sample.int(nrow(y))))
  # Invert label permutations to preserve the original randomization sequence.
  response_permutations <- t(apply(indices, 1, order))
  vegan::adonis2(
    y ~ group,
    method = "euclidean",
    permutations = response_permutations
  )
}

# 4. Test homogeneity of multivariate profile dispersion -------------------

profile_dispersion <- function(y, group, permutations, seed) {
  group <- droplevels(factor(group))
  set.seed(seed)
  permutation_matrix <- permute::shuffleSet(nrow(y), nset = permutations)
  fit <- vegan::betadisper(
    stats::dist(y), group, type = "median", bias.adjust = TRUE
  )
  test <- vegan::permutest(fit, permutations = permutation_matrix)
  list(
    F = unname(stats::anova(fit)[1, "F value"]),
    p = test$tab[1, "Pr(>F)"],
    distances = fit$distances
  )
}

# 5. Combine PERMANOVA, PERMDISP and descriptive profile summaries ---------

run_profile_analysis <- function(data, count_columns, sample_name,
                                 global_seed, dispersion_seed,
                                 permutations = 9999L) {
  y <- as.matrix(data[count_columns]) / data$n_alleles
  colnames(y) <- sub("^n_", "prop_", count_columns)
  group <- droplevels(data$Oligogenic_Effect)
  stopifnot(
    !anyDuplicated(data$OLIDA_ID),
    !anyNA(y),
    all(abs(rowSums(y) - 1) < 1e-12),
    nlevels(group) == 3L
  )

  set.seed(global_seed)
  global <- profile_permanova(y, group, permutations)
  dispersion <- profile_dispersion(y, group, permutations, dispersion_seed)

  profiles <- data.frame(
    sample = sample_name,
    OLIDA_ID = data$OLIDA_ID,
    Oligogenic_Effect = group,
    n_alleles = data$n_alleles,
    y,
    distance_to_spatial_median = dispersion$distances,
    check.names = FALSE
  )

  summary <- profiles |>
    dplyr::group_by(sample, Oligogenic_Effect) |>
    dplyr::summarise(
      n_combinations = dplyr::n(),
      dplyr::across(dplyr::starts_with("prop_"), mean),
      mean_distance = mean(distance_to_spatial_median),
      .groups = "drop"
    )

  test <- data.frame(
    sample = sample_name,
    n_combinations = nrow(data),
    pseudo_F = global$F[1],
    R_squared = global$R2[1],
    p_PERMANOVA = global$`Pr(>F)`[1],
    F_PERMDISP = dispersion$F,
    p_PERMDISP = dispersion$p,
    permutations = permutations,
    global_seed = global_seed,
    dispersion_seed = dispersion_seed
  )

  list(test = test, summary = summary, profiles = profiles)
}

# 6. Fit the grouped binomial regression -----------------------------------

fit_grouped_binomial <- function(data) {
  glmmTMB::glmmTMB(
    cbind(n_success, n_failure) ~ Oligogenic_Effect,
    family = stats::binomial(link = "logit"),
    data = data
  )
}

# 7. Calculate all pairwise odds ratios with Holm correction ---------------

effect_pairwise_holm <- function(model) {
  means <- emmeans::emmeans(model, ~ Oligogenic_Effect)
  contrasts <- emmeans::contrast(means, method = list(
    "Monogenic+Modifier vs True Digenic" = c(-1, 1, 0),
    "Dual Molecular Diagnosis vs True Digenic" = c(-1, 0, 1),
    "Monogenic+Modifier vs Dual Molecular Diagnosis" = c(0, 1, -1)
  ))
  estimates <- as.data.frame(summary(
    contrasts, infer = c(TRUE, TRUE), adjust = "none", type = "response"
  ))
  adjusted <- as.data.frame(summary(contrasts, adjust = "holm"))
  data.frame(
    contrast = estimates$contrast,
    OR = estimates$odds.ratio,
    lower = estimates$asymp.LCL,
    upper = estimates$asymp.UCL,
    p_raw = estimates$p.value,
    p_Holm = adjusted$p.value
  )
}

# 8. Summarize the global regression and pairwise comparisons --------------

summarize_binomial_fit <- function(data, sample_name) {
  full <- fit_grouped_binomial(data)
  null <- glmmTMB::glmmTMB(
    cbind(n_success, n_failure) ~ 1,
    family = stats::binomial(link = "logit"),
    data = data
  )
  stopifnot(
    full$fit$convergence == 0L, isTRUE(full$sdr$pdHess),
    null$fit$convergence == 0L, isTRUE(null$sdr$pdHess)
  )
  comparison <- stats::anova(null, full)

  list(
    model = full,
    global = data.frame(
      sample = sample_name, n_combinations = nrow(data),
      LR_chisq = comparison$Chisq[2], df = comparison$`Chi Df`[2],
      p = comparison$`Pr(>Chisq)`[2]
    ),
    pairwise = dplyr::mutate(
      effect_pairwise_holm(full), sample = sample_name, .before = 1
    )
  )
}

# 9. Evaluate residual deviance by parametric bootstrap --------------------

deviance_bootstrap <- function(data, simulations, seed) {
  model <- stats::glm(
    cbind(n_success, n_failure) ~ Oligogenic_Effect,
    family = stats::binomial(link = "logit"), data = data
  )
  observed <- stats::deviance(model)
  if (!isTRUE(model$converged) || !is.finite(observed)) {
    stop("The observed binomial model did not converge or has non-finite deviance.")
  }
  pearson_statistic <- sum(stats::residuals(model, type = "pearson")^2)
  residual_df <- stats::df.residual(model)
  probability <- stats::fitted(model)
  set.seed(seed)
  simulated <- numeric(simulations)
  simulated_data <- data
  for (i in seq_len(simulations)) {
    simulated_data$n_success <- stats::rbinom(nrow(data), data$n_alleles, probability)
    simulated_data$n_failure <- data$n_alleles - simulated_data$n_success
    refit <- stats::glm(
      cbind(n_success, n_failure) ~ Oligogenic_Effect,
      family = stats::binomial(link = "logit"), data = simulated_data
    )
    simulated[i] <- stats::deviance(refit)
    if (!isTRUE(refit$converged) || !is.finite(simulated[i])) {
      stop("Invalid binomial fit in bootstrap replication ", i, ".")
    }
  }
  data.frame(
    observed_deviance = observed,
    residual_df = residual_df,
    pearson_statistic = pearson_statistic,
    pearson_ratio = pearson_statistic / residual_df,
    pearson_asymptotic_p = stats::pchisq(pearson_statistic, residual_df, lower.tail = FALSE),
    asymptotic_p = stats::pchisq(observed, stats::df.residual(model), lower.tail = FALSE),
    bootstrap_p = (1 + sum(simulated >= observed)) / (simulations + 1),
    simulated_q025 = unname(stats::quantile(simulated, .025)),
    simulated_q975 = unname(stats::quantile(simulated, .975)),
    simulations = simulations,
    seed = seed
  )
}

# 10. Compare binomial and beta-binomial model families --------------------

compare_model_families <- function(data) {
  binomial <- fit_grouped_binomial(data)
  beta_binomial <- glmmTMB::glmmTMB(
    cbind(n_success, n_failure) ~ Oligogenic_Effect,
    family = glmmTMB::betabinomial(link = "logit"), data = data
  )
  for (model in list(binomial = binomial, beta_binomial = beta_binomial)) {
    if (!isTRUE(model$fit$convergence == 0L) || !isTRUE(model$sdr$pdHess) ||
        !is.finite(stats::AIC(model))) {
      stop("Model-family comparison requires converged fits with a positive-definite Hessian and finite AIC.")
    }
  }
  data.frame(
    model = c("Binomial", "Beta-binomial"),
    AIC = c(stats::AIC(binomial), stats::AIC(beta_binomial)),
    convergence = c(binomial$fit$convergence, beta_binomial$fit$convergence),
    positive_definite_hessian = c(
      isTRUE(binomial$sdr$pdHess), isTRUE(beta_binomial$sdr$pdHess)
    )
  ) |>
    dplyr::mutate(delta_AIC_vs_binomial = AIC - AIC[1])
}

# 11. Check whether observed boundary profiles are binomially plausible ----

boundary_simulation <- function(data, simulations = 5000L, seed = 20260908L) {
  model <- fit_grouped_binomial(data)
  set.seed(seed)
  simulated <- stats::simulate(model, nsim = simulations)
  totals <- vapply(simulated, function(x) {
    success <- x[, 1]
    sum(success == 0L | success == data$n_alleles)
  }, numeric(1))
  observed <- sum(data$n_success == 0L | data$n_success == data$n_alleles)
  data.frame(
    observed = observed,
    simulated_mean = mean(totals),
    simulated_q025 = unname(stats::quantile(totals, .025)),
    simulated_q975 = unname(stats::quantile(totals, .975)),
    simulations = simulations,
    seed = seed
  )
}

# 12. Remove combinations containing recurrent alleles --------------------

recurrent_allele_free <- function(sample, memberships) {
  members <- memberships |>
    dplyr::semi_join(sample, by = "OLIDA_ID") |>
    dplyr::filter(!is.na(Allele_ID))
  repeated <- members |>
    dplyr::count(Allele_ID) |>
    dplyr::filter(n > 1L)
  affected <- members |>
    dplyr::semi_join(repeated, by = "Allele_ID") |>
    dplyr::distinct(OLIDA_ID)
  dplyr::anti_join(sample, affected, by = "OLIDA_ID")
}

# 13. Write a result table to the standard output directory ----------------

write_result <- function(x, filename) {
  dir.create(file.path("output", "tables"), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file.path("output", "tables", filename), na = "")
}

