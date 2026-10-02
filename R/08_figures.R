# Figures reported in the thesis.

# 1. Setup, data and consistent visual style -------------------------------

library(tidyverse)
library(readxl)

dir.create(file.path("output", "figures"), recursive = TRUE, showWarnings = FALSE)

effect_levels_all <- c(
  "True Digenic", "Monogenic+Modifier", "Dual Molecular Diagnosis", "Unknown"
)
effect_colours <- c(
  "True Digenic" = "#F8766D", "Monogenic+Modifier" = "#7CAE00",
  "Dual Molecular Diagnosis" = "#00BFC4", "Unknown" = "#C77CFF"
)

combinations <- readr::read_csv(
  file.path("data", "processed", "combination_analysis.csv"),
  show_col_types = FALSE
)

figure_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(colour = "#222222"),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# 2. Define reusable figure functions --------------------------------------

group_labels <- function(data, levels = effect_levels_all) {
  sizes <- data |>
    distinct(OLIDA_ID, Oligogenic_Effect) |>
    count(Oligogenic_Effect)
  setNames(
    paste0(
      c("True\nDigenic", "Monogenic\n+ Modifier", "Dual\nMolecular\nDiagnosis", "Unknown")[seq_along(levels)],
      "\n(n = ", sizes$n[match(levels, sizes$Oligogenic_Effect)], ")"
    ),
    levels
  )
}

save_png <- function(plot, filename, width, height) {
  ggsave(
    file.path("output", "figures", filename), plot,
    width = width, height = height, dpi = 300, bg = "white"
  )
}

profile_plot <- function(data, count_columns, category_labels) {
  long <- data |>
    select(OLIDA_ID, Oligogenic_Effect, n_alleles, all_of(count_columns)) |>
    pivot_longer(all_of(count_columns), names_to = "category", values_to = "count") |>
    mutate(
      category = factor(category, levels = count_columns, labels = category_labels),
      proportion = count / n_alleles,
      Oligogenic_Effect = factor(Oligogenic_Effect, levels = effect_levels_all)
    )

  stopifnot(
    !anyDuplicated(long[c("OLIDA_ID", "category")]),
    all(abs((long |> group_by(OLIDA_ID) |> summarise(x = sum(proportion)))$x - 1) < 1e-10)
  )

  frequencies <- long |>
    count(Oligogenic_Effect, category, proportion, name = "combinations")

  ggplot(frequencies, aes(Oligogenic_Effect, proportion)) +
    geom_point(
      aes(size = combinations, fill = Oligogenic_Effect),
      shape = 21, colour = "#333333", alpha = .75, stroke = .25
    ) +
    facet_wrap(~category, ncol = 2) +
    scale_fill_manual(values = effect_colours, guide = "none") +
    scale_size_area(max_size = 10, name = "Combinations") +
    scale_x_discrete(labels = group_labels(data)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, .25)) +
    labs(x = "Oligogenic effect", y = "Category proportion per combination") +
    figure_theme
}

endpoint_boxplot <- function(data, proportion, y_label) {
  defined_levels <- effect_levels_all[1:3]
  x <- data |>
    mutate(Oligogenic_Effect = factor(Oligogenic_Effect, levels = defined_levels))

  ggplot(x, aes(Oligogenic_Effect, {{ proportion }}, fill = Oligogenic_Effect)) +
    geom_boxplot(width = .58, outlier.shape = NA, alpha = .8) +
    geom_point(
      position = position_jitter(width = .14, height = 0, seed = 18),
      size = 1.5, alpha = .45
    ) +
    scale_fill_manual(values = effect_colours) +
    scale_x_discrete(labels = group_labels(x, defined_levels)) +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, .25)) +
    labs(x = "Oligogenic effect", y = y_label) +
    figure_theme +
    theme(legend.position = "none")
}

# 3. Figure 1: matching workflow -------------------------------------------
submitted <- nrow(readxl::read_excel("data/raw/SMALLVARIANT.xlsx", skip = 1))
resolution <- readr::read_csv(
  "data/processed/franklin_match_resolution.csv", show_col_types = FALSE
)
stopifnot(
  !anyNA(resolution$match_accepted_primary),
  is.logical(resolution$match_accepted_primary)
)
matching_counts <- resolution |> count(resolution_method) |> deframe()
exported <- nrow(resolution)
matched <- sum(resolution$match_accepted_primary)
excluded <- exported - matched

matching_totals <- tibble(
  stage = c("Submitted", "Exported", "Matched", "Excluded"),
  total = c(submitted, exported, matched, excluded),
  y = c(4, 3, 2, 1)
)
matching_segments <- tribble(
  ~stage, ~category, ~count, ~colour,
  "Submitted", "Exported", exported, "#00BFC4",
  "Submitted", "Not returned", submitted - exported, "#DEDEDE",
  "Exported", "Matched", matched, "#7CAE00",
  "Exported", "Excluded", excluded, "#F8766D",
  "Matched", "Single exact candidate",
  matching_counts[["unique_exact_genomic_allele"]], "#7CAE00",
  "Matched", "cDNA resolution",
  matching_counts[["resolved_by_exact_cdna_within_allele"]], "#BED780",
  "Matched", "Protein resolution",
  matching_counts[["resolved_by_exact_protein_within_allele"]], "#E5EFCC",
  "Excluded", "Indels",
  matching_counts[["normalization_required"]], "#F8766D",
  "Excluded", "Unresolved multiple IDs",
  matching_counts[["exact_allele_multiple_variant_ids"]], "#FBC5C1"
) |> left_join(matching_totals, by = "stage")

# Each bar must account for its complete parent category.
matching_checks <- matching_segments |>
  group_by(stage) |>
  summarise(accounted_for = sum(count), total = first(total), .groups = "drop")
stopifnot(
  submitted >= exported, !anyNA(matching_segments$count),
  all(matching_segments$count >= 0),
  all(matching_checks$accounted_for == matching_checks$total),
  all(matching_checks$total > 0)
)
matching_segments <- matching_segments |>
  group_by(stage) |>
  mutate(
    proportion = count / total,
    xmax = cumsum(proportion), xmin = xmax - proportion,
    label_y = y + (n() + 1 - 2 * row_number()) * 0.13,
    label = paste0(category, ": ", scales::comma(count, accuracy = 1),
                   " (", sprintf("%.2f", 100 * proportion), "%)")
  ) |> ungroup()
matching_ticks <- tidyr::crossing(y = matching_totals$y, x = seq(0, 1, by = 0.25))
matching_plot <- ggplot() +
  geom_rect(
    data = matching_segments,
    aes(xmin = xmin, xmax = xmax, ymin = y - 0.20, ymax = y + 0.20, fill = colour)
  ) +
  geom_segment(
    data = matching_ticks,
    aes(x = x, xend = x, y = y - 0.23, yend = y - 0.28),
    colour = "#999999", linewidth = 0.25
  ) +
  geom_text(
    data = matching_totals,
    aes(x = -0.45, y = y,
        label = paste0(stage, "\n", scales::comma(total, accuracy = 1))),
    hjust = 0, size = 3.5, lineheight = 1.05
  ) +
  geom_point(
    data = matching_segments, aes(x = 1.10, y = label_y, fill = colour),
    shape = 22, size = 2.5, colour = "#777777", stroke = 0.3
  ) +
  geom_text(
    data = matching_segments, aes(x = 1.15, y = label_y, label = label),
    hjust = 0, size = 3.5
  ) +
  annotate("text", x = c(0, 0.5, 1), y = 0.49,
           label = c("0%", "50%", "100%"), size = 3.5, colour = "#444444") +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-0.47, 2.03), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(0.35, 4.35), expand = expansion(mult = 0)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(4, 4, 4, 4)
  )
# 4. Figures 2 and 3: Genoox Classification -------------------------------
genoox <- combinations |>
  filter(
    Oligogenic_Effect %in% effect_levels_all,
    n_alleles >= 2,
    coverage == "Complete coverage"
  )
stopifnot(sum(genoox$Oligogenic_Effect != "Unknown") == 147L)

genoox_profile <- profile_plot(
  genoox,
  c("n_PLP", "n_possibly_pathogenic", "n_VUS", "n_possibly_benign", "n_BLB"),
  c("P/LP", "Possibly pathogenic", "VUS", "Possibly benign", "B/LB")
)
genoox_box <- endpoint_boxplot(
  filter(genoox, Oligogenic_Effect != "Unknown"),
  prop_PLP,
  "Proportion of P/LP-classified alleles"
)


# 5. Figures 4 and 5: Aggregated Prediction -------------------------------
ap <- combinations |>
  mutate(n_predictions = n_deleterious + n_uncertain + n_benign) |>
  filter(
    Oligogenic_Effect %in% effect_levels_all,
    n_alleles >= 2,
    n_predictions == n_alleles
  )
stopifnot(sum(ap$Oligogenic_Effect != "Unknown") == 113L)

ap_profile <- profile_plot(
  ap,
  c("n_deleterious", "n_uncertain", "n_benign"),
  c("Deleterious", "Uncertain", "Benign")
)
ap_box <- endpoint_boxplot(
  filter(ap, Oligogenic_Effect != "Unknown"),
  prop_deleterious,
  "Proportion of deleterious alleles"
)


# 6. Export the validated figures
save_png(matching_plot, "01_matching_workflow.png", 10, 4)
save_png(genoox_profile, "02_genoox_profiles.png", 10, 10.5)
save_png(genoox_box, "03_genoox_plp_proportions.png", 8, 5.5)
save_png(ap_profile, "04_aggregated_prediction_profiles.png", 10, 7.5)
save_png(ap_box, "05_ap_deleterious_proportions.png", 8, 5.5)

# 7. Report completion ------------------------------------------------------

cat("Five thesis figures created in output/figures.\n")

