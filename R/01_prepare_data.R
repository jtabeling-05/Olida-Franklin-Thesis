# ==============================================================================
# 01_prepare_data.R
#
# Purpose: Harmonize OLIDA and Franklin data, perform allele-specific matching,
# construct traceable membership datasets, and create combination-level data.
#
# Raw inputs:
#   data/raw/Combination.xlsx
#   data/raw/SMALLVARIANT.xlsx
#   data/raw/Small Variants Franklin complete.csv
# ==============================================================================

# 0. Setup, paths and reusable cleaning functions --------------------------

library(readxl)
library(tidyverse)
library(writexl)

# Canonical OLIDA-Franklin pipeline.
raw_dir <- file.path("data", "raw")
processed_dir <- file.path("data", "processed")
supplementary_matching_dir <- file.path("supplementary", "matching")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplementary_matching_dir, recursive = TRUE, showWarnings = FALSE)

clean_missing <- function(x) {
  x <- str_trim(as.character(x))
  x[str_to_upper(x) %in% c("", "NA", "N.A.", "NONE")] <- NA_character_
  x
}
clean_chr <- function(x) str_remove(str_to_upper(clean_missing(x)), "^CHR")
clean_position <- function(x) {
  position <- as.numeric(clean_missing(x))
  if (any(!is.na(position) &
          (!is.finite(position) | position < 1 | position != floor(position)))) {
    stop("Genomic positions must be positive, finite whole numbers.")
  }
  if (any(!is.na(clean_missing(x)) & is.na(position))) {
    stop("Unrecognized genomic position.")
  }
  ifelse(is.na(position), NA_character_,
         format(position, scientific = FALSE, trim = TRUE, digits = 15))
}
clean_allele <- function(x) str_to_upper(clean_missing(x))
clean_hgvs <- function(x) str_replace_all(str_remove_all(str_to_upper(clean_missing(x)), "\\s+"), "TER", "*")
single_value <- function(x) {
  x <- unique(na.omit(as.character(x)))
  if (length(x) == 1) x else NA_character_
}
make_allele_id <- function(chr, pos, ref, alt) {
  if_else(!is.na(chr) & !is.na(pos) & !is.na(ref) & !is.na(alt),
          str_c(chr, pos, ref, alt, sep = ":"), NA_character_)
}
classify_genoox <- function(x) {
  case_when(
    x %in% c("PATHOGENIC", "LIKELY_PATHOGENIC") ~ "P/LP",
    x %in% c("POSSIBLY_PATHOGENIC_LOW", "POSSIBLY_PATHOGENIC_MODERATE") ~ "Possibly pathogenic",
    x == "UNCERTAIN_SIGNIFICANCE" ~ "VUS",
    x == "POSSIBLY_BENIGN" ~ "Possibly benign",
    x %in% c("BENIGN", "LIKELY_BENIGN") ~ "B/LB",
    TRUE ~ NA_character_
  )
}

franklin_raw <- read.delim(file.path(raw_dir, "Small Variants Franklin complete.csv"), sep = ";", check.names = FALSE)
small_raw <- read_excel(file.path(raw_dir, "SMALLVARIANT.xlsx"), skip = 1)
combination_raw <- read_excel(file.path(raw_dir, "Combination.xlsx"), skip = 1)
expected_genoox_classes <- c(
  "PATHOGENIC", "LIKELY_PATHOGENIC",
  "POSSIBLY_PATHOGENIC_LOW", "POSSIBLY_PATHOGENIC_MODERATE",
  "UNCERTAIN_SIGNIFICANCE", "POSSIBLY_BENIGN", "BENIGN", "LIKELY_BENIGN"
)
stopifnot(
  setequal(unique(na.omit(franklin_raw$`Genoox Classification`)), expected_genoox_classes),
  !anyNA(classify_genoox(franklin_raw$`Genoox Classification`))
)

# 1. Create unique OLIDA variants and complete allele keys -----------------
variant_master_base <- small_raw %>%
  transmute(
    Variant_ID = as.character(`Entry Id`),
    Gene = str_to_upper(clean_missing(Gene)), Chr = clean_chr(Chromosome),
    Position_hg38 = clean_position(`Genomic Position Hg38`),
    Ref = clean_allele(`Ref Allele`), Alt = clean_allele(`Alt Allele`),
    Cdna_Change = clean_missing(`Cdna Change`), Protein_Change = clean_missing(`Protein Change`),
    Transcript_ID = clean_missing(`Transcript Id`), dbSNP_ID = clean_missing(`Dbsnp Id`),
    Cdna_match = clean_hgvs(`Cdna Change`), Protein_match = clean_hgvs(`Protein Change`)
  ) %>%
  mutate(
    has_complete_allele_key = !is.na(Chr) & !is.na(Position_hg38) & !is.na(Ref) & !is.na(Alt),
    Allele_ID = make_allele_id(Chr, Position_hg38, Ref, Alt)
  )
stopifnot(nrow(variant_master_base) == n_distinct(variant_master_base$Variant_ID))

# 2. Prepare Franklin records and allele keys ------------------------------
# Missing key components are never eligible for matching.
franklin_records <- franklin_raw %>%
  mutate(
    Franklin_record_ID = sprintf("FR%05d", row_number()),
    Gene_match = str_to_upper(clean_missing(Gene)), Chr_match = clean_chr(Chr),
    Position_match = clean_position(`Start Position`),
    Ref_match = clean_allele(Ref), Alt_match = clean_allele(Alt),
    Cdna_match = clean_hgvs(Nucleotide), Protein_match = clean_hgvs(`AA Change`),
    has_complete_allele_key = !is.na(Chr_match) & !is.na(Position_match) &
      !is.na(Ref_match) & !is.na(Alt_match),
    Allele_ID = make_allele_id(Chr_match, Position_match, Ref_match, Alt_match),
    .before = 1
  )
stopifnot(nrow(franklin_records) == n_distinct(franklin_records$Franklin_record_ID))

# 3. Identify exact genomic-allele candidates ------------------------------
# Matching key: chromosome + hg38 position + REF + ALT.
# na_matches="never" prevents absent allele values from matching one another.
franklin_match_candidates <- franklin_records %>%
  filter(has_complete_allele_key) %>%
  select(Franklin_record_ID, Gene_franklin = Gene_match, Chr = Chr_match,
         Position_hg38 = Position_match, Ref = Ref_match, Alt = Alt_match,
         Allele_ID, Cdna_franklin = Cdna_match, Protein_franklin = Protein_match) %>%
  inner_join(
    variant_master_base %>% filter(has_complete_allele_key) %>%
      select(Variant_ID, Gene_OLIDA = Gene, Chr, Position_hg38, Ref, Alt,
             Cdna_OLIDA = Cdna_match, Protein_OLIDA = Protein_match),
    by = c("Chr", "Position_hg38", "Ref", "Alt"),
    relationship = "many-to-many", na_matches = "never"
  ) %>%
  mutate(
    Gene_concordant = !is.na(Gene_franklin) & Gene_franklin == Gene_OLIDA,
    Cdna_exact = !is.na(Cdna_franklin) & !is.na(Cdna_OLIDA) & Cdna_franklin == Cdna_OLIDA,
    Protein_exact = !is.na(Protein_franklin) & !is.na(Protein_OLIDA) &
      Protein_franklin == Protein_OLIDA
  ) %>%
  distinct(Franklin_record_ID, Variant_ID, .keep_all = TRUE)

# 4. Resolve exact allele matches and multiple-ID cases --------------------
# Accept unique exact-allele matches directly.
# Resolve multiple Variant_ID candidates by unique exact cDNA agreement,
# followed by unique exact protein agreement.
# Unresolved multiple-ID cases are excluded from the primary matching.
candidate_resolution <- franklin_match_candidates %>%
  group_by(Franklin_record_ID) %>%
  summarise(
    n_exact_allele_candidates = n_distinct(Variant_ID),
    n_cdna_supported = n_distinct(Variant_ID[Cdna_exact]),
    n_protein_supported = n_distinct(Variant_ID[Protein_exact]),
    sole_candidate = if_else(n_exact_allele_candidates == 1L, first(Variant_ID), NA_character_),
    cdna_candidate = if_else(n_cdna_supported == 1L, first(Variant_ID[Cdna_exact]), NA_character_),
    protein_candidate = if_else(n_protein_supported == 1L, first(Variant_ID[Protein_exact]), NA_character_),
    .groups = "drop"
  )

franklin_match_resolution <- franklin_records %>%
  select(Franklin_record_ID, `Variation Type`, has_complete_allele_key, Allele_ID) %>%
  left_join(candidate_resolution, by = "Franklin_record_ID") %>%
  mutate(
    across(c(n_exact_allele_candidates, n_cdna_supported, n_protein_supported), ~replace_na(.x, 0L)),
    resolution_method = case_when(
      !has_complete_allele_key ~ "incomplete_allele_key",
      n_exact_allele_candidates == 0L & str_to_upper(`Variation Type`) == "INDEL" ~ "normalization_required",
      n_exact_allele_candidates == 0L ~ "no_exact_allele_candidate",
      n_exact_allele_candidates == 1L ~ "unique_exact_genomic_allele",
      n_cdna_supported == 1L ~ "resolved_by_exact_cdna_within_allele",
      n_cdna_supported != 1L & n_protein_supported == 1L ~ "resolved_by_exact_protein_within_allele",
      TRUE ~ "exact_allele_multiple_variant_ids"
    ),
    selected_Variant_ID = case_when(
      resolution_method == "unique_exact_genomic_allele" ~ sole_candidate,
      resolution_method == "resolved_by_exact_cdna_within_allele" ~ cdna_candidate,
      resolution_method == "resolved_by_exact_protein_within_allele" ~ protein_candidate,
      TRUE ~ NA_character_
    ),
    match_accepted_primary = resolution_method != "exact_allele_multiple_variant_ids" &
      n_exact_allele_candidates > 0L
  )

franklin_resolved_matches <- franklin_match_candidates %>%
  left_join(
    franklin_match_resolution %>%
      select(Franklin_record_ID, resolution_method, match_accepted_primary,
             selected_Variant_ID),
    by = "Franklin_record_ID", relationship = "many-to-one"
  ) %>%
  filter(
    resolution_method == "exact_allele_multiple_variant_ids" |
      Variant_ID == selected_Variant_ID
  ) %>%
  transmute(
    Franklin_record_ID, Variant_ID, Allele_ID,
    match_method = resolution_method, match_accepted_primary,
    Gene_concordant, Cdna_exact,
    Protein_support = Protein_exact
  )
stopifnot(!anyDuplicated(franklin_resolved_matches[c("Franklin_record_ID", "Variant_ID")]),
          !anyNA(franklin_resolved_matches$Allele_ID))

# 5. Attach both Franklin outputs at variant and allele levels -------------
variant_franklin_long <- franklin_resolved_matches %>%
  left_join(
    franklin_records %>% mutate(
      Franklin_simple = classify_genoox(`Genoox Classification`)) %>%
      select(Franklin_record_ID, `Aggregated Prediction`, `Genoox Classification`, Franklin_simple),
    by = "Franklin_record_ID", relationship = "many-to-one"
  )

variant_franklin_summary <- variant_franklin_long %>%
  group_by(Variant_ID) %>%
  summarise(
    n_franklin_records = n_distinct(Franklin_record_ID),
    Franklin_match_method = single_value(match_method), Gene_concordant = all(Gene_concordant),
    included_primary = any(match_accepted_primary),
    Aggregated_Prediction = single_value(`Aggregated Prediction`),
    Genoox_Classification = single_value(`Genoox Classification`),
    Franklin_simple = single_value(Franklin_simple), .groups = "drop"
  )

variant_master <- variant_master_base %>%
  left_join(variant_franklin_summary, by = "Variant_ID") %>%
  mutate(n_franklin_records = replace_na(n_franklin_records, 0L),
         has_franklin_match = replace_na(included_primary, FALSE))
stopifnot(nrow(variant_master) == n_distinct(variant_master$Variant_ID))

# Allele-level descriptive dataset: duplicated OLIDA record IDs never become independent observations.
# Excluded multiple-ID cases remain documented but do not contribute Franklin classifications.
allele_master <- variant_master %>%
  filter(!is.na(Allele_ID)) %>%
  group_by(Allele_ID) %>%
  summarise(
    n_olida_variant_ids = n_distinct(Variant_ID),
    has_franklin_match = any(has_franklin_match),
    Franklin_match_method = single_value(Franklin_match_method),
    Aggregated_Prediction = single_value(Aggregated_Prediction),
    Genoox_Classification = single_value(Genoox_Classification),
    Franklin_simple = single_value(Franklin_simple),
    .groups = "drop"
  ) %>%
  mutate(
    Aggregated_Prediction = if_else(has_franklin_match, Aggregated_Prediction, NA_character_),
    Genoox_Classification = if_else(has_franklin_match, Genoox_Classification, NA_character_),
    Franklin_simple = if_else(has_franklin_match, Franklin_simple, NA_character_)
  )
stopifnot(!anyDuplicated(allele_master$Allele_ID))

# 6. Construct variant- and allele-combination memberships -----------------
# Preserve the Variant_ID x OLIDA_ID mapping for traceability, then create the primary Allele_ID x OLIDA_ID analysis unit.
# Incomplete OLIDA allele keys remain in a separate QC dataset and do not enter allele-level denominators.
variant_membership_long <- combination_raw %>%
  rename(OLIDA_ID = `OLIDA ID`, Oligogenic_Effect = `Oligogenic Effect`, Disease = Diseases) %>%
  separate_rows(`Associated Variants`, sep = ";") %>%
  mutate(Variant_ID = str_extract(str_trim(`Associated Variants`), "^[0-9]+")) %>%
  select(Variant_ID, OLIDA_ID, Oligogenic_Effect, Disease) %>% distinct() %>%
  left_join(
    variant_master %>% select(Variant_ID, Allele_ID, has_franklin_match,
      Franklin_match_method, Gene_concordant, Aggregated_Prediction,
      Genoox_Classification, Franklin_simple),
    by = "Variant_ID", relationship = "many-to-one"
  ) %>%
  mutate(has_smallvariant_record = Variant_ID %in% variant_master$Variant_ID,
         has_franklin_match = replace_na(has_franklin_match, FALSE))
stopifnot(!anyNA(variant_membership_long$Variant_ID),
          !anyDuplicated(variant_membership_long[c("Variant_ID", "OLIDA_ID")]),
          !any(str_detect(variant_membership_long$OLIDA_ID, fixed("|"))),
          !any(str_detect(variant_membership_long$Oligogenic_Effect, fixed("|"))))

membership_missing_allele_id <- variant_membership_long %>% filter(is.na(Allele_ID))

membership_long <- variant_membership_long %>%
  filter(!is.na(Allele_ID)) %>%
  group_by(Allele_ID, OLIDA_ID) %>%
  summarise(
    Oligogenic_Effect = single_value(Oligogenic_Effect),
    Disease = single_value(Disease),
    n_olida_variant_ids = n_distinct(Variant_ID),
    has_franklin_match = any(has_franklin_match),
    Franklin_match_method = single_value(Franklin_match_method),
    Gene_concordant = if_else(any(has_franklin_match), all(Gene_concordant[has_franklin_match]), NA),
    Aggregated_Prediction = single_value(Aggregated_Prediction),
    Genoox_Classification = single_value(Genoox_Classification),
    Franklin_simple = single_value(Franklin_simple),
    .groups = "drop"
  ) %>%
  mutate(
    Aggregated_Prediction = if_else(has_franklin_match, Aggregated_Prediction, NA_character_),
    Genoox_Classification = if_else(has_franklin_match, Genoox_Classification, NA_character_),
    Franklin_simple = if_else(has_franklin_match, Franklin_simple, NA_character_)
  )
stopifnot(!anyDuplicated(membership_long[c("Allele_ID", "OLIDA_ID")]))

# 7. Construct one analytical row per OLIDA combination -------------------
combination_profile <- membership_long %>%
  group_by(OLIDA_ID) %>%
  summarise(
    n_alleles = n_distinct(Allele_ID),
    n_alleles_with_franklin = n_distinct(Allele_ID[has_franklin_match]),
    n_deleterious = n_distinct(Allele_ID[Aggregated_Prediction == "Deleterious"], na.rm = TRUE),
    n_uncertain = n_distinct(Allele_ID[Aggregated_Prediction == "Uncertain"], na.rm = TRUE),
    n_benign = n_distinct(Allele_ID[Aggregated_Prediction == "Benign"], na.rm = TRUE),
    n_PLP = n_distinct(Allele_ID[Franklin_simple == "P/LP"], na.rm = TRUE),
    n_VUS = n_distinct(Allele_ID[Franklin_simple == "VUS"], na.rm = TRUE),
    n_BLB = n_distinct(Allele_ID[Franklin_simple == "B/LB"], na.rm = TRUE),
    n_possibly_pathogenic = n_distinct(Allele_ID[Franklin_simple == "Possibly pathogenic"], na.rm = TRUE),
    n_possibly_benign = n_distinct(Allele_ID[Franklin_simple == "Possibly benign"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    match_rate = n_alleles_with_franklin / n_alleles,
    coverage = case_when(n_alleles_with_franklin == 0 ~ "No coverage",
      n_alleles_with_franklin == n_alleles ~ "Complete coverage", TRUE ~ "Partial coverage"),
    prop_deleterious = if_else(n_alleles_with_franklin > 0, n_deleterious/n_alleles_with_franklin, NA_real_),
    prop_uncertain = if_else(n_alleles_with_franklin > 0, n_uncertain/n_alleles_with_franklin, NA_real_),
    prop_benign = if_else(n_alleles_with_franklin > 0, n_benign/n_alleles_with_franklin, NA_real_),
    prop_PLP = if_else(n_alleles_with_franklin > 0, n_PLP/n_alleles_with_franklin, NA_real_),
    prop_VUS = if_else(n_alleles_with_franklin > 0, n_VUS/n_alleles_with_franklin, NA_real_),
    prop_BLB = if_else(n_alleles_with_franklin > 0, n_BLB/n_alleles_with_franklin, NA_real_),
    prop_possibly_pathogenic = if_else(n_alleles_with_franklin > 0,
      n_possibly_pathogenic/n_alleles_with_franklin, NA_real_),
    prop_possibly_benign = if_else(n_alleles_with_franklin > 0,
      n_possibly_benign/n_alleles_with_franklin, NA_real_)
  )

combination_analysis <- combination_raw %>%
  rename(OLIDA_ID = `OLIDA ID`, Oligogenic_Effect = `Oligogenic Effect`, Disease = Diseases) %>%
  select(-`Associated Variants`) %>%
  left_join(combination_profile, by = "OLIDA_ID", relationship = "one-to-one") %>%
  mutate(
    across(c(n_alleles, n_alleles_with_franklin, n_deleterious, n_uncertain,
      n_benign, n_PLP, n_VUS, n_BLB, n_possibly_pathogenic, n_possibly_benign), ~replace_na(.x, 0L)),
    coverage = if_else(n_alleles == 0L, "No analyzable allele ID", coverage),
    match_rate = if_else(n_alleles == 0L, NA_real_, match_rate)
  )
stopifnot(nrow(combination_analysis) == n_distinct(combination_analysis$OLIDA_ID),
          all(combination_analysis$n_alleles_with_franklin <= combination_analysis$n_alleles),
          all(with(combination_analysis,
            n_PLP + n_possibly_pathogenic + n_VUS + n_possibly_benign + n_BLB ==
              n_alleles_with_franklin)),
          any(combination_analysis$Oligogenic_Effect == "Unknown"))

# 8. Calculate the required matching quality-control summary ---------------
position_qc <- variant_master_base %>% filter(!is.na(Chr), !is.na(Position_hg38)) %>%
  group_by(Chr, Position_hg38) %>%
  summarise(n_variant_ids = n_distinct(Variant_ID), n_alleles = n_distinct(Allele_ID, na.rm = TRUE), .groups = "drop")
duplicate_allele_count <- variant_master_base %>% filter(!is.na(Allele_ID)) %>%
  count(Allele_ID) %>% summarise(value = sum(n > 1)) %>% pull(value)
matching_qc_summary <- tibble(
  metric = c("franklin_records_total", "franklin_records_complete_allele_key",
    "franklin_records_accepted_primary", "franklin_records_without_accepted_match", "franklin_match_rate",
    "unique_variant_ids_matched_primary",
    "ambiguous_genomic_positions",
    "variant_ids_at_ambiguous_positions", "exact_alleles_with_multiple_variant_ids",
    "primary_excluded_exact_alleles",
    "variant_combination_memberships", "allele_combination_memberships",
    "olida_memberships_missing_allele_id", "unique_olida_combinations_represented",
    "combinations_complete_coverage", "combinations_partial_coverage",
    "combinations_no_coverage", "combinations_no_analyzable_allele_id"),
  value = c(nrow(franklin_records), sum(franklin_records$has_complete_allele_key),
    sum(franklin_match_resolution$match_accepted_primary), sum(!franklin_match_resolution$match_accepted_primary),
    mean(franklin_match_resolution$match_accepted_primary),
    n_distinct(franklin_resolved_matches$Variant_ID[franklin_resolved_matches$match_accepted_primary]),
    sum(position_qc$n_alleles > 1), sum(position_qc$n_variant_ids[position_qc$n_alleles > 1]),
    duplicate_allele_count,
    sum(franklin_match_resolution$resolution_method == "exact_allele_multiple_variant_ids"),
    nrow(variant_membership_long), nrow(membership_long),
    nrow(membership_missing_allele_id),
    n_distinct(membership_long$OLIDA_ID[membership_long$has_franklin_match]),
    sum(combination_analysis$coverage == "Complete coverage"),
    sum(combination_analysis$coverage == "Partial coverage"),
    sum(combination_analysis$coverage == "No coverage"),
    sum(combination_analysis$coverage == "No analyzable allele ID"))
)

# 8a. Checking validate matching and combination --------------------------

# Every Franklin input record must have exactly one matching outcome.
stopifnot(
  !anyNA(franklin_records$Franklin_record_ID),
  !anyDuplicated(franklin_records$Franklin_record_ID),
  !anyNA(franklin_match_resolution$Franklin_record_ID),
  !anyDuplicated(franklin_match_resolution$Franklin_record_ID),
  setequal(
    franklin_records$Franklin_record_ID,
    franklin_match_resolution$Franklin_record_ID
  )
)

# Matching methods must agree with the accepted/excluded decision.
accepted_methods <- c(
  "unique_exact_genomic_allele",
  "resolved_by_exact_cdna_within_allele",
  "resolved_by_exact_protein_within_allele"
)
excluded_methods <- c(
  "incomplete_allele_key",
  "normalization_required",
  "no_exact_allele_candidate",
  "exact_allele_multiple_variant_ids"
)

stopifnot(
  all(franklin_match_resolution$resolution_method %in%
        c(accepted_methods, excluded_methods)),
  all(franklin_match_resolution$match_accepted_primary ==
        (franklin_match_resolution$resolution_method %in% accepted_methods))
)

# Every accepted record must retain exactly one selected candidate.
primary_matches <- franklin_resolved_matches %>%
  filter(match_accepted_primary)

primary_decisions <- franklin_match_resolution %>%
  filter(match_accepted_primary)

stopifnot(
  !anyNA(primary_matches$Variant_ID),
  !anyNA(primary_matches$Allele_ID),
  !anyDuplicated(primary_matches$Franklin_record_ID),
  setequal(
    primary_matches$Franklin_record_ID,
    primary_decisions$Franklin_record_ID
  ),
  nrow(anti_join(
    primary_matches, primary_decisions,
    by = c(
      "Franklin_record_ID",
      "Variant_ID" = "selected_Variant_ID",
      "Allele_ID"
    )
  )) == 0L,
  nrow(anti_join(
    primary_matches, franklin_match_candidates,
    by = c("Franklin_record_ID", "Variant_ID", "Allele_ID")
  )) == 0L
)

# Preserve every combination exactly once, with a valid coverage category.
coverage_levels <- c(
  "Complete coverage",
  "Partial coverage",
  "No coverage",
  "No analyzable allele ID"
)

stopifnot(
  !anyNA(combination_analysis$OLIDA_ID),
  !anyDuplicated(combination_analysis$OLIDA_ID),
  setequal(combination_analysis$OLIDA_ID, combination_raw$`OLIDA ID`),
  all(combination_analysis$coverage %in% coverage_levels)
)


# 9. Export canonical processed datasets and supplementary files -----------
outputs <- list(
  franklin_records = franklin_records, franklin_match_candidates = franklin_match_candidates,
  franklin_match_resolution = franklin_match_resolution,
  franklin_resolved_matches = franklin_resolved_matches,
  variant_franklin_long = variant_franklin_long, variant_master = variant_master,
  allele_master = allele_master,
  variant_membership_long = variant_membership_long,
  membership_missing_allele_id = membership_missing_allele_id,
  membership_long = membership_long,
  combination_analysis = combination_analysis,
  matching_qc_summary = matching_qc_summary
)
indels_for_normalization <- franklin_match_resolution %>%
  filter(resolution_method == "normalization_required") %>%
  select(Franklin_record_ID, match_status = resolution_method,
         n_exact_allele_candidates, Allele_ID) %>%
  left_join(franklin_records, by = c("Franklin_record_ID", "Allele_ID"),
            relationship = "one-to-one")
indel_qc <- tibble(
  Kennzahl = c("Franklin records total", "Franklin INDEL records",
    "Unmatched INDELs requiring normalization", "Matching criterion"),
  Wert = c(nrow(franklin_records),
    sum(str_to_upper(clean_missing(franklin_records$`Variation Type`)) == "INDEL"),
    nrow(indels_for_normalization), "Exact Chr + Start Position + REF + ALT match")
)
# Review export for the 31 exact alleles with several unresolved OLIDA IDs.
multiple_id_records <- franklin_match_resolution %>%
  filter(resolution_method == "exact_allele_multiple_variant_ids") %>%
  select(
    Franklin_record_ID, match_status = resolution_method,
    n_exact_allele_candidates, Allele_ID
  ) %>%
  left_join(
    franklin_records,
    by = c("Franklin_record_ID", "Allele_ID"),
    relationship = "one-to-one"
  )

# Export only after all matching checks have passed.
iwalk(outputs, ~write.csv(
  .x,
  file = file.path(processed_dir, str_c(.y, ".csv")),
  row.names = FALSE,
  na = "NA",
  fileEncoding = "UTF-8"
))
indel_xlsx_written <- tryCatch({
  write_xlsx(
    list(Indels_zu_normalisieren = indels_for_normalization,
         QC_Zusammenfassung = indel_qc),
    file.path(supplementary_matching_dir, "indels.xlsx")
  )
  TRUE
}, error = function(e) {
  warning("INDEL workbook could not be overwritten (possibly open in Excel): ",
          conditionMessage(e), call. = FALSE)
  FALSE
})

write_xlsx(
  list(
    Mehrfach_ID_Allele = multiple_id_records
  ),
  file.path(supplementary_matching_dir, "multiple_id_alleles.xlsx")
)

cat("Canonical matching key: Chr + Position_hg38 + REF + ALT\n")
cat("Franklin records:", nrow(franklin_records), "\n")
cat("Accepted primary:", sum(franklin_match_resolution$match_accepted_primary),
    " Without accepted match:", sum(!franklin_match_resolution$match_accepted_primary), "\n")
cat("Variant_IDs represented among resolved and multiple-ID candidate matches:",
    n_distinct(franklin_resolved_matches$Variant_ID), "\n")
cat("Allele x combination memberships:", nrow(membership_long),
    " Combinations:", nrow(combination_analysis), "\n")
print(count(franklin_match_resolution, resolution_method, sort = TRUE))

