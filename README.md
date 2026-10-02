# OLIDA-Franklin Thesis Analysis

This repository contains the reproducible R workflow and the supporting data for the Bachelor thesis "From Single Variants to Oligogenic Combi-nations: A Statistical Evaluation of Frank-lin/Genoox Classifications Against OLIDA Categories for Rare Genetic Diseases". Results describe associations and do not establish that either framework is diagnostically correct or incorrect.

## Reproduce the analysis

Run from the repository root:

```r
source("run_all.R")
```

Required packages: `tidyverse`, `readxl`, `writexl`, `vegan`, `permute`, `glmmTMB`, `emmeans`, and `scales`.

A completed run records the R and package versions in `output/sessionInfo.txt`. Bootstrap refits and models used for AIC comparisons are checked before their results are exported. The diagnostic tables include the Pearson statistic, dispersion ratio and asymptotic Pearson p-value.

## Workflow

1. `01_prepare_data.R`: harmonization, allele-specific matching, and processed datasets
2. `02_descriptive_statistics.R`: primary sample composition
3. `03_genoox_permutation.R`: five-category Genoox PERMANOVA and PERMDISP
4. `04_genoox_regression.R`: grouped binomial P/LP regression and model checks
5. `05_aggregated_prediction_permutation.R`: three-category AP PERMANOVA and PERMDISP
6. `06_aggregated_prediction_regression.R`: grouped binomial Deleterious regression and model checks
7. `07_sensitivity_analyses.R`: two-allele and recurrent-allele-free restrictions for both outputs, plus exploratory Genoox regression on the complete AP sample
8. `08_figures.R`: the five empirical figures reported in the thesis

Shared statistical functions are stored in `standard_statistical_tests.R`.

The global profile test tables include `R_squared`: the fraction of total Euclidean profile variation associated with OLIDA effect category (multiply by 100 for a percentage). It is not a measure of prediction accuracy or a causal effect.

The common-sample sensitivity analysis retains the Genoox P/LP outcome and selects the same OLIDA combinations as the primary AP analysis. All three restrictions are evaluated independently, retaining the primary model specifications. Every sensitivity regression checks convergence and positive-definite Hessians for the full and intercept-only models. Bootstrap, AIC and boundary diagnostics are performed for the primary analyses in scripts 04 and 06, not repeated in script 07. These numerical estimation checks do not establish goodness of fit.

In `R/07_sensitivity_analyses.R`, section 2 constructs samples, section 3 runs profile tests for the two structural restrictions, and section 4 fits all sensitivity regressions through `summarize_binomial_fit()` in `R/standard_statistical_tests.R`. Sections 5–7 assemble the complete console overview, export the detailed tables and display the overview. The overview includes all five output/sample combinations and both M+M contrasts. Missing profile p-values in the AP-complete Genoox row indicate tests not performed. This common-sample analysis examines sample composition and is not a formal test of differences between the two outputs.

Sensitivity outputs in `output/tables/`:
- `07_sensitivity_profile_tests.csv`: PERMANOVA and PERMDISP for the two structural restrictions and both outputs.
- `07_sensitivity_regression_tests.csv`: global likelihood-ratio tests for all five output/sample combinations.
- `07_sensitivity_pairwise_holm.csv`: all three pairwise contrasts, confidence intervals and Holm-adjusted p-values for each regression.

## Repository contents

- `data/raw/`: the three source files read by the matching pipeline
- `data/processed/`: the canonical combination and allele-membership datasets used by the analyses
- `supplementary/matching/`: the exclusion lists reproduced in Appendix Tables A.1 and A.2
- `supplementary/franklin_upload/Formatted_FranklinUploadTemplate.xlsx`: the original formatted workbook submitted to Franklin, preserved unchanged as a reference for the annotation input. The corresponding Franklin export is `data/raw/Small Variants Franklin complete.csv`; the upload workbook is not required to run the R analysis.
- `output/tables/`: selected tables supporting reported descriptive, diagnostic, inferential, and sensitivity results
- `output/figures/`: the five empirical figures used in the thesis

One observation in the main analyses is one `OLIDA_ID`. Exact genomic alleles are counted once within each combination. Unknown combinations remain in the processed data but are excluded from inference comparing defined mechanisms. Genoox Classification and Aggregated Prediction are analysed separately.

The committed tables in output/tables serve as reference results. The pipeline checks sample sizes and data/model validity rather than requiring exact agreement with historical test statistics or p-values. Random seeds remain fixed for reproducibility; R or package version changes can still affect numerical results.
