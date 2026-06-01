## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(mmfpx)

data("simulated_longitudinal_data")

str(simulated_longitudinal_data)

## -----------------------------------------------------------------------------
fp_results <- mmfp(
  data = simulated_longitudinal_data,
  outcome_vars = "volumeA",
  age_var = "age",
  id_var = "subj_id",
  visit_var = "visit",
  powers = c("log", 1),
  random_intercept_only = TRUE
)

fp_results$results$volumeA$summary

fp_results$results$volumeA$selection_summary

## -----------------------------------------------------------------------------
fp_results_fp1_only <- mmfp(
  data = simulated_longitudinal_data,
  outcome_vars = "volumeA",
  age_var = "age",
  id_var = "subj_id",
  visit_var = "visit",
  powers = c("log", 1),
  fp_models = "fp1",
  random_intercept_only = TRUE
)

fp_results_fp1_only$results$volumeA$summary

## -----------------------------------------------------------------------------
fp_results_fp2_only <- mmfp(
  data = simulated_longitudinal_data,
  outcome_vars = "volumeA",
  age_var = "age",
  id_var = "subj_id",
  visit_var = "visit",
  powers = c("log", 1),
  fp_models = "fp2",
  random_intercept_only = TRUE
)

fp_results_fp2_only$results$volumeA$summary

## -----------------------------------------------------------------------------
fp_results_interact <- mmfp(
  data = simulated_longitudinal_data,
  outcome_vars = "volumeA",
  age_var = "age",
  id_var = "subj_id",
  visit_var = "visit",
  static_formula = "*diagnosis",  # expands to .fp1*diagnosis (+ .fp2*diagnosis)
  powers = c("log", 1),
  random_intercept_only = TRUE
)

fp_results_interact$results$volumeA$summary

## -----------------------------------------------------------------------------
baseline_data <- simulated_longitudinal_data[
  simulated_longitudinal_data$visit == 1,
  ,
  drop = FALSE
]

norm_results <- normbytcv(
  data = baseline_data,
  TBV_var = "volumeB",
  ROI_vars = "volumeA",
  subject_id_vars = "subj_id",
  powers_set = c("log", 0.5, 1),
  proportion = 0.75,
  num_cores = 1,
  plot_models = FALSE,
  save_results = FALSE
)

purrr::map_dfr(
  norm_results$results$volumeA$filtered_models,
  ~ tibble::tibble(
    model_type = .x$model_type,
    p1 = .x$powers$p1,
    p2 = .x$powers$p2,
    AIC = .x$model_aic
  )
)

head(norm_results$adjusted_volumes)
