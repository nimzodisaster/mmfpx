#' Evaluate Fractional Polynomial Mixed-Effects Models
#'
#' Fits first- and second-order fractional polynomial mixed-effects models for one or more
#' outcome variables. By default both FP1 and FP2 model classes are evaluated, but you can
#' restrict fitting to either FP1 or FP2 models. The function returns fit statistics that support
#' selecting the most appropriate model structure while offering optional scaling of the
#' age variable.
#'
#' @param data A data frame containing the variables required for the models.
#' @param outcome_vars A character vector with the names of the outcome variables to model.
#' @param age_var The name of the strictly positive predictor that will be transformed via
#'   fractional polynomials. Defaults to `"age"`.
#' @param static_formula Optional additional fixed-effect terms supplied as a one-sided
#'   formula (for example `~ sex + site`) or a character string (for example
#'   `"sex + site"`). These predictors are added to every model. When supplied as a
#'   character, a shorthand of the form `"*var"` is supported to interact the
#'   fractional polynomial terms with `var`. For example, `"*sex"` yields
#'   `~ .fp1*sex` for FP1 models and `~ .fp1*sex + .fp2*sex` for FP2 models.
#' @param id_var Subject identifier used for the random-effects structure. Defaults to
#'   `"subj_id"`. Not used with `engine = "lm"`; supplying it explicitly with
#'   `engine = "lm"` is an error.
#' @param visit_var Visit (or time) variable used when estimating within-subject
#'   correlation structures. Defaults to `"visit"`. Not used with `engine = "lm"`;
#'   supplying it explicitly with `engine = "lm"` is an error.
#' @param engine Model-fitting backend to use: `"nlme"` (default), `"lme4"`,
#'   or `"lm"`. `engine = "lm"` fits ordinary linear models via `stats::lm`
#'   for cross-sectional or non-clustered data with no random effects. It treats
#'   all rows as independent and is not appropriate for longitudinal or repeated-
#'   measures data where repeated observations per subject violate independence.
#' @param correlation Logical; if `TRUE`, a symmetric correlation structure is included
#'   using `nlme::corSymm` with `visit_var` as the ordering variable. This option
#'   is only supported with `engine = "nlme"` and errors for `engine = "lme4"` or
#'   `engine = "lm"`. Defaults to `FALSE`.
#' @param powers Candidate powers to consider for the fractional polynomial terms.
#'   Defaults to `c(-3, -2, -1, -0.5, "log", 0.5, 1, 2, 3)`.
#' @param fp_models Which fractional polynomial model classes to fit:
#'   `"both"` (default) fits both FP1 and FP2 models, `"fp1"` fits only
#'   one-term FP1 models, and `"fp2"` fits only two-term FP2 models.
#' @param random_intercept_only Logical; if `TRUE`, fits models with random intercepts
#'   only. If `FALSE`, random slopes for the fractional polynomial terms are also fitted.
#'   Not meaningful with `engine = "lm"` because ordinary linear models have no
#'   random effects; `random_intercept_only = TRUE` errors with `engine = "lm"`.
#' @param centering The scaling strategy for the age variable. Despite the argument
#'   name, these options scale but do not mean-center `age_var`. Use `"none"` to keep
#'   raw values (default), `"ratio"` to divide `age_var` by `centering_value`, or
#'   `"sd"` to divide by the sample standard deviation of `age_var`.
#' @param centering_value Optional positive value used for ratio scaling. When `NULL`,
#'   the mean of `age_var` is used. Only relevant when `centering = "ratio"` and
#'   ignored otherwise.
#' @param control Optional control object passed to the selected engine. For
#'   `engine = "nlme"`, supply an object from `nlme::lmeControl`; for
#'   `engine = "lme4"`, supply an object from `lme4::lmerControl`. When `NULL`,
#'   engine-specific defaults are used. `engine = "lm"` takes no control object,
#'   so non-NULL `control` is an error. User-supplied control objects are passed
#'   through unvalidated for mixed-model engines.
#' @param verbose Logical; if `TRUE`, progress messages are printed for models that
#'   fail to converge.
#'
#' @return A list containing:
#'   \describe{
#'     \item{engine}{The backend used to fit the models ("nlme", "lme4", or "lm").}
#'     \item{centering\_mode}{The scaling strategy that was applied ("none", "ratio", or "sd").}
#'     \item{centering\_value}{The divisor used for ratio or SD scaling (or `NA\_real\_` when raw ages are used).}
#'     \item{fp\_models}{The fitted fractional polynomial model classes ("both", "fp1", or "fp2").}
#'     \item{predictor\_variable}{The name of the transformed age column used in the models.}
#'     \item{results}{For each outcome variable, a list containing `summary` and
#'       `selection_summary`. `summary` is a tidy tibble of fit statistics with
#'       `delta_AIC`, `evidence_ratio`, Akaike `weight`, and convergence and
#'       singular-fit flags. Weights are relative AIC support within the fitted
#'       candidate set for that outcome; they are not probabilities that a model
#'       is correct and do not measure out-of-sample performance. For
#'       `engine = "lme4"`, `converged` is `FALSE` when the optimizer reports a
#'       nonzero return code or lme4 reports post-fit convergence messages, and
#'       `singular` is computed with `lme4::isSingular()`. For `engine = "nlme"`,
#'       returned fits are marked `converged = TRUE` because `nlme::lme` failures
#'       are handled as errors; `singular` is `NA` because nlme does not provide
#'       an equivalent singular-fit diagnostic. For `engine = "lm"`, `converged`
#'       is always `TRUE` because ordinary least squares has no iterative
#'       convergence check, and `singular` is `NA`; rank deficiency in `lm` fits is
#'       not separately checked. `selection_summary` reports the best model, its
#'       weight, the number of models with `delta_AIC <= 2`, and the summed weight
#'       of those near-top models. The `delta_AIC <= 2` cutoff is a convention for
#'       substantial support, not a rule.}
#'   }
#'
#' @examples
#' \dontrun{
#' mmfp\_results <- mmfp(
#'   data = demo\_data,
#'   outcome\_vars = c("roi\_volume"),
#'   age\_var = "age",
#'   id\_var = "subj\_id",
#'   visit\_var = "visit"
#' )
#'
#' # Fit only one-term fractional polynomial (FP1) models
#' mmfp\_results\_fp1 <- mmfp(
#'   data = demo\_data,
#'   outcome\_vars = c("roi\_volume"),
#'   age\_var = "age",
#'   fp\_models = "fp1",
#'   id\_var = "subj\_id",
#'   visit\_var = "visit"
#' )
#'
#' # Fit only two-term fractional polynomial (FP2) models
#' mmfp\_results\_fp2 <- mmfp(
#'   data = demo\_data,
#'   outcome\_vars = c("roi\_volume"),
#'   age\_var = "age",
#'   fp\_models = "fp2",
#'   id\_var = "subj\_id",
#'   visit\_var = "visit"
#' )
#'
#' # Interactions with FP terms using shorthand
#' # If your data contain a grouping variable like 'diagnosis', you can request
#' # interactions with the FP terms by supplying a character shorthand such as
#' # "*diagnosis". FP1 models will include .fp1*diagnosis, while FP2 models will
#' # include both .fp1*diagnosis and .fp2*diagnosis.
#' mmfp\_results\_int <- mmfp(
#'   data = demo\_data,
#'   outcome\_vars = c("roi\_volume"),
#'   age\_var = "age",
#'   static\_formula = "*diagnosis",
#'   id\_var = "subj\_id",
#'   visit\_var = "visit"
#' )
#'
#' # Alternatively, you can specify explicit interactions in a one-sided formula
#' # using the internal .fp1/.fp2 placeholders, for example:
#' #   static_formula = ~ .fp1*diagnosis + .fp2*diagnosis
#' # The shorthand is recommended for convenience.
#' }
#'
#' @import dplyr
#' @importFrom tidyr drop_na
#' @importFrom purrr map
#' @importFrom nlme lme lmeControl corSymm
#' @importFrom lme4 isSingular lmer lmerControl
#' @importFrom stats AIC BIC as.formula lm logLik terms sd
#' @importFrom utils combn
#' @importFrom rlang .data sym :=
#' @export
mmfp <- function (data, outcome_vars, age_var = "age", static_formula = NULL,
                  id_var = "subj_id", visit_var = "visit",
                  engine = c("nlme", "lme4", "lm"), correlation = FALSE,
                  powers = c(-3, -2, -1, -0.5, "log", 0.5, 1, 2, 3),
                  fp_models = c("both", "fp1", "fp2"), random_intercept_only = FALSE,
                  centering = c("none", "ratio", "sd"), centering_value = NULL,
                  control = NULL, verbose = FALSE)
{
  engine <- match.arg(engine)
  if (engine == "lm") {
    if (!missing(id_var)) {
      stop("'id_var' is not used when engine = 'lm'. engine = 'lm' fits ordinary linear models with no subject grouping or random effects; use engine = 'nlme' or engine = 'lme4' for subject-level random effects.")
    }
    if (!missing(visit_var)) {
      stop("'visit_var' is not used when engine = 'lm'. engine = 'lm' fits ordinary linear models with no subject grouping or random effects; use engine = 'nlme' or engine = 'lme4' for subject-level random effects.")
    }
    if (isTRUE(random_intercept_only)) {
      stop("'random_intercept_only' is not meaningful when engine = 'lm'. engine = 'lm' fits ordinary linear models with no subject grouping or random effects; use engine = 'nlme' or engine = 'lme4' for random effects.")
    }
    if (isTRUE(correlation)) {
      stop("'correlation' is not supported when engine = 'lm'. Correlation structures require engine = 'nlme'; engine = 'lm' fits ordinary linear models with no subject grouping or random effects.")
    }
    if (!is.null(control)) {
      stop("'control' is not supported when engine = 'lm'. engine = 'lm' fits ordinary linear models and stats::lm takes no control object; remove 'control' or use engine = 'nlme' or engine = 'lme4'.")
    }
  }
  centering <- match.arg(centering)
  fp_models <- match.arg(fp_models)
  if (engine != "nlme" && isTRUE(correlation)) {
    stop("correlation structures require engine = 'nlme'.")
  }
  if (missing(data) || is.null(data)) {
    stop("'data' must be supplied.")
  }
  if (missing(outcome_vars) || length(outcome_vars) == 0) {
    stop("Supply at least one outcome variable name in 'outcome_vars'.")
  }
  if (!all(outcome_vars %in% names(data))) {
    missing_outcomes <- setdiff(outcome_vars, names(data))
    stop("The following outcome variables are missing from 'data': ",
         paste(missing_outcomes, collapse = ", "))
  }
  if (!age_var %in% names(data)) {
    stop("The specified 'age_var' is not found in 'data'.")
  }
  if (engine != "lm" && !id_var %in% names(data)) {
    stop("The specified 'id_var' is not found in 'data'.")
  }
  if (engine != "lm" && !visit_var %in% names(data)) {
    stop("The specified 'visit_var' is not found in 'data'.")
  }
  age_values <- data[[age_var]]
  if (any(is.na(age_values))) {
    warning("'age_var' contains missing values. Rows with missing values will be dropped.")
  }
  if (any(age_values <= 0, na.rm = TRUE)) {
    stop("All values of 'age_var' must be strictly positive to support fractional-polynomial transforms.")
  }
  if (centering == "ratio") {
    if (is.null(centering_value)) {
      centering_value <- mean(age_values, na.rm = TRUE)
    }
    if (!is.numeric(centering_value) || length(centering_value) !=
        1 || is.na(centering_value)) {
      stop("'centering_value' must be a single numeric value.")
    }
    if (centering_value <= 0) {
      stop("'centering_value' must be strictly positive.")
    }
    centering_value_used <- centering_value
  }
  else if (centering == "sd") {
    if (!is.null(centering_value)) {
      warning("'centering_value' is ignored when centering = 'sd'.")
    }
    centering_value_used <- stats::sd(age_values, na.rm = TRUE)
    if (!is.finite(centering_value_used) || centering_value_used <=
        0) {
      stop("The standard deviation of 'age_var' must be positive to support SD scaling.")
    }
  }
  else {
    if (!is.null(centering_value)) {
      warning("'centering_value' is ignored when centering = 'none'.")
    }
    centering_value_used <- NA_real_
  }
  if (is.null(control)) {
    control <- if (engine == "nlme") {
      nlme::lmeControl(opt = "optim", maxIter = 2000,
                       msMaxIter = 2000, niterEM = 100, msVerbose = verbose)
    } else if (engine == "lme4") {
      lme4::lmerControl()
    } else {
      NULL
    }
  }
  parse_static_formula <- function(x) {
    if (is.null(x) || (is.character(x) && all(trimws(x) ==
                                              ""))) {
      return(list(expr = "", vars = character(0), star_terms = character(0),
                  plain_terms = character(0), has_star = FALSE))
    }
    if (inherits(x, "formula")) {
      rhs <- as.character(x)
      rhs <- rhs[length(rhs)]
      rhs <- trimws(rhs)
      if (rhs == "") {
        return(list(expr = "", vars = character(0), star_terms = character(0),
                    plain_terms = character(0), has_star = FALSE))
      }
      parsed_formula <- stats::as.formula(paste("~", rhs))
      vars <- all.vars(parsed_formula)
      plain_terms <- attr(stats::terms(parsed_formula), "term.labels")
      return(list(expr = rhs, vars = vars, star_terms = character(0),
                  plain_terms = plain_terms, has_star = FALSE))
    }
    if (is.character(x) && length(x) == 1) {
      rhs <- trimws(x)
      if (rhs == "") {
        return(list(expr = "", vars = character(0), star_terms = character(0),
                    plain_terms = character(0), has_star = FALSE))
      }
      pieces <- strsplit(rhs, "\\+")[[1]]
      pieces <- trimws(pieces)
      is_star <- grepl("^\\*", pieces)
      star_raw <- gsub("^\\*", "", pieces[is_star])
      star_raw <- trimws(star_raw)
      star_raw <- star_raw[nzchar(star_raw)]
      non_star <- pieces[!is_star]
      non_star <- non_star[nzchar(non_star)]
      plain_terms <- character(0)
      plain_vars <- character(0)
      if (length(non_star) > 0) {
        plain_formula <- stats::as.formula(paste("~", paste(non_star,
                                                            collapse = "+")))
        plain_terms <- attr(stats::terms(plain_formula), "term.labels")
        plain_vars <- all.vars(plain_formula)
      }
      star_vars <- if (length(star_raw) > 0) {
        all.vars(stats::as.formula(paste("~", paste(star_raw,
                                                    collapse = "+"))))
      } else {
        character(0)
      }
      vars <- unique(c(plain_vars, star_vars))
      return(list(expr = rhs, vars = vars, star_terms = star_raw,
                  plain_terms = plain_terms, has_star = any(is_star)))
    }
    stop("'static_formula' must be NULL, a one-sided formula, or a single character string.")
  }

  static_info <- parse_static_formula(static_formula)

  additional_vars <- setdiff(unique(static_info$vars), c(".fp1", ".fp2"))

  required_vars <- unique(c(outcome_vars, age_var, additional_vars))
  if (engine != "lm") {
    required_vars <- unique(c(required_vars, id_var, visit_var))
  }
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following variables required for modelling are missing from 'data': ",
         paste(missing_vars, collapse = ", "))
  }
  predictor_var <- switch(centering, ratio = paste0(age_var,
                                                    "_cr"), sd = paste0(age_var, "_sd"), age_var)
  predictor_sym <- rlang::sym(predictor_var)
  select_vars <- unique(c(required_vars, predictor_var))
  data_model <- data %>% dplyr::mutate(`:=`(!!predictor_sym,
                                            if (centering == "none") {
                                              .data[[age_var]]
                                            }
                                            else {
                                              .data[[age_var]]/centering_value_used
                                            })) %>% dplyr::select(dplyr::all_of(select_vars)) %>%
    tidyr::drop_na()
  if (nrow(data_model) == 0) {
    stop("No observations remain after removing rows with missing values in required variables.")
  }
  format_power <- function(p) {
    if (is.character(p)) {
      p_lower <- tolower(p)
      numeric_power <- suppressWarnings(as.numeric(p_lower))
      if (p_lower == "log" || isTRUE(all.equal(numeric_power, 0))) {
        return("log")
      }
      return(p_lower)
    }
    if (isTRUE(all.equal(p, 0))) {
      return("log")
    }
    as.character(p)
  }
  powers <- powers[!duplicated(vapply(powers, format_power, character(1)))]
  if (length(powers) == 0) {
    stop("Provide at least one candidate power in 'powers'.")
  }
  is_same_power <- function(p1, p2) {
    format_power(p1) == format_power(p2)
  }
  transform_power <- function(x, p) {
    if (is.character(p)) {
      p_lower <- tolower(p)
      numeric_power <- suppressWarnings(as.numeric(p_lower))
      if (p_lower == "log" || isTRUE(all.equal(numeric_power, 0))) {
        return(log(x))
      }
      if (is.na(numeric_power)) {
        stop("Unrecognised power: ", p)
      }
      return(x^numeric_power)
    }
    if (isTRUE(all.equal(p, 0))) {
      return(log(x))
    }
    x^as.numeric(p)
  }
  build_model_data <- function(df, p1, p2 = NULL) {
    predictor_values <- df[[predictor_var]]
    fp1 <- transform_power(predictor_values, p1)
    df$.fp1 <- fp1
    if (is.null(p2)) {
      df$.fp2 <- NULL
      fp_terms <- c(".fp1")
    }
    else {
      if (is_same_power(p1, p2)) {
        # Repeated-power FP2 second term is x^p * log(x), per Royston-Altman.
        df$.fp2 <- fp1 * log(predictor_values)
      }
      else {
        df$.fp2 <- transform_power(predictor_values,
                                   p2)
      }
      fp_terms <- c(".fp1", ".fp2")
    }
    list(data = df, fp_terms = fp_terms)
  }
  build_fixed_formula <- function(outcome, fp_terms) {
    rhs_terms <- fp_terms
    if (static_info$expr != "") {
      if (isTRUE(static_info$has_star) && length(static_info$star_terms) >
          0) {
        star_expanded <- unlist(lapply(static_info$star_terms,
                                       function(v) paste(fp_terms, v,
                                                         sep = ":")))
        rhs_terms <- c(rhs_terms, static_info$star_terms,
                       star_expanded, static_info$plain_terms)
      }
      else {
        rhs_terms <- c(rhs_terms, static_info$plain_terms)
      }
    }
    stats::reformulate(unique(rhs_terms), response = outcome)
  }
  build_random_formula <- function(fp_terms) {
    if (random_intercept_only) {
      return(stats::as.formula(paste("~ 1 |", id_var)))
    }
    rhs <- paste(fp_terms, collapse = " + ")
    stats::as.formula(paste("~", rhs, "|", id_var))
  }
  fit_single_model <- function(model_data, outcome, p1, p2 = NULL) {
    fp_info <- build_model_data(model_data, p1, p2)
    fixed_formula <- build_fixed_formula(outcome, fp_info$fp_terms)
    tryCatch({
      if (engine == "nlme") {
        random_formula <- build_random_formula(fp_info$fp_terms)
        if (correlation) {
          nlme::lme(fixed = fixed_formula, data = fp_info$data,
                    random = random_formula, correlation = nlme::corSymm(form = stats::as.formula(paste("~",
                                                                                                        visit_var))), control = control, method = "ML")
        }
        else {
          nlme::lme(fixed = fixed_formula, data = fp_info$data,
                    random = random_formula, control = control,
                    method = "ML")
        }
      }
      else if (engine == "lme4") {
        re <- if (random_intercept_only) {
          "1"
        } else {
          paste(fp_info$fp_terms, collapse = " + ")
        }
        random_term <- sprintf("(%s | %s)", re, id_var)
        combined_formula <- stats::update(fixed_formula,
                                          paste(". ~ . +", random_term))
        # lme4 convergence warnings are informative and do not mark a fit as failed.
        lme4::lmer(combined_formula, data = fp_info$data, REML = FALSE,
                   control = control)
      }
      else {
        stats::lm(fixed_formula, data = fp_info$data)
      }
    }, error = function(e) {
      if (verbose) {
        message("Model failed for ", outcome, ": p1=",
                format_power(p1), if (!is.null(p2))
                  paste(", p2=", format_power(p2), sep = ""),
                ". Error: ", conditionMessage(e))
      }
      NULL
    })
  }
  model_diagnostics <- function(model) {
    if (engine == "lme4") {
      convergence_messages <- model@optinfo$conv$lme4$messages
      optimizer_code <- model@optinfo$conv$opt
      converged <- is.null(convergence_messages) &&
        (is.null(optimizer_code) || identical(optimizer_code, 0L) ||
           identical(optimizer_code, 0))
      return(list(converged = converged,
                  singular = lme4::isSingular(model)))
    }
    if (engine == "lm") {
      # Rank deficiency in lm fits is not separately checked here.
      return(list(converged = TRUE, singular = NA))
    }
    list(converged = TRUE, singular = NA)
  }
  add_aic_weights <- function(outcome_summary) {
    outcome_summary$delta_AIC <- outcome_summary$AIC -
      min(outcome_summary$AIC, na.rm = TRUE)
    outcome_summary$evidence_ratio <- exp(-0.5 * outcome_summary$delta_AIC)
    outcome_summary$weight <- outcome_summary$evidence_ratio/
      sum(outcome_summary$evidence_ratio)
    outcome_summary
  }
  empty_selection_summary <- function() {
    dplyr::tibble(best_model_id = NA_character_,
                  best_weight = NA_real_,
                  n_models_within_2 = 0L,
                  weight_concentration = NA_real_)
  }
  build_selection_summary <- function(outcome_summary) {
    if (nrow(outcome_summary) == 0) {
      return(empty_selection_summary())
    }
    best_idx <- which.min(outcome_summary$AIC)
    within_2 <- outcome_summary$delta_AIC <= 2
    dplyr::tibble(best_model_id = outcome_summary$model_id[[best_idx]],
                  best_weight = outcome_summary$weight[[best_idx]],
                  n_models_within_2 = sum(within_2),
                  weight_concentration = sum(outcome_summary$weight[within_2]))
  }
  fp2_combinations <- if (length(powers) >= 2) {
    utils::combn(powers, 2, simplify = FALSE)
  } else {
    list()
  }
  fp2_combinations <- c(fp2_combinations, lapply(powers, function(p) c(p,
                                                                       p)))
  outcome_results <- purrr::map(outcome_vars, function(outcome) {
    summaries <- list()
    if (fp_models %in% c("both", "fp1")) {
      for (p in powers) {
        model <- fit_single_model(data_model, outcome, p)
        if (is.null(model)) {
          next
        }
        model_id <- paste0("FP1_", format_power(p))
        diagnostics <- model_diagnostics(model)
        model_summary <- dplyr::tibble(outcome = outcome,
                                       model_id = model_id, model_type = "FP1", p1 = format_power(p),
                                       p2 = NA_character_, logLik = as.numeric(stats::logLik(model)),
                                       AIC = stats::AIC(model), BIC = stats::BIC(model),
                                       converged = diagnostics$converged,
                                       singular = diagnostics$singular)
        summaries[[length(summaries) + 1]] <- model_summary
      }
    }
    if (fp_models %in% c("both", "fp2")) {
      for (pair in fp2_combinations) {
        p1 <- pair[[1]]
        p2 <- pair[[2]]
        model <- fit_single_model(data_model, outcome, p1,
                                  p2)
        if (is.null(model)) {
          next
        }
        model_id <- paste0("FP2_", format_power(p1), "_",
                           format_power(p2))
        diagnostics <- model_diagnostics(model)
        model_summary <- dplyr::tibble(outcome = outcome,
                                       model_id = model_id, model_type = "FP2", p1 = format_power(p1),
                                       p2 = format_power(p2), logLik = as.numeric(stats::logLik(model)),
                                       AIC = stats::AIC(model), BIC = stats::BIC(model),
                                       converged = diagnostics$converged,
                                       singular = diagnostics$singular)
        summaries[[length(summaries) + 1]] <- model_summary
      }
    }
    if (length(summaries) == 0) {
      warning("No models could be fit for outcome '", outcome,
              "'.")
      outcome_summary <- dplyr::tibble(outcome = character(0),
                                       model_id = character(0), model_type = character(0),
                                       p1 = character(0), p2 = character(0), logLik = numeric(0),
                                       AIC = numeric(0), BIC = numeric(0),
                                       delta_AIC = numeric(0), evidence_ratio = numeric(0),
                                       weight = numeric(0),
                                       converged = logical(0), singular = logical(0))
      selection_summary <- empty_selection_summary()
    }
    else {
      outcome_summary <- dplyr::bind_rows(summaries)
      outcome_summary <- add_aic_weights(outcome_summary)
      selection_summary <- build_selection_summary(outcome_summary)
      outcome_summary <- outcome_summary %>%
        dplyr::arrange(.data$model_type, .data$AIC)
    }
    list(summary = outcome_summary, selection_summary = selection_summary)
  })
  names(outcome_results) <- outcome_vars
  list(engine = engine,
       centering_mode = centering, centering_value = centering_value_used,
       fp_models = fp_models,
       predictor_variable = predictor_var, results = outcome_results)
}
