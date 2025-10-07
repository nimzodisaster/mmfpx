#' Multivariate BRMS Fractional-Polynomial Model Selection
#'
#' Fits joint multivariate Bayesian mixed-effects models over 2+ outcomes using
#' fractional-polynomial (FP) time terms, comparing candidate FP1 and/or FP2 bases.
#' Users supply templated fixed/random parts that include the token `{fp_terms}` which
#' this function expands to "fp" (FP1) or "fp + fp2" (FP2) per candidate.
#'
#' @param data Data frame with all variables.
#' @param outcomes Character vector of >=2 outcome column names.
#' @param age_var Strictly positive predictor to transform (default "age").
#' @param fixed_template One-sided RHS string/formula *containing* `{fp_terms}`.
#'   Example: "~ MAR*Sex + ICV + Motion + {fp_terms} + MAR:{fp_terms} + Sex:{fp_terms}".
#' @param rand_template Random-effects string *containing* `{fp_terms}`.
#'   Example: "(1 + {fp_terms} | ID)" or "(1 + {fp_terms} | gr(ID, by = MAR))".
#' @param powers Numeric/character vector of FP powers to try. Use 0 or "log" for log.
#'   Default: c(-2,-1,-0.5,0,0.5,1,2).
#' @param order Integer 1 or 2. If 2, tries all FP2 pairs incl. repeats (Royston–Sauerbrei).
#' @param rescor Logical; `set_rescor(TRUE)` for residual correlation across outcomes.
#' @param center One of c("none","ratio","sd"). If "ratio", divide by \code{center_value}
#'   (default = mean(age)); if "sd", divide by sd(age).
#' @param center_value Optional positive divisor when \code{center="ratio"}.
#' @param family brms family (default gaussian()).
#' @param prior brms prior specification (vector). Mild regularization is recommended.
#' @param chains,cores,iter,seed,control brms sampling controls.
#' @param criteria Character vector among c("loo","waic"). Default c("loo","waic").
#' @param select_by One of c("looic","waic"). Primary selector. Default "looic".
#' @param drop_bad Logical; if TRUE, models with max Rhat>1.01 or min ESS ratio<0.1
#'   are excluded from selection (kept in the table but not eligible to win).
#' @param keep_models Logical; if TRUE, return fitted model objects in result.
#'
#' @return A list with:
#'   \item{config}{Echo of key args and the FP grid considered.}
#'   \item{table}{Tibble of candidates with p1/p2, order, LOOIC/WAIC, diagnostics, winner flag.}
#'   \item{best}{The winning brmsfit (if any).}
#'   \item{fits}{Named list of brmsfit objects if keep_models=TRUE.}
#'
#' @details
#' - Your templates must include the token `{fp_terms}` in both fixed and random parts.
#' - For FP2 with repeated powers, the second term is \eqn{fp1 * log(x)} per standard FP2.
#' - Non-nested comparisons are evaluated by PSIS-LOO (primary) and WAIC.
#' - Convergence issues don't always ruin LOO/WAIC; `drop_bad=TRUE` helps guard selection.
#'
#' @examples
#' \dontrun{
#' best <- brms_mmfp(
#'   data=df,
#'   outcomes=c("vol_frontal","vol_striatum"),
#'   age_var="age",
#'   fixed_template="~ MAR*Sex + ICV + Motion + Scanner + {fp_terms} + MAR:{fp_terms} + Sex:{fp_terms}",
#'   rand_template="(1 + {fp_terms} | gr(ID, by = MAR))",
#'   powers=c(-2,-1,-0.5,0,0.5,1,2),
#'   order=1,
#'   prior=c(prior(normal(0,1), class="b"),
#'           prior(exponential(1), class="sd"),
#'           prior(exponential(1), class="sigma"),
#'           prior(lkj(2), class="cor")),
#'   chains=4, cores=4, iter=4000, seed=2025
#' )
#' best$table
#' summary(best$best)
#' }
#' @import brms dplyr purrr tibble
#' @importFrom stats as.formula sd terms model.matrix
#' @export
brms_mmfp <- function(
  data,
  outcomes,
  age_var = "age",
  fixed_template,
  rand_template,
  powers = c(-2,-1,-0.5,0,0.5,1,2),
  order = 1L,
  rescor = TRUE,
  center = c("none","ratio","sd"),
  center_value = NULL,
  family = gaussian(),
  prior = c(prior(normal(0,1), class="b"),
            prior(exponential(1), class="sd"),
            prior(exponential(1), class="sigma"),
            prior(lkj(2), class="cor")),
  chains = 4, cores = 4, iter = 4000, seed = 2025,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  criteria = c("loo","waic"),
  select_by = c("looic","waic"),
  drop_bad = TRUE,
  keep_models = FALSE
){
  stopifnot(length(outcomes) >= 2, is.character(outcomes))
  if (!all(outcomes %in% names(data))) {
    miss <- setdiff(outcomes, names(data))
    stop("Missing outcome columns: ", paste(miss, collapse=", "))
  }
  if (!age_var %in% names(data)) stop("age_var not found in data.")

  # age must be >0 for log/negative powers
  a <- data[[age_var]]
  if (any(is.na(a))) warning("age_var contains NA; rows with NA in needed columns will be dropped by brms.")
  if (any(a <= 0, na.rm=TRUE)) stop("All age values must be strictly positive for FP transforms.")

  center <- match.arg(center)
  select_by <- match.arg(select_by)
  criteria <- intersect(criteria, c("loo","waic"))

  # scale age if requested
  if (center == "ratio") {
    if (is.null(center_value)) center_value <- mean(a, na.rm=TRUE)
    if (!is.numeric(center_value) || center_value <= 0) stop("center_value must be positive numeric.")
    age_scaled <- a / center_value
    center_val_used <- center_value
  } else if (center == "sd") {
    s <- sd(a, na.rm=TRUE)
    if (!is.finite(s) || s <= 0) stop("sd(age) must be positive for center='sd'.")
    age_scaled <- a / s
    center_val_used <- s
  } else {
    age_scaled <- a
    center_val_used <- NA_real_
  }

  df0 <- data
  df0[[".age_"]] <- age_scaled

  # FP helpers
  to_num_power <- function(p) {
    if (is.character(p) && tolower(p) == "log") 0 else as.numeric(p)
  }
  fp1_transform <- function(x, p) if (isTRUE(all.equal(p,0))) log(x) else x^p
  # FP2: repeated powers use fp2 = fp1 * log(x)
  fp2_transform <- function(x, p1, p2, fp1) {
    if (isTRUE(all.equal(p1, p2))) fp1 * log(x) else fp1_transform(x, p2)
  }

  make_fp_grid <- function(powers, order=1L){
    P <- unique(powers)
    if (order == 1L) {
      tibble::tibble(order=1L, p1 = P, p2 = NA_real_)
    } else {
      # all ordered pairs incl. repeats
      pairs <- as.matrix(expand.grid(P, P, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
      tibble::tibble(order=2L, p1 = pairs[,1], p2 = pairs[,2])
    }
  }
  grid <- make_fp_grid(powers, order)

  # Build brms formulas: paste fixed+random for each outcome after substituting {fp_terms}
  build_bf_list <- function(fp_terms_str){
    # ensure {fp_terms} is present
    if (!grepl("\\{fp_terms\\}", fixed_template, fixed=FALSE))
      stop("fixed_template must contain the token {fp_terms}")
    if (!grepl("\\{fp_terms\\}", rand_template, fixed=FALSE))
      stop("rand_template must contain the token {fp_terms}")

    fixed_rhs <- sub("^~", "", if (inherits(fixed_template,"formula")) deparse(fixed_template) else fixed_template)
    fixed_rhs <- gsub("\\{fp_terms\\}", fp_terms_str, fixed_rhs, fixed=FALSE)
    rand_rhs  <- gsub("\\{fp_terms\\}", fp_terms_str, rand_template, fixed=FALSE)

    bfl <- lapply(outcomes, function(y) {
      as.formula(paste0(y, " ~ ", fixed_rhs, " + ", rand_rhs)) |> brms::bf()
    })
    names(bfl) <- outcomes
    bfl
  }

  # Fit one candidate
  fit_one <- function(df, ord, p1, p2) {
    # build fp columns
    p1n <- to_num_power(p1)
    df$fp  <- fp1_transform(df$.age_, p1n)
    fp_terms <- "fp"

    if (ord == 2L) {
      p2n <- to_num_power(p2)
      df$fp2 <- fp2_transform(df$.age_, p1n, p2n, df$fp)
      fp_terms <- "fp + fp2"
    } else {
      df$fp2 <- NULL
    }

    bfl <- build_bf_list(fp_terms)
    # combine with set_rescor
    mv_formula <- Reduce(`+`, bfl) + set_rescor(rescor)

    fit <- brm(
      mv_formula,
      data = df,
      family = family,
      prior = prior,
      chains = chains, cores = cores, iter = iter, seed = seed,
      backend = "cmdstanr",
      control = control
    )

    # add criteria
    if ("loo" %in% criteria) fit <- add_criterion(fit, "loo", moment_match = TRUE)
    if ("waic" %in% criteria) fit <- add_criterion(fit, "waic")

    fit
  }

  # Loop over grid
  fits <- vector("list", nrow(grid))
  tab  <- tibble::tibble(
    idx = seq_len(nrow(grid)),
    order = grid$order,
    p1 = grid$p1, p2 = grid$p2,
    looic = NA_real_, elpd_loo = NA_real_, loo_warn = NA_character_,
    waic = NA_real_, elpd_waic = NA_real_,
    max_rhat = NA_real_, min_neff_ratio = NA_real_,
    eligible = TRUE
  )

  for (i in seq_len(nrow(grid))) {
    ord <- grid$order[i]; p1 <- grid$p1[i]; p2 <- grid$p2[i]
    df <- df0
    fit <- try(fit_one(df, ord, p1, p2), silent = TRUE)
    if (inherits(fit, "try-error")) {
      tab$eligible[i] <- FALSE
      next
    }
    fits[[i]] <- fit

    # diagnostics
    rhat_max <- suppressWarnings(max(brms::rhat(fit), na.rm = TRUE))
    neff_min <- suppressWarnings(min(brms::neff_ratio(fit), na.rm = TRUE))
    tab$max_rhat[i]      <- rhat_max
    tab$min_neff_ratio[i] <- neff_min

    # criteria
    if (!is.null(fit$criteria$loo)) {
      tab$looic[i]   <- fit$criteria$loo$est["looic"]
      tab$elpd_loo[i] <- fit$criteria$loo$est["elpd_loo"]
      # summarize any Pareto-k problems
      bad_k <- fit$criteria$loo$diagnostics[["Pareto k"]]
      if (!is.null(bad_k) && any(bad_k > 0.7, na.rm=TRUE)) {
        tab$loo_warn[i] <- sprintf("k>0.7 in %d points", sum(bad_k > 0.7, na.rm=TRUE))
      } else tab$loo_warn[i] <- ""
    }
    if (!is.null(fit$criteria$waic)) {
      tab$waic[i]    <- fit$criteria$waic$est["waic"]
      tab$elpd_waic[i] <- fit$criteria$waic$est["elpd_waic"]
    }

    # drop_bad option
    if (drop_bad && (is.finite(rhat_max) && rhat_max > 1.01 || is.finite(neff_min) && neff_min < 0.1)) {
      tab$eligible[i] <- FALSE
    }
  }

  # Choose winner
  if (select_by == "looic" && any(is.finite(tab$looic))) {
    pool <- ifelse(tab$eligible, tab$looic, Inf)
    winner_idx <- which.min(pool)
  } else if (any(is.finite(tab$waic))) {
    pool <- ifelse(tab$eligible, tab$waic, Inf)
    winner_idx <- which.min(pool)
  } else {
    winner_idx <- NA_integer_
  }

  tab$winner <- FALSE
  best_fit <- NULL
  if (is.finite(winner_idx)) {
    tab$winner[winner_idx] <- TRUE
    best_fit <- fits[[winner_idx]]
  }

  out <- list(
    config = list(
      outcomes = outcomes,
      age_var = age_var,
      center = center,
      center_value = center_val_used,
      rescor = rescor,
      powers = powers,
      order = order,
      select_by = select_by,
      drop_bad = drop_bad
    ),
    table = tab |> dplyr::arrange(dplyr::desc(winner), !!dplyr::sym(select_by)),
    best  = best_fit,
    fits  = if (keep_models) stats::setNames(fits, paste0("ord",grid$order,"_p1",grid$p1,"_p2",grid$p2)) else NULL
  )
  class(out) <- c("brms_mmfp_result","list")
  out
}
