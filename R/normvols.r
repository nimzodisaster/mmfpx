#' Normalize ROI Volumes by Total Brain Volume
#'
#' This function normalizes ROI volumes by total brain volume using bagged fractional polynomial models.
#' Adjusts volumetric regions of interest for their association with overall brain volume using fractional polynomials.
#'
#' @param data A data frame containing the ROI and TBV data.
#' @param TBV_var The name of the total brain volume variable as a string.
#' @param ROI_vars A vector of ROI variable names to normalize.
#' @param subject_id_vars A vector of identification variable names (e.g., subject ID, session ID).
#' @param proportion Define threshold for model bagging. The proportion (probability 0-1) that a model fits as well as the best model (default is 0.66).
#' @param powers_set A vector of powers for fractional polynomial transformations.
#'        Default is c(-2, -1.5, -1, -0.5, "log", 0.5, 1, 1.5, 2, 3).
#' @param num_cores Number of cores to use for parallel processing (default 1).
#' @param plot_models Logical, if TRUE, plots the selected model functions for each ROI variable.
#' @param save_results Logical, if TRUE, saves the results object to an RDS file.
#' @param results_name A string specifying the name of the results object to save (without file extension).
#'
#' @return A list containing:
#' \describe{
#'   \item{adjusted\_volumes}{A data frame with the adjusted ROI volumes.}
#'   \item{results}{A list containing the models and adjusted ROI values for each ROI variable.}
#' }
#'
#' @examples
#' \dontrun{
#' norm\_results <- normbytcv(data, "TBV", c("ROI1", "ROI2"), c("ID"), 0.66, plot_models = TRUE)
#' }
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom MASS rlm
#' @importFrom dplyr all_of select mutate group_by slice_head slice_sample ungroup
#' @importFrom tidyr drop_na
#' @importFrom purrr map
#' @importFrom stats AIC predict as.formula coef
#' @importFrom utils combn
#' @importFrom graphics lines legend
#' @importFrom grDevices rainbow
#' @export
normbytcv <- function(data, TBV_var, ROI_vars, subject_id_vars,
                      proportion = 0.66,
                      powers_set = c(-2, -1.5, -1, -0.5, "log", 0.5, 1, 1.5, 2, 3),
                      num_cores = 1,
                      plot_models = FALSE,
                      save_results = FALSE,
                      results_name = "normbytcv_results") {
  # Check that the data contains the necessary variables
  required_vars <- c(subject_id_vars, TBV_var, ROI_vars)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop("The following required variables are missing from the data: ", paste(missing_vars, collapse = ", "))
  }

  # Prepare the data (select only the necessary columns)
  data_filtered <- data %>% dplyr::select(dplyr::all_of(required_vars))

  # Remove rows with NA in any of the required variables
  data_filtered <- data_filtered %>% tidyr::drop_na()

  # Define required functions
  calculate_AIC_difference <- function(x) {
    -2 * log(x)
  }

  transform_TBV <- function(TBV, p) {
    if (p == "log") {
      log(TBV)
    } else {
      TBV^as.numeric(p)
    }
  }

  generate_model_data <- function(data, ROI_var, TBV_var, powers) {
    data <- data %>%
      dplyr::mutate(TBV = .data[[TBV_var]], ROI = .data[[ROI_var]])
    if (length(powers) == 1) {
      p <- powers$p1
      data <- data %>% dplyr::mutate(TBV_p = transform_TBV(.data$TBV, p))
      formula <- stats::as.formula("ROI ~ TBV_p")
    } else if (length(powers) == 2) {
      p1 <- powers$p1
      p2 <- powers$p2
      if (p1 == p2) {
        data <- data %>%
          dplyr::mutate(
            TBV_p = transform_TBV(.data$TBV, p1),
            TBV_p_log = .data$TBV_p * log(.data$TBV)
          )
        formula <- stats::as.formula("ROI ~ TBV_p + TBV_p_log")
      } else {
        data <- data %>%
          dplyr::mutate(
            TBV_p1 = transform_TBV(.data$TBV, p1),
            TBV_p2 = transform_TBV(.data$TBV, p2)
          )
        formula <- stats::as.formula("ROI ~ TBV_p1 + TBV_p2")
      }
    } else {
      stop("Invalid number of powers specified in 'powers' list.")
    }
    list(data = data, formula = formula)
  }

  adjust_ROI_FP <- function(data, ROI_var, TBV_var, powers_set) {
    model_list <- list()
    # FP1 models
    for (p in powers_set) {
      model_info <- generate_model_data(data, ROI_var, TBV_var, list(p1 = p))
      model <- tryCatch({
        MASS::rlm(model_info$formula, data = model_info$data)
      }, error = function(e) NULL)
      if (!is.null(model)) {
        model_aic <- stats::AIC(model)
        predicted_ROI <- stats::predict(model, newdata = model_info$data)
        mean_TBV <- mean(data[[TBV_var]], na.rm = TRUE)
        TBV_p_mean <- transform_TBV(mean_TBV, p)
        mean_predicted_ROI <- stats::coef(model)[1] + stats::coef(model)[2] * TBV_p_mean
        adjusted_ROI <- data[[ROI_var]] - (predicted_ROI - mean_predicted_ROI)
        model_list[[length(model_list) + 1]] <- list(
          model_type = "FP1",
          powers = list(p1 = p, p2 = NA),
          ROI_var = ROI_var,
          model_aic = model_aic,
          adjusted_ROI = adjusted_ROI,
          model = model
        )
      }
    }
    # FP2 models
    combinations <- utils::combn(powers_set, 2, simplify = FALSE)
    combinations <- c(combinations, lapply(powers_set, function(p) c(p, p)))
    for (pair in combinations) {
      p1 <- pair[1]
      p2 <- pair[2]
      model_info <- generate_model_data(data, ROI_var, TBV_var, list(p1 = p1, p2 = p2))
      model <- tryCatch({
        MASS::rlm(model_info$formula, data = model_info$data)
      }, error = function(e) NULL)
      if (!is.null(model)) {
        predicted_ROI <- stats::predict(model, newdata = model_info$data)
        mean_TBV <- mean(data[[TBV_var]], na.rm = TRUE)
        if (p1 == p2) {
          TBV_p_mean <- transform_TBV(mean_TBV, p1)
          TBV_p_log_mean <- TBV_p_mean * log(mean_TBV)
          mean_predicted_ROI <- coef(model)[1] + coef(model)[2] * TBV_p_mean + coef(model)[3] * TBV_p_log_mean
        } else {
          TBV_p1_mean <- transform_TBV(mean_TBV, p1)
          TBV_p2_mean <- transform_TBV(mean_TBV, p2)
          mean_predicted_ROI <- coef(model)[1] + coef(model)[2] * TBV_p1_mean + coef(model)[3] * TBV_p2_mean
        }
        adjusted_ROI <- data[[ROI_var]] - (predicted_ROI - mean_predicted_ROI)
        model_list[[length(model_list) + 1]] <- list(
          model_type = "FP2",
          powers = list(p1 = p1, p2 = p2),
          ROI_var = ROI_var,
          model_aic = stats::AIC(model),
          adjusted_ROI = adjusted_ROI,
          model = model
        )
      }
    }
    model_list
  }

  fit_roi_models <- function(ROI_var) {
    vars_to_keep <- c(subject_id_vars, TBV_var, ROI_var)
    data_for_model <- data_filtered %>% dplyr::select(dplyr::all_of(vars_to_keep))
    data_for_model <- data_for_model %>% tidyr::drop_na()
    result <- adjust_ROI_FP(data_for_model, ROI_var, TBV_var, powers_set)
    result
  }

  if (num_cores > 1) {
    cl <- parallel::makeCluster(num_cores)

    parallel::clusterExport(cl, varlist = c("data_filtered", "TBV_var", "powers_set", "subject_id_vars",
                                            "transform_TBV", "generate_model_data", "adjust_ROI_FP",
                                            "fit_roi_models"), envir = environment())

    parallel::clusterEvalQ(cl, {
      library(dplyr)
      library(MASS)
      library(purrr)
      library(tidyr)
    })

    results <- parallel::parLapply(cl, ROI_vars, fit_roi_models)

    parallel::stopCluster(cl)
  } else {
    results <- lapply(ROI_vars, fit_roi_models)
  }

  # Calculate the AIC criterion
  AICcriterion <- calculate_AIC_difference(proportion)

  # Initialize the adjusted_volumes data.frame with the subject_id_vars
  adjusted_volumes <- data_filtered %>% dplyr::select(dplyr::all_of(subject_id_vars))

  # Initialize the results list
  results_list <- list()

  # For each result (per ROI variable), select the best models and calculate the adjusted ROI values
  for (i in seq_along(results)) {
    result <- results[[i]]
    ROI_var <- ROI_vars[i]

    # Select models within AICcriterion of min_aic
    min_aic <- min(sapply(result, function(x) x$model_aic), na.rm = TRUE)
    filtered_models <- Filter(function(x) (x$model_aic <= (min_aic + AICcriterion)), result)

    # Extract adjusted_ROI values and calculate element-wise mean across lists
    adjusted_ROI_values <- purrr::map(filtered_models, ~ .x$adjusted_ROI) %>%
      simplify2array() %>%
      rowMeans(na.rm = TRUE)

    # Save the element-wise mean array as a new column in adjusted_volumes
    adjusted_ROI_var <- paste0(ROI_var, "_adj")
    adjusted_volumes[[adjusted_ROI_var]] <- adjusted_ROI_values

    # Save the filtered models and adjusted values to results_list
    results_list[[ROI_var]] <- list(
      filtered_models = filtered_models,
      adjusted_ROI_values = adjusted_ROI_values
    )
  }

  # Plot the selected models if requested
  if (plot_models) {
    for (ROI_var in ROI_vars) {
      data_for_model <- data_filtered %>% dplyr::select(dplyr::all_of(subject_id_vars), TBV_var, ROI_var)
      data_for_model <- data_for_model %>% tidyr::drop_na()

      # Get the filtered models
      filtered_models <- results_list[[ROI_var]]$filtered_models

      # Plot observed ROI vs. TBV
      plot(data_for_model[[TBV_var]], data_for_model[[ROI_var]], main = paste("ROI:", ROI_var),
           xlab = TBV_var, ylab = ROI_var, pch = 16, col = "blue")

      # Create a sequence of TBV values for plotting the model functions
      TBV_seq <- seq(min(data_for_model[[TBV_var]], na.rm = TRUE), max(data_for_model[[TBV_var]], na.rm = TRUE), length.out = 100)

      # For each model, compute the predicted ROI values over TBV_seq
      colors <- rainbow(length(filtered_models))
      for (j in seq_along(filtered_models)) {
        model_info <- filtered_models[[j]]
        model <- model_info$model
        powers <- model_info$powers

        # Create a data frame with TBV_seq transformed according to the model
        if (model_info$model_type == "FP1") {
          p <- powers$p1
          TBV_p <- transform_TBV(TBV_seq, p)
          newdata <- data.frame(TBV_p = TBV_p)
        } else if (model_info$model_type == "FP2") {
          p1 <- powers$p1
          p2 <- powers$p2
          if (p1 == p2) {
            TBV_p <- transform_TBV(TBV_seq, p1)
            TBV_p_log <- TBV_p * log(TBV_seq)
            newdata <- data.frame(TBV_p = TBV_p, TBV_p_log = TBV_p_log)
          } else {
            TBV_p1 <- transform_TBV(TBV_seq, p1)
            TBV_p2 <- transform_TBV(TBV_seq, p2)
            newdata <- data.frame(TBV_p1 = TBV_p1, TBV_p2 = TBV_p2)
          }
        }

        # Predict ROI values
        predicted_ROI <- predict(model, newdata = newdata)

        # Add the predicted line to the plot
        lines(TBV_seq, predicted_ROI, col = colors[j], lwd = 2)
      }

      legend("topright", legend = paste0("Model ", seq_along(filtered_models)), col = colors, lwd = 2)
    }
  }

  # Save the results object to a file if requested
  if (save_results) {
    saveRDS(results_list, file = paste0(results_name, ".rds"))
  }

  # Return the adjusted_volumes and results_list
  return(list(adjusted_volumes = adjusted_volumes, results = results_list))
}
