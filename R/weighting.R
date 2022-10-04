#' Estimate Weights for Generalizing Average Treatment Effect
#'
#' This function is designed for use within 'covariate_table()' and 'assess()'.
#'
#' @param data data frame comprised of "stacked" sample and target population data
#' @param sample_indicator variable name denoting binary sample membership (1 = in sample, 0 = out of sample)
#' @param treatment_indicator variable name denoting binary treatment assignment (ok if only available in sample, not population)
#' @param outcome variable name denoting outcome
#' @param covariates vector of covariate names in data set that predict sample membership
#' @param estimation_method method to estimate the probability of sample membership. Default is logistic regression ("lr").Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param disjoint_data logical. If TRUE, then sample and population data are considered disjoint. This affects calculation of the weights - see details for more information.
#' @export
#' @importFrom stats quantile
#' @importFrom crayon bold blue

weighting <- function(data,
                      sample_indicator,
                      treatment_indicator = NULL,
                      outcome = NULL,
                      covariates,
                      estimation_method = "lr",
                      disjoint_data = TRUE) {

  estimation_method <- tolower(estimation_method)

  data_names <- names(data)

  # Check whether data object is of type 'data.frame'
  assertthat::on_failure(is.data.frame) <- function(call, env) {

    "You must pass an object of type 'data.frame' to the 'data' argument."
  }

  assertthat::assert_that(is.data.frame(data))

  # Check whether sample indicator variable is binary
  is_sample_indicator_binary <- function(sample_indicator) {

    all(dplyr::pull(data, tidyselect::all_of(sample_indicator)) %in% c(0, 1))
  }

  assertthat::on_failure(is_sample_indicator_binary) <- function(call, env) {

    "Sample indicator variable must be coded as 0 (not in sample) or 1 (in sample)."
  }

  assertthat::assert_that(is_sample_indicator_binary(sample_indicator))

  # Check whether treatment indicator variable is a binary column in dataframe
  is_treatment_indicator_valid <- function(treatment_indicator) {

    if (is.null(treatment_indicator)) {return(TRUE)}

    if (!treatment_indicator %in% data_names) {return(FALSE)}

    valid_in_sample <- data %>%
      dplyr::filter(!!rlang::sym(sample_indicator) == 1) %>%
      dplyr::pull(!!rlang::sym(treatment_indicator)) %>%
      {. %in% c(0, 1)} %>%
      all()

    valid_out_of_sample <- data %>%
      dplyr::filter(!!rlang::sym(sample_indicator) == 0) %>%
      dplyr::pull(!!rlang::sym(treatment_indicator)) %>%
      is.na() %>%
      all()

    return(valid_in_sample && valid_out_of_sample)
  }

  assertthat::on_failure(is_treatment_indicator_valid) <- function(call, env) {

    if (!treatment_indicator %in% data_names) {

      return("If you wish to specify a treatment indicator variable, it must be one of the columns in the dataframe you have provided.")
    }

    valid_in_sample <- data %>%
      dplyr::filter(!!rlang::sym(sample_indicator) == 1) %>%
      dplyr::pull(!!rlang::sym(treatment_indicator)) %>%
      {. %in% c(0, 1)} %>%
      all()

    if (!valid_in_sample) {

      return("Treatment indicator variable must be coded as 0 (did not receive treatment) or 1 (received treatment) for all observations in the sample.")
    }

    valid_out_of_sample <- data %>%
      dplyr::filter(!!rlang::sym(sample_indicator) == 0) %>%
      dplyr::pull(!!rlang::sym(treatment_indicator)) %>%
      is.na() %>%
      all()

    if (!valid_out_of_sample) {

      return("Treatment indicator variable must be coded as NA for all observations not in the sample.")
    }
  }

  assertthat::assert_that(is_treatment_indicator_valid(treatment_indicator))

  # Check whether outcome variable is one of columns in dataframe
  is_outcome_valid <- function(outcome) {

    if (!is.null(outcome)) {

      outcome %in% data_names
    } else{

      return(TRUE)
    }
  }

  assertthat::on_failure(is_outcome_valid) <- function(call, env) {

    "If you wish to specify an outcome variable, it must be one of the columns in the dataframe you have provided."
  }

  assertthat::assert_that(is_outcome_valid(outcome))

  # Check whether selection covariates are variables in the dataframe provided
  invalid_covariates <- covariates %>% setdiff(data_names)

  assertthat::on_failure(is_empty) <- function(call, env) {

    paste("The following covariates are not variables in the data provided:\n",
          paste(crayon::blue$bold(invalid_covariates),
                collapse = ", ")
    )
  }

  assertthat::assert_that(is_empty(invalid_covariates))

  # Check whether estimation method is valid
  is_estimation_method_valid <- function(estimation_method){

    estimation_method %in% c("lr","rf","lasso")
  }

  assertthat::on_failure(is_estimation_method_valid) <- function(call, env) {

    "Invalid estimation method. Please choose one of 'lr' (logistic regression), 'rf' (random forest), or 'lasso' as your estimation method."
  }

  # Omit rows in which the sample indicator variable or the covariates contain missing values from the data
  data <- data %>%
    tidyr::drop_na(c(tidyselect::all_of(sample_indicator), tidyselect::all_of(covariates)))

  # Generate propensity scores
  data$ps <- .generate.ps(data, sample_indicator, covariates, estimation_method)

  participation_probs <- list(population = data$ps[which(data[,sample_indicator] == 0)],
                              sample = data$ps[which(data[,sample_indicator] == 1)])

  # Generate Weights #
  if(disjoint_data) {

    data <- data %>%
      dplyr::mutate(weights = ifelse(!!sym(sample_indicator) == 0,
                                     0,
                                     (1-ps)/ps))
  }

  else {

    data <- data %>%
      dplyr::mutate(weights = ifelse(!!sym(sample_indicator) == 0,
                                     0,
                                     1/ps))
  }

  # Trim any of the weights if necessary
  data$weights[which(data$weights == 0 & data[,sample_indicator] == 1)] <- quantile(data$weights[which(data[,sample_indicator] == 1)], 0.01, na.rm = TRUE)

  # Make weighted covariate table
  covariate_table_output <- .make.covariate.table(data = data,
                                                  sample_indicator = sample_indicator,
                                                  covariates = covariates,
                                                  weighted_table = TRUE,
                                                  estimation_method = estimation_method,
                                                  disjoint_data = disjoint_data)

  # Make histogram of weights
  weights_hist <- data %>%
    dplyr::filter(!!sym(sample_indicator) == 1) %>%
    ggplot(aes(x = weights)) +
    geom_histogram(bins = 20,
                   fill = viridis::viridis(1, alpha = 0.7),
                   color = "black") +
    theme_minimal() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Weights",
         y = "Frequency",
         title = "Histogram of Sample Weights")

  # Items to return out
  out <- list(participation_probs = participation_probs,
              weights = data$weights,
              covariate_table = covariate_table_output$covariate_table,
              covariate_kable = covariate_table_output$covariate_kable,
              hist = weights_hist)

  if(!is.null(outcome) & !is.null(treatment_indicator)) {

    # ESTIMATE TOTAL AVERAGE TREATMENT EFFECT

    # ADJUSTED TATE

    # Make weighted regression model predicting outcome with treatment
    TATE_model <- paste(outcome, treatment_indicator, sep = "~") %>%
      as.formula() %>%
      lm(data = data, weights = weights)

    # Extract total average treatment effect and standard error from weighted model
    TATE <- TATE_model %>%
      broom::tidy() %>%
      dplyr::filter(term == treatment_indicator) %>%
      dplyr::pull(estimate)

    TATE_se <- TATE_model %>%
      broom::tidy() %>%
      dplyr::filter(term == treatment_indicator) %>%
      dplyr::pull(std.error)

    # Calculate 95% confidence interval for total average treatment effect
    TATE_CI <- TATE + 2.262*TATE_se*c(-1, 1)

    TATE <- list(estimate = TATE,
                 SE = TATE_se,
                 CI = TATE_CI)

    # UNADJUSTED TATE

    # Make unweighted regression model predicting outcome with treatment
    TATE_model_unadj <- paste(outcome, treatment_indicator, sep = "~") %>%
      as.formula() %>%
      lm(data = data)

    # Extract total average treatment effect and standard error from unweighted model
    TATE_unadj <- TATE_model_unadj %>%
      broom::tidy() %>%
      dplyr::filter(term == treatment_indicator) %>%
      dplyr::pull(estimate)

    TATE_se_unadj <- TATE_model_unadj %>%
      broom::tidy() %>%
      dplyr::filter(term == treatment_indicator) %>%
      dplyr::pull(std.error)

    # Calculate 95% confidence interval for unweighted total average treatment effect
    TATE_CI_unadj <- TATE_unadj + 2.262*TATE_se_unadj*c(-1, 1)

    TATE_unadj <- list(estimate = TATE_unadj,
                       SE = TATE_se_unadj,
                       CI = TATE_CI_unadj)

    # Append TATE and TATE_unadj to list of items to return
    out <- c(out, list(TATE = TATE, TATE_unadj = TATE_unadj))
  }

  class(out) <- "generalizer_weighting"

  return(invisible(out))
}

# print.generalize_weighting <- function(x,...) {
#   cat("A generalizer_weighting object: \n")
#   cat(paste0(" - Outcome variable: ", x$outcome, "\n"))
#   cat(paste0(" - Sample indicator variable: ", x$sample_indicator, "\n"))
#   cat(paste0(" - Treatment indicator variable: ", x$treatment, "\n"))
#   cat(paste0(" - Covariates included: ", paste(x$covariates, collapse = ", "), "\n"))
#   cat(paste0(" - Probability of sample membership estimation method: ", x$estimation_method, "\n"))
#   cat(paste0(" - Are sample and population data considered disjoint?: ", ifelse(x$disjoint_data, "Yes", "No"), "\n"))
#   cat(paste0(" - Sample size: ", x$n_sample, "\n"))
#   cat(paste0(" - Population size : ", x$n_pop, "\n"))
#
#   invisible(x)
# }
#
# summary.generalize_weighting <- function(object,...) {
#   estimation_method_name = c("Logistic Regression", "Random Forest", "Lasso")
#   estimation_method = c("lr", "rf", "lasso")
#   prob_dist_table = rbind(summary(object$participation_probs$sample_indicator),
#                           summary(object$participation_probs$population))
#   row.names(prob_dist_table) = paste0(c("Sample","Population"), " (n = ", c(object$n_sample, object$n_pop),")")
#
#   selection_formula = paste0(object$sample_indicator," ~ ", paste(object$covariates, collapse = " + "))
#
#   out = list(
#     selection_formula = selection_formula,
#     estimation_method = estimation_method_name[estimation_method == object$estimation_method],
#     gen_index = object$gen_index,
#     prob_dist_table = prob_dist_table,
#     covariate_table = round(object$covariate_table, 4),
#     trim_pop = object$trim_pop,
#     n_excluded = object$n_excluded
#   )
#
#   class(out) <- "summary.generalize_weighting"
#   return(out)
# }
#
# print.summary.generalize_weighting <- function(x,...){
#   cat("Probability of Sample Participation: \n \n")
#   cat(paste0("Selection Model: ", x$selection_formula," \n \n"))
#   print(x$prob_dist_table)
#   cat("\n")
#   cat(paste0("Estimated by ", x$estimation_method, "\n"))
#   cat(paste0("Generalizability Index: ", gen_index, "\n"))
#   cat("============================================ \n")
#   if(x$trim_pop){
#     cat("Population data were trimmed for covariates to not exceed sample covariate bounds \n")
#     cat(paste0("Number excluded from population: ", x$n_excluded , "\n \n"))
#   }
#   cat("Covariate Distributions: \n \n")
#   print(round(x$covariate_table, 4))
#   invisible(x)
# }



