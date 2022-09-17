#' Estimate Weights for Generalizing Average Treatment Effect
#'
#' This function is designed for use within 'covariate_table()' and 'assess()'.
#'
#' @param outcome variable name denoting outcome
#' @param treatment variable name denoting binary treatment assignment (ok if only available in trial, not population)
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param estimation_method method to estimate the probability of trial participation. Default is logistic regression ("lr").Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent. This affects calculation of the weights - see details for more information.
#' @export
#' @importFrom stats quantile
#' @importFrom crayon bold blue

weighting <- function(data,
                      sample_indicator,
                      treatment_indicator = NULL,
                      outcome = NULL,
                      covariates,
                      estimation_method = "lr",
                      is_data_disjoint = TRUE) {

  weights <- NULL

  ### Make input method lower case ###
  estimation_method <- tolower(estimation_method)

  ### Store the column names ###
  data_names <- names(data)

  # Check whether data object is of type 'data.frame'
  assertthat::on_failure(is.data.frame) <- function(call, env) {
    "You must pass an object of type 'data.frame' to the 'data' argument."
  }

  assertthat::assert_that(is.data.frame(data))

  # Check whether outcome variable is one of columns in dataframe
  is_outcome_valid <- function(outcome) {

    !is.null(outcome) && outcome %in% data_names
  }

  assertthat::on_failure(is_outcome_valid) <- function(call, env) {

      "Your outcome variable must be one of the columns in the dataframe you have provided."
  }

  # Check whether treatment indicator variable is one of columns in dataframe
  is_treatment_indicator_valid <- function(treatment_indicator) {

    !is.null(treatment_indicator) && treatment_indicator %in% data_names
  }

  assertthat::on_failure(is_treatment_indicator_valid) <- function(call, env) {

    "Your treatment indicator variable must be one of the columns in the dataframe you have provided."
  }

  assertthat::assert_that(is_outcome_valid(outcome))

  invalid_covariates <- covariates %>% setdiff(data_names)

  ##### Ensure selection covariates are variables in the dataframe provided #####
  assertthat::on_failure(is_empty) <- function(call, env) {

    paste("The following covariates are not variables in the data provided:\n",
          paste(crayon::blue$bold(invalid_covariates),
                collapse = ", ")
    )
  }

  assertthat::assert_that(is_empty(invalid_covariates))

  ##### Ensure estimation method is valid #####
  is_estimation_method_valid <- function(estimation_method){

    estimation_method %in% c("lr","rf","lasso")
  }

  assertthat::on_failure(is_estimation_method_valid) <- function(call, env) {

    "Invalid estimation method. Please choose one of 'lr' (logistic regression), 'rf' (random forest), or 'lasso' as your estimation method."
  }

  ### Omit rows in which the sample indicator variable or the covariates contain missing values from the data ###
  data <- data %>%
    tidyr::drop_na(c(tidyselect::all_of(sample_indicator), tidyselect::all_of(covariates)))

  ### Generate Participation Probabilities ###
  ps <- .generate.ps(data, sample_indicator, covariates, estimation_method)

  participation_probs <- list(population = ps[which(data[,sample_indicator] == 0)],
                              sample = ps[which(data[,sample_indicator] == 1)])

  ### Generate Weights ###
  if(is_data_disjoint == TRUE) {

    data <- data %>%
      dplyr::mutate(weights = ifelse(sample_indicator == 0,
                                     0,
                                     (1-ps)/ps))
  }

  else {

    data <- data %>%
      dplyr::mutate(weights = ifelse(sample_indicator == 0,
                                     0,
                                     1/ps))
  }

  # Trim any of the weights if necessary
  data$weights[which(data$weights == 0 & data[,sample_indicator] == 1)] <- quantile(data$weights[which(data[,sample_indicator] == 1)], 0.01, na.rm = TRUE)

  # Add histogram of weights

  if(is.null(outcome) & is.null(treatment_indicator)) {TATE <- NULL}

  # SPLIT EVERYTHING BELOW INTO NEW FUNCTION

  else {

    ##### ESTIMATE POPULATION AVERAGE TREATMENT EFFECT #####

    # Model with weights (should outperform model without weights)
    TATE_model <- paste(outcome, treatment_indicator, sep = "~") %>%
      as.formula() %>%
      lm(data = data, weights = weights)

    # Model without weights
    TATE_model_null <- paste(outcome, treatment_indicator, sep = "~") %>%
      as.formula() %>%
      lm(data = data)

    # Total average treatment effect for model with weights
    TATE <- TATE_model %>%
      broom::tidy() %>%
      dplyr::filter(term == "treatment") %>%
      dplyr::pull(estimate)

    # Total average treatment effect for model without weights
    TATE_null <- TATE_model_null %>%
      broom::tidy() %>%
      dplyr::filter(term == "treatment") %>%
      dplyr::pull(estimate)

    # Standard error of total average treatment effect for model with weights
    TATE_se <- TATE_model %>%
      broom::tidy() %>%
      dplyr::filter(term == "treatment") %>%
      dplyr::pull(std.error)

    # Standard error of total average treatment effect for model without weights
    TATE_se_null <- TATE_model_null %>%
      broom::tidy() %>%
      dplyr::filter(term == "treatment") %>%
      dplyr::pull(std.error)

    # 95% confidence interval for total average treatment effect of model with weights
    TATE_CI <- TATE + 2.262*TATE_se*c(-1, 1)

    # 95% confidence interval for total average treatment effect of model without weights
    TATE_CI_null <- TATE_null + 2.262*TATE_se_null*c(-1, 1)

    TATE <- list(estimate = TATE,
                estimate_null = TATE_null,
                se = TATE_se,
                se_null = TATE_se_null,
                CI = TATE_CI,
                CI_null = TATE_CI_null)
  }

  ##### Items to return out #####
  out <- list(participation_probs = participation_probs,
             weights = data$weights,
             TATE = TATE)

  ### Add weighted covariate table

  return(out)
}



