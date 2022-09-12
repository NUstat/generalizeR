#' Estimate Weights for Generalizing Average Treatment Effect
#'
#' This function is designed for use within 'covariate_table()' and 'assess()'.
#'
#' @param outcome variable name denoting outcome
#' @param treatment variable name denoting binary treatment assignment (ok if only available in trial, not population)
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param selection_method method to estimate the probability of trial participation. Default is logistic regression ("lr").Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent. This affects calculation of the weights - see details for more information.
#' @export
#' @importFrom stats quantile
#' @importFrom crayon bold blue

weighting <- function(data,
                      trial,
                      treatment,
                      outcome,
                      selection_covariates,
                      selection_method = "lr",
                      is_data_disjoint = TRUE) {

  weights <- NULL

  ### Make input method lower case ###
  selection_method <- selection_method %>% tolower()

  ### Store the column names ###
  data_names <- data %>% names()

  ### Checks ###
  if (!is.data.frame(data)) {
    stop("Data must be an object of type 'data.frame'.", call. = FALSE)
    }

  if(!(outcome %in% data_names || is.null(outcome))) {
    stop(paste("The outcome variable", blue$bold(outcome), "is not a variable in the data provided."), call. = FALSE)
    }

  if(!(treatment %in% data_names || is.null(treatment))) {
    stop(paste("The treatment variable", blue$bold(treatment), "is not a variable in the data provided."), call. = FALSE)
    }

  invalid_selection_covariates <- selection_covariates %>% setdiff(data_names)

  if(!is_empty(invalid_selection_covariates)) {
    stop(paste("The following covariates are not variables in the data provided:\n", paste(blue$bold(invalid_selection_covariates), collapse = ", ")), call. = FALSE)
    }

  if(!selection_method %in% c("lr","rf","lasso")) {
    stop("Invalid selection method. Please choose one of 'lr' (logistic regression), 'rf' (random forest), or 'lasso' as your selection method.", call. = FALSE)
    }

  ### Omit rows in which the trial variable or the selection covariates contain missing values from the data ###
  data <- data %>% drop_na(c(trial, selection_covariates))

  ### Generate Participation Probabilities ###
  ps <- generate_ps(data, trial, selection_covariates, selection_method)

  participation_probs <- list(population = ps[which(data[,trial] == 0)],
                              trial = ps[which(data[,trial] == 1)])

  ### Generate Weights ###
  if(is_data_disjoint == TRUE) {

    data <- data %>%
      mutate(weights = ifelse(trial == 0,
                              0,
                              (1-ps)/ps
                              )
      )
  }

  else {

    data <- data %>%
      mutate(weights = ifelse(trial == 0,
                              0,
                              1/ps
                              )
      )
  }

  # Trim any of the weights if necessary
  data$weights[which(data$weights == 0 & data[,trial] == 1)] <- quantile(data$weights[which(data[,trial] == 1)], 0.01, na.rm = TRUE)

  # Add histogram of weights

  if(is.null(outcome) & is.null(treatment)) {TATE <- NULL}

  # SPLIT EVERYTHING BELOW INTO NEW FUNCTION

  else {

    ##### ESTIMATE POPULATION AVERAGE TREATMENT EFFECT #####

    # Model with weights (should outperform model without weights)
    TATE_model <- paste(outcome, treatment, sep = "~") %>%
      as.formula() %>%
      lm(data = data, weights = weights)

    # Model without weights
    TATE_model_null <- paste(outcome, treatment, sep = "~") %>%
      as.formula() %>%
      lm(data = data)

    # Total average treatment effect for model with weights
    TATE <- TATE_model %>%
      broom::tidy() %>%
      filter(term == "treatment") %>%
      pull(estimate)

    # Total average treatment effect for model without weights
    TATE_null <- TATE_model_null %>%
      broom::tidy() %>%
      filter(term == "treatment") %>%
      pull(estimate)

    # Standard error of total average treatment effect for model with weights
    TATE_se <- TATE_model %>%
      broom::tidy() %>%
      filter(term == "treatment") %>%
      pull(std.error)

    # Standard error of total average treatment effect for model without weights
    TATE_se_null <- TATE_model_null %>%
      broom::tidy() %>%
      filter(term == "treatment") %>%
      pull(std.error)

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



