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
#' @param seed numeric. By default, the seed is set to 7835, otherwise can be specified (such as for simulation purposes).
#' @export
#' @importFrom glmnet cv.glmnet
#' @importFrom randomForest randomForest
#' @importFrom stats as.formula glm lm predict quantile
#' @importFrom crayon bold blue

weighting <- function(outcome, treatment, trial, selection_covariates, data,
                     selection_method = "lr", is_data_disjoint = TRUE, seed = 7835){

  weights <- NULL

  ##### set the seed #####
  set.seed(seed)

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
  # Logistic Regression
  if(selection_method == "lr") {

    formula <- paste(trial,
                   paste(selection_covariates, collapse = "+"),
                   sep = "~") %>%
      as.formula()

    ps <- formula %>%
      glm(data = data,
          family = "quasibinomial") %>%
      predict(type = "response")
  }

  # Random Forest
  if(selection_method == "rf") {

    formula <- paste(
      paste("as.factor(", trial, ")"),
      paste(selection_covariates, collapse = "+"),
      sep = "~") %>%
      as.formula()

    ps <- randomForest::randomForest(formula,
                                     data = data,
                                     na.action = na.omit) %>%
                                     #sampsize = 454,
                                     #ntree = 1500) %>%
      predict(type = "prob") %>%
      as.data.frame() %>%
      pull("1")
  }

  # Lasso
  if(selection_method == "lasso") {

    test.x <- model.matrix(~ -1 + .,
                           data = data %>% select(selection_covariates))

    test.y <- data %>% pull(trial)

    ps <- glmnet::cv.glmnet(x = test.x,
                            y = test.y,
                            family = "binomial") %>%
      predict(newx = test.x,
              s = "lambda.1se",
              type = "response") %>%
      as.numeric()
  }

  ### Set any participation probabilities of 0 in the trial to the minimum non-zero value ###

  if(0 %in% ps) {
    ps[which(data[,trial] == 1 & ps == 0)] <- min(ps[which(data[,trial] == 1 & ps != 0)], na.rm = TRUE)
  }

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

  participation_probs <- list(population = ps[which(data[,trial] == 0)],
                             trial = ps[which(data[,trial] == 1)])

  if(is.null(outcome) & is.null(treatment)) {TATE <- NULL}

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

  return(out)
}



