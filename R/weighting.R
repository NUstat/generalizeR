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
    stop("Data must be a data.frame.", call. = FALSE)
    }

  if(!outcome %in% data_names){
    stop(paste("The outcome variable", blue$bold(outcome), "is not a variable in the data provided!"), call. = FALSE)
    }

  if(!treatment %in% data_names){
    stop(paste("The treatment variable", blue$bold(treatment), "is not a variable in the data provided!"), call. = FALSE)
    }

  invalid_selection_covariates <- selection_covariates %>% setdiff(data_names)

  if(!is_empty(invalid_selection_covariates)){
    stop(paste("The following covariates are not variables in the data provided:\n", paste(blue$bold(invalid_selection_covariates), collapse = ", ")), call. = FALSE)
    }

  if(!selection_method %in% c("lr","rf","lasso")){
    stop("Invalid selection method!", call. = FALSE)
    }

  ### Omit rows in which the trial variable or the selection covariates contain missing values from the data ###
  data <- data %>% drop_na(c(trial, selection_covariates))

  ### Generate Participation Probabilities ###
  # Logistic Regression
  if(selection_method == "lr"){

    formula <- cat(trial,
                   paste(selection_covariates, collapse = "+"),
                   sep = "~") %>%
      as.formula()
    ps <- formula %>%
      glm(data = data,
          family = "quasibinomial") %>%
      predict(type = "response")

  }

  # Random Forests
  if(selection_method == "rf"){

    formula <- cat(
      paste("as.factor(", trial, ")"),
      paste(selection_covariates, collapse = "+"),
      sep = "~") %>%
      as.formula()
    ps <- randomForest::randomForest(formula,
                                     data = data,
                                     na.action = na.omit,
                                     sampsize = 454,
                                     ntree = 1500) %>%
      predict(type = "prob") %>%
      extract2(2)

  }

  # Lasso
  if(selection_method == "lasso"){

    test.x <- model.matrix(~ -1 + .,
                           data = data[,selection_covariates])
    test.y <- data[,trial]
    ps <- glmnet::cv.glmnet(x = test.x,
                           y = test.y,
                           family = "binomial") %>%
      predict(newx = test.x,
              s = "lambda.1se",
              type = "response") %>%
      as.numeric()

  }

  ### Set any participation probabilities of 0 in the trial to the minimum non-zero value ###
  if(any(ps[which(data[,trial]==1)] == 0)){
    ps[which(data[,trial] == 1 & ps == 0)] = min(ps[which(data[,trial] == 1 & ps != 0)], na.rm=TRUE)
  }

  ### Generate Weights ###
  if(is_data_disjoint == TRUE){
    data$weights = ifelse(data[,trial]==0,0,(1-ps)/ps)
  }

  if(is_data_disjoint == FALSE){
    data$weights = ifelse(data[,trial]==0,0,1/ps)
  }

  # Trim any of the weights if necessary
  data$weights[which(data$weights == 0 & data[,trial] == 1)] = quantile(data$weights[which(data[,trial]==1)],0.01,na.rm=TRUE)

  participation_probs = list(population = ps[which(data[,trial]==0)],
                             trial = ps[which(data[,trial]==1)])

  if(is.null(outcome) & is.null(treatment)){TATE = NULL}
  else{

    ##### ESTIMATE POPULATION AVERAGE TREATMENT EFFECT #####
    TATE_model = lm(as.formula(paste(outcome,treatment,sep="~")),data = data, weights = weights)

    TATE = summary(TATE_model)$coefficients[treatment,"Estimate"]
    TATE_se = summary(TATE_model)$coefficients[treatment,"Std. Error"]

    TATE_CI_l = TATE - 1.96*TATE_se
    TATE_CI_u = TATE + 1.96*TATE_se

    TATE = list(estimate = TATE, se = TATE_se, CI_l = TATE_CI_l, CI_u = TATE_CI_u)
  }

  ##### Items to return out #####
  out = list(participation_probs = participation_probs,
             weights = data$weights,
             TATE = TATE)

  return(out)
}
