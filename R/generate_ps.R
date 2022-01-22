#' Generate Sample Participation Probabilities
#'
#' This function is designed for use within 'weighting()' and 'assess()'.
#'
#' @param data data frame comprised of "stacked" trial and target population data
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param selection_method method to estimate the probability of trial participation. Default is logistic regression ("lr").Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @return sample participation probabilities for each unit in the data frame
#' @export
#' @importFrom glmnet cv.glmnet
#' @importFrom randomForest randomForest
#' @importFrom stats as.formula glm lm predict quantile

generate_ps <- function(data, trial, selection_covariates, selection_method) {

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

  return(ps)
}


