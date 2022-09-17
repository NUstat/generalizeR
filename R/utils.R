#' Generate Sample Participation Probabilities
#'
#' This function is designed for use within 'weighting()' and 'assess()'.
#'
#' @param data data frame comprised of "stacked" sample and target population data
#' @param sample_var variable name denoting sample membership (1 = in sample, 0 = out of sample)
#' @param covariates vector of covariate names in data set that predict sample membership
#' @param estimation_method method to estimate the probability of sample membership. Default is logistic regression ("lr").Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @return sample participation probabilities for each unit in the data frame
#' @export
#' @importFrom glmnet cv.glmnet
#' @importFrom randomForest randomForest
#' @importFrom stats as.formula glm lm predict quantile

.generate.ps <- function(data,
                         sample_var,
                         covariates,
                         estimation_method) {

  # Logistic Regression
  if(estimation_method == "lr") {

    formula <- paste(sample_var,
                     paste(covariates, collapse = "+"),
                     sep = "~") %>%
      as.formula()

    ps <- formula %>%
      glm(data = data,
          family = "quasibinomial") %>%
      predict(type = "response")
  }

  # Random Forest
  if(estimation_method == "rf") {

    formula <- paste(
      paste("as.factor(", sample_var, ")"),
      paste(covariates, collapse = "+"),
      sep = "~") %>%
      as.formula()

    ps <- randomForest::randomForest(formula,
                                     data = data,
                                     na.action = na.omit) %>%
      #sampsize = 454,
      #ntree = 1500) %>%
      predict(type = "prob") %>%
      as.data.frame() %>%
      dplyr::pull("1")
  }

  # Lasso
  if(estimation_method == "lasso") {

    test.x <- model.matrix(~ -1 + .,
                           data = data %>% dplyr::select(tidyselect::all_of(covariates)))

    test.y <- data %>% dplyr::pull(sample_var)

    ps <- glmnet::cv.glmnet(x = test.x,
                            y = test.y,
                            family = "binomial") %>%
      predict(newx = test.x,
              s = "lambda.1se",
              type = "response") %>%
      as.numeric()
  }

  ### Set any participation probabilities of 0 in the sample to the minimum non-zero value ###

  if(0 %in% ps) {

    ps[which(data[,sample_var] == 1 & ps == 0)] <- min(ps[which(data[,sample_var] == 1 & ps != 0)], na.rm = TRUE)
  }

  return(ps)
}
