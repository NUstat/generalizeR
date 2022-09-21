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
  if (estimation_method == "lr") {
    formula <- paste(sample_var,
      paste(covariates, collapse = "+"),
      sep = "~"
    ) %>%
      as.formula()

    ps <- formula %>%
      glm(
        data = data,
        family = "quasibinomial"
      ) %>%
      predict(type = "response")
  }

  # Random Forest
  if (estimation_method == "rf") {
    formula <- paste(
      paste("as.factor(", sample_var, ")"),
      paste(covariates, collapse = "+"),
      sep = "~"
    ) %>%
      as.formula()

    ps <- randomForest::randomForest(formula,
      data = data,
      na.action = na.omit
    ) %>%
      # sampsize = 454,
      # ntree = 1500) %>%
      predict(type = "prob") %>%
      as.data.frame() %>%
      dplyr::pull("1")
  }

  # Lasso
  if (estimation_method == "lasso") {
    test.x <- model.matrix(~ -1 + .,
      data = data %>% dplyr::select(tidyselect::all_of(covariates))
    )

    test.y <- data %>% dplyr::pull(sample_var)

    ps <- glmnet::cv.glmnet(
      x = test.x,
      y = test.y,
      family = "binomial"
    ) %>%
      predict(
        newx = test.x,
        s = "lambda.1se",
        type = "response"
      ) %>%
      as.numeric()
  }

  ### Set any participation probabilities of 0 in the sample to the minimum non-zero value ###

  if (0 %in% ps) {
    ps[which(data[, sample_var] == 1 & ps == 0)] <- min(ps[which(data[, sample_var] == 1 & ps != 0)], na.rm = TRUE)
  }

  return(ps)
}

#' Create Covariate Balance Table
#'
#' This function is designed for use within \code{weighting()} and \code{assess()}.'
#'
#' @param sample_var variable name denoting sample membership (1 = in sample, 0 = out of sample)
#' @param covariates vector of covariate names in data set that predict sample membership
#' @param data data frame comprised of "stacked" sample and target population data
#' @param weighted_table defaults to FALSE; whether weights are already included and do not need to be estimated
#' @param estimation_method method to estimate the probability of sample membership.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint defaults to TRUE. If TRUE, then sample and population data are considered independent.  This affects calculation of the weights
#' @export
#' @importFrom stats model.matrix weighted.mean
#' @importFrom dplyr funs

.make.covariate.table <- function(data,
                                  sample_var,
                                  covariates,
                                  weighted_table = FALSE,
                                  estimation_method = "lr",
                                  is_data_disjoint = TRUE) {

  data <- data %>%
    tidyr::drop_na(tidyselect::all_of(covariates)) %>%
    as.data.frame()

  if (!weighted_table) {

    expanded.data <- data.frame(
      data[, sample_var],
      model.matrix(~ -1 + ., data = data[, covariates])
    )

    names(expanded.data)[1] <- sample_var

    means.tab <- expanded.data %>%
      dplyr::group_by(!!sym(sample_var)) %>%
      dplyr::summarise(dplyr::across(tidyselect::all_of(covariates), mean)) %>%
      t() %>%
      as.data.frame() %>%
      .[-1, ] %>%
      dplyr::select(2,1)

    names(means.tab) <- c("sample", "population")

    n_sample <- as.numeric(table(expanded.data[, sample_var]))[2]
    n_pop <- as.numeric(table(expanded.data[, sample_var]))[1]
    sd.tab <- expanded.data %>%
      dplyr::group_by(!!sym(sample_var)) %>%
      dplyr::summarise(dplyr::across(everything(), var)) %>%
      t() %>%
      as.data.frame() %>%
      .[-1, ] %>%
      dplyr::select(2,1) %>%
      dplyr::mutate(pooled_sd = sqrt(((n_sample - 1) * V1 + (n_pop - 1) * V2) / (n_sample + n_pop - 2)))
    names(sd.tab) <- c("sample_var", "population_var", "pooled_sd")
    } else {

    data$weights <- ifelse(data[, sample_var] == 0, 1, data$weights)

    expanded.data <- data.frame(data[, sample_var],
                                model.matrix(~ -1 + ., data = data[, c(covariates, "weights")]))

    names(expanded.data)[1] <- sample_var

    means.tab <- expanded.data %>%
      dplyr::group_by(!!sym(sample_var)) %>%
      dplyr::summarise(dplyr::across(tidyselect::all_of(covariates), ~weighted.mean(., weights))) %>%
      t() %>%
      as.data.frame() %>%
      .[-1, ] %>%
      dplyr::select(2,1)

    names(means.tab) <- c("sample", "population")

    n_sample <- as.numeric(table(expanded.data[, sample_var]))[2]
    n_pop <- as.numeric(table(expanded.data[, sample_var]))[1]

    weighted.var <- function(x, w) {sum(w * (x - weighted.mean(x, w))^2) / (sum(w)-1)}

    sd.tab <- expanded.data %>%
      dplyr::group_by(!!sym(sample_var)) %>%
      dplyr::summarise(dplyr::across(tidyselect::all_of(covariates), ~weighted.var(., weights))) %>%
      t() %>%
      as.data.frame() %>%
      .[-1, ] %>%
      dplyr::select(2,1) %>%
      dplyr::mutate(pooled_sd = sqrt(((n_sample - 1) * V1 + (n_pop - 1) * V2) / (n_sample + n_pop - 2)))

    names(sd.tab) <- c("sample_var", "population_var", "pooled_sd")
  }

  if (!is_data_disjoint) {

    means.tab <- means.tab %>%
      dplyr::mutate(population = expanded.data %>%
                      dplyr::summarise(dplyr::across(tidyselect::all_of(covariates), mean)) %>%
                      t() %>%
                      as.data.frame() %>%
                      pull())

    n_pop <- nrow(data)

    sd.tab <- sd.tab %>%
      dplyr::mutate(population_var = expanded.data %>%
                      dplyr::summarise(dplyr::across(tidyselect::all_of(covariates), var)) %>%
                      t() %>%
                      as.data.frame() %>%
                      pull(),
                    pooled_sd = sqrt(((n_sample - 1) * sample_var + (n_pop - 1) * population_var) / (n_sample + n_pop - 2)))
  }

  covariate_table <- means.tab %>%
    dplyr::bind_cols(sd.tab) %>%
    dplyr::mutate(ASMD = round(abs((sample - population) / pooled_sd), 3),
                  sample_sd = sqrt(sample_var),
                  population_sd = sqrt(population_var)) %>%
    dplyr::select(
      `Sample Mean` = sample,
      `Population Mean` = population,
      `Sample SD` = sample_sd,
      `Population SD` = population_sd,
      ASMD
    ) %>%
    dplyr::arrange(desc(ASMD)) %>%
    round(digits = 3) %>%
    tibble::rownames_to_column("Covariate")

  if (weighted_table) {

    names(covariate_table)[c(2,4)] <- c("Weighted Sample Mean", "Weighted Sample SD")
  }

  covariate_kable <- covariate_table %>%
    kableExtra::kbl(caption = "Covariate Table",
                    align = "l") %>%
    kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE)

  return(list(covariate_table = covariate_table,
              covariate_kable = covariate_kable))
}
