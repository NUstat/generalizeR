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

  if (weighted_table) {

    data$weights <- ifelse(data[, sample_var] == 0, 1, data$weights)
  } else {

    data$weights <- 1
  }

  expanded.data <- data.frame(data[, sample_var],
                              model.matrix(~ -1 + ., data = data[, c(covariates, "weights")]))

  names(expanded.data)[1] <- sample_var

  if (!is_data_disjoint) {

    expanded.data <- expanded.data %>%
      dplyr::filter(in_study == 1) %>%
      rbind(expanded.data %>%
              dplyr::mutate(in_study = 0,
                            weights = 1))
  }

  get_covariate <- function(name) {

    covariate <- case_when(
      str_detect(name, "_mean_weighted$") ~  str_remove(name, "_mean_weighted$"),
      str_detect(name, "_mean$") ~ str_remove(name, "_mean$"),
      str_detect(name, "_var_weighted$") ~ str_remove(name, "_var_weighted$"),
      str_detect(name, "_var$") ~ str_remove(name, "_var$"),
    )

    return(covariate)
  }

  get_statistic <- function(name) {

    statistic <- case_when(
      str_detect(name, "_mean_weighted$") ~ "mean_weighted",
      str_detect(name, "_mean$") ~ "mean",
      str_detect(name, "_var_weighted$") ~ "var_weighted",
      str_detect(name, "_var$") ~ "var")

    return(statistic)
  }

  tab <- expanded.data %>%
    dplyr::group_by(!!sym(sample_var)) %>%
    dplyr::summarise(dplyr::across(tidyselect::all_of(covariates),
                                   list(mean = mean,
                                        mean_weighted = ~weighted.mean(., weights),
                                        var = var)))

  tab_pop <- tab %>%
    dplyr::filter(in_study == 0) %>%
    tidyr::pivot_longer(cols = -in_study) %>%
    dplyr::mutate(covariate = get_covariate(name),
                  statistic = get_statistic(name)
    ) %>%
    tidyr::pivot_wider(names_from = statistic,
                       values_from = value) %>%
    dplyr::group_by(covariate) %>%
    dplyr::summarise(dplyr::across(c("mean", "var"), na.omit)) %>%
    `colnames<-`(c("covariate", "pop_mean", "pop_var"))

  tab_sample <- tab %>%
    dplyr::filter(in_study == 1) %>%
    tidyr::pivot_longer(cols = -in_study) %>%
    dplyr::mutate(covariate = get_covariate(name),
                  statistic = get_statistic(name)
    ) %>%
    tidyr::pivot_wider(names_from = statistic,
                       values_from = value) %>%
    dplyr::group_by(covariate) %>%
    dplyr::summarise(across(c("mean", "mean_weighted"), na.omit)) %>%
    `colnames<-`(c("covariate", "sample_mean_unweighted", "sample_mean_weighted"))

  tab_merged <- merge(tab_pop, tab_sample, by = "covariate")  %>%
    dplyr::mutate(pop_sd = sqrt(pop_var),
                  ASMD_unweighted = abs((sample_mean_unweighted - pop_mean) / pop_sd),
                  ASMD_weighted = abs((sample_mean_weighted - pop_mean) / pop_sd)) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), round, digits = 3))



  if (weighted_table) {

    covariate_table <- tab_merged %>%
      dplyr::select(covariate, sample_mean_unweighted, sample_mean_weighted,
                    pop_mean, pop_sd, ASMD_unweighted, ASMD_weighted)


    covariate_kable <- covariate_table %>%
      dplyr::mutate(sample_mean = paste0(sample_mean_weighted,
                                         " [",
                                         sample_mean_unweighted,
                                         "]"),
                    ASMD = paste0(ASMD_weighted,
                                  " [",
                                  ASMD_unweighted,
                                  "]")) %>%
      dplyr::select(covariate, sample_mean, pop_mean, pop_sd, ASMD) %>%
      dplyr::rename(Covariate = covariate,
                    `Sample Mean` = sample_mean,
                    `Population Mean` = pop_mean,
                    `Population SD` = pop_sd) %>%
      kableExtra::kbl(caption = "Covariate Table",
                      align = "l") %>%
      kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE)

  } else {

    covariate_table <- tab_merged %>%
      dplyr::select(covariate, sample_mean_unweighted, pop_mean, pop_sd, ASMD_unweighted) %>%
      dplyr::rename(sample_mean = sample_mean_unweighted,
                    ASMD = ASMD_unweighted)

    covariate_kable <- covariate_table %>%
      dplyr::rename(Covariate = covariate,
                    `Sample Mean` = sample_mean,
                    `Population Mean` = pop_mean,
                    `Population SD` = pop_sd) %>%
      kableExtra::kbl(caption = "Covariate Table",
                      align = "l") %>%
      kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE)
  }

  return(list(covariate_table = covariate_table,
              covariate_kable = covariate_kable))
}
