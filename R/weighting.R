#' Estimate Weights for Generalizing Average Treatment Effect
#'
#' @param data data frame comprised of "stacked" sample and target population data
#' @param sample_indicator variable name denoting binary sample membership (1 = in sample, 0 = out of sample)
#' @param treatment_indicator variable name denoting binary treatment assignment (ok if only available in sample, not population)
#' @param outcome variable name denoting outcome
#' @param covariates vector of covariate names in data set that predict sample membership
#' @param estimation_method method to estimate the probability of sample membership (propensity scores). Default is logistic regression ("lr").Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param disjoint_data logical. If TRUE, then sample and population data are considered disjoint. This affects calculation of the weights - see details for more information.
#' @return A summary of propensity scores, covariates, and ASMD for both weighted and unweighted data, as well as a summary of the weights. Also weighted and unweighted TATE if outcome and treatment are given
#' @examples
#' library(tidyverse)
#'
#' # creating a stratified sample and recruiting from the sample to prepare for assess
#' selection_covariates <- c("total", "pct_black_or_african_american", "pct_white",
#'                           "pct_female", "pct_free_and_reduced_lunch")
#'
#' strat_output <- stratify(generalizeR:::inference_pop, guided = FALSE, n_strata = 4,
#'                          variables = selection_covariates, idvar = "ncessch")
#' rec_output <- recruit(strat_output, guided = FALSE, sample_size = 40)
#'
#' # creating the sample dataset from the output of recruit
#' sample_list <- c(rec_output$recruitment_lists[[1]]$ncessch[1:5],
#'                   rec_output$recruitment_lists[[2]]$ncessch[1:20],
#'                   rec_output$recruitment_lists[[3]]$ncessch[1:11],
#'                   rec_output$recruitment_lists[[4]]$ncessch[1:4])
#' inference_pop_sample <- mutate(generalizeR:::inference_pop,
#'                                sample = if_else(ncessch %in% sample_list, 1, 0))
#'
#' # weighting the sample with the given covariates
#' weighting_output <- weighting(inference_pop_sample, sample_indicator = "sample",
#'                               covariates = selection_covariates, disjoint_data = FALSE)
#'
#'
#'
#' @export
#' @importFrom stats quantile
#' @importFrom crayon bold blue
#' @importFrom rlang is_empty


weighting <- function(data,
                      sample_indicator,
                      treatment_indicator = NULL,
                      outcome = NULL,
                      covariates,
                      estimation_method = "lr",
                      disjoint_data = TRUE) {

  estimation_method <- tolower(estimation_method)

  # Check whether data object is of type 'data.frame'
  assertthat::on_failure(is.data.frame) <- function(call, env) {

    "You must pass an object of type 'data.frame' to the 'data' argument."
  }

  assertthat::assert_that(is.data.frame(data))

  data_name <- data %>%
    lazyeval::expr_text()

  data_names <- names(data)

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
  is_estimation_method_valid <- function(estimation_method) {

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

  propensity_scores <- list(not_in_sample = data$ps[which(data[,sample_indicator] == 0)],
                            in_sample = data$ps[which(data[,sample_indicator] == 1)],
                            population = data$ps)

  # Generate Weights #
  if(disjoint_data) {

    data <- data %>%
      dplyr::mutate(sample_weights = ifelse(!!sym(sample_indicator) == 0,
                                            0,
                                            (1-ps)/ps))
  }

  else {

    data <- data %>%
      dplyr::mutate(sample_weights = ifelse(!!sym(sample_indicator) == 0,
                                            0,
                                            1/ps))
  }

  # Trim any of the weights if necessary
  data$sample_weights[which(data$sample_weights == 0 & data[,sample_indicator] == 1)] <- quantile(data$sample_weights[which(data[,sample_indicator] == 1)], 0.01, na.rm = TRUE)

  # Make weighted covariate table
  covariate_table_output <- .make.covariate.table(data = data,
                                                  sample_indicator = sample_indicator,
                                                  covariates = covariates,
                                                  sample_weights = "sample_weights",
                                                  estimation_method = estimation_method,
                                                  disjoint_data = disjoint_data)

  # Make histogram of weights
  weights_hist <- data %>%
    dplyr::filter(!!sym(sample_indicator) == 1) %>%
    ggplot(aes(x = sample_weights)) +
    geom_histogram(bins = 20,
                   fill = viridis::viridis(1, alpha = 0.7),
                   color = "black") +
    theme_minimal() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = "Weights",
         y = "Frequency",
         title = "Histogram of Sample Weights")

  # Calculate sample and population sizes

  n_sample <- data %>%
    dplyr::filter(!!rlang::sym(sample_indicator) == 1) %>%
    nrow()

  n_pop <- data %>%
    nrow()

  if (disjoint_data) {n_pop <- n_pop - n_sample}

  # Items to return out
  out <- list(propensity_scores = propensity_scores,
              sample_weights = data$sample_weights,
              covariate_table = covariate_table_output$covariate_table,
              covariate_kable = covariate_table_output$covariate_kable,
              cov_dist_facet_plot = covariate_table_output$cov_dist_facet_plot,
              cov_dist_plots = covariate_table_output$cov_dist_plots,
              weights_hist = weights_hist,
              data_name = data_name,
              estimation_method = estimation_method,
              disjoint_data = disjoint_data,
              n_sample = n_sample,
              n_pop = n_pop)

  if(!is.null(outcome) & !is.null(treatment_indicator)) {

    # ESTIMATE TOTAL AVERAGE TREATMENT EFFECT

    # ADJUSTED TATE

    # Make weighted regression model predicting outcome with treatment
    TATE_model <- paste(outcome, treatment_indicator, sep = "~") %>%
      as.formula() %>%
      lm(data = data, weights = sample_weights)

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
    TATE_CI <- round(TATE + 2.262*TATE_se*c(-1, 1), 3)

    TATE <- list(model = TATE_model,
                 estimate = TATE,
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
    TATE_CI_unadj <- round(TATE_unadj + 2.262*TATE_se_unadj*c(-1, 1), 3)

    TATE_unadj <- list(model = TATE_model_unadj,
                       estimate = TATE_unadj,
                       SE = TATE_se_unadj,
                       CI = TATE_CI_unadj)

    # Summarize results in table
    TATE_table <- data.frame(TATE = c(TATE$estimate, TATE_unadj$estimate),
                            SE = c(TATE$SE, TATE_unadj$SE),
                            "95% CI" = c(paste0("(", TATE$CI[1], ", ", TATE$CI[2], ")"),
                                         paste0("(", TATE_unadj$CI[1], ", ", TATE_unadj$CI[2], ")")),
                            row.names = c("Weighted", "Unweighted"),
                            check.names = FALSE) %>%
      colorDF::colorDF(theme = "dark")

    print(TATE_table)

    # Append outcome variable, treatment variable, TATE, TATE_unadj, and TATE_table to list of items to return
    out <- c(out, list(outcome = outcome,
                       treatment = treatment_indicator,
                       TATE = TATE,
                       TATE_unadj = TATE_unadj,
                       TATE_table = TATE_table))
  }

  class(out) <- "generalizeR_weighting"

  return(invisible(out))
}

#' Print method for "generalizeR_weighting" class
#'
#' @param x An object of class "generalizeR_weighting"
#' @param ... Other arguments passed to or from other methods
#' @return A summary of the weighted dataset
#'
#' @export print.generalizeR_weighting
#' @export

print.generalizeR_weighting <- function(x, ...) {

  cat("\nA generalizeR_weighting object: \n\n")

  cat(" - Dataset name:", crayon::cyan$bold(x$data_name), "\n\n")

  if (!is.null(x$outcome) & !is.null(x$treatment)) {

    cat(" - Outcome variable:", crayon::cyan$bold(x$outcome), "\n\n")

    cat(" - Treatment variable:", crayon::cyan$bold(x$treatment), "\n\n")
  }

  covariate_names <- x$covariate_table %>%
    pull(covariate) %>%
    crayon::cyan$bold() %>%
    paste(collapse = ", ") %>%
    gsub('(.{200})\\s(,*)', '\\1\n   \\2', .)

  cat(" - Covariates selected:\n\n  ", covariate_names, "\n\n")

  cat(" - Method used to estimate propensity scores:",
      switch(x$estimation_method,
             "lr" = "Logistic Regression",
             "rf" = "Random Forest",
             "lasso" = "LASSO") %>%
        crayon::cyan$bold(),
      "\n\n")

  if (x$disjoint_data) {cat(paste0(" - The sample and the population were considered ", crayon::cyan$bold("disjoint"), " from one another.\n\n"))}

  else {cat(paste0(" - The sample was considered a ", crayon::cyan$bold("subset"), " of the population.\n\n"))}

  cat(" - Sample size:", crayon::cyan$bold(x$n_sample), "\n\n")

  cat(" - Population size:", crayon::cyan$bold(x$n_pop), "\n\n")
}

#' Summary method for "generalizeR_weighting" class
#'
#' @param object An object of class "generalizeR_weighting"
#' @param ... Other arguments passed to or from other methods
#' @return A summary of the weighted and unweighted dataset for propensity scores, covariates and (optionally) TATE
#'
#' @export summary.generalizeR_weighting
#' @export

summary.generalizeR_weighting <- function(object, ...) {

  estimation_method <- switch(object$estimation_method,
                              "lr" = "Logistic Regression",
                              "rf" = "Random Forest",
                              "lasso" = "Lasso")

  if (object$disjoint_data) {

    prop_score_dist_table <- rbind(summary(object$propensity_scores$in_sample),
                                   summary(object$propensity_scores$population))

    row.names(prop_score_dist_table) <- paste0(c("Sample","Population"), " (n = ", c(object$n_sample, object$n_pop),")")

    sample_prop_scores <- data.frame(prop_scores = object$propensity_scores$in_sample) %>%
      mutate(sample_indicator = 1)

    pop_prop_scores <- data.frame(prop_scores = object$propensity_scores$population) %>%
      mutate(sample_indicator = 0)

    prop_scores <- rbind(sample_prop_scores, pop_prop_scores)

    prop_score_dist_plot <- prop_scores %>%
      ggplot2::ggplot() +
      geom_density(aes(x = gtools::logit(prop_scores), fill = factor(sample_indicator)),
                   alpha = 0.7) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_discrete(name = NULL,
                          labels = c("Population", "Sample")) +
      labs(x = "Propensity Score Logits",
           y = "Density",
           title = "Distribution of Propensity Score Logits") +
      theme_minimal() +
      theme(axis.ticks.x = element_line(),
            axis.text.y = element_blank(),
            axis.line = element_line(),
            plot.title = element_text(size = 12))
  }

  else {

    prop_score_dist_table <- rbind(summary(object$propensity_scores$in_sample),
                                   summary(object$propensity_scores$not_in_sample))

    row.names(prop_score_dist_table) <- paste0(c("In Sample","Not In Sample"), " (n = ", c(object$n_sample, object$n_pop),")")

    in_sample_prop_scores <- data.frame(prop_scores = object$propensity_scores$in_sample) %>%
      mutate(sample_indicator = 1)

    not_in_sample_prop_scores <- data.frame(prop_scores = object$propensity_scores$not_in_sample) %>%
      mutate(sample_indicator = 0)

    prop_scores <- rbind(in_sample_prop_scores, not_in_sample_prop_scores)

    prop_score_dist_plot <- prop_scores %>%
      ggplot2::ggplot() +
      geom_density(aes(x = prop_scores, fill = factor(sample_indicator)),
                   alpha = 0.7) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_discrete(name = NULL,
                          labels = c("Not in Sample", "In Sample")) +
      labs(x = "Probability",
           y = "Density",
           title = "Distribution of Propensity Scores") +
      theme_minimal() +
      theme(axis.ticks.x = element_line(),
            axis.text.y = element_blank(),
            axis.line = element_line(),
            axis.title = element_blank(),
            plot.title = element_text(size = 12))
  }

  colnames(prop_score_dist_table) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")

  prop_score_dist_table <- prop_score_dist_table %>%
    data.frame() %>%
    colorDF::colorDF(theme = "dark")

  out <- list(estimation_method = estimation_method,
              prop_score_dist_table = prop_score_dist_table,
              prop_score_dist_plot = prop_score_dist_plot,
              covariate_table = object$covariate_table,
              weights_hist = object$weights_hist,
              weighted_model = object$TATE$model,
              unweighted_model = object$TATE_unadj$model,
              TATE_table = object$TATE_table)

  class(out) <- "summary.generalizeR_weighting"

  return(out)
}

#' Print method for "summary.generalizeR_weighting" class
#'
#' @param x An object of class "summary.generalizeR_weighting"
#' @param ... Other arguments passed to or from other methods
#' @return A summary of the weighted and unweighted dataset for propensity scores, covariates and (optionally) TATE
#'
#' @export print.summary.generalizeR_weighting
#' @export

print.summary.generalizeR_weighting <- function(x, ...) {

  cat("\nSummary of Estimated Propensity Scores: \n\n")

  print(x$prop_score_dist_table)

  print(x$prop_score_dist_plot)

  cat("\n\nEstimation Method:", crayon::cyan$bold(x$estimation_method), "\n\n")

  covariate_names <- x$covariate_table %>%
    pull(covariate) %>%
    crayon::cyan$bold() %>%
    paste(collapse = ", ") %>%
    gsub('(.{200})\\s(,*)', '\\1\n   \\2', .)

  cat("Covariates Used:\n\n  ", covariate_names, "\n\n")

  cat("Covariate Table: \n\n")

  print(x$covariate_table)

  print(x$weights_hist)

  cat("\nSummary of Unweighted Regression Model: \n")

  summary(x$unweighted_model) %>% print()

  cat("\nSummary of Weighted Regression Model: \n")

  summary(x$weighted_model) %>% print()

  cat("TATE Table: \n\n")

  print(x$TATE_table)
}

if(getRversion() >= "2.15.1") utils::globalVariables(c("covariate", ".", "geom_density",
                                                       "sample_indicator", "scale_x_continuous", "scale_y_continuous",
                                                       "scale_fill_discrete", "theme_minimal", "element_line", ".", "sym", "ps",
                                                       "sample_weights", "term", "estimate", "std.error"))


