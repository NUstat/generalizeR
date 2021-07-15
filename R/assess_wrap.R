#' Assess the Generalizability of a Sample to a Population
#'
#' This function is a wrapper for the function 'assess()' which may make assessment more user-friendly and flexible.
#'
#' Given a data frame or tibble consisting of a column of some ID variable (the sample) and a data frame or tibble of the same ID variable and stratifying variables of interest, this function will assess the generalizability of the sample to the population.
#'
#' The user can also specify a grouping variable to assess the generalizability of their sample across population groups -- for example, across US states. This will produce multiple generalizability index values, one per level of the grouping variable.
#'
#' @param sample a data frame with a column of IDs
#' @param population a data frame with a column of IDs, columns for any stratifying variables of interest, and (if applicable) a grouping variable
#' @param join_var defaults to NULL; name of the ID variable in double quotes
#' @param grouping_var defaults to NULL; name of a grouping variable in double quotes (not required)
#' @export
#' @importFrom dplyr full_join group_by_at group_map
#' @importFrom tidyr replace_na
#' @importFrom purrr map

assess_wrap <- function(sample, population, join_var = NULL, grouping_var = NULL){

  trial <- NULL

  # Ideally, want this to work if there are no ids common to both; an option to check where the stratifying variables just have ot be the same. add this.
  sample <- sample %>%
    clean_names() %>%
    mutate(trial = rep(1))

  sample_and_pop <- sample %>%
    full_join((population %>% clean_names()), by = join_var) %>%
    replace_na(list(trial = 0)) %>%
    unique() %>%
    na.omit()

  sample_overall <- sample_and_pop %>%
    filter(trial == 1)

  population_overall <- sample_and_pop %>%
    filter(trial == 0)

  if(is.null(grouping_var)){
    selection_vars <- colnames(
      (sample_and_pop %>% select(-all_of(join_var), -trial))
    )

    full_data <- add_row(sample_overall, population_overall)

    output <- assess(data = full_data, trial = "trial",
              is_data_disjoint = TRUE, selection_covariates = selection_vars)
  }
  if(!is.null(grouping_var)){
    selection_vars <- colnames(
      (sample_and_pop %>% select(-all_of(join_var), -trial, -all_of(grouping_var)))
    )

    output <- population_overall %>%
      group_by_at(grouping_var) %>%
      group_map(~ add_row(sample_overall, .x)) %>% # For each group, add the sample data
      map(~ suppressWarnings(tryCatch(assess(data = data.frame(.x),
                                                   trial = "trial",
                                                   selection_covariates = selection_vars),
                                            error = function(err) NA)))

    g_indexes <- unlist(map(output, function(x){x["g_index"][[1]]}))
    cov_matrices <- map(output, function(x){x["covariate_table"][[1]]})

    output <- list(output, g_indexes, cov_matrices)

  }

  # class(output) <- "generalize_assess"
  return(output)

}

print.generalize_assess <- function(x,...){
  cat("A generalize_assess object: \n")
  cat(paste0(" - probability of trial participation method: ", x$selection_method, "\n"))
  cat(paste0(" - common covariates included: ", paste(x$selection_covariates, collapse = ", "), "\n"))
  cat(paste0(" - sample size of trial: ", x$n_trial, "\n"))
  cat(paste0(" - size of population: ", x$n_pop, "\n"))
  cat(paste0(" - was population trimmed according to trial covariate bounds?: ", ifelse(x$trim_pop == TRUE, "Yes", "No"), "\n"))
  if(x$trim_pop == TRUE){
    cat(paste0("    - number excluded from population data: ", x$n_excluded, "\n"))
  }

  invisible(x)
}

summary.generalize_assess <- function(object,...){
  selection_method_name = c("Logistic Regression","Random Forests","Lasso")
  selection_method = c("lr","rf","lasso")
  prob_dist_table = rbind(summary(object$participation_probs$trial),
                          summary(object$participation_probs$population))
  row.names(prob_dist_table) = paste0(c("Trial","Population"), " (n = ", c(object$n_trial,object$n_pop),")")

  selection_formula = paste0(object$trial_name," ~ ",paste(object$selection_covariates, collapse = " + "))

  out = list(
    selection_formula = selection_formula,
    selection_method = selection_method_name[selection_method == object$selection_method],
    g_index = round(object$g_index,3),
    prob_dist_table = prob_dist_table,
    covariate_table = round(object$covariate_table, 4),
    weighted_covariate_table = round(object$weighted_covariate_table,4),
    trim_pop = object$trim_pop,
    n_excluded = object$n_excluded
  )

  class(out) = "summary.generalize_assess"
  return(out)
}

print.summary.generalize_assess <- function(x,...){
  cat("Probability of Trial Participation: \n \n")
  cat(paste0("Selection Model: ",x$selection_formula," \n \n"))
  print(x$prob_dist_table)
  cat("\n")
  cat(paste0("Estimated by ",x$selection_method, "\n"))
  cat(paste0("Generalizability Index: ", round(x$g_index,3), "\n"))
  cat("============================================ \n")
  cat("Covariate Distributions: \n \n")
  if(x$trim_pop == TRUE){
    cat("Population data were trimmed for covariates to not exceed trial covariate bounds \n")
    cat(paste0("Number excluded from population: ", x$n_excluded ,"\n \n"))
  }
  print(round(x$covariate_table,4))
  invisible(x)
}

