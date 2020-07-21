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

  # class(output) <- "generalizeAssess"
  return(output)

}
