#' Unguided wrapper function for \code{stratify()}
#'
#' The function \code{stratify_unguided()} provides the unguided version of \code{stratify()}. It is intended to be called only within \code{stratify()} when users specify 'guided == FALSE', never as a standalone function.
#' Its primary purpose is to check the validity of the various inputs the user passed to \code{stratify()}.
#' @param data data.frame object containing the population data to be stratified.
#' @return
#' @details
#'
#' @md

stratify_unguided = function(data,
                             n_strata = NULL,
                             variables = NULL,
                             idvar = NULL) {

  # Check user arguments for invalid specifications -------------------------

  # n_strata
  if(is.null(n_strata)) {

    stop(simpleError("argument 'n_strata' is missing. You must specify the number of strata in which to divide your population if you are running the unguided version of this function."))
  }

  if(!is.numeric(n_strata) || n_strata <= 1 || n_strata %% 1 != 0) {

    stop(simpleError("invalid 'n_strata' argument. The number of strata must be an integer greater \nthan 1."))
  }

  # variables
  if(is.null(variables)) {

    stop(simpleError("argument 'variables' is missing. You must specify your stratification variables if you are running the unguided version of this function."))
  }

  if(!is.character(variables)) {

    stop(simpleError("invalid 'variables' argument. You must provide a character vector consisting of the names of the variables in your dataframe by which you wish to stratify your inference population."))
  }

  if(length(variables) <= 1) {

    stop(simpleError("invalid 'variables' argument. You must choose at least 2 variables from your dataframe for stratification."))
  }

  # idvar
  if(is.null(idvar)) {

    stop(simpleError("argument 'idvar' is missing. You must specify your ID variable if you are running the unguided version of this function."))
  }

  if(!is.character(idvar)) {

    stop(simpleError("invalid 'idvar' argument. The name of your ID variable must be a character string."))
  }

  # Check for invalid variable names ----------------------------------------
  invalid_vars = c(variables, idvar) %>%
    setdiff(names(data))

  if(!is_empty(invalid_vars)) {

    stop(simpleError(paste("The following variables are not columns in the dataframe you specified:\n",
                           paste(crayon::blue$bold(invalid_vars),
                                 collapse = ", "))))
  }

  # Check if there are any factor variables with more than 4 levels ---------

  source("R/stratify/helper functions/check_factor_levels.R")

  invalid_factors = data %>%
    select(all_of(variables)) %>%
    select_if(is.factor) %>%
    check_factor_levels()

  if(!is_empty(invalid_factors)) {

    stop(simpleError(

      paste0("This function will not allow a factor variable to have more than 4 levels.\n",
             "The following factor variables have more than 4 levels:\n",
             paste(crayon::blue$bold(invalid_factors), collapse = ", "),
             "\nPlease re-code your desired factor levels from these variables as dummy variables (see \nthe package 'fastDummies')."))

      )
  }
}
