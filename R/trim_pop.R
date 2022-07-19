#' Subset Population so Population Covariates are within bounds of Trial Covariates
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @return \code{trim_pop} returns a data frame, where the target population covariates do not exceed the bounds of the trial covariates

trim_pop <- function(trial, selection_covariates, data){

  ##### CHECKS #####
  if (!is.data.frame(data)) {
    stop("Data must be of type 'data.frame'.", call. = FALSE)
  }

  invalid_selection_covariates <- selection_covariates %>% setdiff(names(data))

  if(!is_empty(invalid_selection_covariates)) {
    stop(paste("The following covariates are not variables in the data provided:\n", paste(blue$bold(invalid_selection_covariates), collapse = ", ")), call. = FALSE)
  }

  trial_valid <- data %>%
    pull(trial) %>%
    na.omit() %>%
    table() %>%
    names() %>%
    setequal(c("0", "1"))

  if(!trial_valid) {
    stop("Trial membership variable must be binary and coded as `0` (not in trial) or `1` (in trial)", call. = FALSE)
  }

  ##### Subset trial data covariates #####
  trial_dat <- data %>%
    filter(trial == 1) %>%
    select(all_of(selection_covariates))

  if(length(selection_covariates) == 1) {

    trial_dat <- trial_dat %>% data.frame()
    names(trial_dat) <- selection_covariates
  }

  ##### Find covariate bounds in the trial #####
  covariate_bounds <- function(covariate) {

    # Convert quoted expression contained in function argument to symbol so it can be evaluated inside dplyr functions
    covariate <- covariate %>% rlang::sym()

    # Make covariate vector but only for observations selected to be part of the trial
    trial_covariate <- trial_dat %>% pull(covariate)

    if(trial_covariate %>% is.factor()) {

      trial_levels <- trial_covariate %>% droplevels() %>% levels()

      return(
        data %>%
          mutate(test = !(!!covariate %in% trial_levels)) %>% # The !!-operator (bang-bang) evaluates covariate first so the expression it contains is what gets passed to mutate()
          pull(test) %>%
          which()
      )
    }

    if(trial_covariate %>% is.numeric()) {

      trial_bounds <- c(trial_covariate %>% min(na.rm = TRUE), trial_covariate %>% max(na.rm = TRUE))

      return(
        data %>%
          mutate(test = !(!!covariate %>% between(trial_bounds[1], trial_bounds[2]))) %>%
          pull(test) %>%
          which()
      )
    }
  }

  ##### Find and remove rows of population data that violate bounds #####
  bound_violations <- purrr::map(selection_covariates, covariate_bounds)

  missing_rows <- bound_violations %>% unlist() %>% unique()

  trimmed_data <- data %>%
    filter(!row.names(data) %in% missing_rows) %>%
    droplevels() # Get rid of unused levels from factors

  ##### Number of rows in population data excluded #####
  n_excluded <- missing_rows %>% length()

  out <- list(n_excluded = n_excluded,
              trimmed_data = trimmed_data,
              untrimmed_data = data)

  return(out)
}
