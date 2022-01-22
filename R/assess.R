#' Assess Generalizability of Randomized Trial to Population
#'
#' This function, given a stacked data frame containing both sample and population data, assesses the generalizability of the sample to the population on given covariates.
#'
#' 'assess_wrap()' is a wrapper for this function that allows assessment over levels of a grouping variable.
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param selection_method method to estimate the probability of trial participation. Default is logistic regression ("lr"). Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint logical. If TRUE, then trial and population data are considered independent. This affects calculation of the weights - see details for more information.
#' @param trim_pop logical. If TRUE, then population data are subset to exclude individuals with covariates outside bounds of trial covariates.
#' @param seed numeric. By default, the seed is set to 1996, otherwise can be specified (such as for simulation purposes).
#' @export


assess <- function(data, guided = TRUE, trial, selection_covariates, selection_method = "lr",
                   is_data_disjoint = TRUE, trim_pop = FALSE, seed = 1996) {

  ##### Set the seed #####
  set.seed(seed)

  ##### Ensure 'data' is actually a dataframe #####
  if(!is.data.frame(data)) {
    stop("Data must be an object of type 'data.frame'.", call. = FALSE)
  }

  ##### GUIDED VERSION PART 1 #####

  if(guided == TRUE) {

    ##### Store the variable names #####
    var_names <- names(data)

    satisfied <- FALSE

    while(satisfied == FALSE) {

      cat("Please select the variable coding binary trial participation in your dataframe: \n")
      trial <- select.list(choices = var_names,
                           graphics = FALSE,
                           multiple = FALSE)

      if(anyNA(match(names(table(data[,trial])), c("0","1")))) {

        cat(red("The trial variable must be coded as '0' (not in trial) or '1' (in trial).\n\n"))
        next
      }

      satisfied <- TRUE
    }

    cat("Please select the covariates that will be used to predict trial participation: \n")
    selection_covariates <- select.list(choices = var_names %>% setdiff(trial),
                                        graphics = FALSE,
                                        multiple = TRUE)

    cat("Are the trial data and population data disjoint? See the vignette for more details.\n")
    is_data_disjoint <- menu(choices = c("Yes", "No"),
                     graphics = FALSE) %>%
      switch("1" = TRUE,
             "2" = FALSE)

    cat("Should the population data be trimmed to exclude units with covariate values outside \nthe bounds of the trial covariates? See the vignette for more details. \n")
    trim_pop <- menu(choices = c("Yes", "No"),
                            graphics = FALSE) %>%
      switch("1" = TRUE,
             "2" = FALSE)

    cat("Please select the method that will be used to estimate the probability of trial participation: \n")
    selection_method <- select.list(choices = c("Logistic Regression", "Random Forest", "Lasso"),
                                    graphics = FALSE,
                                    multiple = FALSE) %>%
      switch("Logistic Regression" = "lr",
             "Random Forest" = "rf",
             "Lasso" = "lasso")
  }

  ##### NON-GUIDED VERSION PART 1 #####
  else {

    ##### Make methods lower case #####
    selection_method = tolower(selection_method)

    ##### Ensure selection covariates are part of dataframe #####
    invalid_selection_covariates <- selection_covariates %>% setdiff(names(data))

    if(!is_empty(invalid_selection_covariates)) {

      stop(paste("The following covariates are not variables in the data provided:\n", paste(blue$bold(invalid_selection_covariates), collapse = ", ")), call. = FALSE)
    }

    ##### Ensure trial variable is binary #####
    if(anyNA(match(names(table(data[,trial])), c("0","1")))) {

      stop("Sample membership variable must be coded as `0` (not in trial) or `1` (in trial).", call. = FALSE)
    }

    ##### Ensure selection method is valid #####
    if(!selection_method %in% c("lr","rf","lasso")) {

      stop("Invalid selection method. Please choose one of 'lr' (logistic regression), 'rf' (random forest), or 'lasso' as your selection method.", call. = FALSE)
    }
  }

  ##### Drop missing values from trial and selection covariate columns in data #####
  data <- data %>%
    drop_na(trial, all_of(selection_covariates))

  if(trim_pop == FALSE) {n_excluded <- 0}

  else {

    trim_pop_output <- trim_pop(trial, selection_covariates, data)

    n_excluded <- trim_pop_output$n_excluded

    data <- trim_pop_output$trimmed_data
  }

  ##### Generate Participation Probabilities #####
  ps <- generate_ps(data, trial, selection_covariates, selection_method)

  participation_probs <- list(nottrial = ps[which(data[,trial] == 0)],
                              trial = ps[which(data[,trial] == 1)],
                              population = ps)

  ##### Calculate Generalizability Index  #####

  ## If data is not disjoint, compare trial to (nottrial + trial)
  if(is_data_disjoint == FALSE) {

    g_index <- gen_index(participation_probs$trial, participation_probs$population) %>% round(4)
  }

  ## If data is disjoint, compare trial to nottrial
  else {

    g_index <- gen_index(participation_probs$trial, participation_probs$nottrial) %>% round(4)
  }


  cat(paste0("The generalizability index of the sample on the selected covariates is ", g_index, ".\n\n"))

  cat("Covariate Distributions: \n \n")
  cov_tab <- covariate_table(trial, selection_covariates, data)

  print(cov_tab)

  n_trial <- data %>%
    filter(trial == 1) %>%
    nrow()

  n_pop <- data %>%
    filter(trial == 0) %>%
    nrow()

  data_output <- data %>%
    select(trial, selection_covariates)

  out <- list(
    g_index = g_index,
    selection_method = selection_method,
    selection_covariates = selection_covariates,
    trial_name = trial,
    n_trial = n_trial,
    n_pop = n_pop,
    trim_pop = trim_pop,
    n_excluded = n_excluded,
    participation_probs = participation_probs,
    covariate_table = cov_tab,
    data = data_output
  )

  class(out) <- "generalizer_assess"

  return(invisible(out))
}

print.generalize_assess <- function(x,...) {
  cat("A generalizer_assess object: \n")
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
  selection_method_name = c("Logistic Regression", "Random Forest", "Lasso")
  selection_method = c("lr", "rf", "lasso")
  prob_dist_table = rbind(summary(object$participation_probs$trial),
                          summary(object$participation_probs$population))
  row.names(prob_dist_table) = paste0(c("Trial","Population"), " (n = ", c(object$n_trial, object$n_pop),")")

  selection_formula = paste0(object$trial_name," ~ ", paste(object$selection_covariates, collapse = " + "))

  out = list(
    selection_formula = selection_formula,
    selection_method = selection_method_name[selection_method == object$selection_method],
    g_index = object$g_index,
    prob_dist_table = prob_dist_table,
    covariate_table = round(object$covariate_table, 4),
    trim_pop = object$trim_pop,
    n_excluded = object$n_excluded
  )

  class(out) <- "summary.generalize_assess"
  return(out)
}

print.summary.generalize_assess <- function(x,...){
  cat("Probability of Trial Participation: \n \n")
  cat(paste0("Selection Model: ", x$selection_formula," \n \n"))
  print(x$prob_dist_table)
  cat("\n")
  cat(paste0("Estimated by ", x$selection_method, "\n"))
  cat(paste0("Generalizability Index: ", g_index, "\n"))
  cat("============================================ \n")
  if(x$trim_pop == TRUE){
    cat("Population data were trimmed for covariates to not exceed trial covariate bounds \n")
    cat(paste0("Number excluded from population: ", x$n_excluded , "\n \n"))
  }
  cat("Covariate Distributions: \n \n")
  print(round(x$covariate_table, 4))
  invisible(x)
}
