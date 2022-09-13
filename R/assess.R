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
#' @export


assess <- function(data,
                   guided = TRUE,
                   trial,
                   selection_covariates,
                   selection_method = "lr",
                   is_data_disjoint = TRUE,
                   trim_pop = FALSE) {

  # Check whether data object is of type 'data.frame'
  assertthat::on_failure(is.data.frame) <- function(call, env) {
    "You must pass an object of type 'data.frame' to the 'data' argument."
  }

  assertthat::assert_that(is.data.frame(data))

  ##### GUIDED VERSION #####

  if(guided == TRUE) {

    user_choices <- .assess.guided(data)

    trial <- user_choices$trial

    selection_covariates <- user_choices$selection_covariates

    selection_method <- user_choices$selection_method

    is_data_disjoint <- user_choices$is_data_disjoint

    trim_pop <- user_choices$trim_pop
  }

  ##### NON-GUIDED VERSION #####
  else {

    selection_method <- tolower(selection_method)

    invalid_selection_covariates <- selection_covariates %>% setdiff(names(data))

    ##### Ensure selection covariates are variables in the dataframe provided #####
    assertthat::on_failure(is_empty) <- function(call, env) {

      paste("The following covariates are not variables in the data provided:\n",
            paste(crayon::blue$bold(invalid_selection_covariates),
                  collapse = ", ")
            )
    }

    assertthat::assert_that(is_empty(invalid_selection_covariates))

    ##### Ensure trial variable is binary #####
    is_trial_binary <- function(trial) {

      !anyNA(match(names(table(data[,trial])), c("0","1")))
    }

    assertthat::on_failure(is_trial_binary) <- function(call, env) {

      "Sample membership variable must be coded as `0` (not in trial) or `1` (in trial)."
    }

    assertthat::assert_that(is_trial_binary(trial))

    ##### Ensure selection method is valid #####
    is_selection_method_valid <- function(selection_method){

      selection_method %in% c("lr","rf","lasso")
    }

    assertthat::on_failure(is_selection_method_valid) <- function(call, env) {

      "Invalid selection method. Please choose one of 'lr' (logistic regression), 'rf' (random forest), or 'lasso' as your selection method."
    }

    assertthat::assert_that(is_selection_method_valid(selection_method))
  }

  ##### Drop missing values from trial and selection covariate columns in data #####
  data <- data %>%
    tidyr::drop_na(trial, tidyselect::all_of(selection_covariates))

  if(trim_pop == FALSE) {

    n_excluded <- 0
  }

  else {

    trim_pop_output <- trim_pop(data, trial, selection_covariates)

    n_excluded <- trim_pop_output$n_excluded

    data <- trim_pop_output$trimmed_data
  }

  ##### Generate Participation Probabilities #####
  ps <- generate_ps(data,
                    trial,
                    selection_covariates,
                    selection_method)

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

  cat(paste0("\nThe generalizability index of the sample on the selected covariates is ", g_index, ".\n\n"))

  cat(crayon::blue$bold("Covariate Distributions:\n"))

  cov_tab <- covariate_table(data,
                             trial,
                             selection_covariates)

  print(cov_tab)

  n_trial <- data %>%
    dplyr::filter(trial == 1) %>%
    nrow()

  n_pop <- data %>%
    dplyr::filter(trial == 0) %>%
    nrow()

  data_output <- data %>%
    dplyr::select(trial, selection_covariates)

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

.assess.guided <- function(data) {

  var_names <- names(data)

  cat(crayon::bold("\nWelcome to assess()! \n\n"))

  cat("Given a data frame containing both sample and population data, this function will \nassess the generalizability of your sample to the population based on your selected \ncovariates.\n\n")

  repeat {

    cat("Please select the variable coding binary trial participation in your dataframe: \n")
    trial <- .select.list(choices = var_names,
                            graphics = FALSE,
                            multiple = FALSE)

    if(!anyNA(match(names(table(data[,trial])), c("0","1")))) {

      break
    }

    cat(crayon::red("\nThe trial variable must be coded as '0' (not in trial) or '1' (in trial).\n\n"))
  }

  cat("Please select the covariates that will be used to predict trial participation: \n")
  selection_covariates <- utils::select.list(choices = var_names %>% setdiff(trial),
                                             graphics = FALSE,
                                             multiple = TRUE)

  cat("Are the trial data and population data disjoint? See the vignette for more details.\n")
  is_data_disjoint <- utils::menu(choices = c("Yes", "No"),
                                  graphics = FALSE) %>%
    switch("1" = TRUE,
           "2" = FALSE)

  cat("Should the population data be trimmed to exclude units with covariate values outside \nthe bounds of the trial covariates? See the vignette for more details. \n")
  trim_pop <- utils::menu(choices = c("Yes", "No"),
                          graphics = FALSE) %>%
    switch("1" = TRUE,
           "2" = FALSE)

  cat("Please select the method that will be used to estimate the probability of trial participation: \n")
  selection_method <- utils::select.list(choices = c("Logistic Regression", "Random Forest", "Lasso"),
                                         graphics = FALSE,
                                         multiple = FALSE) %>%
    switch("Logistic Regression" = "lr",
           "Random Forest" = "rf",
           "Lasso" = "lasso")

  output <- list(trial = trial,
                 selection_covariates = selection_covariates,
                 is_data_disjoint = is_data_disjoint,
                 trim_pop = trim_pop,
                 selection_method = selection_method)

  return(invisible(output))
}

#' Subset Population so Population Covariates are within bounds of Trial Covariates
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @return \code{trim_pop} returns a data frame, where the target population covariates do not exceed the bounds of the trial covariates

trim_pop <- function(data,
                     trial,
                     selection_covariates) {

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

generate_ps <- function(data,
                        trial,
                        selection_covariates,
                        selection_method) {

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
                           data = data %>% select(tidyselect::all_of(selection_covariates)))

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


#' Calculate Generalizability Index
#'
#' This function is easiest to use through 'assess()' but can also be used independently.
#'
#' It calculates the generalizability index, a value between 0 and 1, that represents how generalizable a given sample is to a given population on specified covariates. For more information on calculation and interpretation, please see Tipton (2014).
#'
#' @param trial_ps vector of probabilities of sample participation among individuals in the trial
#' @param pop_ps vector of probabilities of sample participation among individuals in the population
#' @return the generalizability index, a value between 0 and 1, where a higher score indicates greater similarity
#' @export
#' @importFrom stats dnorm integrate
#' @references
#' Tipton, E. (2014). How generalizable is your experiment? An index for comparing experimental samples and populations. *Journal of Educational and Behavioral Statistics*, *39*(6), 478-501.
#' @md

gen_index <- function(trial_ps,
                      pop_ps) {
  ##Baklizi and Eidous (2006) estimator
  # bandwidth

  if(var(trial_ps) == 0 & var(pop_ps) == 0) {return(1)}

  else {

    h = function(x){

      n = length(x)
      optim_binwidth = (4*sqrt(var(x))^5/(3*n))^(1/5)

      if(is.na(optim_binwidth) | is.nan(optim_binwidth)){

        optim_binwidth = 0
      }
      if(optim_binwidth < 0.001) { # this yielded a b index of 0.9999501 for (at least one specific case of) "perfectly" stratified data

        optim_binwidth = 0.001
      }

      return(optim_binwidth)
    }

    # kernel estimators of the density and the distribution
    kg = function(x, data){

      hb = h(data) #bin width
      k = r = length(x)
      for(i in 1:k) r[i] = mean(dnorm((x[i] - data)/hb))/hb # we divide by bin width, which is a problem when bin width goes to zero
      return(r)
    }

    return(as.numeric(integrate(function(x) sqrt(kg(x, trial_ps)*kg(x, pop_ps)), -Inf, Inf)$value))
  }
}


#' Create Covariate Balance Table
#'
#' This function is designed for use within 'assess().'
#'
#' @param trial variable name denoting binary trial participation (1 = trial participant, 0 = not trial participant)
#' @param selection_covariates vector of covariate names in data set that predict trial participation
#' @param data data frame comprised of "stacked" trial and target population data
#' @param weighted_table defaults to FALSE; whether weights are already included and do not need to be estimated
#' @param selection_method method to estimate the probability of trial participation.  Default is logistic regression ("lr").  Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param is_data_disjoint defaults to TRUE. If TRUE, then trial and population data are considered independent.  This affects calculation of the weights
#' @export
#' @importFrom stats model.matrix weighted.mean
#' @importFrom dplyr funs

covariate_table <- function(data,
                            trial,
                            selection_covariates,
                            weighted_table = FALSE,
                            selection_method = "lr",
                            is_data_disjoint = TRUE){

  V1 <- V2 <- weights <- population <- pooled_sd <- ASMD <- . <- NULL

  if(weighted_table == FALSE){
    data = data %>%
      tidyr::drop_na(selection_covariates) %>%
      as.data.frame()

    expanded.data = data.frame(trial = data[, trial],
                               model.matrix(~-1 + ., data = data[, selection_covariates]))

    means.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[-1], mean) %>%
      t() %>% as.data.frame()

    means.tab = means.tab[-1,]

    names(means.tab) = c("trial", "population")
    n_trial = as.numeric(table(expanded.data[,"trial"]))[2]
    n_pop = as.numeric(table(expanded.data[,"trial"]))[1]
    sd.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_all(var) %>% t() %>% as.data.frame() %>% .[-1,] %>%
      mutate(pooled_sd = sqrt(((n_trial - 1) * V1 + (n_pop - 1) * V2)/(n_trial + n_pop - 2)))
    names(sd.tab) = c("trial_var", "population_var", "pooled_sd")
  }

  if(weighted_table == TRUE){
    data = data %>%
      tidyr::drop_na(selection_covariates) %>%
      as.data.frame()

    data$weights = weighting(outcome = NULL, treatment = NULL, trial = trial,
                             selection_covariates = selection_covariates, data = data,
                             selection_method = selection_method, is_data_disjoint = is_data_disjoint)$weights
    data$weights = ifelse(data[,trial] == 0, 1, data$weights)

    expanded.data = data.frame(trial = data[,trial], model.matrix(~ -1 + ., data = data[,c(selection_covariates,"weights")]))

    means.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[-1], funs(weighted.mean(., weights))) %>%
      dplyr::select(-`weights`) %>%
      t() %>%
      as.data.frame() %>% .[-1,]

    names(means.tab) = c("trial","population")

    n_trial = as.numeric(table(expanded.data$trial))[2]
    n_pop = as.numeric(table(expanded.data$trial))[1]

    sd.tab = expanded.data %>%
      dplyr::group_by(trial) %>%
      dplyr::summarise_at(names(expanded.data)[-1],
                          funs(sum(weights * (. - weighted.mean(.,weights))^2)/sum(weights))) %>%
      dplyr::select(-`weights`) %>%
      t() %>%
      as.data.frame() %>% .[-1,] %>%
      mutate(pooled_sd = sqrt(((n_trial - 1)*V1 + (n_pop - 1)*V2)/(n_trial + n_pop - 2)))

    names(sd.tab) = c("trial_var","population_var","pooled_sd")
  }

  covariate_table = means.tab %>%
    dplyr::bind_cols(sd.tab) %>%
    dplyr::mutate(ASMD = round(abs((trial - population)/pooled_sd),3)) %>%
    dplyr::select(trial, population, ASMD)

  if(weighted_table == FALSE){
    row.names(covariate_table) = setdiff(names(expanded.data),c("trial"))
  }

  if(weighted_table == TRUE){
    names(covariate_table)[1] = "trial (weighted)"
    row.names(covariate_table) = setdiff(names(expanded.data),c("trial","weights"))
  }

  return(covariate_table)
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
