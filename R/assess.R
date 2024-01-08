#' Assess Generalizability of Randomized Sample to Population
#'
#' This function, given a stacked data frame containing both sample and population data, assesses the generalizability of the sample to the population on given covariates.
#'
#' @param sample_indicator variable name denoting sample membership (1 = in sample, 0 = out of sample)
#' @param guided logical. Default is TRUE. If FALSE, then user must enter all arguments in function to bypass guided mode
#' @param covariates vector of covariate names in data set that predict sample membership
#' @param data data frame comprised of "stacked" sample and target population data
#' @param estimation_method method to estimate the probability of sample membership (propensity scores). Default is logistic regression ("lr"). Other methods supported are Random Forests ("rf") and Lasso ("lasso")
#' @param disjoint_data logical. If TRUE, then sample and population data are considered disjoint. This affects calculation of the weights - see details for more information.
#' @param trim_pop logical. If TRUE, then population data are subset to exclude individuals with covariates outside bounds of sample covariates.
#' @export
#' @importFrom rlang is_empty
#' @importFrom dplyr pull

assess <- function(data,
                   guided = TRUE,
                   sample_indicator,
                   covariates,
                   estimation_method = "lr",
                   disjoint_data = TRUE,
                   trim_pop = FALSE) {

  # Check whether data object is of type 'data.frame'
  assertthat::on_failure(is.data.frame) <- function(call, env) {
    "You must pass an object of type 'data.frame' to the 'data' argument."
  }

  assertthat::assert_that(is.data.frame(data))

  data_name <<- data %>%
    lazyeval::expr_text()

  # GUIDED VERSION

  if (guided) {
    user_choices <- .assess.guided(data)

    sample_indicator <- user_choices$sample_indicator

    covariates <- user_choices$covariates

    disjoint_data <- user_choices$disjoint_data

    trim_pop <- user_choices$trim_pop

    estimation_method <- user_choices$estimation_method
  }

  # NON-GUIDED VERSION
  else {
    estimation_method <- tolower(estimation_method)

    invalid_covariates <- covariates %>% setdiff(names(data))

    # Ensure selection covariates are variables in the dataframe provided
    assertthat::on_failure(is_empty) <- function(call, env) {
      paste(
        "The following covariates are not variables in the data provided:\n",
        paste(crayon::blue$bold(invalid_covariates),
          collapse = ", "
        )
      )
    }

    assertthat::assert_that(is_empty(invalid_covariates))

    # Ensure sample variable is binary
    is_sample_indicator_binary <- function(sample_indicator) {
      all(dplyr::pull(data, tidyselect::all_of(sample_indicator)) %in% c(0, 1))
    }

    assertthat::on_failure(is_sample_indicator_binary) <- function(call, env) {
      "Sample membership variable must be coded as `0` (out of sample) or `1` (in sample)."
    }

    assertthat::assert_that(is_sample_indicator_binary(sample_indicator))

    # Ensure estimation method is valid
    is_estimation_method_valid <- function(estimation_method) {
      estimation_method %in% c("lr", "rf", "lasso")
    }

    assertthat::on_failure(is_estimation_method_valid) <- function(call, env) {
      "Invalid estimation method. Please choose one of 'lr' (logistic regression), 'rf' (random forest), or 'lasso' as your estimation method."
    }

    assertthat::assert_that(is_estimation_method_valid(estimation_method))
  }

  # Keep only sample_indicator and covariate columns in data and drop missing values
  data <- data %>%
    dplyr::select(tidyselect::all_of(sample_indicator), tidyselect::all_of(covariates)) %>%
    tidyr::drop_na()

  if (!trim_pop) {
    n_excluded <- 0
  } else {
    trim_pop_output <- .trim.pop(data, sample_indicator, covariates)

    n_excluded <- trim_pop_output$n_excluded

    data <- trim_pop_output$trimmed_data
  }

  # Generate Participation Probabilities
  ps <- .generate.ps(
    data,
    sample_indicator,
    covariates,
    estimation_method
  )

  propensity_scores <- list(
    not_in_sample = ps[which(data[, sample_indicator] == 0)],
    in_sample = ps[which(data[, sample_indicator] == 1)],
    population = ps
  )

  # Calculate Generalizability Index and sample and population sizes

  n_sample <- data %>%
    dplyr::filter(!!rlang::sym(sample_indicator) == 1) %>%
    nrow()

  # If data is not disjoint, compare in_sample to (not_in_sample + in_sample)
  if (!disjoint_data) {
    gen_index <- .get.gen.index(propensity_scores$in_sample, propensity_scores$population) %>%
      round(4)

    n_pop <- data %>%
      nrow()
  }

  # If data is disjoint, compare in_sample to not_in_sample
  else {
    gen_index <- .get.gen.index(propensity_scores$in_sample, propensity_scores$not_in_sample) %>% round(4)

    n_pop <- data %>%
      dplyr::filter(!!rlang::sym(sample_indicator) == 0) %>%
      nrow()
  }

  cat(paste0("The generalizability index of the sample to the target population based on the \nselected covariates is ", crayon::cyan$bold(gen_index), ".\n\n"))

  cov_tab_out <- .make.covariate.table(data,
    sample_indicator = sample_indicator,
    covariates = covariates,
    sample_weights = NULL,
    estimation_method = estimation_method,
    disjoint_data = disjoint_data
  )

  data_output <- data %>%
    dplyr::select(tidyselect::all_of(sample_indicator), tidyselect::all_of(covariates))

  out <- list(
    gen_index = gen_index,
    estimation_method = estimation_method,
    disjoint_data = disjoint_data,
    sample_indicator = sample_indicator,
    n_sample = n_sample,
    n_pop = n_pop,
    trim_pop = trim_pop,
    n_excluded = n_excluded,
    propensity_scores = propensity_scores,
    covariate_table = cov_tab_out$covariate_table,
    covariate_kable = cov_tab_out$covariate_kable,
    cov_dist_facet_plot = cov_tab_out$cov_dist_facet_plot,
    cov_dist_plots = cov_tab_out$cov_dist_plots,
    data = data_output,
    data_name = data_name
  )

  rm(
    list = deparse(substitute(data_name)),
    envir = .GlobalEnv
  ) # delete object data_name from global environment

  class(out) <- "generalizeR_assess"

  return(invisible(out))
}

.assess.guided <- function(data) {
  var_names <- names(data)

  cat(crayon::bold("\nWelcome to assess()! \n\n"))

  cat("Given a data frame containing both sample and population data, this function will \nassess the generalizability of your sample to the population based on your selected \ncovariates.\n\n")

  cat("You may exit out of this function at any time by pressing <Esc>.\n\n")

  sample_indicator <- .select.sample.indicator(data)

  covariates <- data %>%
    dplyr::select(-tidyselect::all_of(sample_indicator)) %>%
    .select.covariates()

  disjoint_data <- .yes.no("Are the sample data and population data disjoint? See the vignette for more details.\n")

  trim_pop <- .yes.no("Should the population data be trimmed to exclude units with covariate values outside \nthe bounds of the sample covariates? See the vignette for more details.\n")

  estimation_method <- .select.method()

  output <- list(
    sample_indicator = sample_indicator,
    covariates = covariates,
    disjoint_data = disjoint_data,
    trim_pop = trim_pop,
    estimation_method = estimation_method
  )

  return(invisible(output))
}

.select.sample.indicator <- function(data) {
  choices <- names(data)

  num_choices <- length(choices)

  options <- paste0(
    format(1:num_choices),
    ":  ",
    choices
  )

  repeat {
    cat(paste0(
      "Here is a list of the variables in the '",
      data_name,
      "' dataframe.\n"
    ))

    if (num_choices > 10L) {
      formatted_options <- format(options)
      nw <- nchar(formatted_options[1L], "w") + 2L
      ncol <- getOption("width") %/% nw

      if (ncol > 1L) {
        options <- paste0(formatted_options,
          c(rep.int("  ", ncol - 1L), "\n"),
          collapse = ""
        )
      }

      cat("", options, sep = "\n")
    } else {
      cat("", options, "", sep = "\n")
    }

    cat("Please type the number corresponding to the binary variable encoding sample \nmembership in your dataframe and then hit <Return> to continue.\n\n")

    selection <- tryCatch(

      scan("",
        what = 0,
        quiet = TRUE,
        nlines = 1
      ),
      error = identity
    )

    # Verify that user's selection did not throw an error in tryCatch()
    if (!inherits(selection, "error")) {
      selection <- selection %>%
        unique()

      # Verify that user chose exactly one variable
      if (length(selection) == 1L) {

        # Verify that user only input an integer between 1 and the number of variables in the dataframe
        if (selection %in% 1:num_choices) {

          # Verify user chose a binary variable
          if (all(dplyr::pull(data, choices[selection]) %in% c(0, 1))) {
            return(choices[selection])
          }

          cat(crayon::red("\nERROR: Invalid selection. The variable denoting sample membership must be coded as `0` \n(out of sample) or `1` (in sample).\n\n"))
          next
        }

        cat(paste0(
          crayon::red("\nERROR: Invalid selection. Your input must be a single integer between 1 and "),
          crayon::red(num_choices),
          crayon::red(".\n\n")
        ))
        next
      }

      cat(crayon::red("\nERROR: Invalid selection. You must select exactly one variable.\n\n"))
      next
    }

    cat(crayon::red("\nERROR: Invalid selection. Please try again.\n\n"))
  }
}

.select.covariates <- function(data) {
  choices <- names(data)

  num_choices <- length(choices)

  options <- paste0(
    format(1:num_choices),
    ":  ",
    choices
  )

  cat("\n")

  repeat {
    cat(paste0(
      "Here are the remaining variables in the '",
      data_name,
      "' dataframe.\n"
    ))

    if (num_choices > 10L) {
      formatted_options <- format(options)
      nw <- nchar(formatted_options[1L], "w") + 2L
      ncol <- getOption("width") %/% nw

      if (ncol > 1L) {
        options <- paste0(formatted_options,
          c(rep.int("  ", ncol - 1L), "\n"),
          collapse = ""
        )
      }

      cat("", options, sep = "\n")
    } else {
      cat("", options, "", sep = "\n")
    }

    cat("Please select the covariates that will be used to predict sample membership.\n\n")

    selection <- tryCatch(

      scan("",
        what = 0,
        quiet = TRUE,
        nlines = 1
      ),
      error = identity
    )

    # Verify that user's selection did not throw an error in tryCatch()
    if (!inherits(selection, "error")) {
      selection <- selection %>%
        unique()

      # Verify that user chose at least one variable
      if (length(selection) >= 1L) {

        # Verify that user only input integers between 1 and the number of variables in the dataframe
        if (all(selection %in% 1:num_choices)) {
          return(choices[selection])
        }

        cat(paste0(
          crayon::red("\nERROR: Invalid selection. Each input must be a single integer between 1 and "),
          crayon::red(num_choices),
          crayon::red(".\n\n")
        ))
        next
      }

      cat(crayon::red("\nERROR: Invalid selection. You must select at least one variable.\n\n"))
      next
    }

    cat(crayon::red("\nERROR: Invalid selection. Please try again.\n\n"))
  }
}

.yes.no <- function(message = "") {
  options <- paste0(
    format(1:2),
    ":  ",
    c("Yes", "No")
  )

  cat("\n")

  repeat {
    cat(message)

    cat("", options, "", sep = "\n")

    selection <- tryCatch(

      scan("",
        what = 0,
        quiet = TRUE,
        nlines = 1
      ),
      error = identity
    )

    # Verify that user's selection did not throw an error in tryCatch()
    if (!inherits(selection, "error")) {
      selection <- selection %>%
        unique()

      # Verify that user input a single number
      if (length(selection) == 1) {

        # Verify that user only input either a 1 or a 2
        if (selection %in% 1:2) {
          disjoint_data <- switch(selection,
            "1" = TRUE,
            "2" = FALSE
          )

          return(disjoint_data)
        }
      }
    }

    cat(crayon::red("\nERROR: Invalid selection. You may only input either a 1 or a 2.\n\n"))
  }
}

.select.method <- function() {
  choices <- c("Logistic Regression", "Random Forest", "Lasso")

  options <- paste0(
    format(1:3),
    ":  ",
    c("Logistic Regression", "Random Forest", "Lasso")
  )

  cat("\n")

  repeat {
    cat("Please select the method that will be used to estimate the probability of sample membership (the default is Logistic Regression).\n")

    cat("", options, "", sep = "\n")

    selection <- tryCatch(

      scan("",
        what = 0,
        quiet = TRUE,
        nlines = 1
      ),
      error = identity
    )

    # Verify that user's selection did not throw an error in tryCatch()
    if (!inherits(selection, "error")) {
      selection <- selection %>%
        unique()

      # Verify that user chose exactly one method
      if (length(selection) == 1L) {

        # Verify that user only input an integer between 1 and the number of variables in the dataframe
        if (selection %in% 1:3) {
          estimation_method <- switch(choices[selection],
            "Logistic Regression" = "lr",
            "Random Forest" = "rf",
            "Lasso" = "lasso"
          )

          return(estimation_method)
        }

        cat(crayon::red("\nERROR: Invalid selection. Your input must be a single integer between 1 and 3.\n\n"))
        next
      }

      cat(crayon::red("\nERROR: Invalid selection. You must select exactly one method.\n\n"))
      next
    }

    cat(crayon::red("\nERROR: Invalid selection. Please try again.\n\n"))
  }
}


#' Find covariate bounds in the sample data
#'
#' @param covariate covariate in data set that predicts sample membership
#' @param sample_indicator variable denoting sample membership (1 = in sample, 0 = out of sample)
#' @param data data frame comprised of "stacked" sample and target population data
#' @return \code{covariate_bounds} Returns a dataframe

.get.covariate.bounds <- function(covariate,
                                  sample_indicator,
                                  data) {

  # Make covariate vector but only for observations selected to be part of the sample
  sample_covariate <- data %>%
    dplyr::filter(!!rlang::sym(sample_indicator) == 1) %>%
    dplyr::pull(covariate)

  return(sample_covariate)

  if (sample_covariate %>% is.factor()) {
    sample_covariate_levels <- sample_covariate %>%
      droplevels() %>%
      levels()

    return(
      data %>%
        dplyr::mutate(test = !(!!covariate %in% sample_covariate_levels)) %>% # The !!-operator (bang-bang) evaluates covariate first so the expression it contains is what gets passed to mutate()
        dplyr::pull(test) %>%
        which()
    )
  }

  if (sample_covariate %>% is.numeric()) {
    sample_bounds <- c(sample_covariate %>% min(na.rm = TRUE), sample_covariate %>% max(na.rm = TRUE))

    return(
      data %>%
        dplyr::mutate(test = !(!!covariate %>% dplyr::between(sample_bounds[1], sample_bounds[2]))) %>%
        dplyr::pull(test) %>%
        which()
    )
  }
}

#' Subset Population so Population Covariates are within bounds of Sample Covariates
#'
#' @param data data frame comprised of "stacked" sample and target population data
#' @param sample_indicator variable name denoting sample membership (1 = in sample, 0 = out of sample)
#' @param covariates vector of covariate names in data set that predict sample membership
#' @return \code{trim_pop} returns a data frame, where the target population covariates do not exceed the bounds of the sample covariates

.trim.pop <- function(data,
                      sample_indicator,
                      covariates) {

  ##### CHECKS #####
  if (!is.data.frame(data)) {
    stop("Data must be of type 'data.frame'.")
  }

  invalid_covariates <- covariates %>% setdiff(names(data))

  if (!is_empty(invalid_covariates)) {
    stop(paste("The following covariates are not variables in the data provided:\n", paste(blue$bold(invalid_covariates), collapse = ", ")))
  }

  sample_indicator_valid <- data %>%
    dplyr::pull(tidyselect::all_of(sample_indicator)) %>%
    na.omit() %>%
    table() %>%
    names() %>%
    setequal(c("0", "1"))

  if (!sample_indicator_valid) {
    stop("Sample membership variable must be binary and coded as `0` (out of sample) or `1` (in sample)")
  }

  ##### Find and remove rows of population data that violate bounds #####
  bound_violations <- purrr::map(covariates, .get.covariate.bounds, sample_indicator, data)

  missing_rows <- bound_violations %>%
    unlist() %>%
    unique()

  trimmed_data <- data %>%
    dplyr::filter(!row.names(data) %in% missing_rows) %>%
    droplevels() # Get rid of unused levels from factors

  ##### Number of rows in population data excluded #####
  n_excluded <- missing_rows %>%
    length()

  out <- list(
    n_excluded = n_excluded,
    trimmed_data = trimmed_data,
    untrimmed_data = data
  )

  return(invisible(out))
}

#' Calculate Generalizability Index
#'
#' This function is easiest to use through 'assess()' but can also be used independently.
#'
#' It calculates the generalizability index, a value between 0 and 1, that represents how generalizable a given sample is to a given population on specified covariates. For more information on calculation and interpretation, please see Tipton (2014).
#'
#' @param sample_ps vector of probabilities of sample membership among individuals in the sample
#' @param pop_ps vector of probabilities of sample membership among individuals in the population
#' @return the generalizability index, a value between 0 and 1, where a higher score indicates greater similarity
#' @importFrom stats dnorm integrate
#' @references
#' Tipton, E. (2014). How generalizable is your experiment? An index for comparing experimental samples and populations. *Journal of Educational and Behavioral Statistics*, *39*(6), 478-501.
#' @md

.get.gen.index <- function(sample_ps,
                           pop_ps) {
  ## Baklizi and Eidous (2006) estimator
  # bandwidth

  if (var(sample_ps) == 0 & var(pop_ps) == 0) {
    return(1)
  } else {
    h <- function(x) {
      n <- length(x)
      optim_binwidth <- (4 * sqrt(var(x))^5 / (3 * n))^(1 / 5)

      if (is.na(optim_binwidth) | is.nan(optim_binwidth)) {
        optim_binwidth <- 0
      }
      if (optim_binwidth < 0.001) { # this yielded a b index of 0.9999501 for (at least one specific case of) "perfectly" stratified data

        optim_binwidth <- 0.001
      }

      return(optim_binwidth)
    }

    # Kernel estimators of the density and the distribution
    kg <- function(x, data) {
      hb <- h(data) # Bin width
      k <- r <- length(x)
      for (i in 1:k) r[i] <- mean(dnorm((x[i] - data) / hb)) / hb # We divide by bin width, which is a problem when bin width goes to zero
      return(r)
    }

    return(as.numeric(integrate(function(x) sqrt(kg(x, sample_ps) * kg(x, pop_ps)), -Inf, Inf)$value))
  }
}

#' Print method for "generalizeR_assess" class
#'
#' @param x An object of class "generalizeR_assess"
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.generalizeR_assess
#' @export

print.generalizeR_assess <- function(x, ...) {
  cat("\nA generalizeR_assess object: \n\n")

  cat(" - Dataset name:", crayon::cyan$bold(x$data_name), "\n\n")

  covariate_names <- x$covariate_table %>%
    pull(covariate) %>%
    crayon::cyan$bold() %>%
    paste(collapse = ", ") %>%
    gsub("(.{200})\\s(,*)", "\\1\n   \\2", .)

  cat(" - Covariates selected:\n\n  ", covariate_names, "\n\n")

  cat(
    " - Method used to estimate probability of sample membership:",
    switch(x$estimation_method,
      "lr" = "Logistic Regression",
      "rf" = "Random Forest",
      "lasso" = "LASSO"
    ) %>%
      crayon::cyan$bold(),
    "\n\n"
  )

  if (x$disjoint_data) {
    cat(paste0(" - The sample and the population were considered ", crayon::cyan$bold("disjoint"), " from one another.\n\n"))
  } else {
    cat(paste0(" - The sample was considered a ", crayon::cyan$bold("subset"), " of the population.\n\n"))
  }

  cat(" - Sample size:", crayon::cyan$bold(x$n_sample), "\n\n")

  cat(" - Population size:", crayon::cyan$bold(x$n_pop), "\n\n")

  if (x$n_excluded > 0) {
    cat(" - Number of observations trimmed from population: ", crayon::cyan$bold(x$n_excluded), "\n\n")
  }

  cat(paste0(" - The generalizability index of the sample to the target population based on the selected covariates is ", crayon::cyan$bold(x$gen_index), "."))
}

#' Summary method for "generalizeR_assess" class
#'
#' @param object An object of class "generalizeR_assess"
#' @param ... Other arguments passed to or from other methods
#'
#' @export summary.generalizeR_assess
#' @export

summary.generalizeR_assess <- function(object, ...) {
  estimation_method <- switch(object$estimation_method,
    "lr" = "Logistic Regression",
    "rf" = "Random Forest",
    "lasso" = "Lasso"
  )

  if (object$disjoint_data) {
    prop_score_dist_table <- rbind(
      summary(object$propensity_scores$in_sample),
      summary(object$propensity_scores$population)
    )

    row.names(prop_score_dist_table) <- paste0(c("Sample", "Population"), " (n = ", c(object$n_sample, object$n_pop), ")")

    sample_prop_scores <- data.frame(prop_scores = object$propensity_scores$in_sample) %>%
      mutate(sample_indicator = 1)

    pop_prop_scores <- data.frame(prop_scores = object$propensity_scores$population) %>%
      mutate(sample_indicator = 0)

    prop_scores <- rbind(sample_prop_scores, pop_prop_scores)

    prop_score_dist_plot <- prop_scores %>%
      ggplot2::ggplot() +
      geom_density(aes(x = gtools::logit(prop_scores), fill = factor(sample_indicator)),
        alpha = 0.7
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_discrete(
        name = NULL,
        labels = c("Population", "Sample")
      ) +
      labs(
        x = "Propensity Score Logits",
        y = "Density",
        title = "Distribution of Propensity Score Logits"
      ) +
      theme_minimal() +
      theme(
        axis.ticks.x = element_line(),
        axis.text.y = element_blank(),
        axis.line = element_line(),
        plot.title = element_text(size = 12)
      )
  } else {
    prop_score_dist_table <- rbind(
      summary(object$propensity_scores$in_sample),
      summary(object$propensity_scores$not_in_sample)
    )

    row.names(prop_score_dist_table) <- paste0(c("In Sample", "Not In Sample"), " (n = ", c(object$n_sample, object$n_pop), ")")

    in_sample_prop_scores <- data.frame(prop_scores = object$propensity_scores$in_sample) %>%
      mutate(sample_indicator = 1)

    not_in_sample_prop_scores <- data.frame(prop_scores = object$propensity_scores$not_in_sample) %>%
      mutate(sample_indicator = 0)

    prop_scores <- rbind(in_sample_prop_scores, not_in_sample_prop_scores)

    prop_score_dist_plot <- prop_scores %>%
      ggplot2::ggplot() +
      geom_density(aes(x = gtools::logit(prop_scores), fill = factor(sample_indicator)),
        alpha = 0.7
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_discrete(
        name = NULL,
        labels = c("Not in Sample", "In Sample")
      ) +
      labs(
        x = "Propensity Score Logits",
        y = "Density",
        title = "Distribution of Propensity Score Logits"
      ) +
      theme_minimal() +
      theme(
        axis.ticks.x = element_line(),
        axis.text.y = element_blank(),
        axis.line = element_line(),
        axis.title = element_blank(),
        plot.title = element_text(size = 12)
      )
  }

  colnames(prop_score_dist_table) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")

  prop_score_dist_table <- prop_score_dist_table %>%
    data.frame() %>%
    colorDF::colorDF(theme = "dark")

  out <- list(
    estimation_method = estimation_method,
    gen_index = object$gen_index,
    prop_score_dist_table = prop_score_dist_table,
    prop_score_dist_plot = prop_score_dist_plot,
    covariate_table = object$covariate_table,
    n_excluded = object$n_excluded
  )

  class(out) <- "summary.generalizeR_assess"

  return(out)
}

#' Print method for "summary.generalizeR_assess" class
#'
#' @param x An object of class "summary.generalizeR_assess"
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.summary.generalizeR_assess
#' @export



print.summary.generalizeR_assess <- function(x, ...) {
  cat("\nSummary of Estimated Propensity Scores: \n\n")

  print(x$prop_score_dist_table)

  print(x$prop_score_dist_plot)

  cat("\n\nEstimation Method:", crayon::cyan$bold(x$estimation_method), "\n\n")

  covariate_names <- x$covariate_table %>%
    dplyr::pull(covariate) %>%
    crayon::cyan$bold() %>%
    paste(collapse = ", ") %>%
    gsub("(.{200})\\s(,*)", "\\1\n   \\2", .)

  cat("Covariates Used:\n\n  ", covariate_names, "\n\n")

  cat("Generalizability Index: ", crayon::cyan$bold(x$gen_index), "\n\n")

  if (x$n_excluded > 0) {
    cat("The dataset was trimmed to ensure population covariates do not exceed sample covariate bounds.\n\n")

    cat("Number of observations trimmed from population: ", crayon::cyan$bold(x$n_excluded), "\n\n")
  }

  cat("Covariate Table: \n\n")

  print(x$covariate_table)
}

# if(getRversion() >= "2.15.1") utils::globalVariables(c("test", "data_name", "covariate", ".", "object", "geom_density",
#                                                          "sample_indicator", "scale_x_continuous", "scale_y_continuous",
#                                                          "scale_fill_discrete", "theme_minimal", "element_line"))
