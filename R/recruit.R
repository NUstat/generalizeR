#' Recruit Units from a Population for Sampling
#'
#' This function works with the output of 'stratify()'. The user provides the number of units they wish to sample from their population dataset. The function tells the user how many observations to sample from each stratum and generates recruitment lists, one per stratum, which can either be saved to .csv files in any given directory or accessed later on.
#'
#' This function, and the others in this package, are designed to mimic the website https://www.thegeneralizer.org/ based on the papers referenced below.
#'
#' @param stratify_output output from 'stratify()', of S3 class 'generalizeR_stratify'
#' @param guided logical; defaults to TRUE. Whether the function should be guided (ask questions and behave interactively throughout) or not. If set to FALSE, must provide values for other arguments below
#' @param sample_size defaults to NULL. If guided is set to FALSE, must provide a number of units to sample
#' @param save_as_csv defaults to NULL. If guided is set to FALSE, specify whether or not to save recruitment lists to working directory; TRUE or FALSE
#' @return A three-element list containing the recruitment lists and the recruitment table/kable
#' @export
#' @importFrom readr write_csv
#' @importFrom easycsv choose_dir
#' @importFrom utils tail
#' @references
#' Tipton, E. (2014). Stratified sampling using cluster analysis: A sample selection strategy for improved generalizations from experiments. *Evaluation Review*, *37*(2), 109-139.
#'
#' Tipton, E. (2014). How generalizable is your experiment? An index for comparing experimental samples and populations. *Journal of Educational and Behavioral Statistics*, *39*(6), 478-501.
#' @examples
#' library(tidyverse)
#'
#' selection_covariates <- c("total", "pct_black_or_african_american", "pct_white",
#'                           "pct_female", "pct_free_and_reduced_lunch")
#' strat_output <- stratify(generalizeR:::inference_pop, guided = FALSE, n_strata = 4,
#'                          variables = selection_covariates, idvar = "ncessch")
#'
#' recruit(strat_output, guided = FALSE, sample_size = 72, save_as_csv = FALSE)
#'
#' @md

recruit <- function(stratify_output,
                    guided = TRUE,
                    sample_size = NULL,
                    save_as_csv = FALSE) {

  #check that the stratify_output object passed was generated by the generalizeR::stratify() function
  is_generalizeR_stratify <- function(stratify_output) {
    inherits(stratify_output, "generalizeR_stratify")
  }

  assertthat::on_failure(is_generalizeR_stratify) <- function(call, env) {
    "Argument 'stratify_output' must be an object of class 'generalizeR_stratify', created by assigning the output of stratify() to an object."
  }

  assertthat::assert_that(is_generalizeR_stratify(stratify_output))

  pop_size <- stratify_output$pop_data_by_stratum %>% nrow()

  n_strata <- stratify_output$n_strata

  # checking validity of requested sample size w/ the unguided version of the function
  if (!guided) {
    is_valid_sample_size <- function(sample_size) {
      assertthat::is.count(sample_size) && sample_size <= pop_size
    }

    assertthat::on_failure(is_valid_sample_size) <- function(call, env) {
      paste("You must enter a sample size when using the unguided version of this function. The number of units you wish to recruit must be a single integer ranging from 1 to the total number of units in your population (",
        pop_size,
        ").",
        sep = ""
      )
    }

    assertthat::assert_that(is_valid_sample_size(sample_size))
  } else {
    cat(crayon::bold("\nWelcome to recruit()!\n\n"))

    cat("This function generates recruitment lists for each stratum in your population and \ntells you how many units to recruit from each given the total number of units that you \nwish to recruit (your desired sample size) in order to maximize generalizability.\n\n")

    cat("To store your results, make sure you have assigned this function to an object.\n\n")

    readline(prompt = "Press <Return> to continue or <Esc> to exit.")
  }

  cat("\nThe 'generalizeR_stratify' object you've supplied consists of ", paste(crayon::bold(pop_size)), " population units \ndivided into ", paste(crayon::bold(n_strata)), " strata along these variables:\n\n", paste(crayon::blue$bold(stratify_output$variables), collapse = ", "), ".", sep = "")

  recruitment_lists <- .make.recruitment.lists(stratify_output)

  cat("\n\n")
  cat(paste(n_strata), "recruitment lists have been generated, one per stratum.")
  cat(" Each list contains the ID \ninformation for the units, which have been ranked in order of desirability.\n\n")

  if (guided) {
    cat("The recruitment lists have been printed in the Viewer pane to the right.\n\n")

    for (i in 1:n_strata) {
      recruitment_lists[[i]] %>%
        data.frame() %>%
        kableExtra::kbl(
          caption = paste("Recruitment List for Stratum", i),
          align = "c"
        ) %>%
        kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
        print()
    }

    sample_size <- .get.sample.size(1:pop_size)
    cat("\n")
  }

  table_output <- .make.recruitment.table(stratify_output, sample_size)

  cat("The following table (also shown in the Viewer pane to the right) displays the stratum \nsizes, their proportion relative to the total population size, and consequent \nrecruitment number for each stratum. ")
  cat("Ideally, units should be recruited across strata \naccording to these numbers.")
  cat(" Doing so will lead to the least amount of bias and no \nincrease in standard errors. ")
  cat("Note that the recruitment numbers have been rounded to \nintegers in such a way as to ensure their sum equals the desired total sample size.\n\n")

  cat(crayon::blue$bold("Recruitment Table\n"))

  print(table_output$recruit_table, row.names = FALSE)
  print(table_output$recruit_kable)

  cat("\nAttempt to recruit units starting from the top of each recruitment list. If you are \nunsuccessful in recruiting a particular unit, move on to the next one in the list and \ncontinue until you have reached the ideal recruitment number in each stratum.\n\n", sep = "")

  if (guided) {

    repeat{

      save_as_csv <- utils::menu(
        choices = c("Yes", "No"),
        title = cat("Would you like to save the recruitment lists as .csv files?")
      )

      if (!(save_as_csv %in% 1:2)) {

        cat(crayon::red("\nERROR: You must input either '1' for 'Yes' or '2' for 'No'.\n\n"))
        next
      }

      break
    }

    save_as_csv <- save_as_csv %>%
      switch("1" = TRUE,
             "2" = FALSE)

    cat("\n")

    if (save_as_csv) {
      cat("The lists will be saved as 'Recruitment_list_#', one for each stratum, in the current working directory.\n")
      #cat("Where should they be saved?\n\n")
    }
  }

  # for the unguided version of the function, save the recruitment output as a csv if requested
  if (save_as_csv) {
    filepath <- "" #ifelse(guided, easycsv::choose_dir(), "")

    for (i in 1:(n_strata)) {
      filename <- paste(filepath, "Recruitment_list_", i, ".csv", sep = "")
      readr::write_csv(recruitment_lists[[i]], file = filename)
    }

    #if (guided) {
     # cat("Lists saved successfully!\n\n")
    #} else {
     # cat("You've chosen to save your recruitment lists as .csv files. The lists have been saved \nas 'Recruitment_list_#', one for each stratum, inside your current working directory.\n\n")
    #}
  }

  cat("If you have stored the output of 'recruit()' in an object, you can use it to access \nthese lists by typing the name of the object followed by '$recruitment_lists'.")

  out <- list(
    recruitment_lists = recruitment_lists,
    recruitment_table = table_output$recruit_table,
    recruitment_kable = table_output$recruit_kable,
    recruitment_numbers = table_output$recruit_numbers
  )

  class(out) <- "generalizeR_recruit"

  return(invisible(out))
}

#' Print method for "generalizeR_recruit" class
#'
#' @param x An object of class "generalizeR_recruit"
#' @param ... Other arguments passed to or from other methods
#' @return A summary of the recruitment tables
#'
#' @export print.generalizeR_recruit
#' @export

print.generalizeR_recruit <- function(x, ...) {

  cat("\nA generalizeR_recruit object: \n\n")

  print(x$recruitment_table)

}


#' Internal function that performs the actual recruitment calculations for \code{recruit()}
#'
#' Intended only to be called within \code{recruit()}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 1
#'
#' @param stratify_output list of class generalizeR_stratify
#' @return The function returns a list of dataframes, each corresponding to a recruitment list for a different stratum
#'
#' @md

.make.recruitment.lists <- function(stratify_output) {
  recruitment_lists <- list(NULL)

  for (i in 1:stratify_output$n_strata) {
    dat1 <- stratify_output$pop_data_by_stratum %>%
      dplyr::filter(Stratum == i)
    idvar_values <- dat1 %>% dplyr::select(all_of(stratify_output$idvar))

    dat2 <- dat1 %>%
      dplyr::select(-c(all_of(stratify_output$idvar), Stratum)) %>%
      dplyr::mutate_all(as.numeric)

    mu <- dat2 %>% purrr::map_dbl(mean)
    v <- var(dat2)
    a <- diag(v)

    if (anyNA(a) || any(a == 0)) {
      a[which(is.na(a) | a == 0)] <- 0.00000001 # Make it say NA
    }

    cov.dat <- diag(a)

    tryCatch( ma.s <- stats::mahalanobis(dat2, mu, cov.dat, tol=1e-20),
             error = function(e){
               stop(crayon::yellow("\nERROR: It appears that the covariance matrix associated with your selected variables is linearly dependent. That is, some combination(s) of selected variables are perfect predictors of one another. Please make sure you haven't included duplicate information across variables (columns), choosing fewer variables if necessary, and re-run stratify with a different set of variables."))
               })

    dat3 <- data.frame(idvar_values,
      dat2,
      distance = ma.s,
      Stratum = dat1$Stratum
    ) %>%
      dplyr::tibble()

    recruitment_lists[[i]] <- dat3 %>%
      dplyr::arrange(distance) %>%
      dplyr::mutate(Rank = seq.int(nrow(dat3))) %>%
      dplyr::select(Rank, all_of(stratify_output$idvar))
  }

  return(recruitment_lists)
}

#' Internal function that formats a summary table of stata-based recruitment results for \code{recruit()}
#'
#' Intended only to be called within \code{recruit()}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 1
#'
#' @param stratify_output list of class generalizeR_stratify
#' @param sample_size integer, total sample size selected to guide recruitment
#' @return The function returns a list of formatted recruitment results for user interpretation: a table of results, a corresponding kable, and a the number of units selected from each stratum
#'
#' @md

.make.recruitment.table <- function(stratify_output,
                                    sample_size) {
  n_strata <- stratify_output$n_strata

  recruit_table <- stratify_output$heat_data %>%
    dplyr::select(Stratum, n) %>%
    dplyr::distinct(Stratum, .keep_all = TRUE) %>%
    dplyr::mutate(
      Population_Units = n,
      Proportion = round(n / (dim(stratify_output$pop_data_by_stratum)[1]), digits = 3),
      Recruit_Number = .round.preserve.sum(sample_size * Proportion)
    ) %>%
    dplyr::filter(Stratum != "Population") %>%
    dplyr::select(Stratum, Population_Units, Proportion, Recruit_Number) %>%
    tidyr::pivot_longer(names_to = " ", cols = c(Population_Units, Proportion, Recruit_Number)) %>%
    tidyr::pivot_wider(names_from = Stratum, names_prefix = "Stratum ") %>%
    dplyr::mutate(" " = c("Population Units", "Sampling Proportion", "Recruitment Number")) %>%
    as.data.frame()

  recruit_numbers <- recruit_table %>%
    dplyr::filter(` ` == "Recruitment Number") %>%
    dplyr::select(-c(` `)) %>%
    as.numeric()

  recruit_header <- c(1, n_strata)
  names(recruit_header) <- c(" ", "Stratum")

  recruit_kable <- recruit_table %>%
    kableExtra::kbl(
      caption = "Recruitment Table",
      align = "c",
      col.names = c(" ", 1:n_strata)
    ) %>%
    kableExtra::column_spec(1, bold = TRUE, border_right = TRUE) %>%
    kableExtra::row_spec(3, background = "#5CC863FF") %>%
    kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    kableExtra::add_header_above(recruit_header)

  output <- list(
    recruit_table = recruit_table,
    recruit_kable = recruit_kable,
    recruit_numbers = recruit_numbers
  )

  return(invisible(output))
}

#' Internal function that asks user to specify sample size for guided version of \code{recruit()} and validates user-specified inputs
#'
#' Intended only to be called within \code{recruit()}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 1
#'
#' @param valid_inputs vector, enumerates the possible number of units to recruit, based on the population size
#' @return The function returns a list of formatted recruitment results for user interpretation: a table of results, a corresponding kable, and a the number of units selected from each stratum
#'
#' @md

.get.sample.size <- function(valid_inputs) {
  repeat {
    sample_size <- suppressWarnings(readline(prompt = "How many units would you like to recruit in total? ") %>% as.numeric())

    if (sample_size %in% valid_inputs) {
      break
    } else {
      cat(crayon::red("\nInvalid input. The number of units you wish to recruit must be an integer \nbetween 1 and the total number of units in your population ("),
          crayon::red(length(valid_inputs)),
          crayon::red(").\n\n"),
          sep = ""
          )
    }
  }

  return(invisible(sample_size))
}

#' Internal function used in the calculations for creating recruitment lists, under \code{.make.recruitment.lists()}. Rounds number of observations to be recruited from each stratum to a whole number, while adjusting rounded results to preserve the total.
#'
#' Intended only to be called within \code{recruit()}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 1
#'
#' @param x list of numeric values; enumerates the possible number of units to recruit, based on the population size
#' @param digits integer; signals number of 0 digits to place after each numeric in x (e.g. digits = 2 will multiply x by 100)
#' @return The function returns a list of integers, which will represent the number of observations to be sampled from each stratum
#'
#' @md

.round.preserve.sum <- function(x,
                                digits = 0) {
  up <- 10^digits
  x <- x * up
  y <- floor(x)
  indices <- utils::tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  return(y / up)
}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("Stratum", "distance", "Rank", "n", "Proportion", "Population_Units",
                                                         "Recruit_Number", " "))
