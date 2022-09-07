#' Recruit Units from a Population for Sampling
#'
#' This function works with the output of 'stratify()'. The user provides the number of units they wish to sample from their population dataset. The function tells the user how many to sample from each stratum and generates recruitment lists, one per stratum, which can either be saved to .csv files in any given directory or accessed later on.
#'
#' This function, and the others in this package, are designed to mimic the website https://www.thegeneralizer.org/.
#'
#' @param x output from 'stratify()', of S3 class 'generalizer_output'
#' @param guided logical; defaults to TRUE. Whether the function should be guided (ask questions and behave interactively throughout) or not. If set to FALSE, must provide values for other arguments below
#' @param sample_size defaults to NULL. If guided is set to FALSE, must provide a number of units to sample
#' @param save_as_csv defaults to NULL. If guided is set to FALSE, specify whether or not to save recruitment lists to working directory; TRUE or FALSE
#' @return A three-element list containing the recruitment lists and the recruitment table/kable
#' @export
#' @importFrom readr write_csv
#' @importFrom easycsv choose_dir

recruit <- function(x,
                    guided = TRUE,
                    sample_size = NULL,
                    save_as_csv = FALSE) {

  stopifnot("Argument 'x' must be an object of class \"generalizer_stratify\", created by running stratify()." =
              inherits(x, "generalizer_stratify"))

  pop_data_by_stratum <- x$pop_data_by_stratum

  pop_size <- pop_data_by_stratum %>% nrow()

  valid_inputs <- 1:pop_size

  idvar <- x$idvar

  n_strata <- x$n_strata

  if(guided == TRUE) {

    cat(crayon::bold("\nWelcome to recruit()! \n\n"))
  }

  else {

    #### NON-GUIDED VERSION PART 1 ####

    stopifnot("You must specify the number of units that you wish to recruit." = !is.null(sample_size))

    recruit_num_error_msg <- paste("The number of units you wish to recruit must be an integer between 1 and \nthe total number of units in your population (",
                                   pop_size,
                                   ").",
                                   sep = "")
    stopifnot(recruit_num_error_msg = (sample_size %in% valid_inputs))
  }

  cat("The 'generalizer_stratify' object you've supplied consists of ", paste(crayon::bold(pop_size)), " population units \ndivided into ", paste(crayon::bold(n_strata)), " strata along these variables:\n\n", paste(blue$bold(x$variables), collapse = ", "), ".", sep = "")

  recruitment_lists <- .make.recruitment.lists(x)

  cat("\n")
  cat(paste(n_strata), "recruitment lists have been generated, one per stratum.")
  cat(" Each list contains \nthe ID information for the units, ranked in order of desirability.\n\n")

  if(guided == TRUE) {

    cat("The top 6 rows of the first recruitment list are shown below.\n\n")
    recruitment_lists[[1]] %>%
      utils::head() %>%
      data.frame() %>%
      print(row.names = FALSE)

    cat("\n\nGiven the number of units that you wish to recruit (your desired sample size), this \nfunction can tell you how many units to recruit from each stratum.\n")

    repeat {

      cat("\n")
      sample_size <- readline(prompt = "How many units would you like to recruit? ") %>% as.numeric()

      if(sample_size %in% valid_inputs) {
        break
      }

      else {
        stop(cat(crayon::red("Invalid input. The number of units you wish to recruit must be an integer \nbetween 1 and the total number of units in your population ("),
                 crayon::red(pop_size),
                 crayon::red(")."),
                 sep = ""))
      }
    }
  }

  table_output <- .make.recruitment.table()

  #### CREATE RECRUITMENT TABLE ####

  recruit_table <- x$heat_data %>% dplyr::select(Stratum, n) %>%
    dplyr::distinct(Stratum, .keep_all = TRUE) %>%
    dplyr::mutate(Population_Units = n,
           Proportion = round(n/(dim(pop_data_by_stratum)[1]), digits = 3),
           Recruit_Number = .round.preserve.sum(sample_size * Proportion)) %>%
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

  recruit_kable <- recruit_table %>% kableExtra::kbl(caption = "Recruitment Table",
                                         align = "c",
                                         col.names = c(" ", 1:n_strata)) %>%
    kableExtra::column_spec(1, bold = TRUE, border_right = TRUE) %>%
    kableExtra::row_spec(3, background = "#5CC863FF") %>%
    kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    kableExtra::add_header_above(recruit_header)

  cat("\nThe following table (also shown in the Viewer pane to the right) displays \nthe stratum sizes, their proportion relative to the total population size, \nand consequent recruitment number for each stratum. ")
  cat("Ideally, units should be \nrecruited across strata according to these numbers.")
  cat(" Doing so will lead to the \nleast amount of bias and no increase in standard errors. ")
  cat("Note that the \nrecruitment numbers have been rounded to integers in such a way as to ensure \ntheir sum equals the desired total sample size.\n\n")

  cat(blue$bold("Recruitment Table\n"))

  print(recruit_table, row.names = FALSE)
  print(recruit_kable)

  cat("\nAttempt to recruit units starting from the top of each recruitment list. If you \nare unsuccessful in recruiting a particular unit, move on to the next one in the \nlist and continue until you have reached the ideal recruitment number in each \nstratum.\n\n", sep = "")

  #### GUIDED VERSION PART 2 ####

  if(guided == TRUE) {

    save_as_csv = utils::menu(choices = c("Yes", "No"),
                              title = cat("Would you like to save the recruitment lists as .csv files?"))

    if(save_as_csv == TRUE) {

      cat("\nThe lists will be saved as 'Recruitment_list_#', one for each stratum.\n")
      cat("Where should they be saved?\n\n")
    }
  }

  #### SAVE RECRUITMENT LISTS

  if(save_as_csv == TRUE) {

    filepath <- ifelse(guided, easycsv::choose_dir(), "")

    for(i in 1:(n_strata)) {

      filename <- paste(filepath, "Recruitment_list_", i, ".csv", sep = "")
      readr::write_csv(recruitment_lists[[i]], file = filename)
    }

    #### GUIDED VERSION PART 3 ####

    if(guided == TRUE) {

      cat("Lists saved successfully!\n\n")
    }

    #### NON-GUIDED VERSION PART 2 ####

    if(guided == FALSE) {

      cat("You've chosen to save your recruitment lists as .csv files. The lists have been \nsaved as 'Recruitment_list_#', one for each stratum, inside your current working \ndirectory.\n\n")
    }
  }

  cat("If you have stored the output of 'recruit()' in an object, you can use it to \naccess these lists by typing the name of the object followed by \n'$recruitment_lists'.")

  out <- list(recruitment_lists = recruitment_lists,
                 recruitment_table = recruit_table,
                 recruitment_kable = recruit_kable,
                 recruitment_numbers = recruit_numbers)

  class(out) <- "generalizer_recruit"

  return(invisible(out))
}

#'
#' Internal function that performs recruitment calculations
#'

.make.recruitment.lists <- function(x) {

  recruitment_lists <- list(NULL)

  for(i in 1:x$n_strata) {

    dat1 <- x$pop_data_by_stratum %>%
      dplyr::filter(Stratum == i)
    idvar_values <- dat1 %>% dplyr::select(all_of(idvar))

    dat2 <- dat1 %>%
      dplyr::select(-c(all_of(idvar), Stratum)) %>%
      dplyr::mutate_all(as.numeric)

    mu <- dat2 %>% purrr::map_dbl(mean)
    v <- var(dat2)
    a <- diag(v)

    if(any(a == 0)) { a[which(a == 0)] <- 0.00000001 }
    cov.dat <- diag(a)

    ma.s <- stats::mahalanobis(dat2, mu, cov.dat)

    dat3 <- data.frame(idvar_values,
                       dat2,
                       distance = ma.s,
                       Stratum = dat1$Stratum) %>%
      dplyr::tibble()

    # Produce a list of data frames, one per stratum, sorted by
    # distance (so the top N schools in each data frame are the "best," etc.)
    recruitment_lists[[i]] <- dat3 %>%
      dplyr::arrange(distance) %>%
      dplyr::mutate(rank = seq.int(nrow(dat3))) %>%
      dplyr::select(rank, all_of(idvar))
  }

  return(recruitment_lists)
}

  # cat("\n")
  # cat(paste(n_strata), "recruitment lists have been generated, one per stratum.")
  # cat(" Each list contains \nthe ID information for the units, ranked in order of desirability.\n")
  #
  # if(guided == TRUE) {
  #
  #   cat("\nThe top 6 rows of the first recruitment list are shown below.\n\n")
  #   recruitment_lists[[1]] %>%
  #     dplyr::head() %>%
  #     dplyr::data.frame() %>%
  #     dplyr::print(row.names = FALSE)
  # }

  #### CREATE RECRUITMENT TABLE ####

.make.recruitment.table <- function(x,
                                    sample_size) {

  recruit_table <- x$heat_data %>% dplyr::select(Stratum, n) %>%
    dplyr::distinct(Stratum, .keep_all = TRUE) %>%
    dplyr::mutate(Population_Units = n,
                  Proportion = round(n/(dim(x$pop_data_by_stratum)[1]), digits = 3),
                  Recruit_Number = .round.preserve.sum(sample_size * Proportion)) %>%
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

  recruit_kable <- recruit_table %>% kableExtra::kbl(caption = "Recruitment Table",
                                                     align = "c",
                                                     col.names = c(" ", 1:n_strata)) %>%
    kableExtra::column_spec(1, bold = TRUE, border_right = TRUE) %>%
    kableExtra::row_spec(3, background = "#5CC863FF") %>%
    kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    kableExtra::add_header_above(recruit_header)

  output = list(recruit_table, recruit_kable)

  return(invisible(output))
}

  # cat("\nThe following table (also shown in the Viewer pane to the right) displays \nthe stratum sizes, their proportion relative to the total population size, \nand consequent recruitment number for each stratum. ")
  # cat("Ideally, units should be \nrecruited across strata according to these numbers.")
  # cat(" Doing so will lead to the \nleast amount of bias and no increase in standard errors. ")
  # cat("Note that the \nrecruitment numbers have been rounded to integers in such a way as to ensure \ntheir sum equals the desired total sample size.\n\n")
  #
  # cat(blue$bold("Recruitment Table\n"))
  #
  # print(recruit_table, row.names = FALSE)
  # print(recruit_kable)
  #
  # cat("\nAttempt to recruit units starting from the top of each recruitment list. If you \nare unsuccessful in recruiting a particular unit, move on to the next one in the \nlist and continue until you have reached the ideal recruitment number in each \nstratum.\n\n", sep = "")


#'
#' Helper function that asks user to specify sample size for guided version
#'

.get.sample.size <- function(x) {

  cat("\n\nGiven the number of units that you wish to recruit (your desired sample size), this \nfunction tell you how many units to recruit from each stratum.\n")

  repeat {

    cat("\n")
    sample_size <- readline(prompt = "How many units would you like to recruit? ") %>% as.numeric()

    if(sample_size %in% valid_inputs) {
      break
    }

    else {
      stop(cat(crayon::red("Invalid input. The number of units you wish to recruit must be an integer \nbetween 1 and the total number of units in your population ("),
               crayon::red(pop_size),
               crayon::red(")."),
               sep = ""))
    }
  }

  return(invisible(sample_size))
}

#'
#' Unguided helper function
#'

.recruit.unguided <- function(x) {


  cat("The 'generalizer_stratify' object you've supplied consists of ", paste(crayon::bold(pop_size)), " population units \ndivided into ", paste(crayon::bold(x$n_strata)), " strata along these variables:\n\n", paste(blue$bold(x$variables), collapse = ", "), ".", sep = "")

  recruitment_lists <- .make.recruitment.lists(x)

  cat("\n")
  cat(paste(n_strata), "recruitment lists have been generated, one per stratum.")
  cat(" Each list contains \nthe ID information for the units, ranked in order of desirability.\n\n")

  return(recruitment_lists)

}
