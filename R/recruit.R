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

recruit <- function(x, guided = TRUE, sample_size = NULL, save_as_csv = FALSE) {

  if(!inherits(x, "generalizer_output")) {

    stop("Argument 'x' must be an object of class \"generalizer_output\", \ncreated by running stratify().")
  }

  n <- NULL

  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width") - 1L), collapse = " "))

  pop_data_by_stratum <- x$pop_data_by_stratum

  pop_size <- pop_data_by_stratum %>% nrow()

  valid_inputs <- 1:pop_size

  idnum <- x$idnum

  n_strata <- x$n_strata

  cat("The generalizer output you've supplied consists of ", paste(bold(pop_size)), " population units \ndivided into ", paste(bold(x$n_strata)), " strata along these variables:\n", paste(blue$bold(x$variables), collapse = ", "), ".", sep = "")

  cat("\n\nGiven the number of units that you wish to recruit (your desired sample size), \nthis function can tell you how many units to recruit from each stratum and \ngenerate recruitment lists.\n")

  #### GUIDED VERSION PART 1 ####

  if(guided == TRUE) {

    satisfied <- 0

    while(satisfied == 0) {

      cat("\n")
      sample_size <- readline(prompt = "Number of units to recruit: ") %>% as.numeric()

      if(!(sample_size %in% valid_inputs)) {

        cat(red("Invalid input. The number of units you wish to recruit must be an integer \nbetween 1 and the total number of units in your population ("),
            red(pop_size),
            red(")."),
            sep = "")
      }

      else {

        satisfied <- 1
      }
    }
  }

  #### NON-GUIDED VERSION PART 1 ####

  else {

    if(is.null(sample_size)) {

      stop("You must specify the number of units that you wish to recruit.")
    }

    if(!(sample_size %in% valid_inputs)) {

      errorMsg = paste("The number of units you wish to recruit must be an integer between 1 and \nthe total number of units in your population (",
                       pop_size,
                       ").",
                       sep = "")

      stop(errorMsg)
    }

    if(!inherits(x, "generalizer_output")) {

      stop("Argument 'x' must be an object of class \"generalizer_output\", \ncreated by running stratify().")
    }
  }

  #### CREATE RECRUITMENT LISTS ####

  recruitment_lists <- list(NULL)

  for(i in 1:n_strata) {
    dat1 <- pop_data_by_stratum %>%
      dplyr::filter(Stratum == i)
    idvar <- dat1 %>% select(all_of(idnum))
    dat2 <- dat1 %>% select(-c(all_of(idnum), Stratum)) %>%
      mutate_all(as.numeric)

    mu <- dat2 %>% map_dbl(mean)
    v <- var(dat2)
    a <- diag(v)

    if(any(a == 0)) { a[which(a == 0)] <- 0.00000001 }
    cov.dat <- diag(a)
    ma.s <- mahalanobis(dat2, mu, cov.dat)
    dat3 <- data.frame(idvar, dat2, distance = ma.s, Stratum = dat1$Stratum) %>% tibble()

    recruitment_lists[[i]] <- dat3 %>% # Produces a list of data frames, one per stratum, sorted by
      # distance (so the top N schools in each data frame are the "best," etc.)
      arrange(distance) %>%
      mutate(rank = seq.int(nrow(dat3))) %>%
      select(rank, all_of(idnum))
  }

  cat("\n")
  cat(paste(n_strata), "recruitment lists have been generated, one per stratum.")
  cat(" Each list contains \nthe ID information for the units, ranked in order of desirability.\n")

  if(guided == TRUE) {

    cat("\nThe top 6 rows of the first recruitment list are shown below.\n\n")
    recruitment_lists[[1]] %>% head() %>% print()
  }

  #### CREATE RECRUITMENT TABLE ####

  recruit_table <- x$heat_data %>% select(Stratum, n) %>%
    distinct(Stratum, .keep_all = TRUE) %>%
    mutate(Population_Units = n,
           Proportion = round(n/(dim(pop_data_by_stratum)[1]), digits = 3),
           Recruit_Number = round_preserve_sum(sample_size * Proportion)) %>%
    filter(Stratum != "Population") %>%
    select(Stratum, Population_Units, Proportion, Recruit_Number) %>%
    pivot_longer(names_to = " ", cols = c(Population_Units, Proportion, Recruit_Number)) %>%
    pivot_wider(names_from = Stratum, names_prefix = "Stratum ") %>%
    mutate(" " = c("Population Units", "Sampling Proportion", "Recruitment Number")) %>%
    as.data.frame()

  recruit_header <- c(1, n_strata)
  names(recruit_header) <- c(" ", "Stratum")

  recruit_kable <- recruit_table %>% kbl(caption = "Recruitment Table",
                                         align = "c",
                                         col.names = c(" ", 1:n_strata)) %>%
    column_spec(1, bold = TRUE, border_right = TRUE) %>%
    row_spec(3, background = "#5CC863FF") %>%
    kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    add_header_above(recruit_header)


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

    save_as_csv = menu(choices = c("Yes", "No"), title = cat("Would you like to save the recruitment lists as .csv files?"))

    if(save_as_csv == TRUE) {

      cat("\nThe lists will be saved as 'recruitment_list_#', one for each stratum.\n")
      cat("Where should they be saved?\n\n")
    }
  }

  #### SAVE RECRUITMENT LISTS

  if(save_as_csv == TRUE) {

    filepath <- ifelse(guided, easycsv::choose_dir(), "")

    for(i in 1:(n_strata)) {

      filename <- paste(filepath, "Recruitment_list_", i, ".csv", sep = "")
      write_csv(recruitment_lists[[i]], file = filename)
    }

    #### GUIDED VERSION PART 3 ####

    if(guided == TRUE) {

      cat("Lists saved successfully.\n\n")
    }

    #### NON-GUIDED VERSION PART 2 ####

    if(guided == FALSE) {

      cat("You've chosen to save your recruitment lists as .csv files. The lists have been \nsaved as 'recruitment_list_#', one for each stratum, inside your current working \ndirectory.\n\n")
    }
  }

  cat("If you have stored the output of 'recruit()' in an object, you can use it to \naccess these lists by typing the name of the object followed by \n'$recruitment_lists'.")

  output <- list(recruitment_lists = recruitment_lists,
                 recruitment_table = recruit_table,
                 recruitment_kable = recruit_kable)

  return(invisible(output))
}
