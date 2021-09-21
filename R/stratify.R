#' Stratify a Population Data Frame
#'
#' This function takes as input any data frame that you want to stratify into clusters. Typically, the goal of such stratification is sampling for generalizability. This function, and the others in this package, are designed to mimic the website https://www.thegeneralizer.org/.
#'
#' @param data The R object containing your population data frame
#' @param guided logical; defaults to TRUE. Whether the function should be guided (ask questions and behave interactively throughout) or not. If set to FALSE, must provide values for other arguments below
#' @param n_strata defaults to NULL. If guided is set to FALSE, must provide a number of strata to cluster population into
#' @param variables defaults to NULL. If guided is set to FALSE, must provide a character vector of the names of stratifying variables (from population data frame)
#' @param idnum defaults to NULL. If guided is set to FALSE, must provide a character vector of the name of the ID variable (from population data frame)
#' @return The function returns a list of class "generalizer_output" that can be provided as input to \code{recruit()}. More information on the components of this list can be found above under "Details."
#' @details The list contains 11 components: \code{x2}, \code{solution}, \code{n_strata}, \code{recruitment_lists}, \code{population_summary_stats2}, \code{summary_stats}, \code{summary_stats2}, \code{heat_data}, \code{heat_plot_final}, \code{idnum}, and \code{variables}.
#'
#' \itemize{
#' \item{\code{x2}: }{a tibble with number of rows equal to the number of rows in the inference population (\code{data}) and number of columns equal to the number of stratifying variables (dummy-coded if applicable) plus the ID column (\code{idnum}) and a column representing stratum membership, \code{Stratum}}
#' }
#' @export
#' @importFrom graphics par
#' @importFrom stats mahalanobis median na.omit sd var
#' @importFrom utils menu select.list
#' @importFrom crayon red yellow blue bold
#' @importFrom janitor clean_names
#' @importFrom ggplot2 ggplot aes geom_bar xlab labs geom_histogram geom_text geom_label geom_hline scale_fill_gradientn scale_x_discrete expand_limits geom_tile element_blank element_text theme
#' @importFrom ggthemes theme_base
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggnewscale new_scale
#' @importFrom viridis viridis turbo plasma
#' @importFrom kableExtra kbl kable_styling
#' @importFrom tidyr pivot_longer unite_
#' @importFrom dplyr count arrange filter mutate summarise_all summarize_if left_join group_by select select_if all_of mutate_all case_when bind_rows bind_cols distinct everything
#' @importFrom tibble tibble_row add_row tibble
#' @importFrom ClusterR KMeans_rcpp
#' @importFrom purrr map_dbl map_df negate
#' @importFrom magrittr %>%
#' @importFrom cluster daisy
#' @importFrom stringr str_sub
#' @importFrom grid unit
#' @importFrom tidyselect contains
#' @references
#' Tipton, E. (2014). Stratified sampling using cluster analysis: A sample selection strategy for improved generalizations from experiments. *Evaluation Review*, *37*(2), 109-139.
#'
#' Tipton, E. (2014). How generalizable is your experiment? An index for comparing experimental samples and populations. *Journal of Educational and Behavioral Statistics*, *39*(6), 478-501.
#' @examples
#' \donttest{
#' \dontrun{
#' # To get sample data; must first be installed using install_github("katiecoburn/generalizeRdata")
#' library(generalizeRdata)
#'
#' # Guided:
#' stratify(ipeds)
#'
#' # Not guided:
#' stratify(ipeds, guided = FALSE, n_strata = 4,
#'    variables = c("pct_female", "pct_white"), idnum = "unitid")
#' }
#' }
#' @md

stratify <- function(data, guided = TRUE, n_strata = NULL, variables = NULL,
                     idnum = NULL, seed = 7835, verbose = TRUE){

  skim_variable <- skim_type <- variable <- NULL
  type <- Stratum <- n <- mn <- deviation <- NULL
  data_name <<- data %>% expr_text() # Store name of dataset as global variable so it can be accessed by stratify_basic(). Must be done before function argument 'data' is evaluated for the first time.

  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width") - 1L), collapse = " "));

  data <- data %>%
    mutate(across(where(is_character) & !idnum, as_factor)) # immediately convert all character variables (except the id variable if it is one) to factors

  # Here begins the guided wrapper for the function -------------------------

  if(guided == TRUE) {

    ## Check ##
    if (!is.null(n_strata) | !is.null(variables) | !is.null(idnum)) {
      stop(simpleError("Don't specify n_strata, variables, or idnum as arguments if you are running the guided version of this function."))
    }


    # 1) Introduction: instructions to store object ---------------------------

    cat(bold("\nWelcome to stratify! \n"))

    cat("\nIf you want to adjust or restrict your inference population \n(e.g., if you are interested in only one location, etc.), \nmake sure that you have altered the data frame appropriately. \nIf you need to alter your data frame, you can exit this \nfunction, use ", blue$bold("dplyr::filter()"), ", and then return.\n", sep = "")

    cat(bold("\nTo store your results, make sure you assign \nthis function to an object.\n\n"))

    if(menu(choices = c("Yes", "No"), title = cat("I have assigned this function to an object and wish to proceed:")) == 1) {

    }
    else {
      stop(simpleError(blankMsg))
    }
    cat("Your chosen inference population is the '",
        data_name, "' dataset.", sep = "")

    cat("\n")
    cat("\n")


    # 2) Selection of idnum ---------------------------------------------------

    is_valid_variable_name <- FALSE

    while(is_valid_variable_name == FALSE) {
      idnum <- readline(prompt = "Enter the name of the ID Variable in your dataset: ")

      ## Check ##
      if(!idnum %in% names(data)) {
        cat(red("ERROR: We could not find that variable. Please make sure your \ndataset contains an ID variable."))
        next
      }

      is_valid_variable_name <- TRUE
    }

    id <- data %>% select(all_of(idnum))

    # 3) Selection of variables (for stratification) -----------------------------------------

    make_var_overview(data %>% select(-all_of(idnum)))

    cat("\nIn the Viewer pane to the right you will find a table that displays each \nvariable in your dataset along with its object type and number of levels \n(only relevant for factor variables). ",
        yellow$bold("Please note that any character \nvariables that may have been present in your dataset have been \nautomatically converted to factor variables.\n"),
        sep = "")

    cat("\nYou're now ready to select your stratification variables. The following \nare the variables available in your dataset. Which key variables do you \nthink may explain variation in your treatment effect? Typically, studies \ninclude 4-6 variables for stratification.", yellow$bold("You must choose at least 2 \nvariables and you may not choose any factor variables with more than 4 \nlevels.\n"))

    variables_are_correct <- 0

    while(variables_are_correct != 1){

      names <- names(data %>% select(-all_of(idnum)))
      variables <- select.list_CUSTOMIZED(choices = names,
                                          graphics = FALSE,
                                          multiple = TRUE)

      data_subset <- data %>%
        select(all_of(variables))

      ## Check to see there are no cat vars with > 4 factors ##
      factor_levels_over_4 <- (data_subset %>% select_if(is.factor) %>% sapply(nlevels) > 4L) %>%
        which() %>% names()

      if(!is_empty(factor_levels_over_4)) {
        cat(red("ERROR: The following factor variables have more than 4 levels:\n"),
            paste(blue$bold(factor_levels_over_4), collapse = ", "),
            red("\n4 is the maximum number of levels this function will allow a factor to have."),
            red("\nPlease exit out of this function (Press 'Esc') and re-code your desired factor"),
            red("\nlevels from these variables as dummy variables (see the package 'fastDummies').\n"), sep = "")
        next
      }

      cat("You have selected the following stratifying variables:\n",
          paste(blue$bold(variables), collapse = ", "),
          "\n\n",
          sep = "")

      make_var_overview(data_subset, print_to_console = TRUE)

      if(menu(choices = c("Yes", "No"), title = cat("\nIs this correct?")) == 1) {
        variables_are_correct <- 1
      }
      else {
        variables_are_correct <- 0
      }
    }


# 3b) Let user know about missing observations ----------------------------------

    missing_obs_table <- data.frame(id, data_subset) %>%
      sapply(function(x) sum(is.na(x))) %>%
      as.data.frame() %>%
      `colnames<-`("n_missing") %>%
      rownames_to_column("variable")

    n_missing <- data.frame(id, data_subset) %>% filter(if_any(everything(), is.na)) %>% nrow()

    data_subset <- data.frame(id, data_subset) %>% na.omit() %>% select(-all_of(idnum))

    cat(paste0("The following table shows how many observations are missing for the variables you have chosen. \n",
    yellow$bold("In total,", n_missing, "rows will be dropped from the inference population before stratification due to \nmissing observations. "),
    "Note that there may be multiple missing observations per row, meaning that \nthe total number of dropped observations won't necessarily be the sum of the number of \nmissing observations in each variable.\n\n "))

    missing_obs_table %>%
      kbl(caption = "Missing Observations by Variable",
          align = "l",
          col.names = c("Variable", "Number Missing")) %>%
      kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
      print()

    missing_obs_table %>% print()

    cat("\n")
    readline(prompt = "Press [enter] to proceed once you have viewed the missing observations.")



    # insert a section where we say:
    # no. of observations in population based on variables
    # alert them of no. of observations omitted due to missing information in those variables

    # 4) Overview of categorical variables ------------------------------------
    ## NOTE: These cat variables get converted into factors anyway

    cat_data <- data_subset %>% select_if(is.factor)
    cat_data_vars <- names(cat_data)
    if(dim(cat_data)[2] >= 1) {
      cat_data_plot <- data.frame(cat_data)
      cat("Please review the descriptive statistics of your categorical variables (factors).\n",
          "Bar charts and tables for each variable will also be printed in the Plots and \nViewer panes to the right.\n", sep = "")

      n_cat_vars <- ncol(cat_data_plot)
      fill_colors_cat <- plasma(n_cat_vars, alpha = 0.7, direction = sample(c(-1, 1), size = 1)) %>%
        sample()
      outline_colors_cat <- turbo(n_cat_vars) %>% sample()

      for(i in 1:n_cat_vars) {
        var_name <- cat_data_vars[i]
        levels(cat_data_plot[[var_name]]) <- str_wrap(levels(cat_data_plot[[var_name]]), width = 10)
        barfig <- cat_data_plot %>%
          group_by(across(all_of(var_name))) %>%
          summarise(count = n()) %>%
          mutate(ordered_factor = fct_reorder(.[[var_name]], count)) %>%
          ggplot(aes(x = ordered_factor, y = count)) +
          geom_col(fill = fill_colors_cat[i],
                   color = outline_colors_cat[i]) +
          theme_minimal() +
          xlab(var_name) +
          labs(title = paste("Bar Chart of", var_name))
        print(barfig)
        par(ask = TRUE)
        cat("\nNumber of Observations in Levels of Factor ",
            paste(blue$bold(var_name)),
            ":\n",
            sep = "")
        cat_data_table <- table(cat_data_plot[,i])
        cat_data_table %>% print()
        cat_data_table %>%
          kbl(col.names = c("Level", "Frequency"),
              caption = paste("Number of Observations in Levels of Factor ", var_name),
              align = "l") %>%
          kable_styling(c("striped", "hover")) %>%
          print()
      }
    }



    # 5) Overview of cont data ------------------------------------------------

    cont_data <- data_subset %>%
      select_if(negate(is.factor))
    cont_data_vars <- names(cont_data)

    if(dim(cont_data)[2] >= 1L) {

      cat("\nPlease review the descriptive statistics of your continuous variables. Histograms \nand tables for each variable will also be printed in the Plots and Viewer panes \nto the right. \n\n")

      n_cont_vars <- ncol(cont_data)
      fill_colors_cont <- viridis(n_cont_vars, alpha = 0.7, direction = sample(c(-1, 1), size = 1)) %>%
        sample()
      outline_colors_cont <- turbo(n_cont_vars) %>% sample()

      for(i in 1:n_cont_vars) {
        cont_data_plot <- cont_data %>% data.frame()
        suppressWarnings(
          suppressMessages(
            hist <- cont_data_plot %>% ggplot(aes(x = cont_data_plot[,i])) +
              geom_histogram(bins = 30,
                             fill = fill_colors_cont[i],
                             color = outline_colors_cont[i]) +
              theme_minimal() +
              xlab(cont_data_vars[i]) +
              labs(title = paste("Histogram of", cont_data_vars[i]))
          )
        )
        print(hist)
        par(ask = TRUE)
      }

      sumstats <- cont_data %>%
        map_df(function(x) {
          tibble(min = min(x), pct50 = median(x), max = max(x), mean = mean(x), sd = sd(x))
        }) %>%
        mutate_all(round, digits = 3) %>%
        mutate(variable = cont_data_vars) %>%
        select(variable, everything()) %>%
        clean_names() %>%
        data.frame()

      sumstats %>% print(row.names = FALSE)

      sumstats %>%
        kbl(col.names = c("Variable", "Min", "Pct50", "Max", "Mean", "SD")) %>%
        kable_styling(c("striped", "hover")) %>%
        print()
    }
    par(ask = FALSE)


    # 6) Selection of no. of strata -------------------------------------------

    cat("\nStratification will help you develop a recruitment plan so that your study will \nresult in an unbiased estimate of the ", bold("average treatment effect (ATE)"), ". Without \nusing strata, it is easy to end up with a sample that is very different from your \ninference population. \n\nGeneralization works best when strata are ", bold("homogeneous"), ". That means units within \neach stratum are almost identical in terms of relevant variables.\n\n", sep = "")

    cat("Enter the number of strata in which you wish to divide your population. Typically, ",
        bold("\nthe more strata"),
        ",",
        bold("the better"),
        "; with fewer strata, units in each stratum are no longer \nidentical. However, increasing ",
        "the number of strata uses more resources, because \nyou must sample a given number of units ",
        "from each stratum. Choosing 4-6 strata is \ncommon. \n\nTry a few numbers and choose the 'best' one for you.",
        sep = "")


    ## Catch ##

    n_strata_correct <- FALSE

    while(!n_strata_correct) {

      n_strata <- suppressWarnings(as.numeric(readline(prompt = "# of strata: ")))

      if(is.na(n_strata) || n_strata <= 1) {
        cat(red("ERROR: The number of strata must be a single positive integer greater than 1.\n"))
        next
      }

      if(n_strata%%1==0) {
        n_strata <- round(n_strata)
      }
      n_strata_correct <- TRUE
    }
  }

  # here is where guided_loop ENDS.


  # 7) Stratify_basic call  -------------------------------------------------

    overall_output <- stratify_basic(data = data, n_strata = n_strata, variables = variables,
                                     idnum = idnum, seed = 7835, verbose = verbose)

    rm(list = deparse(substitute(data_name)), envir=.GlobalEnv) # delete object data_name from global environment

  # 8) Final message for guided version -------------------------------------

  if(guided == TRUE) {
    cat(blue$bold("Congratulations, you have successfully grouped your data into", n_strata, "strata!\n"))
    cat("You can pull up the results anytime by passing your stratify_object into summary().\n\n")

    readline(prompt = "Press [enter] to view the results")

    print(summary(overall_output))

    # if(menu(choices = c("Yes", "No"), title = cat("\nWould you like to go back and specify a different number of strata? If you specify \n'No' the stratification process will end and you can proceed to use the output in \n'recruit()' provided that it has been assigned to an object.")) == 2){
    #
    #   satisfied <- 1
    #
    # }else{
    #
    #   satisfied <- 0
    #
    # }

  }

  class(overall_output) <- c("generalizer_output")

  return(invisible(overall_output))

}

