#' Stratify a Population Data Frame
#'
#' The function \code{stratify()} takes as input any data frame with observations (rows) that you wish to stratify into clusters. Typically, the goal of such stratification is developing a sampling design for maximizing generalizability. This function, and the others in this package, are designed to mimic the website https://www.thegeneralizer.org/.
#'
#' @order 4
#'
#' @param data data.frame object containing the population data to be stratified (observations as rows); must include a unique id variable for each observation, as well as covariates.
#' @param guided logical, defaults to TRUE. Whether the function should be guided (ask questions and behave interactively throughout) or not. If set to FALSE, the user must provide values for other arguments below
#' @param n_strata integer, defaults to NULL. If guided is set to FALSE, must provide a number of strata in which to divide to cluster population
#' @param variables character, defaults to NULL. If guided is set to FALSE, must provide a character vector of the names of stratifying variables (from population data frame)
#' @param idvar character, defaults to NULL. If guided is set to FALSE, must provide a character vector of the name of the ID variable (from population data frame)
#' @param verbose logical, defaults to TRUE.
#' @return The function returns a list of class "generalizeR_stratify" that can be provided as input to \code{recruit()}. More information on the components of this list can be found above under "Details."
#' @details The list contains 14 components: \code{idvar}, \code{variables}, \code{dataset}, \code{n_strata}, \code{solution}, \code{pop_data_by_stratum}, \code{summary_stats}, \code{data_omitted}, \code{cont_data_stats}, \code{cat_data_levels}, \code{heat_data}, \code{heat_data_simple}, \code{heat_data_kable}, and \code{heat_plot}.
#'
#' \itemize{
#' \item{\code{pop_data_by_stratum}: }{a tibble with number of rows equal to the number of rows in the inference population (\code{data}) and number of columns equal to the number of stratifying variables (dummy-coded if applicable) plus the ID column (\code{idvar}) and a column representing stratum membership, \code{Stratum}}
#' }
#' @export
#' @importFrom graphics par
#' @importFrom assertthat assert_that not_empty on_failure is.count
#' @importFrom stats mahalanobis median na.omit sd var
#' @importFrom utils menu
#' @importFrom crayon red yellow blue bold
#' @importFrom janitor clean_names
#' @importFrom ggplot2 ggplot aes geom_bar xlab labs geom_histogram geom_text geom_col geom_label geom_hline scale_fill_gradientn scale_x_discrete expand_limits geom_tile element_blank element_text theme
#' @importFrom ggthemes theme_base
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggnewscale new_scale
#' @importFrom viridis viridis turbo plasma
#' @importFrom kableExtra kbl kable_styling add_header_above
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
#' @importFrom forcats fct_reorder as_factor
#' @importFrom rlang is_empty
#' @references
#' Tipton, E. (2014). Stratified sampling using cluster analysis: A sample selection strategy for improved generalizations from experiments. *Evaluation Review*, *37*(2), 109-139.
#'
#' Tipton, E. (2014). How generalizable is your experiment? An index for comparing experimental samples and populations. *Journal of Educational and Behavioral Statistics*, *39*(6), 478-501.
#' @examples
#' \donttest{
#' \dontrun{
#' # To get sample data; must first be installed using devtools::install_github("NUstat/generalizeRdata")
#' library(generalizeRdata)
#'
#' # Guided:
#' stratify(ipeds)
#'
#' # Not guided:
#' stratify(ipeds, guided = FALSE, n_strata = 4, variables = c("pct_female", "pct_white"), idvar= "unitid")
#' }
#' }
#' @md

stratify <- function(data = NULL,
                     guided = TRUE,
                     n_strata = NULL,
                     variables = NULL,
                     idvar = NULL,
                     verbose = TRUE) {

  # Check whether data object is of type 'data.frame'
  assertthat::on_failure(is.data.frame) <- function(call, env) {
    "You must pass an object of type 'data.frame' to the 'data' argument."
  }

  assertthat::assert_that(is.data.frame(data))

  # Store name of data as global variable so it can be accessed by .stratify.calculate(). Must be done before function argument 'data' is evaluated for the first time.
  data_name <<- data %>%
    lazyeval::expr_text()

  # Immediately convert all character variables to factors
  data <- data %>%
    dplyr::mutate(dplyr::across(where(rlang::is_character), forcats::as_factor))

  # Call guided version ------------------------------------------------------
  if (guided) {
    output <- .stratify.guided(
      data,
      verbose
    )
  }

  # Call unguided version ----------------------------------------------------
  else {
    output <- .stratify.unguided(
      data = data,
      n_strata = n_strata,
      variables = variables,
      idvar = idvar,
      verbose = verbose
    )
  }

  rm(
    list = deparse(substitute(data_name)),
    envir = .GlobalEnv
  ) # delete object data_name from global environment

  class(output) <- "generalizeR_stratify"

  return(invisible(output))
}

#' Internal function that performs the actual stratification calculations for \code{stratify()}
#'
#' Intended only to be called within \code{stratify_guided} and \code{stratify_unguided}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 1
#'
#' @param data data.frame object containing the population data to be stratified (observations as rows); must include a unique id variable for each observation, as well as covariates.
#' @param n_strata integer, number of strata in which to divide to cluster population
#' @param variables character, character vector of the names of stratifying variables (from population data frame)
#' @param idvar character, haracter vector of the name of the ID variable (from population data frame)
#' @param verbose logical, defaults to TRUE.
#' @return The function returns a list of class "generalizeR_stratify" that can be provided as input to \code{recruit()}. More information on the components of this list can be found above under "Details."
#' @details The list contains 14 components: \code{idvar}, \code{variables}, \code{dataset}, \code{n_strata}, \code{solution}, \code{pop_data_by_stratum}, \code{summary_stats}, \code{data_omitted}, \code{cont_data_stats}, \code{cat_data_levels}, \code{heat_data}, \code{heat_data_simple}, \code{heat_data_kable}, and \code{heat_plot}.
#'
#' @md

.stratify.calculate <- function(data,
                                n_strata,
                                variables,
                                idvar,
                                verbose) {

  # 1) Store all information given by user

  # Extract id variable
  id <- data %>%
    dplyr::select(tidyselect::all_of(variables), tidyselect::all_of(idvar)) %>%
    tidyr::drop_na() %>%
    dplyr::select(tidyselect::all_of(idvar))

  # Create a dataframe that includes only the variables selected by the user
  data_subset <- data %>%
    dplyr::select(tidyselect::all_of(variables))

  # Save all rows with missing observations in a new dataframe
  data_omitted <- data_subset %>%
    dplyr::filter(dplyr::if_any(tidyselect::everything(), is.na))

  # Drop rows with missing observations in the subsetted dataframe
  data_subset <- data_subset %>%
    tidyr::drop_na()

  # Store categorical and continuous variables in separate dataframes
  cat_data <- data_subset %>%
    dplyr::select_if(is.factor)

  cont_data <- data_subset %>%
    dplyr::select_if(purrr::negate(is.factor))

  # Make dummy variables out of factors
  if (dim(cat_data)[2] >= 1L) {
    cat_data_with_dummies <- cat_data %>%
      fastDummies::dummy_cols(remove_first_dummy = TRUE) %>%
      dplyr::select_if(purrr::negate(is.factor))

    data_full <- cont_data %>%
      cbind(cat_data_with_dummies) %>%
      janitor::clean_names()
  } else {
    data_full <- cont_data %>%
      cbind(cat_data) %>%
      janitor::clean_names()
  }

  var_names <- data_full %>% names()

  # 2) Save summaries of variables in dataset
  cont_data_stats <- cont_data %>%
    .make.cont.data.tbl()

  cat_data_levels <- cat_data %>%
    purrr::map_df(
      function(x) {
        tibble::tibble(
          Type = class(x),
          Levels = nlevels(x)
        )
      }
    ) %>%
    dplyr::mutate(Variable = names(cat_data)) %>%
    dplyr::select(Variable, tidyselect::everything()) %>%
    data.frame()

  # 3) Perform stratification - cluster with KMeans

  if (verbose == TRUE) {
    cat("\nThis might take a little while. Please bear with us.")

    distance <- cluster::daisy(data_full, metric = "gower") %>%
      suppressWarnings()

    cat("\n\nCalculated distance matrix.\n")

    solution <- ClusterR::KMeans_rcpp(as.matrix(distance),
      clusters = n_strata,
      verbose = TRUE
    )
  } else {
    distance <- cluster::daisy(data_full, metric = "gower") %>%
      suppressWarnings()

    solution <- ClusterR::KMeans_rcpp(as.matrix(distance),
      clusters = n_strata,
      verbose = FALSE
    )
  }

  # 4) Create various datasets

  pop_data_by_stratum <- data.frame(id,
    data_full,
    Stratum = solution$clusters
  ) %>%
    dplyr::select(Stratum, tidyselect::everything()) %>%
    dplyr::arrange(Stratum) %>%
    tibble::tibble()

  population_summary_stats2 <- pop_data_by_stratum %>%
    dplyr::select(tidyselect::all_of(var_names)) %>%
    dplyr::summarise_all(list(mean, sd)) %>%
    dplyr::mutate_all(round, digits = 3)

  population_summary_stats <- population_summary_stats2 %>%
    names() %>%
    stringr::str_sub(end = -5) %>%
    unique() %>%
    lapply(
      function(x) {
        tidyr::unite(population_summary_stats2,
          {{ x }},
          grep(x, names(population_summary_stats2), value = TRUE),
          sep = " / ",
          remove = TRUE
        ) %>%
          dplyr::select(tidyselect::all_of(x))
      }
    ) %>%
    dplyr::bind_cols()

  summary_stats <- pop_data_by_stratum %>%
    dplyr::select(tidyselect::all_of(var_names), Stratum) %>%
    dplyr::group_by(Stratum) %>%
    dplyr::summarize_if(is.numeric, mean) %>%
    dplyr::left_join(

      (pop_data_by_stratum %>%
         dplyr::select(tidyselect::all_of(var_names), Stratum) %>%
         dplyr::group_by(Stratum) %>%
         dplyr::summarize_if(is.numeric, sd)),
      by = "Stratum",
      suffix = c("_fn1", "_fn2")
    ) %>%
    dplyr::mutate_all(round, digits = 3)

  summary_stats2 <- summary_stats %>%
    dplyr::select(-Stratum) %>%
    names() %>%
    stringr::str_sub(end = -5) %>%
    unique() %>%
    lapply(
      function(x) {
        tidyr::unite(summary_stats,
          {{ x }},
          grep(x, names(summary_stats), value = TRUE),
          sep = " / ",
          remove = TRUE
        ) %>%
          dplyr::select(tidyselect::all_of(x))
      }
    ) %>%
    dplyr::bind_cols() %>%
    dplyr::mutate(Stratum = summary_stats$Stratum) %>%
    dplyr::select(Stratum, tidyselect::everything()) %>%
    dplyr::left_join(

      (pop_data_by_stratum %>%
         dplyr::group_by(Stratum) %>%
         dplyr::count()),
      by = "Stratum"
    ) %>%
    dplyr::mutate(Stratum = as.character(Stratum)) %>%
    tibble::add_row(tibble::tibble_row(
      Stratum = "Population",
      population_summary_stats,
      n = dim(pop_data_by_stratum)[1]
    ))

  simtab_m <- population_summary_stats2 %>%
    dplyr::select(contains("fn1"))

  names(simtab_m) <- names(simtab_m) %>%
    stringr::str_sub(end = -5)

  sd_tab <- summary_stats %>%
    dplyr::select(contains("fn2")) %>%
    tibble::add_row(tibble::tibble_row((population_summary_stats2 %>% dplyr::select(tidyselect::contains("fn2")))))

  names(sd_tab) <- names(sd_tab) %>%
    stringr::str_sub(end = -5)

  sd_tab <- sd_tab %>%
    dplyr::mutate(Stratum = summary_stats2$Stratum) %>%
    tidyr::pivot_longer(-Stratum,
      names_to = "Variable",
      values_to = "sd"
    )

  mean_tab <- summary_stats %>%
    dplyr::select(contains("fn1")) %>%
    tibble::add_row(tibble::tibble_row((population_summary_stats2 %>% dplyr::select(tidyselect::contains("fn1")))))

  names(mean_tab) <- names(mean_tab) %>%
    stringr::str_sub(end = -5)

  mean_tab <- mean_tab %>%
    dplyr::mutate(Stratum = summary_stats2$Stratum) %>%
    tidyr::pivot_longer(-Stratum, names_to = "Variable", values_to = "mn")

  counts_tab <- summary_stats2 %>%
    dplyr::select(Stratum, n)

  heat_data <- dplyr::left_join(mean_tab, sd_tab, by = c("Stratum", "Variable")) %>%
    dplyr::filter(Variable != "rank")

  heat_data <- heat_data %>%
    dplyr::left_join(counts_tab, by = "Stratum")

  temporary_df <- data.frame(
    Variable = unique(heat_data$Variable),
    pop_mean = (heat_data %>% dplyr::filter(Stratum == "Population") %>% dplyr::select(mn)),
    pop_sd = (heat_data %>% dplyr::filter(Stratum == "Population") %>% dplyr::select(sd)),
    pop_n = (heat_data %>% dplyr::filter(Stratum == "Population") %>% dplyr::select(n))
  ) %>%
    dplyr::mutate(
      pop_mean = mn,
      pop_sd = sd,
      pop_n = n
    ) %>%
    dplyr::select(-mn, -sd, -n)

  heat_data <- heat_data %>%
    dplyr::left_join(temporary_df,
      by = "Variable"
    ) %>%
    dplyr::mutate(deviation = dplyr::case_when(
      (mn - pop_mean) / pop_mean >= 0.7 ~ 0.7,
      (mn - pop_mean) / pop_mean <= -0.7 ~ -0.7,
      TRUE ~ (mn - pop_mean) / pop_mean
    ))

  heat_data_simple <- heat_data %>%
    dplyr::select(Stratum, Variable, mn, sd) %>%
    tidyr::pivot_wider(
      names_from = Stratum,
      values_from = c(mn, sd),
      names_glue = "{Stratum}_{.value}"
    )

  heat_data_simple <- heat_data_simple %>%
    dplyr::select(order(colnames(heat_data_simple))) %>%
    dplyr::select(Variable, tidyselect::everything())

  header <- c(1, rep(2, n_strata + 1))
  header_names <- " "

  for (i in 1:n_strata) {
    header_names <- header_names %>%
      append(paste0(
        "Stratum ",
        i,
        "\nn = ",
        counts_tab$n[i]
      ))
  }

  names(header) <- header_names %>%
    append(paste0(
      "Population\n",
      "n = ",
      counts_tab$n %>% tail(n = 1)
    ))

  heat_data_kable <- heat_data_simple %>%
    kableExtra::kbl(
      caption = "Covariate Statistics by Stratum",
      align = "c",
      col.names = c("Variable", rep(c("Mean", "Standard Deviation"), n_strata + 1))
    ) %>%
    kableExtra::kable_styling(c("striped", "hover"),
      fixed_thead = TRUE
    ) %>%
    kableExtra::add_header_above(header)

  stratum_labels <- "Population"

  for (i in 2:(n_strata + 1)) {
    stratum_labels[i] <- paste("Stratum", (i - 1))
  }

  heat_plot <- ggplot(data = heat_data) +
    geom_tile(aes(x = Stratum, y = Variable, fill = deviation), color = "black", width = 0.95) +
    geom_text(aes(x = Stratum, y = ((ncol(summary_stats) + 1) / 2 - 0.15), label = paste(n, "\nunits")), size = 3.4) +
    scale_fill_gradientn(
      name = NULL,
      breaks = c(-0.5, 0, 0.5),
      labels = c(
        "50% \nBelow Mean",
        "Population\nMean",
        "50% \nAbove Mean"
      ),
      colours = c(
        "#990000", "#CC0000",
        "white", "#3D85C6",
        "#0B5294"
      ),
      limits = c(-0.7, 0.7)
    ) +
    scale_x_discrete(
      position = "top",
      expand = c(0, 0),
      labels = c(stratum_labels[-1], "Population")
    ) +
    expand_limits(y = c(0, (ncol(summary_stats) + 1) / 2 + 0.1)) +
    labs(y = NULL, x = NULL) +
    theme(
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 10, colour = "grey15"),
      legend.key.height = unit(1, "cm"),
      legend.text = element_text(size = 10),
      legend.position = "right"
    )

  longer_heat_data <- heat_data %>%
    tidyr::pivot_longer(
      cols = c("mn", "sd"),
      names_to = "mn_or_sd",
      values_to = "values"
    ) %>%
    dplyr::mutate(values = values %>% round(1))

  suppressWarnings(
    heat_plot_final <- heat_plot +
      ggnewscale::new_scale("fill") +
      ggrepel::geom_label_repel(
        data = longer_heat_data %>% dplyr::filter(mn_or_sd == "mn"),
        aes(
          x = Stratum,
          y = Variable,
          label = values,
          fill = mn_or_sd
        ),
        min.segment.length = 1,
        direction = "y",
        nudge_y = 0.05,
        alpha = 0.7,
        size = ifelse((length(levels(heat_data$Variable %>% factor())) + 1) > 7, 2, 3.5)
      ) +
      ggrepel::geom_label_repel(
        data = longer_heat_data %>% dplyr::filter(mn_or_sd == "sd"),
        aes(
          x = Stratum,
          y = Variable,
          label = values,
          fill = mn_or_sd
        ),
        min.segment.length = 1,
        direction = "y",
        nudge_y = -0.05,
        alpha = 0.7,
        size = ifelse((length(levels(heat_data$Variable %>% factor())) + 1) > 7, 2, 3.5)
      ) +
      viridis::scale_fill_viridis(
        labels = c("Mean", "Standard Deviation"),
        begin = 0.4,
        end = 0.8,
        direction = 1,
        discrete = TRUE,
        option = "D"
      ) +
      theme(legend.title = element_blank()) +
      guides(fill = guide_legend(override.aes = aes(label = "")))
  )

  # 5) Save output

  out <- list(
    idvar = idvar,
    variables = var_names,
    dataset = data_name,
    n_strata = n_strata,
    solution = solution,
    pop_data_by_stratum = pop_data_by_stratum,
    summary_stats = summary_stats2,
    data_omitted = data_omitted,
    cont_data_stats = cont_data_stats,
    cat_data_levels = cat_data_levels,
    heat_data = heat_data,
    heat_data_simple = heat_data_simple,
    heat_data_kable = heat_data_kable,
    heat_plot = heat_plot_final
  )

  class(out) <- "generalizeR_stratify"

  return(invisible(out))
}

#' Print method for "generalizeR_stratify" class
#'
#' @param x An object of class "generalizeR_stratify"
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.generalizeR_stratify
#' @export

print.generalizeR_stratify <- function(x, ...) {

  cat("\nA generalizeR_stratify object: \n\n")

  cat(" - Dataset used:", crayon::bold(x$dataset), "\n\n")

  cat(" - Stratification variables:", paste0(crayon::bold$blue(x$variables), collapse = ", "))

  cat("\n\n - Number of strata:", crayon::bold(x$n_strata))

  cat("\n\n - Number of observations dropped due to missing data:",
      crayon::bold(nrow(x$data_omitted)),
      "(see $data_omitted for dropped observations) \n\n")
}

#' Summary method for "generalizeR_stratify" class
#'
#' @param x An object of class "generalizeR_stratify"
#' @param ... Other arguments passed to or from other methods
#'
#' @export summary.generalizeR_stratify
#' @export

summary.generalizeR_stratify <- function(object, ...) {

  out <- object

  class(out) <- "summary.generalizeR_stratify"

  return(out)
}

#' Print method for "summary.generalizeR_stratify" class
#'
#' @param x An object of class "generalizeR_stratify"
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.summary.generalizeR_stratify
#' @export

print.summary.generalizeR_stratify <- function(x, ...) {

  cat(paste0(rep("=", 80), collapse = ""))

  cat(paste0("\nSummary of stratification performed with '", x$dataset, "' dataset:", "\n\n"))

  cat("Stratification variables:", paste0(crayon::bold$blue(x$variables), collapse = ", "))

  cat(paste0("\nObservations dropped due to missing data: ", crayon::bold(nrow(x$data_omitted)), " (see $data_omitted)\n"))

  cat(paste0("Population size: ", crayon::bold(nrow(x$pop_data_by_stratum)), "\n"))

  cat(paste0("Number of strata specified: ", crayon::bold(x$n_strata), "\n"))

  cat(paste0("Proportion of variation in population explained by strata: "))

  cat(crayon::bold(paste(100 * round(x$solution$between.SS_DIV_total.SS, 4), "%", sep = "")))

  cat("\n")

  cat(paste0(rep("=", 80), collapse = ""))

  cat("\nCovariate Statistics by Stratum: \n\n")

  x$heat_data_simple %>%
    as.data.frame() %>%
    print()

  x$heat_data_kable %>%
    print()

  suppressWarnings(
    tryCatch(
      {
        x$heat_plot %>%
          print()
      },
      error = function(cond) {
        message("Your Plots pane is too small for the heat map to be displayed. If you still want to \nview the plot, try resizing the pane and then running 'x$heat_plot' where 'x' is the \nname assigned to your 'generalizeR_stratify' object.")
        return(NA)
      }
    )
  )
  invisible(x)
}

#' Internal function that provides the guided version of \code{stratify()}. Walks users through the provision and visualization of parameters needed to stratify observations, including id variables, covariates of interest, and the number of strata
#'
#' Intended only to be called within \code{stratify}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 3
#'
#' @param data data.frame object containing the population data to be stratified (rows are observations).
#'
#' @md


.stratify.guided <- function(data,
                             verbose) {

  # 1) Introduction -----------------------------------------------------------

  cat(crayon::bold("\nWelcome to stratify()! \n\n"))

  cat("Given a sampling frame, this function can tell you how to stratify your population of \ninterest so as to draw samples from it with maximal generalizability.\n\n")

  cat("You have chosen the '",
    data_name,
    "' dataframe to represent your population of interest.\n\n",
    sep = ""
  )

  cat("To store your results, make sure you have assigned this function to an object.\n\n")

  cat("To exit this function, press <Esc>. \n\n")

  # 2) Selection of id variable -----------------------------------------------

  repeat {
    idvar <- readline(prompt = "Enter the name of the ID Variable in your dataframe: ")

    ## Check ##
    if (idvar %in% names(data)) {

      break
    }

    cat(crayon::red("\nERROR: We could not find that variable. Your ID variable must be one of the variables \nin the '"),
        crayon::red(data_name),
        crayon::red("' dataframe.\n\n"),
        sep = ""
    )
  }

  id <- data %>%
    dplyr::select(tidyselect::all_of(idvar))

  # 3) Selection of variables for stratification ------------------------------

  data_subset <- data %>%
    dplyr::select(-tidyselect::all_of(idvar))

  data_subset %>%
    .make.var.overview()

  cat("\nIn the Viewer pane to the right you will find a table that displays each variable in \nthe '",
    data_name,
    "' dataframe along with its data type and number of levels (only \nrelevant for factor variables).\n\n",
    sep = ""
  )

  cat(crayon::yellow$bold("Please note that any character variables that may have been present in your \ndataframe have been automatically converted to factor variables.\n\n"))

  readline(prompt = "Press <Return> to continue.")

  cat("\nYou are now ready to select your stratification variables. The following are the \nvariables available in your dataset. Which key variables do you think may explain \nvariation in your treatment effect? Typically, studies include 4-6 variables for \nstratification.\n\n")

  cat(crayon::yellow$bold("You must choose at least 2 variables and you may not choose any factor variables \nwith more than 4 levels.\n"))

  variables <- .select.list.stratify(data_subset)

  data_subset <- data_subset %>%
    dplyr::select(tidyselect::all_of(variables))

  # 4) Let user know about missing observations -------------------------------

  missing_obs_table <- data.frame(id, data_subset) %>%
    sapply(function(x) sum(is.na(x))) %>%
    as.data.frame() %>%
    `colnames<-`("n_missing") %>%
    tibble::rownames_to_column("variable")

  n_missing <- data.frame(id, data_subset) %>%
    dplyr::filter(dplyr::if_any(tidyselect::everything(), is.na)) %>%
    nrow()

  data_subset <- data.frame(id, data_subset) %>%
    tidyr::drop_na() %>%
    dplyr::select(-tidyselect::all_of(idvar))

  cat("\nThe following table shows how many observations are missing for the variables you have \nchosen.\n\n")

  missing_obs_table %>%
    kableExtra::kbl(
      caption = "Missing Observations by Variable",
      align = "l",
      col.names = c("Variable", "Number Missing")
    ) %>%
    kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    print()

  missing_obs_table %>%
    print(row.names = FALSE)

  cat("\nNote that there may be multiple missing observations per row, meaning that the total \nnumber of dropped observations won't necessarily be the sum of the number of missing \nobservations in each variable.\n\n")

  cat(crayon::yellow$bold("In total,", n_missing, "rows will be dropped from the inference population before \nstratification due to missing observations.\n\n"))

  readline(prompt = "Hit <Return> to proceed once you have viewed the missing observations.")

  # 5) Overview of categorical variables --------------------------------------

  cat_data <- data_subset %>%
    dplyr::select_if(is.factor)

  cat_data_vars <- names(cat_data)

  if (dim(cat_data)[2] >= 1L) {
    cat_data_plot <- data.frame(cat_data)

    cat("\nPlease review the descriptive statistics of your categorical variables (factors). ",
      "Bar \ncharts and tables for each variable will also be printed in the Plots and Viewer \npanes to the right.\n",
      sep = ""
    )

    n_cat_vars <- ncol(cat_data_plot)

    fill_colors_cat <- viridis::plasma(n_cat_vars,
      alpha = 0.7,
      direction = sample(c(-1, 1),
        size = 1
      )
    ) %>%
      sample()

    outline_colors_cat <- viridis::turbo(n_cat_vars) %>%
      sample()

    for (i in 1:n_cat_vars) {
      var_name <- cat_data_vars[i]

      levels(cat_data_plot[[var_name]]) <- stringr::str_wrap(levels(cat_data_plot[[var_name]]), width = 10)

      barfig <- cat_data_plot %>%
        dplyr::group_by(dplyr::across(tidyselect::all_of(var_name))) %>%
        dplyr::summarise(count = dplyr::n()) %>%
        dplyr::mutate(ordered_factor = forcats::fct_reorder(.[[var_name]], count)) %>%
        ggplot2::ggplot(aes(x = ordered_factor, y = count)) +
        geom_col(
          fill = fill_colors_cat[i],
          color = outline_colors_cat[i]
        ) +
        theme_minimal() +
        xlab(var_name) +
        labs(title = paste("Bar Chart of", var_name))

      cat("\n")

      print(barfig)

      par(ask = TRUE)

      cat("\nNumber of Observations in Levels of Factor ",
        paste(crayon::blue$bold(var_name)),
        ":\n\n",
        sep = ""
      )

      cat_data_table <- cat_data_plot[, i] %>%
        table() %>%
        data.frame() %>%
        `colnames<-`(c("Level", "Frequency")) %>%
        print(row.names = FALSE)

      cat_data_table %>%
        kableExtra::kbl(
          col.names = c("Level", "Frequency"),
          caption = paste("Number of Observations in Levels of Factor ", var_name),
          align = "l"
        ) %>%
        kableExtra::kable_styling(c("striped", "hover")) %>%
        print()
    }

    cat("\n")

    readline(prompt = "Hit <Return> to proceed once you have reviewed the categorical variables.")
  }

  # 6) Overview of continuous variables ---------------------------------------

  cont_data <- data_subset %>%
    dplyr::select_if(purrr::negate(is.factor))

  cont_data_vars <- cont_data %>%
    names()

  if (dim(cont_data)[2] >= 1L) {
    cat("\nPlease review the descriptive statistics of your continuous variables. Histograms and \ntables for each variable will also be printed in the Plots and Viewer panes to the \nright. \n\n")

    n_cont_vars <- ncol(cont_data)

    fill_colors_cont <- viridis::viridis(n_cont_vars,
      alpha = 0.7,
      direction = sample(c(-1, 1),
        size = 1
      )
    ) %>%
      sample()

    for (i in 1:n_cont_vars) {
      hist <- cont_data %>%
        ggplot2::ggplot(aes(x = cont_data[, i])) +
        geom_histogram(
          bins = 30,
          fill = fill_colors_cont[i],
          color = "black"
        ) +
        theme_minimal() +
        xlab(cont_data_vars[i]) +
        labs(title = paste("Histogram of", cont_data_vars[i]))

      print(hist)

      par(ask = TRUE)
    }

    readline(prompt = "Hit <Return> to view a table of summary statistics for the continuous variables:\n")

    sumstats <- cont_data %>%
      .make.cont.data.tbl()

    sumstats %>%
      print(row.names = FALSE)

    sumstats %>%
      kableExtra::kbl(col.names = c("Variable", "Min", "Pct50", "Max", "Mean", "SD")) %>%
      kableExtra::kable_styling(c("striped", "hover")) %>%
      print()
  }

  par(ask = FALSE)

  cat("\n")

  readline(prompt = "Hit <Return> to proceed.")

  # 7) Selection of number of strata ------------------------------------------

  cat("\nStratification will help you develop a recruitment plan so that your study will result \nin an unbiased estimate of the ", crayon::bold("average treatment effect (ATE)"), ". Without using strata, \nit is easy to end up with a sample that is very different from your inference \npopulation. \n\nGeneralization works best when strata are ", crayon::bold("homogeneous"), ". That means units within each \nstratum are almost identical in terms of relevant variables.\n\n", sep = "")

  cat("Enter the number of strata in which you wish to divide your population. Typically, ",
    crayon::bold("\nthe more strata"),
    ", ",
    crayon::bold("the better"),
    "; with fewer strata, units in each stratum are no longer \nidentical. However, increasing ",
    "the number of strata uses more resources, because you \nmust sample a given number of units ",
    "from each stratum. Choosing 4-6 strata is common. \n\nTry a few numbers and choose the 'best' one for you.\n\n",
    sep = ""
  )

  n_strata_correct <- FALSE

  try_again <- TRUE

  while (try_again) {

    while (!n_strata_correct) {
      n_strata <- suppressWarnings(as.numeric(readline(prompt = "Number of strata: ")))

      if (is.na(n_strata) || !is.numeric(n_strata) || n_strata <= 1 || n_strata %% 1 != 0 || n_strata >= (nrow(data_subset)) - 1) {
        cat(crayon::red("\nERROR: The number of strata must be a single integer greater than 1 and at most 2 less \nthan the the number of observations.\n\n"))

        next
      }

      n_strata_correct <- TRUE
    }

    # 9) Pass arguments to .stratify.calculate() where actual stratification is performed

    output <- .stratify.calculate(
      data = data,
      n_strata = n_strata,
      variables = variables,
      idvar = idvar,
      verbose = verbose
    )

    cat(crayon::blue$bold("Congratulations, you have successfully grouped your data into", n_strata, "strata!\n\n"))

    cat("You can pull up the results at any time by passing your 'generalizeR_stratify' \nobject into summary().\n\n")

    readline(prompt = "Hit <Return> to view the results.")

    output %>%
      summary() %>%
      print()

    try <- "1"

    while (T) {
      cat("\n\nDo you want to try a different number of strata?\n1: Yes\n2: No\n\n")

      try <- readline(prompt = "Enter your selection: ")

      if (try == "1" | try == "2") {
        break
      }
      cat("Invalid Input. Please enter 1 or 2\n\n")
    }

    if (try == "1") {
      try_again <- TRUE
      n_strata_correct <- FALSE
    } else {
      try_again <- FALSE
    }

  }

  # 8) Return output

  return(invisible(output))
}


#' Internal function that provides the unguided version of \code{stratify()}. Given the dataset of interest, the number of strata, covriates of interest, and an id variable, performs stratification
#'
#' Intended only to be called within \code{stratify}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 3
#'
#' @param data data.frame object containing the population data to be stratified (rows are observations).
#' @param n_strata integer, a number of strata in which to divide to cluster population
#' @param variables character, provide a character vector of the names of stratifying variables (from population data frame)
#' @param idvar character, provide a character vector of the name of the ID variable (from population data frame)
#' @param verbose logical, defaults to TRUE.
#' @return The function returns a list of class "generalizeR_stratify" that can be provided as input to \code{recruit()}. More information on the components of this list can be found above under "Details."
#' @details The list contains 14 components: \code{idvar}, \code{variables}, \code{dataset}, \code{n_strata}, \code{solution}, \code{pop_data_by_stratum}, \code{summary_stats}, \code{data_omitted}, \code{cont_data_stats}, \code{cat_data_levels}, \code{heat_data}, \code{heat_data_simple}, \code{heat_data_kable}, and \code{heat_plot}.
#'
#' @md

.stratify.unguided <- function(data,
                               n_strata = NULL,
                               variables = NULL,
                               idvar = NULL,
                               verbose = TRUE) {

  # Check user arguments for invalid specifications -------------------------

  # n_strata
  is_n_strata_valid <- function(n_strata, data, variables, idvar) {

   (n_strata %% 1 == 0) && (n_strata > 1) && .n.strata.less.than.max(n_strata, data, variables, idvar)
  }

  assertthat::on_failure(is_n_strata_valid) <- function(call, env) {

    "You must pass a single integer greater than 1 and at least 2 less than the number of non-missing observations to the 'n_strata' argument if you are running the unguided version of this function."
  }

  assertthat::assert_that(is_n_strata_valid(n_strata, data, variables, idvar))

  # variables
  are_variables_valid <- function(variables) {

    is.character(variables) && length(variables) > 1
  }

  assertthat::on_failure(are_variables_valid) <- function(call, env) {

    "You must pass a character vector consisting of two or more column names from your dataframe to the 'variables' argument if you are running the unguided version of this function."
  }

  assertthat::assert_that(are_variables_valid(variables))

  # idvar
  is_idvar_valid <- function(idvar) {

    is.character(idvar)
  }

  assertthat::on_failure(is_idvar_valid) <- function(call, env) {

    "You must pass one of the column names from your dataframe as a character string to the 'idvar' argument if you are running the unguided version of this function."
  }

  assertthat::assert_that(is_idvar_valid(idvar))

  # Check for invalid variable names ----------------------------------------
  invalid_vars <- c(variables, idvar) %>%
    setdiff(names(data))

  assertthat::on_failure(is_empty) <- function(call, env) {

    paste(
      "The following variables are not columns in the dataframe you specified:\n",
      paste(crayon::blue$bold(invalid_vars),
            collapse = ", "
      )
    )
  }

  assertthat::assert_that(is_empty(invalid_vars))

  # Check if there are any factor variables with more than 4 levels ---------

  invalid_factors <- data %>%
    dplyr::select(tidyselect::all_of(variables)) %>%
    .check.factor.levels()


  assertthat::on_failure(is_empty) <- function(call, env) {

    paste0(
      "This function will not allow a factor variable to have more than 4 levels.\n",
      "The following factor variables have more than 4 levels:\n",
      paste(crayon::blue$bold(invalid_factors), collapse = ", "),
      "\nPlease re-code your desired factor levels from these variables as dummy variables (see \nthe package 'fastDummies')."
    )
  }

  assertthat::assert_that(is_empty(invalid_factors))

  # Verify that there are no variables with all obs. missing --------
  na_variables <- data %>%
    dplyr::select(tidyselect::all_of(variables)) %>%
    sapply(function(x) sum(!is.na(x))) %>%
    data.frame() %>%
    dplyr::filter(.==0) %>%
    row.names()

  assertthat::on_failure(is_empty) <- function(call, env) {

    paste0(
      crayon::red("This function will not allow a variable with all observations missing.\n\n"),
      crayon::red("All observations are missing for the following variables:\n\n"),
      paste(crayon::blue$bold(na_variables), collapse = ", "),
      crayon::red("\n\nPlease choose a different set of variables.\n"),
      sep = ""
    )
  }

  assertthat::assert_that(is_empty(na_variables))

  data_subset <- data %>%
    dplyr::select(tidyselect::all_of(variables))

  list_single <- c()

  for (name in colnames(data_subset)) {
    if (length(unique(data_subset[[name]])) == 1) {
      list_single <- append(list_single, name)
    }
  }


  assertthat::on_failure(is_empty) <- function(call, env) {

    paste0(
      crayon::red("The following variables have only one value:\n\n"),
      paste(crayon::blue$bold(list_single), collapse = ", "),
      crayon::red("\n\nPlease choose covariates with more than one value\n\n"),
      sep = ""
        )
  }

  assertthat::assert_that(is_empty(list_single))

  output <- .stratify.calculate(
    data = data,
    n_strata = n_strata,
    variables = variables,
    idvar = idvar,
    verbose = verbose
  )

  return(invisible(output))
}

.select.list.stratify <- function(data) {

  choices <- names(data)

  num_choices <- length(choices)

  options <- paste0(format(1:num_choices),
                    ":  ",
                    choices)

  repeat {

    if (num_choices > 10L) {

      formatted_options <- format(options)
      nw <- nchar(formatted_options[1L], "w") + 2L
      ncol <- getOption("width") %/% nw

      if (ncol > 1L) {

        options <- paste0(formatted_options,
                          c(rep.int("  ", ncol - 1L), "\n"),
                          collapse = "")
      }

      cat("", options, sep = "\n")
    } else {

      cat("", options, "", sep = "\n")
    }

    cat("Type two or more numbers separated by spaces and then hit <Return> to continue. \n\n")

    selection <- tryCatch(

      scan("",
           what = 0,
           quiet = TRUE,
           nlines = 1
      ),
      error = identity
    )

    # Verify that user's selection did not throw an error in tryCatch()
    if(!inherits(selection, "error")) {

      # Verify that user chose at least two variables
      if (length(selection) >= 2L) {

        # Verify that user only input integers between 1 and the number of variables in the dataframe
        if (all(selection %in% 1:num_choices)) {

          selection <- selection %>%
            unique() %>%
            sort()

          variables <- choices[selection]

          invalid_factors <- data %>%
            dplyr::select(tidyselect::all_of(variables)) %>%
            .check.factor.levels()

          # Verify that user didn't choose any factor variables with more than 4 levels
          if (rlang::is_empty(invalid_factors)) {

            na_variables <- data %>%
              dplyr::select(tidyselect::all_of(variables)) %>%
              sapply(function(x) sum(!is.na(x))) %>%
              data.frame() %>%
              dplyr::filter(. == 0) %>%
              row.names()

            # Verify that user didn't choose any variables with all observations missing
            if (rlang::is_empty(na_variables)) {

              data_subset <- data %>%
                dplyr::select(tidyselect::all_of(variables))

              #lin_dep <- plm::detect.lindep(data_subset)
              lin_dep <- NULL

              # Verify that user didn't choose any linearly dependent variables - currently not checking for this due to difficulty finding a linear matrix
              if (is.null(lin_dep)) {

                num_single <- 0
                list_single <- c()

                for (name in colnames(data_subset)) {
                  if (length(unique(data_subset[[name]])) == 1) {
                    num_single <- num_single + 1
                    list_single <- append(list_single, name)
                  }
                }

                # Verify that users know if they have columns with one value
                if (num_single == 0) {

                  cat("\nYou have selected the following stratifying variables:\n\n",
                      paste(crayon::blue$bold(variables), collapse = ", "),
                      "\n\n",
                      sep = "")

                  data %>%
                    dplyr::select(tidyselect::all_of(variables)) %>%
                    .make.var.overview(print_to_console = TRUE)

                  # Verify that user is happy with their selection
                  if (utils::menu(choices = c("Yes", "No"), title = cat("\nIs this correct?")) == 1) {

                    return(choices[selection])
                  }



                cat("\n")

                data %>%
                  .make.var.overview()
                next

                }

                cat(crayon::red("ERROR: The variable(s),"),
                    paste(crayon::blue$bold(list_single), collapse = ", "),
                    crayon::red("have only one value. Please make sure none of your chosen covariates have only one value\n\n"))

                next
              }

              # will never run as of now
              cat(crayon::red("\nERROR: The variables:"))
              for (name in names(lin_dep)) {
                cat(names)
                cat(" ")
              }
              cat("are linearly dependent. Please make sure none of your chosen covariates are independent \n\n")

              next
            }

            cat(crayon::red("\nERROR: Invalid selection. This function will not allow a selected variable to have only missing observations.\n\n"),
                crayon::red("All observations are missing for the following variables you have chosen:\n\n"),
                paste(crayon::blue$bold(na_variables), collapse = ", "),
                sep = ""
            )

            next
          }

          cat(crayon::red("\nERROR: Invalid selection. This function will not allow a factor variable to have more \nthan 4 levels.\n\n"),
              crayon::red("The following factor variables have more than 4 levels:\n\n"),
              paste(crayon::blue$bold(invalid_factors), collapse = ", "),
              crayon::red("\n\nIf you still wish to use these variables to stratify your population, please exit out of \nstratify() by hitting <Esc> and re-code your desired factor levels from these variables \nas dummy variables (see the package 'fastDummies').\n\n"),
              sep = "")

          next
        }

        cat(paste0(crayon::red("\nERROR: Invalid selection. Each input must be a single integer between 1 and "),
                   crayon::red(num_choices),
                   crayon::red(".\n\n")))
        next
      }

      cat(crayon::red("\nERROR: Invalid selection. You must select at least 2 stratification variables.\n\n"))
      next
    }

    cat(crayon::red("\nERROR: Invalid selection. Please try again.\n\n"))
  }
}

#' Internal function that checks whether any column in a dataset is a factor variable with more than a given number of levels
#'
#' Intended only to be called within \code{stratify}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 3
#'
#' @param data data.frame, a dataframe of factor variables of interest
#' @param maxlevels integer, a maximum permissible number of levels
#' @return invalid_factors, list of variable names with more than the permitted number of levels
#'
#' @md
#'
.check.factor.levels <- function(data,
                                 maxlevels = 4L) {
  invalid_factors <- data %>%
    dplyr::select_if(is.factor) %>%
    sapply(nlevels) %>%
    {. > maxlevels} %>%
    which() %>%
    names()

  return(invalid_factors)
}


#' Internal function that checks whether a list provided is empty (avoiding conflict w/ check on >4 level factor categorical vars)
#'
#' Intended only to be called within \code{stratify}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 3
#'
#' @param na_col_list list, a list of char corresponding to dataframe columns where all rows are NA
#' @return logical, whether the inputted list is empty
#'
#' @md

.check.no.na.cols <- function(na_col_list){
  is_empty(na_col_list)
}

#' Internal function that creates and formats a table of the variables in the dataset provided
#'
#' Intended only to be called within \code{stratify}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 3
#'
#' @param dataset data.frame, a dataframe of continous variables
#' @param print_to_console logical, whether the variable summary table should be printed to the console, defaults to false
#' @return None
#'
#' @md

.make.var.overview <- function(dataset,
                               print_to_console = FALSE) {
  var_overview <- dataset %>%
    purrr::map_df(function(x) {
      tibble::tibble(
        Type = class(x),
        Levels = nlevels(x)
      )
    }) %>%
    dplyr::mutate(Variable = names(dataset)) %>%
    dplyr::select(Variable, everything()) %>%
    data.frame()

  var_overview %>%
    kableExtra::kbl(
      caption = "Variable Overview",
      align = "l",
      row.names = TRUE
    ) %>%
    kableExtra::kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    print()

  if (print_to_console == TRUE) {
    print(var_overview, row.names = FALSE)
  }
}

#' Internal function that creates and formats a summary table of the continuous variables within the provided dataset
#'
#' Intended only to be called within \code{stratify}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 3
#'
#' @param cont_data data.frame, a dataframe of continous variables
#' @return The function returns a dataframe with summary statistics of the inputted data
#'
#' @md

.make.cont.data.tbl <- function(cont_data) {
  cont_data_tbl <- cont_data %>%
    purrr::map_df(
      function(x) {
        tibble::tibble(
          min = min(x),
          pct50 = median(x),
          max = max(x),
          mean = mean(x),
          sd = sd(x)
        )
      }
    ) %>%
    dplyr::mutate_all(round, digits = 3) %>%
    dplyr::mutate(variable = names(cont_data)) %>%
    dplyr::select(variable, tidyselect::everything()) %>%
    data.frame() %>%
    janitor::clean_names()

  return(cont_data_tbl)
}

#' Internal function that ensures that the number of strata requested does not exceed the maximum possible (maximum n-2 where n is number of observations w/ no missing variables of interest)
#'
#' Intended only to be called within \code{stratify}, not as a standalone function
#'
#' @keywords internal
#'
#' @order 3
#'
#' @param n_strata integer, a number of strata in which to divide to cluster population
#' @param data data.frame, a data.frame object with the population of interest (rows are observations)
#' @param variables character, provide a character vector of the names of stratifying variables (from population data frame)
#' @param idvar character, provide a character vector of the name of the ID variable (from population data frame)
#' @return The function returns a boolean value indicating whether the number of observations with non-missing variables of interests exceeds 1 plus the number of strata requested
#'
#' @md

.n.strata.less.than.max <- function(n_strata, data_interest, variables, idvar){
  data_interest %>%
    dplyr::select(tidyselect::all_of(variables), tidyselect::all_of(idvar)) %>%
    tidyr::drop_na() %>%
    nrow() - 1 > n_strata
}
