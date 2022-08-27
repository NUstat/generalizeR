#' Stratify a Population Data Frame
#'
#' The function \code{stratify()} takes as input any data frame that you want to stratify into clusters. Typically, the goal of such stratification is sampling for generalizability. This function, and the others in this package, are designed to mimic the website https://www.thegeneralizer.org/.
#'
#' @param data data.frame object containing the population data to be stratified.
#' @param guided logical, defaults to TRUE. Whether the function should be guided (ask questions and behave interactively throughout) or not. If set to FALSE, the user must provide values for other arguments below
#' @param n_strata integer, defaults to NULL. If guided is set to FALSE, must provide a number of strata in which to divide to cluster population
#' @param variables character, defaults to NULL. If guided is set to FALSE, must provide a character vector of the names of stratifying variables (from population data frame)
#' @param idvar integer, defaults to NULL. If guided is set to FALSE, must provide a character vector of the name of the ID variable (from population data frame)
#' @param verbose logical, defaults to TRUE.
#' @return The function returns a list of class "generalizer_output" that can be provided as input to \code{recruit()}. More information on the components of this list can be found above under "Details."
#' @details The list contains 14 components: \code{idvar}, \code{variables}, \code{dataset}, \code{n_strata}, \code{solution}, \code{pop_data_by_stratum}, \code{summary_stats}, \code{data_omitted}, \code{cont_data_stats}, \code{cat_data_levels}, \code{heat_data}, \code{heat_data_simple}, \code{heat_data_kable}, and \code{heat_plot}.
#'
#' \itemize{
#' \item{\code{pop_data_by_stratum}: }{a tibble with number of rows equal to the number of rows in the inference population (\code{data}) and number of columns equal to the number of stratifying variables (dummy-coded if applicable) plus the ID column (\code{idvar}) and a column representing stratum membership, \code{Stratum}}
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
#'    variables = c("pct_female", "pct_white"), idvar= "unitid")
#' }
#' }
#' @md

stratify = function(data = NULL,
                    guided = TRUE,
                    n_strata = NULL,
                    variables = NULL,
                    idvar = NULL,
                    verbose = TRUE)  {

  # Check whether 'data' argument has been specified
  if(is.null(data)) {

    stop(simpleError("argument 'data' is missing. You must specify the dataframe containing the \ninference population you wish to stratify."))
  }

  # Check whether data object is of type 'data.frame'
  if(!is.data.frame(data)) {

    stop(simpleError("Your data object must be of type 'data.frame'."))
  }

  # Store name of data as global variable so it can be accessed by stratify_basic(). Must be done before function argument 'data' is evaluated for the first time.
  data_name <<- data %>%
    lazyeval::expr_text()

  # Immediately convert all character variables to factors
  data = data %>%
    mutate(across(where(is_character), as_factor))

  # Begin guided version ------------------------------------------------------
  if(guided == TRUE) {

    source("R/stratify/stratify_guided.R")

    inputs = stratify_guided(data)

    n_strata = inputs$n_strata

    variables = inputs$variables

    idvar = inputs$idvar
  }

  # Begin unguided version ----------------------------------------------------
  else {

    source("R/stratify/stratify_unguided.R")

    stratify_unguided(data = data,
                      n_strata = n_strata,
                      variables = variables,
                      idvar = idvar)
  }

  # Pass arguments to stratify_basic() where actual stratification is performed

  source("R/stratify/stratify_basic.R")

  output = stratify_basic(data = data,
                          n_strata = n_strata,
                          variables = variables,
                          idvar = idvar,
                          verbose = verbose)

  rm(list = deparse(substitute(data_name)), envir=.GlobalEnv) # delete object data_name from global environment

  # Final message for guided version ------------------------------------------
  if(guided == TRUE) {

    cat(crayon::blue$bold("Congratulations, you have successfully grouped your data into", n_strata, "strata!\n\n"))

    cat("You can pull up the results at any time by passing your generalizer_stratify \nobject into summary().\n\n")

    readline(prompt = "Hit <Return> to view the results.")

    output %>%
      summary() %>%
      print()
  }

  class(output) = "generalizer_stratify"

  return(invisible(output))
}



