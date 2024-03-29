% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stratify.R
\name{stratify}
\alias{stratify}
\title{Stratify a Population Data Frame}
\usage{
stratify(
  data = NULL,
  guided = TRUE,
  n_strata = NULL,
  variables = NULL,
  idvar = NULL,
  verbose = TRUE
)
}
\arguments{
\item{data}{data.frame object containing the population data to be stratified (observations as rows); must include a unique id variable for each observation, as well as covariates.}

\item{guided}{logical, defaults to TRUE. Whether the function should be guided (ask questions and behave interactively throughout) or not. If set to FALSE, the user must provide values for other arguments below}

\item{n_strata}{integer, defaults to NULL. If guided is set to FALSE, must provide a number of strata in which to divide to cluster population}

\item{variables}{character, defaults to NULL. If guided is set to FALSE, must provide a character vector of the names of stratifying variables (from population data frame)}

\item{idvar}{character, defaults to NULL. If guided is set to FALSE, must provide a character vector of the name of the ID variable (from population data frame)}

\item{verbose}{logical, defaults to TRUE.}
}
\value{
The function returns a list of class "generalizeR_stratify" that can be provided as input to \code{recruit()}. More information on the components of this list can be found above under "Details."
}
\description{
The function \code{stratify()} takes as input any data frame with observations (rows) that you wish to stratify into clusters. Typically, the goal of such stratification is developing a sampling design for maximizing generalizability. This function, and the others in this package, are designed to mimic the website https://www.thegeneralizer.org/.
}
\details{
The list contains 14 components: \code{idvar}, \code{variables}, \code{dataset}, \code{n_strata}, \code{solution}, \code{pop_data_by_stratum}, \code{summary_stats}, \code{data_omitted}, \code{cont_data_stats}, \code{cat_data_levels}, \code{heat_data}, \code{heat_data_simple}, \code{heat_data_kable}, and \code{heat_plot}.

\describe{
\item{\code{pop_data_by_stratum}: }{a tibble with number of rows equal to the number of rows in the inference population (\code{data}) and number of columns equal to the number of stratifying variables (dummy-coded if applicable) plus the ID column (\code{idvar}) and a column representing stratum membership, \code{Stratum}}
}
}
\examples{
library(tidyverse)

selection_covariates <- c("total", "pct_black_or_african_american",
                          "pct_white", "pct_female", "pct_free_and_reduced_lunch")
stratify(generalizeR:::inference_pop, guided = FALSE, n_strata = 4,
         variables = selection_covariates, idvar= "ncessch")

}
\references{
Tipton, E. (2014). Stratified sampling using cluster analysis: A sample selection strategy for improved generalizations from experiments. \emph{Evaluation Review}, \emph{37}(2), 109-139.

Tipton, E. (2014). How generalizable is your experiment? An index for comparing experimental samples and populations. \emph{Journal of Educational and Behavioral Statistics}, \emph{39}(6), 478-501.
}
