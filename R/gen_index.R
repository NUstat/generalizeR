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

gen_index <- function(trial_ps, pop_ps) {
  ##Baklizi and Eidous (2006) estimator
  # bandwidth
  h = function(x){
    n = length(x)
    optim_binwidth = (4*sqrt(var(x))^5/(3*n))^(1/5)
    if(is.na(optim_binwidth) | is.nan(optim_binwidth)){
      optim_binwidth = 0
    }
    if(optim_binwidth < 0.001){ # this yielded a b index of 0.9999501 for (at least one specific case of) "perfectly" stratified data
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

  return(as.numeric(integrate(function(x) sqrt(kg(x, dat1B)*kg(x, dat2B)), -Inf, Inf)$value))
}
