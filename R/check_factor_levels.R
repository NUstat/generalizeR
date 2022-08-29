
check_factor_levels = function(data, maxlevels = 4L) {

  invalid_factors = data %>%
    select_if(is.factor) %>%
    sapply(nlevels) %>%
    `>` (maxlevels) %>%
    which() %>%
    names()

  return(invalid_factors)
}
