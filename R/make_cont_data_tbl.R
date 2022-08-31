
make_cont_data_tbl = function(cont_data) {

  cont_data_tbl = cont_data %>%
    map_df(

      function(x) {

        tibble(
          min = min(x),
          pct50 = median(x),
          max = max(x),
          mean = mean(x),
          sd = sd(x)
        )
      }

    ) %>%
    mutate_all(round, digits = 3) %>%
    mutate(variable = names(cont_data)) %>%
    select(variable, everything()) %>%
    data.frame() %>%
    janitor::clean_names()

  return(cont_data_tbl)
}
