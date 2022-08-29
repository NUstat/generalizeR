
make_var_overview = function(dataset, print_to_console = FALSE) {

  var_overview = dataset %>%
    map_df(function(x) {
      tibble(
        Type = class(x),
        Levels = nlevels(x))
    }) %>%
    mutate(Variable = names(dataset)) %>%
    select(Variable, everything()) %>%
    data.frame()

  var_overview %>%
    kbl(caption = "Variable Overview",
        align = "l",
        row.names = TRUE) %>%
    kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    print()

  if(print_to_console == TRUE){
    print(var_overview, row.names = FALSE)
  }
}
