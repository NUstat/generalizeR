
stratify_basic = function(data,
                          n_strata = NULL,
                          variables = NULL,
                          idvar = NULL,
                          verbose = TRUE) {

  #source("R/stratify/helper functions/make_cont_data_tbl.R")

  # 1) Store all information given by user

  # Extract id variable
  id = data %>%
    select(all_of(variables), all_of(idvar)) %>%
    drop_na() %>%
    select(all_of(idvar))

  # Create a dataframe that includes only the variables selected by the user
  data_subset = data %>%
    select(all_of(variables))

  # Save all rows with missing observations in a new dataframe
  data_omitted = data_subset %>%
    filter(if_any(everything(), is.na))

  # Drop rows with missing observations in the subsetted dataframe
  data_subset = data_subset %>%
    drop_na()

  # Store categorical and continuous variables in separate dataframes
  cat_data = data_subset %>%
    select_if(is.factor)

  cont_data = data_subset %>%
    select_if(negate(is.factor))

  # Make dummy variables out of factors
  if(dim(cat_data)[2] >= 1L) {

    cat_data_with_dummies = cat_data %>%
      fastDummies::dummy_cols(remove_first_dummy = TRUE) %>%
      select_if(negate(is.factor))

    data_full = cont_data %>%
      cbind(cat_data_with_dummies) %>%
      clean_names()
  }

  else {

    data_full = cont_data %>%
      cbind(cat_data) %>%
      clean_names()
  }

  var_names = data_full %>% names()

  # 2) Save summaries of variables in dataset
  cont_data_stats = cont_data %>%
    make_cont_data_tbl()

  cat_data_levels = cat_data %>%
    map_df(

      function(x) {

        tibble(
          Type = class(x),
          Levels = nlevels(x))
        }

      ) %>%
    mutate(Variable = names(cat_data)) %>%
    select(Variable, everything()) %>%
    data.frame()

  # 3) Perform stratification - cluster with KMeans

  if(verbose == TRUE) {

    cat("\nThis might take a little while. Please bear with us.")

    distance = daisy(data_full, metric = "gower") %>%
      suppressWarnings()

    cat("\n\nCalculated distance matrix.\n")

    solution = KMeans_rcpp(as.matrix(distance),
                           clusters = n_strata,
                           verbose = TRUE)
  }

  else {

    distance = daisy(data_full, metric = "gower") %>%
      suppressWarnings()

    solution = KMeans_rcpp(as.matrix(distance),
                           clusters = n_strata,
                           verbose = FALSE)
  }

  # 4) Create various datasets

  pop_data_by_stratum = data.frame(id,
                                   data_full,
                                   Stratum = solution$clusters) %>%
    select(Stratum, everything()) %>%
    arrange(Stratum) %>%
    tibble()

  population_summary_stats2 = pop_data_by_stratum %>%
    select(all_of(var_names)) %>%
    summarise_all(list(mean, sd)) %>%
    mutate_all(round, digits = 3)

  population_summary_stats = population_summary_stats2 %>%
    names() %>% str_sub(end = -5) %>%
    unique() %>%
    lapply(

      function(x) {

        unite(population_summary_stats2,
               {{ x }},
               grep(x, names(population_summary_stats2), value = TRUE),
               sep = ' / ',
               remove = TRUE) %>%
        select(all_of(x))
        }

      ) %>%
    bind_cols()

  summary_stats = pop_data_by_stratum %>%
    select(all_of(var_names), Stratum) %>%
    group_by(Stratum) %>%
    summarize_if(is.numeric, mean) %>%
    left_join(

      (pop_data_by_stratum %>%
         select(all_of(var_names), Stratum) %>%
         group_by(Stratum) %>%
         summarize_if(is.numeric, sd)),

      by = "Stratum",
      suffix = c("_fn1", "_fn2")

      ) %>%
    mutate_all(round, digits = 3)

  summary_stats2 = summary_stats %>%
    select(-Stratum) %>%
    names() %>%
    str_sub(end = -5) %>%
    unique() %>%
    lapply(

      function(x) {

        unite(summary_stats,
               {{ x }},
               grep(x, names(summary_stats), value = TRUE),
               sep = ' / ',
               remove = TRUE) %>%
        select(all_of(x))
        }

      ) %>%
    bind_cols() %>%
    mutate(Stratum = summary_stats$Stratum) %>%
    select(Stratum, everything()) %>%
    left_join(

      (pop_data_by_stratum %>%
         group_by(Stratum) %>%
         count()),

      by = "Stratum") %>%
    mutate(Stratum = as.character(Stratum)) %>%
    add_row(tibble_row(Stratum = "Population",
                       population_summary_stats,
                       n = dim(pop_data_by_stratum)[1]))

  simtab_m = population_summary_stats2 %>%
    select(contains("fn1"))

  names(simtab_m) = names(simtab_m) %>%
    str_sub(end = -5)

  sd_tab = summary_stats %>%
    select(contains("fn2")) %>%
    add_row(tibble_row((population_summary_stats2 %>% select(contains("fn2")))))

  names(sd_tab) = names(sd_tab) %>%
    str_sub(end = -5)

  sd_tab = sd_tab %>%
    mutate(Stratum = summary_stats2$Stratum) %>%
    pivot_longer(-Stratum,
                 names_to = "Variable",
                 values_to = "sd")

  mean_tab = summary_stats %>%
    select(contains("fn1")) %>%
    add_row(tibble_row((population_summary_stats2 %>% select(contains("fn1")))))

  names(mean_tab) = names(mean_tab) %>%
    str_sub(end = -5)

  mean_tab = mean_tab %>%
    mutate(Stratum = summary_stats2$Stratum) %>%
    pivot_longer(-Stratum, names_to = "Variable", values_to = "mn")

  counts_tab = summary_stats2 %>%
    select(Stratum, n)

  heat_data = left_join(mean_tab, sd_tab, by = c("Stratum", "Variable")) %>%
    filter(Variable != "rank")

  heat_data = heat_data %>%
    left_join(counts_tab, by = "Stratum")

  temporary_df = data.frame(Variable = unique(heat_data$Variable),
                             pop_mean = (heat_data %>% filter(Stratum == "Population") %>% select(mn)),
                             pop_sd = (heat_data %>% filter(Stratum == "Population") %>% select(sd)),
                             pop_n = (heat_data %>% filter(Stratum == "Population") %>% select(n))) %>%
    mutate(pop_mean = mn,
           pop_sd = sd,
           pop_n = n) %>%
    select(-mn,-sd,-n)

  heat_data = heat_data %>%
    left_join(temporary_df,
              by = "Variable") %>%
    mutate(deviation = case_when((mn - pop_mean) / pop_mean >= 0.7 ~ 0.7,
                                 (mn - pop_mean) / pop_mean <= -0.7 ~ -0.7,
                                 TRUE ~ (mn - pop_mean) / pop_mean)
    )

  heat_data_simple = heat_data %>%
    select(Stratum, Variable, mn, sd) %>%
    pivot_wider(names_from = Stratum,
                values_from = c(mn, sd),
                names_glue = "{Stratum}_{.value}")

  heat_data_simple = heat_data_simple %>%
    select(order(colnames(heat_data_simple))) %>%
    select(Variable, everything())

  header = c(1, rep(2, n_strata+1))
  header_names = " "

  for(i in 1:n_strata) {

    header_names = header_names %>%
      append(paste0("Stratum ",
                    i,
                    "\nn = ",
                    counts_tab$n[i]))
  }

  names(header) = header_names %>%
    append(paste0("Population\n",
                  "n = ",
                  counts_tab$n %>% tail(n = 1)))

  heat_data_kable = heat_data_simple %>%
    kbl(caption = "Covariate Statistics by Stratum",
        align = "c",
        col.names = c("Variable", rep(c("Mean", "Standard Deviation"), n_strata + 1))) %>%
    kable_styling(c("striped", "hover"),
                  fixed_thead = TRUE) %>%
    add_header_above(header)

  stratum_labels = "Population"

  for(i in 2:(n_strata + 1)) {

    stratum_labels[i] = paste("Stratum", (i - 1))
  }

  heat_plot = ggplot(data = heat_data) +
    geom_tile(aes(x = Stratum, y = Variable, fill = deviation), color = "black", width = 0.95) +
    geom_text(aes(x = Stratum, y = ((ncol(summary_stats) + 1)/2 - 0.15), label = paste(n, "\nunits")), size = 3.4) +
    scale_fill_gradientn(name = NULL,
                         breaks=c(-0.5, 0, 0.5),
                         labels = c("50% \nBelow Mean",
                                    "Population\nMean",
                                    "50% \nAbove Mean"),
                         colours = c("#990000", "#CC0000",
                                     "white", "#3D85C6",
                                     "#0B5294"),
                         limits = c(-0.7, 0.7)) +
    scale_x_discrete(position = "top",
                     expand = c(0, 0),
                     labels = c(stratum_labels[-1], "Population")) +
    expand_limits(y = c(0, (ncol(summary_stats) + 1)/2 + 0.1)) +
    labs(y = NULL, x = NULL) +
    theme(panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 10, colour = "grey15"),
          legend.key.height = unit(1, "cm"),
          legend.text = element_text(size = 10),
          legend.position = "right")

  longer_heat_data = heat_data %>%
    pivot_longer(cols = c("mn", "sd"),
                 names_to = "mn_or_sd",
                 values_to = "values") %>%
    mutate(values = values %>% round(1))

  heat_plot_final = heat_plot +
    new_scale("fill") +
    geom_label_repel(data = longer_heat_data %>% filter(mn_or_sd == "mn"),
                     aes(x = Stratum,
                         y = Variable,
                         label = values,
                         fill = mn_or_sd),
                     min.segment.length = 1,
                     direction = "y",
                     nudge_y = 0.05,
                     alpha = 0.7,
                     size = ifelse((length(levels(heat_data$Variable %>% factor())) + 1) > 7, 2, 3.5)) +
    geom_label_repel(data = longer_heat_data %>% filter(mn_or_sd == "sd"),
                     aes(x = Stratum,
                         y = Variable,
                         label = values,
                         fill = mn_or_sd),
                     min.segment.length = 1,
                     direction = "y",
                     nudge_y = -0.05,
                     alpha = 0.7,
                     size = ifelse((length(levels(heat_data$Variable %>% factor())) + 1) > 7, 2, 3.5)) +
    scale_fill_viridis(labels = c("Mean", "Standard Deviation"),
                       begin = 0.4,
                       end = 0.8,
                       direction = 1,
                       discrete = TRUE,
                       option = "D") +
    theme(legend.title = element_blank()) +
    guides(fill = guide_legend(override.aes = aes(label = "")))

  # 5) Save output

  out = list(idvar = idvar,
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

  class(out) = "generalizer_stratify"

  return(invisible(out))
}

print.generalizer_stratify = function(x,...) {

  cat("A generalizer_stratify object: \n")

  cat(paste0(" - Dataset used: ", crayon::bold(x$dataset),"\n"))

  cat(paste0(" - Stratification variables: "))

  cat(paste0(crayon::bold$blue(x$variables)), sep = ", ")

  cat(paste0("\n"))

  cat(paste0(" - Number of strata: ", crayon::bold(x$n_strata), "\n"))

  cat(paste0(" - Number of observations dropped due to missing data: ", crayon::bold(nrow(x$data_omitted)), "\n"))

  cat(paste0("   (see $data_omitted for dropped observations)"))

  invisible(x)
}

summary.generalizer_stratify = function(object,...) {

  out = object

  class(out) = "summary.generalizer_stratify"

  return(out)
}

print.summary.generalizer_stratify = function(x,...) {

  cat("============================================ \n")

  cat(paste0("Summary of stratification performed with '", x$dataset,"' dataset:", "\n", "\n"))

  cat(paste0("Stratification Variables: "))

  cat(paste0(crayon::blue$bold(x$variables)), sep = ", ")

  cat("\n")

  cat(paste0("Observations dropped due to missing data: ", crayon::bold(nrow(x$data_omitted)), " (see $data_omitted)\n"))

  cat(paste0("Population size: ", crayon::bold(nrow(x$pop_data_by_stratum)),"\n"))

  cat(paste0("Number of strata specified: ", crayon::bold(x$n_strata), "\n"))

  cat(paste0("Proportion of variation in population explained by strata: "))

  cat(bold(paste(100 * round(x$solution$between.SS_DIV_total.SS, 4), "%", sep = "")))

  cat("\n")

  cat("============================================ \n")

  cat("Covariate Statistics by Stratum: \n\n")

  x$heat_data_simple %>%
    as.data.frame() %>%
    print()

  x$heat_data_kable %>%
    print()

  x$heat_plot %>%
    print()


  # tryCatch(
  #
  #   {
  #     x$heat_plot %>%
  #       print()
  #   },
  #
  #   error = function(cond) {
  #
  #     message("Your Plots pane is too small for the heat map to be displayed. \nIf you still want to view the plot, try resizing the pane and \nthen running 'x$heat_plot' where 'x' is the name assigned to your generalizer_stratify object.")
  #     message("Here's the original error message:")
  #     message(cond)
  #     return(NA)
  #   }
  # )

  invisible(x)
}
