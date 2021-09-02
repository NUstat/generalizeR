stratify_basic <- function(data, n_strata = NULL, variables = NULL,
                           idnum = NULL, seed = 7835, verbose = TRUE){

  skim_variable <- skim_type <- variable <- NULL
  type <- clusterID <- n <- mn <- deviation <- NULL

  if(is.null(n_strata) | is.null(variables) | is.null(idnum)){
    stop(simpleError("You must specify n_strata, variables, and idnum as arguments if you are running the non-guided version of this function."))
  }

  if(!is.numeric(n_strata)){
    stop(simpleError("The number of strata must be a number."))
  }

  if((length(n_strata) > 1)){
    stop(simpleError("Only specify one number of strata."))
  }

  if(n_strata <= 1){
    stop(simpleError("The number of strata must be a positive number greater than 1."))
  }

  if(n_strata%%1==0){
    n_strata <- round(n_strata)
  }

  if(!is.character(variables) | (anyNA(match(variables, names(data))))){
    stop(simpleError("You must provide a character vector consisting of the names of stratifying variables in your inference population."))
  }

  if(!is.character(idnum) | is.na(match(idnum, names(data)))){
    stop(simpleError("idnum should be the name of the identifying variable in your inference population -- e.x.: 'id'."))
  }

  # 1) Store all information given by user

  data_name <- deparse(substitute(data)) #store name of dataset

  data <- data %>% select(all_of(variables), all_of(idnum)) # subset data to relevant columns only
  data_full <- data #store a full version of the dataset
  data <- data %>% na.omit() #subset strat data to only covariates
  data_omitted <- data_full %>% anti_join(data, by = all_of(idnum)) #save dropped observations
  id <- data %>% select(all_of(idnum)) #extract id variable
  data <- data %>% select(all_of(variables))

  # 2) Check if there is categorical data.
  #    If there is, and there are not more than 4 factors,
  #    convert categorical data into factors

  cat_data <- data %>%
    select_if(is.factor)

  factor_levels_over_4 <- (cat_data %>% sapply(nlevels) > 4L) %>%
    which() %>% names()

  if(!is_empty(factor_levels_over_4)){
    stop(paste0("The following factor variables have more than 4 levels: ",
               factor_levels_over_4,
               "\n4 is the maximum number of levels this function will allow a factor to have.",
               "\nPlease exit out of this function (Press 'Esc') and re-code your desired factor",
               "\nlevels from these variables as dummy variables (see the package 'fastDummies').\n"))
  }

  if(dim(cat_data)[2] >= 1L & dim(cat_data)[2] <= 4L){
    cat_data <- cat_data %>%
      fastDummies::dummy_cols(remove_first_dummy = TRUE) %>%
      select_if(negate(is.factor))
  }

  data_full <- data %>%
    select_if(negate(is.factor)) %>%
    cbind(cat_data) %>%
    clean_names()

  # 3) Save summaries of variables in dataset
  pop_stats <- data %>%
    map_df(function(x){
      tibble(min = min(x), pct50 = median(x), max = max(x), mean = mean(x), sd = sd(x))
    }) %>%
    mutate_all(round, digits = 3) %>%
    mutate(variable = names(data)) %>%
    select(variable, everything()) %>%
    data.frame() %>%
    clean_names()

  # 4) Perform stratification - cluster with KMeans

  if(verbose == TRUE) {
    cat("\n\nThis might take a little while. Please bear with us.")
    suppressWarnings(distance <- daisy(data_full, metric = "gower"))
    cat("\n1: Calculated distance matrix.")

    solution <- KMeans_rcpp(as.matrix(distance), clusters = n_strata, verbose = TRUE)
  } else {
    suppressWarnings(distance <- daisy(data_full, metric = "gower"))
    solution <- KMeans_rcpp(as.matrix(distance), clusters = n_strata, verbose = FALSE)
  }

  x2 <- data.frame(id, data_full, clusterID = solution$clusters) %>%
    tibble()

  # 5) Create recruitment lists

  recruitment_lists <- list(NULL)

  for(i in 1:n_strata){
    dat3 <- x2 %>%
      dplyr::filter(clusterID == i)
    idvar <- dat3 %>% select(all_of(idnum))
    dat4 <- dat3 %>% select(-c(all_of(idnum), clusterID)) %>%
      mutate_all(as.numeric)

    mu <- dat4 %>% map_dbl(mean)
    v <- var(dat4)
    a <- diag(v)

    if(any(a == 0)){ a[which(a == 0)] <- 0.00000001 }
    cov.dat <- diag(a)
    ma.s <- mahalanobis(dat4, mu, cov.dat)
    final_dat4 <- data.frame(idvar, dat4, distance = ma.s, clusterID = dat3$clusterID) %>% tibble()
    recruitment_lists[[i]] <- final_dat4 %>% # Produces a list of data frames, one per stratum, sorted by
      # distance (so the top N schools in each data frame are the "best," etc.)
      arrange(distance) %>%
      mutate(rank = seq.int(nrow(final_dat4))) %>%
      select(rank, all_of(idnum))
  }

  all_lists <- do.call(rbind, recruitment_lists)

  # 6) Save full dataset and heat data

  x2 <- x2 %>% left_join(all_lists, by = all_of(idnum))

  population_summary_stats2 <- x2 %>% select(all_of(variables)) %>%
    summarise_all(list(mean, sd)) %>%
    mutate_all(round, digits = 3)

  population_summary_stats <- population_summary_stats2 %>%
    names() %>% str_sub(end = -5) %>% unique() %>%
    lapply(function(x){
      unite_(population_summary_stats2, x, grep(x, names(population_summary_stats2), value = TRUE),
             sep = ' / ', remove = TRUE) %>% select(x)
    }) %>%
    bind_cols()

  summary_stats <- x2 %>%
    select(all_of(variables), clusterID) %>%
    group_by(clusterID) %>%
    summarize_if(is.numeric, mean) %>%
    left_join((x2 %>% select(all_of(variables), clusterID) %>% group_by(clusterID) %>% summarize_if(is.numeric, sd)),
              by = "clusterID", suffix = c("_fn1", "_fn2")) %>%
    mutate_all(round, digits = 3)

  summary_stats2 <- summary_stats %>%
    select(-clusterID) %>%
    names() %>%
    str_sub(end = -5) %>%
    unique() %>%
    lapply(function(x){
      unite_(summary_stats, x, grep(x, names(summary_stats), value = TRUE),
             sep = ' / ', remove = TRUE) %>% select(x)
    }) %>%
    bind_cols() %>% mutate(clusterID = summary_stats$clusterID) %>%
    select(clusterID, everything()) %>%
    left_join((x2 %>% group_by(clusterID) %>% count()), by = "clusterID") %>%
    mutate(clusterID = as.character(clusterID)) %>%
    add_row(tibble_row(clusterID = "Population", population_summary_stats, n = dim(x2)[1]))

  simtab_m <- population_summary_stats2 %>%
    select(contains("fn1"))
  names(simtab_m) <- names(simtab_m) %>% str_sub(end = -5)
  sd_tab <- summary_stats %>%
    select(contains("fn2")) %>%
    add_row(tibble_row((population_summary_stats2 %>% select(contains("fn2")))))
  names(sd_tab) <- names(sd_tab) %>% str_sub(end = -5)
  sd_tab <- sd_tab %>%
    mutate(clusterID = summary_stats2$clusterID) %>%
    pivot_longer(-clusterID, names_to = "variable", values_to = "sd")
  mean_tab <- summary_stats %>%
    select(contains("fn1")) %>%
    add_row(tibble_row((population_summary_stats2 %>% select(contains("fn1")))))
  names(mean_tab) <- names(mean_tab) %>% str_sub(end = -5)
  mean_tab <- mean_tab %>%
    mutate(clusterID = summary_stats2$clusterID) %>%
    pivot_longer(-clusterID, names_to = "variable", values_to = "mn")
  counts_tab <- summary_stats2 %>%
    select(clusterID, n)

  heat_data <- left_join(mean_tab, sd_tab, by = c("clusterID", "variable")) %>%
    filter(variable != "rank")

  heat_data <- heat_data %>%
    left_join(counts_tab, by = "clusterID")
  temporary_df <- data.frame(variable = unique(heat_data$variable),
                             pop_mean = (heat_data %>% filter(clusterID == "Population") %>% select(mn)),
                             pop_sd = (heat_data %>% filter(clusterID == "Population") %>% select(sd)),
                             pop_n = (heat_data %>% filter(clusterID == "Population") %>% select(n))) %>%
    mutate(pop_mean = mn,
           pop_sd = sd,
           pop_n = n) %>%
    select(-mn,-sd,-n)
  heat_data <- heat_data %>% left_join(temporary_df, by = "variable") %>%
    mutate(deviation = case_when((mn - pop_mean)/pop_mean >= 0.7 ~ 0.7,
                                 (mn - pop_mean)/pop_mean <= -0.7 ~ -0.7,
                                 TRUE ~ (mn - pop_mean)/pop_mean)
    )

  heat_data_simple <- heat_data %>% select(clusterID, variable, mn, sd) %>%
    filter(clusterID != "Population") %>%
    pivot_wider(names_from = clusterID, values_from = c(mn, sd), names_glue = "{clusterID}_{.value}")

  heat_data_simple <- heat_data_simple %>%
    select(order(colnames(heat_data_simple))) %>%
    select(variable, everything())

  cluster_labels <- "Population"
  for(i in 2:(n_strata + 1)){
    cluster_labels[i] <- paste("Stratum", (i - 1))
  }

  heat_plot_final <- ggplot(data = heat_data) +
    geom_tile(aes(x = clusterID, y = variable, fill = deviation), width = 0.95) +
    geom_text(aes(x = clusterID, y = ((ncol(summary_stats) + 1)/2 - 0.15),
                  label = paste(n, "\nunits")), size = 3.4) +
    geom_label(aes(x = clusterID, y = variable,
                   label = paste0(round(mn, 1), "\n(", round(sd, 1), ")")),
               colour = "black", alpha = 0.7,
               size = ifelse((length(levels(heat_data$variable %>% factor())) + 1) > 7, 2, 3.5)) +
    geom_hline(yintercept = seq(1.5, (ncol(summary_stats) - 1), 1),
               linetype = "dotted",
               colour = "white") +
    scale_fill_gradientn(name = NULL, breaks=c(-0.5, 0, 0.5),
                         labels = c("50% \nBelow Mean",
                                    "Population\nMean",
                                    "50% \nAbove Mean"),
                         colours = c("#990000", "#CC0000",
                                     "white", "#3D85C6",
                                     "#0B5294"),
                         limits = c(-0.7, 0.7)) +
    scale_x_discrete(position = "top", expand = c(0, 0), labels = c(cluster_labels[-1], "Population")) +
    expand_limits(y = c(0, (ncol(summary_stats) + 1)/2 + 0.1)) +
    labs(y = NULL, x = NULL) +
    theme(panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_text(size = 10, colour = "grey15"),
          legend.key.height = unit(1, "cm"),
          legend.text = element_text(size = 10),
          legend.position = "right")

  # 7) Recruit proportions table

  recruit_table <- heat_data %>% select(clusterID, n) %>%
    distinct(clusterID, .keep_all = TRUE) %>%
    mutate(proportion = round(n/(dim(x2)[1]), digits = 3)) %>%
    filter(clusterID != "Population") %>%
    select(clusterID, n, proportion) %>%
    data.frame() %>%
    pivot_longer(names_to = "variable", cols = c(n,proportion)) %>%
    pivot_wider(names_from = clusterID, names_prefix = "Strata_")


  # 8) Save output

  overall_output <- list(x2 = x2, solution = solution, n_strata = n_strata,
                         data_omitted = data_omitted,
                         recruitment_lists = recruitment_lists,
                         recruit_table = recruit_table,
                         pop_stats = pop_stats,
                         heat_data = heat_data,
                         heat_data_simple = heat_data_simple,
                         heat_plot_final = heat_plot_final,
                         idnum = idnum, variables = variables, dataset = data_name
  )

  class(overall_output) <- c("generalizer_output")

  return(invisible(overall_output))

}


print.generalizer_output <- function(x,...){


  cat("A stratify() object: \n")
  cat(paste0(" - dataset used: ", x$dataset,"\n"))
  cat(paste0(" - stratification variables included: "))
  cat(paste0(x$variables), sep = ", ")
  cat(paste0("\n"))
  cat(paste0(" - no. of strata chosen: ", x$n_strata, "\n"))
  cat(paste0(" - no. of observations dropped due to missing data: ", nrow(x$data_omitted), "\n"))
  cat(paste0("   (see $data_omitted for dropped observations)"))

  invisible(x)
}

summary.generalizer_output <- function(object,...){
  out <- object

  class(out) = "summary.generalizer_output"
  return(out)
}

print.summary.generalizer_output <- function(x,...){

  cat(paste0("Summary of Stratification performed with '", x$dataset,"' dataset:", "\n", "\n"))

  cat(paste0("Observations dropped due to missing data: ", nrow(x$data_omitted), " (see $data_omitted)\n"))

  cat(paste0("Stratification Variables: "))
  cat(paste0(x$variables), sep = ", ")
  cat(paste0("\n\n"))
  cat(paste0("Variation in population: \n \n"))

  print(x$pop_stats)

  cat(paste0("\n"))

  cat(paste0("No. in population: ", bold(nrow(x$x2)),"\n"))
  cat(paste0("Number of strata specified: ", bold(x$n_strata), "\n"))
  cat(paste0("Proportion of variation in population explained by strata: "))
  cat(bold(paste(100 * round(x$solution$between.SS_DIV_total.SS, 4), "%", sep = "")))
  cat("\n")

  cat("============================================ \n")

  cat("Covariate Distributions: \n \n")

  print(x$heat_data_simple %>% as.data.frame())

  print(x$heat_plot_final)

  cat("============================================ \n")
  cat(paste0("Recruitment plan: \n \n"))
  cat(paste0(x$n_strata, " recruitment lists have been generated, one per stratum. \nThey can be accessed at $recruitment_lists.\n"))
  cat(paste0("\nEach unit is ranked in ranked in order of desirability. \nIdeally, units should be recruited across strata according to the proportions below.\nDoing so leads to the least bias and no increase in standard errors.\n\n"))
  print(x$recruit_table %>% as.data.frame())


  invisible(x)


}

