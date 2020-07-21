#' Stratify a Population Data Frame
#'
#' This function takes as input any data frame that you want to stratify into clusters. Typically, the goal of such stratification is sampling for generalizability. This function, and the others in this package, are designed to mimic the website https://www.thegeneralizer.org/.
#'
#' @param data The R object containing your population data frame
#' @param guided logical; defaults to TRUE. Whether the function should be guided (ask questions and behave interactively throughout) or not. If set to FALSE, must provide values for other arguments below
#' @param n_strata defaults to NULL. If guided is set to FALSE, must provide a number of strata to cluster population into
#' @param variables defaults to NULL. If guided is set to FALSE, must provide a character vector of the names of stratifying variables (from population data frame)
#' @param idnum defaults to NULL. If guided is set to FALSE, must provide a character vector of the name of the ID variable (from population data frame)
#' @return A list of class "generalizer_output" that can be given to other functions in the package
#' @export
#' @importFrom graphics par
#' @importFrom stats mahalanobis median na.omit sd var
#' @importFrom utils menu select.list
#' @importFrom crayon %+% blue bold
#' @importFrom janitor clean_names
#' @importFrom ggplot2 ggplot aes geom_bar xlab labs geom_histogram geom_text geom_label geom_hline scale_fill_gradientn scale_x_discrete expand_limits geom_tile element_blank element_text theme
#' @importFrom ggthemes theme_base
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

stratify <- function(data, guided = TRUE, n_strata = NULL, variables = NULL,
                     idnum = NULL){

  skim_variable <- skim_type <- variable <- NULL
  type <- clusterID <- n <- mn <- deviation <- NULL

  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width") - 1L), collapse = " "));

  if(guided == TRUE){
    cat(bold("If you want to store your results, make sure you assign \nthis function to an object.\n\n"))
    cat("Your chosen inference population is the '",
        deparse(substitute(data)), "' dataset.", sep = "")

    cat("\n")
    cat("\n")

    idnum <- readline(prompt = "Enter the name of the ID Variable in your dataset: ")

    if(!idnum %in% names(data))
    stop(simpleError("We could not find that variable. Please make sure your \ndataset contains an ID variable."))

    cat("\nIf you want to adjust or restrict your inference population \n(e.g., if you are interested in only one location, etc.), \nmake sure that you have altered the data frame appropriately. \nIf you need to alter your data frame, you can exit this \nfunction, use " %+% blue$bold("dplyr::filter()") %+% " or " %+% blue$bold("set_percentile_limits()") %+% ", and \nreturn.\n")

    if(menu(choices = c("Yes", "No"), title = cat("\nDo you wish to proceed?")) == 1){

    }else{
      stop(simpleError(blankMsg))
    }

    id <- data %>% select(all_of(idnum))
    data <- data %>% select(-all_of(idnum))

    cat("\nYou're now ready to select your stratification variables. \nThe following are the variables available in your dataset.")

    names <- names(data)
    variables <- select.list(choices = names,
                     title = cat("\nWhich key variables do you think may explain variation \nin your treatment effect?",
                                 "Typically, studies include \nup to 10 variables for stratification.\n"),
                     graphics = FALSE, multiple = TRUE)

    if(length(variables) >= 1){
      data <- data %>%
        select(all_of(variables))
    }else{
      stop("You have to select some stratifying variables.")
    }

    var_overview <- skimr::skim(data) %>% tibble() %>%
      distinct(skim_variable, skim_type) %>%
      mutate(variable = skim_variable, type = skim_type) %>%
      select(variable, type) %>%
      data.frame()

    colnames(var_overview) <- c("Variable", "Type")

    cat("\nYou have selected the following stratifying variables: \n")
    cat(paste(blue$bold(colnames(data)), collapse = ", "), ".\n\n", sep = "")
    print(var_overview, row.names = FALSE)

    if(menu(choices = c("Yes", "No"), title = cat("\nIs this correct?")) == 1){

    }else{
      stop(simpleError(blankMsg))
    }

    cat_data <- data %>%
      select_if(is.factor)
    cat_data_vars <- names(cat_data)
    if(dim(cat_data)[2] >= 1){
      cat_data_plot <- data.frame(cat_data)
      cat("Please review the descriptive statistics of your \ncategorical variables (factors). Note that these will \nautomatically be converted to dummy variables for analysis.\n")
      for(i in 1:(ncol(cat_data_plot))){
        var_name <- cat_data_vars[i]
        cat("\nNumber of Observations in Levels of Factor ", paste(blue$bold(var_name)), ":\n", sep = "")
        print(table(cat_data_plot[,i]))
        barfig <- ggplot(data = cat_data_plot, aes(x = cat_data_plot[,i])) +
          geom_bar() +
          theme_base() +
          xlab(var_name) +
          labs(title = paste("Bar Chart of", var_name))
        print(barfig)
        cat("\n")
        par(ask = TRUE)
      }
    }

    cont_data <- data %>%
      select_if(negate(is.factor))
    cont_data_vars <- names(cont_data)
    if(dim(cont_data)[2] >= 1){
      sumstats <- cont_data %>%
        map_df(function(x){
          tibble(min = min(x), pct50 = median(x), max = max(x), mean = mean(x), sd = sd(x))
        }) %>%
        mutate(variable = cont_data_vars) %>%
        select(variable, everything()) %>%
        clean_names() %>%
        data.frame()

      cat("Please review the descriptive statistics of your \ncontinuous variables.\n\n")
      print(sumstats, row.names = FALSE)
      for(i in 1:ncol(cont_data)){
        cont_data_plot <- cont_data %>% data.frame()
        suppressWarnings(
          suppressMessages(
            hist <- ggplot(data = cont_data_plot, aes(x = cont_data_plot[,i])) +
              geom_histogram(bins = 30) +
              theme_base() +
              xlab(cont_data_vars[i]) +
              labs(title = paste("Histogram of", cont_data_vars[i]))
          )
        )
        print(hist)
        par(ask = TRUE)
      }
    }
    par(ask = FALSE)

    satisfied <- 0

    while(satisfied != 1){
      cat("Enter a number of strata to divide your population into. Typically, \nthe more strata, the better. However, increasing the number of strata \nuses more resources, because you must sample a given number of units \nfrom each stratum. Therefore, choose a larger number if possible, and \nonly if you have the resources to accommodate it. Otherwise, \nchoose a smaller number.")
      n_strata <- as.numeric(readline(prompt = "# of strata: "))
      ## Add a catch here, similar to before: MUST enter a number

      cat("This might take a little while. Please bear with us.")

      if(dim(cat_data)[2] >= 1){
        cat_data <- fastDummies::dummy_cols(cat_data, remove_first_dummy = TRUE)
        data_full <- cbind(cat_data, cont_data) %>%
          na.omit()
      }else{
        data_full <- cont_data %>%
          na.omit()
      }

      suppressWarnings(distance <- daisy(data_full, metric = "gower"))
      cat("\n1: Calculated distance matrix.")
      solution <- KMeans_rcpp(as.matrix(distance), clusters = n_strata, verbose = TRUE)

      x2 <- data.frame(id, data_full, clusterID = solution$clusters)
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

        ## Note: pct_female and pct_white work, but other combos throw an error
        ## at the next line. Hunch is that pct_female and pct_white might not
        ## have NA values and NAs might be what are messing it up. To investigate.

        if(any(a == 0)){ a[which(a == 0)] <- 0.00000001 }
        cov.dat <- diag(a)
        ma.s <- mahalanobis(dat4, mu, cov.dat)
        final_dat4 <- data.frame(idvar, dat4, distance = ma.s, clusterID = dat3$clusterID) %>% tibble()
        recruitment_lists[[i]] <- final_dat4 %>% # Produces a list of data frames, one per stratum, sorted by
          # distance (so the top N schools in each data frame are the "best," etc.)
          arrange(distance) %>%
          mutate(rank = seq.int(nrow(final_dat4))) %>%
          select(rank, idnum)
  }

      cat(blue$bold("Congratulations, you have successfully grouped your data into", n_strata, "strata!\n"))

      readline(prompt = "Press [enter] to view the results")

      cat("\nYou have specified ")
      cat(bold(n_strata))
      cat(" strata, which explain ")
      cat(paste(bold(100 * round(solution$between.SS_DIV_total.SS, 4), "%", sep = "")))
      cat(" of the total \nvariation in the population.")

      cat("\n\nThe following table presents the mean and standard deviation \n(mean / sd) of each stratifying variable for each stratum. \nThe bottom row, 'Population,' presents the average values for \nthe entire inference population. The last column, 'n,' lists the \ntotal number of units in the inference population that fall \nwithin each stratum.\n\n")

      x2 <- data.frame(id, data_full, clusterID = solution$clusters) %>% tibble()

      population_summary_stats2 <- x2 %>% select(-c(all_of(idnum), clusterID)) %>%
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
        group_by(clusterID) %>%
        summarize_if(is.numeric, mean) %>%
        left_join((x2 %>% group_by(clusterID) %>% summarize_if(is.numeric, sd)),
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
        left_join((x2 %>% group_by(clusterID) %>% count())) %>%
        mutate(clusterID = as.character(clusterID)) %>%
        add_row(tibble_row(clusterID = "Population", population_summary_stats, n = dim(x2)[1])) %>%
        data.frame()

      print(summary_stats2)

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

      heat_data <- left_join(mean_tab, sd_tab) %>%
        left_join(counts_tab)
      temporary_df <- data.frame(variable = unique(heat_data$variable),
                        population_mean = (heat_data %>% filter(clusterID == "Population") %>% select(mn))) %>%
        mutate(population_mean = mn) %>%
        select(-mn)
      heat_data <- heat_data %>% left_join(temporary_df) %>%
        mutate(deviation = case_when((mn - population_mean)/population_mean >= 0.7 ~ 0.7,
                                     (mn - population_mean)/population_mean <= -0.7 ~ -0.7,
                                     TRUE ~ (mn - population_mean)/population_mean))
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


      readline(prompt = "Press [enter] to continue:")

      print(heat_plot_final)

      if(menu(choices = c("Yes", "No"), title = cat("\nWould you like to go back and specify a different number of strata?")) == 2){

        satisfied <- 1

      }else{

        satisfied <- 0

      }

    }

    ## To add:
    ## If you want to exit, do you want to save your output
    ## (and if so, as what)?

  }else{
    par(ask = FALSE)

    # This is where all the non-guided stuff goes

    cat("Your chosen inference population is the '",
        deparse(substitute(data)), "' dataset.", sep = "")
    cat("\n")

    id <- data %>% select(all_of(idnum))
    data <- data %>% select(-all_of(idnum))

    if(length(variables) >= 1){
      data <- data %>%
        select(all_of(variables))
    }else{
      stop("You have to select some stratifying variables.")
    }

    cat_data <- data %>%
      select_if(is.factor)
    cat_data_vars <- names(cat_data)
    if(dim(cat_data)[2] >= 1){
      cat_data_plot <- data.frame(cat_data)
      cat("Please review the descriptive statistics of your \ncategorical variables (factors). Note that these will \nautomatically be converted to dummy variables for analysis.\n")
      for(i in 1:(ncol(cat_data_plot))){
        var_name <- cat_data_vars[i]
        cat("\nNumber of Observations in Levels of Factor ", paste(blue$bold(var_name)), ":\n", sep = "")
        print(table(cat_data_plot[,i]))
        barfig <- ggplot(data = cat_data_plot, aes(x = cat_data_plot[,i])) +
          geom_bar() +
          theme_base() +
          xlab(var_name) +
          labs(title = paste("Bar Chart of", var_name))
        print(barfig)
        cat("\n")
      }
    }

    cont_data <- data %>%
      select_if(negate(is.factor))
    cont_data_vars <- names(cont_data)
    if(dim(cont_data)[2] >= 1){
      sumstats <- cont_data %>%
        map_df(function(x){
          tibble(min = min(x), pct50 = median(x), max = max(x), mean = mean(x), sd = sd(x))
        }) %>%
        mutate(variable = cont_data_vars) %>%
        select(variable, everything()) %>%
        clean_names() %>%
        data.frame()

      cat("Please review the descriptive statistics of your \ncontinuous variables.\n\n")
      print(sumstats, row.names = FALSE)
      for(i in 1:ncol(cont_data)){
        cont_data_plot <- cont_data %>% data.frame()
        suppressWarnings(
          suppressMessages(
            hist <- ggplot(data = cont_data_plot, aes(x = cont_data_plot[,i])) +
              geom_histogram(bins = 30) +
              theme_base() +
              xlab(cont_data_vars[i]) +
              labs(title = paste("Histogram of", cont_data_vars[i]))
          )
        )
        print(hist)
      }
    }

    cat("This might take a little while. Please bear with us.")

    if(dim(cat_data)[2] >= 1){
      cat_data <- fastDummies::dummy_cols(cat_data, remove_first_dummy = TRUE)
      data_full <- cbind(cat_data, cont_data) %>%
        na.omit()
    }else{
      data_full <- cont_data %>%
        na.omit()
    }

    suppressWarnings(distance <- daisy(data_full, metric = "gower"))
    cat("\n1: Calculated distance matrix.")
    solution <- KMeans_rcpp(as.matrix(distance), clusters = n_strata, verbose = TRUE)

    x2 <- data.frame(id, data_full, clusterID = solution$clusters)
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

      ## Note: pct_female and pct_white work, but other combos throw an error
      ## at the next line. Hunch is that pct_female and pct_white might not
      ## have NA values and NAs might be what are messing it up. To investigate.

      if(any(a == 0)){ a[which(a == 0)] <- 0.00000001 }
      cov.dat <- diag(a)
      ma.s <- mahalanobis(dat4, mu, cov.dat)
      final_dat4 <- data.frame(idvar, dat4, distance = ma.s, clusterID = dat3$clusterID) %>% tibble()
      recruitment_lists[[i]] <- final_dat4 %>% # Produces a list of data frames, one per stratum, sorted by
        # distance (so the top N schools in each data frame are the "best," etc.)
        arrange(distance) %>%
        mutate(rank = seq.int(nrow(final_dat4))) %>%
        select(rank, idnum)
    }

    cat(blue$bold("Congratulations, you have successfully grouped your data into", n_strata, "strata!\n"))

    cat("\nYou have specified ")
    cat(bold(n_strata))
    cat(" strata, which explain ")
    cat(paste(bold(100 * round(solution$between.SS_DIV_total.SS, 4), "%", sep = "")))
    cat(" of the total \nvariation in the population.")

    cat("\n\nThe following table presents the mean and standard deviation \n(mean / sd) of each stratifying variable for each stratum. \nThe bottom row, 'Population,' presents the average values for \nthe entire inference population. The last column, 'n,' lists the \ntotal number of units in the inference population that fall \nwithin each stratum.\n\n")

    x2 <- data.frame(id, data_full, clusterID = solution$clusters) %>% tibble()

    population_summary_stats2 <- x2 %>% select(-c(all_of(idnum), clusterID)) %>%
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
      group_by(clusterID) %>%
      summarize_if(is.numeric, mean) %>%
      left_join((x2 %>% group_by(clusterID) %>% summarize_if(is.numeric, sd)),
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
      left_join((x2 %>% group_by(clusterID) %>% count())) %>%
      mutate(clusterID = as.character(clusterID)) %>%
      add_row(tibble_row(clusterID = "Population", population_summary_stats, n = dim(x2)[1])) %>%
      data.frame()

    print(summary_stats2)

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

    heat_data <- left_join(mean_tab, sd_tab) %>%
      left_join(counts_tab)
    temporary_df <- data.frame(variable = unique(heat_data$variable),
                               population_mean = (heat_data %>% filter(clusterID == "Population") %>% select(mn))) %>%
      mutate(population_mean = mn) %>%
      select(-mn)
    heat_data <- heat_data %>% left_join(temporary_df) %>%
      mutate(deviation = case_when((mn - population_mean)/population_mean >= 0.7 ~ 0.7,
                                   (mn - population_mean)/population_mean <= -0.7 ~ -0.7,
                                   TRUE ~ (mn - population_mean)/population_mean))
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

    print(heat_plot_final)

  }

  overall_output <- list(x2 = x2, solution = solution, n_strata = n_strata,
                         recruitment_lists = recruitment_lists,
                         population_summary_stats2 = population_summary_stats2,
                         summary_stats = summary_stats,
                         summary_stats2 = summary_stats2,
                         heat_data = heat_data, heat_plot_final = heat_plot_final,
                         idnum = idnum, variables = variables)

  class(overall_output) <- c("generalizer_output")

  return(invisible(overall_output))
  # Note to self; you CAN store plots in list of output!

}

