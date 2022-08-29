

stratify_guided = function(data) {

  source("R/make_var_overview.R")
  source("R/select_list.R")
  source("R/round_preserve_sum.R")
  source("R/check_factor_levels.R")
  source("R/make_cont_data_tbl.R")

  # 1) Introduction -----------------------------------------------------------

  cat(bold("\nWelcome to stratify()! \n\n"))

  cat("Your chosen inference population is the '",
      data_name,
      "' dataframe. \n\n",
      sep = "")

  cat("To store your results, make sure you have assigned this function to an object.\n\n")

  cat("To exit this function, press <Esc>. \n\n")

  # 2) Selection of id variable -----------------------------------------------

  is_valid_id_name = FALSE

  while(is_valid_id_name == FALSE) {

    idvar = readline(prompt = "Enter the name of the ID Variable in your dataframe: ")

    ## Check ##
    if(!idvar %in% names(data)) {

      cat(crayon::red("\nERROR: We could not find that variable. Your ID variable must be one of the variables \nin the '"),
          crayon::red(data_name),
          crayon::red("' dataframe.\n\n"),
          sep = "")
      next
    }

    is_valid_id_name = TRUE
  }

  id = data %>%
    select(all_of(idvar))

  # 3) Selection of variables for stratification ------------------------------

  data_subset = data %>%
    select(-all_of(idvar))

  data_subset %>%
    make_var_overview()

  cat("\nIn the Viewer pane to the right you will find a table that displays each variable in \nthe '",
      data_name,
      "' dataframe along with its data type and number of levels (only relevant \nfor factor variables).\n\n",
      sep = "")

  cat(crayon::yellow$bold("Please note that any character variables that may have been present in your \ndataframe have been automatically converted to factor variables.\n\n"))

  readline(prompt = "Press <Return> to continue.")

  cat("\nYou are now ready to select your stratification variables. The following are the \nvariables available in your dataset. Which key variables do you think may explain \nvariation in your treatment effect? Typically, studies include 4-6 variables for \nstratification.\n\n")

  cat(yellow$bold("You must choose at least 2 variables and you may not choose any factor variables \nwith more than 4 levels.\n"))

  var_names = data_subset %>%
    names()

  variables_are_correct = FALSE

  while(variables_are_correct == FALSE){

    variables = select_list(choices = var_names,
                            graphics = FALSE,
                            multiple = TRUE)

    # Verify that there are no factor variables with more than 4 levels
    invalid_factors = data_subset %>%
      select(all_of(variables)) %>%
      check_factor_levels()

    if(!is_empty(invalid_factors)) {

      cat(crayon::red("ERROR: This function will not allow a factor variable to have more than 4 levels.\n\n"),
          crayon::red("The following factor variables have more than 4 levels:\n\n"),
          paste(crayon::blue$bold(invalid_factors), collapse = ", "),
          crayon::red("\n\nPlease exit out of stratify() by pressing <Esc> and re-code your desired factor levels \nfrom these variables as dummy variables (see the package 'fastDummies') or press \n<Return> to choose a different set of variables.\n"),
          sep = "")

      readline(prompt = "")

      next
    }

    cat("\nYou have selected the following stratifying variables:\n\n",
        paste(crayon::blue$bold(variables), collapse = ", "),
        "\n\n",
        sep = "")

    data_subset %>%
      select(all_of(variables)) %>%
      make_var_overview(print_to_console = TRUE)

    if(menu(choices = c("Yes", "No"), title = cat("\nIs this correct?")) == 1) {

      variables_are_correct = TRUE
    }

    else {

      data_subset %>%
        make_var_overview()
    }
  }

  data_subset = data_subset %>%
    select(all_of(variables))

  # 4) Let user know about missing observations -------------------------------

  missing_obs_table = data.frame(id, data_subset) %>%
    sapply(function(x) sum(is.na(x))) %>%
    as.data.frame() %>%
    `colnames<-`("n_missing") %>%
    rownames_to_column("variable")

  n_missing = data.frame(id, data_subset) %>%
    filter(if_any(everything(), is.na)) %>%
    nrow()

  data_subset = data.frame(id, data_subset) %>%
    drop_na() %>%
    select(-all_of(idvar))

  cat("\nThe following table shows how many observations are missing for the variables you have \nchosen.\n\n")

  missing_obs_table %>%
    kbl(caption = "Missing Observations by Variable",
        align = "l",
        col.names = c("Variable", "Number Missing")) %>%
    kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    print()

  missing_obs_table %>%
    print(row.names = FALSE)

  cat("\nNote that there may be multiple missing observations per row, meaning that the total \nnumber of dropped observations won't necessarily be the sum of the number of missing \nobservations in each variable.\n\n")

  cat(crayon::yellow$bold("In total,", n_missing, "rows will be dropped from the inference population before \nstratification due to missing observations.\n\n"))

  readline(prompt = "Hit <Return> to proceed once you have viewed the missing observations.")

  # 5) Overview of categorical variables --------------------------------------

  cat_data = data_subset %>%
    select_if(is.factor)

  cat_data_vars = names(cat_data)

  if(dim(cat_data)[2] >= 1L) {

    cat_data_plot = data.frame(cat_data)

    cat("\nPlease review the descriptive statistics of your categorical variables (factors). ",
        "Bar \ncharts and tables for each variable will also be printed in the Plots and Viewer \npanes to the right.\n", sep = "")

    n_cat_vars = ncol(cat_data_plot)

    fill_colors_cat = plasma(n_cat_vars,
                             alpha = 0.7,
                             direction = sample(c(-1, 1),
                                                size = 1)) %>%
      sample()

    outline_colors_cat = turbo(n_cat_vars) %>%
      sample()

    for(i in 1:n_cat_vars) {

      var_name = cat_data_vars[i]

      levels(cat_data_plot[[var_name]]) = str_wrap(levels(cat_data_plot[[var_name]]), width = 10)

      barfig = cat_data_plot %>%
        group_by(across(all_of(var_name))) %>%
        summarise(count = n()) %>%
        mutate(ordered_factor = fct_reorder(.[[var_name]], count)) %>%
        ggplot(aes(x = ordered_factor, y = count)) +
        geom_col(fill = fill_colors_cat[i],
                 color = outline_colors_cat[i]) +
        theme_minimal() +
        xlab(var_name) +
        labs(title = paste("Bar Chart of", var_name))

      cat("\n")

      print(barfig)

      par(ask = TRUE)

      cat("\nNumber of Observations in Levels of Factor ",
          paste(blue$bold(var_name)),
          ":\n\n",
          sep = "")

      cat_data_table = cat_data_plot[,i] %>%
        table() %>%
        data.frame() %>%
        `colnames<-` (c("Level", "Frequency")) %>%
        print(row.names = FALSE)

      cat_data_table %>%
        kbl(col.names = c("Level", "Frequency"),
            caption = paste("Number of Observations in Levels of Factor ", var_name),
            align = "l") %>%
        kable_styling(c("striped", "hover")) %>%
        print()

    }

    cat("\n")

    readline(prompt = "Hit <Return> to proceed once you have reviewed the categorical variables.")
  }

  # 6) Overview of continuous variables ---------------------------------------

  cont_data = data_subset %>%
    select_if(negate(is.factor))

  cont_data_vars = cont_data %>%
    names()

  if(dim(cont_data)[2] >= 1L) {

    cat("\nPlease review the descriptive statistics of your continuous variables. Histograms and \ntables for each variable will also be printed in the Plots and Viewer panes to the \nright. \n\n")

    n_cont_vars = ncol(cont_data)

    fill_colors_cont = viridis(n_cont_vars,
                               alpha = 0.7,
                               direction = sample(c(-1, 1),
                                                  size = 1)) %>%
      sample()

    outline_colors_cont = turbo(n_cont_vars) %>%
      sample()

    for(i in 1:n_cont_vars) {

      hist = cont_data %>%
        ggplot(aes(x = cont_data[,i])) +
        geom_histogram(bins = 30,
                       fill = fill_colors_cont[i],
                       color = outline_colors_cont[i]) +
        theme_minimal() +
        xlab(cont_data_vars[i]) +
        labs(title = paste("Histogram of", cont_data_vars[i]))

      hist %>%
        print()

      par(ask = TRUE)
    }

    readline(prompt = "Hit <Return> to view a table of summary statistics for the continuous variables:\n")

    sumstats = cont_data %>%
      make_cont_data_tbl()

    sumstats %>%
      print(row.names = FALSE)

    sumstats %>%
      kbl(col.names = c("Variable", "Min", "Pct50", "Max", "Mean", "SD")) %>%
      kable_styling(c("striped", "hover")) %>%
      print()
  }

  par(ask = FALSE)

  cat("\n")

  readline(prompt = "Hit <Return> to proceed.")

  # 7) Selection of number of strata ------------------------------------------

  cat("\nStratification will help you develop a recruitment plan so that your study will result \nin an unbiased estimate of the ", bold("average treatment effect (ATE)"), ". Without using strata, \nit is easy to end up with a sample that is very different from your inference \npopulation. \n\nGeneralization works best when strata are ", bold("homogeneous"), ". That means units within each \nstratum are almost identical in terms of relevant variables.\n\n", sep = "")

  cat("Enter the number of strata in which you wish to divide your population. Typically, ",
      bold("\nthe more strata"),
      ", ",
      bold("the better"),
      "; with fewer strata, units in each stratum are no longer \nidentical. However, increasing ",
      "the number of strata uses more resources, because you \nmust sample a given number of units ",
      "from each stratum. Choosing 4-6 strata is common. \n\nTry a few numbers and choose the 'best' one for you.\n\n",
      sep = "")

  n_strata_correct = FALSE

  while(!n_strata_correct) {

    n_strata = suppressWarnings(as.numeric(readline(prompt = "Number of strata: ")))

    if(is.na(n_strata) || !is.numeric(n_strata) || n_strata <= 1 || n_strata %% 1 != 0) {

      cat(red("\nERROR: The number of strata must be a single positive integer greater than 1.\n\n"))

      next
    }

    n_strata_correct = TRUE
  }

  # 8) Save output

  out = list(n_strata = n_strata,
              variables = variables,
              idvar = idvar)

  return(invisible(out))
}






