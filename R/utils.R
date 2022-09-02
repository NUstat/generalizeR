.check.factor.levels <- function(data,
                                 maxlevels = 4L) {
  invalid_factors <- data %>%
    select_if(is.factor) %>%
    sapply(nlevels) %>%
    `>`(maxlevels) %>%
    which() %>%
    names()

  return(invalid_factors)
}

.select.list <- function(choices,
                         preselect = NULL,
                         multiple = FALSE,
                         title = NULL,
                         graphics = getOption("menu.graphics")) {
  if (!interactive()) {
    stop("select_list() cannot be used non-interactively")
  }


  if (!is.null(title) && (!is.character(title) || length(title) != 1)) {
    stop("'title' must be NULL or a length-1 character vector")
  }

  if (isTRUE(graphics)) {
    if (.Platform$OS.type == "windows" || .Platform$GUI == "AQUA") {
      return(.External2(C_selectlist, choices, preselect, multiple, title))
    } else if (graphics && capabilities("tcltk") && capabilities("X11") && suppressWarnings(tcltk::.TkUp)) {
      return(tcltk::tk_select.list(choices, preselect, multiple, title))
    }
  }

  if (!multiple) {
    res <- menu(choices, FALSE, title)

    if (res < 1L || res > length(choices)) {
      return("")
    } else {
      return(choices[res])
    }
  } else {
    nc <- length(choices)

    if (length(title) && nzchar(title[1L])) {
      cat(title, "\n", sep = "")
    }

    def <- if (is.null(preselect)) {
      rep.int(FALSE, nc)
    } else {
      choices %in% preselect
    }

    op <- paste0(
      format(seq_len(nc)),
      ": ",
      ifelse(def, "+", " "),
      " ",
      choices
    )

    if (nc > 10L) {
      fop <- format(op)
      nw <- nchar(fop[1L], "w") + 2L
      ncol <- getOption("width") %/% nw

      if (ncol > 1L) {
        op <- paste0(fop,
          c(rep.int("  ", ncol - 1L), "\n"),
          collapse = ""
        )
      }

      cat("", op, sep = "\n")
    } else {
      cat("", op, "", sep = "\n")
    }

    cat(gettext("\nType two or more numbers separated by spaces and then hit <Return> to continue. \n\n"))

    repeat {
      res <- tryCatch(

        scan("",
          what = 0,
          quiet = TRUE,
          nlines = 1
        ),
        error = identity
      )

      if (!inherits(res, "error") && length(res) >= 2L) {
        break
      }

      cat(red("\nERROR: Invalid selection. You must select at least 2 stratification variables.\n\n"))

      cat(gettext("Type two or more numbers separated by spaces and then hit <Return> to continue.\n\n"))
    }

    if (any(res == 0)) {
      return(character())
    }


    if (!is.null(preselect)) {
      res <- c(which(def), res)
    }

    res <- unique(res)
    res <- sort(res[1 <= res & res <= nc])

    return(choices[res])
  }
}

.round.preserve.sum <- function(x,
                                digits = 0) {
  up <- 10^digits
  x <- x * up
  y <- floor(x)
  indices <- tail(order(x - y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  return(y / up)
}

.make.var.overview <- function(dataset,
                               print_to_console = FALSE) {
  var_overview <- dataset %>%
    map_df(function(x) {
      tibble(
        Type = class(x),
        Levels = nlevels(x)
      )
    }) %>%
    mutate(Variable = names(dataset)) %>%
    select(Variable, everything()) %>%
    data.frame()

  var_overview %>%
    kbl(
      caption = "Variable Overview",
      align = "l",
      row.names = TRUE
    ) %>%
    kable_styling(c("striped", "hover"), fixed_thead = TRUE) %>%
    print()

  if (print_to_console == TRUE) {
    print(var_overview, row.names = FALSE)
  }
}

.make.cont.data.tbl <- function(cont_data) {
  cont_data_tbl <- cont_data %>%
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
