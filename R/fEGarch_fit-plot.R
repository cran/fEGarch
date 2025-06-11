

set_default <- function(obj, repl) {
  arg_names <- names(repl)
  for (i in seq_along(arg_names)) {
    if (is.null(obj[[arg_names[[i]]]])) {
      obj[[arg_names[[i]]]] <- repl[[i]]
    }
  }
  obj
}

check_which <- function(which) {

    if (is.null(which)) {
    text_prompt <- data.frame(
      c("(1) returns:", "(2) means:", "(3) stand_deviations:", "(4) residuals:"),
      c(
        "Input training series",
        "Fitted conditional means",
        "Fitted conditional standard deviations",
        "Standardized residuals"
      )
    )
    colnames(text_prompt) <- NULL

    cat("\nSelect one of the following plots via the keyword or the position number (exit with 0):\n")

    print.data.frame(text_prompt, row.names = FALSE, right = FALSE)

    cat("\n")

    which <- ""

    while(!(which %in% c(0, 1, 2, 3, 4, "returns", "means", "stand_deviations", "residuals"))) {

        which <- readline("Keyword or position number: ")

    }
    cat("\n")
    }

    which
}


plot_returns <- function(obj, ...) {
  dots <- list(...)
  defaults <- list(
    xlab = "Time",
    ylab = "Value",
    main = paste0('Input series'),
    type = "l"
  )
  dots <- set_default(dots, defaults)
  rt <- obj@rt
  if (inherits(obj@rt, "zoo") || inherits(obj@rt, "ts")) {
    tp <- c(time(rt))
    rt <- zoo::coredata(rt)
  } else {
    tp <- seq_along(rt)
  }
  dots[["y"]] <- rt
  dots[["x"]] <- tp

  do.call(plot, args = dots)
}

plot_cmeans <- function(obj, ...) {
  dots <- list(...)
  defaults <- list(
    xlab = "Time",
    ylab = "Value",
    main = paste0('Estimated conditional means'),
    type = "l"
  )
  dots <- set_default(dots, defaults)
  rt <- obj@cmeans
  if (inherits(obj@rt, "zoo") || inherits(obj@rt, "ts")) {
    tp <- c(time(rt))
    rt <- zoo::coredata(rt)
  } else {
    tp <- seq_along(rt)
  }
  dots[["y"]] <- rt
  dots[["x"]] <- tp

  do.call(plot, args = dots)
}

plot_sigt <- function(obj, ...) {
  dots <- list(...)
  defaults <- list(
    xlab = "Time",
    ylab = "Value",
    main = paste0('Estimated conditional standard deviations'),
    type = "l"
  )
  dots <- set_default(dots, defaults)
  rt <- obj@sigt
  if (inherits(obj@rt, "zoo") || inherits(obj@rt, "ts")) {
    tp <- c(time(rt))
    rt <- zoo::coredata(rt)
  } else {
    tp <- seq_along(rt)
  }
  dots[["y"]] <- rt
  dots[["x"]] <- tp

  do.call(plot, args = dots)
}

plot_etat <- function(obj, ...) {
  dots <- list(...)
  defaults <- list(
    xlab = "Time",
    ylab = "Value",
    main = paste0('Residuals'),
    type = "l"
  )
  dots <- set_default(dots, defaults)
  rt <- obj@etat
  if (inherits(obj@rt, "zoo") || inherits(obj@rt, "ts")) {
    tp <- c(time(rt))
    rt <- zoo::coredata(rt)
  } else {
    tp <- seq_along(rt)
  }
  dots[["y"]] <- rt
  dots[["x"]] <- tp

  do.call(plot, args = dots)
}

#'S4 Plot Generic
#'
#'Imported from base R.
#'
#'@param x see base R \code{plot} documentation.
#'@param y see base R \code{plot} documentation.
#'@param ... see base R \code{plot} documentation.
#'
#'@return
#'Returns nothing to the console but creates a plot in the plot window.
#'
#'@export
#'
setGeneric("plot")

### Base R plot method

#'Plot Method for Showing Fitting Step Results
#'
#'This is method for producing various plots of the estimation results
#'returned by this package.
#'
#'@param x an object returned by the fitting functions of this package,
#'for example by \code{\link{fEGarch}}.
#'@param y for compatibility but without use.
#'@param which various plots can be selected either via a keyword or a number;
#'enter \code{"returns"} or \code{1} to show a plot of the input training series;
#'enter \code{"means"} or \code{2} to show the
#'fitted conditional means; enter \code{"stand_deviations"} or \code{3} to show
#'the fitted conditional standard deviations; use
#'\code{"residuals"} or \code{4} to show the standardized residuals
#'following the fitted model; the
#'default is \code{which = NULL} which then lets you select a plot
#'interactively in the R console.
#'@param ... further arguments to pass to \code{\link[base]{plot}}.
#'
#'@details
#'Create predefined standard plots of the estimation objects returned by the
#'\code{fEGarch} package.
#'Plots are created in the base R plot style. The type of plot can be chosen
#'either interactively from the console, or the argument \code{which} can be
#'used to directly  select the kind of plot to create (see also the description
#'of the argument \code{which}) within the function call.
#'
#'@aliases plot,fEGarch_fit-method
#'
#'@export
#'
#'@return
#'A graphic is created in the plots windows, the function itself, however,
#'returns \code{NULL}.
#'
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Author and Package Creator
#'}
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'# Pure conditional volatility model
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'plot(model, which = 3)
#'

setMethod("plot", "fEGarch_fit", function(x, y = NULL, which = NULL, ...) {

  which <- check_which(which)

  if (which == "0") {
    return(invisible(NULL))
  }

  plot_fun <- switch(
    which,
    "returns" = plot_returns,
    "means" = plot_cmeans,
    "standard_deviations" = plot_sigt,
    "residuals" = plot_etat,
    "1" = plot_returns,
    "2" = plot_cmeans,
    "3" = plot_sigt,
    "4" = plot_etat
  )

  plot_fun(x, ...)

})

#############################################################

### ggplot2 plot method

plot_returns_gg <- function(obj, ...) {

  rt <- obj@rt

  df <- data.frame(
    Time = c(stats::time(rt)),
    Value = zoo::coredata(rt)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data[["Time"]], y = .data[["Value"]])) +
    ggplot2::geom_line() +
    ggplot2::ggtitle('Input series')
}

plot_cmeans_gg <- function(obj, ...) {

  rt <- obj@cmeans

  df <- data.frame(
    Time = c(stats::time(rt)),
    Value = zoo::coredata(rt)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data[["Time"]], y = .data[["Value"]])) +
    ggplot2::geom_line() +
    ggplot2::ggtitle('Estimated conditional means')
}

plot_sigt_gg <- function(obj, ...) {

  rt <- obj@sigt

  df <- data.frame(
    Time = c(stats::time(rt)),
    Value = zoo::coredata(rt)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data[["Time"]], y = .data[["Value"]])) +
    ggplot2::geom_line() +
    ggplot2::ggtitle('Estimated conditional standard deviations')
}

plot_etat_gg <- function(obj, ...) {

  rt <- obj@etat

  df <- data.frame(
    Time = c(stats::time(rt)),
    Value = zoo::coredata(rt)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data[["Time"]], y = .data[["Value"]])) +
    ggplot2::geom_line() +
    ggplot2::ggtitle('Residuals')
}

setGeneric(
  "autoplot",
  function(object, ...) standardGeneric("autoplot"),
  useAsDefault = ggplot2::autoplot
)

#'Plot Method for Fitting Step Results in the Style of ggplot2
#'
#'This is method for producing various plots of the decomposition results
#'returned by this package.
#'
#'@param object an object returned by the fitting functions of this package,
#'for example by \code{\link{fEGarch}}.
#'@param which various plots can be selected either via a keyword or a number;
#'enter \code{"returns"} or \code{1} to show a plot of the input training series;
#'enter \code{"means"} or \code{2} to show the
#'fitted conditional means; enter \code{"stand_deviations"} or \code{3} to show
#'the fitted conditional standard deviations;
#'use \code{"residuals"} or \code{4} to show the standardized residuals
#'following the fitted model;
#'the default is \code{which = NULL} which then lets you select a plot
#'interactively in the R console.
#'@param ... no purpose and only implemented for compatibility.
#'
#'@details
#'Create predefined standard plots of the estimation objects returned by the
#'\code{fEGarch} package.
#'Plots are created in the ggplot2 plot style. The type of plot can be chosen
#'either interactively from the console, or the argument \code{which} can be
#'used to directly  select the kind of plot to create (see also the description
#'of the argument \code{which}) within the function call.
#'
#'@export
#'@importFrom ggplot2 autoplot
#'
#'@aliases autoplot,fEGarch_fit-method
#'
#'@return
#'A ggplot2-graphic object is returned, i.e. an object of classes
#'\code{"gg"} and \code{"ggplot"}.
#'
#'@importFrom rlang .data
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'Author and Package Creator
#'}
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'# Pure conditional volatility model
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'autoplot(model, which = 3)
#'

setMethod("autoplot", "fEGarch_fit", function(object, which = NULL, ...) {

  which <- check_which(which)

  if (which == "0") {
    return(invisible(NULL))
  }

  plot_fun <- switch(
    which,
    "returns" = plot_returns_gg,
    "means" = plot_cmeans_gg,
    "standard_deviations" = plot_sigt_gg,
    "residuals" = plot_etat_gg,
    "1" = plot_returns_gg,
    "2" = plot_cmeans_gg,
    "3" = plot_sigt_gg,
    "4" = plot_etat_gg
  )

  P <- plot_fun(object, ...)
  P

})

