
setClass("fEGarch_risk",
  slots = c(
    measures = "list",
    observations = "ANY",
    cmeans = "ANY",
    sigt = "ANY",
    model = "ANY"
  ),
  prototype = list(
    measures = list(1),
    observations = numeric(1),
    cmeans = numeric(1),
    sigt = numeric(1),
    model = numeric(1)
  )
)

fEGarch_risk <- function(
    measures,
    observations,
    cmeans,
    sigt,
    model
) {

  new("fEGarch_risk",
    measures = measures,
    observations = observations,
    cmeans = cmeans,
    sigt = sigt,
    model = model
  )

}

PoT_calc <- function(object, level = NULL) {
  VaR <- object@measures$VaR
  VaR_names <- if (is.null(level)) {
    names(VaR)
  } else {
    paste0("VaR", level)
  }
  out <- vector(mode = "list", length = length(VaR_names))
  names_adj <- paste0("PoT_", VaR_names)
  names(out) <- names_adj
  for (i in seq_along(VaR_names)) {
    selector <- object@observations < VaR[[VaR_names[[i]]]]
    out[[names_adj[[i]]]] <- object@observations[selector]
  }
  out
}

check_which_risk <- function(which, object) {

    if (is.null(which)) {

    VaRs <- if (!is.null(object@measures$VaR)) {
      paste0(substring(names(object@measures$VaR), 4), collapse = ", ")
    } else {
      ""
    }

    cat("\nThe following VaR levels were found (exit with 0):\n")
    cat(paste0("\n", VaRs, "\n"))

    cat("\n")

    which <- readline("Enter level to display: ")

    cat("\n")
    }

    which

}

#'Plotting of Risk Measure Results (Base R)
#'
#'Plot risk measure results returned by \code{measure_risk}
#'as a points-over-threshold plot in style of base R plots.
#'
#'@param x an object returned by \code{measure_risk}.
#'@param y for compatibility but without use.
#'@param which one of the levels of VaR and ES saved in
#'\code{object}, usually either \code{0.975} or \code{0.99}
#'by default.
#'@param ... without use.
#'
#'@return
#'Returns nothing but produces a base R plot in the plot window.
#'
#'@aliases plot,fEGarch_risk-method
#'
#'@export
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2003-12-31")
#'
#'egarch_spec() %>%
#'  fEGarch(rt = rt, n_test = 250) %>%
#'  predict_roll() %>%
#'  measure_risk() %>%
#'  plot(which = 0.99)
#'
#'
setMethod("plot", "fEGarch_risk", function(x, y = NULL, which = NULL, ...) {

  which <- check_which_risk(which, x)

  VaR_name <- paste0("VaR", which)
  ES_name <- paste0("ES", which)

  VaR <- x@measures$VaR[[VaR_name]]
  ES <- x@measures$ES[[ES_name]]
  obs <- x@observations

  ES_check <- is.null(ES)

  if (which == "0") {
    return(invisible(NULL))
  }

  all_series <- cbind(
    obs, VaR, ES
  )

  tp <- c(time(obs))
  dim_s <- dim(all_series)
  all_series <- matrix(
    as.numeric(all_series),
    nrow = dim_s[[1]],
    ncol = dim_s[[2]]
  )

  input <- list(...)
  if (is.null(input[["type"]])) {
    input[["type"]] <- list("hll", "hl")[[ES_check + 1]]
  }
  if (is.null(input[["lty"]])) {
    input[["lty"]] <- list(c(1, 1, 1), c(1, 1))[[ES_check + 1]]
  }
  if (is.null(input[["col"]])) {
    input[["col"]] <- list(c(1, "gold2", 4), c(1, "gold2"))[[ES_check + 1]]
  }
  if (is.null(input[["xlab"]])) {
    input[["xlab"]] <- "Time"
  }
  if (is.null(input[["ylab"]])) {
    input[["ylab"]] <- "Return, VaR and ES"
  }
  if (is.null(input[["main"]])) {
    input[["main"]] <- "Backtesting results"
  }
  if (is.null(input[["xaxt"]]) && inherits(obs, "zoo")) {
    setting_o <- NULL
    input[["xaxt"]] <- "n"
  } else {
    setting_o <- input[["xaxt"]]
  }

  input[["x"]] <- tp
  input[["y"]] <- all_series

  do.call(graphics::matplot, args = input)

  if (is.null(setting_o) && inherits(obs, "zoo")) {
    dates <- time(obs)
    n <- length(dates)
    month_start <- seq(
      from = as.Date(paste0(format(dates[[1]], "%Y-%m"), "-01")),
      by = "-1 month",
      length.out = 2
    )[[2]]
    month_last <- seq(
      from = as.Date(paste0(format(dates[[n]], "%Y-%m"), "-01")),
      by = "1 month",
      length.out = 2
    )[[2]]
    months <- seq(
      from = month_start,
      to = month_last,
      by = "month"
    )
    graphics::axis(side = 1, at = months, labels = format(months, "%b %y"))
  }

  graphics::abline(h = 0, lty = 3, col = "lightgray")

  idx <- obs < VaR

  obs_s1 <- as.numeric(obs)[idx]
  tp2 <- c(time(VaR))[idx]
  if (sum(idx) >= 1) {
    graphics::lines(tp2, obs_s1, type = "h", col = "red", lwd = 1.5)
  }

  idx <- obs < ES

  obs_s2 <- as.numeric(obs)[idx]
  tp2 <- c(time(ES))[idx]
  if (sum(idx) >= 1) {
    graphics::lines(tp2, obs_s2, type = "h", col = "purple", lwd = 1.5)
  }


})

#'Plotting of Risk Measure Results (\code{ggplot2})
#'
#'Plot risk measure results returned by \code{measure_risk}
#'as a points-over-threshold plot in style of \code{ggplot2}.
#'
#'@param object an object returned by \code{measure_risk}.
#'@param which one of the levels of VaR and ES saved in
#'\code{object}, usually either \code{0.975} or \code{0.99}
#'by default.
#'@param ... without use.
#'
#'@return
#'Returns a \code{ggplot2} plot object.
#'
#'@export
#'
#'@aliases autoplot,fEGarch_risk-method
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2003-12-31")
#'
#'egarch_spec() %>%
#'  fEGarch(rt = rt, n_test = 250) %>%
#'  predict_roll() %>%
#'  measure_risk() %>%
#'  autoplot(which = 0.99)
#'
#'
setMethod("autoplot", "fEGarch_risk", function(object, which = NULL, ...) {

  which <- check_which_risk(which, object)

  VaR_name <- paste0("VaR", which)
  ES_name <- paste0("ES", which)

  VaR <- object@measures$VaR[[VaR_name]]
  ES <- object@measures$ES[[ES_name]]
  obs <- object@observations

  ES_check <- is.null(ES)

  if (which == "0") {
    return(invisible(NULL))
  }

  all_series <- cbind(
    Obs = zoo::coredata(obs),
    VaR = zoo::coredata(VaR),
    ES = zoo::coredata(ES)
  ) %>%
    as.data.frame()

  all_series[["Time"]] <- c(time(obs))

  out <- ggplot2::ggplot(all_series, ggplot2::aes(x = .data[["Time"]], y = .data[["Obs"]])) +
    ggplot2::geom_hline(yintercept = 0, color = "lightgray", linetype = 3) +
    ggplot2::geom_segment(ggplot2::aes(x = .data[["Time"]], y = .data[["Obs"]], yend = 0, color = "1")) +
    ggplot2::geom_line(inherit.aes = FALSE, ggplot2::aes(x = .data[["Time"]], y = .data[["VaR"]], color = "2")) +
    ggplot2::geom_line(inherit.aes = FALSE, ggplot2::aes(x = .data[["Time"]], y = .data[["ES"]], color = "3")) +
    ggplot2::scale_color_manual(name = "Series", values = c(
      "1" = "black", "2" = "gold2", "3" = 4
    ), labels = c(
      "1" = "Returns",
      "2" = paste0(100 * as.numeric(which), "%-VaR"),
      "3" = paste0(100 * as.numeric(which), "%-ES")
    )) +
    ggplot2::ylab("Return, VaR and ES") +
    ggplot2::ggtitle("Backtesting results")

if (inherits(obs, what = "zoo")) {
    dates <- all_series$Time
    n <- length(dates)
    month_start <- seq(
      from = as.Date(paste0(format(dates[[1]], "%Y-%m"), "-01")),
      by = "-1 month",
      length.out = 2
    )[[2]]
    month_last <- seq(
      from = as.Date(paste0(format(dates[[n]], "%Y-%m"), "-01")),
      by = "1 month",
      length.out = 2
    )[[2]]
    months <- seq(
      from = month_start,
      to = month_last,
      by = "month"
    )
    out <- out +
      ggplot2::scale_x_date(date_labels = "%b %y",
                          breaks = months)
  }

  idx <- obs < VaR

  obs2 <- as.numeric(obs)[idx]
  tp2 <- c(time(VaR))[idx]

  if (length(idx) > 0) {

    df <- data.frame(
      Time = tp2,
      Obs = obs2
    )

    out <- out +
      ggplot2::geom_segment(
        data = df,
        mapping = ggplot2::aes(x = .data[["Time"]], y = .data[["Obs"]], yend = 0),
        color = "red", linewidth = 0.6,
        inherit.aes = FALSE
      )

  }

  idx <- obs < ES

  obs3 <- as.numeric(obs)[idx]
  tp2 <- c(time(ES))[idx]

  if (length(idx) > 0) {

    df <- data.frame(
      Time = tp2,
      Obs = obs3
    )

    out <- out +
      ggplot2::geom_segment(
        data = df,
        mapping = ggplot2::aes(x = .data[["Time"]], y = .data[["Obs"]], yend = 0),
        color = "purple", linewidth = 0.6,
        inherit.aes = FALSE
      )

  }

  out

})
