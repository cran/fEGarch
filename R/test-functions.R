#'Goodness-of-Fit Test Generics
#'
#'Generics for different goodness-of-fit tests. Currently
#'without use.
#'
#'@param object currently without use.
#'@param m_max currently without use.
#'@param weight_f currently without use.
#'@param silent currently without use.
#'@param n_bins currently without use.
#'@param adj_df currently without use.
#'@param args_lbt currently without use.
#'@param args_goft currently without use.
#'@param args_sbt currently without use.
#'@param ... currently without use.
#'
#'@return
#'The generics are currently without use and do not
#'return anything.
#'
#'@export
#'@rdname fit-test-generics
#'
setGeneric("ljung_box_test", function(object, m_max = 20, weight_f = function(lag) {m <- max(lag); (m + 1 - lag) / m}, adj_df = NULL, silent = FALSE, ...) {standardGeneric("ljung_box_test")})


#'Weighted Ljung-Box Test for Autocorrelation
#'
#'Apply a (weighted) Ljung-Box test (through the
#'Gamma approximation) to check
#'the standardized residuals of a fitted model
#'from this package for remaining autocorrelation.
#'Two different options allow to check either the
#'simple residuals or the squared residuals.
#'
#'@param object an object \code{"fEGarch_fit"} as returned
#'by the fitting functions of this package, for example
#'by \code{\link{fEGarch}}.
#'@param m_max the maximum lag; tests will be conducted for 1 up to
#'\code{m_max}.
#'@param weight_f a function with argument \code{lag} stating how
#'weights should be calculated.
#'@param silent a logical value reflecting whether or not test results
#'should be printed in a well-formatted manner to the console.
#'@param type either \code{"simple"} or \code{"squared"} for applying
#'the test to simple or squared residuals.
#'@param adj_df degrees of freedom to adjust for as a number or the default
#'\code{NULL}, which uses automatic values from the fitted object; for
#'squared residuals \code{adj_df = 0} is the default and for simple
#'returns, it is the sum of ARMA-parameters.
#'@param ... currently without use.
#'
#'@export
#'
#'@return
#'Returns a numeric matrix invisibly.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'ljung_box_test(model)
#'
setMethod("ljung_box_test", "fEGarch_fit", function(object, m_max = 20, weight_f = function(lag) {m <- max(lag); (m + 1 - lag) / m}, adj_df = NULL, silent = FALSE, type = c("simple", "squared"), ...) {

  type <- match.arg(type)
  et <- zoo::coredata(object@etat)

  fun <- list(
    "simple" = function(et) {et},
    "squared" = function(et) {et^2}
  )[[type]]

  adj_df <- if (!is.null(adj_df)) {
    adj_df
  } else if (is.null(adj_df) && type == "squared") {
    0
  } else if (is.null(adj_df) && type == "simple") {
    sum(object@meanspec@orders)
  }

  et_transf <- fun(et)

  cors <- c(stats::acf(et_transf, plot = FALSE, type = "correlation", lag.max = m_max)$acf)[-1]
  n <- length(et)

  stats <- rep(NA, m_max)

  lagg <- 1:m_max

  for (i in 1:m_max) {
    w <- weight_f(lag = 1:i)
    stats[[i]] <- (n * (n + 2)) * sum((cors[1:i]^2 * w) / (n - (1:i)))
  }

  shape <- (3/4) * (lagg + 1)^2 * lagg / (2 * lagg^2 + 3 * lagg +
            1 - 6 * lagg * adj_df)
  scale <- (2/3) * (2 * lagg^2 + 3 * lagg + 1 - 6 * lagg *
            adj_df) / lagg / (lagg + 1)

  p_vals <- suppressWarnings(1 - stats::pgamma(stats, shape = shape, scale = scale))
  if (any(is.nan(p_vals))) {
    p_vals[is.nan(p_vals)] <- NA
  }

  m <- 1:m_max

  df <- matrix(NA_real_, nrow = m_max, ncol = 2)
  df[, 1] <- m
  df[, 2] <- suppressWarnings(p_vals)
  rownames(df) <- rep("", m_max)
  colnames(df) <- c("m", "pval")

  if (!silent) {
    out1 <- paste0(
      "\n************************************************************\n",
      "* Weighted Ljung-Box test for autocorrelation in residuals *\n",
      "************************************************************\n\n",
      "H0: autocorrelations for lags 1 to m simultaneously equal to zero\n",
      "Mode: ", type, "\n"
    )
    cat(out1, fill = TRUE)
    df2 <- df
    df2[, 2] <- suppressWarnings(as.numeric(sprintf("%.4f", p_vals)))
    print(df2)
  }
  invisible(df)

})


#'@export
#'@rdname fit-test-generics
setGeneric("sign_bias_test", function(object, silent = FALSE, ...) {standardGeneric("sign_bias_test")})

#'Sign Bias Test
#'
#'Apply a sign bias test to check
#'the standardized residuals of a fitted model
#'from this package for remaining significant
#'sign effects.
#'
#'@param object an object \code{"fEGarch_fit"} as returned
#'by the fitting functions of this package, for example
#'by \code{\link{fEGarch}}.
#'@param silent a logical value reflecting whether or not test results
#'should be printed in a well-formatted manner to the console.
#'@param ... currently without use.
#'
#'@export
#'
#'@return
#'Returns a numeric matrix invisibly.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'sign_bias_test(model)
#'
setMethod("sign_bias_test", "fEGarch_fit", function(object, silent = FALSE, ...) {

  scale_f <- if (is.null(object@scale_fun)) {
    1
  } else {
    object@scale_fun
  }

  dm <- zoo::coredata((object@rt - object@cmeans) / scale_f)   # obs adjusted for unconditional and cond. mean
  et <- tail(zoo::coredata(object@etat), -1)    # standardized model residuals (without first)
  dm_s <- utils::head(dm, -1)      # dm without last

  # Indicators
  I_small <- as.numeric(dm_s < 0)
  I_large <- as.numeric(dm_s >= 0)

  V1 <- I_small
  V2 <- I_small * dm_s
  V3 <- I_large * dm_s
  Y <- et^2

  n <- length(Y)

  reg <- stats::lm(Y ~ V1 + V2 + V3)

  sum_stats <- summary(reg)
  R2 <- sum_stats$r.squared

  pvals_reg <- unname(sum_stats$coefficients[, 4, drop = TRUE][-1])

  p_joint <- 1 - stats::pchisq(n * R2, df = 3)
  p_all <- c(pvals_reg, p_joint)

  mat <- matrix(NA, ncol = 1, nrow = 4)
  rownames(mat) <- c(
    "Sign bias",
    "Negative sign bias",
    "Positive sign bias",
    "Joint"
  )
  mat[, 1] <- p_all
  colnames(mat) <- "pval"

  if (!silent) {
    out1 <- paste0(
      "\n******************\n",
      "* Sign bias test *\n",
      "******************\n\n",
      "H0 (sign bias):      negative shocks have no significant effect\n",
      "                     on volatility beyond symmetric model\n",
      "H0 (neg. sign bias): magnitude of negative shocks has no\n",
      "                     additional impact on volatility\n",
      "H0 (pos. sign bias): magnitude of positive shocks has no\n",
      "                     additional impact on volatility\n",
      "H0 (joint):          all three of the above simultaneously\n"
    )
    cat(out1, fill = TRUE)
    mat2 <- mat
    mat2[, 1] <- as.numeric(sprintf("%.4f", p_all))
    print(mat2)
  }
  invisible(mat)

})

#'@export
#'@rdname fit-test-generics
setGeneric("goodn_of_fit_test", function(object, n_bins = c(20, 30, 40, 50), silent = FALSE, ...) {standardGeneric("goodn_of_fit_test")})


#'Adjusted Pearson Goodness-of-Fit Test for Standardized Model Residuals
#'
#'Consider a probability integral transform on the
#'standardized residuals of a fitted model and apply a
#'Pearson goodness-of-fit test to the binned data using a
#'selection of predefined number of bins.
#'
#'@param object an object of class \code{"fEGarch_fit"} as returned
#'by the fitting functions of this package like for example
#'\code{\link{fEGarch}}.
#'@param n_bins a numeric vector giving the number of bins to use.
#'@param silent a logical indicating whether or not to print
#'the test results in a well-formatted manner to the console.
#'@param ... currently without purpose.
#'
#'@details
#'Use a probability integral transform on the
#'standardized residuals of a fitted model. This is then the
#'basis to conduct a Pearson goodness-of-fit chi-square test.
#'
#'@return
#'Returns a numeric matrix invisibly.
#'
#'@export
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'goodn_of_fit_test(model)
#'
setMethod("goodn_of_fit_test", "fEGarch_fit", function(object, n_bins = c(20, 30, 40, 50), silent = FALSE, ...) {

  model <- object
  pars <- model@pars
  cond_d <- model@cond_dist
  n_par <- length(pars)

  dfun <- fun1_selector(cond_d)

  args <- list(
    dfun = dfun
  )

  args[["shape"]] <- if (!(cond_d %in% c("norm", "snorm"))) {
    par_name <- switch(
      cond_d,
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    pars[[par_name]]
  } else {
    0
  }

  args[["skew"]] <- if (cond_d %in% c("snorm", "sstd", "sged", "sald")) {
    pars[["skew"]]
  } else {
    0
  }

  lb <- length(n_bins)
  out_mat <- matrix(NA, nrow = lb, ncol = 3)
  rownames(out_mat) <- rep("", lb)
  colnames(out_mat) <- c("n_bins", "df", "pval")
  et <- zoo::coredata(model@etat)
  n <- length(et)
  #adjust <- n_par * adj

  for (i in 1:lb) {

    nb <- n_bins[[i]]
    one_div_nb <- 1 / nb
    p <- seq(from = one_div_nb, to = 1 - one_div_nb, by = one_div_nb)

    args[["p"]] <- p

    quantiles <- do.call(inv_cdf_fun, args = args)         # Get quantiles corresponding to class bounds

    np1 <- length(p) + 1
    frequ_et <- rep(NA, np1)
    q_adj <- c(-Inf, quantiles, Inf)

    for (j in seq_along(frequ_et)) {
      frequ_et[[j]] <- sum(et >= q_adj[[j]] & et < q_adj[[j + 1]])
    }

    expec <- n * one_div_nb   # expected frequencies in each class

    stat <- sum((frequ_et - expec)^2 / expec)   # test statistic

    df <- nb - 1# - adjust

    pval <- if (df >= 1) {
      1 - stats::pchisq(stat, df = df)
    } else {
      NA
    }

    out_mat[i, ] <- c(
      nb,
      df,
      pval
    )

  }



  if (!silent) {
    out1 <- paste0(
    "\n************************************\n",
      "* Adjusted Pearson goodness-of-fit *\n",
      "*         test on residuals        *\n",
      "************************************\n\n",
      "H0: standardized model residuals align with\n",
      "    assumed innovation distribution\n"
    )
    cat(out1, fill = TRUE)
    mat2 <- out_mat
    mat2[, 3] <- as.numeric(sprintf("%.4f", out_mat[, 3]))
    print(mat2)
  }

  invisible(out_mat)
})

#'@export
#'@rdname fit-test-generics
setGeneric("fit_test_suite", function(object, args_lbt = list(), args_sbt = list(), args_goft = list(), silent = FALSE, ...) {standardGeneric("fit_test_suite")})


#'Post-Estimation Fit-Tests
#'
#'Apply a collection of fit-tests, including
#'a weighted Ljung-Box test for the simple and the squared
#'standardized residuals, a sign-bias test, and an
#'adjusted Pearson goodness-of-fit test.
#'
#'@param object an object of class \code{"fEGarch_fit"} as returned
#'by the fitting functions of this package like for example
#'\code{\link{fEGarch}}.
#'@param args_lbt a list of changes to make to the default argument settings in
#'\code{ljung_box_test}.
#'@param args_sbt a list of changes to make to the default argument settings in
#'\code{sign_bias_test}.
#'@param args_goft a list of changes to make to the default argument settings in
#'\code{goodn_of_fit_test}.
#'@param silent a logical indicating whether or not to print
#'the test results in a well-formatted manner to the console.
#'@param ... currently without purpose.
#'
#'@export
#'
#'@return
#'Returns a list with the four test results invisibly.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'fit_test_suite(model)
#'

setMethod("fit_test_suite", "fEGarch_fit", function(object, args_lbt = list(),
                                                    args_sbt = list(),
                                                    args_goft = list(),
                                                    silent = FALSE, ...) {

  fun1 <- methods::selectMethod("ljung_box_test", "fEGarch_fit")
  fun2 <- methods::selectMethod("sign_bias_test", "fEGarch_fit")
  fun3 <- methods::selectMethod("goodn_of_fit_test", "fEGarch_fit")

  default_lbt = list(
    m_max = 20,
    weight_f = function(lag) {m <- max(lag); (m + 1 - lag) / m},
    adj_df = NULL
  )

  default_sbt = list()

  default_goft = list(
    n_bins = c(20, 30, 40, 50)
  )


  n1 <- length(args_lbt)
  if (n1 > 0) {
    names_lbt <- names(args_lbt)
    for(i in names_lbt) {
      default_lbt[[i]] <- args_lbt[[i]]
    }
  }
  default_lbt[["object"]] <- object
  default_lbt[["silent"]] <- silent

  n2 <- length(args_sbt)
  if (n2 > 0) {
    names_sbt <- names(args_sbt)
    for(i in names_sbt) {
      default_sbt[[i]] <- args_sbt[[i]]
    }
  }
  default_sbt[["object"]] <- object
  default_sbt[["silent"]] <- silent

  n3 <- length(args_goft)
  if (n3 > 0) {
    names_goft <- names(args_goft)
    for(i in names_goft) {
      default_goft[[i]] <- args_goft[[i]]
    }
  }
  default_goft[["object"]] <- object
  default_goft[["silent"]] <- silent

  df1.2 <- default_lbt
  default_lbt[["type"]] <- "simple"
  df1.2[["type"]] <- "squared"

  result1 <- do.call(fun1, args = default_lbt)
  result1.2 <- do.call(fun1, args = df1.2)
  result2 <- do.call(fun2, args = default_sbt)
  result3 <- do.call(fun3, args = default_goft)

  invisible(list(
    "LBT_simple" = result1,
    "LBT_squared" = result1.2,
    "SBT" = result2,
    "GOFT" = result3
  ))

})
