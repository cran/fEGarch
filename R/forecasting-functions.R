fetch_arma_pars_fit <- function(object) {

  lm_arma <-  object@meanspec@long_memo
  arma_orders <- object@meanspec@orders

  p_ar <- arma_orders[[1]]
  q_ma <- arma_orders[[2]]

  pg0 <- p_ar > 0
  qg0 <- q_ma > 0
  pqg0 <- p_ar + q_ma > 0

  pars <- object@pars

  # Get all estimated parameters from output

  ar <- if (pg0) {
    names_ar <- paste0("ar", 1:p_ar)
    pars[names_ar]
  } else {
    numeric(0)
  }

  ma <- if (qg0) {
    names_ma <- paste0("ma", 1:q_ma)
    pars[names_ma]
  } else {
    numeric(0)
  }

  D <- if (lm_arma) pars[["D"]] else numeric(0)

  list(
   ar = ar,
   ma = ma,
   D = D,
   pqg0 = pqg0,
   lm_arma = lm_arma
  )

}

arma_farima_forecasts <- function(ar, ma, D, pqg0, lm_arma, mu, x, innov, n.ahead, trunc) {

  mean_fc <- if (pqg0 || lm_arma) {

    if (!lm_arma) {
      c(forecast_arma_Cpp(
        x = x,
        et = innov,
        mu = mu,
        ma_e = ma,
        ar_e = ar,
        horizon = n.ahead
      ))
    } else {
      # Possibly truncate infinite coefficient series
      coef_inf_arma <- ar_infty(ar = ar, ma = ma, d = D, max_i = trunc)[-1]
      c(forecast_farima_Cpp(
        x = x,
        mu = mu,
        horizon = n.ahead,
        coef_inf = coef_inf_arma
      ))
    }

  } else {

    rep(mu, n.ahead)

  }

  mean_fc
}

# Forecasting methods (multistep)

# For EGARCH-type
#Prediction Methods for Package's Models
#
#Produces forecasts of the conditional standard deviation
#(and of the conditional mean) following package's models.
#These methods are not being exported and are used internally
#only.
#
#@param object these methods are not being exported.
#@param n.ahead these methods are not being exported.
#@param trunc these methods are not being exported.
#@param ... these methods are not being exported.
#
#@return
#They return lists with two numeric vectors as elements named
#\code{sigt} and \code{cmeans}.
#
#@rdname forecasting-model-methods
#
methods::setMethod("fEGarch_predict", "fEGarch_fit_egarch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

  lm_garch <- object@long_memo
  garch_orders <- object@orders
  p <- garch_orders[[1]]
  q <- garch_orders[[2]]

  model_pars <- object@pars
  names_phi <- paste0("phi", 1:p)
  phi <- model_pars[names_phi]

  psi <- if (q > 1) {
    names_psi <- paste0("psi", 1:(q - 1))
    model_pars[names_psi]
  } else {
    numeric(0)
  }
  psi_2 <- psi
  psi <- c(1, psi)
  kappa <- model_pars[["kappa"]]
  gamma <- model_pars[["gamma"]]

  d <- if (lm_garch) model_pars[["d"]] else numeric(0)

  cond_dist <- object@cond_dist

  skew <- if (cond_dist %in% c("snorm", "sstd", "sged", "sald")) {
    model_pars[["skew"]]
  } else {
    numeric(0)
  }

  shape <- if (cond_dist %in% c("std", "sstd", "ged", "sged", "ald", "sald")) {
    sel_word <- switch(
      cond_dist,
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    model_pars[[sel_word]]
  } else  {
    numeric(0)
  }

  mu <- if (exists("mu", where = as.list(model_pars))) {
    if (!is.null(object@nonpar_model)) {
      0
    } else {
      model_pars[["mu"]]
    }
  } else {
    0
  }

  # ARMA / FARIMA specifications
  mean_specs <- fetch_arma_pars_fit(object)
  ar <- mean_specs$ar
  ma <- mean_specs$ma
  D <- mean_specs$D
  pqg0 <- mean_specs$pqg0
  lm_arma <- mean_specs$lm_arma

  # EGARCH settings
  powers <- object@powers
  modulus <- object@modulus

  dfun <- fun1_selector(cond_dist)

  pow1 <- powers[[1]]
  pow1g0 <- pow1 > 0
  mod1 <- modulus[[1]]
  fun_asy <- exp_funs[["asy"]][[mod1 + 1]][[pow1g0 + 1]]
  E_asy <- stats::integrate(fun_asy, lower = -Inf, upper = Inf,
                            subdivisions = 1000L, rel.tol = 1e-6,
                            abs.tol = 1e-6,
                            pow = pow1, shape = shape, skew = skew, dfun = dfun)[[1]]

  pow2 <- powers[[2]]
  pow2g0 <- pow2 > 0
  mod2 <- modulus[[2]]
  fun_mag <- exp_funs[["mag"]][[mod2 + 1]][[pow2g0 + 1]]
  E_mag <- stats::integrate(fun_mag, lower = -Inf, upper = Inf,
                            subdivisions = 1000L, rel.tol = 1e-6,
                            abs.tol = 1e-6,
                            pow = pow2, shape = shape, skew = skew, dfun = dfun)[[1]]

  Elnsig2 <- model_pars[["omega_sig"]]

  # ARMA / FARIMA point forecasts of cond. mean
  mean_fc <- arma_farima_forecasts(ar = ar, ma = ma, D = D, pqg0 = pqg0,
                                   lm_arma = lm_arma, mu = mu,
                                   x = zoo::coredata(object@rt),
                                   innov = zoo::coredata(object@rt - object@cmeans),
                                   n.ahead = n.ahead, trunc = trunc)

  # EGARCH point forecasts of cond. volatility

  mode <- mode_mat[pow1g0 + 1, pow2g0 + 1]

  sigt_fc <- if (!lm_garch) {
    alpha <- psi * kappa
    beta <- psi * gamma
    omega <- Elnsig2 * (1 - sum(phi))
    c(sigt_egarch_forecast_shortCpp(
      et = zoo::coredata(object@etat),
      sigt = zoo::coredata(object@sigt),
      omega = omega,
      phi = phi,
      alpha = alpha,
      beta = beta,
      E_asy = E_asy,
      E_mag = E_mag,
      powers = powers,
      modulus = modulus,
      mode = mode,
      horizon = n.ahead
    ))
  } else {
    # Use MA-representation
    coef_inf_garch <- ma_infty(
      ar = phi, ma = psi_2, d = d, max_i = trunc - 1
    )

    c(sigt_egarch_forecast_longCpp(
      et = zoo::coredata(object@etat),
      coef_inf = coef_inf_garch,
      kappa = kappa,
      gamma = gamma,
      E_asy = E_asy,
      E_mag = E_mag,
      Elnsig2 = Elnsig2,
      powers = powers,
      modulus = modulus,
      mode = mode,
      horizon = n.ahead
    ))
  }

  list(
    cmeans = mean_fc,
    sigt = sigt_fc
  )

  }
)

# For Log-GARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_loggarch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

  lm_garch <- object@long_memo
  garch_orders <- object@orders
  p <- garch_orders[[1]]
  q <- garch_orders[[2]]
  n <- length(object@rt)

  model_pars <- object@pars
  names_phi <- paste0("phi", 1:p)
  phi <- model_pars[names_phi]

  names_psi <- paste0("psi", 1:q)
  psi <- model_pars[names_psi]

  d <- if (lm_garch) model_pars[["d"]] else numeric(0)

  cond_dist <- object@cond_dist

  skew <- if (cond_dist %in% c("snorm", "sstd", "sged", "sald")) {
    model_pars[["skew"]]
  } else {
    numeric(0)
  }

  shape <- if (cond_dist %in% c("std", "sstd", "ged", "sged", "ald", "sald")) {
    sel_word <- switch(
      cond_dist,
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    model_pars[[sel_word]]
  } else  {
    numeric(0)
  }

  mu <- if (exists("mu", where = as.list(model_pars))) {
    if (!is.null(object@nonpar_model)) {
      0
    } else {
      model_pars[["mu"]]
    }
  } else {
    0
  }

  # ARMA / FARIMA specifications
  mean_specs <- fetch_arma_pars_fit(object)
  ar <- mean_specs$ar
  ma <- mean_specs$ma
  D <- mean_specs$D
  pqg0 <- mean_specs$pqg0
  lm_arma <- mean_specs$lm_arma

  dfun <- fun1_selector(cond_dist)

  int_fun <- function(x, shape, skew, dfun) {
    log(x^2) * dfun(x, shape = shape, skew = skew)
  }

  Elneta2 <- stats::integrate(
    int_fun, lower = -Inf, upper = Inf,
    subdivisions = 1000L, rel.tol = 1e-6, abs.tol = 1e-6,
    shape = shape, skew = skew, dfun = dfun
  )[[1]]

  Elnsig2 <- model_pars[["omega_sig"]]

  # ARMA / FARIMA point forecasts of cond. mean
  mean_fc <- arma_farima_forecasts(ar = ar, ma = ma, D = D, pqg0 = pqg0,
                                   lm_arma = lm_arma, mu = mu,
                                   x = zoo::coredata(object@rt),
                                   innov = zoo::coredata(object@rt - object@cmeans),
                                   n.ahead = n.ahead, trunc = trunc)

  # Log-GARCH point forecasts of cond. volatility

  sigt_fc <- if (!lm_garch) {
    omega <- Elnsig2 * (1 - sum(phi))
    l <- max(p, q)
    phi_s <- psi_s <- rep(0, l)
    phi_s[1:p] <- phi
    psi_s[1:q] <- psi
    alpha <- phi_s + psi_s
    c(sigt_loggarch_forecast_short(
      et = zoo::coredata(object@etat),
      sigt = zoo::coredata(object@sigt),
      omega = omega,
      phi = phi,
      psi = alpha,
      Elneta2 = Elneta2,
      horizon = n.ahead
    ))
  } else {
    coef_inf_garch <- ma_infty(ar = phi, ma = psi, d = d,
                                     max_i = trunc)[-1]
    c(sigt_loggarch_forecast_long(
      et = zoo::coredata(object@etat),
      coef_inf = coef_inf_garch,
      Elneta2 = Elneta2,
      Elnsig2 = Elnsig2,
      horizon = n.ahead
    ))
  }

  list(
    cmeans = mean_fc,
    sigt = sigt_fc
  )

  }
)

# For APARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_aparch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

  lm_garch <- object@long_memo
  garch_orders <- object@orders
  p <- garch_orders[[1]]
  q <- garch_orders[[2]]
  n <- length(object@rt)

  model_pars <- object@pars
  names_phi <- paste0("phi", 1:p)
  phi <- model_pars[names_phi]

  names_beta <- paste0("beta", 1:q)
  beta <- model_pars[names_beta]

  names_gamma <- paste0("gamma", 1:p)
  gamma <- model_pars[names_gamma]

  delta <- model_pars[["delta"]]

  omega <- model_pars[["omega"]]

  #d <- if (lm_garch) model_pars[["d"]] else numeric(0)   # no long-memory GARCH-type in this method

  cond_dist <- object@cond_dist

  skew <- if (cond_dist %in% c("snorm", "sstd", "sged", "sald")) {
    model_pars[["skew"]]
  } else {
    numeric(0)
  }

  shape <- if (cond_dist %in% c("std", "sstd", "ged", "sged", "ald", "sald")) {
    sel_word <- switch(
      cond_dist,
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    model_pars[[sel_word]]
  } else  {
    numeric(0)
  }

  mu <- if (exists("mu", where = as.list(model_pars))) {
    if (!is.null(object@nonpar_model)) {
      0
    } else {
      model_pars[["mu"]]
    }
  } else {
    0
  }

  # ARMA / FARIMA specifications
  mean_specs <- fetch_arma_pars_fit(object)
  ar <- mean_specs$ar
  ma <- mean_specs$ma
  D <- mean_specs$D
  pqg0 <- mean_specs$pqg0
  lm_arma <- mean_specs$lm_arma

  dfun <- fun1_selector(cond_dist)

  int_fun <- function(x, gamma_i, delta, shape, skew, dfun) {
    (abs(x) - gamma_i * x)^delta * dfun(x, shape = shape, skew = skew)
  }

  len_gamma <- length(gamma)

  # Store expectations of transformed innovations;
  # one required for each gamma-value
  E_e_transf <- rep(NA, len_gamma)

  for (i in 1:len_gamma) {
    E_e_transf[[i]] <- stats::integrate(
      f = int_fun, lower = -Inf, upper = Inf,
      subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
      gamma_i = gamma[[i]], delta = delta, shape = shape,
      skew = skew, dfun = dfun
    )[[1]]
  }

  # ARMA / FARIMA point forecasts of cond. mean
  mean_fc <- arma_farima_forecasts(ar = ar, ma = ma, D = D, pqg0 = pqg0,
                                   lm_arma = lm_arma, mu = mu,
                                   x = zoo::coredata(object@rt),
                                   innov = zoo::coredata(object@rt - object@cmeans),
                                   n.ahead = n.ahead, trunc = trunc)

  # APARCH point forecasts of cond. volatility

  sigt_fc <- c(sigt_aparch_forecast_short(
      et = zoo::coredata(object@etat),
      sigt = zoo::coredata(object@sigt),
      omega = omega,
      phi = phi,
      beta = beta,
      gamma = gamma,
      delta = delta,
      E_e = E_e_transf,
      horizon = n.ahead
    ))

  list(
    cmeans = mean_fc,
    sigt = sigt_fc
  )

  }
)

# For FIAPARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_fiaparch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

  garch_orders <- object@orders
  p <- garch_orders[[1]]
  q <- garch_orders[[2]]

  n <- length(object@rt)

  model_pars <- object@pars
  names_phi <- paste0("phi", 1:p)
  phi <- model_pars[names_phi]

  names_beta <- paste0("beta", 1:q)
  beta <- -model_pars[names_beta]   # change sign to use esemifar::farima_to_ar function correctly later

  omega <- model_pars[["omega"]]
  gamma <- model_pars[["gamma"]]
  delta <- model_pars[["delta"]]

  d <- model_pars[["d"]]

  mu <- if (exists("mu", where = as.list(model_pars))) {
    if (!is.null(object@nonpar_model)) {
      0
    } else {
      model_pars[["mu"]]
    }
  } else {
    0
  }

  # ARMA / FARIMA specifications
  mean_specs <- fetch_arma_pars_fit(object)
  ar <- mean_specs$ar
  ma <- mean_specs$ma
  D <- mean_specs$D
  pqg0 <- mean_specs$pqg0
  lm_arma <- mean_specs$lm_arma

  cond_dist <- object@cond_dist

  skew <- if (cond_dist %in% c("snorm", "sstd", "sged", "sald")) {
    model_pars[["skew"]]
  } else {
    numeric(0)
  }

  shape <- if (cond_dist %in% c("std", "sstd", "ged", "sged", "ald", "sald")) {
    sel_word <- switch(
      cond_dist,
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    model_pars[[sel_word]]
  } else  {
    numeric(0)
  }

  dfun <- fun1_selector(cond_dist)
  int_fun <- function(x, skew, shape, gamma, delta, dfun) {
    (abs(x) - gamma * x)^delta * dfun(x, shape = shape, skew = skew)
  }
  E_const <- if (n.ahead >= 2) {
    stats::integrate(
      int_fun, lower = -Inf, upper = Inf,
      abs.tol = 1e-6, rel.tol = 1e-6,
      subdivisions = 1000L, skew = skew, shape = shape,
      dfun = dfun, gamma = gamma, delta = delta
    )[[1]]
  } else {
    0
  }


  # ARMA / FARIMA point forecasts of cond. mean
  mean_fc <- arma_farima_forecasts(ar = ar, ma = ma, D = D, pqg0 = pqg0,
                                   lm_arma = lm_arma, mu = mu,
                                   x = zoo::coredata(object@rt),
                                   innov = zoo::coredata(object@rt - object@cmeans),
                                   n.ahead = n.ahead, trunc = trunc)

  # FIAPARCH point forecasts of cond. volatility

  coef_inf_garch <- ar_infty(ar = phi, ma = beta, d = d,
                                     max_i = trunc)[-1]

  sigt_fc <- c(sigt_fiaparch_forecast(
      rt = zoo::coredata(object@rt - object@cmeans),
      coef_inf = coef_inf_garch,
      omega = omega,
      gamma = gamma,
      delta = delta,
      E_const = E_const,
      horizon = n.ahead))

  list(
    cmeans = mean_fc,
    sigt = sigt_fc
  )

  }
)

# For GJR-GARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_gjrgarch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

    object@pars[["delta"]] <- 2

    fun <- methods::selectMethod("fEGarch_predict", "fEGarch_fit_aparch")
    fun(object, n.ahead = n.ahead, trunc = trunc, ...)

  }
)

# For FIGJR-GARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_figjrgarch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

    object@pars[["delta"]] <- 2

    fun <- methods::selectMethod("fEGarch_predict", "fEGarch_fit_fiaparch")
    fun(object, n.ahead = n.ahead, trunc = trunc, ...)

  }
)

# For TGARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_tgarch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

    object@pars[["delta"]] <- 1

    fun <- methods::selectMethod("fEGarch_predict", "fEGarch_fit_aparch")
    fun(object, n.ahead = n.ahead, trunc = trunc, ...)

  }
)

# For FITGARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_fitgarch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

    object@pars[["delta"]] <- 1

    fun <- methods::selectMethod("fEGarch_predict", "fEGarch_fit_fiaparch")
    fun(object, n.ahead = n.ahead, trunc = trunc, ...)

  }
)

# For GARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_garch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

    p <- object@orders[[1]]

    object@pars[["delta"]] <- 2
    l_sub <- rep(0, p)
    names(l_sub) <- paste0("gamma", 1:p)
    object@pars <- c(
      object@pars,
      l_sub
    )

    fun <- methods::selectMethod("fEGarch_predict", "fEGarch_fit_aparch")
    fun(object, n.ahead = n.ahead, trunc = trunc, ...)

  }
)

# For FIGARCH-type
#@rdname forecasting-model-methods
methods::setMethod("fEGarch_predict", "fEGarch_fit_figarch",
  function(object, n.ahead = 10, trunc = NULL, ...) {

    object@pars[["delta"]] <- 2
    object@pars[["gamma"]] <- 0

    fun <- methods::selectMethod("fEGarch_predict", "fEGarch_fit_fiaparch")
    fun(object, n.ahead = n.ahead, trunc = trunc, ...)

  }
)

# Returns a simple list and not the entire input object
#@rdname forecasting-model-methods
methods::setMethod("predict_internal", "fEGarch_fit",
  function(object, n.ahead = 10, trunc = NULL, ...) {

    trunc_obj <- object@trunc
    check <- is.character(trunc_obj) && trunc_obj == "none"

    # Get truncation for long-memory models from fitted object
    if (is.null(trunc) || check) {
      trunc_obj <- object@trunc
      if (check) {
        trunc <- length(object@rt) - 1 + n.ahead
      } else {
        trunc <- trunc_obj
      }
    }

    fcasts <- fEGarch_predict(object = object, n.ahead = n.ahead, trunc = trunc, ...)

    rt <- object@rt
    if (inherits(rt, "zoo")) {
      tp <- stats::time(rt)
      tp_f <- as.Date(utils::tail(tp, 1)) + (1:n.ahead)
      fcasts$cmeans <- zoo::zoo(fcasts$cmeans, order.by = tp_f)
      fcasts$sigt <- zoo::zoo(fcasts$sigt, order.by = tp_f)
    } else if (inherits(rt, "ts")) {
      tp <- time(rt)
      frequ <- stats::frequency(rt)
      start <- utils::tail(tp, 1) + 1 / frequ
      fcasts$cmeans <- stats::ts(fcasts$cmeans, start = start, frequency = frequ)
      fcasts$sigt <- stats::ts(fcasts$sigt, start = start, frequency = frequ)
    }
    fcasts

  }
)

#'Multistep and Rolling Point Forecasts
#'
#'Given a fitted model object from this package, conduct
#'either multistep point forecasts of the conditional
#'means and the conditional standard deviations into
#'the future or rolling point forecasts of arbitrary
#'step size of these quantities for a future test set.
#'
#'@param object an object of class \code{"fEGarch_fit"},
#'i.e. an object returned by either \code{\link{fEGarch}},
#'\code{\link{fiaparch}} or \code{\link{figarch}}, etc.; for
#'\code{predict_roll}, the slot \code{@test_obs} of the
#'fitted model object should not be \code{NULL}.
#'@param n.ahead a single numeric value indicating how
#'far into the future the multistep point forecasts should be
#'produced.
#'@param step_size the step size of the rolling point
#'forecasts; by default, \code{step_size = 1} is employed,
#'i.e. for the immediately subsequent observation time
#'point for the entire test set.
#'@param trunc the truncation setting for the infinite-order
#'polynomial of long-memory model parts; the default uses
#'the setting from the fitted input object \code{object}.
#'@param ... currently without use and included for compatibility
#'with generics.
#'
#'@details
#'Use \code{predict} to compute multistep point forecasts
#'(of the conditional mean and of the conditional standard deviation)
#'into the future. Let \code{n} be the number of observations
#'of the data, to which a model was fitted. Then multistep
#'point forecasts are produced for all future time points
#'from \code{n + 1} to \code{n + n.ahead}.
#'
#'Otherwise, if data was reserved for testing when creating
#'\code{object}, e.g. through the use of the argument
#'\code{n_test} in the corresponding functions, compute
#'rolling point forecasts over the test set using \code{predict_roll}.
#'\code{step_size} then determines the forecasting horizon for
#'the rolling point forecasts. For example, \code{step_size = 1}, i.e.
#'the default, computes one-step rolling point forecasts, whereas for
#'example \code{step_size = 10} computes ten-step rolling point
#'forecasts (starting at the tenth test time point).
#'
#'Refitting of models during the rolling point forecast procedure
#'is currently not yet available.
#'
#'@return
#'Returns an object of class \code{"fEGarch_forecast"} that has the
#'two slots \code{@sigt} and \code{@cmeans} representing the
#'forecasted conditional standard deviations and conditional
#'means, respectively. If the training series saved in \code{object}
#'has a special time series formatting like \code{"zoo"} or \code{"ts"},
#'the formatting is adopted accordingly to these numeric
#'output series. A third slot \code{@model} is the fitted model
#'input object \code{object}.
#'
#'@export
#'
#'@rdname forecasting-methods
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'
#'# Multistep forecasting (EGARCH with cond. normal distr.)
#'model1 <- fEGarch(spec = egarch_spec(), rt)
#'fcast1 <- predict(model1, n.ahead = 10)
#'fcast1
#'
#'# Rolling one-step forecasts (EGARCH with cond. normal distr.)
#'model2 <- fEGarch(spec = egarch_spec(), rt, n_test = 250)
#'fcast2 <- predict_roll(model2, step_size = 1)
#'fcast2
#'
methods::setMethod("predict", "fEGarch_fit",
  function(object, n.ahead = 10, trunc = NULL, ...) {

    np_check <- !is.null(object@nonpar_model)
    object_orig <- object

    if (np_check) {

      mu_h <- utils::tail(zoo::coredata(object@cmeans), 1)
      scale_h <- utils::tail(zoo::coredata(object@scale_fun), 1)

      object@sigt <- object@sigt / object@scale_fun
      object@rt <- (object@rt - mu_h) / object@scale_fun

    }

    fcasts <- predict_internal(object = object, n.ahead = n.ahead, trunc = trunc, ...)

    if (np_check) {
      fcasts$cmeans <- fcasts$cmeans + mu_h
      fcasts$sigt <- fcasts$sigt * scale_h
    }

    fEGarch_forecast(
      cmeans = fcasts$cmeans,
      sigt = fcasts$sigt,
      model = object_orig
    )

  }
)

#'@export
#'
#'@rdname forecasting-methods
#'
methods::setMethod("predict_roll", "fEGarch_fit",
  function(object, step_size = 1, trunc = NULL, ...) {

    # Get truncation for long-memory models from fitted object
    if (is.null(trunc)) {
      trunc_obj <- object@trunc
      if (is.character(trunc_obj) && trunc_obj == "none") {
        trunc <- length(object@rt) - 1 + length(object@test_obs)
      } else {
        trunc <- trunc_obj
      }
    }

    n_test <- length(object@test_obs)
    object2 <- object

    np_check <- !is.null(object@nonpar_model)
    step_size_check <- step_size > 1

    sigt_out <- cmeans_out <- rep(NA, n_test)
    if (np_check) {
      mu_h <- utils::tail(zoo::coredata(object2@cmeans), 1)
      scale_h <- utils::tail(zoo::coredata(object2@scale_fun), 1)
      test_obs <- (object@test_obs - mu_h) / scale_h
      object2@rt <- zoo::coredata((object2@rt - object2@cmeans) / object2@scale_fun)
      object2@sigt <- zoo::coredata(object2@sigt / object2@scale_fun)
      check_mu <- "mu" %in% names(object2@pars)
      object2@cmeans <- if (!check_mu) {
        rep(0, length(object2@rt))
      } else {
        rep(object2@pars[["mu"]], length(object2@rt))
      }
      object2@scale_fun <- NULL
      object2@nonpar_model <- NULL

      if (step_size_check) {
        object3 <- object
        object3@rt <- object2@rt
        object3@sigt <- object2@sigt
        object3@cmeans <- object3@cmeans
        object3@scale_fun <- NULL
        object3@nonpar_model <- NULL
      }

    } else {
      test_obs <- object@test_obs
      object2@rt <- zoo::coredata(object2@rt)
      object2@sigt <- zoo::coredata(object2@sigt)
      object2@cmeans <- zoo::coredata(object2@cmeans)

      if (step_size_check) {
        object3 <- object
        object3@rt <- object2@rt
        object3@sigt <- object2@sigt
        object3@cmeans <- object2@cmeans
      }

    }

    object2@etat <- zoo::coredata(object2@etat)

    for (i in seq_len(n_test)) {
      pred <- predict_internal(object2, n.ahead = 1, trunc = trunc)
      sigt_out[[i]] <- pred[["sigt"]]
      cmeans_out[[i]] <- pred[["cmeans"]]
      object2@rt <- c(object2@rt, test_obs[[i]])
      object2@sigt <- c(object2@sigt, sigt_out[[i]])
      object2@cmeans <- c(object2@cmeans, cmeans_out[[i]])
      object2@etat <- c(object2@etat, (test_obs[[i]] - cmeans_out[[i]]) / sigt_out[[i]])
    }

    if (step_size > 1) {

      object3@etat <- zoo::coredata(object3@etat)
      sigt_out_s <- cmeans_out_s <- rep(NA, n_test)
      for (i in step_size:n_test) {
        pred <- predict_internal(object3, n.ahead = step_size, trunc = trunc)
        sigt_out_s[[i]] <- utils::tail(pred[["sigt"]], 1)
        cmeans_out_s[[i]] <- utils::tail(pred[["cmeans"]], 1)
        imss <- i - step_size + 1
        object3@rt <- c(object3@rt, test_obs[[imss]])
        object3@sigt <- c(object3@sigt, sigt_out[[imss]])
        object3@cmeans <- c(object3@cmeans, cmeans_out[[imss]])
        object3@etat <- c(object3@etat, (test_obs[[imss]] - cmeans_out[[imss]]) / sigt_out[[imss]])
      }

      sigt_out <- sigt_out_s
      cmeans_out <- cmeans_out_s

    }

    rt <- object@rt
    if (inherits(rt, "zoo")) {
      tp_f <- as.Date(time(object@test_obs))
      cmeans_out <- zoo::zoo(cmeans_out, order.by = tp_f)
      sigt_out <- zoo::zoo(sigt_out, order.by = tp_f)
    } else if (inherits(rt, "ts")) {
      tp <- time(object@test_obs)
      frequ <- stats::frequency(object@test_obs)
      start <- tp[[1]]
      cmeans_out <- stats::ts(cmeans_out, start = start, frequency = frequ)
      sigt_out <- stats::ts(sigt_out, start = start, frequency = frequ)
    }

    if (np_check) {
      sigt_out <- sigt_out * scale_h
      cmeans_out <- cmeans_out - cmeans_out + mu_h
    }

    fEGarch_forecast(
      cmeans = cmeans_out,
      sigt = sigt_out,
      model = object
    )

  }
)

