
# For short-memory: coef_inf and presample as arguments irrelevant
arma_fit <- function(x, mu, ar, ma, coef_inf, presample) {
  c(fitted_arma_Cpp(x = x, mu = mu, ma = ma, ar = ar))
}

# Using convolution techniques via FFT
farima_fit <- function(x, mu, ar, ma, coef_inf, presample) {
  x_dm <- x - mu
  pres_val <- if (mu == 0) {     # If mean is fixed at 0, use mean as starting values
    mean(x)
  } else {
    0
  }
  x_upd <- c(rep(pres_val, presample), x_dm)
  n0 <- length(x_dm)
  n <- length(x_upd)
  np2 <- nextpow2(2 * n - 1)
  cl <- length(coef_inf)
  c_adj <- c(coef_inf, rep(0, np2 - cl))
  x_adj <- c(x_upd, rep(0, np2 - n))
  mu + Re(utils::tail(utils::head(stats::fft(stats::fft(c_adj) * stats::fft(x_adj), inverse = TRUE), n), n0)) / np2
}

# Fit a preliminary FARIMA using (conditional) normal QMLE
# (used to determine better starting values for conditional volatility
#  equations in dual models with long memory in the mean; EGF not affected
#  because there the expectation of g(.) is zero)
prelim_farima <- function(x, p, q, include_mean, presample, trunc, Drange) {
  stopifnot("x must be a numeric vector without NAs." = (is.numeric(x) && all(!is.na(x))))
  stopifnot("FARIMA order p must be numeric and of length one." = (is.numeric(p) && length(p) == 1))
  stopifnot("FARIMA order q must be numeric and of length one." = (is.numeric(q) && length(q) == 1))
  stopifnot("include_mean must be logical of length one." = (is.logical(include_mean) && length(include_mean) == 1))
  stopifnot("presample must be a single numerical value." = (is.numeric(presample) && length(presample) == 1))
  stopifnot("trunc must be either NULL or a single numerical value." = (is.null(trunc) || (is.numeric(trunc) && length(trunc) == 1)))
  stopifnot("Drange must be a length two numeric vector." = (is.numeric(Drange) && length(Drange) == 2 && Drange[[1]] < Drange[[2]]))

  if (is.null(trunc)) {trunc <- length(x) + presample - 1}
  mu <- if (include_mean) {mean(x)} else {0}
  n <- length(x)

  obj_fun <- function(theta, x) {

    # Grab values from parameter vector
    idx <- 1
    sigma <- theta[[idx]]
    idx <- idx + 1
    ar <- if (p > 0) {
      theta[idx:(idx + p - 1)]
    } else {
      numeric(0)
    }
    idx <- idx + p
    ma <- if (q > 0) {
      theta[idx:(idx + q - 1)]
    } else {
      numeric(0)
    }
    idx <- idx + q
    D <- theta[[idx]]

    # Compute (truncated) infinite-length coef. series
    coef_inf <- ar_infty(ar = ar, ma = ma, d = D, max_i = trunc)
    coef_inf[[1]] <- 0
    # Compute residuals (obs - fitted cond. means)
    et <- x - farima_fit(x = x, mu = mu, ar = ar, ma = ma, coef_inf = coef_inf, presample = presample)

    # Compute neg. log-likelihood (under normality assumption)
    2 * n * log(sigma) + sum((et / sigma)^2)


  }


  start_pars <- c(
    0.9 * stats::sd(x),
    rep(0.05 / p, p),
    rep(0.05 / q, q),
    0.25
  )

  low_pars <- c(
    1e-6,
    rep(-1.5, p),
    rep(-1.5, q),
    Drange[[1]]
  )

  up_pars <- c(
    stats::sd(x),
    rep(1.5, p),
    rep(1.5, q),
    Drange[[2]]
  )

  fit <- suppressWarnings(stats::nlminb(
    start = start_pars,
    objective = obj_fun,
    lower = low_pars,
    upper = up_pars,
    x = x
  ))

  fit_pars <- c(mu, fit$par)
  names_pars <- c(
    "mu",
    "sigma",
    if (p > 0) {paste0("ar", 1:p)} else {character(0)},
    if (q > 0) {paste0("ma", 1:q)} else {character(0)},
    "D"
  )

  names(fit_pars) <- names_pars


  farima_resid_comp <- function(theta, x) {
    idx <- 1
    sigma <- theta[[idx]]
    idx <- idx + 1
    ar <- if (p > 0) {
      theta[idx:(idx + p - 1)]
    } else {
      numeric(0)
    }
    idx <- idx + p
    ma <- if (q > 0) {
      theta[idx:(idx + q - 1)]
    } else {
      numeric(0)
    }
    idx <- idx + q
    D <- theta[[idx]]

    # Compute (truncated) infinite-length coef. series
    coef_inf <- ar_infty(ar = ar, ma = ma, d = D, max_i = trunc)
    coef_inf[[1]] <- 0
    # Compute residuals (obs - fitted cond. means)
    et <- x - farima_fit(x = x, mu = mu, ar = ar, ma = ma, coef_inf = coef_inf, presample = presample)
    et
  }

  et <- farima_resid_comp(theta = unname(fit_pars)[-1], x = x)

  list(pars = fit_pars, residuals = et)

}
