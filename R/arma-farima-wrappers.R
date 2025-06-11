
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

