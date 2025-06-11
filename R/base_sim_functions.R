cdf_fun <- function(x, dfun, ...) {
  vapply(
    x,
    FUN = function(x, dfun, ...) {
      integrate(f = dfun, lower = -Inf, upper = x,
                ..., subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6)[[1]]
    },
    FUN.VALUE = numeric(1), dfun = dfun, ...
  )
}

inv_cdf_fun_helper <- function(x, p, dfun, ...) {
  cdf_fun(x, dfun = dfun, ...) - p
}

inv_cdf_fun <- function(p, dfun, ...) {
  vapply(
    p,
    FUN = function(p, dfun, ...) {
       uniroot(inv_cdf_fun_helper, interval = c(-5, 5), dfun = dfun, p = p,
               ..., extendInt = "yes")[[1]]
    },
    FUN.VALUE = numeric(1), dfun = dfun, ...
  )

}

#'Sampling Functions for Innovations
#'
#'Draw random samples from a normal distribution, \eqn{t}-distribution,
#'generalized error distribution, or their skewed variants (all standardized
#'to have mean zero and variance one).
#'
#'@param n the number of observations to draw.
#'@param df the degrees of freedom for a (skewed) \eqn{t}-distribution.
#'@param shape the shape parameter for a (skewed) generalized error
#'distribution.
#'@param P the number of Laplace distributions (minus 1) to derive the arithmetic
#'mean from as the basis for a (skewed) average Laplace (AL) distribution
#'distribution.
#'@param skew the skewness parameter in the skewed distributions.
#'
#'@details
#'Draw random samples from a normal distribution, \eqn{t}-distribution,
#'generalized error distribution, an average Laplace distribution,
#'or their skewed variants (all standardized to have mean zero and
#'variance one).
#'
#'@return
#'These functions return a numeric vector of length \code{n}.
#'
#'@name sim_functions
#'
#'@export
#'
#'@examples
#'rnorm_s(10)
#'
#'rsstd_s(10, df = 7, skew = 0.9)

rnorm_s <- function(n) {
  rnorm(n = n, mean = 0, sd = 1)
}

#'@export
#'@rdname sim_functions

rstd_s <- function(n, df = 10000) {
  sdev <- sqrt(df / (df - 2))
  rt(n = n, df = df) / sdev
}

qged_s <- function(p, shape) {
  a <- sqrt(gamma(1 / shape) / gamma(3 / shape))
  intermed <- stats::qgamma(2 * abs(p - 0.5), shape = 1 / shape, scale = 1)
  sign(p - 0.5) * (a^shape * intermed)^(1 / shape)

}

#'@export
#'@rdname sim_functions
rged_s <- function(n, shape = 2) {
  p <- runif(n)
  qged_s(p = p, shape = shape)
}

#'@export
#'@rdname sim_functions
rsged_s <- function(n, shape = 2, skew = 1) {
	weight <- skew / (skew + 1 / skew)
	z <- -runif(n, min = -weight, max = 1.0 - weight)
	fac <- ifelse(z > 0, yes = skew, no = 1 / skew)
	intmed <- fac * abs(rged_s(n, shape = shape)) * sign(z)
	M1 <- gamma(2 / shape) / sqrt(gamma(1 / shape) * gamma(3 / shape))
	mu <- M1 * (skew - 1 / skew)
	sigma <- sqrt((1 - M1^2) * (skew^2 + 1 / skew^2) + 2 * M1^2 - 1)
	out <- (intmed - mu) / sigma
	out
}

qlap <- function(p, scale) {
  idx <- p <= 0.5
  ifelse(
    idx,
    yes = scale * log(2 * p[idx]),
    no = -scale * log(2 - 2 * p[!idx])
  )
}

#'@export
#'@rdname sim_functions
rald_s <- function(n, P = 8) {
  nl <- P + 1  # number of Laplace distributed RVs
  m <- nl * n  # number of Laplace values to simulate

  p <- runif(m)
  a <- sqrt(nl / 2)
  mat <- matrix(
    qlap(p, scale = a),
    nrow = n, ncol = nl
  )

  apply(mat, MARGIN = 1, FUN = mean, simplify = TRUE)

}

#'@export
#'@rdname sim_functions
rsnorm_s <- function(n, skew = 1) {
	weight <- skew / (skew + 1 / skew)
	z <- -runif(n, min = -weight, max = 1.0 - weight)
	fac <- ifelse(z > 0, yes = skew, no = 1 / skew)
	intmed <- fac * abs(rnorm_s(n)) * sign(z)
	M1 <- sqrt(2 / pi)
	mu <- M1 * (skew - 1 / skew)
	sigma <- sqrt((1 - M1^2) * (skew^2 + 1 / skew^2) + 2 * M1^2 - 1)
	out <- (intmed - mu) / sigma
	out
}

#'@export
#'@rdname sim_functions
rsstd_s <- function(n, df = 10000, skew = 1) {
	weight <- skew / (skew + 1 / skew)
	z <- -runif(n, min = -weight, max = 1.0 - weight)
	fac <- ifelse(z > 0, yes = skew, no = 1 / skew)
	intmed <- fac * abs(rstd_s(n, df = df)) * sign(z)
	M1 <- 2 * sqrt(df - 2) * gamma((df + 1) / 2) / ((df - 1) * sqrt(pi) * gamma(df / 2))
	mu <- M1 * (skew - 1 / skew)
	sigma <- sqrt((1 - M1^2) * (skew^2 + 1 / skew^2) + 2 * M1^2 - 1)
	out <- (intmed - mu) / sigma
	out
}


#' rsged_s <- function(n, shape = 2, skew = 1) {
#'   p <- runif(n)
#'   inv_cdf_fun(p = p, dfun = pdf_skew_sged_s, shape = shape, skew = skew)
#' }

#'@export
#'@rdname sim_functions
rsald_s <- function(n, P = 8, skew = 1) {
	weight <- skew / (skew + 1 / skew)
	z <- -runif(n, min = -weight, max = 1.0 - weight)
	fac <- ifelse(z > 0, yes = skew, no = 1 / skew)
	intmed <- fac * abs(rald_s(n, P = P)) * sign(z)
	M1 <- ald_first_abs_mom(P = P)
	mu <- M1 * (skew - 1 / skew)
	sigma <- sqrt((1 - M1^2) * (skew^2 + 1 / skew^2) + 2 * M1^2 - 1)
	out <- (intmed - mu) / sigma
	out
}
