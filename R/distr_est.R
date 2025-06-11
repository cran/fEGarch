#=============================================================#

# Setup for negative log-likelihood

nllhood_setup <- function(theta, y, pdf_fun, extra_par, skew_par, mean_par, sdev_par) {

  idx <- 1
  logic_extra <- is.logical(extra_par)
  if (logic_extra && extra_par) {
    shape <- theta[[idx]]
    idx <- idx + 1
  } else if (logic_extra && !extra_par) {
    shape <- 0
  } else if (!logic_extra && is.numeric(extra_par)) {
    shape <- extra_par
  }

  if (skew_par) {
    skew <- theta[[idx]]
    idx <- idx + 1
  } else {
    skew <- 1
  }

  if (is.null(mean_par)) {
    mean <- theta[[idx]]
    idx <- idx + 1
  } else {
    mean <- mean_par
  }

  if (is.null(sdev_par)) {
    sdev <- theta[[idx]]
  } else {
    sdev <- sdev_par
  }

  # Return negative log-likelihood
  -sum(log(pdf_fun(x = y, shape = shape, skew = skew, mean = mean, sdev = sdev)))

}

#==============================================#

# Estimation functions

### General estimation function

#'MLE for Distribution Fitting
#'
#'Given a vector of values assumed to stem from
#'independent and identically distributed (iid) random variables, fit a selection
#'of distributions, from the normal distribution, the \eqn{t}-distribution, the
#'generalized error distribution (GED), the average Laplace distribution (ALD), and
#'their skewed variants, to the data using maximum-likelihood estimation (MLE).
#'
#'@param x a numeric vector with the data.
#'@param dist a character value that specifies the distribution to consider;
#'available are a normal distribution (\code{"norm"}), a \eqn{t}-distribution
#'(\code{"std"}), a GED (\code{"ged"}), an ALD (\code{"ald"}), and their
#'skewed variants (\code{"snorm"}, \code{"sstd"}, \code{"sged"}, \code{"sald"}).
#'@param fix_mean optional; for the default \code{NULL}, a location parameter representing
#'the (unconditional) mean of the distribution is also being estimated; for any
#'numerical value, however, the mean will be fixed to the corresponding value and
#'therefore excluded from the estimation itself.
#'@param fix_sdev optional; for the default \code{NULL}, a scale parameter representing
#'the (unconditional) standard deviation of the distribution is also being estimated; for any
#'numerical value, however, the standard deviation will be fixed to the corresponding value and
#'therefore excluded from the estimation itself.
#'@param Prange a two-element numeric vector, giving the boundaries of the search space
#'for the shape parameter \eqn{P} in an ALD or its skewed variant.
#'
#'@details
#'
#'Let \eqn{x} be an individual observation. Let \eqn{\mu} a real-valued location
#'parameter, representing the unconditional mean of the distribution, and
#'\eqn{\sigma} a real-valued scale parameter, representing the unconditional
#'standard deviation. Generally, let \eqn{\theta} be a vector with all parameters
#'of the underlying distribution. The likelihood of \eqn{x} is given through
#'\deqn{L_x^{\text{norm}}(\theta)=\frac{\sigma^{-1}}{\sqrt{2\pi}}\exp\left(-\frac{1}{2}\left(\frac{x-\mu}{\sigma}\right)^2\right)}
#'for a normal distribution,
#'\deqn{L_x^{\text{std}}(\theta)=\frac{\sigma^{-1}\Gamma\left(\frac{\nu + 1}{2}\right)}{\Gamma\left(\frac{\nu}{2}\right)\sqrt{\pi(\nu-2)}}\left[1+\frac{1}{\nu - 2}\left(\frac{x-\mu}{\sigma}\right)^2\right]^{-\frac{\nu + 1}{2}}}
#'for a \eqn{t}-distribution with \eqn{\nu} as the degrees of freedom and \eqn{\Gamma}
#'as the gamma function,
#'\deqn{L_x^{\text{ged}}(\theta)=\frac{\sigma^{-1}\beta}{2}\sqrt{\frac{C_{\Gamma,3}}{C_{\Gamma,1}^3}} \exp\left\{-\left|\frac{x-\mu}{\sigma}\right|^{\beta}\left(\frac{C_{\Gamma,3}}{C_{\Gamma,1}}\right)^{\frac{\beta}{2}}\right\}}
#'for a GED with \eqn{\beta} as its real-valued shape and with \eqn{C_{\Gamma,i}=\Gamma\left(\frac{i}{\beta}\right)}, \eqn{i\in\left\{1,3\right\}},
#'and in
#'\deqn{L_x^{\text{ald}}(\theta)=\frac{\sigma^{-1}sB}{2}\exp\left(-s\left|\frac{x-\mu}{\sigma}\right|\right)\sum_{j=0}^{P}c_j\left(s\left|\frac{x-\mu}{\sigma}\right|\right)^j}
#'for an ALD with \eqn{P} as its discrete shape, where \eqn{s = \sqrt{2(P+1)}},
#'\deqn{B=2^{-2P} {{2P}\choose{P}}, \hspace{4mm} P \geq 0,}
#'and
#'\deqn{c_{j}=\frac{2(P-j + 1)}{j(2P-j+1)}c_{j-1}, \hspace{4mm} j =2,3,\dots,P,}
#'with \eqn{c_0 = c_1 = 1}.
#'The individual-observation likelihoods for the skewed variants are derived analogously
#'from the idea by Fernandez and Steel (1998). The log-likelihoods to maximize
#'over are then just the sum of the log-transformed likelihoods for each observation.
#'
#'\code{distr_est} is a general purpose distribution fitting function, where the
#'distribution can be selected through the argument \code{dist}. \code{norm_est},
#'\code{std_est}, \code{ged_est}, \code{ald_est}, \code{snorm_est},
#'\code{sstd_est}, \code{sged_est}, and \code{sald_est} are wrappers around
#'\code{distr_est} in order to directly provide fitting functions for the
#'different distributions available in this package.
#'
#'@return
#'Returns a list with the following elements.
#'
#'@references
#'\itemize{
#'\item{Fernandez, C., & Steel, M. F. J. (1998). Bayesian Modeling of Fat Tails and Skewness.
#'Journal of the American Statistical Association,
#'93(441), 359–371. DOI: 10.1080/01621459.1998.10474117.}
#'}
#'@export
#'
#'@rdname distribution_estimation
#'
#'@examples
#'# Draw obs. from GED and implement standard deviation 1.2
#'# and mean 3.1
#'x <- rged_s(4000, shape = 1.5) * 1.2 + 3.1
#'# Fit GED
#'ged_est(x)
#'# Fit GED differently using distr_est()
#'distr_est(x, dist = "ged")
#'# Fit GED while fixing mean and standard deviation
#'ged_est(x, fix_mean = 3.1, fix_sdev = 1.2)
#'# Fit another distribution
#'sstd_est(x)
#'

distr_est <- function(x, dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                      fix_mean = NULL, fix_sdev = NULL, Prange = c(1, 5)) {

  x_orig <- x
  x <- zoo::coredata(x)

  dist <- match.arg(dist)

  dfun_s <- fun1_selector(dist)
  dfun <- function(x, shape, skew, mean, sdev) {
    dfun_s(x = (x - mean) / sdev, shape = shape, skew = skew) / sdev
  }

  extra_par <- switch(
    dist,
    "norm" = FALSE,
    "std" = TRUE,
    "ged" = TRUE,
    "ald" = TRUE,
    "snorm" = FALSE,
    "sstd" = TRUE,
    "sged" = TRUE,
    "sald" = TRUE
  )

  skew_par <- switch(
    dist,
    "norm" = FALSE,
    "std" = FALSE,
    "ged" = FALSE,
    "ald" = FALSE,
    "snorm" = TRUE,
    "sstd" = TRUE,
    "sged" = TRUE,
    "sald" = TRUE
  )

  # Fix the negative log-likelihood
  nllhood_s <- function(theta, y, extra_par) {
    nllhood_setup(
      theta = theta,
      y = y,
      pdf_fun = dfun,
      extra_par = extra_par,
      skew_par = skew_par,
      mean_par = fix_mean,
      sdev_par = fix_sdev
    )
  }

  if (is.logical(extra_par) && extra_par) {
    start_extra_par <- switch(
      dist,
      "std" = 10,
      "ged" = 2,
      "sstd" = 10,
      "sged" = 2
    )
    LB_extra_par <- switch(
      dist,
      "std" = 2,
      "ged" = 1,
      "sstd" = 2,
      "sged" = 1
    )
    UB_extra_par <- switch(
      dist,
      "std" = Inf,
      "ged" = 5,
      "sstd" = Inf,
      "sged" = 5
    )
    name_extra_par <- switch(
      dist,
      "norm" = character(0),
      "std" = "df",
      "ged" = "shape",
      "ald" = "P",
      "snorm" = character(0),
      "sstd" = "df",
      "sged" = "shape",
      "sald" = "P"
    )
  } else {
    start_extra_par <- LB_extra_par <- UB_extra_par <- numeric(0)
    name_extra_par <- character(0)
  }

  if (skew_par) {
    start_skew <- 1
    LB_skew <- 0.3
    UB_skew <- 3
    name_skew <- "skew"
  } else {
    start_skew <- LB_skew <- UB_skew <- numeric(0)
    name_skew <- character(0)
  }

  if (is.null(fix_mean)) {
    start_mean <- mean(x)
    LB_mean <- min(x)
    UB_mean <- max(x)
    name_mean <- "mean"
  } else {
    start_mean <- LB_mean <- UB_mean <- numeric(0)
    name_mean <- character(0)
  }

  if (is.null(fix_sdev)) {
    SD_x <- stats::sd(x)
    start_sdev <- SD_x
    LB_sdev <- SD_x / 4
    UB_sdev <- SD_x * 4
    name_sdev <- "sdev"
  } else {
    start_sdev <- LB_sdev <- UB_sdev <- numeric(0)
    name_sdev <- character(0)
  }

  start <- c(
    start_extra_par,
    start_skew,
    start_mean,
    start_sdev
  )
  LB <- c(
    LB_extra_par,
    LB_skew,
    LB_mean,
    LB_sdev
  )
  UB <- c(
    UB_extra_par,
    UB_skew,
    UB_mean,
    UB_sdev
  )

  check_lstart <- length(start) > 0
  check_ald <- dist %in% c("ald", "sald")

  if (check_ald) {

    P_all <- seq(Prange[[1]], Prange[[2]], 1)
    i <- 1
    lP <- length(P_all)
    ests <- vector(mode = "list", length = lP)
    nllh <- rep(NA, lP)

    for (P0 in P_all) {

      nllhood <- function(theta, y) {
        nllhood_s(
          theta = theta,
          y = y,
          extra_par = P0
        )
      }

      ests[[i]] <- if (check_lstart) {
        suppressWarnings(Rsolnp::solnp(
          pars = start, fun = nllhood, LB = LB, UB = UB, y = x,
          control = list(trace = FALSE)
        ))
      } else {
        list(
          pars = numeric(0),
          values = nllhood(theta = numeric(0), y = x)
        )
      }

      nllh[[i]] <- utils::tail(ests[[i]]$values, 1)

      i <- i + 1
    }

    ind <- which.min(nllh)

    est <- ests[[ind]]
    P <- P_all[[ind]]

  } else {

    nllhood <- function(theta, y) {
      nllhood_s(
        theta = theta,
        y = y,
        extra_par = extra_par
      )
    }

    est <- if (check_lstart) {

      suppressWarnings(Rsolnp::solnp(
        pars = start, fun = nllhood, LB = LB, UB = UB, y = x,
        control = list(trace = FALSE)
      ))

    } else {

      list(
        pars = numeric(0),
        values = nllhood(theta = numeric(0), y = x)
      )

    }

  }

  if (!is.null(est$convergence) && est$convergence != 0){
    warning("Convergence failed.")
  }

  par <- est$pars
  par_names <- c(
    name_extra_par,
    name_skew,
    name_mean,
    name_sdev
  )

  if (check_lstart) {
    if (check_ald) {
      nllhood <- function(theta, y) {
       nllhood_s(
         theta = theta,
         y = y,
         extra_par = P
       )
      }
    }
    hess <- numDeriv::hessian(nllhood, x = par, y = x)
    vcov <- solve(hess)
    serror <- sqrt(diag(vcov))
    if (any(is.nan(serror))) {
      warning("Unable to compute standard errors.")
    }
    if (check_ald) {
      par <- c(P, par)
      serror <- c(NA_real_, serror)
      vcov <- rbind(NA_real_, cbind(NA_real_, vcov))
    }
    names(par) <- par_names
    names(serror) <- par_names
    rownames(vcov) <- par_names
    colnames(vcov) <- par_names
  } else {
    if (check_ald) {
      par <- P
      serror <- NA_real_
      vcov <- matrix(NA_real_, ncol = 1, nrow = 1)
      names(par) <- "P"
      names(serror) <- "P"
      rownames(vcov) <- "P"
      colnames(vcov) <- "P"
    } else {
      vcov <- serror <- numeric(0)
    }

  }

  k <- length(par)
  llhood <- -utils::tail(est$values, 1)
  n <- length(x)
  m2llhood <- -2 * llhood

  aic <- (m2llhood + 2 * k) / n
  bic <- (m2llhood + k * log(n)) / n

  if (!is.null(fix_mean) || !is.null(fix_sdev)) {
    fm <- if (!is.null(fix_mean)) {
      name_m <- "mean"
      fix_mean
    } else {
      name_m <- character(0)
      numeric(0)
    }

    sm <- if (!is.null(fix_sdev)) {
      name_s <- "sdev"
      fix_sdev
    } else {
      name_s <- character(0)
      numeric(0)
    }

    fixed <- c(fm, sm)
    names(fixed) <- c(name_m, name_s)

  } else {
    fixed <- numeric(0)
  }

  fEGarch_distr_est(
    pars = par,
    se = serror,
    vcov_mat = vcov,
    x = x_orig,
    inf_criteria = c("aic" = aic, "bic" = bic),
    dist = dist,
    fixed = fixed,
    llhood = llhood
  )

}

# Special normal distribution estimation wrapper

#'@export
#'@rdname distribution_estimation
norm_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "norm", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special t-distribution estimation wrapper

#'@export
#'@rdname distribution_estimation
std_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "std", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special GED estimation wrapper

#'@export
#'@rdname distribution_estimation
ged_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "ged", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special ALD estimation wrapper

#'@export
#'@rdname distribution_estimation
ald_est <- function(x, fix_mean = NULL, fix_sdev = NULL, Prange = c(1, 5)) {

  distr_est(x = x, dist = "ald", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = Prange)

}

# Special skewed normal distribution estimation wrapper

#'@export
#'@rdname distribution_estimation
snorm_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "snorm", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special skewed t-distribution estimation wrapper

#'@export
#'@rdname distribution_estimation
sstd_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "sstd", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special skewed GED estimation wrapper

#'@export
#'@rdname distribution_estimation
sged_est <- function(x, fix_mean = NULL, fix_sdev = NULL) {

  distr_est(x = x, dist = "sged", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = c(1, 5))

}

# Special skewed ALD estimation wrapper

#'@export
#'@rdname distribution_estimation
sald_est <- function(x, fix_mean = NULL, fix_sdev = NULL, Prange = c(1, 5)) {

  distr_est(x = x, dist = "sald", fix_mean = fix_mean, fix_sdev = fix_sdev,
            Prange = Prange)

}

#'Optimal Distribution Fitting to IID Data
#'
#'Given a series supposed to be from independent and identically distributed (iid)
#'random variables, fit all eight distributions of this package to the data
#'using maximum-likelihood estimation (MLE) and select the best one following
#'either the BIC (the default) or the AIC.
#'
#'@param x the vector of iid values to fit distributions to.
#'@param dists a vector with all the distribution abbreviations which
#'should be considered in the selection process; by default, all eight
#'distributions of this package are considered.
#'@param fix_mean a value to fix the unconditional mean of the distribution to;
#'with the default \code{NULL}, the unconditional mean is estimated as an extra parameter.
#'@param fix_sdev a value to fix the unconditional standard deviation of the distribution to;
#'with the default \code{NULL}, the unconditional standard deviation is estimated as an extra parameter.
#'@param Prange a two-element vector giving the search range for the shape parameter
#'\eqn{P} of the (skewed) average Laplace distribution.
#'@param criterion either \code{"bic"} or \code{"aic"} to use BIC or AIC as
#'the final selection criterion; by default \code{"bic"} is implemented.
#'
#'@details
#'For information on the method and distributions, we refer the reader to
#'\code{\link{distr_est}}.
#'
#'@return
#'Returns an object of class \code{"fEGarch_distr_est"} with various slots
#'representing the estimation results of the selected fitted distribution.
#'
#'@export
#'
#'@examples
#'x <- rnorm(2000) * 2.1 + 10.5
#'find_dist(x)
#'
find_dist <- function(x, dists = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"), fix_mean = NULL, fix_sdev = NULL, Prange = c(1, 5), criterion = c("bic", "aic")) {

  criterion <- match.arg(criterion)

  ests <- lapply(
    X = dists,
    FUN = function(.y, x, fix_mean, fix_sdev, Prange) {
      suppressWarnings(distr_est(x = x, dist = .y, fix_mean = fix_mean, fix_sdev = fix_sdev, Prange = Prange))
    },
    x = x, fix_mean = fix_mean, fix_sdev = fix_sdev, Prange = Prange
  )

  bics <- vapply(
    X = ests,
    FUN = function(.x, criterion) {
      inf_criteria(.x)[[criterion]]
    },
    FUN.VALUE = numeric(1), criterion = criterion
  )

  idx <- which.min(bics)

  ests[[idx]]

}


#'@rdname measure_risk
#'@export

setMethod("measure_risk", "fEGarch_distr_est",
  function(object, measure = c("VaR", "ES"), level = c(0.975, 0.99), test_obs, sigt, cmeans, ...) {

  # check inheritance of time series formatting
  inp_list <- list(
    test_obs, sigt, cmeans
  )
  checks_zoo <- vapply(
    X = inp_list,
    FUN = function(.x) {
      inherits(.x, what = "zoo")
    },
    FUN.VALUE = logical(1)
  )
  check_zoo <- any(checks_zoo)

  checks_ts <- vapply(
    X = inp_list,
    FUN = function(.x) {
      inherits(.x, what = "ts")
    },
    FUN.VALUE = logical(1)
  )
  check_ts <- any(checks_ts)

  format_giver <- if (check_zoo) {
    inp_list[checks_zoo][[1]]
  } else if (check_ts) {
    inp_list[checks_ts][[1]]
  }

  test_obs <- zoo::coredata(test_obs)
  sigt <- zoo::coredata(sigt)
  cmeans <- zoo::coredata(cmeans)

  if (check_zoo || check_ts) {
    series_out <- format_applier_ts(
      rt = format_giver,
      list_of_ts = list(
        "test_obs" = test_obs,
        "sigt" = sigt,
        "cmeans" = cmeans
      )
    )

  test_obs <- series_out$test_obs
  sigt <- series_out$sigt
  cmeans <- series_out$cmeans

  }

  dist_name <- object@dist
  par <- pars(object)

  args <- list("dist" = dist_name)

  if (dist_name %in% c("std", "ged", "ald", "sstd", "sged", "sald")) {
    par_name <- switch(
      dist_name,
      "std" = "df",
      "ged" = "shape",
      "ald" = "P",
      "sstd" = "df",
      "sged" = "shape",
      "sald" = "P"
    )
    par_val <- par[[par_name]]
    args[[par_name]] <- par_val
  }

  if(dist_name %in% c("snorm", "sstd", "sged", "sald")) {
    args[["skew"]] <- par[["skew"]]
  }

  VaR_out <- ES_out <- vector(mode = "list", length = length(level))



  if ("VaR" %in% measure) {

    i <- 0

    for (lev in level) {
      i <- i + 1
      args[["level"]] <- lev
      VaR <- do.call(what = VaR_calc, args = args)

      VaR_out[[i]] <- cmeans + sigt * VaR

    }

    names(VaR_out) <- paste0("VaR", level)

  } else {
    VaR_out <- list()
  }

  if ("ES" %in% measure) {

    i <- 0

    for (lev in level) {
      i <- i + 1
      args[["level"]] <- lev
      ES <- do.call(what = ES_calc, args = args)

      ES_out[[i]] <- cmeans + sigt * ES

    }

    names(ES_out) <- paste0("ES", level)

  } else {
    ES_out <- list()
  }

  fEGarch_risk(
    measures = list("VaR" = VaR_out, "ES" = ES_out),
    observations = test_obs,
    model = object,
    sigt = sigt,
    cmeans = cmeans
  )

})
