#'General GARCH-Type Model Estimation
#'
#'Fit any of the additional short- or long-memory GARCH-type
#'models from the \code{fEGarch} package aside from those
#'of the extended EGARCH family.
#'
#'@param rt the observed series ordered from past to present; can be
#'a numeric vector, a \code{"zoo"} class time series object, or a
#'\code{"ts"} class time series object.
#'@param model any character object among \code{"garch"}, \code{"gjrgarch"},
#'\code{"aparch"}, \code{"tgarch"}, \code{"figarch"}, \code{"figjrgarch"},
#'\code{"fitgarch"}and \code{"fiaparch"}.
#'@param orders a two-element numeric vector containing the two model
#'orders \eqn{p} and \eqn{q} (see Details for more information); currently,
#'only the default \code{orders = c(1, 1)} is supported for long-memory
#'models; other specifications
#'of a two-element numeric vector will lead to \code{orders = c(1, 1)} being
#'run and a warning message being returned for long-memory models.
#'@param cond_dist the conditional distribution to consider as a
#'character object; the default is a conditional normal distribution
#'\code{"norm"}; available are also, however, a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution (\code{"ged"}), an
#'average Laplace distribution (\code{"ald"}),
#'and their four skewed variants (\code{"snorm"}, \code{"sstd"},
#'\code{"sged"}, \code{"sald"}).
#'@param start_pars the starting parameters for the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param drange a two-element numeric vector that gives the boundaries of the
#'search interval for the fractional differencing parameter \eqn{d} in the
#'conditional volatility model part of a long-memory model; is
#'overwritten by the settings of the arguments \code{LB} and \code{UB}.
#'@param meanspec an object of class "mean_spec"; indicates the
#'specifications for the model in the conditional mean.
#'@param Drange a two-element numeric vector that indicates the boundaries
#'of the interval over which to search for the fractional differencing
#'parameter \eqn{D} in a long-memory ARMA-type model in the conditional
#'mean model part; by default,
#'\eqn{D} being searched for on the
#'interval from 0 to \eqn{0.5 - 1\times 10^{-6}}; note that specific
#'settings in the arguments
#'\code{LB} and \code{UB} overwrite this argument.
#'@param nonparspec an object of class \code{"locpol_spec"} returned
#'by \code{\link{locpol_spec}}; defines the settings of the nonparametric
#'smoothing technique for \code{use_nonpar = TRUE}.
#'@param use_nonpar a logical indicating whether or not to implement a
#'semiparametric extension of the volatility model defined through \code{spec};
#'see "Details" for more information.
#'@param n_test a single numerical value indicating, how many observations
#'at the end of \code{rt} not to include in the fitting process and to
#'reserve for backtesting.
#'@param LB the lower boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param UB the upper boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param control a list that is passed to \code{control} of the
#'function \code{solnp} of the package \code{Rsolnp}.
#'@param control_nonpar a list containing changes to the arguments
#'for the hyperparameter estimation algorithm in the nonparametric
#'scale function estimation for
#'\code{use_nonpar = TRUE}; see "Details" for more information.
#'@param mean_after_nonpar only for \code{use_nonpar = TRUE}; considers the unconditional mean
#'of the parametric model part in the QMLE step in a semiparametric model; by default, a zero-mean
#'model is considered for the parametric part in a semiparametric model.
#'@param parallel only relevant for a (skewed) average Laplace (AL)
#'distribution, i.e.
#'if \code{cond_dist} in \code{spec} is set to \code{cond_dist = "ald"} or
#'\code{cond_dist = "sald"}; \code{parallel} is a logical value indicating whether
#'or not the slices for the positive integer-valued parameter of the SM
#'distribution should be fitted in parallel for a speed boost.
#'@param ncores only relevant for a (skewed) average Laplace (AL)
#'distribution, i.e.
#'if \code{cond_dist} in \code{spec} is set to \code{cond_dist = "ald"} or
#'\code{cond_dist = "sald"}, and if simultaneously \code{parallel = TRUE};
#'\code{ncores} is a single numeric value indicating the number of cores to
#'use for parallel computations.
#'@param trunc a positive integer indicating the finite truncation length of the
#'infinite-order polynomials of the infinite-order representations of the
#'long-memory model parts; the character \code{"none"} is an optional input
#'that specifies that truncation should always be applied back to the first (presample) observation
#'time point, i.e. that maximum length filters should be applied at all times.
#'@param presample the presample length for initialization (for extended EGARCH- / Log-GARCH-type
#'models only relevant for the FARIMA-part, as series in log-transformed
#'conditional variance are initialized by zero).
#'@param Prange a two-element vector that indicates the search boundaries for
#'the parameter \eqn{P} in a (skewed) average Laplace distribution.
#'
#'@details
#'See the documentation on \code{\link{garch}}, \code{\link{gjrgarch}},
#'\code{\link{tgarch}}, \code{\link{aparch}}, \code{\link{figarch}},
#'\code{\link{figjrgarch}}, \code{\link{fitgarch}} and \code{\link{fiaparch}} for more detailed
#'information on the corresponding models and functions selectable through this
#'wrapper function.
#'
#'@export
#'
#'@return
#'An object of S4-class \code{"fEGarch_fit_garch"}, \code{"fEGarch_fit_gjrgarch"},
#'\code{"fEGarch_fit_tgarch"},
#'\code{"fEGarch_fit_aparch"}, \code{"fEGarch_fit_figarch"}, \code{"fEGarch_fit_figjrgarch"},
#'\code{"fEGarch_fit_fitgarch"} or \code{"fEGarch_fit_fiaparch"}
#'is returned depending on the selected input for the argument
#'\code{model}. The object then contains the following elements.
#'\describe{
#'\item{\code{pars}:}{a named numeric vector with the parameter estimates.}
#'\item{\code{se}:}{a named numeric vector with the obtained standard errors in accordance with the parameter estimates.}
#'\item{\code{vcov_mat}:}{the variance-covariance matrix of the parameter estimates with named columns and rows.}
#'\item{\code{rt}:}{the input object \code{rt} (or at least the training data, if \code{n_test} is greater than zero);
#'if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is kept.}
#'\item{\code{cmeans}:}{the estimated conditional means; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{cmeans}.}
#'\item{\code{sigt}:}{the estimated conditional standard deviations (or for \code{use_nonpar = TRUE} the estimated total
#'volatilities, i.e. scale function value times conditional standard deviation); if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{sigt}.}
#'\item{\code{etat}:}{the obtained residuals; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{etat}.}
#'\item{\code{orders}:}{a two-element numeric vector stating the considered model orders.}
#'\item{\code{cond_dist}:}{a character value stating the conditional distribution considered in the model fitting.}
#'\item{\code{long_memo}:}{a logical value stating whether or not long memory was considered in the model fitting.}
#'\item{\code{llhood}:}{the log-likelihood value obtained at the optimal parameter combination.}
#'\item{\code{inf_criteria}:}{a named two-element numeric vector with the corresponding AIC (first element) and BIC (second element) of the fitted parametric model part; for purely parametric models, these criteria are valid for the entire model; for semiparametric models, they are only valid for the parametric step and are not valid for the entire model.}
#'\item{\code{meanspec}:}{the settings for the model in the conditional mean; is an object
#'of class \code{"mean_spec"} that is identical to the object passed to the input argument
#'\code{meanspec}.}
#'\item{\code{test_obs}:}{the observations at the end up the input \code{rt} reserved for
#'testing following \code{n_test}.}
#'\item{\code{scale_fun}:}{the estimated scale function values, if \code{use_nonpar = TRUE}, otherwise
#'\code{NULL}; formatting of \code{rt} is reused.}
#'\item{\code{nonpar_model}:}{the estimation object returned by \code{\link[smoots]{tsmooth}} or \code{\link[esemifar]{tsmoothlm}} for
#'\code{use_nonpar = TRUE}.}
#'\item{\code{trunc}:}{the input argument \code{trunc}.}
#'}
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- garchm_estim(rt, model = "garch")
#'model
#'

garchm_estim <- function(rt,
                        model = c("garch", "gjrgarch", "tgarch", "aparch", "figarch", "figjrgarch", "fitgarch", "fiaparch"),
                        orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        drange = c(0, 1),
                        meanspec = mean_spec(), Drange = c(0, 1),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        control_nonpar = list(), mean_after_nonpar = FALSE,
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1),
                        trunc = "none", presample = 50, Prange = c(1, 5)) {


  model <- match.arg(model)

  fun <- switch(
    model,
    "garch" = garch,
    "gjrgarch" = gjrgarch,
    "tgarch" = tgarch,
    "aparch" = aparch,
    "figarch" = figarch,
    "figjrgarch" = figjrgarch,
    "fitgarch" = fitgarch,
    "fiaparch" = fiaparch
  )

  args <- list(
    rt = rt,
    orders = orders,
    cond_dist = cond_dist,
    drange = drange,
    meanspec = meanspec,
    nonparspec = nonparspec,
    use_nonpar = use_nonpar,
    n_test = n_test,
    start_pars = start_pars,
    LB = LB, UB = UB,
    control = control,
    control_nonpar = control_nonpar,
    mean_after_nonpar = mean_after_nonpar,
    parallel = parallel,
    ncores = ncores,
    trunc = trunc,
    presample = presample,
    Prange = Prange
  )

  if (model %in% c("garch", "gjrgarch", "tgarch", "aparch")) {
    args[["drange"]] <- NULL
  }

  do.call(what = fun, args = args)

}
