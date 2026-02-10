nonpar_est <- function(rt, lm, nonparspec, n_test, control_nonpar) {

    split_data <- split_ts(rt, n_test = n_test)
    rt_train <- split_data$train
    rt_test <- split_data$test

    bw <- nonparspec@bwidth
    if (is.null(bw)) {
      args <- list(smoots_default, esemifar_default)[[lm + 1]]
      smooth_fun <- list(smoots::tsmooth, esemifar::tsmoothlm)[[lm + 1]]
      if (length(control_nonpar) > 0) {
        names_l <- names(control_nonpar)
        for (i in names_l) {
          args[[i]] <- control_nonpar[[i]]
        }
      }
    } else {
      args <- list()
      smooth_fun <- smoots::gsmooth
      args[["b"]] <- bw
      args[["v"]] <- 0
    }
    args[["p"]] <- nonparspec@poly_order
    args[["mu"]] <- nonparspec@kernel_order
    args[["bb"]] <- list("shorten" = 0, "extend" = 1)[[nonparspec@boundary_method]]
    mu <- mean(rt_train)
    rt_dm <- rt_train - mu
    yt <- zoo::coredata(log(rt_dm^2))
    args[["y"]] <- yt
    est_nonpar <- do.call(smooth_fun, args = args)
    res <- est_nonpar$res
    C <- -log(mean(exp(res)))
    scale_fun <- exp((est_nonpar$ye - C) / 2)
    rt_orig <- rt_train
    et <- rt_dm / scale_fun
    list(
      et = et,
      scale_fun = scale_fun,
      test_obs = if (n_test > 0) {rt_test} else if (n_test == 0) {NULL},
      rt_train = rt_train,
      est_nonpar = est_nonpar,
      mu = mu
    )
}

#'Fitting Function for Models of the Broader EGARCH Family
#'
#'Use quasi-maximum-likelihood estimation to fit a model from the
#'broader EGARCH family to some observed time series.
#'
#'@param spec an S4 object of class \code{"egarch_type_spec"}
#'or \code{"loggarch_type_spec"} as returned by the various spec
#'functions of this package, for example \link[fEGarch]{fEGarch_spec}
#'or its wrappers like \link[fEGarch]{megarch_spec} or
#'\link[fEGarch]{mloggarch_spec}, among others.
#'@param rt the input time series to fit the model to ordered from
#'past to present; can also be a \code{"zoo"} class object or a
#'\code{"ts"} class object.
#'@param drange a two-element numeric vector that indicates the boundaries
#'of the interval over which to search for the fractional differencing
#'parameter \eqn{d} in a long-memory GARCH-type model in the conditional
#'volatility model part; by default,
#'\eqn{d} being searched for on the
#'interval from 0 to 1; note that specific
#'settings in the arguments
#'\code{LB} and \code{UB} overwrite this argument.
#'@param meanspec an object of class "mean_spec"; indicates the
#'specifications for the model in the conditional mean.
#'@param Drange a two-element numeric vector that indicates the boundaries
#'of the interval over which to search for the fractional differencing
#'parameter \eqn{d} in a long-memory ARMA-type model in the conditional
#'mean model part; by default,
#'\eqn{D} being searched for on the
#'interval from 0 to 1; note that specific
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
#'@param start_pars a vector with starting parameters for the
#'optimization; must be of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of starting values; see "Details" for the
#'order of elements.
#'@param LB a vector with lower boundaries for parameters; must be
#'of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of lower boundary values; see "Details" for the
#'order of elements.
#'@param UB a vector with upper boundaries for parameters; must be
#'of the same length as the output vector of
#'parameters; the default NULL uses internally saved default sets
#'of upper boundary values; see "Details" for the
#'order of elements.
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
#'@param skip_vcov a logical indicating whether or not to skip the computation
#'of the variance-covariance matrix of the parameter estimators and therefore
#'also standard error computation.
#'
#'@export
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib fEGarch
#'
#'@return
#'An object of S4-class \code{"fEGarch_fit_egarch"} \code{"fEGarch_fit_loggarch"}
#'is returned depending on the selected model type in the model
#'specification. It contains the following elements.
#'\describe{
#'\item{\code{pars}:}{a named numeric vector with the parameter estimates.}
#'\item{\code{se}:}{a named numeric vector with the obtained standard errors in accordance with the parameter estimates.}
#'\item{\code{vcov_mat}:}{the variance-covariance matrix of the parameter estimates with named columns and rows.}
#'\item{\code{rt}:}{the input object \code{rt} (or at least the training data, if \code{n_test} is greater than zero);
#'if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is kept.}
#'\item{\code{sigt}:}{the estimated conditional standard deviations (or for \code{use_nonpar = TRUE} the estimated total
#'volatilities, i.e. scale function value times conditional standard deviation); if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{sigt}.}
#'\item{\code{cmeans}:}{the estimated conditional means; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{cmeans}.}
#'\item{\code{etat}:}{the obtained residuals; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{etat}.}
#'\item{\code{orders}:}{a two-element numeric vector stating the considered model orders.}
#'\item{\code{cond_dist}:}{a character value stating the conditional distribution considered in the model fitting.}
#'\item{\code{long_memo}:}{a logical value stating whether or not long memory was considered in the model fitting.}
#'\item{\code{llhood}:}{the log-likelihood value obtained at the optimal parameter combination.}
#'\item{\code{inf_criteria}:}{a named two-element numeric vector with the corresponding AIC (first element) and BIC (second element) of the fitted parametric model part; for purely parametric models, these criteria are valid for the entire model; for semiparametric models, they are only valid for the parametric step and are not valid for the entire model.}
#'\item{\code{powers}:}{a two-element numeric vector stating the powers considered for the asymmetry term (first element) and the magnitude term (second element); only exists, if a type I model was fitted.}
#'\item{\code{modulus}:}{a two-element logical vector stating the modulus transformations (\code{TRUE}, otherwise \code{FALSE}) considered for the asymmetry term (first element) and the magnitude term (second element); only exists, if a type I model was fitted.}
#'\item{\code{meanspec}:}{the settings for the model in the conditional mean; is an object
#'of class \code{"mean_spec"} that is identical to the object passed to the input argument
#'\code{meanspec}.}
#'\item{\code{test_obs}:}{the observations at the end up the input \code{rt} reserved for
#'testing following \code{n_test}.}
#'\item{\code{scale_fun}:}{the estimated scale function values, if \code{use_nonpar = TRUE}, otherwise
#'\code{NULL}; formatting of \code{rt} is reused.}
#'\item{\code{nonpar_model}:}{the estimation object returned by either
#'\code{\link[smoots]{tsmooth}} or \code{\link[esemifar]{tsmoothlm}} for
#'\code{use_nonpar = TRUE}.}
#'\item{\code{trunc}:}{the input argument \code{trunc}.}
#'}
#'
#'@details
#'For details on the models in the conditional variance, see \link[fEGarch]{fEGarch_spec}.
#'For details on the models in the conditional mean, see \link[fEGarch]{mean_spec}.
#'The combined model defined through \link[fEGarch]{mean_spec} and
#'\link[fEGarch]{fEGarch_spec} is the specified model. It can be thought of
#'as the model described in \link[fEGarch]{mean_spec} with \eqn{\left\{r_t\right\}}
#'therein being governed by a model from the EGARCH family (see for example Feng et al., 2025) as described in
#'\link[fEGarch]{fEGarch_spec}, however with mean of \eqn{\left\{r_t\right\}} fixed
#'to zero.
#'
#'The specified model is then fitted using quasi maximum likelihood estimation,
#'where pre-sample values of \eqn{g\left(\eta_{t-1}\right)} are
#'filled in through its theoretical expectation of zero, which is analogous to not
#'setting pre-sample values in the long-memory case. In addition, in
#'short-memory models, pre-sample values of \eqn{ln(\sigma_t^2)} are roughly
#'approximated through \eqn{\ln\left(\hat{\sigma}_{t}^{2}\right)}, where
#'\eqn{\hat{\sigma}_{t}^{2}} is the sample variance of the observations.
#'
#'See the references section for sources on the EGARCH (Nelson, 1991),
#'FIEGARCH (Bollerslev and Mikkelsen, 1996),
#'Log-GARCH (Geweke, 1986; Pantula, 1986; Milhoj, 1987)
#'and FILog-GARCH (Feng et al., 2020) models. For information on the
#'FIMLog-GARCH, see Feng et al. (2023).
#'
#'In the current package version, standard errors of parameter estimates are
#'computed from the Hessian at the optimum of the log-likelihood using
#'\code{\link[numDeriv]{hessian}}. To ensure numerical stability and
#'applicability to a huge variety of differently scaled data, parametric
#'models are first fitted to data that is scaled to have sample variance
#'1. Parameter estimates and other quantities are then either
#'retransformed or recalculated afterwards for the original data.
#'
#'For a conditional average Laplace distribution, an optimal model for each
#'distribution parameter \eqn{P} from 1 to 5 is estimated (assuming that
#'\eqn{P} is then fixed to the corresponding value). Afterwards, \eqn{P} is then
#'estimated by selecting the estimated model among the five fitted models that
#'has the largest log-likelihood. The five models are, by default, fitted
#'simultaneously using parallel programming techniques (see also the arguments
#'\code{parallel} and \code{ncores}, which are only relevant for a conditional
#'average Laplace distribution). After the optimal model (including
#'the estimate of \eqn{P} called \eqn{\hat{P}}) has been determined, \eqn{P=\hat{P}}
#'is seen as fixed to obtain the standard errors via the Hessian matrix for the
#'estimates of the continuous parameters. A standard error for \eqn{\hat{P}} is therefore
#'not obtained and the ones obtained for the remaining estimates do not account
#'for \eqn{\hat{P}}.
#'
#'As an alternative, a semiparametric extension of the pure models
#'in the conditional variance can be implemented. If \code{use_nonpar = TRUE},
#'\code{meanspec} is omitted and before fitting a zero-mean model in the
#'conditional volatility following \code{spec}, a smooth scale function,
#'i.e. a function representing the unconditional standard deviation over time,
#'is being estimated following the specifications in \code{nonparspec} and
#'\code{control_nonpar}. This preliminary step stabilizes the input
#'series \code{rt}, as long-term changes in the unconditional variance
#'are being estimated and removed before the parametric step using
#'either \code{\link[smoots]{tsmooth}} or \code{\link[esemifar]{tsmoothlm}}
#'depending on whether \code{spec} specifies a model with short memory
#'or with long memory. \code{control_nonpar} can be adjusted following
#'the arguments of \code{\link[smoots]{tsmooth}} for short-memory
#'specifications of \code{spec}, on the other hand changes to
#'arguments of \code{\link[esemifar]{tsmoothlm}} can be passed to it
#'for long-memory specifications. These arguments specify settings
#'for the automated bandwidth selection algorithms implemented by these
#'two functions. By default, we use the settings
#'\code{Mcf = "NP"}, \code{InfR = "Opt"},
#'\code{bStart = 0.15}, \code{bvc = "Y"}, \code{cb = 0.05},
#'and \code{method = "lpr"} within the call to \code{\link[smoots]{tsmooth}},
#'as well as \code{pmin = 0}, \code{pmax = 1}, \code{qmin = 0},
#'\code{qmax = 1}, \code{InfR = "Opt"},
#'\code{bStart = 0.15}, \code{cb = 0.05}, and
#'\code{method = "lpr"} for \code{\link[esemifar]{tsmoothlm}}.
#'\code{\link{locpol_spec}} passed to \code{nonparspec} handles
#'more direct settings of the local polynomial smoother itself. See
#'the documentation for these functions to get a detailed overview
#'of these settings. Assume \eqn{\{r_t\}} to be the observed series, where
#'\eqn{t = 1, 2, \dots, n},
#'then \eqn{r_t^{*} = r_t - \bar{r}}, with \eqn{\bar{r}} being the arithmetic
#'mean over the observed \eqn{r_t}, is computed and subsequently
#'\eqn{y_t = \ln\left[\left(r_t^{*}\right)^2\right]}. The subtraction of
#'\eqn{\bar{r}} is necessary so that \eqn{r_t^{*}} are all different from zero
#'almost surely. Once \eqn{y_t} are available, its trend \eqn{m(x_t)},
#'with \eqn{x_t} as the rescaled time on the interval \eqn{[0, 1]}, is
#'being estimated using either \code{\link[smoots]{tsmooth}} or
#'\code{\link[esemifar]{tsmoothlm}} and denoted here by
#'\eqn{\hat{m}(x_t)}. Then from \eqn{\hat{\xi}_t = y_t - \hat{m}(x_t)}
#'obtain \eqn{\hat{C} = -\ln\left\{\sum_{t=1}^{n}\exp\left(\hat{\xi}_t\right)\right\}},
#'and obtain the estimated scale function as
#'\eqn{\hat{s}(x_t)=\exp\left[\left(\hat{\mu}(x_t) - \hat{C}\right) / 2\right]}.
#'The stabilized / standardized version of the series \eqn{\left\{r_t\right\}}
#'is then \eqn{\tilde{r}_t = r_t^{*} / \hat{s}(x_t)}, to which
#'a purely parametric volatility model following \code{spec} is then
#'fitted. The estimated volatility at a given time point is then
#'the product of the estimate of the corresponding scale function value
#'and of the estimated conditional standard deviation (following the parametric
#'model part) for that same time point. See for example Feng et al. (2022)
#'or Letmathe et al. (2023) for more information on the semiparametric extension
#'of volatility models. Moreover, if \code{bwidth} in the object passed to \code{nonparspec}
#'is not at its default \code{NULL} but instead a numeric value between \code{0} and \code{0.5},
#'the automated bandwidth selection is skipped and the provided bandwidth in \code{bwidth}
#'is used.
#'
#'The order for manual settings of \code{start_pars}, \code{LB} and \code{UB}
#'is crucial. The correct order is: \eqn{\mu}, \eqn{\text{ar}_1,\dots,\text{ar}_{p^{*}}},
#'\eqn{\text{ma}_1,\dots,\text{ma}_{q^{*}}},\eqn{D},\eqn{\omega_{\sigma}},
#'\eqn{\phi_1,\dots,\phi_p}, \eqn{\psi_1, \dots, \psi_{q-1}}, \eqn{\kappa}, \eqn{\gamma},
#'\eqn{d}, \eqn{\text{shape parameter}},
#'\eqn{\text{skewness parameter}} for Type I models (see \code{\link{fEGarch_spec}}). For Type
#'II models, we have \eqn{\mu}, \eqn{\text{ar}_1,\dots,\text{ar}_{p^{*}}},
#'\eqn{\text{ma}_1,\dots,\text{ma}_{q^{*}}},\eqn{D},\eqn{\omega_{\sigma}},
#'\eqn{\phi_1,\dots,\phi_p}, \eqn{\psi_1, \dots, \psi_{q}}, \eqn{d},
#'shape parameter,
#'skewness parameter. Depending on the exact model specification,
#'parameters irrelevant for the specification at hand should be dropped
#'in \code{start_pars}, \code{LB} and \code{UB}.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'# Pure conditional volatility model
#'spec <- fEGarch_spec()
#'model <- fEGarch(spec, rt)
#'model
#'# Simultaneously model conditional mean
#'spec <- egarch_spec()
#'model2 <- suppressWarnings(
#'  fEGarch(spec, rt, meanspec = mean_spec(orders = c(1, 1)))
#')
#'model2
#'
#'
#'@references
#'\itemize{
#'\item{Bollerslev, T., & Mikkelsen, H. O. (1996). Modeling and pricing long memory in stock market volatility. Journal of Econometrics,
#'73(1), 151–184. DOI: 10.1016/0304-4076(95)01749-6.}
#'\item{Feng, Y., Beran, J., Ghosh, S., & Letmathe, S. (2020). Fractionally integrated Log-GARCH with application to value at risk and expected shortfall.
#'Working Papers CIE No. 137, Paderborn University, Center for International Economics.
#'URL: http://groups.uni-paderborn.de/wp-wiwi/RePEc/pdf/ciepap/WP137.pdf.}
#'\item{Feng, Y., Gries, T., Letmathe, S., & Schulz, D. (2022). The smoots Package in R for Semiparametric Modeling of
#'Trend Stationary Time Series. The R Journal,
#'14(1), 182-195. URL: https://journal.r-project.org/articles/RJ-2022-017/.}
#'\item{Feng, Y., Gries, T., & Letmathe, S. (2023). FIEGARCH, modulus asymmetric FILog-GARCH
#'and trend-stationary dual long memory time series. Working Papers CIE No. 156, Paderborn University.
#'URL: https://econpapers.repec.org/paper/pdnciepap/156.htm.}
#'\item{Feng, Y., Peitz, C., & Siddiqui, S. (2025). A few useful members of the EGARCH-family
#'with short- or long-memory in volatility. Unpublished working paper at Paderborn University.}
#'\item{Geweke, J. (1986). Modeling the persistence of conditional variances: A comment. Econometric Reviews, 5(1),
#'57-61. DOI: 10.1080/07474938608800088.}
#'\item{Letmathe, S., Beran, J., & Feng, Y. (2023). An extended exponential SEMIFAR model with application
#'in R. Communications in Statistics - Theory and Methods,
#'53(22), 7914–7926. DOI: 10.1080/03610926.2023.2276049.}
#'\item{Milhoj, A. (1987). A Multiplicative Parameterization of ARCH Models. University of Copenhagen, Denmark.}
#'\item{Nelson, D. B. (1991). Conditional Heteroskedasticity in Asset Returns: A New Approach. Econometrica,
#'59(2), 347–370. DOI: 10.2307/2938260.}
#'\item{Pantula, S. G. (1986). Modeling the persistence of conditional variances: A comment. Econometric Reviews, 5(1),
#'71-74. DOI: 10.1080/07474938608800089.}
#'}

fEGarch <- function(spec = egarch_spec(),
    rt, drange = c(0, 1),
    meanspec = mean_spec(), Drange = c(0, 1),
    nonparspec = locpol_spec(), use_nonpar = FALSE,
    n_test = 0,
    start_pars = NULL, LB = NULL, UB = NULL, control = list(),
    control_nonpar = list(), mean_after_nonpar = FALSE,
    parallel = TRUE, ncores = max(1, future::availableCores() - 1),
    trunc = "none", presample = 50, Prange = c(1, 5), skip_vcov = FALSE) {

  if (use_nonpar) {

    nonpar_out <- run_nonpar(nonparspec = nonparspec,
               control_nonpar = control_nonpar, rt = rt, n_test = n_test,
               mean_after_nonpar = mean_after_nonpar, lm = spec@long_memo)
    rt <- nonpar_out$rt
    meanspec <- nonpar_out$meanspec
    nonpar_result <- nonpar_out$nonpar_result
    n_test <- nonpar_out$n_test

  }

  out <- fEGarch_fit(
    spec = spec,
    rt = rt,
    drange = drange,
    meanspec = meanspec,
    Drange = Drange,
    n_test = n_test,
    start_pars = start_pars,
    LB = LB,
    UB = UB,
    control = control,
    parallel = parallel,
    ncores = ncores,
    trunc = trunc,
    presample = presample,
    Prange = Prange,
    skip_vcov = skip_vcov
  )

  if (use_nonpar) {
    mu <- rep(nonpar_result$mu, length(rt))
    out@rt <- nonpar_result$rt_train
    format_list <- format_applier_ts(
      rt = nonpar_result$rt_train,
      list_of_ts = list(
        scale_fun = nonpar_result$scale_fun,
        cmeans = mu
      )
    )
    out@scale_fun <- format_list$scale_fun
    out@cmeans <- format_list$cmeans
    out@test_obs <- nonpar_result$test_obs
    out@sigt <- out@sigt * format_list$scale_fun
    out@nonpar_model <- nonpar_result$est_nonpar
  }

  out
}

