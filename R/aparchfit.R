aparch_grab_pars <- function(theta, p, q, incl_mean, extra_par, skew_par) {

  imean <- !incl_mean
  mu <- if (incl_mean) theta[[1]] else 0
  step0 <- 2 - imean
  omega <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  step3 <- step2 + p
  phi <- theta[(step0 + 1):step1]
  beta <- theta[(step1 + 1):step2]
  gamma <- theta[(step2 + 1):step3]
  delta <- theta[[step3 + 1]]
  shape <- if (extra_par) theta[[step3 + 2]] else 0
  skew <- if (skew_par) theta[[step3 + 2 + extra_par]] else 0

  list(
    mu = mu,
    omega = omega,
    phi = phi,
    beta = beta,
    gamma = gamma,
    delta = delta,
    shape = shape,
    skew = skew
  )

}

sigt_aparch_R <- function(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par,
              trunc, presample, ...) {

  list2env(list(...), envir = environment())

  all_pars <- aparch_grab_pars(
    theta = theta,
    p = p,
    q = q,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par
  )
  mu <- all_pars$mu
  omega <- all_pars$omega
  phi <- all_pars$phi
  beta <- all_pars$beta
  gamma <- all_pars$gamma
  delta <- all_pars$delta
  shape <- all_pars$shape
  skew <- all_pars$skew

  sigt <- c(sigt_aparch_Cpp(
    x = x,
    mu = mu,
    phi = phi,
    beta = beta,
    gamma = gamma,
    delta = delta,
    omega = omega,
    sigd_init,
    etransf_init
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape, cmeans = rep(mu, length(sigt)))

}

goal_fun_creator_aparch <- function(
    theta,
    x,
    pdf_fun2,
    p,
    q,
    p_ar,
    q_ma,
    incl_mean,
    extra_par,
    skew_par,
    sig_fun,
    sigd_init,
    etransf_init,
    trunc,
    presample
) {

  res_list <- sig_fun(
    theta = theta,
    x = x,
    p = p,
    q = q,
    p_ar = p_ar,
    q_ma = q_ma,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par,
    trunc = trunc,
    presample = presample,
    sigd_init = sigd_init,
    etransf_init = etransf_init
  )

  nllhood_calc(res_list, x, pdf_fun2)

}


#'APARCH Model Fitting
#'
#'Fit an APARCH model under the six
#'most common and further conditional distributions to observed data
#'using quasi maximum-likelihood estimation.
#'
#'@param rt the observed series ordered from past to present; can be
#'a numeric vector, a \code{"zoo"} class time series object, or a
#'\code{"ts"} class time series object.
#'@param orders a two-element numeric vector containing the two model
#'orders \eqn{p} and \eqn{q} (see Details for more information); currently,
#'only the default \code{orders = c(1, 1)} is supported; other specifications
#'of a two-element numeric vector will lead to \code{orders = c(1, 1)} being
#'run and a warning message being returned.
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
#'an internally saved default set of values is used; see "Details" for the
#'order of elements; elements should be set with respect to a series rescaled
#'to have sample variance one.
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
#'an internally saved default set of values is used; see "Details" for the
#'order of elements; elements should be set with respect to a series rescaled
#'to have sample variance one.
#'@param UB the upper boundaries of the parameters in the numerical optimization
#'routine; should be of the same length as the parameter output vector
#'within the output object (also keeping the same order); for \code{NULL},
#'an internally saved default set of values is used; see "Details" for the
#'order of elements; elements should be set with respect to a series rescaled
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
#'@param fix_delta let the parameter \eqn{\delta} be either free (\code{fix_delta = NA}),
#'fix it to 1 (\code{fix_delta = 1}) or to 2 (\code{fix_delta = 2}); the latter
#'two are specific submodels of a FIAPARCH.
#'@param skip_vcov a logical indicating whether or not to skip the computation
#'of the variance-covariance matrix of the parameter estimators and therefore
#'also standard error computation.
#'
#'@details
#'Let \eqn{\left\{r_t\right\}}, with \eqn{t \in \mathbb{Z}} as the
#'time index, be a theoretical time series that follows
#'\deqn{r_t=\mu+\varepsilon_t \text{ with } \varepsilon_t=\sigma_t \eta_t \text{ and } \eta_t \sim \text{IID}(0,1), \text{ where}}
#'\deqn{\sigma_t^{\delta}=\omega+\sum_{i=1}^{p}\phi_i\left(\left|\varepsilon_{t-i}\right|-\gamma_{i}\varepsilon_{t-i}\right)^{\delta} + \sum_{j=1}^{j} \beta_{j}\sigma_{t-j}^{\delta}.}
#'Here, \eqn{\eta_t\sim\text{IID}(0,1)} means that the innovations
#'\eqn{\eta_t} are independent and identically distributed (iid) with mean zero
#'and variance one, whereas \eqn{\sigma_t > 0} are the conditional standard
#'deviations in \eqn{r_t}.
#'Moreover, \eqn{\beta_j}, \eqn{j=1,2,\dots, q}, \eqn{\phi_i}, \eqn{i=1,2,\dots, p},
#'are real-valued coefficients. \eqn{p}
#'and \eqn{q} are the model orders definable through the argument \code{orders},
#'where \eqn{p} is the first element and \eqn{q} is the second element in the
#'argument. In addition, we have \eqn{\mu = E\left(r_t\right)} as a
#'real-valued parameter, and \eqn{\gamma \in (-1,1)}. \eqn{\omega > 0} is the intercept. It is assumed that
#'all \eqn{\beta_j} and \eqn{\phi_i} are non-negative.
#'
#'See also the reference section for sources on the APARCH (Ding et al., 1993)
#'and the FIAPARCH (Tse, 1998) models.
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
#'An ARMA-APARCH or a FARIMA-APARCH can be fitted by adjusting the
#'argument \code{meanspec} correspondingly.
#'
#'As an alternative, a semiparametric extension of the pure models
#'in the conditional variance can be implemented. If \code{use_nonpar = TRUE},
#'\code{meanspec} is omitted and before fitting a model in the
#'conditional volatility following the remaining function arguments, a smooth scale function,
#'i.e. a function representing the unconditional standard deviation over time,
#'is being estimated following the specifications in \code{nonparspec} and
#'\code{control_nonpar}. This preliminary step stabilizes the input
#'series \code{rt}, as long-term changes in the unconditional variance
#'are being estimated and removed before the parametric step using
#'\code{\link[smoots]{tsmooth}}. \code{control_nonpar} can be adjusted following
#'to make changes to the arguments of \code{\link[smoots]{tsmooth}}
#'for short-memory specifications. These arguments specify settings
#'for the automated bandwidth selection algorithms implemented by this
#'function. By default, we use the settings
#'\code{InfR = "Nai"},
#'\code{bStart = 0.15}, \code{cb = 0.05}, and
#'\code{method = "lpr"} for \code{\link[smoots]{tsmooth}}.
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
#'being estimated using
#'\code{\link[smoots]{tsmooth}} and denoted here by
#'\eqn{\hat{m}(x_t)}. Then from \eqn{\hat{\xi}_t = y_t - \hat{m}(x_t)}
#'obtain \eqn{\hat{C} = -\ln\left\{\sum_{t=1}^{n}\exp\left(\hat{\xi}_t\right)\right\}},
#'and obtain the estimated scale function as
#'\eqn{\hat{s}(x_t)=\exp\left[\left(\hat{\mu}(x_t) - \hat{C}\right) / 2\right]}.
#'The stabilized / standardized version of the series \eqn{\left\{r_t\right\}}
#'is then \eqn{\tilde{r}_t = r_t^{*} / \hat{s}(x_t)}, to which
#'a purely parametric volatility model following the remaining function arguments
#'is then
#'fitted. The estimated volatility at a given time point is then
#'the product of the estimate of the corresponding scale function value
#'and of the estimated conditional standard deviation (following the parametric
#'model part) for that same time point. See for example Feng et al. (2022)
#'or Letmathe et al. (2023) for more information on the semiparametric extension
#'of volatility models.
#'
#'The order for manual settings of \code{start_pars}, \code{LB} and \code{UB}
#'is crucial. The correct order is: \eqn{\mu}, \eqn{\text{ar}_1,\dots,\text{ar}_{p^{*}}},
#'\eqn{\text{ma}_1,\dots,\text{ma}_{q^{*}}},\eqn{D},\eqn{\omega},
#'\eqn{\phi_1,\dots,\phi_p}, \eqn{\beta_1, \dots, \beta_{q}},
#'\eqn{\gamma_1,\dots,\gamma_p}, \eqn{\delta},
#'shape parameter,
#'skewness parameter. Depending on the exact model specification,
#'parameters irrelevant for the specification at hand should be dropped
#'in \code{start_pars}, \code{LB} and \code{UB}.
#'
#'@export
#'
#'@return
#'An object of S4-class \code{"fEGarch_fit_aparch"}
#'is returned. It contains the following elements.
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
#'\item{\code{nonpar_model}:}{the estimation object returned by \code{\link[smoots]{tsmooth}} for
#'\code{use_nonpar = TRUE}.}
#'\item{\code{trunc}:}{the input argument \code{trunc}.}
#'}
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- aparch(rt)
#'model
#'
#'@references
#'\itemize{
#'\item{Ding, Z., Granger, C. W. J., & Engle, R. F. (1993). A long memory property of stock market returns
#'and a new model. Journal of Empirical Finance, 1(1):
#'83-106. DOI: 10.1016/0927-5398(93)90006-D.}
#'\item{Feng, Y., Gries, T., Letmathe, S., & Schulz, D. (2022). The smoots Package in R for Semiparametric Modeling of
#'Trend Stationary Time Series. The R Journal,
#'14(1), 182-195. URL: https://journal.r-project.org/articles/RJ-2022-017/.}
#'\item{Letmathe, S., Beran, J., & Feng, Y. (2023). An extended exponential SEMIFAR model with application
#'in R. Communications in Statistics - Theory and Methods,
#'53(22), 7914–7926. DOI: 10.1080/03610926.2023.2276049.}
#'\item{Tse, Y. K. (1998). The conditional heteroskedasticity of the
#'yen-dollar exchange rate. Journal of Applied Econometrics, 13(1):
#'49-55. DOI: 10.1002/(SICI)1099-1255(199801/02)13:1<49::AID-JAE459>3.0.CO;2-O.}
#'}
#'
aparch <- function(rt, orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
                        meanspec = mean_spec(), Drange = c(0, 1),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL, LB = NULL, UB = NULL, control = list(),
                        control_nonpar = list(), mean_after_nonpar = FALSE,
                        parallel = TRUE, ncores = max(1, future::availableCores() - 1),
                        trunc = "none", presample = 50, Prange = c(1, 5),
                        fix_delta = c(NA, 1, 2), skip_vcov = FALSE) {

  inp_args <- mget(names(formals()), envir = environment())

  if (length(fix_delta) == 3 && (is.na(fix_delta[[1]]) && fix_delta[[2]] == 1 && fix_delta[[3]] == 2)) {
    fix_delta <- NA
  }

  inp_args[["fix_delta"]] <- NULL

  mt <- if (is.na(fix_delta)) {
    "aparch"
  } else if (fix_delta == 1) {
    "tgarch"
  } else if (fix_delta == 2) {
    "gjrgarch"
  }

  inp_args[["model_type"]] <- mt
  do.call(general_garch_fitting, args = inp_args)

}

#'@export
#'@rdname show-methods
#'@aliases show,fEgarch_fit_aparch-method
setMethod(
  "show",
  "fEGarch_fit_aparch",
  function(object) {
  x <- object
  par_names <- names(x@pars)
  se <- unname(x@se)
  pars <- unname(x@pars)

  tvals <- pars / se
  atvals <- abs(tvals)
  pvals <- 2 * (1 - pnorm(atvals))

  df <- data.frame(
    par = sprintf("%.4f", pars),
    se = sprintf("%.4f", se),
    tval = sprintf("%.4f", tvals),
    pval = sprintf("%.4f", pvals)
  )
  row.names(df) <- par_names

  check_nonpar <- !is.null(x@nonpar_model)
  b <- if (is.null(x@nonpar_model$b0)) {
    x@nonpar_model$b
  } else {
    x@nonpar_model$b0
  }
  bwidth <- if (check_nonpar) {
    paste0(
      " (poly_order = ", x@nonpar_model$p,"; bandwidth = ", sprintf("%.4f", b), ")"
    )
  } else {
    ""
  }

  part1 <- paste0(
    "*************************************\n",
    "*        Fitted APARCH Model        *\n",
    "*************************************\n",
    " \n",
    "Type: aparch\n",
    "Orders: (", x@orders[[1]], ", ", x@orders[[2]], ")\n",
    "Long memory: ", x@long_memo, "\n",
    "Cond. distribution: ", x@cond_dist, "\n",
    " \n",
    "ARMA orders (cond. mean): (", x@meanspec@orders[[1]], ", ", x@meanspec@orders[[2]], ")\n",
    "Long memory (cond. mean): ", x@meanspec@long_memo, "\n",
    " \n",
    "Scale estimation: ", paste0(check_nonpar, bwidth), "\n",
    " \n",
    "Fitted parameters:\n"
  )

  part2 <- paste0(
    " \n",
    "Information criteria (parametric part):\n",
    "AIC: ", sprintf("%.4f", unname(x@inf_criteria[[1]])),
    ", BIC: ", sprintf("%.4f", unname(x@inf_criteria[[2]]))
  )

  cat(part1, fill = TRUE)
  print(df)
  cat(part2, fill = TRUE)
  }
)
