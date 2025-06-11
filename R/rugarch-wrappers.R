rugarch_modelfit <- function(rt, model = "apARCH", orders = c(1, 1),
                        cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
                        include_mean = TRUE,
                        meanspec = mean_spec(),
                        nonparspec = locpol_spec(), use_nonpar = FALSE,
                        n_test = 0,
                        start_pars = NULL,
                        control_nonpar = list()) {

  if (use_nonpar) {

    nonpar_result <- nonpar_est(
      rt = rt,
      lm = FALSE,   # For now only under short memory
      nonparspec = nonparspec,
      n_test = n_test,
      control_nonpar = control_nonpar
    )

    meanspec <- mean_spec(
      orders = c(0, 0),
      long_memo = FALSE,
      include_mean = FALSE
    )    # Force to not estimate the mean in parametric part
    n_test <- 0
    rt <- nonpar_result$et

  }

  split_rt <- split_ts(rt, n_test)
  rt_core <- zoo::coredata(split_rt$train)

  cond_dist <- match.arg(cond_dist)

  if (is.null(start_pars)) {
    start.pars <- list()
  } else {
    start.pars <- start_pars
  }

  if (model == "gjrGARCH") {
    fixed.pars <- list(delta = 2)
    model_s <- "apARCH"
  } else {
    fixed.pars <- list()
    model_s <- model
  }

  spec <- rugarch::ugarchspec(
    variance.model = list(model = model_s, garchOrder = orders),
    mean.model = list(armaOrder = meanspec@orders, include.mean = include_mean,
                      arfima = meanspec@long_memo),
    distribution.model = cond_dist, start.pars = start.pars,
    fixed.pars = fixed.pars
  )

  est <- rugarch::ugarchfit(spec = spec, data = rt_core, out.sample = n_test)
  sig <- est@fit$sigma
  cmeans <- est@fit$fitted.values
  etat <- est@fit$residuals
  coef <- est@fit$coef
  if (model == "gjrGARCH") {
    coef <- coef[names(coef) != "delta"]
  }
  serrors <- est@fit$se.coef
  names(serrors) <- names(coef)
  vcov_mat <- est@fit$cvar
  llh <- est@fit$LLH
  inf_crit <- rugarch::infocriteria(est)[c(1, 2)]
  names(inf_crit) <- c("aic", "bic")

  series_out <- format_applier_ts(
    rt = split_rt$train,
    list_of_ts = list(
      "cmeans" = cmeans,
      "sigt" = sig,
      "etat" = etat
    )
  )
  cmeans <- series_out$cmeans
  sigt <- series_out$sigt
  etat <- series_out$etat

  out <- fEGarch_fit_rugarch_wrapper(
    rt = split_rt$train,
    sigt = sig,
    cmeans = cmeans,
    etat = etat,
    pars = coef,
    se = serrors,
    scale_fun = NULL,
    vcov_mat = vcov_mat,
    llhood = llh,
    orders = orders,
    cond_dist = cond_dist,
    inf_criteria = inf_crit,
    rugarch_model = est,
    long_memo = FALSE,
    model = model,
    meanspec = meanspec,
    test_obs = split_rt$test,
    nonpar_model = NULL,
    trunc = NA
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

#'Wrapper Functions for Selected \code{rugarch} GARCH Models
#'
#'Easy to use functions for fitting selected GARCH-type models from
#'the widely-known \code{rugarch} package by Ghalanos (2024). These
#'functions are meant as an easy to use way to compare the main results from
#'\code{rugarch} to the newly established models in the \code{fEGarch} package
#'and are by no means considered to be a replacement of \code{rugarch}.
#'
#'@param rt the input time series to fit the model to ordered from
#'past to present; can also be a \code{"zoo"} class object or a
#'\code{"ts"} class object.
#'@param orders the ARCH and GARCH orders of the model as a two-element
#'numeric vector.
#'@param cond_dist a single-element character vector with the conditional
#'distribution to consider.
#'@param meanspec an object of class "mean_spec"; indicates the
#'specifications for the model in the conditional mean.
#'@param n_test a single numerical value indicating, how many observations
#'at the end of \code{rt} not to include in the fitting process and to
#'reserve for backtesting.
#'@param start_pars a named list with starting parameters.
#'@param nonparspec an object of class \code{"locpol_spec"} returned
#'by \code{\link{locpol_spec}}; defines the settings of the nonparametric
#'smoothing technique for \code{use_nonpar = TRUE}.
#'@param use_nonpar a logical indicating whether or not to implement a
#'semiparametric extension of the volatility model defined through \code{spec};
#'see "Details" for more information.
#'@param control_nonpar a list containing changes to the arguments
#'for the hyperparameter estimation algorithm in the nonparametric
#'scale function estimation for
#'\code{use_nonpar = TRUE}; see "Details" for more information.
#'
#'@details
#'For most details, please see the documentation of the \code{rugarch} package
#'(Ghalanos, 2024).
#'
#'These functions also provide an extension, so that a nonparametric,
#'smooth scale function in the unconditional standard deviation
#'can be estimated before the parametric step.
#'If \code{use_nonpar = TRUE},
#'\code{meanspec} is omitted and before fitting a zero-mean model in the
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
#'\code{Mcf = "NP"}, \code{InfR = "Nai"},
#'\code{bStart = 0.15}, \code{bvc = "Y"}, \code{cb = 0.05}, and
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
#'\code{\link[esemifar]{tsmoothlm}} and denoted here by
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
#'@export
#'
#'@return
#'A list with the following named elements is returned.
#'\describe{
#'\item{\code{pars}:}{a named numeric vector with the parameter estimates.}
#'\item{\code{se}:}{a named numeric vector with the obtained standard errors in accordance with the parameter estimates.}
#'\item{\code{vcov_mat}:}{the variance-covariance matrix of the parameter estimates with named columns and rows.}
#'\item{\code{rt}:}{the input object \code{rt} (or at least the training data, if \code{n_test} is greater than zero); if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is kept.}
#'\item{\code{cmeans}:}{the estimated conditional means; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{cmeans}.}
#'\item{\code{sigt}:}{the estimated conditional standard deviations (or for \code{use_nonpar = TRUE} the estimated total
#'volatilities, i.e. scale function value times conditional standard deviation); if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{sigt}.}
#'\item{\code{etat}:}{the obtained residuals; if \code{rt} was a \code{"zoo"} or \code{"ts"} object, the formatting is also applied to \code{etat}.}
#'\item{\code{orders}:}{a two-element numeric vector stating the considered model orders.}
#'\item{\code{cond_dist}:}{a character value stating the conditional distribution considered in the model fitting.}
#'\item{\code{llhood}:}{the log-likelihood value obtained at the optimal parameter combination.}
#'\item{\code{inf_criteria}:}{a named two-element numeric vector with the corresponding AIC (first element) and BIC (second element) of the fitted model.}
#'\item{\code{rugarch_model}:}{the estimation object returned by \code{ugarchfit()} of the \code{rugarch} package (Ghalanos, 2024).}
#'\item{\code{meanspec}:}{the settings for the model in the conditional mean; is an object
#'of class \code{"mean_spec"} that is identical to the object passed to the input argument
#'\code{meanspec}.}
#'\item{\code{test_obs}:}{the observations at the end up the input \code{rt} reserved for
#'testing following \code{n_test}.}
#'\item{\code{scale_fun}:}{the estimated scale function values, if \code{use_nonpar = TRUE}, otherwise
#'\code{NULL}; formatting of \code{rt} is reused.}
#'\item{\code{nonpar_model}:}{the estimation object returned by \code{\link[esemifar]{tsmoothlm}} for
#'\code{use_nonpar = TRUE}.}
#'}
#'
#'@references
#'\itemize{
#'\item{Feng, Y., Gries, T., Letmathe, S., & Schulz, D. (2022). The smoots Package in R for Semiparametric Modeling of
#'Trend Stationary Time Series. The R Journal,
#'14(1), 182-195. URL: https://journal.r-project.org/articles/RJ-2022-017/.}
#'\item{Ghalanos, A. (2024). \code{rugarch}: Univariate GARCH models. R package
#'version 1.5-3. DOI: 10.32614/CRAN.package.rugarch.}
#'\item{Letmathe, S., Beran, J., & Feng, Y. (2023). An extended exponential SEMIFAR model with application
#'in R. Communications in Statistics - Theory and Methods,
#'53(22), 7914–7926. DOI: 10.1080/03610926.2023.2276049.}
#'}
#'
#'@name rugarch_wrappers
#'
#'@examples
#'est <- gjrgarch_ru(SP500)
#'est@pars
#'est@se
#'plot(est@sigt)
#'
aparch_ru <- function(rt, orders = c(1, 1),
            cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
            meanspec = mean_spec(),
            nonparspec = locpol_spec(), use_nonpar = FALSE,
            n_test = 0,
            start_pars = NULL,
            control_nonpar = list()) {

  rugarch_modelfit(
    rt = rt,
    orders = orders,
    cond_dist = cond_dist,
    model = "apARCH",
    include_mean = meanspec@include_mean, start_pars = start_pars,
    n_test = n_test,
    meanspec = meanspec,
    nonparspec = nonparspec,
    use_nonpar = use_nonpar,
    control_nonpar = control_nonpar
  )

}

#'@export
#'@rdname rugarch_wrappers

gjrgarch_ru <- function(rt, orders = c(1, 1),
            cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
            meanspec = mean_spec(),
            nonparspec = locpol_spec(), use_nonpar = FALSE,
            n_test = 0,
            start_pars = NULL,
            control_nonpar = list()) {

  rugarch_modelfit(
    rt = rt,
    orders = orders,
    cond_dist = cond_dist,
    model = "gjrGARCH",
    include_mean = meanspec@include_mean, start_pars = start_pars,
    n_test = n_test,
    meanspec = meanspec,
    nonparspec = nonparspec,
    use_nonpar = use_nonpar,
    control_nonpar = control_nonpar
  )

}

#'@export
#'@rdname rugarch_wrappers

egarch_ru <- function(rt, orders = c(1, 1),
            cond_dist = c("norm", "std", "ged", "snorm", "sstd", "sged"),
            meanspec = mean_spec(),
            nonparspec = locpol_spec(), use_nonpar = FALSE,
            n_test = 0,
            start_pars = NULL,
            control_nonpar = list()) {

  rugarch_modelfit(
    rt = rt,
    orders = orders,
    cond_dist = cond_dist,
    model = "eGARCH",
    include_mean = meanspec@include_mean, start_pars = start_pars,
    n_test = n_test,
    meanspec = meanspec,
    nonparspec = nonparspec,
    use_nonpar = use_nonpar,
    control_nonpar = control_nonpar
  )

}

#'@export
#'@rdname show-methods
#'@aliases show,fEgarch_fit_rugarch_wrapper-method
setMethod(
  "show",
  "fEGarch_fit_rugarch_wrapper",
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
    "* Fitted GARCH-Type Model (rugarch) *\n",
    "*************************************\n",
    " \n",
    paste0("GARCH-type model (rugarch): ", x@model,"\n"),
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

