
#'@rdname measure_risk
#'@export
#'

setGeneric("measure_risk", function(object, measure = c("VaR", "ES"), level = c(0.975, 0.99), ...) {standardGeneric("measure_risk")})

#'VaR and ES Computation Following Fitted Models or Forecasts
#'
#'Provides easy access to value-at-risk (VaR) and expected shortfall (ES)
#'computation for available models in this package. VaR and ES can either
#'be computed based on (a) fitted conditional means and conditional
#'standard deviations for a training period, or following (b) point forecasts
#'(either multistep or rolling) of the conditional means and conditional
#'standard deviations.
#'
#'@param object either an object of class \code{"fEGarch_fit"} returned by the
#'fitting / estimation functions of this package like returned by for example
#'\code{\link{fEGarch}} among others, an object of class \code{"fEGarch_forecast"}
#'as returned by \code{\link{predict,fEGarch_fit-method}} or
#'\code{\link{predict_roll,fEGarch_fit-method}}, or an object of class \code{"fEGarch_distr_est"} returned by the
#'distribution fitting functions of this package like returned by for example
#'\code{\link{find_dist}} among others.
#'@param measure a character vector with element \code{"VaR"}, \code{"ES"} or both;
#'indicates, what risk measures should be computed; by default, both VaR and ES
#'are computed.
#'@param level a numeric vector of arbitrary length indicating the confidence
#'levels to compute the VaR and the ES at; by default, the levels
#'\code{0.975} and \code{0.99}, i.e. 97.5 percent and 99 percent, are considered.
#'@param test_obs a series of test observations (only required when \code{object} is of class \code{"fEGarch_distr_est"}).
#'@param sigt a series of forecasted conditional standard deviations for the
#'same time points as \code{test_obs} (only required when \code{object} is of class \code{"fEGarch_distr_est"}).
#'@param cmeans a series of forecasted conditional means for the
#'same time points as \code{test_obs} (only required when \code{object} is of class \code{"fEGarch_distr_est"}).
#'@param ... currently without use.
#'
#'@details
#'Given a fitted model with fitted conditional means and conditional standard deviations
#'or given point forecasts of such series based on a fitted model, the risk measures
#'VaR and ES can be computed (at arbitrary confidence levels) following the conditional
#'loss distribution defined through the estimated / forecasted conditional mean value,
#'the estimated / forecasted conditional standard deviation value, and the assumed
#'conditional distribution (including potential estimates of distribution parameters).
#'
#'Let \eqn{\hat{\mu}_t} be the estimated / forecasted conditional mean and \eqn{\hat \sigma_t} be the
#'estimated / forecasted conditional standard deviation at some time point \eqn{t}. Furthermore,
#'define \eqn{\text{VaR}_{\eta,\alpha}} and \eqn{\text{ES}_{\eta,\alpha}} be
#'the time-invariant VaR and ES, respectively, of some identically but independently
#'distributed random variables \eqn{\eta_t} with mean zero and variance one. Given that
#'the relationship \eqn{r_t = \mu_t + \sigma_t\eta_t}, where \eqn{\mu_t} and \eqn{\sigma_t}
#'are the true conditional mean and conditional standard deviation at time \eqn{t}, is assumed
#'for some return series \eqn{\{r_t\}}, the estimated / forecasted conditional VaR and ES of \eqn{r_t} are simply
#'\deqn{\widehat{\text{VaR}}_{r,\alpha}(t) = \hat{\mu}_t + \hat{\sigma}_t \text{VaR}_{\eta,\alpha}  \hspace{3mm} \text{ and } \hspace{3mm} \widehat{\text{ES}}_{r,\alpha}(t) = \hat{\mu}_t + \hat{\sigma}_t \text{ES}_{\eta,\alpha}.}
#'This definition holds, when losses and therefore also
#'\eqn{\text{VaR}_{\eta,\alpha}(t)} and \eqn{\text{ES}_{\eta,\alpha}(t)} (for common \eqn{\alpha} such as
#'\eqn{\alpha = 0.975} or \eqn{\alpha = 0.99}) are considered
#'to be negative in sign.
#'
#'Define
#'\deqn{\text{VaR}_{\eta,\alpha} = f_{\eta}^{-1}(1-\alpha) \hspace{3mm} \text{ and } \hspace{3mm} \text{ES}_{\eta,\alpha} = (1-\alpha)^{-1}\int_{\alpha}^{1} \text{VaR}_{\eta, x} dx,}
#'which also need to be estimated for some distributions, if a distribution parameter needed to be estimated. \eqn{f} in the previous formula
#'is the cumulative distribution function of the random variables \eqn{\eta_t}. Therefore,
#'\eqn{f^{-1}_{\eta}(1-\alpha)} returns the quantile of the innovation distribution at level
#'\eqn{1-\alpha}.
#'
#'In some cases, when rolling
#'one-step forecasts of the conditional standard deviations and the conditional means
#'were obtained following a nonparametric approach, for example through
#'neural networks or similar approaches, VaR and ES are not directly to be calculated
#'because distribution assumptions have not been made. If an \code{object} that
#'is a fitted distribution to the model's standardized in-sample residuals is provided,
#'and if also test observations as well as forecasted conditional standard deviations
#'and conditional means for the test time points are passed to the method, VaR
#'and ES will be computed using the fitted distribution in \code{object}. Note
#'that \code{object} must be of class \code{"fEGarch_distr_est"}. A natural
#'selection of \code{object} is the output of \code{\link{find_dist}}, which returns
#'the best fitted model among a normal distribution, a \eqn{t}-distribution, a
#'generalized error distribution, an average Laplace distribution, and their skewed
#'variants, following either BIC (the default) or AIC. It is recommended to then
#'set \code{fix_mean = 0} and \code{fix_sdev = 1} in the call to
#'\code{\link{find_dist}} to reflect the known property that the residuals are assumed
#'to be estimated from innovations with mean zero and variance one.
#'
#'@rdname measure_risk
#'
#'@export
#'
#'@return
#'The S4 methods all return an object of class \code{"fEGarch_risk"} with elements \code{measures},
#'\code{observations} and \code{model}.
#'\code{observations} is the observed series at the time points, for which the
#'risk measures are calculated. \code{measures} is a list with elements
#'\code{VaR} and \code{ES}, distinguishing between computed VaR and ES values.
#'These elements again are list with named elements representing the various
#'computed series. \code{model} is the fitted model object.
#'
#'@examples
#'
#'# In-sample
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- fEGarch(egarch_spec(), rt)
#'risk <- measure_risk(model, measure = c("VaR", "ES"), level = c(0.95, 0.975, 0.99))
#'risk
#'
#'# Out-of-sample rolling point forecasts
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model2 <- fEGarch(egarch_spec(), rt, n_test = 250)
#'fcast <- predict_roll(model2)
#'risk2 <- measure_risk(fcast, measure = c("VaR", "ES"), level = c(0.95, 0.975, 0.99))
#'risk2
#'
#'# Use some model to obtain rolling point forecasts of
#'# the conditional mean and the conditional standard deviation for
#'# some test period; in practice, this will not be from a GARCH-type
#'# model, because it is parametric and includes a distribution assumption,
#'# but instead from some nonparametric model
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2005-12-31")
#'model <- fEGarch(egarch_spec(), rt, n_test = 250)
#'fc <- model %>% predict_roll()
#'
#'test_obs <- model@test_obs   # Test observations
#'sigt <- fc@sigt              # Conditional volatility forecasts
#'cmeans <- fc@cmeans          # Conditional mean forecasts
#'
#'resids <- model@etat         # In-sample standardized residuals
#'
#'# Given 'test_obs', 'sigt', 'cmeans' and 'resids', we can now
#'# compute the VaR and ES forecasts for the test period
#'
#'dist <- find_dist(resids, fix_mean = 0, fix_sdev = 1)
#'dist
#'
#'risk <- dist %>%
#'  measure_risk(test_obs = test_obs, sigt = sigt, cmeans = cmeans)
#'
#'plot(risk, which = 0.975)
#'
#'
#'
setMethod("measure_risk", "fEGarch_fit",
  function(object, measure = c("VaR", "ES"), level = c(0.975, 0.99), ...) {
    cond_d <- object@cond_dist

    par_name <- switch(
      cond_d,
      "norm" = "df",   # irrelevant
      "snorm" = "df",  # irrelevant
      "std" = "df",
      "sstd" = "df",
      "ged" = "shape",
      "sged" = "shape",
      "ald" = "P",
      "sald" = "P"
    )
    shape <- if (!(cond_d %in% c("norm", "snorm"))) {
      object@pars[[par_name]]
    } else {
      0
    }
    skew <- if (!(cond_d %in% c("norm", "std", "ged", "ald"))) {
      object@pars[["skew"]]
    } else {
      0
    }


    args <- list(
      level = level,
      dist = cond_d,
      skew = skew
    )
    args[[par_name]] <- shape

    if ("VaR" %in% measure) {
      q_VaR <- do.call(what = VaR_calc, args = args)
      names_q <- paste0("VaR", level)
      lVaR <- list()
      for (i in seq_along(q_VaR)) {
        lVaR[[names_q[[i]]]] <- object@cmeans + q_VaR[[i]] * object@sigt
      }
    } else {
      lVaR <- list()
    }

    if ("ES" %in% measure) {
      q_ES <- do.call(what = ES_calc, args = args)
      names_q <- paste0("ES", level)
      lES <- list()
      for (i in seq_along(q_ES)) {
        lES[[names_q[[i]]]] <- object@cmeans + q_ES[[i]] * object@sigt
      }
    } else {
      lES <- list()
    }

    fEGarch_risk(
      measures = list("VaR" = lVaR, "ES" = lES),
      observations = object@rt,
      model = object,
      sigt = object@sigt,
      cmeans = object@cmeans
    )

  }
)

#'@export
#'@rdname measure_risk
setMethod("measure_risk", "fEGarch_forecast",
  function(object, measure = c("VaR", "ES"), level = c(0.975, 0.99), ...) {

    object2 <- object@model
    object2@sigt <- object@sigt
    object2@cmeans <- object@cmeans

    func <- methods::selectMethod("measure_risk", "fEGarch_fit")
    out <- func(object = object2, measure = measure, level = level, ...)
    out@observations <- object2@test_obs
    out@sigt <- object2@sigt
    out@cmeans <- object2@cmeans
    out
  }
)

#'@rdname backtest-generics
#'@export
setGeneric("trafflight_test", function(object, ...) {standardGeneric("trafflight_test")})

zone_selector <- function(p) {
  out <- if (p < 0.95) {
    cli::col_green("Green zone")
  } else if (p >= 0.95 && p < 0.9999) {
    cli::col_yellow("Yellow zone")
  } else if (p >= 0.9999) {
    cli::col_red("Red zone")
  }
  out
}

#'Backtesting VaR and ES
#'
#'Run traffic light tests for value at risk (VaR) and expected
#'shortfall (ES) as well as a selection of coverage and independence
#'tests for VaR.
#'
#'@param object an object of class \code{"fEGarch_risk"}.
#'@param silent a logical value indicating whether or not to
#'print test results in a nicely formatted manner to the console.
#'@param ... currently without use.
#'
#'@details
#'
#'\code{backtest_suite} runs all the other backtesting methods.
#'\code{cov_tests} runs all of \code{uncond_cov_test},
#'\code{indep_test} and \code{cond_cov_test}.
#'
#'Traffic light tests (\code{trafflight_test}):
#'
#'Given an input object \code{object} of class \code{"fEGarch_risk"},
#'traffic light tests for value at risk (VaR) and expected shortfall (ES)
#'are applied to the individual risk measure series in the object. Note
#'that in order for a traffic light test in context of ES being applicable,
#'the corresponding VaR series of the same confidence level must also
#'be present in \code{object}. If this is not fulfilled,
#'messages will be printed to the console, making the user aware of
#'this issue.
#'
#'Let the number of test observations be denoted by
#'\eqn{n\in \mathbb{N}} and let \eqn{\{r_t\}}, \eqn{t=1,\dots,n},
#'be the test returns. \eqn{\{\text{VaR}_t\}} are the (one-step rolling) VaR point forecasts
#'for the same period with confidence level \eqn{\alpha}. Denote by
#'\eqn{I_t} an indicator that equals 1, whenever \eqn{r_t < \text{VaR}_t}, and
#'0 otherwise, and define \eqn{K_1 = \sum_{t=1}^{n}I_t}. \eqn{I_t} are assumed
#'to follow a binomial distribution with probability \eqn{P = \alpha} for any \eqn{I_t = 0}.
#'Then \eqn{C} is computed as the cumulative probability of observing \eqn{K_1}
#'under \eqn{P}. The forecasted VaR series is then classified following \eqn{C}. If
#'\eqn{C < 0.95}, then it is sorted into the green zone, if \eqn{0.95 \leq C < 0.9999},
#'then the series belongs to the yellow zone, and if \eqn{C \geq 0.9999},
#'then the class of the VaR series is the red zone (Basel Committee on Banking Supervision, 1996).
#'
#'The traffic light test for the ES (Costanzino and Curran, 2018) uses a similar classification system
#'based on the severity of breaches
#'\deqn{B = \sum_{t = 1}^{n} \frac{1-F(\hat \eta_t)-\alpha}{1-\alpha}I_t,}
#'where \eqn{F} is the (fitted) cumulative distribution function of the
#'standardized innovations and with \eqn{\hat \eta_t} as the standardized residuals
#'of a fitted GARCH-type model (or of its semiparametric extension).
#'Then \eqn{B \overset{a}{\sim}N(\mu_{\text{ES}}, \sigma^2_{\text{ES}})} with
#'\eqn{\mu_{\text{ES}} = 0.5(1-\alpha)n} and
#'\eqn{\sigma^2_{\text{ES}} = (1-\alpha)[(1+3\alpha) / 12]}. The cumulative
#'probability of observing a severity of breaches of \eqn{B} or less can
#'be computed and classified in the same way as for the VaR traffic light test
#'using this asymptotic distribution.
#'
#'Weighted Absolute Deviation (WAD) (\code{WAD}):
#'Following the standard computation of the 99%-VaR, the 97.5%-VaR and the
#'97.5%-ES for the traffic light tests, the WAD criterion takes all of these
#'into account and summarizes them into one numeric value. Let \eqn{N_1} be the
#'observed breaches for the 99%-VaR for the test set and let \eqn{\mu_1} be the
#'corresponding expected number of breaches. \eqn{N_2} and \eqn{\mu_2} are to
#'understood analogously for the 97.5%-VaR. \eqn{N_3} is then the severity of
#'breaches of the 97.5%-ES (i.e. it is equal to \eqn{B} from before) and \eqn{\mu_3}
#'is \eqn{\mu_{\text{ES}}} from before. Then
#'\deqn{\text{WAD} = \frac{|N_1-\mu_1|}{\mu_1} + \frac{|N_2-\mu_2|}{\mu_2} + \frac{|N_3-\mu_3|}{\mu_3}.}
#'See also Letmathe et al. (2022) for further information.
#'
#'Coverage and independence tests (\code{cov_tests}):
#'
#'Following Christoffersen (1998), the backtesting suite also includes
#'a selection of coverage and independence tests regarding the VaR. Let the number of test
#'observations be denoted by \eqn{n\in \mathbb{N}} and let \eqn{\{r_t\}}, \eqn{t=1,\dots,n},
#'be the test returns. \eqn{\{\text{VaR}_t\}} are the (one-step rolling) VaR point forecasts
#'for the same period with confidence level \eqn{\alpha}. Furthermore, define
#'\eqn{I_t} to be an indicator that equals \eqn{1}, whenever \eqn{r_t < \text{VaR}_t} and
#'zero otherwise. Let \eqn{K_1 = \sum_{t=1}^{n}I_t} and \eqn{K_0 = n - K_1}. Furthermore,
#'\eqn{\hat z_1 = K_1 / (K_0 + K_1)} and \eqn{\hat z_0 = K_0 / (K_0 + K_1)} as well as
#'\deqn{L_{\hat z} = \hat z_0^{K_0} \hat z_1^{K_1}}
#'and
#'\deqn{L_{\alpha} = \alpha^{K_0}(1-\alpha)^{K_1}.}
#'
#'In addition, we require \eqn{I^{*}_{i,j}(t)}, \eqn{t = 2,\dots,n} and \eqn{i,j \in \{0,1\}},
#'to be other indicators that equal 1, whenever \eqn{I_t=j} and
#'simultaneously \eqn{I_{t-1} = i}. Per consequence,
#'\eqn{K_{i,j}=\sum_{t=2}^{n} I^{*}_{i,j}(t)} and \eqn{\hat z_{i,j} = K_{i,j} / (K_{i,0} + K_{i, 1})}.
#'Moreover, \eqn{\hat z_1^{*} = (K_{0,1}+K_{1,1}) / (n - 1)} and
#'\eqn{\hat z_0^{*} = 1-\hat z_1^{*}}. Now,
#'\deqn{L_{\hat z_{0,0}} = \hat z_{0,0}^{K_{0,0}} \hat z_{0,1}^{K_{0,1}} \hat z_{1,0}^{K_{1,0}} \hat z_{1,1}^{K_{1,1}}}
#'and
#'\deqn{L_{\hat z^{*}} = (\hat z_{0}^{*})^{(K_{0,0} + K_{1,0})} (\hat z_{1}^{*})^{(K_{0,1} + K_{1,1})}.}
#'
#'Ultimately,
#'\deqn{L_{\alpha^{*}} = \alpha^{(K_{0,0} + K_{1,0})}(1-\alpha)^{(K_{0,1} + K_{1, 1})}.}
#'
#'The three test statistics following Christoffersen (1998) are then
#'\deqn{S_{\text{uc}} = -2 \ln\left[L_{\alpha} / L_{\hat{z}}\right] \overset{a}{\sim} \chi^2 (1),}
#'\deqn{S_{\text{ind}} = -2 \ln\left[L_{\hat z^{*}} / L_{\hat{z}_{0,0}}\right] \overset{a}{\sim} \chi^2 (1), \hspace{4mm} \text{and}}
#'\deqn{S_{\text{cc}} = -2 \ln\left[L_{\alpha^{*}} / L_{\hat{z}_{0,0}}\right] \overset{a}{\sim} \chi^2 (2),}
#'where \eqn{S_{\text{uc}}} is the test statistic of the unconditional coverage test,
#'\eqn{S_{\text{ind}}} is that of the independence test and
#'\eqn{S_{\text{cc}}} is that of the conditional coverage test.
#'
#'@rdname backtest-tests
#'
#'@return
#'All methods return a list invisibly. The elements of the list differ
#'slightly depending on the method. Moreover, for \code{silent = FALSE},
#'the default, test results are printed to the console.
#'
#'@references
#'\itemize{
#'\item{Basel Committee on Banking Supervision (1996). Supervisory Framework For The Use of "Backtesting" in
#'Conjunction With The Internal Models Approach to Market Risk Capital Requirements.
#'URL: https://www.bis.org/publ/bcbs22.pdf.}
#'\item{Christoffersen, P. F. (1998). Evaluating Interval Forecasts. International Economic Review,
#'39(4): 841-862. DOI: 10.2307/2527341.}
#'\item{Costanzino, N., & Curran, M. (2018). A Simple Traffic Light Approach to Backtesting Expected Shortfall.
#'Risks, 6(1). DOI: 10.3390/risks6010002.}
#'\item{Letmathe, S., Feng, Y., & Uhde, A. (2022). Semiparametric GARCH models with long memory
#'applied to Value at Risk and Expected Shortfall. Journal of Risk,
#'25(2). DOI: 10.21314/JOR.2022.044.}
#'}
#'
#'
#'@export
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- fEGarch(egarch_spec(), rt, n_test = 250)
#'fcast <- predict_roll(model)
#'risk <- measure_risk(fcast, measure = c("VaR", "ES"), level = c(0.95, 0.975, 0.99))
#'trafflight_test(risk)
#'cov_tests(risk)
#'backtest_suite(risk)
#'
setMethod("trafflight_test", "fEGarch_risk",
  function(object, silent = FALSE, ...) {

  obs <- object@observations
  n <- length(obs)

  slotnames <- methods::slotNames(object@model)
  cond_d <- if ("cond_dist" %in% slotnames) {
    object@model@cond_dist
  } else if ("dist" %in% slotnames) {
    object@model@dist
  } else {
    stop("Neither slot @cond_dist nor @dist to be found in model.")
  }

  par_name <- switch(
    cond_d,
    "norm" = "df",    # irrelevant
    "snorm" = "df",   # irrelevant
    "std" = "df",
    "sstd" = "df",
    "ged" = "shape",
    "sged" = "shape",
    "ald" = "P",
    "sald" = "P"
  )
  shape <- if (cond_d %in% c("norm", "snorm")) {
    0
  } else {
    object@model@pars[[par_name]]
  }

  skew <- if (cond_d %in% c("snorm", "sstd", "sged", "sald")) {
    object@model@pars[["skew"]]
  } else {
    0
  }

  VaR_names <- names(object@measures$VaR)
  ES_names <- names(object@measures$ES)

  n_VaR <- length(VaR_names)
  n_ES <- length(ES_names)

  list_VaR <- vector(mode = "list", length = n_VaR)
  list_ES <- vector(mode = "list", length = n_ES)

  if (n_VaR > 0) {

    names(list_VaR) <- VaR_names
    for (i in 1:n_VaR) {
      conf_lvl <- as.numeric(substring(VaR_names[[i]], 4))
      pot <- sum(obs < object@measures$VaR[[i]])
      prob <- stats::pbinom(q = pot, size = n, prob = 1 - conf_lvl)
      list_VaR[[i]] <- list(
        conf_lvl = conf_lvl,
        pot = pot,
        prob = prob,
        zone = zone_selector(prob)
      )
    }

  } else {
    message("Note: No VaR measures to analyze in input object")
  }

  if (n_ES > 0 && n_VaR > 0) {

    cmeans <- object@cmeans
    sigt <- object@sigt

    et <- (obs - cmeans) / sigt     # standardized residuals
    args <- list(
      dfun = fun1_selector(cond_d),
      skew = skew,
      shape = shape
    )


    names(list_ES) <- ES_names
    for (i in 1:n_ES) {
      conf_lvl <- as.numeric(substring(ES_names[[i]], 3))
      check_VaR <- paste0("VaR", conf_lvl)
      if (!(check_VaR %in% VaR_names)) {
        message(paste0(ES_names, "cannot be analyzed because ", check_VaR, "is missing from input object"))
        next
      }

      pot_idx <- obs < object@measures$VaR[[check_VaR]]
      alpha <- conf_lvl
      sum_severity <- if (sum(pot_idx) > 0) {
        et_sub <- zoo::coredata(et[pot_idx])
        args[["x"]] <- et_sub
        nprobs <- do.call(what = cdf_fun, args = args)
        probs <- 1 - nprobs

        sum((probs - alpha) / (1 - alpha))
      } else {
        0
      }
        mean_ES <- 0.5 * (1 - alpha) * n
        var_ES <- n * (1 - alpha) * ((1 + 3 * alpha) / 12)

        prob <- pnorm(q = sum_severity, mean = mean_ES, sd = sqrt(var_ES))

      list_ES[[i]] <- list(
        conf_lvl = conf_lvl,
        severity = sum_severity,
        prob = prob,
        zone = zone_selector(prob)
      )
    }

  } else {
    message("Note: No VaR and / or ES measures to analyze in input object")
  }

  if (!silent) {
    msg_VaR <- paste0(
      "\nVaR results:\n",
      "************\n\n",
      paste0(vapply(
        X = list_VaR,
        FUN = function(.x) {
          paste0(
            "Conf. level: ", .x$conf_lvl, "\n",
            "Breaches: ", .x$pot, "\n",
            "Cumul. prob.: ", sprintf("%.4f", .x$prob), "\n",
            "Zone: ", .x$zone, "\n\n"
          )
        }, FUN.VALUE = character(1)
      ), collapse = "")
    )

    msg_ES <- paste0(
      "\nES results:\n",
      "***********\n\n",
      paste0(vapply(
        X = list_ES,
        FUN = function(.x) {
          paste0(
            "Conf. level: ", .x$conf_lvl, "\n",
            "Severity of breaches: ", sprintf("%.4f", .x$severity), "\n",
            "Cumul. prob.: ", sprintf("%.4f", .x$prob), "\n",
            "Zone: ", .x$zone, "\n\n"
          )
        }, FUN.VALUE = character(1)
      ), collapse = "")
    )

    cli::cat_line(paste0(
      "\n***********************\n",
      "* Traffic light tests *\n",
      "***********************\n",
      msg_VaR,
      msg_ES))

  }

  invisible(list(VaR = list_VaR, ES = list_ES))

})

#'@rdname backtest-generics
#'@export
setGeneric("uncond_cov_test", function(object, ...) {standardGeneric("uncond_cov_test")})

rejecter <- function(pval) {
  if (pval < 0.05) {
    cli::col_red("Reject H0")
  } else if (pval >= 0.05) {
    cli::col_green("Do not reject H0")
  }
}

#'@rdname backtest-tests
#'@export
setMethod("uncond_cov_test", "fEGarch_risk",
  function(object, silent = FALSE, ...) {

    VaR <- object@measures$VaR

    VaR_names <- names(VaR)
    conf_lvls <- as.numeric(substring(VaR_names, 4))

    m <- length(conf_lvls)

    list_out <- vector(mode = "list", length = m)
    names(list_out) <- VaR_names
    obs <- zoo::coredata(object@observations)
    n <- length(obs)

    for (i in 1:m) {
      alpha <- conf_lvls[[i]]   # Covered theoretically
      K1 <- sum(obs < zoo::coredata(VaR[[i]]))
      K0 <- n - K1

      z0 <- K0 / n

      stat <- -2 * log(alpha^K0 * (1 - alpha)^K1 / (z0^K0 * (1 - z0)^K1))
      pval <- 1 - stats::pchisq(stat, df = 1)

      list_out[[i]] <- list(
        conf_lvl = alpha,
        breaches = K1,
        test_statistic = stat,
        p_value = pval,
        decision = rejecter(pval)
      )

    }

    if (!silent) {

      header <- paste0(
        "\n***************************************\n",
        "* Unconditional coverage test for VaR *\n",
        "***************************************\n\n",
        "H0: true share of covered observations = theoretical share of VaR\n"
      )

      tests <- paste0(vapply(
        list_out, FUN = function(.x) {
          paste0(
            "\nConf. level: ", .x$conf_lvl, "\n",
            "Breaches: ", .x$breaches, "\n",
            "Test statistic: ", sprintf("%.4f", .x$test_statistic), "\n",
            "p-value: ", sprintf("%.4f", .x$p_value), "\n",
            "Decision: ", .x$decision, "\n"
          )
        }, FUN.VALUE = character(1)
      ), collapse = "")

      cli::cat_line(paste0(
        header, tests
      ))

    }

    invisible(list_out)

})

#'@rdname backtest-generics
#'@export
setGeneric("indep_test", function(object, ...) {standardGeneric("indep_test")})

#'@rdname backtest-tests
#'@export
setMethod("indep_test", "fEGarch_risk",
  function(object, silent = FALSE, ...) {

    VaR <- object@measures$VaR

    VaR_names <- names(VaR)
    conf_lvls <- as.numeric(substring(VaR_names, 4))

    m <- length(conf_lvls)

    list_out <- vector(mode = "list", length = m)
    names(list_out) <- VaR_names
    obs <- zoo::coredata(object@observations)
    n <- length(obs)
    n_adj <- n - 1

    for (i in 1:m) {
      alpha <- conf_lvls[[i]]   # Covered theoretically
      It <- obs < zoo::coredata(VaR[[i]])

      pre1 <- utils::head(It, -1) == 1
      pre0 <- utils::head(It, -1) == 0
      post1 <- utils::tail(It, n_adj) == 1
      post0 <- utils::tail(It, n_adj) == 0

      K00 <- sum(pre0 & post0)
      K10 <- sum(pre1 & post0)
      K01 <- sum(pre0 & post1)
      K11 <- sum(pre1 & post1)

      K0 <- K00 + K01
      K1 <- K10 + K11

      z00 <- K00 / K0
      z10 <- K10 / K1

      L_z00 <- z00^K00 * (1 - z00)^K01 * z10^K10 * (1 - z10)^K11

      K <- K0 + K1
      Ki0 <- K00 + K10
      Ki1 <- K11 + K01
      z0 <- Ki0 / K
      z1 <- Ki1 / K

      L_z0 <- z0^Ki0 * z1^Ki1

      stat <- -2 * log(L_z0 / L_z00)

      pval <- 1 - stats::pchisq(stat, df = 1)

      list_out[[i]] <- list(
        conf_lvl = alpha,
        breaches = K1,
        test_statistic = stat,
        p_value = pval,
        decision = rejecter(pval)
      )

    }

    if (!silent) {

      header <- paste0(
        "\n*****************************\n",
        "* Independence test for VaR *\n",
        "*****************************\n\n",
        "H0: true share of covered observations independent\n    of breach or no breach at previous time point\n"
      )

      tests <- paste0(vapply(
        list_out, FUN = function(.x) {
          paste0(
            "\nConf. level: ", .x$conf_lvl, "\n",
            "Breaches: ", .x$breaches, "\n",
            "Test statistic: ", sprintf("%.4f", .x$test_statistic), "\n",
            "p-value: ", sprintf("%.4f", .x$p_value), "\n",
            "Decision: ", .x$decision, "\n"
          )
        }, FUN.VALUE = character(1)
      ), collapse = "")

      cli::cat_line(paste0(
        header, tests
      ))

    }

    invisible(list_out)

})

#'Generics for backtests
#'
#'Generic functions to build backtesting methods from.
#'
#'@param object the generics are currently without use.
#'@param ... the generics are currently without use.
#'
#'@details
#'The generics are currently without use.
#'
#'@export
#'
#'@return
#'The generics are currently without use, can therefore not be called
#'and thus don't produce results.
#'
#'@rdname backtest-generics
#'

setGeneric("cond_cov_test", function(object, ...) {standardGeneric("cond_cov_test")})

#'@rdname backtest-tests
#'@export
setMethod("cond_cov_test", "fEGarch_risk",
  function(object, silent = FALSE, ...) {

    VaR <- object@measures$VaR

    VaR_names <- names(VaR)
    conf_lvls <- as.numeric(substring(VaR_names, 4))

    m <- length(conf_lvls)

    list_out <- vector(mode = "list", length = m)
    names(list_out) <- VaR_names
    obs <- zoo::coredata(object@observations)
    n <- length(obs)
    n_adj <- n - 1

    for (i in 1:m) {
      alpha <- conf_lvls[[i]]   # Covered theoretically
      It <- obs < zoo::coredata(VaR[[i]])
      K1 <- sum(It[-1])
      K0 <- n_adj - K1

      L_a <- alpha^K0 * (1 - alpha)^K1

      pre1 <- utils::head(It, -1) == 1
      pre0 <- utils::head(It, -1) == 0
      post1 <- utils::tail(It, n_adj) == 1
      post0 <- utils::tail(It, n_adj) == 0

      K00 <- sum(pre0 & post0)
      K10 <- sum(pre1 & post0)
      K01 <- sum(pre0 & post1)
      K11 <- sum(pre1 & post1)

      z00 <- K00 / (K00 + K01)
      z10 <- K10 / (K11 + K10)

      L_z00 <- z00^K00 * (1 - z00)^K01 * z10^K10 * (1 - z10)^K11

      stat <- -2 * log(L_a / L_z00)

      pval <- 1 - stats::pchisq(stat, df = 2)

      list_out[[i]] <- list(
        conf_lvl = alpha,
        breaches = K1,
        test_statistic = stat,
        p_value = pval,
        decision = rejecter(pval)
      )

    }

    if (!silent) {

      header <- paste0(
        "\n*************************************\n",
        "* Conditional coverage test for VaR *\n",
        "*************************************\n\n",
        "H0: true share of covered observations simultaneously\n    independent of breach or no breach at previous time point\n    and equal to theoretical share of VaR\n"
      )

      tests <- paste0(vapply(
        list_out, FUN = function(.x) {
          paste0(
            "\nConf. level: ", .x$conf_lvl, "\n",
            "Breaches: ", .x$breaches, "\n",
            "Test statistic: ", sprintf("%.4f", .x$test_statistic), "\n",
            "p-value: ", sprintf("%.4f", .x$p_value), "\n",
            "Decision: ", .x$decision, "\n"
          )
        }, FUN.VALUE = character(1)
      ), collapse = "")

      cli::cat_line(paste0(
        header, tests
      ))

    }

    invisible(list_out)

})



#'@rdname backtest-generics
#'@export
setGeneric("cov_tests", function(object, ...) {standardGeneric("cov_tests")})

#'@rdname backtest-tests
#'@export
setMethod("cov_tests", "fEGarch_risk",
  function(object, silent = FALSE, ...) {

    l1 <- uncond_cov_test(object, silent = silent, ...)
    l2 <- indep_test(object, silent = silent, ...)
    l3 <- cond_cov_test(object, silent = silent, ...)

    list_out <- list(
      uncond_cov_test = l1,
      indep_test = l2,
      cond_cov_test = l3
    )

    invisible(list_out)

})

#'@rdname backtest-generics
#'@export
setGeneric("backtest_suite", function(object, ...) {standardGeneric("backtest_suite")})

#'@rdname backtest-tests
#'@export
setMethod("backtest_suite", "fEGarch_risk",
  function(object, silent = FALSE, ...) {

    l0 <- trafflight_test(object, silent = silent, ...)
    l1 <- WAD(object, silent = silent, ...)
    l2 <- uncond_cov_test(object, silent = silent, ...)
    l3 <- indep_test(object, silent = silent, ...)
    l4 <- cond_cov_test(object, silent = silent, ...)

    list_out <- list(
      trafflight_test = l0,
      wad = l1,
      uncond_cov_test = l2,
      indep_test = l3,
      cond_cov_test = l4
    )

    invisible(list_out)

})


########### loss functions ###############

loss0 <- function(object, penalty = 1e-4, ...) {

  obs <- zoo::coredata(object@observations)
  n <- length(obs)

  VaR_names <- names(object@measures$VaR)
  m_VaR <- length(VaR_names)
  ES_names <- names(object@measures$ES)
  m_ES <- length(ES_names)

  list_out <- list(
    VaR = vector(mode = "list", length = m_VaR),
    ES = vector(mode = "list", length = m_ES)
  )
  names(list_out$VaR) <- VaR_names
  names(list_out$ES) <- ES_names

  if (m_VaR > 0) {
    for (i in 1:m_VaR) {
      VaR <- zoo::coredata(object@measures$VaR[[i]])
      idx <- obs < VaR

      l1 <- rep(NA, n)
      l1[idx] <- (VaR[idx] - obs[idx])^2
      l1[!idx] <- 0

      list_out$VaR[[i]] <- sum(l1)
    }
  }

  if (m_ES > 0) {
    for (i in 1:m_ES) {
      ES <- zoo::coredata(object@measures$ES[[i]])
      idx <- obs < ES

      l1 <- rep(NA, n)
      l1[idx] <- (ES[idx] - obs[idx])^2
      l1[!idx] <- 0

      list_out$ES[[i]] <- sum(l1)
    }
  }

  list_out

}

loss1 <- function(object, penalty = 1e-4, ...) {

  obs <- zoo::coredata(object@observations)
  n <- length(obs)

  VaR_names <- names(object@measures$VaR)
  m_VaR <- length(VaR_names)
  ES_names <- names(object@measures$ES)
  m_ES <- length(ES_names)

  list_out <- list(
    VaR = vector(mode = "list", length = m_VaR),
    ES = vector(mode = "list", length = m_ES)
  )
  names(list_out$VaR) <- VaR_names
  names(list_out$ES) <- ES_names

  if (m_VaR > 0) {
    for (i in 1:m_VaR) {
      VaR <- zoo::coredata(object@measures$VaR[[i]])
      idx <- obs < VaR

      l1 <- rep(NA, n)
      l1[idx] <- (VaR[idx] - obs[idx])^2
      l1[!idx] <- penalty * abs(VaR[!idx])

      list_out$VaR[[i]] <- sum(l1)
    }
  }

  if (m_ES > 0) {
    for (i in 1:m_ES) {
      ES <- zoo::coredata(object@measures$ES[[i]])
      idx <- obs < ES

      l1 <- rep(NA, n)
      l1[idx] <- (ES[idx] - obs[idx])^2
      l1[!idx] <- penalty * abs(ES[!idx])

      list_out$ES[[i]] <- sum(l1)
    }
  }

  list_out

}

loss2 <- function(object, penalty = 1e-4, ...) {

  obs <- zoo::coredata(object@observations)
  n <- length(obs)

  VaR_names <- names(object@measures$VaR)
  m_VaR <- length(VaR_names)
  ES_names <- names(object@measures$ES)
  m_ES <- length(ES_names)

  list_out <- list(
    VaR = vector(mode = "list", length = m_VaR),
    ES = vector(mode = "list", length = m_ES)
  )
  names(list_out$VaR) <- VaR_names
  names(list_out$ES) <- ES_names

  if (m_VaR > 0) {
    for (i in 1:m_VaR) {
      VaR <- zoo::coredata(object@measures$VaR[[i]])
      idx <- obs < VaR

      l1 <- rep(NA, n)
      l1[idx] <- (VaR[idx] - obs[idx])^2
      l1[!idx] <- penalty * (abs(VaR[!idx] - obs[!idx]))

      list_out$VaR[[i]] <- sum(l1)
    }
  }

  if (m_ES > 0) {
    for (i in 1:m_ES) {
      ES <- zoo::coredata(object@measures$ES[[i]])
      idx <- obs < ES

      l1 <- rep(NA, n)
      l1[idx] <- (ES[idx] - obs[idx])^2
      l1[!idx] <- penalty * (abs(ES[!idx] - obs[!idx]))

      list_out$ES[[i]] <- sum(l1)
    }
  }

  list_out

}

loss3 <- function(object, penalty = 1e-4, ...) {

  obs <- zoo::coredata(object@observations)
  n <- length(obs)

  VaR_names <- names(object@measures$VaR)
  m_VaR <- length(VaR_names)
  ES_names <- names(object@measures$ES)
  m_ES <- length(ES_names)

  list_out <- list(
    VaR = vector(mode = "list", length = m_VaR),
    ES = vector(mode = "list", length = m_ES)
  )
  names(list_out$VaR) <- VaR_names
  names(list_out$ES) <- ES_names

  if (m_VaR > 0) {
    for (i in 1:m_VaR) {
      VaR <- zoo::coredata(object@measures$VaR[[i]])
      idx <- obs < VaR

      l1 <- rep(NA, n)
      l1[idx] <- (VaR[idx] - obs[idx])^2
      diffe <- abs(VaR[!idx] - obs[!idx])
      VaR_a <- abs(VaR[!idx])
      selector <- diffe > VaR_a
      diffe[selector] <- VaR_a[selector]
      l1[!idx] <- penalty * diffe

      list_out$VaR[[i]] <- sum(l1)
    }
  }

  if (m_ES > 0) {
    for (i in 1:m_ES) {
      ES <- zoo::coredata(object@measures$ES[[i]])
      idx <- obs < ES

      l1 <- rep(NA, n)
      l1[idx] <- (ES[idx] - obs[idx])^2
      diffe <- abs(ES[!idx] - obs[!idx])
      ES_a <- abs(ES[!idx])
      selector <- diffe > ES_a
      diffe[selector] <- ES_a[selector]
      l1[!idx] <- penalty * diffe

      list_out$ES[[i]] <- sum(l1)
    }
  }

  list_out

}

#'Generic for Loss Function Calculation
#'
#'Currently without use. Use the methods derived from
#'this generic.
#'
#'@param object currently without use.
#'@param penalty currently without use.
#'@param ... currently without use.
#'
#'@return
#'The generic itself is currently without use and thus
#'does not return anything.
#'
#'@export
#'
setGeneric("loss_functions", function(object, penalty = 1e-4, ...) {standardGeneric("loss_functions")})

#'Loss Function Calculation
#'
#'Compute loss function values given log-returns and corresponding
#'value at risk (VaR) and expected shortfall (ES) series.
#'
#'@param object an object of class \code{"fEGarch_risk"} as returned
#'by the model fitting functions of this package, for example
#'\code{\link{fEGarch}}.
#'@param penalty the penalty term to use in the opportunity cost terms.
#'@param ... currently without use.
#'
#'@details
#'Let \eqn{n \in \mathbb{N}} be the number of observations of a (log-)return
#'series \eqn{\{r_t\}}, \eqn{t=1,\dots,n}, and let \eqn{\text{VaR}_t} and
#'\eqn{\text{ES}_t} be the
#'estimated or forecasted VaR and ES (at some confidence level \eqn{\alpha}) at time
#'\eqn{t}, respectively. Such series are included in an object of class
#'\code{"fEGarch_risk"}. In the following, a risk measure at time \eqn{t} is
#'simply denoted by \eqn{\text{RM}_t} and can either mean
#'\eqn{\text{VaR}_t} or
#'\eqn{\text{ES}_t}.
#'
#'Based on a calculated VaR and / or expected shortfall (ES), capital needs
#'to be held back following regulatory rules. Commonly, among many models
#'used for forecasting risk measures that fulfill regulatory conditions,
#'loss functions are computed that also consider opportunity costs in to
#'assess, what model that fulfills regulatory rules minimizes such loss
#'functions. Let \eqn{\Omega \geq 0} be the penalty term.
#'
#'For all loss functions we have
#'\deqn{\text{LF}_i = \sum_{t=1}^{n} l_{t,i}, \hspace{3mm} i = 0,1,2,3,}
#'as the loss function with
#'\deqn{l_{t,i} = (\text{RM}_t - r_t)^2, \hspace{3mm} i = 0,1,2,3,}
#'for \eqn{r_t < \text{RM}_t}. They differ in how the case
#'\eqn{r_t \geq \text{RM}_t} is treated.
#'
#'The regulatory loss function (\code{rlf}) uses \eqn{l_{t,0} = 0}.
#'
#'The firm's loss function (Sarma et al., 2003) (\code{flf}) considers
#'\eqn{l_{t,1}=\Omega |\text{RM}_t|}.
#'
#'The adjusted loss function (Abad et al., 2015) (\code{alf}) makes use
#'of \eqn{l_{t,2} = \Omega |\text{RM}_t - r_t|}.
#'
#'The corrected loss function (Feng, forthcoming) (\code{clf}) has
#'\eqn{l_{t,3} = \Omega \text{min}\left(|\text{RM}_t - r_t|, |\text{RM}_t|\right)}.
#'
#'@return
#'Returns a list with the four elements \code{rlf}, \code{flf},
#'\code{alf} and \code{clf}, each lists with numeric vector
#'elements \code{VaR} and \code{ES}. The four elements correspond
#'to the regulatory loss function, the firm's loss function,
#'the adjusted loss function and the corrected loss function.
#'
#'@export
#'
#'@references
#'\itemize{
#'\item{Abad, P., Muela, S. B., & Martín, C. L. (2015). The role of
#'the loss function in value-at-risk comparisons. The Journal of Risk
#'Model Validation, 9(1): 1-19. DOI: 10.21314/JRMV.2015.132.}
#'\item{Sarma, M., Thomas, S., & Shah, A. (2003). Selection of Value-at-Risk
#'models. Journal of Forecasting,
#'22(4): 337-358. DOI: 10.1002/for.868.}
#'}
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- fEGarch(egarch_spec(), rt, n_test = 250)
#'fcast <- predict_roll(model)
#'risk <- measure_risk(fcast, measure = c("VaR", "ES"), level = c(0.95, 0.975, 0.99))
#'loss_functions(risk)
#'
setMethod("loss_functions", "fEGarch_risk",
  function(object, penalty = 1e-4, ...) {


    list_out <- list(
      rlf = loss0(object, penalty = penalty, ...),
      flf = loss1(object, penalty = penalty, ...),
      alf = loss2(object, penalty = penalty, ...),
      clf = loss3(object, penalty = penalty, ...)
    )

    list_out

})

#'@rdname backtest-generics
#'@export
setGeneric("WAD", function(object, ...) {standardGeneric("WAD")})

#'@rdname backtest-tests
#'@export
setMethod("WAD", "fEGarch_risk",
  function(object, silent = FALSE, ...) {

  VaR_names <- names(object@measures$VaR)
  ES_names <- names(object@measures$ES)

  if ("VaR0.975" %in% VaR_names && "VaR0.99" %in% VaR_names && "ES0.975" %in% ES_names) {

    results <- trafflight_test(object, silent = TRUE)

    N3 <- results$ES$ES0.975$severity    # severity of breaches for 97.5%-ES
    N2 <- results$VaR$VaR0.975$pot       # number of breaches for 97.5%-VaR
    N1 <- results$VaR$VaR0.99$pot        # number of breaches for 99%-VaR

    n <- length(object@observations)

    mu1 <- n * 0.01
    mu2 <- n * 0.025
    mu3 <- 0.5 * 0.025 * n

    crit1 <- abs(N1 - mu1) / mu1
    crit2 <- abs(N2 - mu2) / mu2
    crit3 <- abs(N3 - mu3) / mu3

    crit_total <- crit1 + crit2 + crit3

    if (!silent) {

      cat(paste0(
        "\n*******************************\n",
        "* Weighted Absolute Deviation *\n",
        "*******************************\n",
        "\n",
        "Following 99%-VaR, 97.5%-VaR and 97.5%-ES.\n",
        "\n",
        "WAD: ", sprintf("%.4f", crit_total), "\n\n"
      ))

    }

    return(invisible(crit_total))

  } else {

    warning("Not able to compute WAD criterion. 99%-VaR, 97.5%-VaR and 97.5%-ES required.")

    return(NULL)

  }


})
