#'Estimation of a Broad Family of EGARCH Models
#'
#'A library of quasi maximum-likelihood estimation (QMLE) methods for
#'fitting various short- and long-memory models from a broad family of
#'exponential generalized autoregressive conditional heteroskedasticity (EGARCH) models.
#'For the purpose of comparison, a FIAPARCH (fractionally integrated
#'asymmetric power ARCH), a FIGARCH (fractionally integrated GARCH),
#'a FITGARCH, a FIGJR-GARCH and their short-memory variants
#'can be implemented as well.
#'
#'@details
#'\code{fEGarch} is an R package for estimating a broad family of
#'EGARCH models (Feng et al., 2025; Ayensu et al., 2025)
#'including both short- and long-memory as well as a
#'selection of varying transformations for the asymmetry and the
#'magnitude term in such a model, for example in form of the
#'FIMLog-GARCH (Feng et al., 2023). Log-GARCH specifications can be
#'implemented as well as a special case of the broad EGARCH family.
#'The six most common conditional distributions are supported, namely
#'a normal distribution, a \eqn{t}- distribution, a generalized error
#'distribution, as well as the skewed variants of these three
#'distributions. Furthermore, as a novelty, an average Laplace (AL)
#'distribution (see for example Feng et al., 2025) and its skewed
#'version are provided as well.
#'The main functions to implement these models are
#'\code{\link{fEGarch_spec}} in combination with
#'\code{\link{fEGarch}}. Further details on these models can also
#'be found in the documentation of these two functions. For
#'convenience, further specification functions for particular
#'submodels are available as well: \code{\link{egarch_spec}},
#'\code{\link{loggarch_spec}}, \code{\link{megarch_spec}},
#'\code{\link{mloggarch_spec}}, \code{\link{fiegarch_spec}},
#'\code{\link{filoggarch_spec}},
#'\code{\link{fimegarch_spec}} and \code{\link{fimloggarch_spec}}.
#'
#'As a popular alternative for the sake of comparison, a FIAPARCH
#'model can be fitted as well using \code{\link{fiaparch}}. The
#'corresponding documentation page also includes further information
#'on the model. Similarly, \code{\link{figarch}} can be utilized for
#'fitting FIGARCH models, \code{\link{fitgarch}} for fitting
#'FITGARCH models and \code{\link{figjrgarch}} for fitting
#'FIGJR-GARCH models. A general function for the estimation of additional
#'GARCH-type models, including the aforementioned additional models
#'as well as their short-memory variants, is
#'\code{\link{garchm_estim}}.
#'
#'In addition, the package provides functionalities in order to
#'simultaneously model the conditional mean (using either autoregressive
#'moving-average (ARMA) models or fractionally integrated ARMA (FARIMA)
#'models) alongside the
#'conditional variance. For this purpose, the function
#'\code{\link{mean_spec}} can be utilized and its result needs to be
#'passed to \code{\link{fEGarch}} alongside the result of either
#'\code{\link{fEGarch_spec}} or one of its wrappers.
#'
#'Further options include the specification of semiparametric
#'volatility models (see also Ayensu et al., 2025), where a smooth, nonparametric scale function
#'is at first estimated and removed from an observed series, before
#'estimating a parametric model. The scale estimation is currently
#'done through automated local polynomial regression with designated
#'bandwidth selection algorithms under short memory and long memory.
#'
#'@section Main Functions:
#'The main functions of the package are:
#'\describe{
#'\item{\code{\link{fEGarch_spec}}:}{setting the model specifications for
#'a model from the broader EGARCH family,}
#'\item{\code{\link{mean_spec}}:}{setting the model specifications for
#'the conditional mean,}
#'\item{\code{\link{fEGarch}}:}{fitting a model from the broad family of EGARCH
#'models given a model specification and an observation series,}
#'\item{\code{\link{garchm_estim}}:}{fitting a GARCH-type model selectable from a standard GARCH,
#'a GJR-GARCH, a TGARCH, an APARCH, a FIGARCH, a FIGJR-GARCH, a FITGARCH and a FIAPARCH,}
#'\item{\code{\link{fEGarch_sim}}:}{simulating from an EGARCH family model,}
#'\item{\code{\link{fiaparch_sim}}:}{simulating from a FIAPARCH model,}
#'\item{\code{\link{figarch_sim}}:}{simulating from a FIGARCH model,}
#'\item{\code{\link{figjrgarch_sim}}:}{simulating from a FIGJR-GARCH model,}
#'\item{\code{\link{fitgarch_sim}}:}{simulating from a FITGARCH model,}
#'\item{\code{\link{aparch_sim}}:}{simulating from an APARCH model,}
#'\item{\code{\link{garch_sim}}:}{simulating from a GARCH model,}
#'\item{\code{\link{gjrgarch_sim}}:}{simulating from a GJR-GARCH model,}
#'\item{\code{\link{tgarch_sim}}:}{simulating from a TGARCH model,}
#'\item{\code{\link{predict,fEGarch_fit-method}}:}{multistep point forecasts of
#'the conditional mean and the conditional standard deviation,}
#'\item{\code{\link{predict_roll,fEGarch_fit-method}}:}{rolling point forecasts of
#'the conditional mean and the conditional standard deviation over a test set.}
#'\item{\code{\link{measure_risk}}:}{value at risk and expected shortfall computation for various model specifications.}
#'\item{\code{\link{find_dist}}:}{fits all eight distributions considered in this
#'package to a supposed iid series and selects the best fitted distribution
#'following either BIC (the default) or AIC.}
#'\item{\code{\link{backtest_suite,fEGarch_risk-method}}:}{runs a selection of functions
#'for backtesting VaR and ES.}
#'}
#'
#'@section Datasets:
#'The package includes a few datasets. Follow the corresponding links to the
#'documentation of the datasets to find additional information including the
#'sources.
#'\describe{
#'\item{\code{\link{UKinflation}}:}{monthly inflation rate of the UK.}
#'\item{\code{\link{SP500}}:}{daily log-returns of the S&P 500 index.}
#'}
#'
#'@section License:
#'The package is distributed under the General Public License v3
#'([GPL-3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))).
#'
#'@references
#'\itemize{
#'\item{Ayensu, O. K., Feng, Y., & Schulz, D. (2025). Recent Extensions of Exponential GARCH Models:
#'Theory and Application. Forthcoming preprint, Paderborn University.}
#'\item{Baillie, R., Bollerslev, T., & Mikkelsen, H. O. (1996). Fractionally integrated generalized autoregressive conditional heteroskedasticity.
#'Journal of Econometrics,
#'74(1), 3-30. DOI: 10.1016/S0304-4076(95)01749-6.}
#'\item{Bollerslev, T. (1986). Generalized autoregressive conditional heteroskedasticity.
#'Journal of Econometrics, 31(3): 307-327. DOI: 10.1016/0304-4076(86)90063-1.}
#'\item{Bollerslev, T., & Mikkelsen, H. O. (1996). Modeling and pricing long memory in stock market volatility. Journal of Econometrics,
#'73(1), 151–184. DOI: 10.1016/0304-4076(95)01749-6.}
#'\item{Conrad, C., & Haag, B. R. (2006). Inequality constraints in the fractionally
#'integrated GARCH model. Journal of Financial Econometrics, 4(3):
#'413-449. DOI: 10.1093/jjfinec/nbj015.}
#'\item{Conrad, C., & Karanasos, M. (2006). The impulse response function
#'of the long memory GARCH process. Economics Letters, 90(1):
#'34-41. DOI: 10.1016/j.econlet.2005.07.001.}
#'\item{Ding, Z., Granger, C. W. J., & Engle, R. F. (1993). A long memory property of stock market returns
#'and a new model. Journal of Empirical Finance, 1(1):
#'83-106. DOI: 10.1016/0927-5398(93)90006-D.}
#'\item{Engle, R. F. (1982). Autoregressive Conditional Heteroscedasticity with Estimates of the Variance of United Kingdom Inflation.
#'Econometrica, 50(4): 987-1007. DOI: 10.2307/1912773.}
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
#'\item{Glosten, L. R., Jagannathan, R., & Runkle, D. E. (1993). On The Relation between The Expected Value and The
#'Volatility of Nominal Excess Return on stocks. Journal of Finance 48(5), 1779-1801.
#'DOI: 10.1111/j.1540-6261.1993.tb05128.x.}
#'\item{Karanasos, M., Psaradakis, Z., & Sola, M. (2004). On the autocorrelation properties of
#'long-memory GARCH processes. Journal of Time Series Analysis, 25(2):
#'265-281. DOI: 10.1046/j.0143-9782.2003.00349.x.}
#'\item{Letmathe, S., Beran, J., & Feng, Y. (2023). An extended exponential SEMIFAR model with application
#'in R. Communications in Statistics - Theory and Methods,
#'53(22), 7914–7926. DOI: 10.1080/03610926.2023.2276049.}
#'\item{Milhoj, A. (1987). A Multiplicative Parameterization of ARCH Models. University of Copenhagen, Denmark.}
#'\item{Missiakoulis, S. (1983). Sargan Densities: Which One?. Journal of Econometrics, 23(2): 223-233.
#'DOI: 10.1016/0304-4076(93)90078-J}
#'\item{Nelson, D. B. (1991). Conditional Heteroskedasticity in Asset Returns: A New Approach. Econometrica,
#'59(2), 347–370. DOI: 10.2307/2938260.}
#'\item{Nielsen, M. O., & Noel, A. L. (2021). To infinity and beyond: Efficient computation of ARCH(\eqn{\infty})
#'models. Journal of Time Series Analysis,
#'42(3), 338–354. DOI: 10.1111/jtsa.12570.}
#'\item{Pantula, S. G. (1986). Modeling the persistence of conditional variances: A comment. Econometric Reviews, 5(1),
#'71-74. DOI: 10.1080/07474938608800089.}
#'\item{Tse, Y. K. (1987). A Note On Sargan Densities. Journal of Econometrics, 34(3): 349-354.
#'DOI: 10.1016/0304-4076(87)90017-0}
#'\item{Tse, Y. K. (1998). The conditional heteroskedasticity of the
#'yen-dollar exchange rate. Journal of Applied Econometrics, 13(1):
#'49-55. DOI: 10.1002/(SICI)1099-1255(199801/02)13:1<49::AID-JAE459>3.0.CO;2-O.}
#'\item{Zakoian, J.-M. (1994). Threshold heteroskedastic models. Journal of Economic Dynamics and Control, 18(5):
#'931-955. DOI: 10.1016/0165-1889(94)90039-6.}
#'}
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Department of Economics, Paderborn
#'University), \cr
#'Author and Package Creator
#'\item Yuanhua Feng (Department of Economics, Paderborn
#'University), \cr
#'Author
#'\item Christian Peitz (Financial Intelligence Unit, German Government), \cr
#'Author
#'\item Oliver Kojo Ayensu (Department of Economics, Paderborn
#'University), \cr
#'Author
#'}
#'
#'@useDynLib fEGarch
#'@name fEGarch-package
#'@importFrom Rcpp sourceCpp
#'@importFrom methods new validObject
#'@importFrom stats integrate pnorm rnorm rt runif time uniroot
#'@importFrom utils tail
#'@importFrom ggplot2 autoplot
#'@import methods
#'
"_PACKAGE"



utils::globalVariables(c(
  "D_low", "D_up", "D_start", "D_par_name"
))

utils::globalVariables(c(
  "add_low", "add_up", "add_pars_start"
))

utils::globalVariables(c(
  "aic", "bic"
))

utils::globalVariables(c(
  "ar_lb", "ar_ub", "ar_start", "ma_lb", "ma_ub", "ma_start", "arma_grab_fun",
  "names_ar", "names_ma", "n_arma_pars", "lm_arma", "low_ar", "up_ar",
  "pg0", "qg0", "p_ar", "q_ma"
))

utils::globalVariables(c(
  "d_low", "d_name", "d_spar", "d_up"
))

utils::globalVariables(c(
  "cmean_fun", "dfun", "int_fun", "int_fun_asy", "int_fun_mag"
))

utils::globalVariables(c(
  "omega_low", "omega_start", "omega_up"
))

utils::globalVariables(c(
  "phi_low", "phi_start", "phi_up", "psi_low", "psi_start", "psi_up"
))

utils::globalVariables(c(
  "etransf_init", "lnsig2_init", "sig2_init", "sigd_init", "presample_val"
))

utils::globalVariables(c(
  "sval", "lval", "uval"
))

utils::globalVariables(c(
  "test_obs", "train_obs", "rt_core", "rt_core_o"
))

utils::globalVariables(c(
  "incl_mean", "mods", "pows"
))

utils::globalVariables(c(
  "sc_delta", "scale_const"
))

utils::globalVariables(c(
  "serrors"
))

