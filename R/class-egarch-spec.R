setClass("base_garch_spec",
  slots = c(
    orders = "numeric",
    long_memo = "logical",
    cond_dist = "character"
  ),
  prototype = list(
    orders = c(1, 1),
    long_memo = FALSE,
    cond_dist = "norm"
  )
)

setValidity("base_garch_spec", function(object){
  if (order_checker(object@orders)) {
    "@orders must be of length two and both elements must be at least one"
  } else if (lmemo_checker(object@long_memo)) {
    "@long_memo must be of length one"
  } else if (distr_checker(object@cond_dist)) {
    "@cond_dist must be of length one and from the available model options"
  } else {
    TRUE
  }
})

base_garch_spec <- function(
    orders = c(1, 1),
    long_memo = TRUE,
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {

  cond_dist <- tryCatch(
    expr = {match.arg(cond_dist)},
    error = function(e1) {
      stop('"cond_dist" must be kept at default or must be of length 1 and one of the defaults')
    }
  )

  new("base_garch_spec",
    orders = orders,
    long_memo = long_memo,
    cond_dist = cond_dist
  )

}

setClass("loggarch_type_spec",
  contains = "base_garch_spec"
)

setClass("egarch_type_spec",
  contains = "base_garch_spec",
  slots = c(
    powers = "numeric",
    modulus = "logical"
  ),
  prototype = list(
    powers = c(NA_real_, NA_real_),
    modulus = c(NA, NA)
  )
)

setValidity("egarch_type_spec", function(object){
  if (power_checker(object@powers)) {
    "@powers should have length two"
  } else if (modulus_checker(object@modulus)){
    "@modulus must be a logical two-element vector"
  } else {
    TRUE
  }
})

#'Subspecification of EGARCH Family Models
#'
#'Two common subspecifications of the broad EGARCH family,
#'namely for a EGARCH-type and a Log-GARCH-type model.
#'
#'@param orders a two-element numeric vector with the model
#'orders; the first element is the order \eqn{p} for the term based on
#'\eqn{\ln\left(\sigma_t^2\right)}, i.e. the log-transformed
#'conditional variance, while the second element is the order \eqn{q} for
#'the innovation-based term (see Details below for more information).
#'@param long_memo a logical value that indicates whether the
#'long-memory version of the model should be considered or not.
#'@param cond_dist a character value stating the underlying
#'conditional distribution to consider; available are a normal
#'distribution (\code{"norm"}), a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution
#'(\code{"ged"}), an average Laplace distribution (\code{"ald"})
#'and the skewed versions of them
#'(\code{"snorm"}, \code{"sstd"}, \code{"sged"}, \code{"sald"}).
#'@param powers a two-element numeric vector that states the
#'exponents in the power-transformations of the asymmetry and the
#'magnitude terms in that order (see Details for more information).
#'@param modulus a two-element logical vector indicating if the
#'innovations in the asymmetry and the magnitude terms (in that
#'order) should use a modulus transformation (see Details for
#'more information).
#'
#'@details
#'These are wrappers for \code{\link{fEGarch_spec}}.
#'\code{egarch_type_spec()} is when setting \code{model_type} in
#'\code{\link{fEGarch_spec}} to \code{eGARCH}.
#'\code{loggarch_type_spec()} is a shortcut when setting
#'\code{model_type = "loggarch"}. See \code{\link{fEGarch_spec}}
#'for further details.
#'
#'@return
#'An object of class \code{"egarch_type_spec"} or
#'\code{"loggarch_type_spec"} is returned.
#'
#'@export
#'@rdname fEGarch-subspecs
egarch_type_spec <- function(
  orders = c(1, 1),
  long_memo = TRUE,
  cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
  powers = c(1, 1),
  modulus = c(FALSE, FALSE)
) {

  cond_dist <- tryCatch(
    expr = {match.arg(cond_dist)},
    error = function(e1) {
      stop('"cond_dist" must be kept at default or must be of length 1 and one of the defaults')
    }
  )

  new("egarch_type_spec",
    orders = orders,
    long_memo = long_memo,
    cond_dist = cond_dist,
    powers = powers,
    modulus = modulus
  )

}

#'@export
#'@rdname fEGarch-subspecs
loggarch_type_spec <- function(
  orders = c(1, 1),
  long_memo = TRUE,
  cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {

  cond_dist <- tryCatch(
    expr = {match.arg(cond_dist)},
    error = function(e1) {
      stop('"cond_dist" must be kept at default or must be of length 1 and one of the defaults')
    }
  )

  new("loggarch_type_spec",
    orders = orders,
    long_memo = long_memo,
    cond_dist = cond_dist
  )

}

#'General EGARCH Family Model Specification
#'
#'Create an object with specifications for a model from
#'the broader EGARCH family.
#'
#'@param model_type a character value (either \code{"egarch"} or
#'\code{"loggarch"}) that indicates the type of model to
#'implement (see Details for more information).
#'@param orders a two-element numeric vector with the model
#'orders; the first element is the order \eqn{p} for the term based on
#'\eqn{\ln\left(\sigma_t^2\right)}, i.e. the log-transformed
#'conditional variance, while the second element is the order \eqn{q} for
#'the innovation-based term (see Details below for more information).
#'@param long_memo a logical value that indicates whether the
#'long-memory version of the model should be considered or not.
#'@param cond_dist a character value stating the underlying
#'conditional distribution to consider; available are a normal
#'distribution (\code{"norm"}), a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution
#'(\code{"ged"}), an average Laplace distribution (\code{"ald"})
#'and the skewed versions of them
#'(\code{"snorm"}, \code{"sstd"}, \code{"sged"}, \code{"sald"}).
#'@param powers a two-element numeric vector that states the
#'exponents in the power-transformations of the asymmetry and the
#'magnitude terms in that order (see Details for more information).
#'@param modulus a two-element logical vector indicating if the
#'innovations in the asymmetry and the magnitude terms (in that
#'order) should use a modulus transformation (see Details for
#'more information).
#'
#'@export
#'
#'@details
#'Let \eqn{\left\{r_t\right\}}, with \eqn{t \in \mathbb{Z}} as the
#'time index, be a theoretical time series that follows
#'\deqn{r_t=\mu+\sigma_t \eta_t \text{ with } \eta_t \sim \text{IID}(0,1), \text{ where}}
#'\deqn{\ln\left(\sigma_t^2\right)=\omega_{\sigma}+\theta(B)g(\eta_{t-1}).}
#'Here, \eqn{\eta_t\sim\text{IID}(0,1)} means that the innovations
#'\eqn{\eta_t} are independent and identically distributed (iid) with mean zero
#'and variance one, whereas \eqn{\sigma_t > 0} are the conditional standard
#'deviations in \eqn{r_t}. Note that \eqn{\ln\left(\cdot\right)} denotes the natural logarithm.
#'Moreover, \eqn{B} is the backshift operator and
#'\eqn{\theta(B) = 1 +\sum_{i=1}^{\infty}\theta_i B^{i}}, where
#'\eqn{\theta_i}, \eqn{i=1,2,\dots}, are real-valued coefficients.
#'\eqn{g\left(\eta_{t-1}\right)} is a suitable function in \eqn{\eta_{t-1}}.
#'Generally, \eqn{\left\{g\left(\eta_{t}\right)\right\}} should be an iid
#'zero-mean sequence with finite variance.
#'We have \eqn{\mu = E\left(r_t\right)} as a real-valued parameter.
#'The real-valued parameter
#'\eqn{\omega_{\sigma}} is in fact
#'\eqn{\omega_{\sigma}=E\left[\ln\left(\sigma_t^2\right)\right]}.
#'This previous set of
#'equations defines the broader family of EGARCH models (Feng et al., 2025; Ayensu et al., 2025), from which
#'subtypes are described in the following that depend on the choice of
#'\eqn{g}.
#'
#'\eqn{\textbf{Type I}:}
#'
#'We have \eqn{\theta(B) = \phi^{-1}(B)(1-B)^{-d}\psi(B)}, where
#'\deqn{\phi(B) = 1-\sum_{i=1}^{p}\phi_{i}B^{i} \text{ and}}
#'\deqn{\psi(B) = 1+\sum_{j=1}^{q-1}\psi_{j}B^{j},}
#'are characteristic polynomials
#'with real coefficients \eqn{\phi_{0},\dots,\phi_{p},\psi_{0},\dots,\psi_{q-1}},
#'by fixing \eqn{\phi_{0}=\psi_{0}=1}, and without common roots. Furthermore,
#'the
#'fractional differencing parameter is \eqn{d \in [0,1]}.
#'
#'\eqn{g\left(\cdot\right)} can be defined in different ways.
#'Following a type I specification (\code{model_type = "egarch"}), we have
#'\deqn{g\left(\eta_t\right)=\kappa \left\{g_a\left(\eta_t\right) - E\left[g_a\left(\eta_t\right)\right] \right\} + \gamma\left\{g_m\left(\eta_t\right)-E\left[ g_m\left(\eta_t\right)\right]\right\}}
#'with \eqn{g_a\left(\eta_t\right)} and \eqn{g_m\left(\eta_t\right)} being
#'suitable transformations of \eqn{\eta_t} and
#'where \eqn{\kappa} and \eqn{\gamma} are two additional real-valued
#'parameters. In case of a simple (FI)EGARCH, we have
#'\eqn{g_a\left(\eta_t\right) = \eta_t} (and therefore with
#'\eqn{E\left[g_a\left(\eta_t\right)\right] = 0}) and
#'\eqn{g_m\left(\eta_t\right) = \left|\eta_t \right|}.
#'
#'Generally, we consider two cases:
#'\deqn{g_{1,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\left|\eta_t\right|^{p_a}/p_a \text{ and}}
#'\deqn{g_{2,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\left[\left(\left|\eta_t\right| + 1\right)^{p_a} - 1\right] / p_{a}}
#'whereas
#'\deqn{g_{1,m}\left(\eta_t\right)=\left|\eta_t\right|^{p_m}/p_m \text{ and}}
#'\deqn{g_{2,m}\left(\eta_t\right)=\left[\left(\left|\eta_t\right| + 1\right)^{p_m} - 1\right] / p_{m}.}
#'Note that \eqn{\text{sgn}\left(\eta_t\right)} denotes the sign of
#'\eqn{\eta_t}. \eqn{g_{1,\cdot}} incorporates a power transformation
#'and \eqn{g_{2,\cdot}} a modulus transformation together with a power
#'transformation. The choices
#'\eqn{g_{1,a}} and \eqn{g_{2,a}} correspond to
#'setting the first element in \code{modulus} to \code{FALSE} or
#'\code{TRUE}, respectively, under
#'\code{model_type = "egarch"}, where \eqn{p_a} can be selected
#'via the first element in \code{powers}. As a special case, for
#'\eqn{p_a = 0}, a log-transformation is employed and the division
#'through \eqn{p_a} is dropped, i.e.
#'\eqn{g_{1,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\ln\left(\left|\eta_t\right|\right)}
#'and \eqn{g_{2,a}\left(\eta_t\right)=\text{sgn}\left(\eta_t\right)\ln\left(\left|\eta_t\right|+1\right)}
#'are employed for \eqn{p_a=0}. Completely analogous
#'thoughts hold for \eqn{g_{1,m}} and \eqn{g_{2,m}} and the second
#'elements in the arguments \code{modulus} and \code{powers}. The aforementioned
#'model family is a type I model selectable through
#'\code{model_type = "egarch"}. Simple (FI)EGARCH models are given through
#'selection of \eqn{g_{1,a}(\cdot)} and \eqn{g_{1,m}(\cdot)}
#'in combination with \eqn{p_a = p_m = 1}. Another of such type I models
#'is the FIMLog-GARCH (Feng et al., 2023), where instead
#'\eqn{g_{2,a}(\cdot)} and \eqn{g_{2,m}(\cdot)} are selected
#'with \eqn{p_a = p_m = 0}.
#'
#'\eqn{\textbf{Type II}:}
#'
#'As additional specifications of a Log-GARCH and a FILog-GARCH, which belong
#'to the broader EGARCH family, we redefine
#'\deqn{\psi(B) = 1+\sum_{j=1}^{q}\psi_{j}B^{j}}
#'now as a polynomial of order \eqn{q} and
#'\deqn{\ln\left(\sigma_t^2\right)=\omega_{\sigma}+\left[\phi^{-1}(B)(1-B)^{-d}\psi(B) - 1\right]\xi_t,}
#'where
#'\deqn{\xi_t=\ln\left(\eta_t^2\right)-E\left[\ln\left(\eta_t^2\right)\right].}
#'Everything else is defined as before. Since
#'\deqn{\phi^{-1}(B)(1-B)^{-d}\psi(B) - 1=\sum_{i=1}^{\infty}\gamma_i B^{i}=\gamma_1 B\left[\sum_{i=1}^{\infty}(\gamma_i/\gamma_1)B^{i-1}\right]=\gamma_1 B\left[1+\sum_{i=1}^{\infty}\theta_i B^{i}\right]=\theta(B)\gamma_1 B,}
#'where \eqn{\theta_i = \gamma_{i+1}/\gamma_{1}}, \eqn{i=0,1,\dots}, and by defining
#'\deqn{g\left(\eta_{t-1}\right)=\gamma_1\left\{\ln\left(\eta_{t-1}^2\right)-E\left[\ln\left(\eta_{t-1}^2\right)\right]\right\} = 2\gamma_1\left\{\ln\left(\left|\eta_{t-1}\right|\right)-E\left[\ln\left(\left|\eta_{t-1}\right|\right)\right]\right\},}
#'the equation of \eqn{\ln\left(\sigma_t^2\right)} can be stated to be
#'\deqn{\ln\left(\sigma_t^2\right)=\omega_{\sigma}+\theta(B)g(\eta_{t-1})}
#'as in the broad EGARCH family at the very beginning. Therefore, Log-GARCH
#'and FI-Log-GARCH models are equivalent to the type I models, where
#'\eqn{\kappa = 0} with usage of \eqn{g_{1,m}} with \eqn{p_m = 0} and where
#'\eqn{\gamma = 2\gamma_1}. Nonetheless,
#'in this package, the type II models make use of the more common
#'parameterization of \eqn{\ln\left(\sigma_t^2\right)} stated at the
#'beginning of the type II model description.
#'
#'This describes the
#'established Log-GARCH models as part of the broad EGARCH family
#'(type II models; \code{model_type = "loggarch"}).
#'
#'\eqn{\textbf{General information}:}
#'
#'While the arguments \code{powers} and \code{modulus} are only relevant
#'under a type I model, i.e. for \code{model_type = "egarch"}, the arguments
#'\code{orders}, \code{long_memo} and \code{cond_dist} are
#'meaningful for both \code{model_type = "egarch"}
#'and \code{model_type = "loggarch"}, i.e. both under type I and II models.
#'The first element of the two-element vector orders is the order \eqn{p},
#'while the second element is the order \eqn{q}. Furthermore, for
#'\code{long_memo = TRUE}, the mentioned models are kept as they are, while
#'for \code{long_memo = FALSE} the parameter \eqn{d} is set to zero.
#'\code{cond_dist} controls the conditional distribution.
#'The unconditional mean \eqn{\mu} is controlled via the function
#'\code{\link{mean_spec}}; for \code{include_mean = FALSE} therein, \eqn{\mu} is not being estimated and fixed
#'to zero; its default is however \code{include_mean = TRUE}.
#'
#'See also the closely related spec-functions that immediately create
#'specifications of specific submodels of the broad EGARCH family. These
#'functions are \code{egarch_spec()}, \code{fiegarch_spec()},
#'\code{loggarch_spec()},\code{filoggarch_spec()}, \code{megarch_spec()},
#'\code{mloggarch_spec()} and \code{mafiloggarch_spec()}, which are all
#'wrappers for \code{fEGarch_spec()}.
#'
#'See the references section for sources on the EGARCH (Nelson, 1991),
#'FIEGARCH (Bollerslev and Mikkelsen, 1996),
#'Log-GARCH (Geweke, 1986; Pantula, 1986; Milhoj, 1987)
#'and FILog-GARCH (Feng et al., 2020) models.
#'
#'@return
#'An object of either class \code{"egarch_type_spec"} or
#'\code{"loggarch_type_spec"} is returned, depending on the
#'choice for the input argument \code{model_type}.
#'
#'@references
#'\itemize{
#'\item{Ayensu, O. K., Feng, Y., & Schulz, D. (2025). Recent Extensions of Exponential GARCH Models:
#'Theory and Application. Forthcoming preprint, Paderborn University.}
#'\item{Bollerslev, T., & Mikkelsen, H. O. (1996). Modeling and pricing long memory in stock market volatility. Journal of Econometrics,
#'73(1), 151–184. DOI: 10.1016/0304-4076(95)01749-6.}
#'\item{Feng, Y., Beran, J., Ghosh, S., & Letmathe, S. (2020). Fractionally integrated Log-GARCH with application to value at risk and expected shortfall.
#'Working Papers CIE No. 137, Paderborn University, Center for International Economics.
#'URL: http://groups.uni-paderborn.de/wp-wiwi/RePEc/pdf/ciepap/WP137.pdf.}
#'\item{Feng, Y., Gries, T., & Letmathe, S. (2023). FIEGARCH, modulus asymmetric FILog-GARCH
#'and trend-stationary dual long memory time series. Working Papers CIE No. 156, Paderborn University.
#'URL: https://econpapers.repec.org/paper/pdnciepap/156.htm.}
#'\item{Feng, Y., Peitz, C., & Siddiqui, S. (2025). A few useful members of the EGARCH-family
#'with short- or long-memory in volatility. Unpublished working paper at Paderborn University.}
#'\item{Geweke, J. (1986). Modeling the persistence of conditional variances: A comment. Econometric Reviews, 5(1),
#'57-61. DOI: 10.1080/07474938608800088.}
#'\item{Milhoj, A. (1987). A Multiplicative Parameterization of ARCH Models. University of Copenhagen, Denmark.}
#'\item{Nelson, D. B. (1991). Conditional Heteroskedasticity in Asset Returns: A New Approach. Econometrica,
#'59(2), 347–370. DOI: 10.2307/2938260.}
#'\item{Pantula, S. G. (1986). Modeling the persistence of conditional variances: A comment. Econometric Reviews, 5(1),
#'71-74. DOI: 10.1080/07474938608800089.}
#'}
#'
#'@examples
#'# EGARCH(1, 1) with cond. normal distribution
#'spec1 <- fEGarch_spec()
#'# EGARCH(2, 1) with cond. t-distribution
#'spec2 <- fEGarch_spec(orders = c(2, 1), cond_dist = "std")
#'# FIEGARCH(1, 1) with cond. normal distribution
#'spec3 <- fEGarch_spec(long_memo = TRUE)
#'# MEGARCH(1, 1) with cond. generalized error distribution
#'spec4 <- fEGarch_spec(modulus = c(TRUE, FALSE), powers = c(0, 1))
#'# Some unnamed specification
#'spec5 <- fEGarch_spec(
#'  model_type = "egarch",
#'  orders = c(1, 1),
#'  long_memo = TRUE,
#'  cond_dist = "std",
#'  powers = c(0.25, 0.75),
#'  modulus = c(TRUE, FALSE)
#')
#'

fEGarch_spec <- function(
  model_type = c("egarch", "loggarch"),
  orders = c(1, 1),
  long_memo = FALSE,
  cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald"),
  powers = c(1, 1),
  modulus = c(FALSE, FALSE)
) {
  model_type <- tryCatch(
    expr = {match.arg(model_type)},
    error = function(e1) {
      stop('"model_type" must be kept at default or must be of length 1 and one of the defaults')
    }
  )
  if (is.null(model_type) || is.na(model_type)) {
    model_type <- "egarch"
  }
  if (is.null(orders) || all(is.na(orders))) {
    orders <- c(1, 1)
  }
  if (is.null(long_memo) || is.na(long_memo)) {
    long_memo <- FALSE
  }
  if ((is.null(modulus) || all(is.na(modulus))) && (model_type == "egarch")) {
    modulus <- c(FALSE, FALSE)
  }
  cond_dist <- tryCatch(
    expr = {match.arg(cond_dist)},
    error = function(e1) {
      stop('"cond_dist" must be kept at default or must be of length 1 and one of the defaults')
    }
  )
  if ((is.null(powers) || all(is.na(powers))) && (model_type == "egarch")) {
    powers <- c(1, 1)
  }

  if (model_type == "egarch") {
    out <- egarch_type_spec(
      orders = orders,
      long_memo = long_memo,
      cond_dist = cond_dist,
      powers = powers,
      modulus = modulus
    )
  } else if (model_type == "loggarch") {
    out <- loggarch_type_spec(
      orders = orders,
      long_memo = long_memo,
      cond_dist = cond_dist
    )
  }

  out

}



#'Accessors for Classes \code{"base_garch_spec"} and \code{"egarch_spec"}
#'
#'Access and change elements in objects of class \code{"egarch_spec"}.
#'The method names represent the name of the element to access /
#'manipulate.
#'
#'@param x the input object or object to modify.
#'@param value the value to modify the object \code{x} with.
#'
#'@details
#'These methods are intended to be used for accessing or manipulating
#'individual elements of objects of class \code{"egarch_spec"}.
#'
#'@name spec_methods
#'@aliases orders,base_garch_spec-method
#'
#'@return
#'These methods either return an object of class \code{"egarch_type_spec"}
#'or \code{"loggarch_type_spec"} or the corresponding element of
#'object of such class objects.
#'
#'@export
#'
#'@examples
#'test_obj <- egarch_spec()
#'orders(test_obj)
#'orders(test_obj) <- c(2, 1)
#'orders(test_obj)
#'
setMethod("orders", "base_garch_spec", function(x) {x@orders})
#'@export
#'@rdname spec_methods
#'@aliases powers,egarch_type_spec-method
setMethod("powers", "egarch_type_spec", function(x) {x@powers})
#'@export
#'@rdname spec_methods
#'@aliases long_memo,base_garch_spec-method
setMethod("long_memo", "base_garch_spec", function(x) {x@long_memo})
#'@export
#'@rdname spec_methods
#'@aliases modulus,egarch_type_spec-method
setMethod("modulus", "egarch_type_spec", function(x) {x@modulus})
#'@export
#'@rdname spec_methods
#'@aliases cond_dist,base_garch_spec-method
setMethod("cond_dist", "base_garch_spec", function(x) {x@cond_dist})

#'@export
#'@rdname spec_methods
setMethod("orders<-", "base_garch_spec", function(x, value) {x@orders <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("powers<-", "egarch_type_spec", function(x, value) {x@powers <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("long_memo<-", "base_garch_spec", function(x, value) {x@long_memo <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("modulus<-", "egarch_type_spec", function(x, value) {x@modulus <- value; validObject(x); x})
#'@export
#'@rdname spec_methods
setMethod("cond_dist<-", "base_garch_spec", function(x, value) {x@cond_dist <- value; validObject(x); x})

#'EGARCH Family Submodel Specification
#'
#'Wrappers of \code{fEGarch_spec()} that create
#'specifications of specific submodels of the broad
#'EGARCH family.
#'
#'@param orders a two-element numeric vector with the model
#'orders.
#'@param cond_dist a character value stating the underlying
#'conditional distribution to consider; available are a normal
#'distribution (\code{"norm"}), a \eqn{t}-distribution
#'(\code{"std"}), a generalized error distribution
#'(\code{"ged"}), an average Laplace distribution (\code{"ald"})
#'and the skewed versions of them
#'(\code{"snorm"}, \code{"sstd"}, \code{"sged"}, \code{"sald"}).
#'
#'@export
#'
#'@details
#'Available are shortcut specification functions for
#'EGARCH \code{egarch_spec()}, FIEGARCH \code{fiegarch_spec()},
#'MEGARCH \code{megarch_spec()}, Log-GARCH \code{loggarch_spec()},
#'FILog-GARCH \code{filoggarch_spec()}, MLog-GARCH \code{mloggarch_spec()},
#'FIMEGARCH \code{fimegarch_spec()} and
#'FIMLog-GARCH \code{fimloggarch_spec()}.
#'
#'The following descriptions are following the descriptions in the
#'documentation of the more general \code{fEGarch_spec()}. Please go there
#'first to understand the following descriptions on the arguments of
#'\code{fEGarch_spec()} to obtain these wrappers.
#'
#'\eqn{\textbf{EGARCH:}}
#'
#'\code{model_type = "egarch"}, \code{long_memo = FALSE},
#'\code{powers = c(1, 1)}, \code{modulus = c(FALSE, FALSE)}
#'
#'\eqn{\textbf{FIEGARCH:}}
#'
#'\code{model_type = "egarch"}, \code{long_memo = TRUE},
#'\code{powers = c(1, 1)}, \code{modulus = c(FALSE, FALSE)}
#'
#'\eqn{\textbf{MEGARCH:}}
#'
#'\code{model_type = "egarch"}, \code{long_memo = FALSE},
#'\code{powers = c(0, 1)}, \code{modulus = c(TRUE, FALSE)}
#'
#'\eqn{\textbf{Log-GARCH:}}
#'
#'\code{model_type = "loggarch"}, \code{long_memo = FALSE}
#'
#'\eqn{\textbf{FILog-GARCH:}}
#'
#'\code{model_type = "loggarch"}, \code{long_memo = TRUE}
#'
#'\eqn{\textbf{MLog-GARCH:}}
#'
#'\code{model_type = "egarch"}, \code{long_memo = FALSE},
#'\code{powers = c(0, 0)}, \code{modulus = c(TRUE, TRUE)}
#'
#'\eqn{\textbf{FIMLog-GARCH:}}
#'
#'\code{model_type = "egarch"}, \code{long_memo = TRUE},
#'\code{powers = c(0, 0)}, \code{modulus = c(TRUE, TRUE)}
#'
#'\eqn{\textbf{FIMEGARCH:}}
#'
#'\code{model_type = "egarch"}, \code{long_memo = TRUE},
#'\code{powers = c(0, 1)}, \code{modulus = c(TRUE, FALSE)}
#'
#'@name submodel-specs
#'
#'@return
#'Depending on the spec-fun function, either an object of class
#'\code{"egarch-type-spec"} or \code{"loggarch-type-spec"} is returned.
#'
#'@examples
#'spec <- megarch_spec(cond_dist = "std")
#'

megarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {
  fEGarch_spec(
    model_type = "egarch",
    orders = orders,
    long_memo = FALSE,
    cond_dist = cond_dist,
    powers = c(0, 1),
    modulus = c(TRUE, FALSE)
  )
}

#'@export
#'@rdname submodel-specs
fimegarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {
  fEGarch_spec(
    model_type = "egarch",
    orders = orders,
    long_memo = TRUE,
    cond_dist = cond_dist,
    powers = c(0, 1),
    modulus = c(TRUE, FALSE)
  )
}

#'@export
#'@rdname submodel-specs
egarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {
  fEGarch_spec(
    model_type = "egarch",
    orders = orders,
    long_memo = FALSE,
    cond_dist = cond_dist,
    powers = c(1, 1),
    modulus = c(FALSE, FALSE)
  )
}

#'@export
#'@rdname submodel-specs
fiegarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {
  fEGarch_spec(
    model_type = "egarch",
    orders = orders,
    long_memo = TRUE,
    cond_dist = cond_dist,
    powers = c(1, 1),
    modulus = c(FALSE, FALSE)
  )
}

#'@export
#'@rdname submodel-specs
mloggarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {
  fEGarch_spec(
    model_type = "egarch",
    orders = orders,
    long_memo = FALSE,
    cond_dist = cond_dist,
    powers = c(0, 0),
    modulus = c(TRUE, TRUE)
  )
}

#'@export
#'@rdname submodel-specs
fimloggarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {
  fEGarch_spec(
    model_type = "egarch",
    orders = orders,
    long_memo = TRUE,
    cond_dist = cond_dist,
    powers = c(0, 0),
    modulus = c(TRUE, TRUE)
  )
}

#'@export
#'@rdname submodel-specs
loggarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {
  fEGarch_spec(
    model_type = "loggarch",
    orders = orders,
    long_memo = FALSE,
    cond_dist = cond_dist
  )
}

#'@export
#'@rdname submodel-specs
filoggarch_spec <- function(
    orders = c(1, 1),
    cond_dist = c("norm", "std", "ged", "ald", "snorm", "sstd", "sged", "sald")
) {
  fEGarch_spec(
    model_type = "loggarch",
    orders = orders,
    long_memo = TRUE,
    cond_dist = cond_dist
  )
}


# Given a parameter vector "theta" and settings on what parameters are included,
# grab the parameters from the vector (following an EGARCH-specification)
egarch_grab_pars <- function(theta, p, q, incl_mean, d_par, extra_par, skew_par) {

  imean <- !incl_mean

  mu <- if (incl_mean) theta[[1]] else 0
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q - 1

  phi <- theta[(step0 + 1):step1]

  psi <- if (q > 1) c(1, theta[(step1 + 1):step2]) else 1

  kappa <- theta[[(step2 + 1)]]
  gamma <- theta[[step2 + 2]]
  step2p3 <- step2 + 3
  d <- if (d_par) theta[[step2p3]] else 0
  shape <- if (extra_par) theta[[step2p3 + d_par]] else 0
  skew <- if (skew_par) theta[[step2p3 + d_par + extra_par]] else 0

  list(
    "mu" = mu,
    "omega_sig" = omega_sig,
    "phi" = phi,
    "psi" = psi,
    "kappa" = kappa,
    "gamma" = gamma,
    "d" = d,
    "shape" = shape,
    "skew" = skew
  )

}

arma_grab_pars <- function(theta, incl_mean, pg0, qg0, p, q) {

  imean <- !incl_mean
  mu <- if (incl_mean) theta[[1]] else 0
  step0 <- 2 - imean
  ar <- if (pg0) {
    step1 <- step0 + p - 1
    step2 <- step1 + 1
    theta[step0:step1]
  } else {
    step2 <- step0
    numeric(0)
  }
  ma <- if (qg0) {
    step3 <- step2 + q - 1
    theta[step2:step3]
  } else {
    numeric(0)
  }

  list(
    mu = mu,
    ar = ar,
    ma = ma
  )

}

farima_grab_pars <- function(theta, incl_mean, pg0, qg0, p, q) {

  imean <- !incl_mean
  mu <- if (incl_mean) theta[[1]] else 0
  step0 <- 2 - imean
  ar <- if (pg0) {
    step1 <- step0 + p - 1
    step2 <- step1 + 1
    theta[step0:step1]
  } else {
    step2 <- step0
    numeric(0)
  }
  ma <- if (qg0) {
    step3 <- step2 + q - 1
    D <- theta[[step3 + 1]]
    theta[step2:step3]
  } else {
    D <- theta[[step2]]
    numeric(0)
  }

  list(
    mu = mu,
    ar = ar,
    ma = ma,
    D = D
  )

}

# Calculate the two constants (i.e. the expectations) in the EGARCH formula
expect_calc_egarch <- function(int_fun_asy, int_fun_mag, dfun, shape, skew, pows) {
  E_asy <- stats::integrate(f = int_fun_asy, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun, pow = pows[[1]],
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]
  E_mag <- stats::integrate(f = int_fun_mag, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun, pow = pows[[2]],
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]
  list("E_asy" = E_asy, "E_mag" = E_mag)
}

# Given parameters, observations, etc., calculate conditional standard
# deviations following a short-memory EGARCH specification
sigt_egarch_short_R <- function(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par, trunc, presample, ...) {

  # Make ellipsis objects available
  list2env(list(...), envir = environment())

  all_pars <- egarch_grab_pars(
    theta = theta,
    p = p,
    q = q,
    incl_mean = incl_mean,
    d_par = FALSE,
    extra_par = extra_par,
    skew_par = skew_par
  )
  mu <- all_pars$mu
  omega_sig <- all_pars$omega_sig
  phi <- all_pars$phi
  psi <- all_pars$psi
  kappa <- all_pars$kappa
  gamma <- all_pars$gamma
  shape <- all_pars$shape
  skew <- all_pars$skew

  alpha <- psi * kappa
  beta <- psi * gamma

  E_vals <- expect_calc_egarch(
    int_fun_asy = int_fun_asy,
    int_fun_mag = int_fun_mag,
    dfun = dfun,
    shape = shape,
    skew = skew,
    pows = pows)

  E_asy <- E_vals$E_asy
  E_mag <- E_vals$E_mag

  Elnsig2 <- omega_sig
  omega <- Elnsig2 * (1 - sum(phi))

  sigt <- c(sigt_egarch_shortCpp(
    x = x,
    mu = mu,
    omega = omega,
    phi = phi,
    alpha = alpha,
    beta = beta,
    E_asy = E_asy,
    E_mag = E_mag,
    powers = pows, modulus = mods, mode = mode,
    lnsig2_init = lnsig2_init
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape, cmeans = rep(mu, length(sigt)))

}

# Given parameters, observations, etc., calculate conditional standard
# deviations following a long-memory EGARCH specification
sigt_egarch_long_R <- function(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par, trunc, presample, ...) {

  # Make ellipsis objects available
  list2env(list(...), envir = environment())

  all_pars <- egarch_grab_pars(
    theta = theta,
    p = p,
    q = q,
    incl_mean = incl_mean,
    d_par = TRUE,
    extra_par = extra_par,
    skew_par = skew_par
  )
  mu <- all_pars$mu
  omega_sig <- all_pars$omega_sig
  phi <- all_pars$phi
  psi <- all_pars$psi
  kappa <- all_pars$kappa
  gamma <- all_pars$gamma
  d <- all_pars$d
  shape <- all_pars$shape
  skew <- all_pars$skew

  psi_2 <- if (q > 1) psi[-1] else numeric(0)
  n <- length(x)
  coef_inf <- ma_infty(ar = phi, ma = psi_2, d = d,
                                     max_i = trunc - 1 - presample)

  E_vals <- expect_calc_egarch(
    int_fun_asy = int_fun_asy,
    int_fun_mag = int_fun_mag,
    dfun = dfun,
    shape = shape,
    skew = skew,
    pows = pows)

  E_asy <- E_vals$E_asy
  E_mag <- E_vals$E_mag

  Elnsig2 <- omega_sig

  sigt <- c(sigt_egarch_longCpp(
    x = x,
    mu = mu,
    coef_inf = coef_inf,
    kappa = kappa,
    gamma = gamma,
    E_asy = E_asy,
    E_mag = E_mag,
    Elnsig2 = Elnsig2,
    powers = pows,
    modulus = mods,
    mode = mode
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape, cmeans = rep(mu, length(sigt)))

}

goal_fun_creator_egarch <- function(
    theta,
    x,
    dfun1,
    dfun2,
    p,
    q,
    p_ar,
    q_ma,
    incl_mean,
    extra_par,
    skew_par,
    int_fun_asy,
    int_fun_mag,
    sigt_fun,
    pows, mods,
    mode,
    lnsig2_init,
    trunc,
    presample
) {

  res_list <- sigt_fun(
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
    dfun = dfun1,
    int_fun_asy = int_fun_asy,
    int_fun_mag = int_fun_mag,
    pows = pows,
    mods = mods,
    mode = mode,
    lnsig2_init = lnsig2_init
  )

  nllhood_calc(res_list, x, dfun2)

}

#Fitting Method for Type I EGARCH-Family Models
#
#Fits an EGARCH-family model of Type I. The method
#is not being exported.
#
#@param spec the generic is currently without use.
#@param rt the generic is currently without use.
#@param drange the generic is currently without use.
#@param meanspec the generic is currently without use.
#@param Drange the generic is currently without use.
#@param n_test the generic is currently without use.
#@param start_pars the generic is currently without use.
#@param LB the generic is currently without use.
#@param UB the generic is currently without use.
#@param control the generic is currently without use.
#@param parallel the generic is currently without use.
#@param ncores the generic is currently without use.
#@param trunc the generic is currently without use.
#@param presample the generic is currently without use.
#@param Prange the generic is currently without use.
#
#@return
#Returns an object of class \code{"fEGarch_fit_egarch"}.
#
#
setMethod("fEGarch_fit", "egarch_type_spec",
  function(spec = egarch_spec(), rt, drange = c(0, 1), meanspec = mean_spec(), Drange = c(0, 1),
    n_test = 0, start_pars = NULL, LB = NULL, UB = NULL, control = list(), parallel = TRUE,
    ncores = max(1, future::availableCores() - 1), trunc = "none", presample = 50, Prange = c(1, 5),
    skip_vcov = FALSE) {

  # Make rt have sample variance 1; more robust for
  # differently scaled data; results like parameters
  # and the log-likelihood with regard to the original data
  # can be computed easily later,
  # as, if required, the corresponding values returned for the scaled data
  # are just multiplied by a constant
  # or shifted by a constant. This transformation just shifts the
  # entire log-likelihood by a constant and therefore does not affect
  # the shape of the log-likelihood.

  # Makes available
  # - rt_core_o: pure original obs. of the training set,
  # - rt_core: pure rescaled obs.,
  # - scale_const: the scaling constant used,
  # - test_obs: the test obs. with keeping initial formatting,
  # - train_obs: the train obs. with keeping initial formatting
  initial_split(rt, n_test)


  n <- length(rt_core)
  if (Prange[[1]] < 1) {Prange[[1]] <- 1}
  trunc_i <- trunc

  # Obtain model orders
  order <- orders(spec)
  p <- order[[1]]
  q <- order[[2]]
  # Conditional distribution
  cond_d <- cond_dist(spec)
  # Long-memory setting (for GARCH-part)
  lm <- long_memo(spec)

  # Settings for ARMA / FARIMA part in the mean; makes objects available in this
  # environment
  identify_arma_settings(meanspec = meanspec, Drange = Drange)

  if ((lm || lm_arma) && is.character(trunc) && trunc == "none") {
    trunc <- n + presample - 1
  } else if ((lm || lm_arma) && is.numeric(trunc) && trunc >= n + presample) {
    trunc <- n + presample - 1
    message("trunc reduced to n + ", presample, " - 1")
  }


  # Look up internally saved density functions
  # from "data-raw/lookup.R"
  dfun1 <- fun1_selector(cond_d)
  dfun2 <- fun2_selector(cond_d)

  # Look up internally saved settings for additional parameters
  # from "data-raw/lookup.R"
  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  # Adjust (if ARMA relevant)
  model_type <- c("egarch", "fiegarch")[[lm + 1]]
  sigt_fun_adj <- adj_sigt_fun(
    p_ar = p_ar, q_ma = q_ma, pg0 = pg0, qg0 = qg0,
    lm_arma = lm_arma, incl_mean = incl_mean,
    n_arma_pars = n_arma_pars, grab_fun = arma_grab_fun,
    arma_fun = cmean_fun,
    model_type = model_type
  )

  # Modulus settings
  mod <- modulus(spec)
  mod_asy <- mod[[1]]
  mod_mag <- mod[[2]]

  # Power parameter settings
  pow <- powers(spec)
  pow_asy <- pow[[1]]
  pow_mag <- pow[[2]]
  pow_asy_check <- (pow_asy != 0) + 1
  pow_mag_check <- (pow_mag != 0) + 1

  # Look up internally saved functions to integrate over to numerically
  # in order to obtain relevant expectation values for the volatility
  # function (from "data-raw/lookup.R")
  int_fun_asy <- exp_funs[["asy"]][[mod_asy + 1]][[pow_asy_check]]
  int_fun_mag <- exp_funs[["mag"]][[mod_mag + 1]][[pow_mag_check]]

  # Look up the mode for the volatility calculation functions based
  # on power and modulus settings from "data-raw/lookup.R"
  mode <- mode_mat[pow_asy_check, pow_mag_check]

  # Initialize ln(sig_t^2) with some reasonable, fixed value
  lnsig2_init <- log(var(rt_core))

  # Create the goal function to optimize over
  # (this is intermediate goal function, as dfun1 and dfun2
  # must be kept flexible for (skewed) average Laplace distribution
  # in order to be able to fix its parameter at different integer values)
  goal_fun_s <- function(theta, rt, dfun1, dfun2) {
    goal_fun_creator_egarch(
      theta = theta,
      x = rt,
      dfun1 = dfun1,
      dfun2 = dfun2,
      p = p,
      q = q,
      p_ar = p_ar,
      q_ma = q_ma,
      incl_mean = incl_mean,
      extra_par = extra_par,
      skew_par = skew_par,
      int_fun_asy = int_fun_asy,
      int_fun_mag = int_fun_mag,
      sigt_fun = sigt_fun_adj,
      pows = pow,
      mods = mod,
      mode = mode,
      lnsig2_init = lnsig2_init,
      trunc = trunc,
      presample = presample
    )
  }

  # Get cond. distribution restrictions and starting values
  # as objects sval, lval and uval in the environment
  cond_d_vals <- cond_d_par_restr(cond_d)

  # Make further starting parameter values and restrictions available
  # using internally saved table function
  pars_and_restr <- start_pars_and_restr[[model_type]](drange, p, q)
  list2env(pars_and_restr, envir = environment())
  rm(pars_and_restr)

  if (is.null(start_pars) || is.null(LB) || is.null(UB)) {
    mean_rt <- list(NULL, mean(rt_core))[[incl_mean + 1]]
    min_rt <- list(NULL, suppressWarnings(0.2 * min(rt_core - mean_rt) + mean_rt))[[incl_mean + 1]]
    max_rt <- list(NULL, suppressWarnings(0.2 * max(rt_core - mean_rt) + mean_rt))[[incl_mean + 1]]
    skew_start <- list(NULL, 0.98)[[skew_par + 1]]
    skew_low <- list(NULL, 1e-15)[[skew_par + 1]]
    skew_up <- list(NULL, Inf)[[skew_par + 1]]

    start_LB_UB <- set_start_LB_UB(
      mean_val = mean_rt,
      ar_start = ar_start,
      ma_start = ma_start,
      D_start = D_start,
      omega_val = omega_start,
      phi_start = phi_start,
      psi_start = psi_start,
      add_pars_start = add_pars_start,
      d_start = d_spar,
      sval = sval,
      skew_start = skew_start,
      mean_low = min_rt,
      ar_low = ar_lb,
      ma_low = ma_lb,
      D_low = D_low,
      omega_low = omega_low,
      phi_low = phi_low,
      psi_low = psi_low,
      add_low = add_low,
      d_low = d_low,
      lval = lval, skew_low = skew_low,
      mean_up = max_rt,
      ar_up = ar_ub, ma_up = ma_ub,
      D_up = D_up,
      omega_up = omega_up,
      phi_up = phi_up, psi_up = psi_up,
      add_up = add_up,
      d_up = d_up,
      uval = uval, skew_up = skew_up,
      start_pars = start_pars, LB = LB, UB = UB
    )


  }

  low_phi <- 2 + n_arma_pars
  up_phi <- 1 + n_arma_pars + p

  # For stationarity, restrict phi / AR values
  constr_fun <- function(theta, rt) {
    phi <- c(1, -theta[low_phi:up_phi])
    c1 <- abs(polyroot(phi)) - 1
    if (pg0) {
      ar <- c(1, -theta[low_ar:up_ar])
      c2 <- abs(polyroot(ar)) - 1
    } else {
      c2 <- numeric(0)
    }
    c(c1, c2)
  }

  if (is.null(control$trace)) {control$trace <- FALSE}

  p_all <- p + p_ar

  ineqLB <- rep(1e-6, p_all)
  ineqUB <- rep(Inf, p_all)

  if (cond_d %in% c("ald", "sald")) {

    # Makes "P_sel" and "result" available
    ald_fit(
      Prange = Prange,
      parallel = parallel,
      ncores = ncores,
      dfun1 = dfun1,
      dfun2 = dfun2,
      goal_fun_s = goal_fun_s,
      start_pars = start_pars,
      LB = LB, UB = UB,
      ineqfun = constr_fun,
      ineqLB = ineqLB,
      ineqUB = ineqUB,
      rt_core = rt_core,
      control = control
    )

  } else {

    P_sel <- NULL

    # Use dfun1 and dfun2 directly without further adjustments
    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun1 = dfun1, dfun2 = dfun2)
    }

    result <- tryCatch(
      expr = {suppressWarnings(Rsolnp::solnp(
        pars = start_pars,
        fun = goal_fun,
        LB = LB,
        UB = UB,
        ineqfun = constr_fun,
        ineqLB = ineqLB,    # abs. values of roots must be outside of unit circle
        ineqUB = ineqUB,
        control = control,
        rt = rt_core
      ))},
      error = function(e1) {
        stop("Error during optimization. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
      }
    )

  }

  if (result$convergence %in% c(1, 2)) {
    warning("Convergence failed. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
  }

  pars <- result$pars

  # Compute log-likelihood which is shifted from the obtained negative
  # log-likelihood at the optimum; "-tail(result$values, 1)" is the
  # log-likelihood at the optimum for the scaled data, which is too large by
  # "n * log(scale_const)" in comparison to original data
  llhood <- -tail(result$values, 1) - n * log(scale_const)


  psi_names <- if (q > 1) {paste0("psi", 1:(q - 1))} else {NULL}

  extra_par_name <- switch(
    cond_d,
    "std" = "df",
    "sstd" = "df",
    "ged" = "shape",
    "sged" = "shape",
    "ald" = "P",
    "sald" = "P",
    character(0)
  )

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    names_ar,
    names_ma,
    D_par_name,
    "omega_sig",
    paste0("phi", 1:p),
    psi_names,
    "kappa",
    "gamma",
    d_name,
    extra_par_name,
    list(NULL, "skew")[[skew_par + 1]]
  )

  # Adjust goal function one last time for hessian computation
  if (cond_d %in% c("ald", "sald")) {

    # Fix dfun1
    dfun1_P <- function(x, shape, skew) {
      dfun1(x = x, shape = P_sel, skew = skew)
    }

    # Fix dfun2
    dfun2_P <- function(x, mu, sigt, shape, skew) {
      dfun2(x = x, mu = mu, sigt = sigt, shape = P_sel, skew = skew)
    }

    # Fix the final goal function to optimize over
    # using dfun1_P and dfun2_P
    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun1 = dfun1_P, dfun2 = dfun2_P)
    }

    fdfun1 <- dfun1_P

  } else {

    fdfun1 <- dfun1

  }


  # Compute the conditional standard deviations at the optimum

  final_series <- sigt_fun_adj(
    theta = unname(pars),
    x = rt_core,
    dfun = fdfun1,
    p = p,
    q = q,
    p_ar = p_ar,
    q_ma = q_ma,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par,
    trunc = trunc,
    presample = presample,
    int_fun_asy = int_fun_asy,
    int_fun_mag = int_fun_mag,
    pows = pow,
    mods = mod,
    mode = mode,
    lnsig2_init = lnsig2_init
  )

  # Also updates "pars", if ald or sald
  compute_vcov(
    goal_fun = goal_fun,
    pars = pars,
    rt_core = rt_core,
    scale_const = scale_const,
    incl_mean = incl_mean,
    P_sel = P_sel,
    cond_d = cond_d,
    par_names = par_names,
    model_type = model_type,
    skip_vcov = skip_vcov
  )

  # Create rescaled final parameters and results
  rescale_estim(
    pars = pars,
    rt_core_o = rt_core_o,
    sigt = final_series$sigt,
    cmeans = final_series$cmeans,
    scale_const = scale_const,
    n_arma_pars = n_arma_pars,
    incl_mean = incl_mean,
    model_type = model_type
  )

  names(pars) <- par_names

  # Makes "aic" and "bic" available
  crit_calc(pars = pars, rt = rt_core_o, llhood = llhood)

  # Apply time series formatting (if relevant) to all
  # estimated series
  series_out <- format_applier_ts(
    rt = train_obs,
    list_of_ts = list(
      "cmeans" = cmeans,
      "sigt" = sigt,
      "etat" = etat
    )
  )
  cmeans <- series_out$cmeans
  sigt <- series_out$sigt
  etat <- series_out$etat

  list_out <- fEGarch_fit_egarch(
    pars = pars,
    se = serrors,
    vcov_mat = vcov_mat,
    rt = train_obs,
    cmeans = cmeans,
    scale_fun = NULL,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_d,
    orders = order,
    powers = pow,
    modulus = mod,
    long_memo = lm,
    meanspec = meanspec,
    test_obs = test_obs,
    nonpar_model = NULL,
    trunc = trunc_i
  )

  list_out

})

loggarch_grab_pars <- function(theta, p, q, incl_mean, d_par, extra_par, skew_par) {

  imean <- !incl_mean
  mu <- if (incl_mean) theta[[1]] else 0
  step0 <- 2 - imean
  omega_sig <- theta[[step0]]
  step1 <- step0 + p
  step2 <- step1 + q
  phi <- theta[(step0 + 1):step1]
  psi <- theta[(step1 + 1):step2]

  d <- if (d_par) theta[[step2 + 1]] else 0

  shape <- if (extra_par) theta[[step2 + d_par + 1]] else 0
  skew <- if (skew_par) theta[[step2 + d_par + 1 + extra_par]] else 0

  list(
    mu = mu,
    omega_sig = omega_sig,
    phi = phi,
    psi = psi,
    d = d,
    shape = shape,
    skew = skew
  )

}

sigt_loggarch_short_R <- function(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par, trunc, presample, ...) {

  # Make ellipsis objects available
  list2env(list(...), envir = environment())

  all_pars <- loggarch_grab_pars(
    theta = theta,
    p = p,
    q = q,
    incl_mean = incl_mean,
    d_par = FALSE,
    extra_par = extra_par,
    skew_par = skew_par
  )
  mu <- all_pars$mu
  omega_sig <- all_pars$omega_sig
  phi <- all_pars$phi
  psi <- all_pars$psi
  shape <- all_pars$shape
  skew <- all_pars$skew

  Elneta2 <- stats::integrate(f = int_fun, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun,
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]

  omega <- omega_sig * (1 - sum(phi))
  l <- max(p, q)
  phi_s <- psi_s <- rep(0, l)
  phi_s[1:p] <- phi
  psi_s[1:q] <- psi
  alpha <- phi_s + psi_s

  sigt <- c(sigt_loggarch_short(
    x = x,
    mu = mu,
    omega = omega,
    phi = phi,
    psi = alpha,
    Elneta2 = Elneta2,
    lnsig2_init = lnsig2_init
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape, cmeans = rep(mu, length(sigt)))

}

sigt_loggarch_long_R <- function(theta, x, p, q, p_ar, q_ma,
              incl_mean, extra_par, skew_par, trunc, presample, ...) {

  # Make ellipsis objects available
  list2env(list(...), envir = environment())

  all_pars <- loggarch_grab_pars(
    theta = theta,
    p = p,
    q = q,
    incl_mean = incl_mean,
    d_par = TRUE,
    extra_par = extra_par,
    skew_par = skew_par
  )
  mu <- all_pars$mu
  omega_sig <- all_pars$omega_sig
  phi <- all_pars$phi
  psi <- all_pars$psi
  d <- all_pars$d
  shape <- all_pars$shape
  skew <- all_pars$skew

  Elneta2 <- stats::integrate(f = int_fun, lower = -Inf, upper = Inf, shape = shape, skew = skew,
                  dfun = dfun,
                  subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
                  stop.on.error = FALSE)[[1]]

  n <- length(x)
  coef_inf <- ma_infty(ar = phi, ma = psi, d = d,
                                     max_i = trunc - presample)[-1]

  Elnsig2 <- omega_sig

  sigt <- c(sigt_loggarch_long(
    x = x,
    mu = mu,
    coef_inf = coef_inf,
    Elneta2 = Elneta2,
    Elnsig2 = Elnsig2
  ))

  list(sigt = sigt, mu = mu, skew = skew, shape = shape, cmeans = rep(mu, length(sigt)))

}

goal_fun_creator_loggarch <- function(
    theta,
    x,
    dfun1,
    dfun2,
    p,
    q,
    p_ar,
    q_ma,
    incl_mean,
    extra_par,
    skew_par,
    int_fun,
    sigt_fun,
    lnsig2_init,
    trunc,
    presample
) {

  res_list <- sigt_fun(
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
    dfun = dfun1,
    int_fun = int_fun,
    lnsig2_init = lnsig2_init)

  nllhood_calc(res_list, x, dfun2)

}

#Fitting Method for Type II EGARCH-Family Models
#
#Fits an EGARCH-family model of Type II. The method
#is not being exported.
#
#@param spec the generic is currently without use.
#@param rt the generic is currently without use.
#@param drange the generic is currently without use.
#@param meanspec the generic is currently without use.
#@param Drange the generic is currently without use.
#@param n_test the generic is currently without use.
#@param start_pars the generic is currently without use.
#@param LB the generic is currently without use.
#@param UB the generic is currently without use.
#@param control the generic is currently without use.
#@param parallel the generic is currently without use.
#@param ncores the generic is currently without use.
#@param trunc the generic is currently without use.
#@param presample the generic is currently without use.
#@param Prange the generic is currently without use.
#
#@return
#Returns an object of class \code{"fEGarch_fit_loggarch"}.
#
#
setMethod("fEGarch_fit", "loggarch_type_spec",
  function(spec = loggarch_spec(), rt, drange = c(0, 1), meanspec = mean_spec(), Drange = c(0, 1), n_test = 0, start_pars = NULL, LB = NULL, UB = NULL, control = list(), parallel = TRUE, ncores = max(1, future::availableCores() - 1), trunc = "none", presample = 50, Prange = c(1, 5), skip_vcov = skip_vcov) {

  # Make rt have sample variance 1; more robust for
  # differently scaled data; results like parameters
  # and the log-likelihood with regard to the original data
  # can be computed easily later,
  # as, if required, the corresponding values returned for the scaled data
  # are just multiplied by a constant
  # or shifted by a constant. This transformation just shifts the
  # entire log-likelihood by a constant and therefore does not affect
  # the shape of the log-likelihood.

  # Makes available
  # - rt_core_o: pure original obs. of the training set,
  # - rt_core: pure rescaled obs.,
  # - scale_const: the scaling constant used,
  # - test_obs: the test obs. with keeping initial formatting,
  # - train_obs: the train obs. with keeping initial formatting
  initial_split(rt, n_test)


  n <- length(rt_core)
  if (Prange[[1]] < 1) {Prange[[1]] <- 1}
  trunc_i <- trunc

  # Obtain model orders
  order <- orders(spec)

  p <- order[[1]]
  q <- order[[2]]
  # Conditional distribution
  cond_d <- cond_dist(spec)
  # Long-memory setting (for GARCH-part)
  lm <- long_memo(spec)

  # Settings for ARMA / FARIMA part in the mean; makes objects available in this
  # environment
  identify_arma_settings(meanspec = meanspec, Drange = Drange)

  if ((lm || lm_arma) && is.character(trunc) && trunc == "none") {
    trunc <- n + presample - 1
  } else if ((lm || lm_arma) && is.numeric(trunc) && trunc >= n + presample) {
    trunc <- n + presample - 1
    message("trunc reduced to n + ", presample, " - 1")
  }


  # Look up internally saved density functions
  # from "data-raw/lookup.R"
  dfun1 <- fun1_selector(cond_d)
  dfun2 <- fun2_selector(cond_d)

  # Look up internally saved settings for additional parameters
  # from "data-raw/lookup.R"
  skew_par <- lookup_table[["skew_par"]][[cond_d]]
  extra_par <- lookup_table[["extra_par"]][[cond_d]]

  # Adjust (if ARMA relevant)
  model_type <- c("loggarch", "filoggarch")[[lm + 1]]
  sigt_fun_adj <- adj_sigt_fun(
    p_ar = p_ar, q_ma = q_ma, pg0 = pg0, qg0 = qg0,
    lm_arma = lm_arma, incl_mean = incl_mean,
    n_arma_pars = n_arma_pars, grab_fun = arma_grab_fun,
    arma_fun = cmean_fun,
    model_type = model_type
  )

  # Define function to integrate over to compute
  # the expectation of ln(e_t^2) numerically in the
  # numerical optimization; quantity depends on parameter
  # settings for everything but a standard normal distribution
  int_fun <- function(x, shape, skew, dfun) {
    log(x^2) * dfun(x, shape = shape, skew = skew)
  }

  # Initialize ln(sig_t^2) with some reasonable, fixed value
  lnsig2_init <- log(var(rt_core))

  # Create the goal function to optimize over
  # (this is intermediate goal function, as dfun1 and dfun2
  # must be kept flexible for (skewed) average Laplace distribution
  # in order to be able to fix its parameter at different integer values)
  goal_fun_s <- function(theta, rt, dfun1, dfun2) {
    goal_fun_creator_loggarch(
      theta = theta,
      x = rt,
      dfun1 = dfun1,
      dfun2 = dfun2,
      p = p,
      q = q,
      p_ar = p_ar,
      q_ma = q_ma,
      incl_mean = incl_mean,
      extra_par = extra_par,
      skew_par = skew_par,
      int_fun = int_fun,
      sigt_fun = sigt_fun_adj,
      lnsig2_init = lnsig2_init,
      trunc = trunc,
      presample = presample
    )
  }

  # Get cond. distribution restrictions and starting values
  # as objects sval, lval and uval in the environment
  cond_d_vals <- cond_d_par_restr(cond_d)

  # Make further starting parameter values and restrictions available
  # using internally saved table function
  pars_and_restr <- start_pars_and_restr[[model_type]](drange, p, q)
  list2env(pars_and_restr, envir = environment())
  rm(pars_and_restr)

  if (is.null(start_pars) || is.null(LB) || is.null(UB)) {
    mean_rt <- list(NULL, mean(rt_core))[[incl_mean + 1]]
    min_rt <- list(NULL, suppressWarnings(0.2 * min(rt_core - mean_rt) + mean_rt))[[incl_mean + 1]]
    max_rt <- list(NULL, suppressWarnings(0.2 * max(rt_core - mean_rt) + mean_rt))[[incl_mean + 1]]
    skew_start <- list(NULL, 0.98)[[skew_par + 1]]
    skew_low <- list(NULL, 1e-15)[[skew_par + 1]]
    skew_up <- list(NULL, Inf)[[skew_par + 1]]

    start_LB_UB <- set_start_LB_UB(
      mean_val = mean_rt,
      ar_start = ar_start,
      ma_start = ma_start,
      D_start = D_start,
      omega_val = omega_start,
      phi_start = phi_start,
      psi_start = psi_start,
      add_pars_start = add_pars_start,
      d_start = d_spar,
      sval = sval,
      skew_start = skew_start,
      mean_low = min_rt,
      ar_low = ar_lb,
      ma_low = ma_lb,
      D_low = D_low,
      omega_low = omega_low,
      phi_low = phi_low,
      psi_low = psi_low,
      add_low = add_low,
      d_low = d_low,
      lval = lval, skew_low = skew_low,
      mean_up = max_rt,
      ar_up = ar_ub, ma_up = ma_ub,
      D_up = D_up,
      omega_up = omega_up,
      phi_up = phi_up, psi_up = psi_up,
      add_up = add_up,
      d_up = d_up,
      uval = uval, skew_up = skew_up,
      start_pars = start_pars, LB = LB, UB = UB
    )


  }

  low_phi <- 2 + n_arma_pars
  up_phi <- 1 + n_arma_pars + p

  # For stationarity, restrict phi / AR values
  constr_fun <- function(theta, rt) {
    phi <- c(1, -theta[low_phi:up_phi])
    c1 <- abs(polyroot(phi)) - 1
    if (pg0) {
      ar <- c(1, -theta[low_ar:up_ar])
      c2 <- abs(polyroot(ar)) - 1
    } else {
      c2 <- numeric(0)
    }
    c(c1, c2)
  }

  if (is.null(control$trace)) {control$trace <- FALSE}

  p_all <- p + p_ar

  ineqLB <- rep(1e-6, p_all)
  ineqUB <- rep(Inf, p_all)

  if (cond_d %in% c("ald", "sald")) {

    # Makes "P_sel" and "result" available
    ald_fit(
      Prange = Prange,
      parallel = parallel,
      ncores = ncores,
      dfun1 = dfun1,
      dfun2 = dfun2,
      goal_fun_s = goal_fun_s,
      start_pars = start_pars,
      LB = LB, UB = UB,
      ineqfun = constr_fun,
      ineqLB = ineqLB,
      ineqUB = ineqUB,
      rt_core = rt_core,
      control = control
    )

  } else {

    P_sel <- NULL

    # Use dfun1 and dfun2 directly without further adjustments
    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun1 = dfun1, dfun2 = dfun2)
    }

    result <- tryCatch(
      expr = {suppressWarnings(Rsolnp::solnp(
        pars = start_pars,
        fun = goal_fun,
        LB = LB,
        UB = UB,
        ineqfun = constr_fun,
        ineqLB = ineqLB,    # abs. values of roots must be outside of unit circle
        ineqUB = ineqUB,
        control = control,
        rt = rt_core
      ))},
      error = function(e1) {
        stop("Error during optimization. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
      }
    )

  }

  if (result$convergence %in% c(1, 2)) {
    warning("Convergence failed. You may want to try different starting parameter and / or solver settings.", call. = FALSE)
  }

  pars <- result$pars

  # Compute log-likelihood which is shifted from the obtained negative
  # log-likelihood at the optimum; "-tail(result$values, 1)" is the
  # log-likelihood at the optimum for the scaled data, which is too large by
  # "n * log(scale_const)" in comparison to original data
  llhood <- -tail(result$values, 1) - n * log(scale_const)


  psi_names <- paste0("psi", 1:q)

  extra_par_name <- switch(
    cond_d,
    "std" = "df",
    "sstd" = "df",
    "ged" = "shape",
    "sged" = "shape",
    "ald" = "P",
    "sald" = "P",
    character(0)
  )

  par_names <- c(
    list(NULL, "mu")[[incl_mean + 1]],
    names_ar,
    names_ma,
    D_par_name,
    "omega_sig",
    paste0("phi", 1:p),
    psi_names,
    d_name,
    extra_par_name,
    list(NULL, "skew")[[skew_par + 1]]
  )

  # Adjust goal function one last time for hessian computation
  if (cond_d %in% c("ald", "sald")) {

    # Fix dfun1
    dfun1_P <- function(x, shape, skew) {
      dfun1(x = x, shape = P_sel, skew = skew)
    }

    # Fix dfun2
    dfun2_P <- function(x, mu, sigt, shape, skew) {
      dfun2(x = x, mu = mu, sigt = sigt, shape = P_sel, skew = skew)
    }

    # Fix the final goal function to optimize over
    # using dfun1_P and dfun2_P
    goal_fun <- function(theta, rt) {
      goal_fun_s(theta = theta, rt = rt, dfun1 = dfun1_P, dfun2 = dfun2_P)
    }

    fdfun1 <- dfun1_P

  } else {

    fdfun1 <- dfun1

  }

  # Compute the conditional standard deviations at the optimum
  final_series <- sigt_fun_adj(
    theta = unname(pars),
    x = rt_core,
    dfun = fdfun1,
    p = p,
    q = q,
    p_ar = p_ar,
    q_ma = q_ma,
    incl_mean = incl_mean,
    extra_par = extra_par,
    skew_par = skew_par,
    trunc = trunc,
    presample = presample,
    int_fun = int_fun,
    lnsig2_init = lnsig2_init
  )

  # Also updates "pars", if ald or sald
  compute_vcov(
    goal_fun = goal_fun,
    pars = pars,
    rt_core = rt_core,
    scale_const = scale_const,
    incl_mean = incl_mean,
    P_sel = P_sel,
    cond_d = cond_d,
    par_names = par_names,
    model_type = model_type,
    skip_vcov = skip_vcov
  )

  # Create rescaled final parameters and results
  rescale_estim(
    pars = pars,
    rt_core_o = rt_core_o,
    sigt = final_series$sigt,
    cmeans = final_series$cmeans,
    scale_const = scale_const,
    n_arma_pars = n_arma_pars,
    incl_mean = incl_mean,
    model_type = model_type
  )

  names(pars) <- par_names

  # Makes "aic" and "bic" available
  crit_calc(pars = pars, rt = rt_core_o, llhood = llhood)

  # Apply time series formatting (if relevant) to all
  # estimated series
  series_out <- format_applier_ts(
    rt = train_obs,
    list_of_ts = list(
      "cmeans" = cmeans,
      "sigt" = sigt,
      "etat" = etat
    )
  )
  cmeans <- series_out$cmeans
  sigt <- series_out$sigt
  etat <- series_out$etat

  list_out <- fEGarch_fit_loggarch(
    pars = pars,
    se = serrors,
    vcov_mat = vcov_mat,
    rt = train_obs,
    cmeans = cmeans,
    scale_fun = NULL,
    sigt = sigt,
    etat = etat,
    llhood = llhood,
    inf_criteria = c("aic" = aic, "bic" = bic),
    cond_dist = cond_d,
    orders = order,
    long_memo = lm,
    meanspec = meanspec,
    test_obs = test_obs,
    nonpar_model = NULL,
    trunc = trunc_i
  )

  list_out
  })

#'Show Method for Estimation Output
#'
#'Display estimation results from the EGARCH family in
#'a convenient way in the console.
#'
#'@param object an object returned by one of this package's fitting functions,
#'for example \code{\link{fEGarch}},
#'\code{\link{fiaparch}}, \code{\link{figarch}}, etc.
#'
#'@aliases show,fEGarch_fit_egarch-method
#'@export
#'
#'@rdname show-methods
#'
#'@return
#'Nothing is returned. Results are printed in the R-console.
#'
setMethod(
  "show",
  "fEGarch_fit_egarch",
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
    "*     Fitted EGARCH Family Model    *\n",
    "*************************************\n",
    " \n",
    "Type: egarch\n",
    "Orders: (", x@orders[[1]], ", ", x@orders[[2]], ")\n",
    "Modulus: (", x@modulus[[1]], ", ", x@modulus[[2]], ")\n",
    "Powers: (", x@powers[[1]], ", ", x@powers[[2]], ")\n",
    "Long memory: ", x@long_memo, "\n",
    "Cond. distribution: ", x@cond_dist, "\n",
    " \n",
    "ARMA orders (cond. mean): (", x@meanspec@orders[[1]], ", ", x@meanspec@orders[[2]], ")\n",
    "Long memory (cond. mean): ", x@meanspec@long_memo, "\n",
    "\n",
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

#'@export
#'@rdname show-methods
#'@aliases show,fEgarch_fit_loggarch-method
setMethod(
  "show",
  "fEGarch_fit_loggarch",
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
    "*     Fitted EGARCH Family Model    *\n",
    "*************************************\n",
    " \n",
    "Type: loggarch\n",
    "Orders: (", x@orders[[1]], ", ", x@orders[[2]], ")\n",
    "Long memory: ", x@long_memo, "\n",
    "Cond. distribution: ", x@cond_dist, "\n",
    " \n",
    "ARMA orders (cond. mean): (", x@meanspec@orders[[1]], ", ", x@meanspec@orders[[2]], ")\n",
    "Long memory (cond. mean): ", x@meanspec@long_memo, "\n",
    "\n",
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

#'Methods for Accessing Model Estimation and Forecasting Output Elements
#'
#'Accessors to access the elements of the same name in
#'output objects returned by either \code{\link{fEGarch}},
#'\code{\link{garchm_estim}}, \code{predict} or \code{predict_roll}.
#'
#'@param x an object returned by either \code{\link{fEGarch}},
#'\code{\link{garchm_estim}}, \code{predict} (for \code{sigt} and \code{cmeans})
#'or \code{predict_roll} (for \code{sigt} and \code{cmeans}).
#'
#'@details
#'Convenience methods to access the elements of the same name
#'that can otherwise be accessed via the operator \code{@} within
#'objects that inherit from class \code{"fEGarch_fit"}, which covers
#'objects returned by either \code{\link{fEGarch}},
#'\code{\link{garchm_estim}}, \code{predict} (for \code{sigt} and \code{cmeans})
#'or \code{predict_roll} (for \code{sigt} and \code{cmeans}).
#'
#'As alternatives for \code{sigt}, \code{cmeans} and \code{etat}, see also
#'\code{\link{sigma,fEGarch_fit-method}}, \code{\link{fitted,fEGarch_fit-method}}
#'and \code{\link{residuals,fEGarch_fit-method}}.
#'
#'@return
#'The element within the input object of the same name as the method
#'is returned. Depending on the element that can be a numeric vector,
#'an object of class \code{"zoo"} or a numeric matrix.
#'
#'@export
#'
#'@name accessor_methods
#'@aliases sigt,fEGarch_fit-method
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, start = "2010-01-01", end = "2012-12-31")
#'est <- fEGarch(egarch_spec(), rt)
#'
#'# Access estimated conditional standard deviations using
#'# the common operator "@" ...
#'sigt1 <- est@sigt
#'
#'# ... or use the accessor method "sigt()"
#'sigt2 <- sigt(est)
#'
#'zoo::plot.zoo(
#'  cbind("Approach 1" = sigt1, "Approach 2" = sigt2)
#')
#'
#'# Other methods
#'cmeans(est)
#'etat(est)
#'inf_criteria(est)
#'llhood(est)
#'pars(est)
#'se(est)
#'vcov_mat(est)
#'
setMethod("sigt", "fEGarch_fit", function(x) {x@sigt})
#'@export
#'@rdname accessor_methods
#'@aliases cmeans,fEGarch_fit-method
setMethod("cmeans", "fEGarch_fit", function(x) {x@cmeans})
#'@export
#'@rdname accessor_methods
#'@aliases etat,fEGarch_fit-method
setMethod("etat", "fEGarch_fit", function(x) {x@etat})
#'@export
#'@rdname accessor_methods
#'@aliases inf_criteria,fEGarch_fit-method
setMethod("inf_criteria", "fEGarch_fit", function(x) {x@inf_criteria})
#'@export
#'@rdname accessor_methods
#'@aliases llhood,fEGarch_fit-method
setMethod("llhood", "fEGarch_fit", function(x) {x@llhood})
#'@export
#'@rdname accessor_methods
#'@aliases pars,fEGarch_fit-method
setMethod("pars", "fEGarch_fit", function(x) {x@pars})
#'@export
#'@rdname accessor_methods
#'@aliases se,fEGarch_fit-method
setMethod("se", "fEGarch_fit", function(x) {x@se})
#'@export
#'@rdname accessor_methods
#'@aliases vcov_mat,fEGarch_fit-method
setMethod("vcov_mat", "fEGarch_fit", function(x) {x@vcov_mat})
