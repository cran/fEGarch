setClass("locpol_spec",
  slots = c(
    poly_order = "numeric",
    kernel_order = "numeric",
    boundary_method = "character",
    bwidth = "ANY"
  ),
  prototype = list(
    poly_order = 3,
    kernel_order = 1,
    boundary_method = "extend",
    bwidth = NULL
  )
)

setValidity("locpol_spec", function(object){
  if (poly_checker_locpol(object@poly_order)) {
    "@poly_order must be of length one and either 1 or 3"
  } else if (kernel_checker_locpol(object@kernel_order)) {
    "@kernel_order must be of length one and one of the integers from 0, 1, 2, 3"
  } else if (bmethod_checker_locpol(object@boundary_method)) {
    '@boundary_method must be of length one and either "shorten" or "extend"'
  } else if (bwidth_checker_locpol(object@bwidth)) {
    '@bwidth must be either NULL (for automatic bandwidth selection) or a single numeric value between 0 and 0.5'
  } else {
    TRUE
  }
})

#'Specification of Nonparametric Local Polynomial Models
#'
#'Specify the nonparametric local polynomial model part in a
#'semiparametric volatility model.
#'
#'@param poly_order a single numeric value, in detail either \code{1}
#'or \code{3}, that represents the local polynomial order.
#'@param kernel_order a single numeric value representing
#'the smoothness of the underlying kernel function;
#'available are \code{0} (uniform), \code{1} (epanechnikov),
#'\code{2} (bisquare), and \code{3} (triweight).
#'@param boundary_method a single character value indicating
#'the smoothing concept to use at boundary points; for
#'\code{"extend"}, the smoothing window is extended toward the
#'interior by the amount that is lost toward the boundary; for
#'\code{"shorten"}, there is no compensation in the smoothing window
#'toward the interior for the loss of window width toward the boundary,
#'i.e. the total smoothing window width reduces more and more when
#'getting closer to the first and last time points.
#'@param bwidth the smoothing bandwidth; for NULL, i.e. the default, an automated
#'bandwidth selection is employed; otherwise a single numeric value between
#'0 and 0.5 must be provided.
#'
#'@return
#'An object of class \code{"locpol_spec"} is returned.
#'
#'@details
#'Assume that a time series \eqn{\{r_t\}}, \eqn{t=1,\dots,n}, follows
#'\deqn{r_t = \mu + \sigma_t \eta_t,}
#'where \eqn{\mu = E(r_t)} and \eqn{\eta_t} are independent and identically
#'distributed random variables with mean zero and variance one. \eqn{\sigma_t > 0}
#'are total volatilities composed of \eqn{s(x_t)}, a smooth, deterministic
#'scale function in the unconditional variance over time (with \eqn{x_t} being
#'the rescaled time on the interval \eqn{[0, 1]}), and of \eqn{\lambda_t},
#'the conditional standard deviation in \eqn{\zeta_t=\lambda_t\eta_t}, so that
#'\eqn{\sigma_t = s(x_t)\lambda_t}, or alternatively \eqn{r_t = \mu + s(x_t)\zeta_t}.
#'It is assumed that the unconditional variance of the \eqn{\zeta_t} is one.
#'
#'The package's estimation of \eqn{\sigma_t} is based on the following relations:
#'
#'\eqn{r_t^{*} = r_t - \mu},
#'\eqn{y_t=\ln\left[\left(r_t^{*}\right)^2\right]},
#'\eqn{C_{\mu}=E\left[\ln\left(\zeta_t^2\right)\right]},
#'\eqn{m(x_t) = \ln\left[s^2 (x_t)\right] + C_{\mu}},
#'\eqn{\xi_t = \ln\left(\zeta_t^2\right) - C_{\mu}}, so that
#'\deqn{y_t = m(x_t)+\xi_t,}
#'where \eqn{m} describes a smooth, deterministic trend in \eqn{y_t}.
#'Nonparametric estimation of \eqn{m} and subsequent retransformation
#'allows to obtain a suitable estimate of the scale function \eqn{s}
#'in \eqn{r_t}. Following Feng et al. (2022) and Letmathe et al. (2023),
#'we employ local polynomial regression with automatically selected
#'bandwidth (specially for the time-series context). The function
#'\code{locpol_spec} allows to set the basic characteristics of the
#'local polynomial estimator considered, like the order of polynomial
#'used in the local regressions, and the kernel function order. After
#'the scale function has been estimated, a zero-mean GARCH-type model
#'can be fitted to the estimated \eqn{\zeta_t}.
#'
#'Depending on whether \eqn{\zeta_t} is assumed to follow a short-memory
#'or a long-memory model, the bandwidth selection algorithm in the
#'local polynomial regression step differs and follows either
#'Feng et al. (2022) and Letmathe et al. (2023). The algorithm
#'selection is done automatically based on the remaining model
#'specifications in the call to the estimation functions like
#'\code{\link{fEGarch}}.
#'
#'@export
#'
#'@references
#'\itemize{
#'\item{Feng, Y., Gries, T., Letmathe, S., & Schulz, D. (2022). The smoots Package in R for Semiparametric Modeling of
#'Trend Stationary Time Series. The R Journal,
#'14(1), 182-195. URL: https://journal.r-project.org/articles/RJ-2022-017/.}
#'\item{Letmathe, S., Beran, J., & Feng, Y. (2023). An extended exponential SEMIFAR model with application
#'in R. Communications in Statistics - Theory and Methods,
#'53(22), 7914–7926. DOI: 10.1080/03610926.2023.2276049.}
#'}
#'
#'@examples
#'locpol_spec()
#'locpol_spec(poly_order = 1)
#'locpol_spec(kernel_order = 2)
#'
locpol_spec <- function(
    poly_order = c(3, 1),
    kernel_order = c(1, 0, 2, 3),
    boundary_method = c("extend", "shorten"),
    bwidth = NULL
) {

  if (all(poly_order %in% c(3, 1))) {
    poly_order <- poly_order[[1]]
  } else {
    stop("@poly_order must be one from 1 or 3.")
  }

  if (all(kernel_order %in% c(1, 0, 2, 3))) {
    kernel_order <- kernel_order[[1]]
  } else {
    stop("@kernel_order must be one from 0, 1, 2 or 3.")
  }

  boundary_method <- match.arg(boundary_method)

  new("locpol_spec",
    poly_order = poly_order,
    kernel_order = kernel_order,
    boundary_method = boundary_method,
    bwidth = bwidth
  )

}

#'Accessors for Class \code{"locpol_spec"}
#'
#'Access and change elements in objects of class \code{"locpol_spec"}.
#'The method names represent the name of the element to access /
#'manipulate.
#'
#'@param x the input object or object to modify.
#'@param value the value to modify the object \code{x} with.
#'
#'@details
#'These methods are intended to be used for accessing or manipulating
#'individual elements of objects of class \code{"locpol_spec"}.
#'
#'@name locpol_spec_methods
#'@aliases poly_order,locpol_spec-method
#'
#'@return
#'These methods return an object of class \code{"locpol_spec"}.
#'
#'@export
#'
#'@examples
#'test_obj <- locpol_spec()
#'poly_order(test_obj)
#'poly_order(test_obj) <- 1
#'poly_order(test_obj)
#'
setMethod("poly_order", "locpol_spec", function(x) {x@poly_order})
#'@export
#'@rdname locpol_spec_methods
#'@aliases kernel_order,locpol_spec-method
setMethod("kernel_order", "locpol_spec", function(x) {x@kernel_order})
#'@export
#'@rdname locpol_spec_methods
#'@aliases boundary_method,locpol_spec-method
setMethod("boundary_method", "locpol_spec", function(x) {x@boundary_method})
#'@export
#'@rdname locpol_spec_methods
#'@aliases bwidth,locpol_spec-method
setMethod("bwidth", "locpol_spec", function(x) {x@bwidth})
#'@export
#'@rdname locpol_spec_methods
setMethod("poly_order<-", "locpol_spec", function(x, value) {x@poly_order <- value; validObject(x); x})
#'@export
#'@rdname locpol_spec_methods
setMethod("kernel_order<-", "locpol_spec", function(x, value) {x@kernel_order <- value; validObject(x); x})
#'@export
#'@rdname locpol_spec_methods
setMethod("boundary_method<-", "locpol_spec", function(x, value) {x@boundary_method <- value; validObject(x); x})
#'@export
#'@rdname locpol_spec_methods
setMethod("bwidth<-", "locpol_spec", function(x, value) {x@bwidth <- value; validObject(x); x})
