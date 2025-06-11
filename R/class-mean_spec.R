setClass("mean_spec",
  slots = c(
    orders = "numeric",
    long_memo = "logical",
    include_mean = "logical"
  ),
  prototype = list(
    orders = c(0, 0),
    long_memo = FALSE,
    include_mean = TRUE
  )
)

setValidity("mean_spec", function(object){
  if (order_checker_arma(object@orders)) {
    "@orders must be of length two and both elements must be at least zero"
  } else if (lmemo_checker(object@long_memo)) {
    "@long_memo must be of length one"
  } else if (mean_checker(object@include_mean)) {
    "@include_mean must be of length one"
  } else {
    TRUE
  }
})

#'Specification of Conditional Mean Models
#'
#'Specify the model for the conditional mean in a dual model, where the
#'conditional mean is modelled through an ARMA or a FARIMA model and the
#'conditional standard deviations through a GARCH-type model simultaneously.
#'
#'@param orders a two-element numeric vector with the model
#'orders; the first element is the autoregressive order \eqn{p^{*}}, while
#'the second element is the moving-average order \eqn{q^{*}}.
#'@param long_memo a logical value that indicates whether the
#'long-memory version of the model should be considered or not.
#'@param include_mean a logical value indicating whether or
#'not to include the constant unconditional mean in the estimation
#'procedure; for \code{include_mean = FALSE}, the unconditional
#'mean of the series is fixed to zero and not being estimated.
#'@return
#'An object of class \code{"mean_spec"} is returned.
#'
#'@details
#'Let \eqn{\left\{y_t\right\}}, with \eqn{t \in \mathbb{Z}} as the time index,
#'be a theoretical time series that follows
#'\deqn{\beta(B)(1- B)^{D}(y_t - \mu)=\alpha(B)r_t,}
#'where \eqn{\beta(B) = 1 - \sum_{i=1}^{p^{*}}\beta_i B^{i}} and
#'\eqn{\alpha(B) = 1 + \sum_{j=1}^{q^{*}}\alpha_j B^{j}} are the AR- and MA-polynomials
#'of orders \eqn{p^{*}} and \eqn{q^{*}}, respectively, with real coefficients
#'\eqn{\beta_i}, \eqn{i=1,\dots,p^{*}}, and \eqn{\alpha_j}, \eqn{j=1,\dots,q^{*}}.
#'\eqn{B} is the backshift operator. \eqn{\beta(B)} and \eqn{\alpha(B)} are
#'commonly assumed to be without common roots and to have roots outside of
#'the unit circle.
#'Furthermore, \eqn{\mu} is a real-valued coefficient representing the unconditional
#'mean in \eqn{\left\{y_t\right\}}. \eqn{D \in [0, 0.5)} is the fractional
#'differencing parameter. \eqn{\left\{r_t\right\}} is a zero-mean (weak) white
#'noise process, for example a member of the GARCH-models (with mean set to zero)
#'presented in this package (see the
#'descriptions in \link{fEGarch_spec}, \link{fiaparch}, \link{figarch}, etc.).
#'
#'The for \eqn{D=0}, which can be achieved through \code{long_memo = FALSE},
#'the formulas above describe an autoregressive moving-average (ARMA) model.
#'For \eqn{D \in (0, 0.5)}, they describe a fractionally integrated ARMA (FARIMA)
#'model.
#'
#'@export
#'
#'@examples
#'mean_spec()
#'mean_spec(orders = c(1, 1))
#'
mean_spec <- function(
    orders = c(0, 0),
    long_memo = FALSE,
    include_mean = TRUE
) {

  new("mean_spec",
    orders = orders,
    long_memo = long_memo,
    include_mean = include_mean
  )

}

#'Accessors for Class \code{"mean_spec"}
#'
#'Access and change elements in objects of class \code{"mean_spec"}.
#'The method names represent the name of the element to access /
#'manipulate.
#'
#'@param x the input object or object to modify.
#'@param value the value to modify the object \code{x} with.
#'
#'@details
#'These methods are intended to be used for accessing or manipulating
#'individual elements of objects of class \code{"mean_spec"}.
#'
#'@name mean_spec_methods
#'@aliases orders,mean_spec-method
#'
#'@return
#'These methods return an object of class \code{"mean_spec"}.
#'
#'@export
#'
#'@examples
#'test_obj <- mean_spec()
#'orders(test_obj)
#'orders(test_obj) <- c(1, 1)
#'orders(test_obj)
#'
setMethod("orders", "mean_spec", function(x) {x@orders})
#'@export
#'@rdname mean_spec_methods
#'@aliases long_memo,mean_spec-method
setMethod("long_memo", "mean_spec", function(x) {x@long_memo})
#'@export
#'@rdname mean_spec_methods
#'@aliases include_mean,mean_spec-method
setMethod("include_mean", "mean_spec", function(x) {x@include_mean})
#'@export
#'@rdname mean_spec_methods
setMethod("orders<-", "mean_spec", function(x, value) {x@orders <- value; validObject(x); x})
#'@export
#'@rdname mean_spec_methods
setMethod("long_memo<-", "mean_spec", function(x, value) {x@long_memo <- value; validObject(x); x})
#'@export
#'@rdname mean_spec_methods
setMethod("include_mean<-", "mean_spec", function(x, value) {x@include_mean <- value; validObject(x); x})
