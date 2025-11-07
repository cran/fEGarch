#'Extract Fitted Conditional Standard Deviations
#'
#'An alternative to \code{\link{sigt,fEGarch_fit-method}} to extract
#'fitted conditional standard deviations from an estimation object
#'in this package.
#'
#'@param object an object either of class \code{"fEGarch_fit"} or
#'\code{"fEGarch_forecast"}.
#'
#'@importFrom stats sigma
#'@exportMethod sigma
#'@rdname sigma
#'@aliases sigma,fEGarch_fit-method
#'@export
#'
#'@details
#'Extract fitted conditional standard deviations from an estimation object
#'in this package.
#'
#'@return
#'The element within the input object with name \code{sigt} is returned.
#'Depending on the element that can be a numeric vector, an object of
#'class "zoo" or a numeric matrix.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- fEGarch(egarch_spec(), rt, n_test = 250)
#'sigma(model)
#'
setMethod("sigma", signature(object = "fEGarch_fit"), function(object) {
  fun <- methods::selectMethod("sigt", "fEGarch_fit")
  fun(object)
})

#'Extract Fitted Conditional Means
#'
#'An alternative to \code{\link{cmeans,fEGarch_fit-method}} to extract
#'fitted conditional means from an estimation object
#'in this package.
#'
#'@param object an object either of class \code{"fEGarch_fit"} or
#'\code{"fEGarch_forecast"}.
#'
#'@importFrom stats fitted
#'@exportMethod fitted
#'@rdname fitted
#'@aliases fitted,fEGarch_fit-method
#'@export
#'
#'@details
#'Extract fitted conditional means from an estimation object
#'in this package.
#'
#'@return
#'The element within the input object with name \code{cmeans} is returned.
#'Depending on the element that can be a numeric vector, an object of
#'class "zoo" or a numeric matrix.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- fEGarch(egarch_spec(), rt, n_test = 250)
#'fitted(model)
#'
setMethod("fitted", signature(object = "fEGarch_fit"), function(object) {
  fun <- methods::selectMethod("cmeans", "fEGarch_fit")
  fun(object)
})

#'Extract Standardized Residuals
#'
#'An alternative to \code{\link{etat,fEGarch_fit-method}} to extract
#'standardized residuals from an estimation object
#'in this package.
#'
#'@param object an object either of class \code{"fEGarch_fit"} or
#'\code{"fEGarch_forecast"}.
#'
#'@importFrom stats residuals
#'@exportMethod residuals
#'@rdname residuals
#'@aliases residuals,fEGarch_fit-method
#'@export
#'
#'@details
#'Extract fitted standardized residuals from an estimation object
#'in this package.
#'
#'@return
#'The element within the input object with name \code{etat} is returned.
#'Depending on the element that can be a numeric vector, an object of
#'class "zoo" or a numeric matrix.
#'
#'@examples
#'window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
#'rt <- window.zoo(SP500, end = "2002-12-31")
#'model <- fEGarch(egarch_spec(), rt, n_test = 250)
#'residuals(model)
#'
setMethod("residuals", signature(object = "fEGarch_fit"), function(object) {
  fun <- methods::selectMethod("etat", "fEGarch_fit")
  fun(object)
})

#'@export
#'@rdname sigma
#'@aliases sigma,fEGarch_forecast-method
setMethod("sigma", "fEGarch_forecast", function(object) {
  fun <- methods::selectMethod("sigt", "fEGarch_forecast")
  fun(object)
})
#'@export
#'@rdname fitted
#'@aliases fitted,fEGarch_forecast-method
setMethod("fitted", "fEGarch_forecast", function(object) {
  fun <- methods::selectMethod("cmeans", "fEGarch_forecast")
  fun(object)
})
