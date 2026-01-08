#'Generics for Model Specification Accessors
#'
#'A collection of generics for accessors of the model specification
#'objects.
#'
#'@param x the input object or object to modify.
#'@param value the value to modify the object \code{x} with.
#'
#'@details
#'These generics are intended to provide a basis to construct methods
#'for adjusting the the output of model specification objects like
#'\code{egarch_spec()}.
#'
#'@name spec_generics
#'@aliases orders orders,ANY-method
#'
#'@return
#'These generics themselves do not return anything and are just the
#'foundation for more sophisticated methods.
#'
#'@export
#'@rdname spec_generics
setGeneric("orders", function(x) {standardGeneric("orders")})
#'@export
#'@rdname spec_generics
#'@aliases powers powers,ANY-method
setGeneric("powers", function(x) {standardGeneric("powers")})
#'@export
#'@rdname spec_generics
#'@aliases long_memo long_memo,ANY-method
setGeneric("long_memo", function(x) {standardGeneric("long_memo")})
#'@export
#'@rdname spec_generics
#'@aliases modulus modulus,ANY-method
setGeneric("modulus", function(x) {standardGeneric("modulus")})
#'@export
#'@rdname spec_generics
#'@aliases cond_dist cond_dist,ANY-method
setGeneric("cond_dist", function(x) {standardGeneric("cond_dist")})
#'@export
#'@rdname spec_generics
#'@aliases include_mean include_mean,ANY-method
setGeneric("include_mean", function(x) {standardGeneric("include_mean")})

#'@export
#'@rdname spec_generics
setGeneric("orders<-", function(x, value) {standardGeneric("orders<-")})
#'@export
#'@rdname spec_generics
setGeneric("powers<-", function(x, value) {standardGeneric("powers<-")})
#'@export
#'@rdname spec_generics
setGeneric("long_memo<-", function(x, value) {standardGeneric("long_memo<-")})
#'@export
#'@rdname spec_generics
setGeneric("modulus<-", function(x, value) {standardGeneric("modulus<-")})
#'@export
#'@rdname spec_generics
setGeneric("cond_dist<-", function(x, value) {standardGeneric("cond_dist<-")})
#'@export
#'@rdname spec_generics
setGeneric("include_mean<-", function(x, value) {standardGeneric("include_mean<-")})


#====================================================#

#Generic for Fitting EGARCH Family Models
#
#The generic is currently without use.
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
#The generic is currently without use. Nothing is returned.
#
#
#

setGeneric("fEGarch_fit", function(spec, rt, drange = c(0, 1), meanspec = mean_spec(), Drange = c(0, 1), n_test = 0, start_pars = NULL, LB = NULL, UB = NULL, control = list(), parallel = TRUE, ncores = max(1, future::availableCores() - 1), trunc = floor(0.4 * length(rt)), presample = 50, Prange = c(1, 5), skip_vcov = FALSE) {standardGeneric("fEGarch_fit")})

#====================================================#

#@rdname forecasting-generics
setGeneric("fEGarch_predict", function(object, n.ahead = 10, trunc = NULL, ...) {standardGeneric("fEGarch_predict")})

#====================================================#

#@rdname forecasting-generics
setGeneric("predict_internal", function(object, n.ahead = 10, trunc = NULL, ...) {standardGeneric("predict_internal")})
setGeneric("predict_roll_internal", function(object, step_size = 1, trunc = NULL, ...) {standardGeneric("predict_roll_internal")})


#====================================================#

#'Generics for Forecasts
#'
#'The generics are themselves without use.
#'
#'@param object the generics are currently without use.
#'@param step_size the generics are currently without use.
#'@param trunc the generics are currently without use.
#'@param refit_after the generics are currently without use.
#'@param steady_window the generics are currently without use.
#'@param parallel the generics are currently without use.
#'@param ncores the generics are currently without use.
#'@param fitting_args the generics are currently without use.
#'@param ... the generics are currently without use.
#'
#'@return
#'The generics do not work themselves and therefore
#'do not return anything.
#'
#'@export
#'
#'@rdname forecasting-generics
#'
setGeneric("predict_roll", function(object, step_size = 1, trunc = NULL, refit_after = NULL, steady_window = FALSE, parallel = TRUE, ncores = max(1, future::availableCores() - 1), fitting_args = list(), ...) {standardGeneric("predict_roll")})

#====================================================#
#'Generics for Accessing Model Estimation Output Elements
#'
#'Accessors for the output returned by the main fitting functions of
#'the \code{fEGarch} package. The generics themselves are without use.
#'
#'@param x an object returned by either \code{\link{fEGarch}},
#'\code{\link{fiaparch}} or \code{\link{figarch}}, etc.
#'
#'@details
#'These generics are without direct use. Consider the
#'specific methods based on them.
#'
#'@name fitted_object_generics
#'@aliases sigt sigt,ANY-method
#'
#'@return
#'These generics do not return anything. Inspect the methods based on
#'them for practical purpose.
#'
#'@export
setGeneric("sigt", function(x) {standardGeneric("sigt")})
#'@export
#'@rdname fitted_object_generics
#'@aliases cmeans cmeans,ANY-method
setGeneric("cmeans", function(x) {standardGeneric("cmeans")})
#'@rdname fitted_object_generics
#'@export
#'@aliases etat etat,ANY-method
setGeneric("etat", function(x) {standardGeneric("etat")})
#'@rdname fitted_object_generics
#'@export
#'@aliases llhood llhood,ANY-method
setGeneric("llhood", function(x) {standardGeneric("llhood")})
#'@rdname fitted_object_generics
#'@export
#'@aliases inf_criteria inf_criteria,ANY-method
setGeneric("inf_criteria", function(x) {standardGeneric("inf_criteria")})
#'@rdname fitted_object_generics
#'@export
#'@aliases pars pars,ANY-method
setGeneric("pars", function(x) {standardGeneric("pars")})
#'@rdname fitted_object_generics
#'@export
#'@aliases se se,ANY-method
setGeneric("se", function(x) {standardGeneric("se")})
#'@rdname fitted_object_generics
#'@export
#'@aliases vcov_mat vcov_mat,ANY-method
setGeneric("vcov_mat", function(x) {standardGeneric("vcov_mat")})

#########################

#'Generics for Nonparametric Smoothing Setting Adjustments
#'
#'Generics to build methods on for adjusting certain settings
#'in context of nonparametric smoothing. The generics are currently
#'without use.
#'
#'@param x the generics are currently without use.
#'@param value the generics are currently without use.
#'
#'@export
#'
#'@return
#'The generics do not work on their own and thus do not return anything.
#'
#'@rdname nonpar-generics
#'
setGeneric("poly_order", function(x) {standardGeneric("poly_order")})
#'@rdname nonpar-generics
#'@export
setGeneric("kernel_order", function(x) {standardGeneric("kernel_order")})
#'@rdname nonpar-generics
#'@export
setGeneric("boundary_method", function(x) {standardGeneric("boundary_method")})
#'@rdname nonpar-generics
#'@export
setGeneric("bwidth", function(x) {standardGeneric("bwidth")})


#'@rdname nonpar-generics
#'@export
setGeneric("poly_order<-", function(x, value) {standardGeneric("poly_order<-")})
#'@rdname nonpar-generics
#'@export
setGeneric("kernel_order<-", function(x, value) {standardGeneric("kernel_order<-")})
#'@rdname nonpar-generics
#'@export
setGeneric("boundary_method<-", function(x, value) {standardGeneric("boundary_method<-")})
#'@rdname nonpar-generics
#'@export
setGeneric("bwidth<-", function(x, value) {standardGeneric("bwidth<-")})

