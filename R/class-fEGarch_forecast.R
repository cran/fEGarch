setClass("fEGarch_forecast",
  slots = c(
    cmeans = "ANY",
    sigt = "ANY",
    model = "fEGarch_fit"
  ),
  prototype = list(
    cmeans = numeric(1),
    sigt = numeric(1),
    model = new("fEGarch_fit")
  )
)


fEGarch_forecast <- function(
    cmeans,
    sigt,
    model
) {

  new("fEGarch_forecast",
    cmeans = cmeans,
    sigt = sigt,
    model = model
  )

}

#'@export
#'@rdname accessor_methods
#'@aliases sigt,fEGarch_forecast-method
setMethod("sigt", "fEGarch_forecast", function(x) {x@sigt})
#'@export
#'@rdname accessor_methods
#'@aliases cmeans,fEGarch_forecast-method
setMethod("cmeans", "fEGarch_forecast", function(x) {x@cmeans})
