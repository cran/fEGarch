#'Log-Return Calculation From Closing Prices
#'
#'Makes log-returns available from an input closing price series.
#'
#'@param close a closing price series as a numeric vector or some time series object
#'like \code{"ts"} or \code{"zoo"} ordered chronologically.
#'
#'@details
#'Let \eqn{P_t}, \eqn{t=1^,\dots,n}, be an observed closing price series. The function
#'returns
#'\deqn{r_t = \ln{P_t}-\ln{P_{t-1}}},  t = 2,\dots,n.
#'
#'@export
#'
#'@return
#'Returns the log-return series following the input \code{close}. The output object
#'has one observation less than \code{close}, but keeps potential time series
#'formatting.
#'
#'@examples
#'# Assume SP500 + 100 was a closing price series,
#'# which it is not
#'close <- SP500 + 100
#'close_to_lreturn(close)
#'

close_to_lreturn <- function(close) {
  diff(log(close))
}
