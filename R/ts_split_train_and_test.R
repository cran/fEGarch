split_ts <- function(xt, n_test = 0) {
  if (n_test == 0) {
    return(list(train = xt, test = NULL))
  } else {
    if (inherits(xt, "zoo")) {
      tp <- time(xt)
      tp_train <- as.Date(utils::head(tp, -n_test))
      tp_test <- as.Date(utils::tail(tp, n_test))
      xt_core <- zoo::coredata(xt)
      xt_train <- zoo::zoo(
        utils::head(xt_core, -n_test),
        order.by = tp_train
      )
      xt_test <- zoo::zoo(
        utils::tail(xt_core, n_test),
        order.by = tp_test
      )
    } else if (inherits(xt, "ts")) {
      tp <- time(xt)
      start_train <- stats::start(xt)
      start_test <- utils::tail(tp, n_test)[[1]]
      frequ <- stats::frequency(xt)
      xt_core <- zoo::coredata(xt)
      xt_train <- stats::ts(
        utils::head(xt_core, -n_test),
        start = start_train, frequency = frequ
      )
      xt_test <- stats::ts(
        utils::tail(xt_core, n_test),
        start = start_test, frequency = frequ
      )
    } else {
      xt_train <- utils::head(xt, -n_test)
      xt_test <- utils::tail(xt, n_test)
    }
    out <- list(
      train = xt_train,
      test = xt_test
    )
  }
  out
}
