format_applier_ts <- function(rt, list_of_ts) {

  if (inherits(rt, what = "zoo")) {
    time_rt <-  time(rt)
    tp <- as.Date(time_rt)
    out <- lapply(list_of_ts, FUN = function(.x, tp) {
      zoo::zoo(.x, order.by = tp)
    }, tp = tp)
  } else if (inherits(rt, what = "ts")) {
    start_rt <- stats::start(rt)
    frequ_rt <- stats::frequency(rt)
    out <- lapply(list_of_ts, FUN = function(.x, start_rt, frequ_rt) {
      stats::ts(.x, start = start_rt, frequency = frequ_rt)
    }, start_rt = start_rt, frequ_rt = frequ_rt)
  } else {
    out <- list_of_ts
  }
  out
}
