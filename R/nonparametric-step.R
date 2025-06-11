run_nonpar <- function(nonparspec, control_nonpar, rt, n_test, mean_after_nonpar) {

  nonpar_result <- nonpar_est(
    rt = rt,
    lm = TRUE,
    nonparspec = nonparspec,
    n_test = n_test,
    control_nonpar = control_nonpar
  )

  meanspec <- mean_spec(
    orders = c(0, 0),
    long_memo = FALSE,
    include_mean = FALSE
  )

  if (mean_after_nonpar) {
    include_mean(meanspec) <- TRUE
  }

  rt <- nonpar_result$et
  n_test <- 0

  list(nonpar_result = nonpar_result, meanspec = meanspec, rt = rt, n_test = n_test)
}
