test_that("Time series formatting of results works for 'zoo' and 'ts' (and stays neutral otherwise)", {
  expect_s3_class({
    func <- get("zoo", envir = asNamespace("zoo"))
    tp <- as.Date(c("2020-01-01", "2020-01-02", "2020-01-03", "2020-01-04"))
    x <- func(c(1, 2, 3, 4, 5), order.by = tp)
    test_l <- list(
      c(7, 3, 2, 9, 12),
      c(1, 5, 2, 10, 1)
    )
    out <- format_applier_ts(rt = x, list_of_ts = test_l)
    out[[1]]
  }, "zoo")
  expect_s3_class({
    func <- get("zoo", envir = asNamespace("zoo"))
    tp <- as.Date(c("2020-01-01", "2020-01-02", "2020-01-03", "2020-01-04"))
    x <- func(c(1, 2, 3, 4, 5), order.by = tp)
    test_l <- list(
      c(7, 3, 2, 9, 12),
      c(1, 5, 2, 10, 1)
    )
    out <- format_applier_ts(rt = x, list_of_ts = test_l)
    out[[2]]
  }, "zoo")
  expect_s3_class({
    func <- get("ts", envir = asNamespace("stats"))
    x <- func(c(1, 2, 3, 4, 5), start = c(2000, 1), frequency = 4)
    test_l <- list(
      c(7, 3, 2, 9, 12),
      c(1, 5, 2, 10, 1)
    )
    out <- format_applier_ts(rt = x, list_of_ts = test_l)
    out[[1]]
  }, "ts")
  expect_s3_class({
    func <- get("ts", envir = asNamespace("stats"))
    x <- func(c(1, 2, 3, 4, 5), start = c(2000, 1), frequency = 4)
    test_l <- list(
      c(7, 3, 2, 9, 12),
      c(1, 5, 2, 10, 1)
    )
    out <- format_applier_ts(rt = x, list_of_ts = test_l)
    out[[2]]
  }, "ts")
  expect_identical({
    x <- c(1, 2, 3, 4, 5)
    test_l <- list(
      c(7, 3, 2, 9, 12),
      c(1, 5, 2, 10, 1)
    )
    out <- format_applier_ts(rt = x, list_of_ts = test_l)
    out
  }, list(
      c(7, 3, 2, 9, 12),
      c(1, 5, 2, 10, 1)
    ))
})
