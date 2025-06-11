window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")
spec <- egarch_spec()

test_that("fEGarch function works as intended", {
  expect_no_error({
   fEGarch(spec, rt)
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(1, 1)))
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(0, 1)))
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(1, 0)))
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(0, 0), long_memo = TRUE))
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(0, 0), long_memo = TRUE, include_mean = FALSE))
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(0, 0), include_mean = FALSE))
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(1, 1), include_mean = FALSE))
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(1, 0), include_mean = FALSE))
  })

  expect_no_error({
   fEGarch(spec, rt, meanspec = mean_spec(orders = c(0, 1), include_mean = FALSE))
  })

  expect_no_error({
   fEGarch(spec, rt, nonparspec = locpol_spec(), use_nonpar = TRUE)
  })

  expect_no_error({
   fEGarch(spec, rt, nonparspec = locpol_spec(p = 1), use_nonpar = TRUE)
  })

  expect_no_error({
   fEGarch(spec, rt, nonparspec = locpol_spec(p = 3), use_nonpar = TRUE)
  })

  expect_no_error({
   fEGarch(spec, rt, nonparspec = locpol_spec(kernel_order = 2), use_nonpar = TRUE)
  })

  expect_no_error({
   fEGarch(spec, rt, nonparspec = locpol_spec(boundary_method = "shorten"), use_nonpar = TRUE)
  })

  expect_no_error({
   fEGarch(spec, rt, nonparspec = locpol_spec(bwidth = 0.14), use_nonpar = TRUE)
  })

  expect_no_error({
   suppressWarnings(fEGarch(fiegarch_spec(), rt, trunc = 100))
  })

  expect_no_error({
   suppressWarnings(fEGarch(fiegarch_spec(), rt, trunc = "none"))
  })

  expect_no_error({
   fEGarch(egarch_spec(cond_dist = "ald"), rt, parallel = FALSE, Prange = c(1, 3))
  })

  expect_no_error({
   fEGarch(egarch_spec(cond_dist = "std"), rt)
  })

  expect_no_error({
   fEGarch(egarch_spec(cond_dist = "ged"), rt)
  })

  expect_no_error({
   fEGarch(egarch_spec(cond_dist = "snorm"), rt)
  })

  expect_no_error({
   fEGarch(egarch_spec(cond_dist = "sstd"), rt)
  })

  expect_no_error({
   fEGarch(egarch_spec(cond_dist = "sged"), rt)
  })

})
