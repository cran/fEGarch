window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")

test_that("garchm_estim function works as intended", {

  expect_no_error({
   garchm_estim(rt, model = "garch")
  })

  expect_no_error({
   suppressWarnings(garchm_estim(rt, model = "aparch"))
  })

  expect_no_error({
   garchm_estim(rt, model = "gjrgarch")
  })

  expect_no_error({
   garchm_estim(rt, model = "figarch")
  })

  expect_no_error({
   suppressWarnings(garchm_estim(rt, model = "fiaparch"))
  })

  expect_no_error({
   garchm_estim(rt, model = "figjrgarch")
  })

  expect_no_error({
   garchm_estim(rt, model = "figarch", trunc = 50)
  })

  expect_no_error({
   suppressWarnings(garchm_estim(rt, model = "fiaparch", trunc = 50))
  })

  expect_no_error({
   garchm_estim(rt, model = "figjrgarch", trunc = 50)
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", cond_dist = "norm")
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", cond_dist = "std")
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", cond_dist = "ged")
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", cond_dist = "ald", parallel = FALSE, Prange = c(1, 3))
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", cond_dist = "sald", parallel = FALSE, Prange = c(1, 3))
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", cond_dist = "snorm")
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", cond_dist = "sstd")
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", cond_dist = "sged")
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", orders = c(1, 1))
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", orders = c(2, 1))
  })

  expect_no_error({
   garchm_estim(rt, model = "garch", orders = c(1, 2))
  })

})
