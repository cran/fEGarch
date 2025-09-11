window.zoo <- get("window.zoo", envir = asNamespace("zoo"))
rt <- window.zoo(SP500, end = "2002-12-31")
spec <- egarch_spec()
spec2 <- fiegarch_spec()

test_that("Bandwidth selections are employed for short and long memory correctly", {
  expect_equal(
   round(fEGarch(spec, rt, use_nonpar = TRUE)@nonpar_model$b0, 4), 0.1896
  )
  expect_equal(
   round(fEGarch(spec2, rt, use_nonpar = TRUE)@nonpar_model$b0, 4), 0.1894
  )


})
