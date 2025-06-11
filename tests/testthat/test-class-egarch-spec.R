test_that("base_garch_spec is validated correctly", {
  expect_no_error(
    base_garch_spec(
      orders = c(1, 1),
      long_memo = TRUE,
      cond_dist = "norm"
    )
  )
  expect_no_error(
    base_garch_spec()
  )
  expect_error(
    base_garch_spec(orders = c(0, 1)),
    "@orders must be of length two and both elements must be at least one"
  )
  expect_error(
    base_garch_spec(orders = c(0, 0)),
    "@orders must be of length two and both elements must be at least one"
  )
  expect_error(
    base_garch_spec(orders = c(1, 0)),
    "@orders must be of length two and both elements must be at least one"
  )
  expect_error(
    base_garch_spec(orders = numeric(0)),
    "@orders must be of length two and both elements must be at least one"
  )
  expect_error(
    base_garch_spec(orders = c(1, 1, 1)),
    "@orders must be of length two and both elements must be at least one"
  )
  expect_no_error(
    base_garch_spec(orders = c(2, 2))
  )
  expect_no_error(
    base_garch_spec(orders = c(2, 1))
  )
  expect_no_error(
    base_garch_spec(orders = c(1, 2))
  )
  expect_error(base_garch_spec(orders = c("a", "b")), "invalid object for slot \"orders\" in class \"base_garch_spec\": got class \"character\", should be or extend class \"numeric\"")

  expect_no_error(base_garch_spec(long_memo = TRUE))
  expect_no_error(base_garch_spec(long_memo = FALSE))
  expect_error(base_garch_spec(long_memo = logical(0)), "@long_memo must be of length one")
  expect_error(base_garch_spec(long_memo = c(TRUE, FALSE)), "@long_memo must be of length one")
  expect_error(base_garch_spec(long_memo = "a"), "invalid object for slot \"long_memo\" in class \"base_garch_spec\": got class \"character\", should be or extend class \"logical\"")

  expect_no_error(base_garch_spec(cond_dist = "norm"))
  expect_error(base_garch_spec(cond_dist = c("norm", "std")), "\"cond_dist\" must be kept at default or must be of length 1 and one of the defaults")
  expect_error(base_garch_spec(cond_dist = character(0)), "\"cond_dist\" must be kept at default or must be of length 1 and one of the defaults")
  expect_error(base_garch_spec(cond_dist = "smd"), "\"cond_dist\" must be kept at default or must be of length 1 and one of the defaults")

})
