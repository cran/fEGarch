test_that("Standard normal distribution function is correctly specified", {
  expect_equal(
    integrate(function(x) {x * pdf_norm_v1(x)}, lower = -Inf, upper = Inf)[[1]],
    0
  )
  expect_equal(
    integrate(function(x) {x^2 * pdf_norm_v1(x)}, lower = -Inf, upper = Inf)[[1]],
    1
  )
  expect_equal(
    integrate(function(x) {x^3 * pdf_norm_v1(x)}, lower = -Inf, upper = Inf)[[1]],
    0
  )
  expect_equal(
    pdf_norm_v1(c(-1.96, -0.5, 0, 0.5, 1.96)),
    matrix(c(0.05844094, 0.35206533, 0.39894228, 0.35206533, 0.05844094), ncol = 1)
  )
  expect_equal(
    c(pdf_norm_v1(c(-1.96, -0.5, 0, 0.5, 1.96))),
    c(0.05844094, 0.35206533, 0.39894228, 0.35206533, 0.05844094)
  )
  expect_equal(
    pdf_norm_v1(-1.96),
    pdf_norm_v1(1.96)
  )
  expect_equal(
    pdf_norm_v1(-1),
    pdf_norm_v1(1)
  )
  expect_equal(
    pdf_norm_v1(-0.5),
    pdf_norm_v1(0.5)
  )
})

#------------------------------------------------------

test_that("Check likelihood function for normal distr.", {
  expect_equal(
  c(pdf_norm(c(-1.96, -1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5, 1.96),
      mu = rep(0.05, 11), sigt = rep(0.5, 11))),
c(0.000247032326682047, 0.00653363811239984, 0.0879671919608544,
0.435704354065101, 0.666449205783599, 0.793905094954024, 0.736540280606647,
0.53217049979751, 0.131231629549353, 0.0119050648395517, 0.000541054062923042
)
)
})

