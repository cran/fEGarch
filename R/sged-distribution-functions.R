# ## Absolute moments of GEDs.
#
# absMomentGED <- function(order, shape) {
#
#   force(order)
#
#   FUN <- function(x, shape) {
#     x^order * c(pdf_ged_v1(x, shape))
#   }
#
#   2 * integrate(FUN, lower = 0, upper = Inf, shape = shape,
#                 subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
#                 stop.on.error = FALSE)[[1]]
#
# }

## Central moments of skewed standardized GED
central_mom_sged <- function(skew, shape) {
  M1 <- gamma(2 / shape) / sqrt(gamma(1 / shape) * gamma(3 / shape))
  M2 <- 1
  Ez <- M1 * (skew - 1 / skew)
  skew2 <- skew^2
  M1_2 <- M1^2
  Varz <- (M2 - M1_2) * (skew2 + 1 / skew2) + 2 * M1_2 - M2
  c("E" = Ez, "Var" = Varz)
}

# pdf of standardized skewed GED
pdf_skew_sged_s <- function(x, shape, skew) {
  moms <- central_mom_sged(skew, shape)
  mom1 <- moms[[1]]
  mom2 <- sqrt(moms[[2]])
  fac <- mom2
  fac * c(pdf_skew_sged(x * fac + mom1, shape = shape, skew = skew))
}

# pdf of standardized skewed GED with adj. for mean and cond. SD
pdf_skew_sged_final <- function(rt, mu, sigt, shape, skew) {
  pdf_skew_sged_s((rt - mu) / sigt, shape = shape, skew = skew) / sigt
}
