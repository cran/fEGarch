# ## Absolute moments normal distr.
# absMomentNorm <- function(order) {
#
#   force(order)
#
#   FUN <- function(x) {
#     x^order * c(pdf_norm_v1(x))
#   }
#
#   2 * integrate(FUN, lower = 0, upper = Inf,
#                 subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
#                 stop.on.error = FALSE)[[1]]
# }

## Central moments of skewed standard normal
central_mom_snorm <- function(skew) {
  M1 <- sqrt(2 / pi)
  M2 <- 1
  Ez <- M1 * (skew - 1 / skew)
  skew2 <- skew^2
  M1_2 <- M1^2
  Varz <- (M2 - M1_2) * (skew2 + 1 / skew2) + 2 * M1_2 - M2
  c("E" = Ez, "Var" = Varz)
}

# pdf of standardized skewed normal
pdf_skew_norm_s <- function(x, skew) {
  moms <- central_mom_snorm(skew)
  mom1 <- moms[[1]]
  mom2 <- sqrt(moms[[2]])
  fac <- mom2
  fac * c(pdf_skew_snorm(x * fac + mom1, skew = skew))
}

# pdf of standardized skewed normal with adjustments for mean and cond. SD
pdf_skew_norm_final <- function(rt, mu, sigt, skew) {
  pdf_skew_norm_s((rt - mu) / sigt, skew = skew) / sigt
}
