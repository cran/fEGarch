# ## Absolute moments t-distr.
# absMomentStd <- function(order, df) {
#
#   force(order)
#
#   FUN <- function(x, df) {
#     x^order * c(pdf_std_v1(x, df))
#   }
#
#   2 * integrate(FUN, lower = 0, upper = Inf, df = df,
#                 subdivisions = 1000L, abs.tol = 1e-6, rel.tol = 1e-6,
#                 stop.on.error = FALSE)[[1]]
#
# }

## Central moments of skewed standardized t
central_mom_sstd <- function(skew, df) {
  M1 <- 2 * sqrt(df - 2) * gamma((df + 1) / 2) / ((df - 1) * sqrt(pi) * gamma(df / 2))
  M2 <- 1
  Ez <- M1 * (skew - 1 / skew)
  skew2 <- skew^2
  M1_2 <- M1^2
  Varz <- (M2 - M1_2) * (skew2 + 1 / skew2) + 2 * M1_2 - M2
  c("E" = Ez, "Var" = Varz)
}

# pdf of standardized skewed t-distr.
pdf_skew_sstd_s <- function(x, df, skew) {
  moms <- central_mom_sstd(skew, df)
  mom1 <- moms[[1]]
  mom2 <- sqrt(moms[[2]])
  fac <- mom2
  fac * c(pdf_skew_sstd(x * fac + mom1, df = df, skew = skew))
}

# pdf of standardized skewed t-distr. with adj. for mean and cond. SD
pdf_skew_sstd_final <- function(rt, mu, sigt, df, skew) {
  pdf_skew_sstd_s((rt - mu) / sigt, df = df, skew = skew) / sigt
}
