
## Central moments of skewed standardized average Laplace distribution
central_mom_sald <- function(skew, P) {
  M1 <- ald_first_abs_mom(P = P)
  M2 <- 1
  Ez <- M1 * (skew - 1 / skew)
  skew2 <- skew^2
  M1_2 <- M1^2
  Varz <- (M2 - M1_2) * (skew2 + 1 / skew2) + 2 * M1_2 - M2
  c("E" = Ez, "Var" = Varz)
}

# pdf of standardized skewed average Laplace distribution
pdf_skew_sald_s <- function(x, P, skew) {
  moms <- central_mom_sald(skew, P)
  mom1 <- moms[[1]]
  mom2 <- sqrt(moms[[2]])
  fac <- mom2
  fac * c(pdf_skew_sald(x * fac + mom1, P = P, skew = skew))
}

# pdf of standardized skewed average Laplace distribution with adj. for mean and cond. SD
pdf_skew_sald_final <- function(rt, mu, sigt, P, skew) {
  pdf_skew_sald_s((rt - mu) / sigt, P = P, skew = skew) / sigt
}
