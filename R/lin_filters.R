# Find the next power of 2 equal to or above x
# (FFT works best on vectors of length "power of 2")

nextpow2 <- function(x) 2^ceiling(log2(x))

# Compute the infinite-order MA-representation
# from a given FARIMA using FFT; "max_i" gives the truncation
# (see also Nielsen and Noel, 2021)

ma_infty <- function(ar, ma, d, max_i) {

  n <- max_i + 1
  np2 <- nextpow2(2 * n - 1)
  ma <- c(1, ma)
  l <- length(ma)
  ma_adj <- c(ma, rep(0, np2 - l))
  d_adj <- c(d_coefs_inv(d, max_i), rep(0, np2 - n))

  ma_inf <- Re(stats::fft(stats::fft(d_adj) * stats::fft(ma_adj), inverse = TRUE) / np2)

  ar_adj <- c(poly_inv(ar, max_i = max_i), rep(0, np2 - n))

  utils::head(Re(stats::fft(stats::fft(ar_adj) * stats::fft(ma_inf), inverse = TRUE) / np2), n)

}

ar_infty <- function(ar, ma, d, max_i) {

  n <- max_i + 1
  np2 <- nextpow2(2 * n - 1)
  ar <- c(1, -ar)
  l <- length(ar)
  ar_adj <- c(ar, rep(0, np2 - l))

  d_adj <- c(d_coefs(d, max_i), rep(0, np2 - n))

  ar_inf <- Re(stats::fft(stats::fft(d_adj) * stats::fft(ar_adj), inverse = TRUE) / np2)

  ma_adj <- c(poly_inv(-ma, max_i = max_i), rep(0, np2 - n))

  -utils::head(Re(stats::fft(stats::fft(ma_adj) * stats::fft(ar_inf), inverse = TRUE) / np2), n)

}
