#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Author: Dominik Schulz
// Version: March 16, 2025

// Compute coefficients of inverse fractional differencing operator

// [[Rcpp::export]]
arma::vec d_coefs_inv(const double d, const int max_i) {

  arma::vec d_out(max_i + 1);
  d_out(0) = 1.0;
  const arma::vec k = arma::regspace<arma::vec>(1, max_i);
  const arma::vec k_inv = 1.0 / k;
  d_out.subvec(1, max_i) = arma::cumprod((k - 1.0 + d) % k_inv);

  return d_out;
}

// [[Rcpp::export]]
arma::vec d_coefs(const double d, const int max_i) {

  arma::vec d_out(max_i + 1);
  d_out(0) = 1.0;
  const arma::vec k = arma::regspace<arma::vec>(1, max_i);
  const arma::vec k_inv = 1.0 / k;
  d_out.subvec(1, max_i) = arma::cumprod((k - 1.0 - d) % k_inv);

  return d_out;
}

// Compute the inversion of a polynomial of the
// form 1 - ar_1 * B - ar_2 * B^2 - ... - ar_p * B^p
// with result in the form
// c_0 + c_1 * B + c_2 * B^2 + ... with c_0 = 1.

// [[Rcpp::export]]
arma::vec poly_inv(const arma::vec& ar, const int max_i) {

  const int maxip1 = max_i + 1;
  arma::vec ma_out(maxip1);

  const int p = ar.size();
  if (p <= 0) {
    ma_out.fill(0.0);
    ma_out(0) = 1.0;
    return ma_out;
  }

  ma_out(0) = 1.0;

  const arma::vec ar_r = arma::reverse(ar);

  for (int i = 1; i < p; ++i) {
    ma_out(i) = arma::as_scalar(arma::dot(ma_out.subvec(0, i - 1), ar_r.subvec(p - i, p - 1)));
  }

  for (int j = p; j <= max_i; ++j) {
    ma_out(j) = arma::as_scalar(arma::dot(ma_out.subvec(j - p, j - 1), ar_r));
  }

  return ma_out;
}
