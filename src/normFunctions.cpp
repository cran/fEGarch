#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Author: Dominik Schulz
// Version: November 08, 2024

// [[Rcpp::export]]
arma::vec pdf_norm(const arma::vec& rt, const arma::vec& mu, const arma::vec& sigt) {
  const arma::vec one_div_sigt = 1.0 / sigt;
  const arma::vec eta_t2 = arma::pow((rt - mu) % one_div_sigt, 2.0);
  const arma::vec f1 = one_div_sigt / std::sqrt(2.0 * M_PI);
  const arma::vec f2 = -0.5 * eta_t2;
  return f1 % arma::exp(f2);
}

// [[Rcpp::export]]
arma::vec pdf_norm_v1(const arma::vec& rt) {
  const arma::vec eta_t2 = arma::pow(rt, 2.0);
  const double f1 = std::sqrt(2 * M_PI);
  const arma::vec f2 = -0.5 * eta_t2;
  return arma::exp(f2) / f1;
}

// [[Rcpp::export]]
arma::vec pdf_skew_snorm(const arma::vec& x, const double skew) {
  const double inv_skew = 1.0 / skew;
  const arma::uvec Hnz = (-x >= 0);
  const arma::uvec Hz = 1 - Hnz;
  const double f1 = 2.0 / (skew + inv_skew);
  const arma::vec f2 = Hnz % pdf_norm_v1(skew * x);
  const arma::vec f3 = Hz % pdf_norm_v1(inv_skew * x);
  return f1 * (f2 + f3);
}
