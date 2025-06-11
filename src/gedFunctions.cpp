#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Author: Dominik Schulz
// Version: November 08, 2024

// [[Rcpp::export]]
arma::vec pdf_ged(const arma::vec& rt, const arma::vec& mu, const arma::vec& sigt,
                  const double shape) {
  const arma::vec one_div_sigt = 1.0 / sigt;
  const arma::vec abs_eta_t = arma::abs((rt - mu) % one_div_sigt);
  const double g3b = std::tgamma(3.0 / shape);
  const double g1b = std::tgamma(1.0 / shape);
  const arma::vec f1 = one_div_sigt * std::sqrt(g3b / std::pow(g1b, 3.0));
  const double f2 = shape / 2.0;
  const double f3s = std::pow(g3b / g1b, f2);
  const arma::vec f3 = -arma::pow(abs_eta_t, shape) * f3s;
  return f2 * (f1 % arma::exp(f3));
}


// [[Rcpp::export]]
arma::vec pdf_ged_v1(const arma::vec& rt, const double shape) {
  const arma::vec abs_eta_t = arma::abs(rt);
  const double g3b = std::tgamma(3.0 / shape);
  const double g1b = std::tgamma(1.0 / shape);
  const double f1 = std::sqrt(g3b / std::pow(g1b, 3.0));
  const double f2 = shape / 2.0;
  const double f3s = std::pow(g3b / g1b, f2);
  const arma::vec f3 = -arma::pow(abs_eta_t, shape) * f3s;
  return f2 * (f1 * arma::exp(f3));
}

// [[Rcpp::export]]
arma::vec pdf_skew_sged(const arma::vec& x, const double shape, const double skew) {
  const double inv_skew = 1.0 / skew;
  const arma::uvec Hnz = (-x >= 0);
  const arma::uvec Hz = 1 - Hnz;
  const double f1 = 2.0 / (skew + inv_skew);
  const arma::vec f2 = Hnz % pdf_ged_v1(skew * x, shape);
  const arma::vec f3 = Hz % pdf_ged_v1(inv_skew * x, shape);
  return f1 * (f2 + f3);
}
