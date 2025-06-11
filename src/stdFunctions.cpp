#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Author: Dominik Schulz
// Version: November 08, 2024

// [[Rcpp::export]]
arma::vec pdf_std(const arma::vec& rt, const arma::vec& mu, const arma::vec& sigt,
                  const double df) {
  const arma::vec one_div_sigt = 1.0 / sigt;
  const arma::vec eta_t2 = arma::pow((rt - mu) % one_div_sigt, 2.0);
  const double dfp1div2 = (df + 1.0) / 2.0;
  const double dfm2 = df - 2.0;
  const double f1 = std::tgamma(dfp1div2);
  const arma::vec f2 = one_div_sigt / (std::tgamma(df / 2.0) * std::sqrt(M_PI * dfm2));
  const arma::vec f3 = arma::pow(1.0 + eta_t2 / (dfm2), -dfp1div2);
  return f1 * (f2 % f3);
}


// [[Rcpp::export]]
arma::vec pdf_std_v1(const arma::vec& rt, const double df) {
  const arma::vec eta_t2 = arma::pow(rt, 2.0);
  const double dfp1div2 = (df + 1.0) / 2.0;
  const double dfm2 = df - 2.0;
  const double f1 = std::tgamma(dfp1div2);
  const double f2 = std::tgamma(df / 2.0) * std::sqrt(M_PI * dfm2);
  const arma::vec f3 = arma::pow(1.0 + eta_t2 / (dfm2), -dfp1div2);
  return f1 * (f3 / f2);
}

// [[Rcpp::export]]
arma::vec pdf_skew_sstd(const arma::vec& x, const double df, const double skew) {
  const double inv_skew = 1.0 / skew;
  const arma::uvec Hnz = (-x >= 0);
  const arma::uvec Hz = 1 - Hnz;
  const double f1 = 2.0 / (skew + inv_skew);
  const arma::vec f2 = Hnz % pdf_std_v1(skew * x, df);
  const arma::vec f3 = Hz % pdf_std_v1(inv_skew * x, df);
  return f1 * (f2 + f3);
}
