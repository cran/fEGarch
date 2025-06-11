#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Multistep point forecasts using an ARMA (using Box-Cox idea)

// [[Rcpp::export]]
arma::vec forecast_arma_Cpp(const arma::vec& x, const arma::vec& et,
                            const double mu,
                            const arma::vec& ma_e,
                            const arma::vec& ar_e,
                            const int horizon) {

  const int p = ar_e.size();
  const int q = ma_e.size();
  const int l = std::max(p, q);
  const int n = x.size();
  const int Hl = horizon + l;

  arma::vec x_dm_fc = arma::zeros<arma::vec>(Hl);
  x_dm_fc.subvec(0, l - 1) = x.subvec(n - l, n - 1) - mu;
  arma::vec et_s = arma::zeros<arma::vec>(Hl);
  et_s.subvec(0, l - 1) = et.subvec(n - l, n - 1);

  const arma::vec ar_rev = (p > 0) ? arma::reverse(ar_e) : arma::vec();
  const arma::vec ma_rev = (q > 0) ? arma::reverse(ma_e) : arma::vec();

  const bool use_ar = (p > 0);
  const bool use_ma = (q > 0);

  double ar_term = 0.0;
  double ma_term = 0.0;

  int im1;

  for (int i = l; i < Hl; ++i) {

    im1 = i - 1;

    if (use_ar) {
      ar_term = arma::as_scalar(arma::dot(ar_rev, x_dm_fc.subvec(i - p, im1)));
    }
    if (use_ma) {
      ma_term = arma::as_scalar(arma::dot(ma_rev, et_s.subvec(i - q, im1)));
    }

    x_dm_fc(i) = ar_term + ma_term;

  }

  return x_dm_fc.subvec(l, Hl - 1) + mu;

}

// Multistep point forecasts using a FARIMA (approx. via AR-representation)

// [[Rcpp::export]]
arma::vec forecast_farima_Cpp(const arma::vec& x, const double mu,
                            const arma::vec& coef_inf,
                            const int horizon) {

  const int n = x.size();
  const int nH = n + horizon;

  arma::vec x_dm_fc = arma::zeros<arma::vec>(nH);
  x_dm_fc.subvec(0, n - 1) = x - mu;
  const arma::vec coef_rev = arma::reverse(coef_inf);

  const int nc = coef_inf.size();

  const int ncm1 = nc - 1;
  const int ncmn = nc - n;

  int ipn;

  for (int i = 0; i < horizon; ++i) {

  ipn = i + n;

    // Depending on coef. truncation length and series length, make a  distinction
    if (ipn < nc) {

      x_dm_fc(ipn) = arma::as_scalar(arma::dot(coef_rev.subvec(ncmn - i, ncm1), x_dm_fc.subvec(0, ipn - 1)));

    } else {

      x_dm_fc(ipn) = arma::as_scalar(arma::dot(coef_rev, x_dm_fc.subvec(ipn - nc, ipn - 1)));

    }

  }

  return x_dm_fc.subvec(n, nH - 1) + mu;

}


