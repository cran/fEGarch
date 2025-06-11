#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec fitted_arma_Cpp(const arma::vec& x, const double mu,
                            const arma::vec& ma,
                            const arma::vec& ar) {

  const int p = ar.size();
  const int q = ma.size();
  const int l = std::max(p, q);
  const int n = x.size();
  const int nl = n + l;

  arma::vec x_dm(nl);
  if (mu == 0) {
    x_dm = arma::vec(nl).fill(arma::mean(x));
  } else {
    x_dm = arma::zeros<arma::vec>(nl);
  }
  arma::vec fitted = arma::zeros<arma::vec>(nl);
  x_dm.subvec(l, nl - 1) = x - mu;

  const arma::vec ar_rev = (p > 0) ? arma::reverse(ar) : arma::vec();
  const arma::vec ma_rev = (q > 0) ? arma::reverse(ma) : arma::vec();

  const bool use_ar = (p > 0);
  const bool use_ma = (q > 0);

  arma::vec et = arma::zeros<arma::vec>(nl);

  double ar_term = 0.0;
  double ma_term = 0.0;

  int im1;

  for (int i = l; i < nl; ++i) {

    im1 = i - 1;

    if (use_ar) {
      ar_term = arma::as_scalar(arma::dot(ar_rev, x_dm.subvec(i - p, im1)));
    }
    if (use_ma) {
      ma_term = arma::as_scalar(arma::dot(ma_rev, et.subvec(i - q, im1)));
    }

    fitted(i) = ar_term + ma_term;

    et(i) = x_dm(i) - fitted(i);

  }

  return fitted.subvec(l, nl - 1) + mu;

}

