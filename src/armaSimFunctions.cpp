#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List arma_sim_Cpp(const arma::vec& innov, const double mu,
                            const arma::vec& ma,
                            const arma::vec& ar,
                            const int nout) {

  const int p = ar.size();
  const int q = ma.size();
  const int l = std::max(p, q);
  const int N = innov.size();
  const int nl = N + l;

  arma::vec x_dm = arma::zeros<arma::vec>(nl);

  const arma::vec ar_rev = (p > 0) ? arma::reverse(ar) : arma::vec();
  const arma::vec ma_rev = (q > 0) ? arma::reverse(ma) : arma::vec();

  const bool use_ar = (p > 0);
  const bool use_ma = (q > 0);

  arma::vec et = arma::zeros<arma::vec>(nl);
  et.subvec(l, nl - 1) = innov;

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

    x_dm(i) = ar_term + ma_term + et(i);

  }

  const arma::vec rt_mu = x_dm.subvec(nl - nout, nl - 1) + mu;

  Rcpp::List out = Rcpp::List::create(
    _["rt"] = rt_mu,
    _["cmeans"] = rt_mu - et.subvec(nl - nout, nl - 1)
  );

  return out;

}

// [[Rcpp::export]]
Rcpp::List farima_sim_Cpp(const arma::vec& innov, const double mu,
                            const arma::vec& coef_inf,
                            const int nout) {

  const int pre = 1;
  const int N = innov.size();
  const int n_pre = N + pre;

  arma::vec x_dm = arma::zeros<arma::vec>(n_pre);
  const arma::vec coef_rev = arma::reverse(coef_inf);
  const int lc = coef_inf.size();
  const int lc_pre = lc + pre;

  int jm1;
  const int lm1 = lc - 1;

  for (int j = pre; j < lc_pre; ++j) {

    jm1 = j - 1;

    x_dm(j) = arma::as_scalar(arma::dot(coef_rev.subvec(lc - j, lm1), x_dm.subvec(0, jm1))) + innov(jm1);

  }

  for (int i = lc_pre; i < n_pre; ++i) {

    jm1 = i - 1;
    x_dm(i) = arma::as_scalar(arma::dot(coef_rev, x_dm.subvec(i - lc, jm1))) + innov(jm1);

  }

  const arma::vec rt_mu = x_dm.subvec(n_pre - nout, n_pre - 1) + mu;
  Rcpp::List out = Rcpp::List::create(
    _["rt"] = rt_mu,
    _["cmeans"] = rt_mu - innov.subvec(N - nout, N - 1)
  );

  return out;

}
