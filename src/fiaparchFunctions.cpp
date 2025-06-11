#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List fiaparch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::vec& coef_inf,
                          const double gamma, const double delta,
                          const int n_out,
                          const double mu) {

  const int n = innov.size();

  double sigd;
  arma::vec sigt = arma::vec(n);

  const arma::vec coef_rev = arma::reverse(coef_inf);

  arma::vec rt = arma::zeros<arma::vec>(n);
  arma::vec eps2_transf = arma::zeros<arma::vec>(n);

  const int lc = coef_inf.size();

  const int n_coef_m1 = lc - 1;

  const double onedivdelta = 1.0 / delta;

  for (int j = 0; j < n_coef_m1; ++j) {

    sigd = omega +
      arma::as_scalar(arma::dot(coef_rev.subvec(n_coef_m1 - j, n_coef_m1), eps2_transf.subvec(0, j)));

    sigt(j) = std::pow(sigd, onedivdelta);

    rt(j) = innov(j) * sigt(j);
    eps2_transf(j) = std::pow(std::abs(rt(j)) - gamma * rt(j), delta);

  }

  for (int i = n_coef_m1; i < n; ++i) {

    sigd = omega +
      arma::as_scalar(arma::dot(coef_rev, eps2_transf.subvec(0 + i - n_coef_m1, i)));

    sigt(i) = std::pow(sigd, onedivdelta);

    rt(i) = innov(i) * sigt(i);
    eps2_transf(i) = std::pow(std::abs(rt(i)) - gamma * rt(i), delta);

  }

  const arma::vec sigt_out = sigt.subvec(n - n_out, n - 1);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt_out = rt.subvec(n - n_out, n - 1) + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}
