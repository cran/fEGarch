#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_fiaparch_forecast(const arma::vec& rt,
                       const arma::vec& coef_inf,
                       const double omega,
                       const double gamma,
                       const double delta,
                       const double E_const,
                       const int horizon) {

  const int n = rt.size();
  arma::vec sig_delta = arma::vec(horizon);
  arma::vec g_fun = arma::zeros<arma::vec>(n + horizon);
  g_fun.subvec(0, n - 1) = arma::pow(arma::abs(rt) - gamma * rt, delta);

  const arma::vec coef_inf_r = arma::reverse(coef_inf);

  const int nm1 = n - 1;
  const int nc = coef_inf_r.size();
  const int ncm1 = nc - 1;
  int npi;

  for (int i = 0; i < horizon; ++i) {

    npi = n + i;

    if (npi < nc) {

      sig_delta(i) = omega +
        arma::as_scalar(arma::dot(coef_inf_r.subvec(nc - npi, ncm1), g_fun.subvec(0, nm1 + i)));

    } else {

      sig_delta(i) = omega +
        arma::as_scalar(arma::dot(coef_inf_r, g_fun.subvec(npi - nc, npi - 1)));

    }

    g_fun(n + i) = sig_delta(i) * E_const;

  }

  return arma::pow(sig_delta, 1.0 / delta);

}
