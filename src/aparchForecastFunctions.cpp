#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec sigt_aparch_forecast_short(
                               const arma::vec& et,
                               const arma::vec& sigt,
                               const double omega,
                               const arma::vec& phi,
                               const arma::vec& beta,
                               const arma::vec& gamma,
                               const double delta,
                               const arma::vec& E_e,
                               const int horizon) {


  const int n = sigt.size();
  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int Hl = horizon + l;

  int il, il1, ilp, ilq;

  const arma::vec eta_i = et.subvec(n - p, n - 1);
  const arma::vec gamma_r = arma::reverse(gamma);

  arma::vec transf(p);

  arma::vec sigd(Hl);
  sigd.subvec(0, l - 1) = arma::pow(sigt.subvec(n - l, n - 1), delta);

  const arma::vec phi_r = arma::reverse(phi);
  const arma::vec beta_r = arma::reverse(beta);

  int pm1;

  const arma::vec E_e_r = arma::reverse(E_e);

  for (int i = 0; i < horizon; ++i) {
    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    pm1 = p - 1;

    if ((i >= 0) & (i <= pm1)) {
      transf.subvec(0, pm1 - i) = arma::pow(arma::abs(eta_i.subvec(i, pm1)) - gamma_r.subvec(0 , pm1 - i) % eta_i.subvec(i, pm1), delta);
      if (i >= 1) {
        transf.subvec(p - i, pm1) = E_e_r.subvec(p - i, pm1);
      }
    } else {
      transf = E_e_r;
    }

    sigd(il) = omega +
      arma::as_scalar(arma::dot(phi_r, (sigd.subvec(ilp, il1) % transf))) +
      arma::as_scalar(arma::dot(beta_r, sigd.subvec(ilq, il1)));

  }

  return arma::pow(sigd.subvec(l, Hl - 1), 1.0 / delta);
}
