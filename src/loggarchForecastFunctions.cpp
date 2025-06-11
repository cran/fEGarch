#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec sigt_loggarch_forecast_short(
                               const arma::vec& et,
                               const arma::vec& sigt,
                               const double omega,
                               const arma::vec& phi,
                               const arma::vec& psi,
                               const double Elneta2,
                               const int horizon) {

  const int n = sigt.size();
  const int p = phi.size();
  const int q = psi.size();
  const int l = std::max(p, q);
  const int Hl = horizon + l;

  int il, il1, ilp, ilq;

  const arma::vec eta_i = et.subvec(n - l, n - 1);

  const arma::vec transf = arma::log(arma::pow(eta_i, 2.0)) - Elneta2;

  arma::vec et_fun = arma::zeros<arma::vec>(Hl);
  et_fun.subvec(0, l - 1) = transf;

  arma::vec lnsig2 = arma::zeros<arma::vec>(Hl);
  lnsig2.subvec(0, l - 1) = arma::log(arma::pow(sigt.subvec(n - l, n - 1), 2.0));

  const arma::vec psi_r = arma::reverse(psi);
  const arma::vec phi_r = arma::reverse(phi);

  for (int i = 0; i < horizon; ++i) {
    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;

    lnsig2(il) = omega +
      arma::as_scalar(arma::dot(psi_r, et_fun.subvec(ilq, il1))) +
      arma::as_scalar(arma::dot(phi_r, lnsig2.subvec(ilp, il1)));

  }

  return arma::exp(0.5 * lnsig2.subvec(l, Hl - 1));
}

// [[Rcpp::export]]
arma::vec sigt_loggarch_forecast_long(const arma::vec& et,
                       const arma::vec& coef_inf,
                       const double Elneta2,
                       const double Elnsig2,
                       const int horizon) {

  const int n = et.size();
  arma::vec lnsig2 = arma::vec(horizon);
  arma::vec g_fun = arma::zeros<arma::vec>(n + horizon);
  g_fun.subvec(0, n - 1) = arma::log(arma::pow(et, 2.0)) - Elneta2;

  const arma::vec coef_inf_r = arma::reverse(coef_inf);

  const int nm1 = n - 1;
  const int nc = coef_inf_r.size();
  const int ncm1 = nc - 1;
  int npi;

  for (int i = 0; i < horizon; ++i) {

    npi = n + i;

    if (npi < nc) {

      lnsig2(i) = Elnsig2 +
        arma::as_scalar(arma::dot(coef_inf_r.subvec(nc - npi, ncm1), g_fun.subvec(0, nm1 + i)));

    } else {

      lnsig2(i) = Elnsig2 +
        arma::as_scalar(arma::dot(coef_inf_r, g_fun.subvec(npi - nc, npi - 1)));

    }

  }

  return arma::exp(0.5 * lnsig2);

}
