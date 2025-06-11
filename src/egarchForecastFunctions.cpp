#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "signVec.h"
#include "signSingle.h"

// [[Rcpp::export]]
arma::vec sigt_egarch_forecast_shortCpp(
                               const arma::vec& et,
                               const arma::vec& sigt,
                               const double omega,
                               const arma::vec& phi,
                               const arma::vec& alpha,
                               const arma::vec& beta,
                               const double E_asy,
                               const double E_mag,
                               const arma::vec& powers,
                               const arma::vec& modulus,
                               const std::string& mode,
                               const int horizon) {

  const int n = sigt.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int Hl = horizon + l;

  int il, il1, ilp, ilq;

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  arma::vec transf1(l);
  arma::vec transf2(l);
  const arma::vec eta_i = et.subvec(n - l, n - 1);

  if (mode == "pow_pow") {
    transf1 = (signVecCpp(eta_i) % (arma::pow(arma::abs(eta_i) + mod1, pow1) - mod1) - E_asy) / pow1;
    transf2 = (arma::pow(arma::abs(eta_i) + mod2, pow2) - mod2 - E_mag) / pow2;
  } else if (mode == "log_log") {
    transf1 = signVecCpp(eta_i) % arma::log(arma::abs(eta_i) + mod1) - E_asy;
    transf2 = arma::log(arma::abs(eta_i) + mod2) - E_mag;
  } else if (mode == "pow_log") {
    transf1 = (signVecCpp(eta_i) % (arma::pow(arma::abs(eta_i) + mod1, pow1) - mod1) - E_asy) / pow1;
    transf2 = arma::log(arma::abs(eta_i) + mod2) - E_mag;
  } else if (mode == "log_pow") {
    transf1 = signVecCpp(eta_i) % arma::log(arma::abs(eta_i) + mod1) - E_asy;
    transf2 = (arma::pow(arma::abs(eta_i) + mod2, pow2) - mod2 - E_mag) / pow2;
  } else {
    Rcpp::stop("Invalid mode specified");
  }

  arma::vec et_asy = arma::zeros<arma::vec>(Hl);
  et_asy.subvec(0, l - 1) = transf1;

  arma::vec et_mag = arma::zeros<arma::vec>(Hl);
  et_mag.subvec(0, l - 1) = transf2;

  arma::vec lnsig2 = arma::zeros<arma::vec>(Hl);
  lnsig2.subvec(0, l - 1) = arma::log(arma::pow(sigt.subvec(n - l, n - 1), 2.0));

  const arma::vec alpha_r = arma::reverse(alpha);
  const arma::vec phi_r = arma::reverse(phi);
  const arma::vec beta_r = arma::reverse(beta);

  for (int i = 0; i < horizon; ++i) {
    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;

    lnsig2(il) = omega +
      arma::as_scalar(arma::dot(alpha_r, et_asy.subvec(ilq, il1))) +
      arma::as_scalar(arma::dot(beta_r, et_mag.subvec(ilq, il1))) +
      arma::as_scalar(arma::dot(phi_r, lnsig2.subvec(ilp, il1)));

  }

  return arma::exp(0.5 * lnsig2.subvec(l, Hl - 1));
}


// [[Rcpp::export]]
arma::vec sigt_egarch_forecast_longCpp(const arma::vec& et,
                       const arma::vec& coef_inf,
                       const double kappa,
                       const double gamma,
                       const double E_asy,
                       const double E_mag,
                       const double Elnsig2,
                       const arma::vec& powers,
                       const arma::vec& modulus,
                       const std::string& mode,
                       const int horizon) {

  const int n = et.size();
  arma::vec lnsig2 = arma::zeros<arma::vec>(horizon);

  const arma::vec coef_r = arma::reverse(coef_inf);

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  arma::vec et_asy(n);
  arma::vec et_mag(n);

  if (mode == "pow_pow") {
    et_asy = (signVecCpp(et) % (arma::pow(arma::abs(et) + mod1, pow1) - mod1) - E_asy) / pow1;
    et_mag = (arma::pow(arma::abs(et) + mod2, pow2) - mod2 - E_mag) / pow2;
  } else if (mode == "log_log") {
    et_asy = signVecCpp(et) % arma::log(arma::abs(et) + mod1) - E_asy;
    et_mag = arma::log(arma::abs(et) + mod2) - E_mag;
  } else if (mode == "pow_log") {
    et_asy = (signVecCpp(et) % (arma::pow(arma::abs(et) + mod1, pow1) - mod1) - E_asy) / pow1;
    et_mag = arma::log(arma::abs(et) + mod2) - E_mag;
  } else if (mode == "log_pow") {
    et_asy = signVecCpp(et) % arma::log(arma::abs(et) + mod1) - E_asy;
    et_mag = (arma::pow(arma::abs(et) + mod2, pow2) - mod2 - E_mag) / pow2;
  } else {
    Rcpp::stop("Invalid mode specified");
  }

  arma::vec series_g = arma::zeros<arma::vec>(n + horizon);
  series_g.subvec(0, n - 1) = kappa * et_asy + gamma * et_mag;

  const int nm1 = n - 1;
  const int nc = coef_r.size();
  const int ncm1 = nc - 1;
  int npi;

  for (int i = 0; i < horizon; ++i) {

    npi = n + i;

    if (npi < nc) {

      lnsig2(i) = Elnsig2 +
        arma::as_scalar(arma::dot(coef_r.subvec(nc - npi, ncm1), series_g.subvec(0, nm1 + i)));

    } else {

      lnsig2(i) = Elnsig2 +
        arma::as_scalar(arma::dot(coef_r, series_g.subvec(npi - nc, npi - 1)));

    }

  }

  return arma::exp(0.5 * lnsig2);

}

