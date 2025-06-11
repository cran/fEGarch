#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_aparch_Cpp(const arma::vec& x, const double mu,
                       const arma::vec& phi, const arma::vec& beta,
                       const arma::vec& gamma,
                       const double delta,
                       const double omega,
                       const double sigd_init,
                       const double etransf_init) {

  const int n = x.size();
  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int nl = n + l;

  const arma::vec x_dm = x - mu;

  arma::vec r_transf = arma::vec(p).fill(etransf_init);
  arma::vec sig_delta = arma::vec(nl).fill(sigd_init);

  const arma::vec phi_r = arma::reverse(phi);
  const arma::vec gamma_r = arma::reverse(gamma);
  const arma::vec beta_r = arma::reverse(beta);

  const int pm1 = p - 1;

  for (int j = l; j <= l + pm1; ++j) {

    if (j > l) {
      r_transf.subvec(pm1 + l + 1 - j, pm1) = arma::pow(arma::abs(x_dm.subvec(0, j - 1 - l)) - gamma_r.subvec(p - j + l, pm1) % x_dm.subvec(0, j - 1 - l), delta);
    }

    sig_delta(j) = omega +
      arma::as_scalar(arma::dot(phi_r, r_transf)) +
      arma::as_scalar(arma::dot(beta_r, sig_delta.subvec(j - q, j - 1)));

  }

  const int lp = l + p;

  for (int i = lp; i < nl; ++i) {

    sig_delta(i) = omega +
      arma::as_scalar(arma::dot(phi_r, arma::pow(arma::abs(x_dm.subvec(i - lp, i - lp + pm1)) - gamma_r % x_dm.subvec(i - lp, i - lp + pm1), delta))) +
      arma::as_scalar(arma::dot(beta_r, sig_delta.subvec(i - q, i - 1)));

  }

  return arma::pow(sig_delta.subvec(l, nl - 1), 1.0 / delta);

}

// [[Rcpp::export]]
Rcpp::List aparch_sim_Cpp(const arma::vec& innov, const double omega,
                          const arma::rowvec& phi, const arma::rowvec& beta,
                          const arma::vec& gamma, const double delta,
                          const int n_out,
                          const double mu,
                          const double E_sigd) {

  const int n = innov.size();

  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);

  const int n_pre = l;
  const int n_all = n + n_pre;

  arma::vec sig_d(n_all);
  sig_d.fill(E_sigd);

  // Initialize innovs with zeros
  arma::vec innov_a = arma::zeros<arma::vec>(n_all);
  innov_a.subvec(l, n_all - 1) = innov;

  const arma::rowvec phi_rev = arma::reverse(phi);
  const arma::rowvec beta_rev = arma::reverse(beta);
  const arma::vec gamma_r = arma::reverse(gamma);


  int imp;
  int im1;
  int imq;

  arma::vec gfun(p);
  arma::vec innov_sel(p);

  for (int i = l; i < n_all; ++i) {

    imp = i - p;
    imq = i - q;
    im1 = i - 1;

    innov_sel = innov_a.subvec(imp, im1);
    gfun = arma::pow(arma::abs(innov_sel) - gamma_r % innov_sel, delta);

    sig_d(i) = omega +
      arma::as_scalar(arma::dot(phi_rev, (sig_d.subvec(imp, im1) % gfun))) +
      arma::as_scalar(arma::dot(beta_rev, sig_d.subvec(imq, im1)));

  }

  const arma::vec sigt_out = arma::pow(sig_d.subvec(n_all - n_out, n_all - 1), 1.0 / delta);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt_out = (sigt_out % innov_out) + mu;

  return Rcpp::List::create(_["rt"] = rt_out, _["sigt"] = sigt_out, _["etat"] = innov_out);
}
