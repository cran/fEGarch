#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_garch_Cpp(const arma::vec& x, const double mu,
                       const arma::vec& phi, const arma::vec& beta,
                       const double omega, const double sig2_init) {

  const int n = x.size();
  const int p = phi.size();
  const int q = beta.size();
  const int l = std::max(p, q);
  const int nl = n + l;

  const arma::vec x_dm = x - mu;
  const arma::vec Y2 = arma::pow(x_dm, 2.0);

  arma::vec Y_fun = arma::vec(nl).fill(sig2_init);
  arma::vec sig_2 = arma::vec(nl).fill(sig2_init);
  Y_fun.subvec(l, nl - 1) = Y2;

  const arma::vec phi_r = arma::reverse(phi);
  const arma::vec beta_r = arma::reverse(beta);

  for (int i = l; i < nl; ++i) {

    sig_2(i) = omega +
      arma::as_scalar(arma::dot(phi_r, Y_fun.subvec(i - p, i - 1))) +
      arma::as_scalar(arma::dot(beta_r, sig_2.subvec(i - q, i - 1)));

  }

  return arma::sqrt(sig_2.subvec(l, nl - 1));

}

