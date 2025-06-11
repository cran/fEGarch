#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec sigt_loggarch_long(const arma::vec& x, const double mu,
                       const arma::vec& coef_inf,
                       const double Elneta2,
                       const double Elnsig2) {

  const int n = x.size();
  const int n1 = n + 1;
  arma::vec lnsig2 = arma::vec(n1).fill(0);
  arma::vec epst = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::vec coef_inf_r = arma::reverse(coef_inf);
  double eta_i;

  int il;

  const int coef_len = coef_inf.size();
  const int cm1 = coef_len - 1;
  const int m = std::min(coef_len, n);

  for (int j = 0; j < m; ++j) {

    il = j + 1;

    lnsig2(il) = Elnsig2 +
      arma::as_scalar(arma::dot(coef_inf_r.subvec(cm1 - j, cm1), epst.subvec(0, j)));

    sigt(j) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(j) / sigt(j);
    epst(il) = std::log(std::pow(eta_i, 2.0)) - Elneta2;

  }

  if (n > coef_len) {

    for (int i = m; i < n; ++i) {

      il = i + 1;

      lnsig2(il) = Elnsig2 +
        arma::as_scalar(arma::dot(coef_inf_r, epst.subvec(i - m + 1, i)));

      sigt(i) = std::exp(lnsig2(il) / 2.0);
      eta_i = x_cent(i) / sigt(i);
      epst(il) = std::log(std::pow(eta_i, 2.0)) - Elneta2;

    }

  }



  return sigt;

}

// [[Rcpp::export]]
arma::vec sigt_loggarch_short(const arma::vec& x, const double mu,
                       const double omega,
                       const arma::vec& phi,
                       const arma::vec& psi,
                       const double Elneta2,
                       const double lnsig2_init) {

  const int n = x.size();
  const int p = phi.size();
  const int q = psi.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(lnsig2_init);
  double eta_i;
  arma::vec series_innov = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::vec psi_r = arma::reverse(psi);
  const arma::vec phi_r = arma::reverse(phi);

  int il;
  int il1;
  int ilp;
  int ilq;

  for (int i = 0; i < n; ++i) {

    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;
    lnsig2(il) = omega +
      arma::as_scalar(arma::dot(psi_r, series_innov.subvec(ilq, il1))) +
      arma::as_scalar(arma::dot(phi_r, lnsig2.subvec(ilp, il1)));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);
    series_innov(il) = std::log(std::pow(eta_i, 2.0)) - Elneta2;

  }

  return sigt;

}

// [[Rcpp::export]]
Rcpp::List loggarch_shortmemo_type_sim(const arma::vec& innov, const double omega_sig,
                          const arma::vec& ar_coef, const arma::vec& ma_coef,
                          const int n_out,
                          const double mu,
                          const double Econst) {

  const int n = innov.size();

  const int n_ar = ar_coef.size();
  const int n_ma = ma_coef.size();
  const int l = std::max(n_ar, n_ma);

  const int n_pre = l;
  const int n_all = n + n_pre;

  arma::vec lnsig2(n_all);
  lnsig2.fill(omega_sig);

  const arma::vec ar_rev = arma::reverse(ar_coef);
  const arma::vec ma_rev = arma::reverse(ma_coef);

  arma::vec xi = arma::zeros<arma::vec>(n_all);
  xi.subvec(l, n_all - 1) = arma::log(arma::pow(innov, 2.0)) - Econst;

  const double omega = omega_sig * (1 - arma::accu(ar_coef));
  const int lar = l - n_ar;
  const int lma = l - n_ma;
  const int lm1 = l - 1;
  int li;

  for (int i = 0; i < n; ++i) {

    li = lm1 + i;

    lnsig2(l + i) = omega +
      arma::as_scalar(arma::dot(ar_rev, lnsig2.subvec(lar + i, li))) +
      arma::as_scalar(arma::dot(ma_rev, xi.subvec(lma + i, li)));

  }

  const arma::vec sigt = arma::exp(lnsig2.subvec(n_all - n_out, n_all - 1) / 2.0);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt = mu + sigt % innov_out;

  return Rcpp::List::create(_["rt"] = rt, _["sigt"] = sigt, _["etat"] = innov_out);
}

// [[Rcpp::export]]
Rcpp::List loggarch_longmemo_type_sim(const arma::vec& innov, const double omega_sig,
                          const arma::vec& coef_inf,
                          const int n_out,
                          const double mu,
                          const double Econst,
                          const int np2) {

  const int n = innov.size();
  const int l = coef_inf.size();

  const arma::vec i_transf = arma::log(arma::pow(innov, 2.0)) - Econst;

  arma::vec e_adj = arma::zeros<arma::vec>(np2);
  arma::vec c_adj = arma::zeros<arma::vec>(np2);

  e_adj.subvec(0, n - 1) = i_transf;
  c_adj.subvec(0, l - 1) = coef_inf;

  const arma::vec interm = arma::real(arma::ifft(arma::fft(e_adj) % arma::fft(c_adj)));
  const arma::vec sigt = arma::exp((omega_sig + interm.subvec(n - n_out, n - 1)) / 2.0); // no division by np2 required in C++; is implemented by ifft automatically

  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt = mu + sigt % innov_out;

  return Rcpp::List::create(_["rt"] = rt, _["sigt"] = sigt, _["etat"] = innov_out);
}
