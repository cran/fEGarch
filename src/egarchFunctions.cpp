#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
#include "signVec.h"
#include "signSingle.h"

// [[Rcpp::export]]
int signCpp(const double val) {
    return (0.0 < val) - (val < 0.0);
}

// [[Rcpp::export]]
arma::vec signVecCpp(const arma::vec& vals) {
    return arma::conv_to<arma::vec>::from(0.0 < vals) - arma::conv_to<arma::vec>::from(vals < 0.0);
}

// [[Rcpp::export]]
arma::vec sigt_egarch_shortCpp(const arma::vec& x, const double mu,
                               const double omega,
                               const arma::vec& phi,
                               const arma::vec& alpha,
                               const arma::vec& beta,
                               const double E_asy,
                               const double E_mag,
                               const arma::vec& powers,
                               const arma::vec& modulus,
                               const std::string& mode,
                               const double lnsig2_init) {

  const int n = x.size();
  const int p = phi.size();
  const int q = alpha.size();
  const int l = std::max(p, q);
  const int nl = n + l;
  arma::vec lnsig2 = arma::vec(nl).fill(lnsig2_init);
  double eta_i;
  arma::vec series1 = arma::vec(nl).fill(0);
  arma::vec series2 = arma::vec(nl).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::vec alpha_r = arma::reverse(alpha);
  const arma::vec phi_r = arma::reverse(phi);
  const arma::vec beta_r = arma::reverse(beta);

  int il, il1, ilp, ilq;

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  std::function<double(double)> compute_abs_eta_pow1;
  std::function<double(double)> compute_abs_eta_pow2;

  if (mode == "pow_pow") {
    compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
      return (signCpp(eta_i) * (std::pow(std::abs(eta_i) + mod1, pow1) - mod1) - E_asy) / pow1;
    };
    compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
      return (std::pow(std::abs(eta_i) + mod2, pow2) - mod2 - E_mag) / pow2;
    };
  } else if (mode == "log_log") {
    compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
      return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
    };
    compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
      return std::log(std::abs(eta_i) + mod2) - E_mag;
    };
  } else if (mode == "pow_log") {
    compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
      return (signCpp(eta_i) * (std::pow(std::abs(eta_i) + mod1, pow1) - mod1) - E_asy) / pow1;
    };
    compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
      return std::log(std::abs(eta_i) + mod2) - E_mag;
    };
  } else if (mode == "log_pow") {
    compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
      return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
    };
    compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
      return (std::pow(std::abs(eta_i) + mod2, pow2) - mod2 - E_mag) / pow2;
    };
  } else {
    Rcpp::stop("Invalid mode specified");
  }

  for (int i = 0; i < n; ++i) {
    il = i + l;
    il1 = il - 1;
    ilp = il - p;
    ilq = il - q;

    lnsig2(il) = omega +
      arma::as_scalar(arma::dot(alpha_r, series1.subvec(ilq, il1))) +
      arma::as_scalar(arma::dot(beta_r, series2.subvec(ilq, il1))) +
      arma::as_scalar(arma::dot(phi_r, lnsig2.subvec(ilp, il1)));

    sigt(i) = std::exp(lnsig2(il) / 2.0);
    eta_i = x_cent(i) / sigt(i);

    series1(il) = compute_abs_eta_pow1(eta_i);
    series2(il) = compute_abs_eta_pow2(eta_i);
  }

  return sigt;
}


// [[Rcpp::export]]
arma::vec sigt_egarch_longCpp(const arma::vec& x, const double mu,
                       const arma::vec& coef_inf,
                       const double kappa,
                       const double gamma,
                       const double E_asy,
                       const double E_mag,
                       const double Elnsig2,
                       const arma::vec& powers,
                       const arma::vec& modulus,
                       const std::string& mode) {

  const int n = x.size();
  const int n1 = n + 1;
  double lnsig2 = 0;
  arma::vec series_g = arma::vec(n1).fill(0);

  const arma::vec x_cent = x - mu;
  arma::vec sigt = arma::zeros<arma::vec>(n);

  const arma::vec coef_r = arma::reverse(coef_inf);
  double eta_i;

  int il;
  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  std::function<double(double)> compute_abs_eta_pow1;
  std::function<double(double)> compute_abs_eta_pow2;

  if (mode == "pow_pow") {
    compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
      return (signCpp(eta_i) * (std::pow(std::abs(eta_i) + mod1, pow1) - mod1) - E_asy) / pow1;
    };
    compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
      return (std::pow(std::abs(eta_i) + mod2, pow2) - mod2 - E_mag) / pow2;
    };
  } else if (mode == "log_log") {
    compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
      return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
    };
    compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
      return std::log(std::abs(eta_i) + mod2) - E_mag;
    };
  } else if (mode == "pow_log") {
    compute_abs_eta_pow1 = [pow1, mod1, E_asy](double eta_i) {
      return (signCpp(eta_i) * (std::pow(std::abs(eta_i) + mod1, pow1) - mod1) - E_asy) / pow1;
    };
    compute_abs_eta_pow2 = [mod2, E_mag](double eta_i) {
      return std::log(std::abs(eta_i) + mod2) - E_mag;
    };
  } else if (mode == "log_pow") {
    compute_abs_eta_pow1 = [mod1, E_asy](double eta_i) {
      return signCpp(eta_i) * std::log(std::abs(eta_i) + mod1) - E_asy;
    };
    compute_abs_eta_pow2 = [pow2, mod2, E_mag](double eta_i) {
      return (std::pow(std::abs(eta_i) + mod2, pow2) - mod2 - E_mag) / pow2;
    };
  } else {
    Rcpp::stop("Invalid mode specified");
  }

  const int coef_len = coef_inf.size();   // gives truncation length
  const int m = std::min(coef_len, n);
  const int cm1 = coef_len - 1;

  for (int j = 0; j < m; ++j) {

    il = j + 1;

    lnsig2 = Elnsig2 +
      arma::as_scalar(arma::dot(coef_r.subvec(cm1 - j, cm1), series_g.subvec(0, j)));

    sigt(j) = std::exp(lnsig2 / 2.0);
    eta_i = x_cent(j) / sigt(j);

    series_g(il) = kappa * compute_abs_eta_pow1(eta_i) + gamma * compute_abs_eta_pow2(eta_i);

  }

  if (n > coef_len) {

    for (int i = m; i < n; ++i) {

      il = i + 1;

      lnsig2 = Elnsig2 +
        arma::as_scalar(arma::dot(coef_r, series_g.subvec(i - m + 1, i)));

      sigt(i) = std::exp(lnsig2 / 2.0);
      eta_i = x_cent(i) / sigt(i);

      series_g(il) = kappa * compute_abs_eta_pow1(eta_i) + gamma * compute_abs_eta_pow2(eta_i);

    }

  }


  return sigt;

}


// [[Rcpp::export]]
Rcpp::List egarch_longmemo_type_sim(const arma::vec& innov, const double omega_sig,
                          const arma::vec& coef_inf,
                          const double kappa,
                          const double gamma, const int n_out,
                          const double mu, const double E_mag,
                          const double E_asy,
                          const arma::vec& powers,
                          const arma::vec& modulus,
                          const std::string& mode,
                          const int np2) {

  const int n = innov.size();
  const int n_coef = coef_inf.size();

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  arma::vec gfun = arma::zeros<arma::vec>(np2);
  arma::vec asy_term = arma::vec(n);
  arma::vec mag_term = arma::vec(n);

  if (mode == "pow_pow") {
    asy_term = (signVecCpp(innov) % (arma::pow(arma::abs(innov) + mod1, pow1) - mod1) - E_asy) / pow1;
    mag_term = (arma::pow(arma::abs(innov) + mod2, pow2) - mod2 - E_mag) / pow2;
  } else if (mode == "log_log") {
    asy_term = signVecCpp(innov) % arma::log(arma::abs(innov) + mod1) - E_asy;
    mag_term = arma::log(arma::abs(innov) + mod2) - E_mag;
  } else if (mode == "pow_log") {
    asy_term = (signVecCpp(innov) % (arma::pow(arma::abs(innov) + mod1, pow1) - mod1) - E_asy) / pow1;
    mag_term = arma::log(arma::abs(innov) + mod2) - E_mag;
  } else if (mode == "log_pow") {
    asy_term = signVecCpp(innov) % arma::log(arma::abs(innov) + mod1) - E_asy;
    mag_term = (arma::pow(arma::abs(innov) + mod2, pow2) - mod2 - E_mag) / pow2;
  } else {
    Rcpp::stop("Invalid mode specified");
  }

  gfun.subvec(0, n - 1) = kappa * asy_term + gamma * mag_term;
  arma::vec coef_adj = arma::zeros<arma::vec>(np2);
  coef_adj.subvec(0, n_coef - 1) = coef_inf;
  const arma::vec out = arma::real(arma::ifft(arma::fft(coef_adj) % arma::fft(gfun)));
  const arma::vec sigt = arma::exp((omega_sig + out.subvec(n - n_out, n - 1)) / 2.0); // no division by np2 required in C++; is implemented by ifft automatically
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt = mu + sigt % innov_out;

  return Rcpp::List::create(_["rt"] = rt, _["sigt"] = sigt, _["etat"] = innov_out);
}

// [[Rcpp::export]]
Rcpp::List egarch_shortmemo_type_sim(const arma::vec& innov, const double omega_sig,
                          const arma::vec& ar_coef, const arma::vec& ma_coef,
                          const double kappa,
                          const double gamma, const int n_out,
                          const double mu, const double E_mag,
                          const double E_asy,
                          const arma::vec& powers,
                          const arma::vec& modulus,
                          const std::string& mode) {

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

  const double mod1 = modulus(0);
  const double mod2 = modulus(1);
  const double pow1 = powers.size() > 0 ? powers(0) : 1.0; // Default to 1.0 for safety
  const double pow2 = powers.size() > 1 ? powers(1) : 1.0;

  arma::vec gfun = arma::zeros<arma::vec>(n_all);
  arma::vec asy_term = arma::vec(n);
  arma::vec mag_term = arma::vec(n);

  if (mode == "pow_pow") {
    asy_term = (signVecCpp(innov) % (arma::pow(arma::abs(innov) + mod1, pow1) - mod1) - E_asy) / pow1;
    mag_term = (arma::pow(arma::abs(innov) + mod2, pow2) - mod2 - E_mag) / pow2;
  } else if (mode == "log_log") {
    asy_term = signVecCpp(innov) % arma::log(arma::abs(innov) + mod1) - E_asy;
    mag_term = arma::log(arma::abs(innov) + mod2) - E_mag;
  } else if (mode == "pow_log") {
    asy_term = (signVecCpp(innov) % (arma::pow(arma::abs(innov) + mod1, pow1) - mod1) - E_asy) / pow1;
    mag_term = arma::log(arma::abs(innov) + mod2) - E_mag;
  } else if (mode == "log_pow") {
    asy_term = signVecCpp(innov) % arma::log(arma::abs(innov) + mod1) - E_asy;
    mag_term = (arma::pow(arma::abs(innov) + mod2, pow2) - mod2 - E_mag) / pow2;
  } else {
    Rcpp::stop("Invalid mode specified");
  }

  gfun.subvec(l, n_all - 1) = kappa * asy_term + gamma * mag_term;
  const double omega = omega_sig * (1 - arma::accu(ar_coef));
  const int lar = l - n_ar;
  const int lma = l - n_ma;
  const int lm1 = l - 1;
  int li;

  for (int i = 0; i < n; ++i) {

    li = lm1 + i;

    lnsig2(l + i) = omega +
      arma::as_scalar(arma::dot(ar_rev, lnsig2.subvec(lar + i, li))) +
      arma::as_scalar(arma::dot(ma_rev, gfun.subvec(lma + i, li)));

  }

  const arma::vec sigt = arma::exp(lnsig2.subvec(n_all - n_out, n_all - 1) / 2.0);
  const arma::vec innov_out = innov.subvec(n - n_out, n - 1);

  const arma::vec rt = mu + sigt % innov_out;

  return Rcpp::List::create(_["rt"] = rt, _["sigt"] = sigt, _["etat"] = innov_out);
}

