// simulate_count.cpp
// Family generators for the count / multistate families (N-mixture, removal,
// open N-mixture, false-positive occupancy). Each draws ONE replicate at the
// parameter vector `theta` (the fit's draw layout) plus each arm's realised
// offset, the latent abundance / state and the observations coming from R's
// RNG stream via the R:: samplers in a fixed order. simulate_<family>(),
// simulate(fit) and the SBC replicate all call these, so the generative
// process is written once.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "simulate_helpers.h"
using namespace Rcpp;
using tulpaObs::arm_eta;
using tulpaObs::draw_latent_N;

// N-mixture: latent N per site, then a structural-zero draw per site when
// `zi_omega` is finite (zero-inflated mixtures; the plain mixtures draw
// nothing here), then the binomial detections in the long-form observation
// order (site_idx, visit_idx, 1-based).
// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_simulate_nmix(
    Rcpp::NumericMatrix X_lambda, Rcpp::NumericMatrix X_p, Rcpp::NumericVector theta,
    Rcpp::NumericVector off_lambda, Rcpp::NumericVector off_p,
    Rcpp::IntegerVector site_idx, Rcpp::IntegerVector visit_idx,
    int n_sites, int max_visits, int p_lam, int p_p,
    bool is_nb, double r_size, double zi_omega = NA_REAL
) {
  const int n_obs = site_idx.size();
  const bool zi = R_finite(zi_omega);
  Rcpp::RNGScope scope;
  const double* pXl = X_lambda.begin(); const double* pXp = X_p.begin();
  const double* th = theta.begin();
  std::vector<int> N(n_sites);
  for (int i = 0; i < n_sites; ++i)
    N[i] = draw_latent_N(std::exp(arm_eta(pXl, n_sites, i, th, 0, p_lam, off_lambda)),
                         is_nb ? r_size : R_PosInf);
  std::vector<int> zero(n_sites, 0);
  if (zi) {
    for (int i = 0; i < n_sites; ++i)
      if ((zero[i] = (int) R::rbinom(1.0, zi_omega)) == 1) N[i] = 0;
  }
  Rcpp::IntegerMatrix ys(n_sites, max_visits);
  std::fill(ys.begin(), ys.end(), NA_INTEGER);
  for (int k = 0; k < n_obs; ++k) {
    double po = tulpa::math::inv_logit(arm_eta(pXp, n_obs, k, th, p_lam, p_p, off_p));
    ys(site_idx[k] - 1, visit_idx[k] - 1) = (int) R::rbinom((double) N[site_idx[k] - 1], po);
  }
  ys.attr("latent") = Rcpp::List::create(Rcpp::Named("N") = Rcpp::wrap(N),
                                         Rcpp::Named("zero") = Rcpp::wrap(zero));
  return ys;
}

// Removal: latent N per site, then each site's passes in turn, every pass
// removing a binomial share of the individuals not yet removed.
// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_simulate_removal(
    Rcpp::NumericMatrix X_lambda, Rcpp::NumericMatrix X_p, Rcpp::NumericVector theta,
    Rcpp::NumericVector off_lambda, Rcpp::NumericVector off_p,
    Rcpp::IntegerVector site_idx, Rcpp::IntegerVector visit_idx,
    int n_sites, int n_pass, int p_lam, int p_p,
    bool is_nb, double r_size
) {
  const int n_obs = site_idx.size();
  Rcpp::RNGScope scope;
  const double* pXl = X_lambda.begin(); const double* pXp = X_p.begin();
  const double* th = theta.begin();
  std::vector<int> N(n_sites);
  for (int i = 0; i < n_sites; ++i)
    N[i] = draw_latent_N(std::exp(arm_eta(pXl, n_sites, i, th, 0, p_lam, off_lambda)),
                         is_nb ? r_size : R_PosInf);
  // Per-site-pass detection prob (long form scattered to [n_sites x n_pass]).
  std::vector<double> pmat((std::size_t) n_sites * n_pass, NA_REAL);
  for (int k = 0; k < n_obs; ++k)
    pmat[(std::size_t) (visit_idx[k] - 1) * n_sites + (site_idx[k] - 1)] =
      tulpa::math::inv_logit(arm_eta(pXp, n_obs, k, th, p_lam, p_p, off_p));
  Rcpp::IntegerMatrix ys(n_sites, n_pass);
  for (int i = 0; i < n_sites; ++i) {
    int rem = N[i];
    for (int k = 0; k < n_pass; ++k) {
      int yk = (int) R::rbinom((double) rem, pmat[(std::size_t) k * n_sites + i]);
      ys(i, k) = yk; rem -= yk;
    }
  }
  ys.attr("latent") = Rcpp::wrap(N);
  return ys;
}

// Open N-mixture (dyn_abun): per site, initial abundance then a survival +
// recruitment HMM across seasons, observed by binomial detection. The lambda
// and p arms are site-level; omega and gamma are INTERVAL-level -- the
// transition into season t uses interval (t-1)'s rate -- so their designs (and
// offsets) carry either n_sites rows (constant rates) or n_sites * (n_seasons-1)
// rows in the site-major layout .tobs_period_arm_design() builds. Both shapes
// and the row index are read exactly as the likelihood reads them
// (src/dyn_abun_laplace.cpp: `eta_omega[i * nIv + iv]` against `eta_omega[i]`).
// A zero-inflated mixture (`zi` > 0) first draws each site's structural-zero
// indicator; such a site stays at zero and draws nothing further. RNG order is
// then site-major (draw N, then per season the transition draws + the J
// detections), and the interval lookup draws nothing, so a constant-rate design
// in either shape draws the same replicate.
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_simulate_dyn_abun(
    Rcpp::NumericMatrix X_lambda, Rcpp::NumericMatrix X_p,
    Rcpp::NumericMatrix X_omega, Rcpp::NumericMatrix X_gamma,
    Rcpp::NumericVector theta,
    Rcpp::NumericVector off_lambda, Rcpp::NumericVector off_p,
    Rcpp::NumericVector off_omega, Rcpp::NumericVector off_gamma,
    int n_sites, int n_seasons, int J,
    int p_lam, int p_p, int p_om, int p_gm,
    bool is_nb, double r_disp, double zi = 0.0
) {
  const int o_p = p_lam, o_om = p_lam + p_p, o_gm = p_lam + p_p + p_om;
  // The two transition designs are accepted at the two shapes the likelihood
  // accepts and refused at any other, rather than indexed on an assumption:
  // read at the per-site stride, a season-varying design's non-intercept
  // columns land in the wrong memory region and the interval variation is
  // silently dropped.
  const int nIv = n_seasons - 1;
  const int nrow_om = X_omega.nrow(), nrow_gm = X_gamma.nrow();
  const bool om_iv = (nIv > 0 && nrow_om == n_sites * nIv);
  const bool gm_iv = (nIv > 0 && nrow_gm == n_sites * nIv);
  if (!om_iv && nrow_om != n_sites)
    Rcpp::stop("X_omega must have n_sites or n_sites*(n_seasons-1) rows.");
  if (!gm_iv && nrow_gm != n_sites)
    Rcpp::stop("X_gamma must have n_sites or n_sites*(n_seasons-1) rows.");
  Rcpp::RNGScope scope;
  Rcpp::IntegerVector out((std::size_t) n_sites * J * n_seasons);
  const double* pXl = X_lambda.begin(); const double* pXp = X_p.begin();
  const double* pXo = X_omega.begin(); const double* pXg = X_gamma.begin();
  const double* th = theta.begin();
  std::vector<int> zero(n_sites, 0);
  Rcpp::IntegerMatrix Nmat(n_sites, n_seasons);
  if (zi > 0.0) {
    for (int i = 0; i < n_sites; ++i) zero[i] = (int) R::rbinom(1.0, zi);
  }
  for (int i = 0; i < n_sites; ++i) {
    if (zero[i] == 1) continue;
    double lambda = std::exp(arm_eta(pXl, n_sites, i, th, 0, p_lam, off_lambda));
    double pdet = tulpa::math::inv_logit(arm_eta(pXp, n_sites, i, th, o_p, p_p, off_p));
    int N = draw_latent_N(lambda, is_nb ? r_disp : R_PosInf);
    for (int t = 0; t < n_seasons; ++t) {
      if (t > 0) {
        const int iv = t - 1;
        double omega = tulpa::math::inv_logit(arm_eta(
            pXo, nrow_om, om_iv ? i * nIv + iv : i, th, o_om, p_om, off_omega));
        double gamma = std::exp(arm_eta(
            pXg, nrow_gm, gm_iv ? i * nIv + iv : i, th, o_gm, p_gm, off_gamma));
        N = (int) R::rbinom((double) N, omega) + (int) R::rpois(gamma);
      }
      Nmat(i, t) = N;
      for (int j = 0; j < J; ++j)
        out[(std::size_t) i + (std::size_t) j * n_sites + (std::size_t) t * n_sites * J] =
          (int) R::rbinom((double) N, pdet);
    }
  }
  out.attr("dim") = Rcpp::IntegerVector::create(n_sites, J, n_seasons);
  out.attr("latent") = Rcpp::List::create(Rcpp::Named("N") = Nmat,
                                          Rcpp::Named("zero") = Rcpp::wrap(zero));
  return out;
}

// False-positive occupancy: the site state, then per site either the J
// detections and their J certainty flags (occupied: 2 certain, 1 uncertain,
// 0 missed) or the J false-positive detections (unoccupied).
// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_simulate_fp_occu(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p11, Rcpp::NumericMatrix X_p10,
    Rcpp::NumericMatrix X_b, Rcpp::NumericVector theta,
    Rcpp::NumericVector off_psi, Rcpp::NumericVector off_p11,
    Rcpp::NumericVector off_p10, Rcpp::NumericVector off_b,
    int n_sites, int J, int p_psi, int p_p11, int p_p10, int p_b
) {
  const int o_p11 = p_psi, o_p10 = p_psi + p_p11, o_b = p_psi + p_p11 + p_p10;
  Rcpp::RNGScope scope;
  const double* th = theta.begin();
  std::vector<double> psi(n_sites), p11(n_sites), p10(n_sites), b(n_sites);
  for (int i = 0; i < n_sites; ++i) {
    psi[i] = tulpa::math::inv_logit(arm_eta(X_psi.begin(), n_sites, i, th, 0, p_psi, off_psi));
    p11[i] = tulpa::math::inv_logit(arm_eta(X_p11.begin(), n_sites, i, th, o_p11, p_p11, off_p11));
    p10[i] = tulpa::math::inv_logit(arm_eta(X_p10.begin(), n_sites, i, th, o_p10, p_p10, off_p10));
    b[i]   = tulpa::math::inv_logit(arm_eta(X_b.begin(),   n_sites, i, th, o_b,   p_b,   off_b));
  }
  std::vector<int> z(n_sites);
  for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, psi[i]);
  Rcpp::IntegerMatrix ys(n_sites, J);
  for (int i = 0; i < n_sites; ++i) {
    if (z[i] == 1) {
      for (int j = 0; j < J; ++j) ys(i, j) = (int) R::rbinom(1.0, p11[i]);
      for (int j = 0; j < J; ++j) {
        int cert = (int) R::rbinom(1.0, b[i]);
        ys(i, j) = (ys(i, j) == 1) ? (cert == 1 ? 2 : 1) : 0;
      }
    } else {
      for (int j = 0; j < J; ++j) ys(i, j) = (int) R::rbinom(1.0, p10[i]);
    }
  }
  ys.attr("latent") = Rcpp::wrap(z);
  return ys;
}
