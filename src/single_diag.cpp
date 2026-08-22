// single_diag.cpp
// C++ kernels for the single-season occupancy posterior diagnostics whose
// per-draw / per-site loops were pure R:
//   * cpp_single_ppc -- posterior predictive check (ppc(), single). Per
//     selected draw: latent z ~ Bernoulli(full conditional), detection replicate
//     y_rep ~ Bernoulli(z p), Freeman-Tukey / chi-squared discrepancy.
//   * cpp_single_pit -- randomized PIT residuals (pit_residuals()). Per site:
//     the posterior-mean predictive CDF plus a uniform jitter.
// The posterior draw SELECTION (sample.int) stays in R and its indices are
// passed in; the per-draw RNG (z, y_rep, the PIT jitter) is drawn here from R's
// stream via the R:: samplers, in the SAME order as the former R loops, so under
// a fixed seed the results are byte-identical.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "tobs_shape.h"
using namespace Rcpp;
using tulpaObs::stable_plogis;
using tulpaObs::row_draw_dot;
using tulpaObs::ppc_stat;
namespace shape = tulpaObs::shape;


// [[Rcpp::export]]
Rcpp::List cpp_single_ppc(
    Rcpp::NumericMatrix X_occ,    // [n_sites x p_occ]
    Rcpp::NumericMatrix X_det,    // [n_sites x p_det]
    Rcpp::NumericMatrix draws,    // [ndraws x (p_occ + p_det)]
    Rcpp::IntegerVector draw_idx, // [n.samples] 1-based, from sample.int in R
    Rcpp::IntegerMatrix y,        // [n_sites x max_visits], < 0 = NA
    Rcpp::IntegerVector n_valid,  // [n_sites]
    Rcpp::IntegerVector any_det,  // [n_sites] 0/1
    bool freeman
) {
  const int n_sites = X_occ.nrow();
  const int p_occ = X_occ.ncol(), p_det = X_det.ncol();
  const int max_v = y.ncol();
  const int ndr = draws.nrow();
  const int nsamp = draw_idx.size();
  shape::check_nrow(X_det, n_sites, "X_det");
  shape::check_nrow(y, n_sites, "y");
  shape::check_ncol_min(draws, (R_xlen_t) p_occ + p_det, "draws");
  shape::check_index1(draw_idx, ndr, "draw_idx");
  shape::check_len(n_valid, n_sites, "n_valid");
  shape::check_len(any_det, n_sites, "any_det");
  Rcpp::RNGScope scope;
  Rcpp::NumericVector fit_y(nsamp), fit_rep(nsamp);
  const double* pXo = X_occ.begin(); const double* pXd = X_det.begin();
  const double* pdr = draws.begin(); const int* py = y.begin();
  auto stat = [&](double o, double e) { return ppc_stat(o, e, freeman); };
  std::vector<double> psi(n_sites), p(n_sites);
  std::vector<int> z(n_sites);

  for (int s = 0; s < nsamp; ++s) {
    int idx = draw_idx[s] - 1;
    for (int i = 0; i < n_sites; ++i) {
      psi[i] = stable_plogis(row_draw_dot(pXo, n_sites, i, pdr, ndr, idx, 0, p_occ));
      p[i]   = stable_plogis(row_draw_dot(pXd, n_sites, i, pdr, ndr, idx, p_occ, p_det));
    }
    // z_prob (deterministic), then z ~ Bernoulli in site order.
    for (int i = 0; i < n_sites; ++i) {
      double zp;
      if (n_valid[i] == 0) zp = psi[i];
      else if (any_det[i]) zp = 1.0;
      else {
        double a = psi[i] * std::pow(1.0 - p[i], (double) n_valid[i]);
        zp = a / (a + (1.0 - psi[i]));
      }
      z[i] = (int) R::rbinom(1.0, zp);
    }
    // y_rep ~ Bernoulli(z p) in (i outer, j inner) order, matching the R loop.
    double f_obs = 0.0, f_rep = 0.0;
    for (int i = 0; i < n_sites; ++i)
      for (int j = 0; j < max_v; ++j) {
        int yij = py[(std::size_t) j * n_sites + i];
        if (yij < 0) continue;
        double e = z[i] * p[i];
        int yr = (int) R::rbinom(1.0, e);
        f_obs += stat((double) yij, e);
        f_rep += stat((double) yr, e);
      }
    fit_y[s] = f_obs; fit_rep[s] = f_rep;
  }
  return Rcpp::List::create(Rcpp::Named("fit.y") = fit_y,
                            Rcpp::Named("fit.y.rep") = fit_rep);
}

// The site-level occupancy response is the ordered detected/all-zero event D_i
// = 1{>=1 detection}, not the individual visit outcomes, so its randomized-PIT
// interval is [P(D=0), 1] when D_i = 1 (detected) and [0, P(D=0)] when D_i = 0
// (all-zero), with P(D=0) = psi*(1-p)^n_valid + (1-psi) posterior-averaged over
// draws. This mirrors the cdf_lower/cdf_upper construction cpp_cover_pit_cdf /
// cpp_occu_cover_cdf_limits already use; the randomization itself is left to
// tulpa::tulpa_pit() in R, so no RNG runs here.
// [[Rcpp::export]]
Rcpp::List cpp_single_pit_cdf(
    Rcpp::NumericMatrix X_occ, Rcpp::NumericMatrix X_det,
    Rcpp::NumericMatrix draws, Rcpp::IntegerVector draw_idx,
    Rcpp::IntegerMatrix y
) {
  const int n_sites = X_occ.nrow();
  const int p_occ = X_occ.ncol(), p_det = X_det.ncol();
  const int max_v = y.ncol(); const int ndr = draws.nrow();
  const int n_draws = draw_idx.size();
  shape::check_nrow(X_det, n_sites, "X_det");
  shape::check_nrow(y, n_sites, "y");
  shape::check_ncol_min(draws, (R_xlen_t) p_occ + p_det, "draws");
  shape::check_index1(draw_idx, ndr, "draw_idx");
  // The site limit is a posterior mean over the selected draws, so an empty
  // selection has no value to report.
  if (n_draws < 1) {
    Rcpp::stop("draw_idx must select at least one draw; got 0.");
  }
  Rcpp::NumericVector lower(n_sites), upper(n_sites);
  const double* pXo = X_occ.begin(); const double* pXd = X_det.begin();
  const double* pdr = draws.begin(); const int* py = y.begin();

  for (int i = 0; i < n_sites; ++i) {
    int n_valid = 0, n_det = 0;
    for (int j = 0; j < max_v; ++j) {
      int yij = py[(std::size_t) j * n_sites + i];
      if (yij >= 0) { ++n_valid; if (yij == 1) ++n_det; }
    }
    if (n_valid == 0) { lower[i] = 0.0; upper[i] = 1.0; continue; }
    double acc = 0.0;
    for (int s = 0; s < n_draws; ++s) {
      int idx = draw_idx[s] - 1;
      double psi = stable_plogis(row_draw_dot(pXo, n_sites, i, pdr, ndr, idx, 0, p_occ));
      double p   = stable_plogis(row_draw_dot(pXd, n_sites, i, pdr, ndr, idx, p_occ, p_det));
      acc += psi * std::pow(1.0 - p, (double) n_valid) + (1.0 - psi);
    }
    double q = acc / n_draws;  // posterior-mean P(all-zero) at this site
    if (n_det > 0) { lower[i] = q;   upper[i] = 1.0; }
    else           { lower[i] = 0.0; upper[i] = q; }
  }
  return Rcpp::List::create(Rcpp::Named("cdf_lower") = lower,
                            Rcpp::Named("cdf_upper") = upper);
}
