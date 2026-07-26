// occu_cover_diag.cpp
// C++ kernels for the occu_cover() posterior diagnostics whose per-draw loops
// were pure R. Two pieces:
//   * cpp_occu_cover_cdf_limits -- the deterministic per-site detection-summary
//     CDF limits (any-detection vs all-zero), the randomized-PIT / LOO-PIT
//     building block (.occu_cover_pit_cdf_limits). Parallel over draws.
//   * cpp_occu_cover_ppc -- the posterior predictive check (.tobs_ppc_occu_cover,
//     cover_aggregate = "none" path). It draws the latent z, detection replicate
//     y_rep, and cover replicate from R's RNG via the R:: samplers, in the SAME
//     order as the former R loop, so under a fixed seed the result is
//     byte-identical. Serial (the RNG stream is inherently ordered).
// Both take the per-draw linear predictors as [n_sites x S] / [n_sites x max_v]
// blocks the R caller builds by BLAS.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using tulpaObs::stable_plogis;
using tulpaObs::clamp_eta;
using tulpaObs::ppc_stat;

namespace {
// Positive-arm response-scale mean / replicate draw. `positive` follows the
// shared cover scheme (lognormal 0, beta 3, gaussian 4, gcol33/tulpaObs#112).
inline double mean_pos(double eta, double d, int positive) {
  if (positive == 3) return stable_plogis(clamp_eta(eta));                 // beta mean
  if (positive == 4) return eta;                              // gaussian: mu = eta
  return std::exp(clamp_eta(eta) + d * d / 2.0);                 // lognormal mean
}
inline double draw_pos(double eta, double d, int positive) {
  if (positive == 3) { double mu = stable_plogis(clamp_eta(eta)); return R::rbeta(mu * d, (1.0 - mu) * d); }
  if (positive == 4) return R::rnorm(eta, d);                 // gaussian draw
  return std::exp(R::rnorm(eta, d));                          // lognormal draw
}
// Draw-d detection predictor at (site i, visit j): site block + visit block.
inline double eta_p_ij(const double* Xds, int n_sites, int i,
                       const double* bdet, int S, int d, int p_det_site,
                       const double* Xdv, int nvrow, int vrow, int p_det_visit,
                       bool has_dv) {
  double e = 0.0;
  for (int k = 0; k < p_det_site; ++k)
    e += Xds[(std::size_t) k * n_sites + i] * bdet[(std::size_t) k * S + d];
  if (has_dv)
    for (int k = 0; k < p_det_visit; ++k)
      e += Xdv[(std::size_t) k * nvrow + vrow] *
           bdet[(std::size_t) (p_det_site + k) * S + d];
  return e;
}
}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_occu_cover_cdf_limits(
    Rcpp::NumericMatrix X_occ,        // [n_sites x p_occ]
    Rcpp::NumericMatrix X_det_site,   // [n_sites x p_det_site]
    Rcpp::NumericMatrix X_det_visit,  // [n_sites*max_v x p_det_visit] (0 cols none)
    Rcpp::NumericMatrix b_occ,        // [S x p_occ]
    Rcpp::NumericMatrix b_det,        // [S x (p_det_site + p_det_visit)]
    Rcpp::NumericMatrix field_occ,    // [n_sites x S]
    Rcpp::IntegerMatrix valid,        // [n_sites x max_v]
    Rcpp::IntegerVector any_det,      // [n_sites] 0/1
    int n_threads
) {
  const int S = b_occ.nrow();
  const int n_sites = X_occ.nrow();
  const int max_v = valid.ncol();
  const int p_occ = X_occ.ncol();
  const int p_det_site = X_det_site.ncol();
  const int p_det_visit = X_det_visit.ncol();
  const bool has_dv = p_det_visit > 0;
  const int nvrow = n_sites * max_v;

  Rcpp::NumericMatrix Fl(S, n_sites), Fu(S, n_sites);
  const double* pXo = X_occ.begin(); const double* pXds = X_det_site.begin();
  const double* pXdv = X_det_visit.begin();
  const double* pbo = b_occ.begin(); const double* pbd = b_det.begin();
  const double* pfo = field_occ.begin();
  const int* pv = valid.begin(); const int* pad = any_det.begin();
  double* pFl = Fl.begin(); double* pFu = Fu.begin();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads > 0 ? n_threads : 1)
#endif
  for (int d = 0; d < S; ++d) {
    for (int i = 0; i < n_sites; ++i) {
      double eta_psi = pfo[(std::size_t) d * n_sites + i];
      for (int k = 0; k < p_occ; ++k)
        eta_psi += pXo[(std::size_t) k * n_sites + i] * pbo[(std::size_t) k * S + d];
      double psi = stable_plogis(clamp_eta(eta_psi));
      double sum_l1mp = 0.0;
      for (int j = 0; j < max_v; ++j) {
        if (pv[(std::size_t) j * n_sites + i] == 0) continue;
        int vrow = i * max_v + j;
        double p = stable_plogis(clamp_eta(eta_p_ij(pXds, n_sites, i, pbd, S, d, p_det_site,
                                       pXdv, nvrow, vrow, p_det_visit, has_dv)));
        sum_l1mp += std::log(1.0 - p);
      }
      double pdet0 = psi * std::exp(sum_l1mp) + (1.0 - psi);
      std::size_t off = (std::size_t) i * S + d;
      if (pad[i] != 0) { pFl[off] = pdet0; pFu[off] = 1.0; }
      else             { pFl[off] = 0.0;   pFu[off] = pdet0; }
    }
  }
  return Rcpp::List::create(Rcpp::Named("cdf_lower") = Fl,
                            Rcpp::Named("cdf_upper") = Fu);
}

// Posterior predictive check for occu_cover(), cover_aggregate = "none". Serial:
// the RNG draws (R::rbinom for z and detection y_rep, R::rbeta / R::rnorm for the
// cover replicate) follow the exact order of the former R loop, so under a fixed
// seed the returned discrepancies are byte-identical. `freeman` selects the
// Freeman-Tukey (else chi-squared) statistic. Predictors are per-draw [n_sites x
// max_v] blocks: p_mat (detection prob), ep_mat (cover linear predictor);
// psi_all [n_sites x S]. `disp` is per draw.
// [[Rcpp::export]]
Rcpp::List cpp_occu_cover_ppc(
    Rcpp::NumericMatrix psi_all,      // [n_sites x S]
    Rcpp::NumericMatrix p_all,        // [n_sites x (S*max_v)] draw-major visit blocks
    Rcpp::NumericMatrix ep_all,       // [n_sites x (S*max_v)] cover predictor
    Rcpp::IntegerMatrix y,            // [n_sites x max_v]
    Rcpp::NumericMatrix y_pos,        // [n_sites x max_v]
    Rcpp::IntegerMatrix valid,        // [n_sites x max_v]
    Rcpp::IntegerVector any_det,      // [n_sites]
    Rcpp::IntegerVector n_valid,      // [n_sites]
    Rcpp::NumericVector disp,         // [S]
    int positive, bool freeman        // positive: 0 lognormal, 3 beta, 4 gaussian
) {
  const int S = psi_all.ncol();
  const int n_sites = psi_all.nrow();
  const int max_v = valid.ncol();
  Rcpp::RNGScope scope;                      // ties to R's RNG stream

  Rcpp::NumericVector fit_y(S), fit_rep(S);
  const int* pv = valid.begin(); const int* py = y.begin();
  const double* pyp = y_pos.begin();
  const int* pad = any_det.begin(); const int* pnv = n_valid.begin();

  auto stat = [&](double o, double e) { return ppc_stat(o, e, freeman); };
  std::vector<double> zmass(n_sites);

  for (int s = 0; s < S; ++s) {
    const double* psi = &psi_all[(std::size_t) s * n_sites];
    const double* pmat = &p_all[(std::size_t) s * max_v * n_sites];
    const double* epmat = &ep_all[(std::size_t) s * max_v * n_sites];
    double d = disp[s];

    // z_prob per site (deterministic), then z ~ Bernoulli in site order.
    for (int i = 0; i < n_sites; ++i) {
      double prod1mp = 0.0;                 // log
      for (int j = 0; j < max_v; ++j)
        if (pv[(std::size_t) j * n_sites + i])
          prod1mp += std::log(1.0 - pmat[(std::size_t) j * n_sites + i]);
      double pr = std::exp(prod1mp);
      double zp;
      if (pnv[i] == 0) zp = psi[i];
      else if (pad[i]) zp = 1.0;
      else zp = psi[i] * pr / (psi[i] * pr + (1.0 - psi[i]));
      zmass[i] = zp;
    }
    std::vector<int> z(n_sites);
    for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, zmass[i]);

    // Detection replicate y_rep ~ Bernoulli(z p), column-major [n_sites x max_v].
    double det_obs = 0.0, det_rep = 0.0;
    std::vector<double> exp_det((std::size_t) n_sites * max_v);
    for (int j = 0; j < max_v; ++j)
      for (int i = 0; i < n_sites; ++i)
        exp_det[(std::size_t) j * n_sites + i] = z[i] * pmat[(std::size_t) j * n_sites + i];
    std::vector<int> yrep((std::size_t) n_sites * max_v);
    for (std::size_t idx = 0; idx < (std::size_t) n_sites * max_v; ++idx)
      yrep[idx] = (int) R::rbinom(1.0, exp_det[idx]);
    for (int j = 0; j < max_v; ++j)
      for (int i = 0; i < n_sites; ++i) {
        std::size_t idx = (std::size_t) j * n_sites + i;
        if (!pv[idx]) continue;
        det_obs += stat((double) py[idx], exp_det[idx]);
        det_rep += stat((double) yrep[idx], exp_det[idx]);
      }

    // Cover replicate over all cells (draw_pos on as.vector(ep_mat)), then the
    // positive-part discrepancy on detected cells (Epos = mean_pos).
    std::vector<double> yrep_cov((std::size_t) n_sites * max_v);
    for (std::size_t idx = 0; idx < (std::size_t) n_sites * max_v; ++idx) {
      double eta = epmat[idx];
      yrep_cov[idx] = draw_pos(eta, d, positive);
    }
    double cov_obs = 0.0, cov_rep = 0.0;
    for (int j = 0; j < max_v; ++j)
      for (int i = 0; i < n_sites; ++i) {
        std::size_t idx = (std::size_t) j * n_sites + i;
        if (!(pv[idx] && py[idx] == 1)) continue;
        double eta = epmat[idx];
        double Epos = mean_pos(eta, d, positive);
        cov_obs += stat(pyp[idx], Epos);
        cov_rep += stat(yrep_cov[idx], Epos);
      }
    fit_y[s] = det_obs + cov_obs;
    fit_rep[s] = det_rep + cov_rep;
  }
  return Rcpp::List::create(Rcpp::Named("fit.y") = fit_y,
                            Rcpp::Named("fit.y.rep") = fit_rep);
}

// Aggregated / latent-mode PPC for occu_cover() (cover_aggregate = "mean" /
// "median" / "latent"). The detection replicate is identical to the none-mode
// path; the cover replicate is one aggregated cover per detected unit (mode 1,
// against the precomputed observed aggregate `yv`) or the shared cover-RE
// marginal (mode 2: a per-unit RE u ~ N(0, disp) then a per-visit draw at
// dispersion `disp2`). Units are CSR: `pos_site[u]` the cell, `vals_flat` /
// `unit_off` the detected covers. RNG order matches the R loop, so byte-identical.
// [[Rcpp::export]]
Rcpp::List cpp_occu_cover_ppc_agg(
    Rcpp::NumericMatrix psi_all, Rcpp::NumericMatrix p_all,
    Rcpp::NumericMatrix ep_all, Rcpp::IntegerMatrix y, Rcpp::IntegerMatrix valid,
    Rcpp::IntegerVector any_det, Rcpp::IntegerVector n_valid,
    Rcpp::NumericVector disp, int mode_code, Rcpp::IntegerVector pos_site,
    Rcpp::NumericVector yv, Rcpp::NumericVector vals_flat,
    Rcpp::IntegerVector unit_off, double disp2, int positive, bool freeman
) {
  const int S = psi_all.ncol(), n_sites = psi_all.nrow(), max_v = valid.ncol();
  const int n_units = pos_site.size();
  Rcpp::RNGScope scope;
  Rcpp::NumericVector fit_y(S), fit_rep(S);
  const int* pv = valid.begin(); const int* py = y.begin();
  const int* pad = any_det.begin(); const int* pnv = n_valid.begin();
  const int* pps = pos_site.begin();
  auto stat = [&](double o, double e) { return ppc_stat(o, e, freeman); };
  std::vector<double> zmass(n_sites);

  for (int s = 0; s < S; ++s) {
    const double* psi = &psi_all[(std::size_t) s * n_sites];
    const double* pmat = &p_all[(std::size_t) s * max_v * n_sites];
    const double* epmat = &ep_all[(std::size_t) s * max_v * n_sites];
    double d = disp[s];
    // Detection replicate (identical to the none-mode path).
    for (int i = 0; i < n_sites; ++i) {
      double lp = 0.0;
      for (int j = 0; j < max_v; ++j)
        if (pv[(std::size_t) j * n_sites + i]) lp += std::log(1.0 - pmat[(std::size_t) j * n_sites + i]);
      double pr = std::exp(lp), zp;
      if (pnv[i] == 0) zp = psi[i];
      else if (pad[i]) zp = 1.0;
      else zp = psi[i] * pr / (psi[i] * pr + (1.0 - psi[i]));
      zmass[i] = zp;
    }
    std::vector<int> z(n_sites);
    for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, zmass[i]);
    std::vector<double> exp_det((std::size_t) n_sites * max_v);
    for (int j = 0; j < max_v; ++j)
      for (int i = 0; i < n_sites; ++i)
        exp_det[(std::size_t) j * n_sites + i] = z[i] * pmat[(std::size_t) j * n_sites + i];
    std::vector<int> yrep((std::size_t) n_sites * max_v);
    for (std::size_t k = 0; k < (std::size_t) n_sites * max_v; ++k)
      yrep[k] = (int) R::rbinom(1.0, exp_det[k]);
    double det_obs = 0.0, det_rep = 0.0;
    for (int j = 0; j < max_v; ++j)
      for (int i = 0; i < n_sites; ++i) {
        std::size_t k = (std::size_t) j * n_sites + i;
        if (!pv[k]) continue;
        det_obs += stat((double) py[k], exp_det[k]);
        det_rep += stat((double) yrep[k], exp_det[k]);
      }

    // Aggregated / latent cover replicate.
    double cov_obs = 0.0, cov_rep = 0.0;
    if (n_units > 0) {
      if (mode_code == 1) {                        // mean / median
        for (int u = 0; u < n_units; ++u) {
          double eta = epmat[pps[u]];              // ep_mat[pos_site, 1]
          double Ep = mean_pos(eta, d, positive);
          double rp = draw_pos(eta, d, positive);
          cov_obs += stat(yv[u], Ep);
          cov_rep += stat(rp, Ep);
        }
      } else {                                     // latent
        std::vector<double> ure(n_units);
        for (int u = 0; u < n_units; ++u) ure[u] = R::rnorm(0.0, d);
        for (int u = 0; u < n_units; ++u) {
          double eta = epmat[pps[u]];
          for (int t = unit_off[u]; t < unit_off[u + 1]; ++t) {
            double e = mean_pos(eta, disp2, positive);
            cov_obs += stat(vals_flat[t], e);
            cov_rep += stat(draw_pos(eta + ure[u], disp2, positive), e);
          }
        }
      }
    }
    fit_y[s] = det_obs + cov_obs;
    fit_rep[s] = det_rep + cov_rep;
  }
  return Rcpp::List::create(Rcpp::Named("fit.y") = fit_y,
                            Rcpp::Named("fit.y.rep") = fit_rep);
}
