// occu_cover_diag.cpp
// C++ kernels for the occu_cover() posterior diagnostics whose per-draw loops
// were pure R. Three pieces:
//   * cpp_occu_cover_cdf_limits -- the deterministic per-site detection-summary
//     CDF limits (any-detection vs all-zero), the randomized-PIT / LOO-PIT
//     building block (.occu_cover_pit_cdf_limits). Parallel over draws.
//   * cpp_occu_cover_ppc -- the posterior predictive check (.tobs_ppc_occu_cover,
//     cover_aggregate = "none" path). It draws the latent z, detection replicate
//     y_rep, and cover replicate from R's RNG via the R:: samplers. Serial (the
//     RNG stream is inherently ordered).
//   * cpp_occu_cover_ppc_agg -- the same check for the aggregated / latent cover
//     modes, which score one cover term per detected UNIT rather than per visit.
//
// The first two read the compact (one row per valid visit) layout described in
// occu_cover_ragged.h and assemble their predictors from the shared `Arms`
// view, so a dense fit -- flattened to that layout by .occu_cover_visit_view()
// in R -- and a compact fit of the same data give identical results
// (gcol33/tulpaObs#185). The aggregated modes are dense-only (compact input is
// gated to cover_aggregate = "none"), so cpp_occu_cover_ppc_agg keeps the
// padded [n_sites x max_visits] grid and takes its predictors precomputed.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "occu_cover_ragged.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using tulpaObs::stable_plogis;
using tulpaObs::clamp_eta;
using tulpaObs::ppc_stat;
using tulpaObs::occu_cover_ragged::Arms;

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
}  // namespace

// Per-draw CDF limits of the per-site detection summary (any detection vs all
// zero), the latent occupancy state marginalized out. Deterministic, so this
// parallelises over draws.
// [[Rcpp::export]]
Rcpp::List cpp_occu_cover_cdf_limits(
    Rcpp::NumericMatrix X_occ,        // [n_sites x p_occ]
    Rcpp::NumericMatrix X_det_site,   // [n_sites x p_det_site]
    Rcpp::NumericMatrix X_det_visit,  // [V x p_det_visit] (0 cols if absent)
    Rcpp::IntegerVector site_of_visit,// [V], 1-based
    Rcpp::NumericMatrix b_occ,        // [S x p_occ]
    Rcpp::NumericMatrix b_det,        // [S x (p_det_site + p_det_visit)]
    Rcpp::NumericMatrix field_occ,    // [n_sites x S]
    Rcpp::NumericMatrix off_det,      // [V x S] (0 cols if absent)
    Rcpp::IntegerVector any_det,      // [n_sites] 0/1
    double eta_bound,
    int n_threads
) {
  Arms arms = tulpaObs::occu_cover_ragged::make_arms(
      X_occ, X_det_site, X_det_visit, site_of_visit, b_occ, b_det, field_occ,
      eta_bound);
  arms.off_det_visit = tulpaObs::occu_cover_ragged::visit_offset(
      off_det, arms.V, arms.S, "off_det");
  const int S = arms.S, n_sites = arms.n_sites, V = arms.V;
  if (any_det.size() != n_sites) {
    Rcpp::stop("any_det must be length n_sites.");
  }

  Rcpp::NumericMatrix Fl(S, n_sites), Fu(S, n_sites);
  const int* pad = any_det.begin();
  double* pFl = Fl.begin(); double* pFu = Fu.begin();

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
#endif
  {
    std::vector<double> p_site(n_sites), sum_l1mp(n_sites);

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int d = 0; d < S; ++d) {
      for (int i = 0; i < n_sites; ++i) {
        p_site[i] = arms.eta_p_site(i, d);
        sum_l1mp[i] = 0.0;
      }
      // A site's log(1 - p) terms accumulate in its own visit order, which the
      // compact layout keeps ascending, matching the dense per-site sweep.
      for (int v = 0; v < V; ++v) {
        int s = arms.site(v);
        double p = stable_plogis(clamp_eta(p_site[s] + arms.eta_p_visit(v, d),
                                           eta_bound));
        sum_l1mp[s] += std::log(1.0 - p);
      }
      for (int i = 0; i < n_sites; ++i) {
        double psi = arms.psi(i, d);
        // P(no detection at any visit) = psi prod(1 - p) + (1 - psi).
        double pdet0 = psi * std::exp(sum_l1mp[i]) + (1.0 - psi);
        std::size_t off = (std::size_t) i * S + d;
        if (pad[i] != 0) { pFl[off] = pdet0; pFu[off] = 1.0; }
        else             { pFl[off] = 0.0;   pFu[off] = pdet0; }
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("cdf_lower") = Fl,
                            Rcpp::Named("cdf_upper") = Fu);
}

// Posterior predictive check for occu_cover(), cover_aggregate = "none". Serial:
// the RNG draws (R::rbinom for z and the detection replicate, R::rbeta /
// R::rnorm for the cover replicate) run in a fixed order -- z in site order,
// then one detection and one cover replicate per valid visit in visit order --
// so under a fixed seed the discrepancies are reproducible, and identical
// between a dense and a compact build of the same data. `freeman` selects the
// Freeman-Tukey (else chi-squared) statistic.
// [[Rcpp::export]]
Rcpp::List cpp_occu_cover_ppc(
    Rcpp::NumericMatrix X_occ,        // [n_sites x p_occ]
    Rcpp::NumericMatrix X_det_site,   // [n_sites x p_det_site]
    Rcpp::NumericMatrix X_pos_site,   // [n_sites x p_pos_site]
    Rcpp::NumericMatrix X_det_visit,  // [V x p_det_visit] (0 cols if absent)
    Rcpp::NumericMatrix X_pos_visit,  // [V x p_pos_visit] (0 cols if absent)
    Rcpp::IntegerVector site_of_visit,// [V], 1-based
    Rcpp::IntegerVector y_det_visit,  // [V], 0/1
    Rcpp::NumericVector y_pos_visit,  // [V], cover (NA where unobserved)
    Rcpp::NumericMatrix b_occ,        // [S x p_occ]
    Rcpp::NumericMatrix b_det,        // [S x (p_det_site + p_det_visit)]
    Rcpp::NumericMatrix b_pos,        // [S x (p_pos_site + p_pos_visit)]
    Rcpp::NumericVector disp,         // [S]
    Rcpp::NumericMatrix field_occ,    // [n_sites x S]
    Rcpp::NumericMatrix field_pos,    // [n_sites x S]
    Rcpp::NumericMatrix off_det,      // [V x S] (0 cols if absent)
    Rcpp::NumericMatrix off_pos,      // [V x S] (0 cols if absent)
    Rcpp::IntegerVector any_det,      // [n_sites] 0/1
    Rcpp::IntegerVector n_valid,      // [n_sites]
    int positive, double eta_bound, bool freeman
) {
  Arms arms = tulpaObs::occu_cover_ragged::make_arms(
      X_occ, X_det_site, X_det_visit, site_of_visit, b_occ, b_det, field_occ,
      eta_bound);
  tulpaObs::occu_cover_ragged::attach_cover(arms, X_pos_site, X_pos_visit,
                                            b_pos, field_pos);
  arms.off_det_visit = tulpaObs::occu_cover_ragged::visit_offset(
      off_det, arms.V, arms.S, "off_det");
  arms.off_pos_visit = tulpaObs::occu_cover_ragged::visit_offset(
      off_pos, arms.V, arms.S, "off_pos");
  const int S = arms.S, n_sites = arms.n_sites, V = arms.V;
  if (y_det_visit.size() != V || y_pos_visit.size() != V) {
    Rcpp::stop("y_det_visit / y_pos_visit must be length V.");
  }
  if (any_det.size() != n_sites || n_valid.size() != n_sites) {
    Rcpp::stop("any_det / n_valid must be length n_sites.");
  }
  Rcpp::RNGScope scope;                      // ties to R's RNG stream

  Rcpp::NumericVector fit_y(S), fit_rep(S);
  const int* ydet = y_det_visit.begin();
  const double* ypos = y_pos_visit.begin();
  const int* pad = any_det.begin();
  const int* pnv = n_valid.begin();

  auto stat = [&](double o, double e) { return ppc_stat(o, e, freeman); };
  std::vector<double> psi(n_sites), p_site(n_sites), ep_site(n_sites);
  std::vector<double> sum_l1mp(n_sites);
  std::vector<double> p_vis(V), ep_vis(V);
  std::vector<int> z(n_sites);

  for (int s = 0; s < S; ++s) {
    const double d = disp[s];

    // Deterministic per-draw predictors: occupancy prob, per-visit detection
    // prob, per-visit cover predictor.
    for (int i = 0; i < n_sites; ++i) {
      psi[i]     = arms.psi(i, s);
      p_site[i]  = arms.eta_p_site(i, s);
      ep_site[i] = arms.eta_pos_site(i, s);
      sum_l1mp[i] = 0.0;
    }
    for (int v = 0; v < V; ++v) {
      int i = arms.site(v);
      double p = stable_plogis(clamp_eta(p_site[i] + arms.eta_p_visit(v, s),
                                         eta_bound));
      p_vis[v]  = p;
      ep_vis[v] = ep_site[i] + arms.eta_pos_visit(v, s);
      sum_l1mp[i] += std::log(1.0 - p);
    }

    // Latent z from its full conditional given the detection history: certain
    // at a site with a detection, the prior at a site with no valid visit.
    for (int i = 0; i < n_sites; ++i) {
      double zp;
      if (pnv[i] == 0) zp = psi[i];
      else if (pad[i]) zp = 1.0;
      else {
        double pr = std::exp(sum_l1mp[i]);
        zp = psi[i] * pr / (psi[i] * pr + (1.0 - psi[i]));
      }
      z[i] = (int) R::rbinom(1.0, zp);
    }

    // Detection replicate y_rep ~ Bernoulli(z p) at every valid visit.
    double det_obs = 0.0, det_rep = 0.0;
    for (int v = 0; v < V; ++v) {
      double e = z[arms.site(v)] * p_vis[v];
      double yr = R::rbinom(1.0, e);
      det_obs += stat((double) ydet[v], e);
      det_rep += stat(yr, e);
    }

    // Cover replicate at every valid visit (so the RNG stream does not depend on
    // which visits were detected), scored where a cover was observed: detected
    // AND finite. A detected visit with a missing cover carries no cover term,
    // matching the likelihood's missing-at-random gate.
    double cov_obs = 0.0, cov_rep = 0.0;
    for (int v = 0; v < V; ++v) {
      double eta = ep_vis[v];
      double rp = draw_pos(eta, d, positive);
      if (ydet[v] == 1 && std::isfinite(ypos[v])) {
        double Epos = mean_pos(eta, d, positive);
        cov_obs += stat(ypos[v], Epos);
        cov_rep += stat(rp, Epos);
      }
    }
    fit_y[s] = det_obs + cov_obs;
    fit_rep[s] = det_rep + cov_rep;
  }
  return Rcpp::List::create(Rcpp::Named("fit.y") = fit_y,
                            Rcpp::Named("fit.y.rep") = fit_rep);
}

// Aggregated / latent-mode PPC for occu_cover() (cover_aggregate = "mean" /
// "median" / "latent"). The detection replicate is the padded-grid counterpart
// of the none-mode path; the cover replicate is one aggregated cover per
// detected unit (mode 1, against the precomputed observed aggregate `yv`) or the
// shared cover-RE marginal (mode 2: a per-unit RE u ~ N(0, disp) then a
// per-visit draw at dispersion `disp2`). Units are CSR: `pos_site[u]` the cell,
// `vals_flat` / `unit_off` the detected covers. RNG order matches the R loop,
// so byte-identical.
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
