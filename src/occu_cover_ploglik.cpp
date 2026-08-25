// occu_cover_ploglik.cpp
// Parallel pointwise log-likelihood for the compact (ragged) occu_cover() fit,
// the input WAIC / PSIS-LOO / LOO-stacking consume. The R reference is
// .occu_cover_site_ll_ragged (+ .occu_cover_draw_eta_ragged /
// .occu_cover_eta_components) in R/occu_cover.R; this port mirrors it draw for
// draw and is cross-checked byte-close against it (test-ploglik-cpp.R).
//
// The per-site marginal integrates the latent occupancy z in closed form over
// its two states, exactly as the fit's negative-log-posterior does:
//   det > 0 : log psi + sum_v log h(y_v) + sum_{det} log f_pos(cover_v)
//   det = 0 : logsumexp( log psi + sum_v log(1 - p_v),  log(1 - psi) )
// with h(y_v) = p_v if y_v = 1 else (1 - p_v). Each draw's [n_sites] row is
// independent, so the draw loop parallelises with no shared writes -- the one
// axis WAIC scales on (n.draws x total plots). Per-visit sums accumulate in
// visit (index) order, matching R's rowsum(), so the two paths agree to the
// bit modulo libm plogis/log/lgamma rounding.
//
// Memory: per-visit predictors are formed one draw at a time in thread-private
// scratch (never the [V x n_draws] transient the R chunker guards), so the
// draw count is not RAM-bounded here.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_shape.h"
#include "tobs_math.h"
#include "occu_cover_ragged.h"      // Arms -- the shared per-draw predictor view
#include "occu_coupling_shared.h"   // pos_log_density -- the fit-kernel positive density
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using tulpaObs::stable_plogis;
using tulpaObs::clamp_eta;
using tulpaObs::logsumexp2;
using tulpaObs::occu_cover_ragged::Arms;

namespace {

// Positive-arm log-density = the fit kernel src/occu_coupling_shared.h::
// pos_log_density (Beta/Lognormal/Gaussian::log_density; code 0/3/4). Routing
// through it makes WAIC / LOO score the positive arm with the exact density the
// model was fit with: log_safe at the cover boundary (cover exactly 0 or 1) so
// the density is finite rather than -Inf, and no eta clamp.
inline double pos_logdens(double y, double eta, double disp, int positive) {
  return tulpaObs::pos_log_density(positive, y, eta, disp);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_occu_cover_ploglik_ragged(
    Rcpp::NumericMatrix X_occ,        // [n_sites x p_occ]
    Rcpp::NumericMatrix X_det_site,   // [n_sites x p_det_site]
    Rcpp::NumericMatrix X_pos_site,   // [n_sites x p_pos_site]
    Rcpp::NumericMatrix X_det_visit,  // [V x p_det_visit] (0 cols if absent)
    Rcpp::NumericMatrix X_pos_visit,  // [V x p_pos_visit] (0 cols if absent)
    Rcpp::IntegerVector site_of_visit,// [V], 1-based
    Rcpp::IntegerVector y_det_visit,  // [V], 0/1
    Rcpp::NumericVector y_pos_visit,  // [V], cover value (0 where undetected)
    Rcpp::NumericMatrix b_occ,        // [S x p_occ]
    Rcpp::NumericMatrix b_det,        // [S x (p_det_site + p_det_visit)]
    Rcpp::NumericMatrix b_pos,        // [S x (p_pos_site + p_pos_visit)]
    Rcpp::NumericVector disp,         // [S]
    Rcpp::NumericMatrix field_occ,    // [n_sites x S]
    Rcpp::NumericMatrix field_pos,    // [n_sites x S]
    Rcpp::NumericMatrix off_det,      // [V x S] (0 cols if absent)
    Rcpp::NumericMatrix off_pos,      // [V x S] (0 cols if absent)
    int positive,                     // 0 lognormal, 3 beta, 4 gaussian (#112)
    double eta_bound,
    int n_threads
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

  const int n_sites = arms.n_sites;
  const int V       = arms.V;
  const int S       = arms.S;
  namespace sh = tulpaObs::shape;
  sh::check_len(y_det_visit, V, "y_det_visit");
  sh::check_len(y_pos_visit, V, "y_pos_visit");
  sh::check_len(disp, S, "disp");

  Rcpp::NumericMatrix ll(S, n_sites);

  const int*    ydet   = y_det_visit.begin();
  const double* ypos   = y_pos_visit.begin();
  double* pll          = ll.begin();          // column-major [S x n_sites]

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
#endif
  {
    std::vector<double> psi(n_sites), ep_site(n_sites), p_site(n_sites);
    std::vector<double> slh(n_sites), sl1mp(n_sites), cov(n_sites);
    std::vector<int>    ndet(n_sites);

#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int d = 0; d < S; ++d) {
      const double  disp_d = disp[d];

      // Per-site: occupancy prob, site-level detection / cover predictors, and
      // reset the per-site accumulators.
      for (int i = 0; i < n_sites; ++i) {
        psi[i]     = arms.psi(i, d);
        p_site[i]  = arms.eta_p_site(i, d);
        ep_site[i] = arms.eta_pos_site(i, d);
        slh[i] = 0.0; sl1mp[i] = 0.0; cov[i] = 0.0; ndet[i] = 0;
      }

      // Per-visit: fold the visit-level block onto its site, accumulate the
      // detection mixture terms and the cover density at detected visits.
      for (int v = 0; v < V; ++v) {
        int s = arms.site(v);
        double eta_p = p_site[s] + arms.eta_p_visit(v, d);
        double p    = stable_plogis(clamp_eta(eta_p, eta_bound));
        double l1mp = std::log(1.0 - p);
        int    y    = ydet[v];
        slh[s]   += (y == 1) ? std::log(p) : l1mp;
        sl1mp[s] += l1mp;
        ndet[s]  += y;
        // Cover term at detected visits with an observed cover; a missing (NA ->
        // non-finite) cover drops out (missing-at-random cover), the detection
        // mixture above still counts it.
        if (y == 1 && std::isfinite(ypos[v])) {
          double ep = ep_site[s] + arms.eta_pos_visit(v, d);
          cov[s] += pos_logdens(ypos[v], ep, disp_d, positive);
        }
      }

      // Per-site: fold in the closed-form occupancy marginal.
      for (int i = 0; i < n_sites; ++i) {
        double one_m = 1.0 - psi[i];
        double lpsi   = std::log(psi[i] > 1e-300 ? psi[i] : 1e-300);
        double l1mpsi = std::log(one_m  > 1e-300 ? one_m  : 1e-300);
        double val;
        if (ndet[i] > 0) {
          val = lpsi + slh[i] + cov[i];
        } else {
          val = logsumexp2(lpsi + sl1mp[i], l1mpsi);
        }
        pll[(std::size_t) i * S + d] = val;
      }
    }
  }

  return ll;
}
