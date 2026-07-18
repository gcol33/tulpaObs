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
#include "occu_coupling_shared.h"   // pos_log_density -- the fit-kernel positive density
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {

// Stable logistic, matching R's plogis(x) two-branch form.
inline double stable_plogis(double x) {
  if (x >= 0.0) {
    double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }
  double z = std::exp(x);
  return z / (1.0 + z);
}

inline double clamp_eta(double e, double bound) {
  if (e >  bound) return  bound;
  if (e < -bound) return -bound;
  return e;
}

// log(exp(a) + exp(b)), max-shifted -- the C++ analogue of .tobs_logsumexp2.
inline double logsumexp2(double a, double b) {
  double m = a > b ? a : b;
  double s = a > b ? b : a;
  return m + std::log1p(std::exp(s - m));
}

// Positive-arm log-density = the fit kernel src/occu_coupling_shared.h::
// pos_log_density (Beta/Lognormal/Gaussian::log_density; code 0/3/4). Routing
// through it makes WAIC / LOO score the positive arm with the exact density the
// model was fit with: log_safe at the cover boundary (cover exactly 0 or 1) so
// the density is finite rather than -Inf, and no eta clamp (gcol33/tulpaObs#133).
inline double pos_logdens(double y, double eta, double disp, int positive) {
  return tulpaObs::pos_log_density(positive, y, eta, disp);
}

// Dot of design row i (column-major [nrow x p]) with draw-d coefficient row of
// b (column-major [S x p]); ncol columns of b are offset by `boff`.
inline double row_dot(const double* Xcol, int nrow, int i,
                      const double* bcol, int S, int d, int boff, int p) {
  double acc = 0.0;
  for (int j = 0; j < p; ++j) {
    acc += Xcol[(std::size_t) j * nrow + i] *
           bcol[(std::size_t) (boff + j) * S + d];
  }
  return acc;
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
    int positive,                     // 0 lognormal, 3 beta, 4 gaussian (#112)
    double eta_bound,
    int n_threads
) {
  const int n_sites = X_occ.nrow();
  const int V       = site_of_visit.size();
  const int S       = b_occ.nrow();

  const int p_occ      = X_occ.ncol();
  const int p_det_site = X_det_site.ncol();
  const int p_pos_site = X_pos_site.ncol();
  const int p_det_vis  = X_det_visit.ncol();
  const int p_pos_vis  = X_pos_visit.ncol();
  const bool has_det_visit = p_det_vis > 0;
  const bool has_pos_visit = p_pos_vis > 0;

  if (field_occ.nrow() != n_sites || field_occ.ncol() != S ||
      field_pos.nrow() != n_sites || field_pos.ncol() != S) {
    Rcpp::stop("field_occ / field_pos must be [n_sites x S].");
  }
  if (y_det_visit.size() != V || y_pos_visit.size() != V) {
    Rcpp::stop("y_det_visit / y_pos_visit must be length V.");
  }

  Rcpp::NumericMatrix ll(S, n_sites);

  const double* pXocc  = X_occ.begin();
  const double* pXds   = X_det_site.begin();
  const double* pXps   = X_pos_site.begin();
  const double* pXdv   = X_det_visit.begin();
  const double* pXpv   = X_pos_visit.begin();
  const double* pBocc  = b_occ.begin();
  const double* pBdet  = b_det.begin();
  const double* pBpos  = b_pos.begin();
  const double* pFocc  = field_occ.begin();
  const double* pFpos  = field_pos.begin();
  const int*    sov    = site_of_visit.begin();
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
      const double* focc_d = pFocc + (std::size_t) d * n_sites;
      const double* fpos_d = pFpos + (std::size_t) d * n_sites;
      const double  disp_d = disp[d];

      // Per-site: occupancy prob, site-level detection / cover predictors, and
      // reset the per-site accumulators.
      for (int i = 0; i < n_sites; ++i) {
        double eta_psi = row_dot(pXocc, n_sites, i, pBocc, S, d, 0, p_occ) +
                         focc_d[i];
        psi[i]     = stable_plogis(clamp_eta(eta_psi, eta_bound));
        p_site[i]  = row_dot(pXds, n_sites, i, pBdet, S, d, 0, p_det_site);
        ep_site[i] = row_dot(pXps, n_sites, i, pBpos, S, d, 0, p_pos_site) +
                     fpos_d[i];
        slh[i] = 0.0; sl1mp[i] = 0.0; cov[i] = 0.0; ndet[i] = 0;
      }

      // Per-visit: fold the visit-level block onto its site, accumulate the
      // detection mixture terms and the cover density at detected visits.
      for (int v = 0; v < V; ++v) {
        int s = sov[v] - 1;
        double eta_p = p_site[s];
        if (has_det_visit) {
          eta_p += row_dot(pXdv, V, v, pBdet, S, d, p_det_site, p_det_vis);
        }
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
          double ep = ep_site[s];
          if (has_pos_visit) {
            ep += row_dot(pXpv, V, v, pBpos, S, d, p_pos_site, p_pos_vis);
          }
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
