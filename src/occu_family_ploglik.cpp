// occu_family_ploglik.cpp
// Parallel pointwise log-likelihood kernels for the draw-matrix occupancy
// families whose per-observation marginal was a pure-R loop in R/diagnostics.R
// (single-season replicated, multi-source integrated, multi-season dynamic).
// Each family's R marginal is the oracle (reproduced in the tests); these ports
// mirror it and parallelise over the observation index (each column of the
// [S x N] output is independent, so there are no shared writes). The linear
// predictors arrive as [S x N] matrices built by the R caller (BLAS), so the
// kernel does only the per-observation latent-state marginal.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using tulpaObs::logsumexp2;

namespace {

inline double log_plogis(double x) {              // log(plogis(x))
  if (x >= 0.0) return -std::log1p(std::exp(-x));
  return x - std::log1p(std::exp(x));
}
inline double log_1m_plogis(double x) { return log_plogis(-x); }

}  // namespace

// Single-season occupancy (.tobs_ploglik_replicated): per replicate row i,
// latent z marginalised. eta_psi / eta_p are [S x N]; y is [N x max_visits]
// with entries < 0 marking an invalid visit. Column i of the output is the
// per-draw log-likelihood of observation i.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_occu_single_ploglik(
    Rcpp::NumericMatrix eta_psi,   // [S x N]
    Rcpp::NumericMatrix eta_p,     // [S x N]
    Rcpp::IntegerMatrix y,         // [N x max_visits], < 0 = NA
    int n_threads
) {
  const int S = eta_psi.nrow();
  const int N = eta_psi.ncol();
  const int mv = y.ncol();
  if (eta_p.nrow() != S || eta_p.ncol() != N)
    Rcpp::stop("eta_p must match eta_psi dims.");
  if (y.nrow() != N) Rcpp::stop("y must have N rows.");

  Rcpp::NumericMatrix ll(S, N);
  const double* ppsi = eta_psi.begin();
  const double* pp   = eta_p.begin();
  const int*    py   = y.begin();               // column-major [N x mv]
  double* pll        = ll.begin();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads > 0 ? n_threads : 1)
#endif
  for (int i = 0; i < N; ++i) {
    int nv = 0, k1 = 0;
    for (int v = 0; v < mv; ++v) {
      int yiv = py[(std::size_t) v * N + i];
      if (yiv >= 0) { ++nv; if (yiv == 1) ++k1; }
    }
    double* col = pll + (std::size_t) i * S;
    if (nv == 0) continue;                        // no data -> 0 contribution
    int k0 = nv - k1;
    const double* ep  = pp   + (std::size_t) i * S;
    const double* eps = ppsi + (std::size_t) i * S;
    for (int s = 0; s < S; ++s) {
      double lp    = log_plogis(ep[s]);
      double l1mp  = log_1m_plogis(ep[s]);
      double lpsi  = log_plogis(eps[s]);
      if (k1 > 0) {
        col[s] = lpsi + k1 * lp + k0 * l1mp;
      } else {
        col[s] = logsumexp2(lpsi + nv * l1mp, log_1m_plogis(eps[s]));
      }
    }
  }
  return ll;
}

// Multi-season dynamic occupancy (.tobs_ploglik_dynamic): per-site forward HMM
// recursion in log space. eta_psi1 / eta_p / eta_gam / eta_eps are [S x n_sites];
// y is the flat [n_sites x max_visits x n_seasons] array (< 0 = NA); n_visits and
// any_detected are length n_sites * n_seasons, indexed site-major
// (idx = i * n_seasons + t). Column i of the output is the per-draw site
// log-likelihood.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_occu_dynamic_ploglik(
    Rcpp::NumericMatrix eta_psi1,   // [S x n_sites]
    Rcpp::NumericMatrix eta_p,
    Rcpp::NumericMatrix eta_gam,
    Rcpp::NumericMatrix eta_eps,
    Rcpp::IntegerVector y,          // flat [n_sites x mv x Tn], < 0 = NA
    Rcpp::IntegerVector n_visits,   // [n_sites * Tn], site-major
    Rcpp::IntegerVector any_detected,
    int n_sites, int max_visits, int n_seasons,
    int n_threads
) {
  const int S = eta_psi1.nrow();
  const int mv = max_visits, Tn = n_seasons;
  const double NEG = -1e10;
  Rcpp::NumericMatrix out(S, n_sites);
  const double* ppsi = eta_psi1.begin();
  const double* pp   = eta_p.begin();
  const double* pg   = eta_gam.begin();
  const double* pe   = eta_eps.begin();
  const int* py = y.begin();
  const int* pnv = n_visits.begin();
  const int* pad = any_detected.begin();
  double* pout = out.begin();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads > 0 ? n_threads : 1)
#endif
  for (int i = 0; i < n_sites; ++i) {
    const double* eps1 = ppsi + (std::size_t) i * S;
    const double* ep   = pp   + (std::size_t) i * S;
    const double* eg   = pg   + (std::size_t) i * S;
    const double* ee   = pe   + (std::size_t) i * S;
    double* col = pout + (std::size_t) i * S;
    for (int s = 0; s < S; ++s) {
      double lp = log_plogis(ep[s]),   l1mp   = log_1m_plogis(ep[s]);
      double lgam = log_plogis(eg[s]), l1mgam = log_1m_plogis(eg[s]);
      double leps = log_plogis(ee[s]), l1meps = log_1m_plogis(ee[s]);
      double a_occ = log_plogis(eps1[s]), a_un = log_1m_plogis(eps1[s]);
      double site_ll = 0.0;
      for (int t = 0; t < Tn; ++t) {
        int idx = i * Tn + t;
        int nv = pnv[idx];
        if (nv > 0) {
          int k1 = 0, k0 = 0;
          for (int v = 0; v < mv; ++v) {
            int yv = py[(std::size_t) i + (std::size_t) v * n_sites +
                        (std::size_t) t * n_sites * mv];
            if (yv >= 0) { if (yv == 1) ++k1; else ++k0; }
          }
          if (pad[idx] != 0) {
            site_ll += a_occ + (k1 * lp + k0 * l1mp);
            a_occ = 0.0; a_un = NEG;               // z_t = 1 known
          } else {
            double term1 = a_occ + nv * l1mp, term2 = a_un;
            double lnorm = logsumexp2(term1, term2);
            site_ll += lnorm;
            a_occ = term1 - lnorm; a_un = term2 - lnorm;
          }
        }
        if (t < Tn - 1) {
          double new_occ = logsumexp2(a_occ + l1meps, a_un + lgam);
          double new_un  = logsumexp2(a_occ + leps,   a_un + l1mgam);
          a_occ = new_occ; a_un = new_un;
        }
      }
      col[s] = site_ll;
    }
  }
  return out;
}

// Integrated multi-source occupancy (.tobs_ploglik_integrated): shared psi,
// detection summed over the sources that observed each site. The per-(site,
// source) detection counts K1 / K0 are draw-invariant, so the R caller
// precomputes them; eta_src is the [S x n_sites] detection predictor stacked per
// source (column-major slab `src`). Column i marginalises z per site.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_occu_integrated_ploglik(
    Rcpp::NumericMatrix eta_psi,    // [S x n_sites]
    Rcpp::NumericVector eta_src,    // flat [S x n_sites x n_sources]
    Rcpp::IntegerMatrix K1,         // [n_sites x n_sources]
    Rcpp::IntegerMatrix K0,         // [n_sites x n_sources]
    int n_sources, int n_threads
) {
  const int S = eta_psi.nrow();
  const int N = eta_psi.ncol();
  Rcpp::NumericMatrix ll(S, N);
  const double* ppsi = eta_psi.begin();
  const double* psrc = eta_src.begin();
  const int* pk1 = K1.begin();
  const int* pk0 = K0.begin();
  double* pll = ll.begin();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads > 0 ? n_threads : 1)
#endif
  for (int i = 0; i < N; ++i) {
    const double* eps = ppsi + (std::size_t) i * S;
    double* col = pll + (std::size_t) i * S;
    bool any_det = false;
    for (int src = 0; src < n_sources; ++src)
      if (pk1[(std::size_t) src * N + i] > 0) { any_det = true; break; }
    for (int s = 0; s < S; ++s) {
      double log_det_occ = 0.0;
      for (int src = 0; src < n_sources; ++src) {
        int k1 = pk1[(std::size_t) src * N + i];
        int k0 = pk0[(std::size_t) src * N + i];
        if (k1 == 0 && k0 == 0) continue;
        // eta_src slab `src`: [S x N], column i, row s.
        double e = psrc[(std::size_t) src * S * N + (std::size_t) i * S + s];
        log_det_occ += k1 * log_plogis(e) + k0 * log_1m_plogis(e);
      }
      double lpsi = log_plogis(eps[s]);
      col[s] = any_det ? (lpsi + log_det_occ)
                       : logsumexp2(lpsi + log_det_occ, log_1m_plogis(eps[s]));
    }
  }
  return ll;
}
