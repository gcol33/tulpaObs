// ms_ocs_ploglik.cpp
// Batched per-(cell, species) pointwise log-likelihood for the spatial-factor
// community occupancy + cover family (ms_occu_cover_spatial). The former R loop
// (.tobs_ploglik_ms_occu_cover_spatial) unpacked the NUTS draw (community mean
// mu, per-species deviation b, shared fields W, occupancy loadings L, optional
// cover loadings Lpos, log-dispersion), assembled each species' occ / detection
// / cover predictors -- with the shared-factor offset W L[s,] on psi (and
// W Lpos[s,] on cover) -- and evaluated the dense occu_cover per-cell marginal
// (.occu_cover_site_ll). All of that now runs in C++, parallel over draws; the
// per-cell z-marginal + cover density mirror the dense occu_cover kernel, so the
// result is byte-close (~1e-13) to the R oracle. Output is [M x (N * S)] with
// species blocks contiguous (column s*N + c), matching as.numeric(LL).

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
using tulpaObs::logsumexp2;

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_ms_ocs_ploglik(
    Rcpp::NumericMatrix draws,        // [M x total_inner]
    Rcpp::NumericMatrix X_occ,        // [N x P_occ]
    Rcpp::NumericMatrix X_det_site,   // [N x P_p_site]
    Rcpp::NumericMatrix X_det_visit,  // [N*J x P_p_visit] (0 cols if none)
    Rcpp::NumericMatrix X_pos_site,   // [N x P_pos_site]
    Rcpp::NumericMatrix X_pos_visit,  // [N*J x P_pos_visit] (0 cols if none)
    Rcpp::IntegerVector y,            // [N*J*S] flat (site-major within species)
    Rcpp::NumericVector y_pos,        // [N*J*S]
    Rcpp::IntegerVector valid,        // [N*J*S] 0/1
    int N, int J, int S, int K,
    int P_occ, int P_p, int P_pos, int P_p_site, int P_pos_site,
    bool cover_factor, bool is_beta, int n_threads
) {
  const int M = draws.nrow();
  const int P = P_occ + P_p + P_pos;
  const int P_p_visit = P_p - P_p_site;
  const int P_pos_visit = P_pos - P_pos_site;
  // Packed-parameter offsets (0-based), matching .ms_ocs_unpack.
  const int off_mu = 0;
  const int off_b  = P;
  const int off_L  = off_b + S * P;
  const int off_Lpos = off_L + S * K;
  const int off_W  = off_Lpos + (cover_factor ? S * K : 0);
  const int off_ld = off_W + N * K;

  Rcpp::NumericMatrix out(M, (std::size_t) N * S);
  const double* pd = draws.begin();
  const int* py = y.begin(); const double* pyp = y_pos.begin();
  const int* pv = valid.begin();
  double* pout = out.begin();
  const std::size_t sp_stride = (std::size_t) N * J;   // per-species y slab

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
#endif
  {
    std::vector<double> thocc(P_occ), thp(P_p), thpos(P_pos);
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int m = 0; m < M; ++m) {
      auto dr = [&](int off) { return pd[(std::size_t) off * M + m]; };
      double ld = dr(off_ld);
      double disp = std::exp(ld);

      for (int s = 0; s < S; ++s) {
        // th = mu + b_s, split into arms.
        for (int k = 0; k < P_occ; ++k)
          thocc[k] = dr(off_mu + k) + dr(off_b + s * P + k);
        for (int k = 0; k < P_p; ++k)
          thp[k] = dr(off_mu + P_occ + k) + dr(off_b + s * P + P_occ + k);
        for (int k = 0; k < P_pos; ++k)
          thpos[k] = dr(off_mu + P_occ + P_p + k) + dr(off_b + s * P + P_occ + P_p + k);

        const int* y_s = py + (std::size_t) s * sp_stride;
        const double* yp_s = pyp + (std::size_t) s * sp_stride;
        const int* v_s = pv + (std::size_t) s * sp_stride;

        for (int c = 0; c < N; ++c) {
          // Occupancy predictor with the shared-factor offset W L[s,].
          double eta_psi = 0.0;
          for (int k = 0; k < P_occ; ++k) eta_psi += X_occ((std::size_t) c, k) * thocc[k];
          for (int f = 0; f < K; ++f)
            eta_psi += dr(off_W + f * N + c) * dr(off_L + f * S + s);
          double psi = stable_plogis(clamp_eta(eta_psi));
          double lpsi = std::log(psi), l1mpsi = std::log(1.0 - psi);

          // Site-level detection / cover predictors (broadcast across visits).
          double ep_site = 0.0, pos_site = 0.0;
          for (int k = 0; k < P_p_site; ++k) ep_site += X_det_site((std::size_t) c, k) * thp[k];
          for (int k = 0; k < P_pos_site; ++k) pos_site += X_pos_site((std::size_t) c, k) * thpos[k];
          double field_pos = 0.0;
          if (cover_factor)
            for (int f = 0; f < K; ++f)
              field_pos += dr(off_W + f * N + c) * dr(off_Lpos + f * S + s);

          double sum_hdet = 0.0, sum_1mp = 0.0, cover = 0.0;
          int ndet = 0;
          for (int j = 0; j < J; ++j) {
            std::size_t cj = (std::size_t) c + (std::size_t) j * N;   // [N x J] col-major
            if (v_s[cj] == 0) continue;
            int vrow = c * J + j;                    // site-major visit design row
            double eta_p = ep_site;
            for (int k = 0; k < P_p_visit; ++k)
              eta_p += X_det_visit((std::size_t) vrow, k) * thp[P_p_site + k];
            double p = stable_plogis(clamp_eta(eta_p));
            double lp = std::log(p), l1mp = std::log(1.0 - p);
            sum_1mp += l1mp;
            int yij = y_s[cj];
            if (yij == 1) {
              ndet++;
              sum_hdet += lp;
              double ep = pos_site + field_pos;
              for (int k = 0; k < P_pos_visit; ++k)
                ep += X_pos_visit((std::size_t) vrow, k) * thpos[P_pos_site + k];
              double cv = yp_s[cj];
              double dens;
              if (is_beta) {
                double mu = stable_plogis(clamp_eta(ep));
                double a = mu * disp, b = (1.0 - mu) * disp;
                dens = std::lgamma(disp) - std::lgamma(a) - std::lgamma(b) +
                       (a - 1.0) * std::log(cv) + (b - 1.0) * std::log(1.0 - cv);
              } else {
                double ly = std::log(cv), r = (ly - ep) / disp;
                dens = -ly - std::log(disp) - 0.5 * std::log(2.0 * M_PI) - 0.5 * r * r;
              }
              cover += dens;
            } else {
              sum_hdet += l1mp;
            }
          }
          double val = (ndet > 0) ? (lpsi + sum_hdet + cover)
                                  : logsumexp2(lpsi + sum_1mp, l1mpsi);
          std::size_t col = (std::size_t) s * N + c;
          pout[col * M + m] = val;
        }
      }
    }
  }
  return out;
}
