// ms_nmix_ploglik.cpp
// Batched pointwise log-likelihood for the community N-mixture (ms_abun) family.
// The former R loop (.tobs_ploglik_ms_nmix) reconstructed each species' deviation
// b = C z from the non-centered NUTS draw (log-Cholesky factor C per arm) and
// called the per-species Royle marginal in R. Both the reconstruction and the
// per-(species, site) marginal now run in C++, parallel over draws, reusing the
// SAME per-site kernel compute_nmix_site (nmix_kernel.h). The log-Cholesky unpack
// mirrors .ms_ocs_chol_unpack; the result is byte-identical to the former loop.
// Output is [M x (n_species * n_sites)] with the per-species blocks contiguous.

#include <Rcpp.h>
#include <vector>
#include "nmix_kernel.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {
// Unpack a log-Cholesky coordinate vector to the lower-triangular factor C
// (column-major storage), matching .ms_ocs_chol_unpack: per column j the
// log-diagonal exp() then the sub-diagonal entries.
inline void chol_unpack(const double* vec, int P, std::vector<double>& C) {
  C.assign((std::size_t) P * P, 0.0);
  int pos = 0;
  for (int j = 0; j < P; ++j) {
    C[(std::size_t) j * P + j] = std::exp(vec[pos++]);
    for (int i = j + 1; i < P; ++i) C[(std::size_t) j * P + i] = vec[pos++];
  }
}
// b = C z for a lower-triangular C (column-major [P x P]).
inline void chol_apply(const std::vector<double>& C, int P, const double* z,
                       double* b) {
  for (int i = 0; i < P; ++i) {
    double acc = 0.0;
    for (int j = 0; j <= i; ++j) acc += C[(std::size_t) j * P + i] * z[j];
    b[i] = acc;
  }
}
}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_ms_nmix_ploglik_batch(
    Rcpp::IntegerVector y,           // [n_obs] longform counts
    Rcpp::IntegerVector species_idx, // [n_obs] 1-based
    Rcpp::IntegerVector site_idx,    // [n_obs] 1-based
    Rcpp::NumericMatrix X_p,         // [n_obs x p_p] per-obs detection design
    Rcpp::NumericMatrix X_lambda,    // [n_sites x p_lam] shared abundance design
    Rcpp::NumericMatrix draws,       // [M x total]
    int mu_off, int b_off, int chol_lam_off, int chol_p_off, int chol_logr_off,
    int p_lam, int p_p, int n_species, int n_sites,
    bool is_nb, int K_max, int n_threads
) {
  const int M = draws.nrow();
  const int n_obs = y.size();
  const int P = p_lam + p_p + (is_nb ? 1 : 0);

  // Group observation rows by (species, site): obs[(s)*n_sites + site].
  std::vector<std::vector<int>> obs_by(
      (std::size_t) n_species * n_sites);
  for (int o = 0; o < n_obs; ++o) {
    int s = species_idx[o] - 1, si = site_idx[o] - 1;
    if (s < 0 || s >= n_species) Rcpp::stop("species_idx out of range.");
    if (si < 0 || si >= n_sites) Rcpp::stop("site_idx out of range.");
    obs_by[(std::size_t) s * n_sites + si].push_back(o);
  }

  Rcpp::NumericMatrix out(M, (std::size_t) n_species * n_sites);
  const double* pd = draws.begin();
  const int* py = y.begin();
  double* pout = out.begin();

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
#endif
  {
    std::vector<double> C_lam, C_p, b(P), coef_lam(p_lam), coef_p(p_p);
    std::vector<int> y_site; std::vector<double> ep_site;
    std::vector<double> eta_lambda(n_sites);
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int m = 0; m < M; ++m) {
      auto dr = [&](int off) { return pd[(std::size_t) off * M + m]; };
      // Cholesky factors for this draw.
      std::vector<double> cl_coord(p_lam * (p_lam + 1) / 2),
                          cp_coord(p_p * (p_p + 1) / 2);
      for (std::size_t k = 0; k < cl_coord.size(); ++k) cl_coord[k] = dr(chol_lam_off + (int) k);
      for (std::size_t k = 0; k < cp_coord.size(); ++k) cp_coord[k] = dr(chol_p_off + (int) k);
      chol_unpack(cl_coord.data(), p_lam, C_lam);
      chol_unpack(cp_coord.data(), p_p, C_p);
      double C_lr = is_nb ? std::exp(dr(chol_logr_off)) : 0.0;

      for (int s = 0; s < n_species; ++s) {
        // Reconstruct b_s = C z_s per arm; add community mean mu.
        int zb = b_off + s * P;
        std::vector<double> z_lam(p_lam), z_p(p_p);
        for (int k = 0; k < p_lam; ++k) z_lam[k] = dr(zb + k);
        for (int k = 0; k < p_p; ++k)   z_p[k]   = dr(zb + p_lam + k);
        std::vector<double> b_lam(p_lam), b_p(p_p);
        chol_apply(C_lam, p_lam, z_lam.data(), b_lam.data());
        chol_apply(C_p, p_p, z_p.data(), b_p.data());
        for (int k = 0; k < p_lam; ++k) coef_lam[k] = dr(mu_off + k) + b_lam[k];
        for (int k = 0; k < p_p; ++k)   coef_p[k]   = dr(mu_off + p_lam + k) + b_p[k];
        double r = is_nb
          ? std::exp(dr(mu_off + P - 1) + C_lr * dr(zb + P - 1))
          : R_PosInf;

        // eta_lambda over all sites (shared design).
        for (int si = 0; si < n_sites; ++si) {
          double e = 0.0;
          for (int k = 0; k < p_lam; ++k) e += X_lambda((std::size_t) si, k) * coef_lam[k];
          eta_lambda[si] = e;
        }
        for (int si = 0; si < n_sites; ++si) {
          const std::vector<int>& obs = obs_by[(std::size_t) s * n_sites + si];
          const int J = (int) obs.size();
          double val = 0.0;
          if (J > 0) {
            y_site.resize(J); ep_site.resize(J);
            for (int j = 0; j < J; ++j) {
              y_site[j] = py[obs[j]];
              double ep = 0.0;
              for (int k = 0; k < p_p; ++k) ep += X_p((std::size_t) obs[j], k) * coef_p[k];
              ep_site[j] = ep;
            }
            val = tulpaObs::compute_nmix_site(
              y_site.data(), ep_site.data(), J, eta_lambda[si], K_max, r).log_lik;
          }
          // out column = s*n_sites + si, row m: column-major index col*M + m.
          std::size_t col = (std::size_t) s * n_sites + si;
          pout[col * M + m] = val;
        }
      }
    }
  }
  return out;
}
