// nmix_ploglik.cpp
// Batched pointwise log-likelihood for the N-mixture family: the per-draw loop
// that R (.tobs_ploglik_nmix) ran around the per-site Royle marginal now runs in
// C++, parallel over draws, reusing the SAME per-site kernel compute_nmix_site
// (nmix_kernel.h) the fit and the single-draw cpp_nmix_total_log_lik use -- one
// source of truth for the marginal. The linear predictors arrive as [S x n_obs]
// / [S x n_sites] matrices built by BLAS in R; the kernel gathers each site's
// visits and evaluates log L_i per draw. Byte-identical to the former R loop
// (same kernel), just without the S R-level calls.

#include <Rcpp.h>
#include <vector>
#include "nmix_kernel.h"
#include "removal_kernel.h"
#include "fp_occu_kernel.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {
// Group observation rows by 1-based site index, input order preserved.
std::vector<std::vector<int>> group_by_site(const int* site_idx, int n_obs,
                                            int n_sites) {
  std::vector<std::vector<int>> obs(n_sites);
  for (int o = 0; o < n_obs; ++o) {
    int s = site_idx[o] - 1;
    if (s < 0 || s >= n_sites) Rcpp::stop("site_idx out of range.");
    obs[s].push_back(o);
  }
  return obs;
}
}  // namespace

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_nmix_ploglik_batch(
    Rcpp::IntegerVector y,          // [n_obs]
    Rcpp::IntegerVector site_idx,   // [n_obs], 1-based
    Rcpp::NumericMatrix eta_p,      // [S x n_obs]
    Rcpp::NumericMatrix eta_lambda, // [S x n_sites]
    int K_max,
    Rcpp::NumericVector r_vec,      // [S] (Inf for Poisson)
    int n_threads
) {
  const int S = eta_lambda.nrow();
  const int n_sites = eta_lambda.ncol();
  const int n_obs = y.size();
  if (eta_p.nrow() != S || eta_p.ncol() != n_obs)
    Rcpp::stop("eta_p must be [S x n_obs].");
  if (r_vec.size() != S) Rcpp::stop("r_vec must be length S.");

  // Group observation rows by site once (input order preserved).
  std::vector<std::vector<int>> obs_by_site =
    group_by_site(site_idx.begin(), n_obs, n_sites);

  Rcpp::NumericMatrix ll(S, n_sites);
  const double* pep = eta_p.begin();       // column-major [S x n_obs]
  const double* pel = eta_lambda.begin();  // column-major [S x n_sites]
  const int* py = y.begin();
  double* pll = ll.begin();

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
#endif
  {
    std::vector<int> y_site;
    std::vector<double> ep_site;
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int d = 0; d < S; ++d) {
      double r = r_vec[d];
      for (int s = 0; s < n_sites; ++s) {
        const std::vector<int>& obs = obs_by_site[s];
        const int J = (int) obs.size();
        double val;
        if (J == 0) {
          val = 0.0;                          // no visits -> marginal 1
        } else {
          y_site.resize(J); ep_site.resize(J);
          for (int j = 0; j < J; ++j) {
            y_site[j]  = py[obs[j]];
            ep_site[j] = pep[(std::size_t) obs[j] * S + d];
          }
          double el = pel[(std::size_t) s * S + d];
          tulpaObs::NMixSiteResult res = tulpaObs::compute_nmix_site(
            y_site.data(), ep_site.data(), J, el, K_max, r);
          val = res.log_lik;
        }
        pll[(std::size_t) s * S + d] = val;
      }
    }
  }
  return ll;
}

// Removal sampling: same site grouping as N-mixture, per-site depleting-binomial
// removal marginal (compute_removal_site, removal_kernel.h). y is the per-pass
// removal counts; K_max here is site-total-based (set by the R caller).
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_removal_ploglik_batch(
    Rcpp::IntegerVector y, Rcpp::IntegerVector site_idx,
    Rcpp::NumericMatrix eta_p, Rcpp::NumericMatrix eta_lambda,
    int K_max, Rcpp::NumericVector r_vec, int n_threads
) {
  const int S = eta_lambda.nrow(), n_sites = eta_lambda.ncol(), n_obs = y.size();
  if (eta_p.nrow() != S || eta_p.ncol() != n_obs)
    Rcpp::stop("eta_p must be [S x n_obs].");
  if (r_vec.size() != S) Rcpp::stop("r_vec must be length S.");
  std::vector<std::vector<int>> obs_by_site =
    group_by_site(site_idx.begin(), n_obs, n_sites);
  Rcpp::NumericMatrix ll(S, n_sites);
  const double* pep = eta_p.begin(); const double* pel = eta_lambda.begin();
  const int* py = y.begin(); double* pll = ll.begin();
#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
#endif
  {
    std::vector<int> y_site; std::vector<double> ep_site;
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int d = 0; d < S; ++d) {
      double r = r_vec[d];
      for (int s = 0; s < n_sites; ++s) {
        const std::vector<int>& obs = obs_by_site[s];
        const int J = (int) obs.size();
        double val = 0.0;
        if (J > 0) {
          y_site.resize(J); ep_site.resize(J);
          for (int j = 0; j < J; ++j) {
            y_site[j]  = py[obs[j]];
            ep_site[j] = pep[(std::size_t) obs[j] * S + d];
          }
          val = tulpaObs::compute_removal_site(
            y_site.data(), ep_site.data(), J,
            pel[(std::size_t) s * S + d], K_max, r).log_lik;
        }
        pll[(std::size_t) s * S + d] = val;
      }
    }
  }
  return ll;
}

// False-positive occupancy (Miller et al. 2011 multistate): per-site marginal
// over the latent z (compute_fp_occu_site, fp_occu_kernel.h). Four site-level
// arms (psi, p11, p10, b), each [S x n_sites]; no latent-count sum, so no K_max.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_fp_occu_ploglik_batch(
    Rcpp::IntegerVector y, Rcpp::IntegerVector site_idx,
    Rcpp::NumericMatrix eta_psi, Rcpp::NumericMatrix eta_p11,
    Rcpp::NumericMatrix eta_p10, Rcpp::NumericMatrix eta_b, int n_threads
) {
  const int S = eta_psi.nrow(), n_sites = eta_psi.ncol(), n_obs = y.size();
  std::vector<std::vector<int>> obs_by_site =
    group_by_site(site_idx.begin(), n_obs, n_sites);
  Rcpp::NumericMatrix ll(S, n_sites);
  const double* ppsi = eta_psi.begin(); const double* pp11 = eta_p11.begin();
  const double* pp10 = eta_p10.begin(); const double* pb = eta_b.begin();
  const int* py = y.begin(); double* pll = ll.begin();
#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
#endif
  {
    std::vector<int> y_site;
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int d = 0; d < S; ++d) {
      for (int s = 0; s < n_sites; ++s) {
        const std::vector<int>& obs = obs_by_site[s];
        const int J = (int) obs.size();
        double val = 0.0;
        if (J > 0) {
          y_site.resize(J);
          for (int j = 0; j < J; ++j) y_site[j] = py[obs[j]];
          std::size_t off = (std::size_t) s * S + d;
          val = tulpaObs::compute_fp_occu_site(
            y_site.data(), J, ppsi[off], pp11[off], pp10[off], pb[off]).log_lik;
        }
        pll[(std::size_t) s * S + d] = val;
      }
    }
  }
  return ll;
}
