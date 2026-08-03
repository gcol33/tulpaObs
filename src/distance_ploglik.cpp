// distance_ploglik.cpp
// Batched pointwise log-likelihood for the binned distance-sampling family: the
// per-draw loop that R (.tobs_ploglik_distance) ran around the per-site binned-
// multinomial-over-N marginal now runs in C++, parallel over draws, reusing the
// SAME per-site kernel compute_distance_site (distance_kernel.h) and the SAME
// quadrature build dist_build_quad the fit uses. The quad is built once; the
// linear predictors arrive as [S x n_sites] matrices from BLAS in R. Byte-
// identical to the former R loop (same kernel + quad).

#include <Rcpp.h>
#include <vector>
#include "distance_quad.h"
#include "distance_kernel.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_distance_ploglik_batch(
    Rcpp::IntegerMatrix y,          // [n_sites x n_bins]
    SEXP quad_xptr,                 // per-fit quadrature (cpp_distance_build_quad)
    int key, int K_max,
    Rcpp::NumericMatrix eta_lambda, // [S x n_sites]
    Rcpp::NumericMatrix eta_sigma,  // [S x n_sites]
    Rcpp::NumericVector eta_b,      // [S] (hazard shape; 0 otherwise)
    Rcpp::NumericVector r_vec,      // [S] (Inf for Poisson)
    int n_threads,
    int headroom = -1               // gcol33/tulpaObs#168: per-site K_hi cap
) {
  const int S = eta_lambda.nrow();
  const int n_sites = eta_lambda.ncol();
  const int n_bins = y.ncol();
  if (y.nrow() != n_sites) Rcpp::stop("y must have n_sites rows.");
  if (eta_sigma.nrow() != S || eta_sigma.ncol() != n_sites)
    Rcpp::stop("eta_sigma must be [S x n_sites].");
  if (eta_b.size() != S || r_vec.size() != S)
    Rcpp::stop("eta_b / r_vec must be length S.");

  Rcpp::XPtr<tulpaObs::DistQuad> quad_ptr = tulpaObs::dist_quad_from_xptr(quad_xptr);
  const tulpaObs::DistQuad& quad = *quad_ptr;
  if (quad.n_bins != n_bins) Rcpp::stop("quad_xptr's bin count must equal ncol(y).");

  // Per-site bin counts (draw-invariant), gathered once.
  std::vector<std::vector<int>> y_by_site(n_sites, std::vector<int>(n_bins));
  for (int s = 0; s < n_sites; ++s)
    for (int b = 0; b < n_bins; ++b) y_by_site[s][b] = y(s, b);

  Rcpp::NumericMatrix ll(S, n_sites);
  const double* pel = eta_lambda.begin();
  const double* pes = eta_sigma.begin();
  const double* peb = eta_b.begin();
  const double* prv = r_vec.begin();
  double* pll = ll.begin();

  // Shared, read-only: the K_max-indexed combinatorial table (gcol33/tulpaObs#167)
  // is the SAME for every (site, draw), so it is built once and reused across
  // the whole S x n_sites sweep instead of every one of those calls repeating
  // its own O(K_max) run of R::lgammafn().
  const std::vector<double> comb_table = tulpaObs::dist_build_comb_table(K_max);

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
  {
    tulpaObs::DistScratch scratch;
    #pragma omp for schedule(static)
    for (int d = 0; d < S; ++d) {
      double eb = peb[d], r = prv[d];
      for (int s = 0; s < n_sites; ++s) {
        std::size_t off = (std::size_t) s * S + d;
        // Only log_lik is read below, so value_only=true skips the detection-arm
        // gradient/Fisher block and its five per-bin derivative vectors entirely
        // (gcol33/tulpaObs#164) -- this call previously computed and discarded them.
        double val = tulpaObs::compute_distance_site(
          y_by_site[s].data(), n_bins, pel[off], pes[off], eb, key, quad,
          K_max, r, /*value_only=*/true, &comb_table, &scratch,
          headroom).log_lik;
        pll[off] = val;
      }
    }
  }
#else
  tulpaObs::DistScratch scratch;
  for (int d = 0; d < S; ++d) {
    double eb = peb[d], r = prv[d];
    for (int s = 0; s < n_sites; ++s) {
      std::size_t off = (std::size_t) s * S + d;
      double val = tulpaObs::compute_distance_site(
        y_by_site[s].data(), n_bins, pel[off], pes[off], eb, key, quad,
        K_max, r, /*value_only=*/true, &comb_table, &scratch,
        headroom).log_lik;
      pll[off] = val;
    }
  }
#endif
  return ll;
}
