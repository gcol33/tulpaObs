// dyn_abun_ploglik.cpp
// Batched pointwise log-likelihood for the open-population N-mixture (dyn_abun)
// family: the per-draw loop that R (.tobs_ploglik_dyn_abun) ran around the
// per-site HMM forward marginal now runs in C++, parallel over draws, reusing the
// SAME per-site kernel compute_dyn_abun_site (dyn_abun_kernel.h) the fit and the
// single-draw cpp_dyn_abun_total_log_lik use. The four arms (lambda, p, omega,
// gamma) are site-level [S x n_sites] predictors from BLAS in R; each site's
// count block is y_flat + i * T * J. eta_logr = 0 mirrors the former R loop
// (which called eval_beta without the log r argument). Byte-identical to it.

#include <Rcpp.h>
#include <vector>
#include "dyn_abun_kernel.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dyn_abun_ploglik_batch(
    Rcpp::IntegerVector y_flat,     // [n_sites * T * J], site-major
    int n_sites, int T, int J, int K,
    Rcpp::NumericMatrix eta_lambda, // [S x n_sites]
    Rcpp::NumericMatrix eta_p,
    Rcpp::NumericMatrix eta_omega,
    Rcpp::NumericMatrix eta_gamma,
    bool use_nb,
    Rcpp::NumericVector eta_logr,   // [S] (0 in the WAIC path, mirroring R)
    int n_threads
) {
  const int S = eta_lambda.nrow();
  if (eta_lambda.ncol() != n_sites) Rcpp::stop("eta_lambda must be [S x n_sites].");
  if (eta_logr.size() != S) Rcpp::stop("eta_logr must be length S.");
  const std::size_t site_stride = (std::size_t) T * J;

  Rcpp::NumericMatrix ll(S, n_sites);
  const int* py = y_flat.begin();
  const double* pel = eta_lambda.begin();
  const double* pep = eta_p.begin();
  const double* peo = eta_omega.begin();
  const double* peg = eta_gamma.begin();
  const double* plr = eta_logr.begin();
  double* pll = ll.begin();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads > 0 ? n_threads : 1)
#endif
  for (int d = 0; d < S; ++d) {
    double elogr = plr[d];
    for (int i = 0; i < n_sites; ++i) {
      std::size_t off = (std::size_t) i * S + d;
      tulpaObs::DynAbunSiteResult r = tulpaObs::compute_dyn_abun_site(
        py + (std::size_t) i * site_stride, T, J, K,
        pel[off], pep[off], peo[off], peg[off], use_nb, elogr);
      pll[off] = r.log_lik;
    }
  }
  return ll;
}
