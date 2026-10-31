// marginal_count_laplace.h
// Shared non-spatial fixed-effects Laplace driver for count-marginal abundance
// families whose per-site likelihood sums a latent abundance N out in closed
// form with a Poisson OR negative-binomial mixing distribution and a logit
// detection arm: the Royle (2004) N-mixture (nmix_kernel.h) and the
// removal-sampling family (removal_kernel.h).
//
// Both families expose the SAME per-site kernel interface -- a callable
//   NMixSiteResult kern(const int* y, const double* eta_p, int n_visits,
//                       double eta_lambda, int K_max, double r)
// filling the abundance score / Fisher, the per-visit detection score / Fisher,
// the abundance score weight (1 Poisson, 1-q NB), the var_N rank-1 coupling, and
// (NB) the theta = log r dispersion row/col. Everything the driver does --
// inner Newton on beta = (beta_lambda, beta_p) against the marginal observed
// Fisher (complete-data Fisher fallback when not PSD), the outer profile-score
// step on theta, and the final joint observed-information Hessian -- is
// identical for every such family, so it is compiled once
// (marginal_count_laplace.cpp) and each family passes its per-site kernel as a
// CountKernelRef. See nmix_laplace.cpp / removal_laplace.cpp for the callers
// and nmix_kernel.h for the curvature derivation.

#ifndef TULPAOBS_MARGINAL_COUNT_LAPLACE_H
#define TULPAOBS_MARGINAL_COUNT_LAPLACE_H

#include "nmix_kernel.h"    // NMixSiteResult
#include "function_ref.h"
#include <Rcpp.h>

namespace tulpaObs {
namespace mcl {

// The family's per-site kernel: (y, eta_p, n_visits, eta_lambda, K_max, r).
using CountKernelRef = FunctionRef<NMixSiteResult(const int*, const double*, int,
                                                  double, int, double)>;

// The full fixed-effects Laplace fit. `kern` is the family's per-site kernel.
Rcpp::List marginal_count_laplace_fixed(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::NumericMatrix X_p_R,
    Rcpp::NumericVector beta_lambda_init,
    Rcpp::NumericVector beta_p_init,
    int K_max,
    int max_iter,
    double tol,
    bool verbose,
    bool nb,
    double log_r_init,
    double theta_max,
    CountKernelRef kern
);

}  // namespace mcl
}  // namespace tulpaObs

#endif  // TULPAOBS_MARGINAL_COUNT_LAPLACE_H
