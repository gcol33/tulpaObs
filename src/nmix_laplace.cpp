// nmix_laplace.cpp
// Non-spatial fixed-effects Laplace fit for the Royle (2004) N-mixture model
// with a Poisson OR negative-binomial abundance mixing distribution.
//
// The fit -- inner Newton on beta = (beta_lambda, beta_p) against the marginal
// observed Fisher information, the outer profile-score step on theta = log r
// under NB, and the final joint observed-information inverse (Louis 1982; see
// nmix_kernel.h) -- is the shared count-marginal Laplace driver
// (marginal_count_laplace.h). This file is the (single-line) instantiation with
// the N-mixture per-site kernel; the removal-sampling family is the same driver
// with its own kernel (removal_laplace.cpp).

#include "nmix_kernel.h"
#include "marginal_count_laplace.h"
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List cpp_nmix_laplace_fixed(
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
    int headroom = -1   // latent-N states above each site's own max(y); <0 = none
) {
    return tulpaObs::mcl::marginal_count_laplace_fixed(
        y, site_idx, X_lambda_R, X_p_R, beta_lambda_init, beta_p_init,
        K_max, max_iter, tol, verbose, nb, log_r_init, theta_max,
        [headroom](const int* yy, const double* ep, int J, double el, int Km,
                   double r) {
            const tulpaObs::NMixSiteCache c =
                tulpaObs::nmix_precompute_site(yy, J, Km, headroom);
            return tulpaObs::compute_nmix_site_cached(c, ep, el, r);
        });
}
