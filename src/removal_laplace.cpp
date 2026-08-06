// removal_laplace.cpp
// Non-spatial fixed-effects Laplace fit for the removal-sampling abundance
// model with a Poisson OR negative-binomial abundance mixing distribution.
//
// Same shared count-marginal Laplace driver (marginal_count_laplace.h) the
// N-mixture uses -- inner Newton on (beta_lambda, beta_p), outer profile step on
// theta = log r under NB, final joint observed-information inverse -- with the
// removal per-site kernel (removal_kernel.h), whose detection arm sees the
// depleting available count A_k = N - sum_{l<k} y_l rather than the full N.
// `site_idx` groups visits by site preserving input (pass) order, so the kernel
// receives each site's removals in pass order; callers must lay the passes out
// in order with no gaps.

#include "removal_kernel.h"
#include "marginal_count_laplace.h"
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List cpp_removal_laplace_fixed(
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
    double theta_max
) {
    return tulpaObs::mcl::marginal_count_laplace_fixed(
        y, site_idx, X_lambda_R, X_p_R, beta_lambda_init, beta_p_init,
        K_max, max_iter, tol, verbose, nb, log_r_init, theta_max,
        [](const int* yy, const double* ep, int J, double el, int Km, double r) {
            return tulpaObs::compute_removal_site(yy, ep, J, el, Km, r);
        });
}

// Total marginal log-likelihood + per-site/per-pass gradients and observed-info
// pieces for the removal model, the per-site-marginal surface the NUTS oracle
// and any composable random-effect integrator read (mirrors
// cpp_nmix_total_log_lik for the N-mixture). Visits are grouped by site
// preserving input (pass) order.
// [[Rcpp::export]]
Rcpp::List cpp_removal_total_log_lik(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::NumericVector eta_p,
    Rcpp::NumericVector eta_lambda,
    int K_max,
    double r
) {
    const int n_obs   = y.size();
    const int n_sites = eta_lambda.size();
    if (site_idx.size() != n_obs)
        Rcpp::stop("site_idx length must match y length (%d vs %d).",
                   (int)site_idx.size(), n_obs);
    if (eta_p.size() != n_obs)
        Rcpp::stop("eta_p length must match y length (%d vs %d).",
                   (int)eta_p.size(), n_obs);
    if (K_max < 0) Rcpp::stop("K_max must be >= 0.");
    if (R_finite(r) && r <= 0.0) Rcpp::stop("r (NB size) must be > 0.");

    const std::vector<std::vector<int>> obs_by_site =
        tulpaObs::count_group_by_site(site_idx, n_sites);

    // Field set, per-site scatter and returned list shared with the N-mixture
    // sweep (nmix_kernel.h, gcol33/tulpaObs#173).
    tulpaObs::CountSweepAccum acc(n_sites, n_obs);

    for (int s = 0; s < n_sites; ++s) {
        const auto& idx = obs_by_site[s];
        const int J = (int)idx.size();
        if (J == 0) {
            acc.empty_site(s, eta_lambda[s]);
            continue;
        }
        std::vector<int>    y_site(J);
        std::vector<double> eta_p_site(J);
        for (int j = 0; j < J; ++j) {
            y_site[j]     = y[idx[j]];
            eta_p_site[j] = eta_p[idx[j]];
        }
        tulpaObs::NMixSiteResult res = tulpaObs::compute_removal_site(
            y_site.data(), eta_p_site.data(), J, eta_lambda[s], K_max, r);
        acc.scatter(s, idx, res);
    }

    return acc.result();
}
