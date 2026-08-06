// nmix_loglik.cpp
// Rcpp surface for the N-mixture marginal log-likelihood kernel
// (see nmix_kernel.h). Provides a thin pass-through used both for
// validation against `unmarked::pcount` and as the per-iteration evaluator
// the Laplace Newton driver will call (forthcoming).
//
// Data layout (long form):
//   y          int[n_obs]   observed counts, visits stacked across sites
//   site_idx   int[n_obs]   1-based site index for each visit
//   eta_p      double[n_obs] detection logit linear predictor per visit
//   eta_lambda double[n_sites] abundance log linear predictor per site
//
// All exported quantities aggregate site-level kernel outputs:
//   log_lik           = sum_i log L_i
//   log_lik_site      = log L_i                 per site  (grouped-marginal sums)
//   grad_eta_lambda   = (E[N|y_i] - lambda_i)   per site
//   grad_eta_p        = (y_ij - E[N|y_i] p_ij)  per visit
//   info_eta_lambda   = lambda_i                per site  (complete-data Fisher)
//   info_eta_p        = E[N|y_i] p_ij (1-p_ij)  per visit (complete-data Fisher)
//   score_wt_lambda   = N-coeff of s_lambda     per site  (1 Poisson / 1-q NB)
//   mean_N, var_N     posterior moments per site (diagnostics, observed-info path)
//   boundary_weight   posterior mass on N = K_max per site (K_max sanity check)
//
// The per-site marginal observed-information block in the eta coordinates
// (eta_lambda_i, eta_p_{i,1..J}) is assembled from these pieces as
//   B_i = diag(info_eta_lambda_i, info_eta_p_{ij}) - var_N_i * v_i v_i^T,
//   v_i = (-score_wt_lambda_i, p_{i1}, ..., p_{iJ}),   p_{ij} = plogis(eta_p_{ij}),
// i.e. the complete-data Fisher diagonal minus the rank-1 shared-latent
// score covariance (Louis 1982; same coupling assemble_beta_obs_info() in
// nmix_laplace.cpp sandwiches with the design matrices). The lambda<->p cross
// term var_N * score_wt_lambda * p_ij is the both-arm coupling a random-effect
// integrator needs. Under NB the theta = log r pieces (info_theta,
// info_lambda_theta, cov_N_stheta, var_stheta, per site) carry the dispersion
// row/col of the joint observed information when log r is treated as a free
// (global) parameter; they are zero on the Poisson path.

#include "nmix_kernel.h"
#include <Rcpp.h>
#include <vector>

// [[Rcpp::export]]
Rcpp::List cpp_nmix_total_log_lik(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::NumericVector eta_p,
    Rcpp::NumericVector eta_lambda,
    int K_max,
    double r,         // NB size; pass Inf for the Poisson kernel
    int headroom = -1 // latent-N states above each site's own max(y); <0 = none
) {
    const int n_obs = y.size();
    const int n_sites = eta_lambda.size();
    if (site_idx.size() != n_obs) {
        Rcpp::stop("site_idx length must match y length (%d vs %d).",
                   site_idx.size(), n_obs);
    }
    if (eta_p.size() != n_obs) {
        Rcpp::stop("eta_p length must match y length (%d vs %d).",
                   eta_p.size(), n_obs);
    }
    if (K_max < 0) Rcpp::stop("K_max must be >= 0.");
    if (R_finite(r) && r <= 0.0) Rcpp::stop("r (NB size) must be > 0.");

    // Group observations by site (preserves input order within each site).
    const std::vector<std::vector<int>> obs_by_site =
        tulpaObs::count_group_by_site(site_idx, n_sites);

    // The output field set, the per-site scatter and the returned list are
    // shared with the removal sweep (nmix_kernel.h, gcol33/tulpaObs#173), so a
    // field added here reaches every count family.
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
        const tulpaObs::NMixSiteCache cache =
            tulpaObs::nmix_precompute_site(y_site.data(), J, K_max, headroom);
        tulpaObs::NMixSiteResult res = tulpaObs::compute_nmix_site_cached(
            cache, eta_p_site.data(), eta_lambda[s], r
        );
        acc.scatter(s, idx, res);
    }

    return acc.result();
}
