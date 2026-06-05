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

    std::vector<std::vector<int>> obs_by_site(n_sites);
    for (int o = 0; o < n_obs; ++o) {
        const int s = site_idx[o] - 1;
        if (s < 0 || s >= n_sites)
            Rcpp::stop("site_idx[%d] = %d out of range [1, %d].",
                       o + 1, site_idx[o], n_sites);
        obs_by_site[s].push_back(o);
    }

    double total_log_lik = 0.0, total_grad_theta = 0.0;
    Rcpp::NumericVector log_lik_site(n_sites), grad_eta_lambda(n_sites),
        info_eta_lambda(n_sites), score_wt_lambda(n_sites), mean_N(n_sites),
        var_N(n_sites), boundary_weight(n_sites), info_theta(n_sites),
        info_lambda_theta(n_sites), cov_N_stheta(n_sites), var_stheta(n_sites);
    Rcpp::NumericVector grad_eta_p(n_obs), info_eta_p(n_obs);
    int n_K_inadmissible = 0;

    for (int s = 0; s < n_sites; ++s) {
        const auto& idx = obs_by_site[s];
        const int J = (int)idx.size();
        if (J == 0) {
            log_lik_site[s] = 0.0; grad_eta_lambda[s] = 0.0;
            info_eta_lambda[s] = 0.0; score_wt_lambda[s] = 1.0;
            mean_N[s] = std::exp(eta_lambda[s]);
            var_N[s]  = std::exp(eta_lambda[s]);
            boundary_weight[s] = 0.0;
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
        if (!R_finite(res.log_lik)) ++n_K_inadmissible;
        total_log_lik    += res.log_lik;
        total_grad_theta += res.grad_theta;
        log_lik_site[s]    = res.log_lik;
        grad_eta_lambda[s] = res.grad_eta_lambda;
        info_eta_lambda[s] = res.info_eta_lambda;
        score_wt_lambda[s] = res.score_wt_lambda;
        mean_N[s] = res.mean_N; var_N[s] = res.var_N;
        boundary_weight[s] = res.boundary_weight;
        info_theta[s]        = res.info_theta;
        info_lambda_theta[s] = res.info_lambda_theta;
        cov_N_stheta[s]      = res.cov_N_stheta;
        var_stheta[s]        = res.var_stheta;
        for (int j = 0; j < J; ++j) {
            grad_eta_p[idx[j]] = res.grad_eta_p[j];
            info_eta_p[idx[j]] = res.info_eta_p[j];
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("log_lik")           = total_log_lik,
        Rcpp::Named("log_lik_site")      = log_lik_site,
        Rcpp::Named("grad_eta_lambda")   = grad_eta_lambda,
        Rcpp::Named("grad_eta_p")        = grad_eta_p,
        Rcpp::Named("grad_theta")        = total_grad_theta,
        Rcpp::Named("info_eta_lambda")   = info_eta_lambda,
        Rcpp::Named("info_eta_p")        = info_eta_p,
        Rcpp::Named("score_wt_lambda")   = score_wt_lambda,
        Rcpp::Named("mean_N")            = mean_N,
        Rcpp::Named("var_N")             = var_N,
        Rcpp::Named("boundary_weight")   = boundary_weight,
        Rcpp::Named("info_theta")        = info_theta,
        Rcpp::Named("info_lambda_theta") = info_lambda_theta,
        Rcpp::Named("cov_N_stheta")      = cov_N_stheta,
        Rcpp::Named("var_stheta")        = var_stheta,
        Rcpp::Named("n_K_inadmissible")  = n_K_inadmissible
    );
}
