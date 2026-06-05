// fp_occu_laplace.cpp
// Total marginal log-likelihood + per-site gradients for the multistate
// false-positive occupancy model (fp_occu_kernel.h). The per-site marginal
// surface the R Laplace optimiser (optim BFGS over the analytic gradient), the
// NUTS oracle, and the WAIC pointwise log-likelihood read. The arms are
// site-level, so the eta inputs are one value per site and the returned
// gradients are per-site eta-gradients the caller sandwiches with each arm's
// design matrix.

#include <Rcpp.h>
#include <vector>
#include "fp_occu_kernel.h"

// [[Rcpp::export]]
Rcpp::List cpp_fp_occu_total_log_lik(
    Rcpp::IntegerVector y,          // valid detection states (0/1/2), long form
    Rcpp::IntegerVector site_idx,   // 1-based site per observation
    Rcpp::NumericVector eta_psi,    // length n_sites
    Rcpp::NumericVector eta_p11,
    Rcpp::NumericVector eta_p10,
    Rcpp::NumericVector eta_b
) {
    const int n_obs = y.size();
    const int n_sites = eta_psi.size();
    if (site_idx.size() != n_obs)
        Rcpp::stop("site_idx length must match y length (%d vs %d).",
                   (int)site_idx.size(), n_obs);
    if (eta_p11.size() != n_sites || eta_p10.size() != n_sites || eta_b.size() != n_sites)
        Rcpp::stop("all eta vectors must have length n_sites.");

    std::vector<std::vector<int>> obs_by_site(n_sites);
    for (int o = 0; o < n_obs; ++o) {
        const int s = site_idx[o] - 1;
        if (s < 0 || s >= n_sites)
            Rcpp::stop("site_idx[%d] = %d out of range [1, %d].", o + 1, site_idx[o], n_sites);
        obs_by_site[s].push_back(o);
    }

    double total = 0.0;
    Rcpp::NumericVector log_lik_site(n_sites), grad_eta_psi(n_sites),
        grad_eta_p11(n_sites), grad_eta_p10(n_sites), grad_eta_b(n_sites),
        w1(n_sites);
    std::vector<int> y_site;
    for (int s = 0; s < n_sites; ++s) {
        const std::vector<int>& idx = obs_by_site[s];
        const int J = (int)idx.size();
        y_site.resize(J);
        for (int j = 0; j < J; ++j) y_site[j] = y[idx[j]];
        tulpaObs::FpOccuSiteResult r = tulpaObs::compute_fp_occu_site(
            y_site.data(), J, eta_psi[s], eta_p11[s], eta_p10[s], eta_b[s]);
        total += r.log_lik;
        log_lik_site[s] = r.log_lik;
        grad_eta_psi[s] = r.grad_eta_psi;
        grad_eta_p11[s] = r.grad_eta_p11;
        grad_eta_p10[s] = r.grad_eta_p10;
        grad_eta_b[s]   = r.grad_eta_b;
        w1[s] = r.w1;
    }
    return Rcpp::List::create(
        Rcpp::Named("log_lik")       = total,
        Rcpp::Named("log_lik_site")  = log_lik_site,
        Rcpp::Named("grad_eta_psi")  = grad_eta_psi,
        Rcpp::Named("grad_eta_p11")  = grad_eta_p11,
        Rcpp::Named("grad_eta_p10")  = grad_eta_p10,
        Rcpp::Named("grad_eta_b")    = grad_eta_b,
        Rcpp::Named("w1")            = w1);
}
