// dyn_abun_laplace.cpp
// Total marginal log-likelihood + per-site gradients for the Dail-Madsen open
// N-mixture (dyn_abun_kernel.h) -- the per-site forward surface the R Laplace
// optimiser (analytic-gradient BFGS), the NUTS oracle, and the WAIC pointwise
// log-likelihood read. The four arms are site-level, so the eta inputs are one
// value per site and the returned gradients are per-site eta-gradients the caller
// sandwiches with each arm's design matrix. `y` is laid out site-major then
// season then visit (index ((i*T)+t)*J + j), with -1 marking a missing visit.

#include <Rcpp.h>
#include <vector>
#include "dyn_abun_kernel.h"

// [[Rcpp::export]]
Rcpp::List cpp_dyn_abun_total_log_lik(
    Rcpp::IntegerVector y, int n_sites, int T, int J, int K,
    Rcpp::NumericVector eta_lambda, Rcpp::NumericVector eta_p,
    Rcpp::NumericVector eta_omega, Rcpp::NumericVector eta_gamma
) {
    if ((int)y.size() != n_sites * T * J)
        Rcpp::stop("y length %d != n_sites*T*J = %d.", (int)y.size(), n_sites * T * J);
    if (eta_lambda.size() != n_sites || eta_p.size() != n_sites ||
        eta_omega.size() != n_sites || eta_gamma.size() != n_sites)
        Rcpp::stop("all eta vectors must have length n_sites.");

    double total = 0.0;
    Rcpp::NumericVector log_lik_site(n_sites), grad_eta_lambda(n_sites),
        grad_eta_p(n_sites), grad_eta_omega(n_sites), grad_eta_gamma(n_sites),
        mean_N1(n_sites);
    int n_inadmissible = 0;
    const int* yp = y.begin();
    for (int i = 0; i < n_sites; ++i) {
        tulpaObs::DynAbunSiteResult r = tulpaObs::compute_dyn_abun_site(
            yp + (std::size_t)i * T * J, T, J, K,
            eta_lambda[i], eta_p[i], eta_omega[i], eta_gamma[i]);
        if (!R_finite(r.log_lik)) ++n_inadmissible;
        total += r.log_lik;
        log_lik_site[i]    = r.log_lik;
        grad_eta_lambda[i] = r.grad_eta_lambda;
        grad_eta_p[i]      = r.grad_eta_p;
        grad_eta_omega[i]  = r.grad_eta_omega;
        grad_eta_gamma[i]  = r.grad_eta_gamma;
        mean_N1[i]         = r.mean_N1;
    }
    return Rcpp::List::create(
        Rcpp::Named("log_lik")          = total,
        Rcpp::Named("log_lik_site")     = log_lik_site,
        Rcpp::Named("grad_eta_lambda")  = grad_eta_lambda,
        Rcpp::Named("grad_eta_p")       = grad_eta_p,
        Rcpp::Named("grad_eta_omega")   = grad_eta_omega,
        Rcpp::Named("grad_eta_gamma")   = grad_eta_gamma,
        Rcpp::Named("mean_N1")          = mean_N1,
        Rcpp::Named("n_inadmissible")   = n_inadmissible);
}
