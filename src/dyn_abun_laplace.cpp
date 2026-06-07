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

// `use_nb` switches the initial abundance to negative-binomial NB(mean=lambda,
// size = exp(eta_logr)); `eta_logr` is a single shared dispersion coordinate.
// `grad_eta_logr` is summed across sites (the score for the one log r parameter).
// [[Rcpp::export]]
Rcpp::List cpp_dyn_abun_total_log_lik(
    Rcpp::IntegerVector y, int n_sites, int T, int J, int K,
    Rcpp::NumericVector eta_lambda, Rcpp::NumericVector eta_p,
    Rcpp::NumericVector eta_omega, Rcpp::NumericVector eta_gamma,
    bool use_nb = false, double eta_logr = 0.0
) {
    if ((int)y.size() != n_sites * T * J)
        Rcpp::stop("y length %d != n_sites*T*J = %d.", (int)y.size(), n_sites * T * J);
    if (eta_lambda.size() != n_sites || eta_p.size() != n_sites ||
        eta_omega.size() != n_sites || eta_gamma.size() != n_sites)
        Rcpp::stop("all eta vectors must have length n_sites.");

    double total = 0.0, grad_eta_logr = 0.0;
    Rcpp::NumericVector log_lik_site(n_sites), grad_eta_lambda(n_sites),
        grad_eta_p(n_sites), grad_eta_omega(n_sites), grad_eta_gamma(n_sites),
        mean_N1(n_sites);
    int n_inadmissible = 0;
    const int* yp = y.begin();
    for (int i = 0; i < n_sites; ++i) {
        tulpaObs::DynAbunSiteResult r = tulpaObs::compute_dyn_abun_site(
            yp + (std::size_t)i * T * J, T, J, K,
            eta_lambda[i], eta_p[i], eta_omega[i], eta_gamma[i], use_nb, eta_logr);
        if (!R_finite(r.log_lik)) ++n_inadmissible;
        total += r.log_lik;
        log_lik_site[i]    = r.log_lik;
        grad_eta_lambda[i] = r.grad_eta_lambda;
        grad_eta_p[i]      = r.grad_eta_p;
        grad_eta_omega[i]  = r.grad_eta_omega;
        grad_eta_gamma[i]  = r.grad_eta_gamma;
        grad_eta_logr     += r.grad_eta_logr;
        mean_N1[i]         = r.mean_N1;
    }
    return Rcpp::List::create(
        Rcpp::Named("log_lik")          = total,
        Rcpp::Named("log_lik_site")     = log_lik_site,
        Rcpp::Named("grad_eta_lambda")  = grad_eta_lambda,
        Rcpp::Named("grad_eta_p")       = grad_eta_p,
        Rcpp::Named("grad_eta_omega")   = grad_eta_omega,
        Rcpp::Named("grad_eta_gamma")   = grad_eta_gamma,
        Rcpp::Named("grad_eta_logr")    = grad_eta_logr,
        Rcpp::Named("mean_N1")          = mean_N1,
        Rcpp::Named("n_inadmissible")   = n_inadmissible);
}

// Per-site conditional-likelihood weights c(n1) = P(all data | N_1 = n1) for the
// grouped random-effect AGHQ path on the initial-abundance arm (tulpaObs#51).
// `site` are 0-based site indices into `y`; eta_p / eta_omega / eta_gamma are the
// (RE-fixed) detection / survival / recruitment predictors aligned with `site`.
// Returns an m x (K+1) matrix of c-weights, one row per requested site. The
// engine precomputes this ONCE per make_site call; the O(K) per-node log marginal
// then comes from cpp_dyn_abun_init_loglik below.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_dyn_abun_init_weights_mat(
    Rcpp::IntegerVector y, int n_sites, int T, int J, int K,
    Rcpp::IntegerVector site, Rcpp::NumericVector eta_p,
    Rcpp::NumericVector eta_omega, Rcpp::NumericVector eta_gamma
) {
    if ((int)y.size() != n_sites * T * J)
        Rcpp::stop("y length %d != n_sites*T*J = %d.", (int)y.size(), n_sites * T * J);
    const int m = site.size();
    if (eta_p.size() != m || eta_omega.size() != m || eta_gamma.size() != m)
        Rcpp::stop("eta vectors must be aligned with `site` (length %d).", m);

    const int S = K + 1;
    Rcpp::NumericMatrix C(m, S);
    std::vector<double> c_row(S);
    const int* yp = y.begin();
    for (int k = 0; k < m; ++k) {
        const int i = site[k];
        if (i < 0 || i >= n_sites) Rcpp::stop("site index out of range.");
        tulpaObs::compute_dyn_abun_init_weights(
            yp + (std::size_t)i * T * J, T, J, K,
            eta_p[k], eta_omega[k], eta_gamma[k], c_row.data());
        for (int n = 0; n < S; ++n) C(k, n) = c_row[n];
    }
    return C;
}

// Per-site log marginal L(eta_lambda) = sum_{n1} pi_{n1}(eta_lambda) c(n1) from
// precomputed c-weights `C` (m x (K+1)) and the per-site initial predictor
// eta_lambda (length m). With `deriv`, also the first and second eta_lambda
// derivatives (for the engine's per-group Newton mode and curvature):
//   pi Poisson: d log pi = (n - lambda),  d2 log pi = -lambda
//   pi NB(mu, r): d log pi = n - mu(n+r)/(r+mu),  d2 log pi = -mu(n+r)r/(r+mu)^2
// and  d1 = L'/L,  d2 = L''/L - (L'/L)^2  with L', L'' the pi'-, pi''-weighted
// sums (pi' = pi * dlogpi, pi'' = pi * (dlogpi^2 + d2logpi)). O(K) per site.
// [[Rcpp::export]]
Rcpp::List cpp_dyn_abun_init_loglik(
    Rcpp::NumericMatrix C, Rcpp::NumericVector eta_lambda,
    bool use_nb = false, double eta_logr = 0.0, bool deriv = true
) {
    const int m = C.nrow(), S = C.ncol();
    if (eta_lambda.size() != m)
        Rcpp::stop("eta_lambda length %d != nrow(C) %d.", (int)eta_lambda.size(), m);
    Rcpp::NumericVector logL(m), d1(m), d2(m);
    const double rr = use_nb ? std::exp(eta_logr) : 0.0;
    for (int k = 0; k < m; ++k) {
        const double el = eta_lambda[k], lam = std::exp(el);
        const double rpm = rr + lam;
        double L = 0.0, L1 = 0.0, L2 = 0.0;
        for (int n = 0; n < S; ++n) {
            double pi_n, g_l = 0.0, d2log = 0.0;
            if (use_nb) {
                const double lpn = R::lgammafn((double)n + rr) - R::lgammafn(rr)
                    - R::lgammafn((double)n + 1.0)
                    + rr * std::log(rr / rpm) + (double)n * std::log(lam / rpm);
                pi_n = std::exp(lpn);
                if (deriv) {
                    g_l = (double)n - lam * ((double)n + rr) / rpm;
                    d2log = -lam * ((double)n + rr) * rr / (rpm * rpm);
                }
            } else {
                const double lpn = -lam + (double)n * el - R::lgammafn((double)n + 1.0);
                pi_n = std::exp(lpn);
                if (deriv) { g_l = (double)n - lam; d2log = -lam; }
            }
            const double a = pi_n * C(k, n);
            L += a;
            if (deriv) { L1 += a * g_l; L2 += a * (g_l * g_l + d2log); }
        }
        if (!(L > 0.0)) {
            logL[k] = -std::numeric_limits<double>::infinity();
            continue;
        }
        logL[k] = std::log(L);
        if (deriv) { const double dd = L1 / L; d1[k] = dd; d2[k] = L2 / L - dd * dd; }
    }
    return Rcpp::List::create(
        Rcpp::Named("logL") = logL,
        Rcpp::Named("d1")   = d1,
        Rcpp::Named("d2")   = d2);
}
