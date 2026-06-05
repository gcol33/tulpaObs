// abun_nuts.cpp
// C++ joint log-posterior + gradient for the single-species N-mixture (abun())
// NUTS target. The R reference .tobs_abun_nuts_logpost (R/abun_nuts.R) is the
// oracle; this port mirrors it and is cross-checked against it before driving
// tulpa's NUTS engine.
//
// theta = (beta_lambda [p_lam], beta_p [p_p], [log_r under NB]).
//   log p(theta|y) = sum_i log m_i(theta)  - 0.5||beta||^2/sigma.beta^2
//                    [ - 0.5 log_r^2 / sigma.logr^2 ]   (NB only)
// where m_i is the Royle (2004) per-site marginal. compute_nmix_site()
// (nmix_kernel.h) already returns log m_i, grad_eta_lambda, grad_eta_p, and (NB)
// grad_theta = d log m_i / d log_r, so the coefficient gradient is the
// design-sandwiched eta-gradient -- the same arithmetic as the R target, no new
// likelihood math. No latent field, no random effect: the parameter vector is the
// flat coefficient block, so the FullGradFn is a plain coefficient loop plus the
// Gaussian priors (far simpler than the spatial-factor community target).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
#include "nmix_kernel.h"

using namespace Rcpp;

namespace tulpaObs {

// Marshalled per-fit data for the N-mixture NUTS target. The long form is the
// same .tobs_build_abun() produces: y / site_idx in site-major order, X_lambda
// site-level (n_sites rows), X_p long form (n_obs rows). Per-site visit row lists
// are built once at construction so the eval loop gathers each site's counts /
// detection etas contiguously for compute_nmix_site().
struct NMixNutsData {
    int n_sites = 0, n_obs = 0, p_lam = 0, p_p = 0, K_max = 0;
    bool is_nb = false;
    std::vector<int> y;                       // length n_obs (long form)
    NumericMatrix X_lambda;                   // n_sites x p_lam
    NumericMatrix X_p;                        // n_obs   x p_p
    std::vector<std::vector<int>> obs_by_site; // per-site row indices into y / X_p
    int total = 0;
};

inline NMixNutsData abun_nuts_build_data(const Rcpp::List& spec) {
    NMixNutsData d;
    IntegerVector y        = spec["y"];
    IntegerVector site_idx = spec["site_idx"];   // 1-based
    d.X_lambda = Rcpp::as<NumericMatrix>(spec["X_lambda"]);
    d.X_p      = Rcpp::as<NumericMatrix>(spec["X_p"]);
    d.n_sites  = Rcpp::as<int>(spec["n_sites"]);
    d.K_max    = Rcpp::as<int>(spec["K_max"]);
    d.is_nb    = Rcpp::as<bool>(spec["is_nb"]);
    d.n_obs    = y.size();
    d.p_lam    = d.X_lambda.ncol();
    d.p_p      = d.X_p.ncol();
    d.y.assign(y.begin(), y.end());
    d.obs_by_site.assign(d.n_sites, std::vector<int>());
    for (int o = 0; o < d.n_obs; ++o) {
        const int i = site_idx[o] - 1;
        if (i < 0 || i >= d.n_sites)
            Rcpp::stop("site_idx out of range in abun_nuts_build_data");
        d.obs_by_site[i].push_back(o);
    }
    d.total = d.p_lam + d.p_p + (d.is_nb ? 1 : 0);
    return d;
}

// Joint log-posterior + gradient over (beta_lambda, beta_p[, log_r]). Returns the
// log-posterior (NUTS maximises, no negation) and writes the gradient into `grad`
// (length d.total). Mirrors .tobs_abun_nuts_logpost byte-for-byte.
inline double abun_nuts_eval(const NMixNutsData& d, const double* theta,
                             double sigma_beta, double sigma_logr, double* grad) {
    const int p_lam = d.p_lam, p_p = d.p_p;
    const double r = d.is_nb
        ? std::exp(theta[p_lam + p_p])
        : std::numeric_limits<double>::infinity();
    for (int j = 0; j < d.total; ++j) grad[j] = 0.0;

    // Per-observation detection eta (long form) at the current beta_p.
    std::vector<double> eta_p_all(d.n_obs, 0.0);
    for (int o = 0; o < d.n_obs; ++o) {
        double e = 0.0;
        for (int k = 0; k < p_p; ++k) e += d.X_p(o, k) * theta[p_lam + k];
        eta_p_all[o] = e;
    }

    double lp = 0.0;
    std::vector<int>    y_site;
    std::vector<double> eta_p_site;
    for (int i = 0; i < d.n_sites; ++i) {
        const std::vector<int>& obs = d.obs_by_site[i];
        const int J = (int) obs.size();
        // eta_lambda_i = X_lambda[i,] . beta_lambda
        double eta_lambda = 0.0;
        for (int k = 0; k < p_lam; ++k) eta_lambda += d.X_lambda(i, k) * theta[k];
        // gather this site's counts + detection etas contiguously
        y_site.resize(J); eta_p_site.resize(J);
        for (int jj = 0; jj < J; ++jj) {
            y_site[jj]     = d.y[obs[jj]];
            eta_p_site[jj] = eta_p_all[obs[jj]];
        }
        const NMixSiteResult res = compute_nmix_site(
            y_site.data(), eta_p_site.data(), J, eta_lambda, d.K_max, r);
        lp += res.log_lik;
        // grad_beta_lambda += grad_eta_lambda * X_lambda[i,]
        for (int k = 0; k < p_lam; ++k)
            grad[k] += res.grad_eta_lambda * d.X_lambda(i, k);
        // grad_beta_p += sum_visits grad_eta_p[jj] * X_p[obs,]
        for (int jj = 0; jj < J; ++jj) {
            const double ge = res.grad_eta_p[jj];
            const int o = obs[jj];
            for (int k = 0; k < p_p; ++k)
                grad[p_lam + k] += ge * d.X_p(o, k);
        }
        // grad_log_r += grad_theta (theta = log r)
        if (d.is_nb) grad[p_lam + p_p] += res.grad_theta;
    }

    // Weak Gaussian coefficient priors N(0, sigma_beta^2).
    const double ib2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < p_lam + p_p; ++k) {
        lp        -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k]   -= ib2 * theta[k];
    }
    // NB: N(0, sigma_logr^2) prior on log_r.
    if (d.is_nb) {
        const double lr = theta[p_lam + p_p];
        const double ilr2 = 1.0 / (sigma_logr * sigma_logr);
        lp -= 0.5 * ilr2 * lr * lr;
        grad[p_lam + p_p] -= ilr2 * lr;
    }
    return lp;
}

// NUTS model carrying the marshalled data + prior scales; the FullGradFn reaches
// it through ModelData.model_response_data.
struct NMixNutsModel {
    NMixNutsData d;
    double sigma_beta = 10.0, sigma_logr = 1.5;
};

// FullGradFn: log-posterior + gradient over the entire parameter vector.
inline void abun_nuts_full_grad(const std::vector<double>& params,
                                const tulpa::ModelData& data,
                                const tulpa::ParamLayout& /*layout*/,
                                std::vector<double>& grad, double* log_post_out) {
    const NMixNutsModel* m =
        static_cast<const NMixNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = abun_nuts_eval(m->d, params.data(),
                                     m->sigma_beta, m->sigma_logr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

} // namespace tulpaObs

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .tobs_abun_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_abun_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                       double sigma_beta, double sigma_logr) {
    tulpaObs::NMixNutsData d = tulpaObs::abun_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::abun_nuts_eval(
        d, theta.begin(), sigma_beta, sigma_logr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the N-mixture target via tulpa's engine and the FullGradFn
// (gradient mode "H"). `theta0` is the warm-start (the Laplace mode); `inv_metric`
// an optional length-n_params inverse-mass diagonal (the Laplace curvature).
// [[Rcpp::export]]
Rcpp::List cpp_abun_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                         double sigma_beta, double sigma_logr,
                         Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                         int n_iter, int n_warmup, int max_treedepth,
                         double adapt_delta, int seed, bool verbose) {
    tulpaObs::NMixNutsModel m;
    m.d = tulpaObs::abun_nuts_build_data(spec);
    m.sigma_beta = sigma_beta; m.sigma_logr = sigma_logr;
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    tulpa::LikelihoodSpec lspec;
    lspec.name = "abun_nmix";
    lspec.n_processes = 1;
    lspec.gradient_fn = &tulpaObs::abun_nuts_full_grad;

    tulpa::ModelData data;
    data.N = m.d.n_sites;
    data.n_processes = 1;
    data.sigma_beta = sigma_beta;
    data.model_response_data = &m;
    data.likelihood_spec = &lspec;
    data.sharing.init(1);
    data.zi_type = tulpa::ZIType::NONE;
    data.p_zi = 0; data.p_oi = 0;

    tulpa::ParamLayout layout;
    layout.total_params = m.d.total;

    tulpa::set_gradient_mode_str("H");

    std::vector<double> init(theta0.begin(), theta0.end());
    std::vector<double> imv;
    const double* im = nullptr;
    if (inv_metric.isNotNull()) {
        Rcpp::NumericVector v(inv_metric);
        imv.assign(v.begin(), v.end()); im = imv.data();
    }

    tulpa::NUTSFn run_nuts = tulpa::get_nuts_fn();
    tulpa::NUTSResult result = {};
    run_nuts(&data, &layout, init.data(), m.d.total, n_iter, n_warmup,
             max_treedepth, adapt_delta, static_cast<unsigned int>(seed),
             verbose ? 1 : 0, im, &result);

    const int n_samples = result.n_sample, np = m.d.total;
    Rcpp::NumericMatrix draws(n_samples, np);
    Rcpp::NumericVector lp(n_samples), ap(n_samples);
    Rcpp::IntegerVector div(n_samples), td(n_samples);
    for (int s = 0; s < n_samples; ++s) {
        for (int j = 0; j < np; ++j) draws(s, j) = result.samples[s * np + j];
        lp[s] = result.log_prob[s]; ap[s] = result.accept_prob[s];
        div[s] = result.divergent[s]; td[s] = result.treedepth[s];
    }
    const double epsilon = result.epsilon;
    result.free_buffers();
    return Rcpp::List::create(
        Rcpp::Named("draws") = draws, Rcpp::Named("log_prob") = lp,
        Rcpp::Named("accept_prob") = ap, Rcpp::Named("divergent") = div,
        Rcpp::Named("treedepth") = td, Rcpp::Named("epsilon") = epsilon,
        Rcpp::Named("n_params") = np);
}
