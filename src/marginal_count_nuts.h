// marginal_count_nuts.h
// Shared NUTS target for the non-spatial count-marginal abundance families
// (the Royle 2004 N-mixture, abun_nuts.cpp; the removal-sampling family,
// removal_nuts.cpp). The flat coefficient vector is
//   theta = (beta_lambda [p_lam], beta_p [p_p], [log_r under NB])
// and the joint log-posterior is
//   log p(theta|y) = sum_i log m_i(theta) - 0.5||beta||^2/sigma.beta^2
//                    [ - 0.5 log_r^2 / sigma.logr^2 ]   (NB only)
// where m_i is the family's per-site marginal. Each family's kernel
// (compute_nmix_site / compute_removal_site) is a free function with the SAME
// signature, so the eval loop, the FullGradFn, and the tulpa NUTS driver are
// shared verbatim and the kernel enters as a plain function pointer -- no new
// likelihood math, no per-family sampler plumbing.

#ifndef TULPAOBS_MARGINAL_COUNT_NUTS_H
#define TULPAOBS_MARGINAL_COUNT_NUTS_H

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
#include "nmix_kernel.h"        // NMixSiteResult
#include "nuts_engine.h"        // run_tulpa_nuts (shared engine plumbing)
#include "nuts_field_block.h"
#include "nuts_re_block.h"   // FieldBlock (shared fixed-hyper areal field)

namespace tulpaObs {

// Per-site marginal kernel as a plain function pointer: both families' kernels
// match this signature, so the NUTS target is kernel-agnostic.
using CountKernelFn = NMixSiteResult (*)(const int*, const double*, int,
                                         double, int, double);

// Marshalled per-fit data: long form (y / site_idx site-major), X_lambda
// site-level, X_p long form, plus per-site visit row lists built once.
struct CountNutsData {
    int n_sites = 0, n_obs = 0, p_lam = 0, p_p = 0, K_max = 0;
    bool is_nb = false;
    std::vector<int> y;
    Rcpp::NumericMatrix X_lambda;
    Rcpp::NumericMatrix X_p;
    std::vector<std::vector<int>> obs_by_site;
    int total = 0;
    // Optional single-grouping intercept random effect, arm 0 = abundance
    // (lambda), 1 = detection (p). The offset is per SITE (uniform over a
    // site's visits); see nuts_re_block.h. The flat vector grows to
    // [beta_lambda, beta_p, (log_r), z_1..z_G, log_sigma_re].
    ReBlock re;
    // Optional fixed-hyper areal field on the abundance arm: the shared
    // non-centered Gaussian field z = Linv %*% raw added to eta_lambda
    // (nuts_field_block.h). The field covariance is fixed at the nested-Laplace
    // estimate; NUTS samples only the whitened raw and the betas. Field XOR RE
    // (gated upstream).
    FieldBlock field;
};

inline CountNutsData count_nuts_build_data(const Rcpp::List& spec) {
    CountNutsData d;
    Rcpp::IntegerVector y        = spec["y"];
    Rcpp::IntegerVector site_idx = spec["site_idx"];   // 1-based
    d.X_lambda = Rcpp::as<Rcpp::NumericMatrix>(spec["X_lambda"]);
    d.X_p      = Rcpp::as<Rcpp::NumericMatrix>(spec["X_p"]);
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
            Rcpp::stop("site_idx out of range in count_nuts_build_data");
        d.obs_by_site[i].push_back(o);
    }
    int base = d.p_lam + d.p_p + (d.is_nb ? 1 : 0);
    d.re = re_block_build(spec, base, d.n_sites, /*max_arm=*/1);
    base += re_block_size(d.re);
    d.field = field_block_build(spec, base, d.n_sites);
    base += field_block_size(d.field);
    d.total = base;
    return d;
}

// Joint log-posterior + gradient over (beta_lambda, beta_p[, log_r]). Returns the
// log-posterior (NUTS maximises) and writes the gradient into `grad` (length
// d.total). The coefficient gradient is the design-sandwiched eta-gradient the
// per-site kernel returns -- identical arithmetic for every count family.
inline double count_nuts_eval(const CountNutsData& d, const double* theta,
                              double sigma_beta, double sigma_logr,
                              double* grad, CountKernelFn kern) {
    const int p_lam = d.p_lam, p_p = d.p_p;
    const double r = d.is_nb
        ? std::exp(theta[p_lam + p_p])
        : std::numeric_limits<double>::infinity();
    for (int j = 0; j < d.total; ++j) grad[j] = 0.0;

    std::vector<double> eta_p_all(d.n_obs, 0.0);
    for (int o = 0; o < d.n_obs; ++o) {
        double e = 0.0;
        for (int k = 0; k < p_p; ++k) e += d.X_p(o, k) * theta[p_lam + k];
        eta_p_all[o] = e;
    }

    // Non-centered per-site intercept offset added to the chosen arm's eta.
    const double sigma_re = re_block_sigma(d.re, theta);
    double grad_logsig = 0.0;

    // Fixed-hyper areal field: z = Linv %*% raw (per unit), added to eta_lambda.
    const bool has_field = d.field.active();
    std::vector<double> zfield, grad_z;
    field_block_forward(d.field, theta, zfield);
    field_block_init_grad(d.field, grad_z);

    double lp = 0.0;
    std::vector<int>    y_site;
    std::vector<double> eta_p_site;
    for (int i = 0; i < d.n_sites; ++i) {
        const std::vector<int>& obs = d.obs_by_site[i];
        const int J = (int) obs.size();
        const double re_off = re_block_offset(d.re, sigma_re, theta, i);
        double eta_lambda = 0.0;
        for (int k = 0; k < p_lam; ++k) eta_lambda += d.X_lambda(i, k) * theta[k];
        if (d.re.arm == 0) eta_lambda += re_off;
        if (has_field) eta_lambda += zfield[d.field.field_map[i]];
        y_site.resize(J); eta_p_site.resize(J);
        for (int jj = 0; jj < J; ++jj) {
            y_site[jj]     = d.y[obs[jj]];
            eta_p_site[jj] = eta_p_all[obs[jj]] + ((d.re.arm == 1) ? re_off : 0.0);
        }
        const NMixSiteResult res = kern(
            y_site.data(), eta_p_site.data(), J, eta_lambda, d.K_max, r);
        lp += res.log_lik;
        for (int k = 0; k < p_lam; ++k)
            grad[k] += res.grad_eta_lambda * d.X_lambda(i, k);
        if (has_field) grad_z[d.field.field_map[i]] += res.grad_eta_lambda;
        double g_off = 0.0;
        for (int jj = 0; jj < J; ++jj) {
            const double ge = res.grad_eta_p[jj];
            const int o = obs[jj];
            for (int k = 0; k < p_p; ++k)
                grad[p_lam + k] += ge * d.X_p(o, k);
            if (d.re.arm == 1) g_off += ge;
        }
        if (d.re.arm == 0) g_off = res.grad_eta_lambda;
        re_block_accumulate(d.re, sigma_re, g_off, i, theta, grad, grad_logsig);
        if (d.is_nb) grad[p_lam + p_p] += res.grad_theta;
    }

    const double ib2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < p_lam + p_p; ++k) {
        lp      -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k] -= ib2 * theta[k];
    }
    if (d.is_nb) {
        const double lr = theta[p_lam + p_p];
        const double ilr2 = 1.0 / (sigma_logr * sigma_logr);
        lp -= 0.5 * ilr2 * lr * lr;
        grad[p_lam + p_p] -= ilr2 * lr;
    }
    // Non-centered RE prior: z ~ N(0, I) + Gaussian hyperprior on log_sigma_re.
    lp += re_block_backward(d.re, theta, grad_logsig, grad);
    // Whitened field prior raw ~ N(0, I); the chain grad_raw = Linv^T grad_z.
    lp += field_block_backward(d.field, theta, grad_z, grad);
    return lp;
}

// NUTS model carrying the data, prior scales, and the family kernel; the
// FullGradFn reaches it through ModelData.model_response_data.
struct CountNutsModel {
    CountNutsData d;
    double sigma_beta = 10.0, sigma_logr = 1.5;
    CountKernelFn kern = nullptr;
};

inline void count_nuts_full_grad(const std::vector<double>& params,
                                 const tulpa::ModelData& data,
                                 const tulpa::ParamLayout& /*layout*/,
                                 std::vector<double>& grad, double* log_post_out) {
    const CountNutsModel* m =
        static_cast<const CountNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = count_nuts_eval(m->d, params.data(),
                                      m->sigma_beta, m->sigma_logr,
                                      grad.data(), m->kern);
    if (log_post_out) *log_post_out = lp;
}

// Run NUTS on a count-marginal target via tulpa's engine and the shared
// FullGradFn. `kern` selects the family. Returns the draws + sampler diagnostics.
inline Rcpp::List count_nuts_run(const Rcpp::List& spec,
                                 const Rcpp::NumericVector& theta0,
                                 double sigma_beta, double sigma_logr,
                                 const std::vector<double>& inv_metric,
                                 int n_iter, int n_warmup, int max_treedepth,
                                 double adapt_delta, int seed, bool verbose,
                                 CountKernelFn kern) {
    CountNutsModel m;
    m.d = count_nuts_build_data(spec);
    m.sigma_beta = sigma_beta; m.sigma_logr = sigma_logr; m.kern = kern;
    return run_tulpa_nuts(&count_nuts_full_grad, &m, m.d.total, theta0, sigma_beta,
                          inv_metric, n_iter, n_warmup, max_treedepth, adapt_delta,
                          seed, verbose);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_MARGINAL_COUNT_NUTS_H
