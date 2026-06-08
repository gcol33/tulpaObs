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
#include "nmix_kernel.h"   // NMixSiteResult
#include "nuts_engine.h"   // run_tulpa_nuts (shared engine plumbing)

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
    // Optional single-grouping intercept random effect (tulpaObs#51). re_arm =
    // -1 none, 0 abundance (lambda), 1 detection (p). The RE is a per-site
    // offset (uniform over a site's visits) sigma_re * z[group], with z ~ N(0,I)
    // non-centered and one variance hyperparameter log_sigma_re. The flat vector
    // grows to [beta_lambda, beta_p, (log_r), z_1..z_G, log_sigma_re].
    int re_arm = -1, n_re_groups = 0;
    std::vector<int> re_group;        // 0-based group per site (length n_sites)
    double sigma_re_lsd = 1.5;        // prior SD on log_sigma_re
    int o_z = 0, o_logsig = 0;        // offsets of the z block / log_sigma_re
    // Optional fixed-hyper areal field on the abundance arm (tulpaObs#51): the
    // non-centered Gaussian field z = Linv %*% raw, raw ~ N(0, I), with Linv the
    // inverse Cholesky of the FIXED field precision tau Q(rho) (a small ridge
    // proper-ises an intrinsic ICAR). The field covariance is fixed at the
    // nested-Laplace estimate; NUTS samples only the whitened raw and the betas.
    // eta_lambda[site] += z[field_map[site]]. Field XOR RE (gated upstream).
    int n_field_units = 0, o_raw = 0;
    std::vector<int> field_map;       // 0-based unit per site (length n_sites)
    std::vector<double> Linv;         // row-major n_field_units x n_field_units
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
    if (spec.containsElementNamed("re_arm")) {
        d.re_arm = Rcpp::as<int>(spec["re_arm"]);
        if (d.re_arm >= 0) {
            Rcpp::IntegerVector rg = spec["re_group"];     // 1-based, length n_sites
            if ((int) rg.size() != d.n_sites)
                Rcpp::stop("re_group must have length n_sites");
            d.re_group.resize(d.n_sites);
            for (int i = 0; i < d.n_sites; ++i) d.re_group[i] = rg[i] - 1;
            d.n_re_groups = Rcpp::as<int>(spec["n_re_groups"]);
            if (spec.containsElementNamed("sigma_re_lsd"))
                d.sigma_re_lsd = Rcpp::as<double>(spec["sigma_re_lsd"]);
            d.o_z      = base;
            d.o_logsig = base + d.n_re_groups;
            base       = d.o_logsig + 1;
        }
    }
    if (spec.containsElementNamed("n_field_units")) {
        d.n_field_units = Rcpp::as<int>(spec["n_field_units"]);
        if (d.n_field_units > 0) {
            Rcpp::IntegerVector fm = spec["field_map"];   // 1-based site -> unit
            if ((int) fm.size() != d.n_sites)
                Rcpp::stop("field_map must have length n_sites");
            d.field_map.resize(d.n_sites);
            for (int i = 0; i < d.n_sites; ++i) d.field_map[i] = fm[i] - 1;
            Rcpp::NumericMatrix Li = spec["field_Linv"];
            if (Li.nrow() != d.n_field_units || Li.ncol() != d.n_field_units)
                Rcpp::stop("field_Linv must be n_field_units x n_field_units");
            d.Linv.resize((std::size_t) d.n_field_units * d.n_field_units);
            for (int u = 0; u < d.n_field_units; ++u)
                for (int v = 0; v < d.n_field_units; ++v)
                    d.Linv[(std::size_t) u * d.n_field_units + v] = Li(u, v);
            d.o_raw = base; base += d.n_field_units;
        }
    }
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

    // Random-effect setup: non-centered per-site intercept offset
    // sigma_re * z[group(site)] added to the chosen arm's eta.
    const bool has_re = d.re_arm >= 0;
    const double sigma_re = has_re ? std::exp(theta[d.o_logsig]) : 0.0;
    double grad_logsig = 0.0;

    // Fixed-hyper areal field: z = Linv %*% raw (per unit), added to eta_lambda.
    const bool has_field = d.n_field_units > 0;
    std::vector<double> zfield, grad_z;
    if (has_field) {
        zfield.assign(d.n_field_units, 0.0);
        grad_z.assign(d.n_field_units, 0.0);
        for (int u = 0; u < d.n_field_units; ++u) {
            double zz = 0.0;
            const double* Lu = &d.Linv[(std::size_t) u * d.n_field_units];
            for (int v = 0; v < d.n_field_units; ++v) zz += Lu[v] * theta[d.o_raw + v];
            zfield[u] = zz;
        }
    }

    double lp = 0.0;
    std::vector<int>    y_site;
    std::vector<double> eta_p_site;
    for (int i = 0; i < d.n_sites; ++i) {
        const std::vector<int>& obs = d.obs_by_site[i];
        const int J = (int) obs.size();
        const double re_off = has_re ? sigma_re * theta[d.o_z + d.re_group[i]] : 0.0;
        double eta_lambda = 0.0;
        for (int k = 0; k < p_lam; ++k) eta_lambda += d.X_lambda(i, k) * theta[k];
        if (has_re && d.re_arm == 0) eta_lambda += re_off;
        if (has_field) eta_lambda += zfield[d.field_map[i]];
        y_site.resize(J); eta_p_site.resize(J);
        for (int jj = 0; jj < J; ++jj) {
            y_site[jj]     = d.y[obs[jj]];
            eta_p_site[jj] = eta_p_all[obs[jj]] + ((has_re && d.re_arm == 1) ? re_off : 0.0);
        }
        const NMixSiteResult res = kern(
            y_site.data(), eta_p_site.data(), J, eta_lambda, d.K_max, r);
        lp += res.log_lik;
        for (int k = 0; k < p_lam; ++k)
            grad[k] += res.grad_eta_lambda * d.X_lambda(i, k);
        if (has_field) grad_z[d.field_map[i]] += res.grad_eta_lambda;
        double g_off = 0.0;
        for (int jj = 0; jj < J; ++jj) {
            const double ge = res.grad_eta_p[jj];
            const int o = obs[jj];
            for (int k = 0; k < p_p; ++k)
                grad[p_lam + k] += ge * d.X_p(o, k);
            if (has_re && d.re_arm == 1) g_off += ge;
        }
        if (has_re && d.re_arm == 0) g_off = res.grad_eta_lambda;
        if (has_re) {
            // d eta / d z_g = sigma_re ; d eta / d log_sigma = sigma_re*z = re_off.
            grad[d.o_z + d.re_group[i]] += sigma_re * g_off;
            grad_logsig += g_off * re_off;
        }
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
    if (has_re) {
        // Non-centered RE prior: z ~ N(0, I); weak Gaussian hyperprior on
        // log_sigma_re.
        for (int g = 0; g < d.n_re_groups; ++g) {
            const double zg = theta[d.o_z + g];
            lp -= 0.5 * zg * zg;
            grad[d.o_z + g] -= zg;
        }
        const double ls   = theta[d.o_logsig];
        const double ils2 = 1.0 / (d.sigma_re_lsd * d.sigma_re_lsd);
        lp -= 0.5 * ils2 * ls * ls;
        grad[d.o_logsig] += grad_logsig - ils2 * ls;
    }
    if (has_field) {
        // Whitened field prior raw ~ N(0, I); the chain grad_raw = Linv^T grad_z.
        for (int v = 0; v < d.n_field_units; ++v) {
            double gr = 0.0;
            for (int u = 0; u < d.n_field_units; ++u)
                gr += d.Linv[(std::size_t) u * d.n_field_units + v] * grad_z[u];
            const double rv = theta[d.o_raw + v];
            lp -= 0.5 * rv * rv;
            grad[d.o_raw + v] += gr - rv;
        }
    }
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
                                 Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
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
