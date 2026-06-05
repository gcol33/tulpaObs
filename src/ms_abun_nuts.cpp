// ms_abun_nuts.cpp
// C++ joint log-posterior + gradient for the community / multispecies N-mixture
// (ms_abun()) NUTS target. The R reference .tobs_ms_abun_nuts_logpost
// (R/ms_abun_nuts.R) is the oracle; this port mirrors it and is cross-checked
// against it before driving tulpa's NUTS engine.
//
// The non-spatial community N-mixture is per-species Royle (2004) with Gaussian
// community hyperpriors on the per-species coefficients:
//   N_{s,i}       ~ Poisson(lambda_{s,i})            (or NegBin(lambda, r_s))
//   y_{s,i,j}|N   ~ Binomial(N_{s,i}, p_{s,i,j})
//   log lambda_{s,i} = X_lambda_i  . (mu_lambda + b_lambda_s)
//   logit p_{s,i,j}  = X_p_{ij}    . (mu_p      + b_p_s)
//   b_lambda_s ~ N(0, Sigma_lambda),  b_p_s ~ N(0, Sigma_p)
// Under NB the dispersion is a per-species random effect log_r_s ~
// N(mu_log_r, sigma_log_r^2), with r_s = exp(mu_log_r + b_logr_s) constant across
// that species' sites. N_{s,i} integrates out per species-site in closed form;
// the Laplace-EM (nmix_laplace_re) profiles the per-species deviations and
// community covariances out, returning a Gaussian community-mean posterior. NUTS
// instead samples EVERYTHING jointly -- the community means, the per-species
// deviations, AND the community covariances -- from the exact joint posterior,
// which gives calibrated (non-Gaussian) community-mean / covariance intervals and
// the per-(species, site) pointwise likelihood WAIC / LOO need.
//
// The target factorises as
//   log p = sum_{s,i} log m_{s,i}(theta)                 # per-species-site marginal
//         - 0.5 ||mu_coef||^2 / sigma.beta^2             # community-mean priors
//         [ - 0.5 mu_log_r^2 / sigma.logr^2 ]            # (NB)
//         - 0.5 sum_s b_{s,arm}' Sigma_arm^{-1} b_{s,arm}   (per arm)  # community RE
//         - 0.5 S sum_arm log|Sigma_arm|                 # MVN normalisers
//         + log p(Sigma coords)                          # log-Cholesky hyperpriors
// where m_{s,i} is the Royle marginal exposed by compute_nmix_site()
// (nmix_kernel.h, the same kernel the Laplace fit and the AGHQ RE path use); it
// returns grad_eta_lambda, grad_eta_p, and (NB) grad_theta = d log m / d log_r, so
// the coefficient gradient is the design-sandwiched eta-gradient -- no new
// likelihood math. The arm community covariances are carried by their
// log-Cholesky factors and use the shared helpers in community_chol.h.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
#include "nmix_kernel.h"
#include "community_chol.h"

using namespace Rcpp;

namespace tulpaObs {

// Log-Cholesky coordinate hyperprior scalars for the community covariances.
struct MsNmixPri {
    double logdiag_mean = 0.0, logdiag_sd = 1.5, offdiag_sd = 1.0;
};

// Marshalled per-fit data for the community N-mixture NUTS target. The long form
// is what .tobs_ms_nmix_longform() produces: y / site_idx / species_idx in the
// stacked per-(species, site, visit) order, X_lambda site-level (n_sites rows,
// shared across species), X_p long form (n_obs rows). Per-(species, site) visit
// row lists are built once at construction so the eval loop gathers each
// species-site's counts / detection etas contiguously for compute_nmix_site().
struct MsNmixNutsData {
    int n_sites = 0, n_obs = 0, n_species = 0, p_lam = 0, p_p = 0, K_max = 0;
    bool is_nb = false;
    std::vector<int> y;                                 // length n_obs (long form)
    NumericMatrix X_lambda;                             // n_sites x p_lam
    NumericMatrix X_p;                                  // n_obs   x p_p
    std::vector<std::vector<std::vector<int>>> obs;     // [s][i] -> rows into y / X_p
    // Per-(species, site) lgamma cache. The combinatorial lgamma terms of the
    // Royle marginal are eta-independent, so they are built once here (NOT on
    // every leapfrog step, which dominated the runtime); the eval loop calls
    // compute_nmix_site_cached with the precomputed cache.
    std::vector<std::vector<NMixSiteCache>> cache;      // [s][i]

    // Packed-coordinate layout (mirrors .tobs_ms_abun_nuts_layout).
    int P_tot = 0, q_lam = 0, q_p = 0, q_logr = 0, total = 0;
    int mu_off = 0, b_off = 0, chol_lam_off = 0, chol_p_off = 0, chol_logr_off = 0;
};

inline void ms_abun_nuts_layout(MsNmixNutsData& d) {
    d.P_tot   = d.p_lam + d.p_p + (d.is_nb ? 1 : 0);
    d.q_lam   = d.p_lam * (d.p_lam + 1) / 2;
    d.q_p     = d.p_p   * (d.p_p   + 1) / 2;
    d.q_logr  = d.is_nb ? 1 : 0;
    d.mu_off  = 0;
    d.b_off   = d.P_tot;
    d.chol_lam_off  = d.P_tot + d.n_species * d.P_tot;
    d.chol_p_off    = d.chol_lam_off + d.q_lam;
    d.chol_logr_off = d.chol_p_off + d.q_p;
    d.total   = d.chol_logr_off + d.q_logr;
}

inline MsNmixNutsData ms_abun_nuts_build_data(const Rcpp::List& spec) {
    MsNmixNutsData d;
    IntegerVector y        = spec["y"];
    IntegerVector site_idx = spec["site_idx"];        // 1-based, into X_lambda rows
    IntegerVector sp_idx   = spec["species_idx"];     // 1-based
    d.X_lambda  = Rcpp::as<NumericMatrix>(spec["X_lambda"]);
    d.X_p       = Rcpp::as<NumericMatrix>(spec["X_p"]);
    d.n_sites   = Rcpp::as<int>(spec["n_sites"]);
    d.n_species = Rcpp::as<int>(spec["n_species"]);
    d.K_max     = Rcpp::as<int>(spec["K_max"]);
    d.is_nb     = Rcpp::as<bool>(spec["is_nb"]);
    d.n_obs     = y.size();
    d.p_lam     = d.X_lambda.ncol();
    d.p_p       = d.X_p.ncol();
    d.y.assign(y.begin(), y.end());
    d.obs.assign(d.n_species,
                 std::vector<std::vector<int>>(d.n_sites, std::vector<int>()));
    for (int o = 0; o < d.n_obs; ++o) {
        const int s = sp_idx[o] - 1, i = site_idx[o] - 1;
        if (s < 0 || s >= d.n_species)
            Rcpp::stop("species_idx out of range in ms_abun_nuts_build_data");
        if (i < 0 || i >= d.n_sites)
            Rcpp::stop("site_idx out of range in ms_abun_nuts_build_data");
        d.obs[s][i].push_back(o);
    }
    // Build the eta-independent lgamma cache for every observed (species, site).
    d.cache.assign(d.n_species, std::vector<NMixSiteCache>(d.n_sites));
    std::vector<int> y_site;
    for (int s = 0; s < d.n_species; ++s) {
        for (int i = 0; i < d.n_sites; ++i) {
            const std::vector<int>& obs = d.obs[s][i];
            const int J = (int) obs.size();
            if (J == 0) continue;
            y_site.resize(J);
            for (int jj = 0; jj < J; ++jj) y_site[jj] = d.y[obs[jj]];
            d.cache[s][i] = nmix_precompute_site(y_site.data(), J, d.K_max);
        }
    }
    ms_abun_nuts_layout(d);
    return d;
}

// Joint log-posterior + gradient over the packed vector
//   theta = (mu, {b_s} species-major, chol_lambda, chol_p [, chol_logr]).
// NUTS maximises, so this returns the log-posterior (no negation) and writes the
// gradient into `g` (length d.total). Mirrors .tobs_ms_abun_nuts_logpost.
inline double ms_abun_nuts_eval(const MsNmixNutsData& d, const double* th,
                                double sigma_beta, double sigma_logr,
                                const MsNmixPri& pr, double* g) {
    const int P = d.P_tot, S = d.n_species, p_lam = d.p_lam, p_p = d.p_p;
    const bool nb = d.is_nb;
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;

    const double* mu = th + d.mu_off;
    const double* b  = th + d.b_off;       // species-major, length P each
    double* g_mu = g + d.mu_off;
    double* g_b  = g + d.b_off;

    // ---- data log-lik + inner gradient (per species, per site) ----
    double lp = 0.0;
    std::vector<double> eta_p_site;
    for (int s = 0; s < S; ++s) {
        const double* b_s = b + s * P;
        const double r = nb
            ? std::exp(mu[p_lam + p_p] + b_s[p_lam + p_p])
            : std::numeric_limits<double>::infinity();
        double* g_b_lam = g_b + s * P;             // coords [0, p_lam)
        double* g_b_p   = g_b + s * P + p_lam;      // coords [p_lam, p_lam + p_p)
        for (int i = 0; i < d.n_sites; ++i) {
            const std::vector<int>& obs = d.obs[s][i];
            const int J = (int) obs.size();
            if (J == 0) continue;
            double eta_lambda = 0.0;
            for (int k = 0; k < p_lam; ++k)
                eta_lambda += d.X_lambda(i, k) * (mu[k] + b_s[k]);
            eta_p_site.resize(J);
            for (int jj = 0; jj < J; ++jj) {
                const int o = obs[jj];
                double e = 0.0;
                for (int k = 0; k < p_p; ++k)
                    e += d.X_p(o, k) * (mu[p_lam + k] + b_s[p_lam + k]);
                eta_p_site[jj] = e;
            }
            const NMixSiteResult res = compute_nmix_site_cached(
                d.cache[s][i], eta_p_site.data(), eta_lambda, r);
            lp += res.log_lik;
            // abundance arm: mu_lambda and b_lambda_s share grad_eta_lambda
            for (int k = 0; k < p_lam; ++k) {
                const double gx = res.grad_eta_lambda * d.X_lambda(i, k);
                g_mu[k]    += gx;
                g_b_lam[k] += gx;
            }
            // detection arm: mu_p and b_p_s share the per-visit grad_eta_p
            for (int jj = 0; jj < J; ++jj) {
                const int o = obs[jj];
                const double ge = res.grad_eta_p[jj];
                for (int k = 0; k < p_p; ++k) {
                    const double gx = ge * d.X_p(o, k);
                    g_mu[p_lam + k] += gx;
                    g_b_p[k]        += gx;
                }
            }
            // dispersion arm: mu_log_r and b_logr_s share grad_theta (theta = log r)
            if (nb) {
                g_mu[p_lam + p_p]        += res.grad_theta;
                g_b[s * P + p_lam + p_p] += res.grad_theta;
            }
        }
    }

    // ---- community covariance: per-arm b-quadratic + log-det normaliser + chol
    //      block gradient; accumulate the b-prior into g_b. ----
    const int P_arm[3]     = {p_lam, p_p, nb ? 1 : 0};
    const int arm_start[3] = {0, p_lam, p_lam + p_p};
    const int chol_off[3]  = {d.chol_lam_off, d.chol_p_off, d.chol_logr_off};
    const int n_arms = nb ? 3 : 2;
    for (int a = 0; a < n_arms; ++a) {
        const int Pa = P_arm[a]; if (Pa == 0) continue;
        std::vector<double> C, Cinv, Si;
        chol_unpack_cpp(th + chol_off[a], Pa, C);
        lower_tri_inv(C, Pa, Cinv);
        sinv_from_cinv(Cinv, Pa, Si);
        double logdet = 0.0;
        for (int j = 0; j < Pa; ++j) logdet += 2.0 * std::log(C[(std::size_t) j * Pa + j]);

        std::vector<double> M((std::size_t) Pa * Pa, 0.0);
        double quad_sum = 0.0;
        for (int s = 0; s < S; ++s) {
            const double* bsa = b + s * P + arm_start[a];
            for (int u = 0; u < Pa; ++u) {
                double sib = 0.0;
                for (int v = 0; v < Pa; ++v) sib += Si[(std::size_t) u * Pa + v] * bsa[v];
                g_b[s * P + arm_start[a] + u] -= sib;     // b-prior gradient
                quad_sum += bsa[u] * sib;
                for (int v = 0; v < Pa; ++v) M[(std::size_t) u * Pa + v] += bsa[u] * bsa[v];
            }
        }
        lp += -0.5 * quad_sum - 0.5 * S * logdet;

        chol_block_grad_cpp(C, Si, M, Pa, (double) S, th + chol_off[a],
                            pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                            g + chol_off[a]);
        // chol coordinate hyperprior contribution to lp
        int pos = chol_off[a];
        for (int j = 0; j < Pa; ++j) {
            const double vd = th[pos++];
            const double zd = (vd - pr.logdiag_mean) / pr.logdiag_sd;
            lp += -0.5 * zd * zd;
            for (int i = j + 1; i < Pa; ++i) {
                const double vo = th[pos++];
                const double zo = vo / pr.offdiag_sd;
                lp += -0.5 * zo * zo;
            }
        }
    }

    // ---- community-mean priors: N(0, sigma.beta^2) on the coefficient means,
    //      N(0, sigma.logr^2) on the log-dispersion community mean (NB). ----
    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < p_lam + p_p; ++k) {
        g_mu[k] -= inv_sb2 * mu[k];
        lp      += -0.5 * inv_sb2 * mu[k] * mu[k];
    }
    if (nb) {
        const double inv_sl2 = 1.0 / (sigma_logr * sigma_logr);
        const double m = mu[p_lam + p_p];
        g_mu[p_lam + p_p] -= inv_sl2 * m;
        lp                += -0.5 * inv_sl2 * m * m;
    }
    return lp;
}

// NUTS model carrying the marshalled data + prior scales; the FullGradFn reaches
// it through ModelData.model_response_data.
struct MsNmixNutsModel {
    MsNmixNutsData d;
    double sigma_beta = 10.0, sigma_logr = 1.5;
    MsNmixPri pr;
};

// FullGradFn: log-posterior + gradient over the entire parameter vector.
inline void ms_abun_nuts_full_grad(const std::vector<double>& params,
                                   const tulpa::ModelData& data,
                                   const tulpa::ParamLayout& /*layout*/,
                                   std::vector<double>& grad, double* log_post_out) {
    const MsNmixNutsModel* m =
        static_cast<const MsNmixNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = ms_abun_nuts_eval(m->d, params.data(), m->sigma_beta,
                                        m->sigma_logr, m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

inline MsNmixPri ms_abun_pri_from_list(const Rcpp::List& pri) {
    MsNmixPri pr;
    pr.logdiag_mean = Rcpp::as<double>(pri["chol_logdiag_mean"]);
    pr.logdiag_sd   = Rcpp::as<double>(pri["chol_logdiag_sd"]);
    pr.offdiag_sd   = Rcpp::as<double>(pri["chol_offdiag_sd"]);
    return pr;
}

} // namespace tulpaObs

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .tobs_ms_abun_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_ms_abun_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                          Rcpp::List pri, double sigma_beta,
                                          double sigma_logr) {
    tulpaObs::MsNmixNutsData d = tulpaObs::ms_abun_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    tulpaObs::MsNmixPri pr = tulpaObs::ms_abun_pri_from_list(pri);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::ms_abun_nuts_eval(
        d, theta.begin(), sigma_beta, sigma_logr, pr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the community N-mixture target via tulpa's engine and the
// FullGradFn (gradient mode "H"). `theta0` is the warm-start (the Laplace-EM
// mode); `inv_metric` an optional length-n_params inverse-mass diagonal (the
// Laplace curvature). Returns draws + diagnostics.
// [[Rcpp::export]]
Rcpp::List cpp_ms_abun_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                            Rcpp::List pri, double sigma_beta, double sigma_logr,
                            Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                            int n_iter, int n_warmup, int max_treedepth,
                            double adapt_delta, int seed, bool verbose) {
    tulpaObs::MsNmixNutsModel m;
    m.d = tulpaObs::ms_abun_nuts_build_data(spec);
    m.sigma_beta = sigma_beta; m.sigma_logr = sigma_logr;
    m.pr = tulpaObs::ms_abun_pri_from_list(pri);
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    tulpa::LikelihoodSpec lspec;
    lspec.name = "ms_abun_nmix";
    lspec.n_processes = 1;
    lspec.gradient_fn = &tulpaObs::ms_abun_nuts_full_grad;

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
