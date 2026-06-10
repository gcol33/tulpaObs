// dyn_abun_nuts.cpp
// NUTS target for the Dail-Madsen open N-mixture family (dyn_abun()). The flat
// coefficient vector is theta = (beta_lambda, beta_p, beta_omega, beta_gamma) and
// the joint log-posterior is the forward marginal (dyn_abun_kernel.h) plus weak
// Gaussian priors. The shared engine (nuts_engine.h) drives tulpa's NUTS; the R
// reference .tobs_dyn_abun_nuts_logpost (R/dyn_abun_nuts.R) is the oracle that
// cpp_dyn_abun_nuts_joint_logpost is cross-checked against.

#include <Rcpp.h>
#include <vector>
#include "dyn_abun_kernel.h"
#include "nuts_engine.h"
#include "nuts_field_block.h"   // FieldBlock (shared fixed-hyper areal field)

namespace tulpaObs {

struct DynNutsModel {
    int n_sites = 0, T = 0, J = 0, K = 0;
    int p_lam = 0, p_p = 0, p_om = 0, p_gm = 0, total = 0;
    bool use_nb = false;                          // NB initial abundance
    int o_logr = -1;                             // trailing log r coord (NB only)
    double sigma_beta = 10.0;
    std::vector<int> y;                          // site-major, season, visit; -1 = NA
    Rcpp::NumericMatrix X_lambda, X_p, X_omega, X_gamma;
    // Optional single intercept RE on the initial-abundance (lambda) arm
    // (tulpaObs#51): per-site offset sigma_re * z[group], non-centered, with one
    // log_sigma_re hyperparameter. re_arm = -1 none, 0 lambda. The RE block
    // [z_1..z_G, log_sigma_re] follows the (optional) log r coord.
    int re_arm = -1, n_re_groups = 0, o_re_z = 0, o_re_logsig = 0, n_pre_re = 0;
    double sigma_re_lsd = 1.5;
    std::vector<int> re_group;        // 0-based group per site
    // Optional fixed-hyper areal field on the initial-abundance (lambda) arm
    // (tulpaObs#72): the shared non-centered field z = Linv %*% raw added to
    // eta_lambda (nuts_field_block.h). Field XOR RE (gated upstream).
    FieldBlock field;
};

inline DynNutsModel dyn_nuts_build(const Rcpp::List& spec) {
    DynNutsModel m;
    Rcpp::IntegerVector y = spec["y"];
    m.X_lambda = Rcpp::as<Rcpp::NumericMatrix>(spec["X_lambda"]);
    m.X_p      = Rcpp::as<Rcpp::NumericMatrix>(spec["X_p"]);
    m.X_omega  = Rcpp::as<Rcpp::NumericMatrix>(spec["X_omega"]);
    m.X_gamma  = Rcpp::as<Rcpp::NumericMatrix>(spec["X_gamma"]);
    m.n_sites = Rcpp::as<int>(spec["n_sites"]);
    m.T = Rcpp::as<int>(spec["T"]);
    m.J = Rcpp::as<int>(spec["J"]);
    m.K = Rcpp::as<int>(spec["K_max"]);
    if (spec.containsElementNamed("use_nb")) m.use_nb = Rcpp::as<bool>(spec["use_nb"]);
    m.p_lam = m.X_lambda.ncol(); m.p_p = m.X_p.ncol();
    m.p_om = m.X_omega.ncol(); m.p_gm = m.X_gamma.ncol();
    m.total = m.p_lam + m.p_p + m.p_om + m.p_gm;
    if (m.use_nb) { m.o_logr = m.total; m.total += 1; }   // log r after the betas
    m.n_pre_re = m.total;                                 // coords under the beta prior
    if (spec.containsElementNamed("re_arm")) m.re_arm = Rcpp::as<int>(spec["re_arm"]);
    if (m.re_arm == 0) {
        Rcpp::IntegerVector rg = spec["re_group"];
        if ((int) rg.size() != m.n_sites) Rcpp::stop("re_group must have length n_sites");
        m.re_group.resize(m.n_sites);
        for (int i = 0; i < m.n_sites; ++i) m.re_group[i] = rg[i] - 1;
        m.n_re_groups = Rcpp::as<int>(spec["n_re_groups"]);
        if (spec.containsElementNamed("sigma_re_lsd"))
            m.sigma_re_lsd = Rcpp::as<double>(spec["sigma_re_lsd"]);
        m.o_re_z = m.total; m.o_re_logsig = m.total + m.n_re_groups;
        m.total = m.o_re_logsig + 1;
    } else { m.re_arm = -1; }
    m.field = field_block_build(spec, m.total, m.n_sites);
    m.total += field_block_size(m.field);
    m.y.assign(y.begin(), y.end());
    return m;
}

inline double dyn_nuts_eval(const DynNutsModel& m, const double* theta, double* grad) {
    const int p_lam = m.p_lam, p_p = m.p_p, p_om = m.p_om, p_gm = m.p_gm;
    const int o_lam = 0, o_p = p_lam, o_om = p_lam + p_p, o_gm = p_lam + p_p + p_om;
    for (int j = 0; j < m.total; ++j) grad[j] = 0.0;
    const double eta_logr = m.use_nb ? theta[m.o_logr] : 0.0;
    const bool has_re = m.re_arm == 0;
    const double sigma_re = has_re ? std::exp(theta[m.o_re_logsig]) : 0.0;
    double grad_logr = 0.0, grad_re_logsig = 0.0;
    const bool has_field = m.field.active();
    std::vector<double> zfield, grad_z;
    field_block_forward(m.field, theta, zfield);
    field_block_init_grad(m.field, grad_z);
    double lp = 0.0;
    for (int i = 0; i < m.n_sites; ++i) {
        double el = 0.0, ep = 0.0, eo = 0.0, eg = 0.0;
        for (int k = 0; k < p_lam; ++k) el += m.X_lambda(i, k) * theta[o_lam + k];
        for (int k = 0; k < p_p;   ++k) ep += m.X_p(i, k)      * theta[o_p + k];
        for (int k = 0; k < p_om;  ++k) eo += m.X_omega(i, k)  * theta[o_om + k];
        for (int k = 0; k < p_gm;  ++k) eg += m.X_gamma(i, k)  * theta[o_gm + k];
        const double re_off = has_re ? sigma_re * theta[m.o_re_z + m.re_group[i]] : 0.0;
        if (has_re) el += re_off;
        if (has_field) el += zfield[m.field.field_map[i]];
        DynAbunSiteResult r = compute_dyn_abun_site(
            m.y.data() + (std::size_t)i * m.T * m.J, m.T, m.J, m.K, el, ep, eo, eg,
            m.use_nb, eta_logr);
        lp += r.log_lik;
        for (int k = 0; k < p_lam; ++k) grad[o_lam + k] += r.grad_eta_lambda * m.X_lambda(i, k);
        for (int k = 0; k < p_p;   ++k) grad[o_p + k]   += r.grad_eta_p      * m.X_p(i, k);
        for (int k = 0; k < p_om;  ++k) grad[o_om + k]  += r.grad_eta_omega  * m.X_omega(i, k);
        for (int k = 0; k < p_gm;  ++k) grad[o_gm + k]  += r.grad_eta_gamma  * m.X_gamma(i, k);
        grad_logr += r.grad_eta_logr;
        if (has_re) {
            grad[m.o_re_z + m.re_group[i]] += sigma_re * r.grad_eta_lambda;
            grad_re_logsig += r.grad_eta_lambda * re_off;
        }
        if (has_field) grad_z[m.field.field_map[i]] += r.grad_eta_lambda;
    }
    if (m.use_nb) grad[m.o_logr] += grad_logr;
    const double ib2 = 1.0 / (m.sigma_beta * m.sigma_beta);
    for (int k = 0; k < m.n_pre_re; ++k) {     // beta (+ log r) prior only
        lp -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k] -= ib2 * theta[k];
    }
    if (has_re) {
        for (int g = 0; g < m.n_re_groups; ++g) {
            const double zg = theta[m.o_re_z + g];
            lp -= 0.5 * zg * zg;
            grad[m.o_re_z + g] -= zg;
        }
        const double ls = theta[m.o_re_logsig], ils2 = 1.0 / (m.sigma_re_lsd * m.sigma_re_lsd);
        lp -= 0.5 * ils2 * ls * ls;
        grad[m.o_re_logsig] += grad_re_logsig - ils2 * ls;
    }
    lp += field_block_backward(m.field, theta, grad_z, grad);
    return lp;
}

inline void dyn_nuts_full_grad(const std::vector<double>& params,
                               const tulpa::ModelData& data,
                               const tulpa::ParamLayout& /*layout*/,
                               std::vector<double>& grad, double* log_post_out) {
    const DynNutsModel* m = static_cast<const DynNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->total, 0.0);
    const double lp = dyn_nuts_eval(*m, params.data(), grad.data());
    if (log_post_out) *log_post_out = lp;
}

}  // namespace tulpaObs

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_dyn_abun_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                           double sigma_beta) {
    tulpaObs::DynNutsModel m = tulpaObs::dyn_nuts_build(spec);
    m.sigma_beta = sigma_beta;
    if ((int) theta.size() != m.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), m.total);
    Rcpp::NumericVector grad(m.total);
    const double lp = tulpaObs::dyn_nuts_eval(m, theta.begin(), grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// [[Rcpp::export]]
Rcpp::List cpp_dyn_abun_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                             double sigma_beta,
                             Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                             int n_iter, int n_warmup, int max_treedepth,
                             double adapt_delta, int seed, bool verbose) {
    tulpaObs::DynNutsModel m = tulpaObs::dyn_nuts_build(spec);
    m.sigma_beta = sigma_beta;
    return tulpaObs::run_tulpa_nuts(&tulpaObs::dyn_nuts_full_grad, &m, m.total,
                                    theta0, sigma_beta, inv_metric, n_iter, n_warmup,
                                    max_treedepth, adapt_delta, seed, verbose);
}
