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

namespace tulpaObs {

struct DynNutsModel {
    int n_sites = 0, T = 0, J = 0, K = 0;
    int p_lam = 0, p_p = 0, p_om = 0, p_gm = 0, total = 0;
    double sigma_beta = 10.0;
    std::vector<int> y;                          // site-major, season, visit; -1 = NA
    Rcpp::NumericMatrix X_lambda, X_p, X_omega, X_gamma;
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
    m.p_lam = m.X_lambda.ncol(); m.p_p = m.X_p.ncol();
    m.p_om = m.X_omega.ncol(); m.p_gm = m.X_gamma.ncol();
    m.total = m.p_lam + m.p_p + m.p_om + m.p_gm;
    m.y.assign(y.begin(), y.end());
    return m;
}

inline double dyn_nuts_eval(const DynNutsModel& m, const double* theta, double* grad) {
    const int p_lam = m.p_lam, p_p = m.p_p, p_om = m.p_om, p_gm = m.p_gm;
    const int o_lam = 0, o_p = p_lam, o_om = p_lam + p_p, o_gm = p_lam + p_p + p_om;
    for (int j = 0; j < m.total; ++j) grad[j] = 0.0;
    double lp = 0.0;
    for (int i = 0; i < m.n_sites; ++i) {
        double el = 0.0, ep = 0.0, eo = 0.0, eg = 0.0;
        for (int k = 0; k < p_lam; ++k) el += m.X_lambda(i, k) * theta[o_lam + k];
        for (int k = 0; k < p_p;   ++k) ep += m.X_p(i, k)      * theta[o_p + k];
        for (int k = 0; k < p_om;  ++k) eo += m.X_omega(i, k)  * theta[o_om + k];
        for (int k = 0; k < p_gm;  ++k) eg += m.X_gamma(i, k)  * theta[o_gm + k];
        DynAbunSiteResult r = compute_dyn_abun_site(
            m.y.data() + (std::size_t)i * m.T * m.J, m.T, m.J, m.K, el, ep, eo, eg);
        lp += r.log_lik;
        for (int k = 0; k < p_lam; ++k) grad[o_lam + k] += r.grad_eta_lambda * m.X_lambda(i, k);
        for (int k = 0; k < p_p;   ++k) grad[o_p + k]   += r.grad_eta_p      * m.X_p(i, k);
        for (int k = 0; k < p_om;  ++k) grad[o_om + k]  += r.grad_eta_omega  * m.X_omega(i, k);
        for (int k = 0; k < p_gm;  ++k) grad[o_gm + k]  += r.grad_eta_gamma  * m.X_gamma(i, k);
    }
    const double ib2 = 1.0 / (m.sigma_beta * m.sigma_beta);
    for (int k = 0; k < m.total; ++k) {
        lp -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k] -= ib2 * theta[k];
    }
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
