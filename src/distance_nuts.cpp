// distance_nuts.cpp
// NUTS target for the binned distance-sampling family (distance()). The flat
// coefficient vector is
//   theta = (beta_lambda [p_lam], beta_sigma [p_sig],
//            eta_b [1, hazard-rate only], log_r [1, NB only])
// and the joint log-posterior is the distance marginal (distance_kernel.h) plus
// weak Gaussian priors. The per-site gradient is the design-sandwiched eta
// gradient the kernel returns. The shared engine (nuts_engine.h) drives tulpa's
// NUTS; the R reference .tobs_distance_nuts_logpost (R/distance_nuts.R) is the
// oracle that cpp_distance_nuts_joint_logpost is cross-checked against.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include "distance_kernel.h"
#include "nuts_engine.h"

namespace tulpaObs {

struct DistNutsModel {
    int n_sites = 0, n_bins = 0, p_lam = 0, p_sig = 0, key = 0, K_max = 0;
    bool is_nb = false, hazard = false;
    double sigma_beta = 10.0, sigma_shape = 1.5, sigma_logr = 1.5;
    Rcpp::IntegerMatrix y;
    Rcpp::NumericMatrix X_lambda, X_sigma;
    DistQuad quad;
    int total = 0;
};

inline DistNutsModel dist_nuts_build(const Rcpp::List& spec) {
    DistNutsModel m;
    m.y        = Rcpp::as<Rcpp::IntegerMatrix>(spec["y"]);
    m.X_lambda = Rcpp::as<Rcpp::NumericMatrix>(spec["X_lambda"]);
    m.X_sigma  = Rcpp::as<Rcpp::NumericMatrix>(spec["X_sigma"]);
    m.n_sites  = m.y.nrow();
    m.n_bins   = m.y.ncol();
    m.p_lam    = m.X_lambda.ncol();
    m.p_sig    = m.X_sigma.ncol();
    m.key      = Rcpp::as<int>(spec["key"]);
    m.K_max    = Rcpp::as<int>(spec["K_max"]);
    m.is_nb    = Rcpp::as<bool>(spec["is_nb"]);
    m.hazard   = (m.key == DIST_HAZARD);
    Rcpp::NumericVector cutpoints = spec["cutpoints"];
    const int transect  = Rcpp::as<int>(spec["transect"]);
    const int quad_order = Rcpp::as<int>(spec["quad_order"]);
    std::vector<double> cut(cutpoints.begin(), cutpoints.end());
    m.quad = dist_build_quad(cut, transect, quad_order);
    m.total = m.p_lam + m.p_sig + (m.hazard ? 1 : 0) + (m.is_nb ? 1 : 0);
    return m;
}

// Joint log-posterior + gradient over the distance coefficient vector.
inline double dist_nuts_eval(const DistNutsModel& m, const double* theta,
                             double* grad) {
    const int p_lam = m.p_lam, p_sig = m.p_sig;
    int idx = p_lam + p_sig;
    const double eta_b = m.hazard ? theta[idx] : 0.0;
    const int b_idx = idx; if (m.hazard) ++idx;
    const int lr_idx = idx;
    const double r = m.is_nb ? std::exp(theta[lr_idx])
                             : std::numeric_limits<double>::infinity();
    for (int j = 0; j < m.total; ++j) grad[j] = 0.0;

    std::vector<int> y_site(m.n_bins);
    double lp = 0.0;
    for (int s = 0; s < m.n_sites; ++s) {
        double eta_lambda = 0.0, eta_sigma = 0.0;
        for (int k = 0; k < p_lam; ++k) eta_lambda += m.X_lambda(s, k) * theta[k];
        for (int k = 0; k < p_sig; ++k) eta_sigma  += m.X_sigma(s, k) * theta[p_lam + k];
        for (int b = 0; b < m.n_bins; ++b) y_site[b] = m.y(s, b);
        const DistSiteResult res = compute_distance_site(
            y_site.data(), m.n_bins, eta_lambda, eta_sigma, eta_b, m.key,
            m.quad, m.K_max, r);
        lp += res.log_lik;
        for (int k = 0; k < p_lam; ++k) grad[k] += res.grad_eta_lambda * m.X_lambda(s, k);
        for (int k = 0; k < p_sig; ++k) grad[p_lam + k] += res.grad_eta_d[0] * m.X_sigma(s, k);
        if (m.hazard) grad[b_idx] += res.grad_eta_d[1];
        if (m.is_nb)  grad[lr_idx] += res.grad_theta;
    }

    const double ib2 = 1.0 / (m.sigma_beta * m.sigma_beta);
    for (int k = 0; k < p_lam + p_sig; ++k) {
        lp -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k] -= ib2 * theta[k];
    }
    if (m.hazard) {
        const double is2 = 1.0 / (m.sigma_shape * m.sigma_shape);
        lp -= 0.5 * is2 * eta_b * eta_b;
        grad[b_idx] -= is2 * eta_b;
    }
    if (m.is_nb) {
        const double lr = theta[lr_idx], il2 = 1.0 / (m.sigma_logr * m.sigma_logr);
        lp -= 0.5 * il2 * lr * lr;
        grad[lr_idx] -= il2 * lr;
    }
    return lp;
}

inline void dist_nuts_full_grad(const std::vector<double>& params,
                                const tulpa::ModelData& data,
                                const tulpa::ParamLayout& /*layout*/,
                                std::vector<double>& grad, double* log_post_out) {
    const DistNutsModel* m =
        static_cast<const DistNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->total, 0.0);
    const double lp = dist_nuts_eval(*m, params.data(), grad.data());
    if (log_post_out) *log_post_out = lp;
}

}  // namespace tulpaObs

using namespace Rcpp;

// Full-vector joint log-posterior + gradient, the cross-check for the R oracle.
// [[Rcpp::export]]
Rcpp::List cpp_distance_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                           double sigma_beta, double sigma_shape,
                                           double sigma_logr) {
    tulpaObs::DistNutsModel m = tulpaObs::dist_nuts_build(spec);
    m.sigma_beta = sigma_beta; m.sigma_shape = sigma_shape; m.sigma_logr = sigma_logr;
    if ((int) theta.size() != m.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), m.total);
    Rcpp::NumericVector grad(m.total);
    const double lp = tulpaObs::dist_nuts_eval(m, theta.begin(), grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the distance target via tulpa's engine and the shared FullGradFn.
// [[Rcpp::export]]
Rcpp::List cpp_distance_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                             double sigma_beta, double sigma_shape, double sigma_logr,
                             Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                             int n_iter, int n_warmup, int max_treedepth,
                             double adapt_delta, int seed, bool verbose) {
    tulpaObs::DistNutsModel m = tulpaObs::dist_nuts_build(spec);
    m.sigma_beta = sigma_beta; m.sigma_shape = sigma_shape; m.sigma_logr = sigma_logr;
    return tulpaObs::run_tulpa_nuts(&tulpaObs::dist_nuts_full_grad, &m, m.total,
                                    theta0, sigma_beta, inv_metric, n_iter, n_warmup,
                                    max_treedepth, adapt_delta, seed, verbose);
}
