// fp_occu_nuts.cpp
// NUTS target for the multistate false-positive occupancy family (fp_occu()).
// The flat coefficient vector is theta = (beta_psi, beta_p11, beta_p10, beta_b)
// and the joint log-posterior is the false-positive occupancy marginal
// (fp_occu_kernel.h) plus weak Gaussian priors. The shared engine (nuts_engine.h)
// drives tulpa's NUTS; the R reference .tobs_fp_occu_nuts_logpost
// (R/fp_occu_nuts.R) is the oracle that cpp_fp_occu_nuts_joint_logpost is
// cross-checked against.

#include <Rcpp.h>
#include <vector>
#include "tobs_shape.h"
#include "fp_occu_kernel.h"
#include "nuts_engine.h"
#include "nuts_field_block.h"   // FieldBlock (shared fixed-hyper areal field)
#include "nuts_re_block.h"      // ReBlock (shared non-centered grouped RE)

namespace tulpaObs {

struct FpNutsModel {
    int n_sites = 0, n_obs = 0;
    int p_psi = 0, p_p11 = 0, p_p10 = 0, p_b = 0, total = 0;
    double sigma_beta = 10.0;
    std::vector<int> y;
    Rcpp::NumericMatrix X_psi, X_p11, X_p10, X_b;
    std::vector<std::vector<int>> obs_by_site;
    // Optional single intercept RE on the occupancy (psi) arm, loaded on
    // eta_psi (nuts_re_block.h). The block [z_1..z_G, log_sigma_re] follows
    // the four coefficient blocks.
    int n_pre_re = 0;
    ReBlock re;
    // Optional fixed-hyper areal field on the occupancy (psi) arm: the shared
    // non-centered field z = Linv %*% raw added to eta_psi (nuts_field_block.h).
    // Field XOR RE (gated upstream).
    FieldBlock field;
};

inline FpNutsModel fp_nuts_build(const Rcpp::List& spec) {
    FpNutsModel m;
    Rcpp::IntegerVector y        = spec["y"];
    Rcpp::IntegerVector site_idx = spec["site_idx"];
    m.X_psi = Rcpp::as<Rcpp::NumericMatrix>(spec["X_psi"]);
    m.X_p11 = Rcpp::as<Rcpp::NumericMatrix>(spec["X_p11"]);
    m.X_p10 = Rcpp::as<Rcpp::NumericMatrix>(spec["X_p10"]);
    m.X_b   = Rcpp::as<Rcpp::NumericMatrix>(spec["X_b"]);
    m.n_sites = Rcpp::as<int>(spec["n_sites"]);
    m.n_obs = y.size();
    m.p_psi = m.X_psi.ncol(); m.p_p11 = m.X_p11.ncol();
    m.p_p10 = m.X_p10.ncol(); m.p_b = m.X_b.ncol();
    m.total = m.p_psi + m.p_p11 + m.p_p10 + m.p_b;
    m.n_pre_re = m.total;                                 // coords under the beta prior
    m.re = re_block_build(spec, m.total, m.n_sites);       // psi arm only
    m.total += re_block_size(m.re);
    m.field = field_block_build(spec, m.total, m.n_sites);
    m.total += field_block_size(m.field);
    m.y.assign(y.begin(), y.end());
    m.obs_by_site.assign(m.n_sites, std::vector<int>());
    for (int o = 0; o < m.n_obs; ++o) {
        const int s = site_idx[o] - 1;
        if (s < 0 || s >= m.n_sites) Rcpp::stop("site_idx out of range in fp_nuts_build");
        m.obs_by_site[s].push_back(o);
    }
    return m;
}

inline double fp_nuts_eval(const FpNutsModel& m, const double* theta, double* grad) {
    const int p_psi = m.p_psi, p_p11 = m.p_p11, p_p10 = m.p_p10, p_b = m.p_b;
    const int o_psi = 0, o_p11 = p_psi, o_p10 = p_psi + p_p11, o_b = p_psi + p_p11 + p_p10;
    for (int j = 0; j < m.total; ++j) grad[j] = 0.0;
    const double sigma_re = re_block_sigma(m.re, theta);
    double grad_re_logsig = 0.0;

    const bool has_field = m.field.active();
    std::vector<double> zfield, grad_z;
    field_block_forward(m.field, theta, zfield);
    field_block_init_grad(m.field, grad_z);

    std::vector<int> y_site;
    double lp = 0.0;
    for (int s = 0; s < m.n_sites; ++s) {
        double eta_psi = 0.0, eta_p11 = 0.0, eta_p10 = 0.0, eta_b = 0.0;
        for (int k = 0; k < p_psi; ++k) eta_psi += m.X_psi(s, k) * theta[o_psi + k];
        for (int k = 0; k < p_p11; ++k) eta_p11 += m.X_p11(s, k) * theta[o_p11 + k];
        for (int k = 0; k < p_p10; ++k) eta_p10 += m.X_p10(s, k) * theta[o_p10 + k];
        for (int k = 0; k < p_b;   ++k) eta_b   += m.X_b(s, k)   * theta[o_b + k];
        eta_psi += re_block_offset(m.re, sigma_re, theta, s);
        if (has_field) eta_psi += zfield[m.field.field_map[s]];
        const std::vector<int>& idx = m.obs_by_site[s];
        const int J = (int)idx.size();
        y_site.resize(J);
        for (int j = 0; j < J; ++j) y_site[j] = m.y[idx[j]];
        FpOccuSiteResult r = compute_fp_occu_site(y_site.data(), J,
                                                  eta_psi, eta_p11, eta_p10, eta_b);
        lp += r.log_lik;
        for (int k = 0; k < p_psi; ++k) grad[o_psi + k] += r.grad_eta_psi * m.X_psi(s, k);
        for (int k = 0; k < p_p11; ++k) grad[o_p11 + k] += r.grad_eta_p11 * m.X_p11(s, k);
        for (int k = 0; k < p_p10; ++k) grad[o_p10 + k] += r.grad_eta_p10 * m.X_p10(s, k);
        for (int k = 0; k < p_b;   ++k) grad[o_b + k]   += r.grad_eta_b   * m.X_b(s, k);
        re_block_accumulate(m.re, sigma_re, r.grad_eta_psi, s, theta, grad,
                            grad_re_logsig);
        if (has_field) grad_z[m.field.field_map[s]] += r.grad_eta_psi;
    }
    const double ib2 = 1.0 / (m.sigma_beta * m.sigma_beta);
    for (int k = 0; k < m.n_pre_re; ++k) {     // beta prior only
        lp -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k] -= ib2 * theta[k];
    }
    lp += re_block_backward(m.re, theta, grad_re_logsig, grad);
    lp += field_block_backward(m.field, theta, grad_z, grad);
    return lp;
}

inline void fp_nuts_full_grad(const std::vector<double>& params,
                              const tulpa::ModelData& data,
                              const tulpa::ParamLayout& /*layout*/,
                              std::vector<double>& grad, double* log_post_out) {
    const FpNutsModel* m = static_cast<const FpNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->total, 0.0);
    const double lp = fp_nuts_eval(*m, params.data(), grad.data());
    if (log_post_out) *log_post_out = lp;
}

}  // namespace tulpaObs

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_fp_occu_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                          double sigma_beta) {
    tulpaObs::FpNutsModel m = tulpaObs::fp_nuts_build(spec);
    m.sigma_beta = sigma_beta;
    if ((int) theta.size() != m.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), m.total);
    Rcpp::NumericVector grad(m.total);
    const double lp = tulpaObs::fp_nuts_eval(m, theta.begin(), grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// [[Rcpp::export]]
Rcpp::List cpp_fp_occu_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                            double sigma_beta,
                            Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                            int n_iter, int n_warmup, int max_treedepth,
                            double adapt_delta, int seed, bool verbose) {
    tulpaObs::FpNutsModel m = tulpaObs::fp_nuts_build(spec);
    m.sigma_beta = sigma_beta;
    return tulpaObs::run_tulpa_nuts(&tulpaObs::fp_nuts_full_grad, &m, m.total,
                                    theta0, sigma_beta, tulpaObs::shape::optional_numeric(inv_metric.get(), "inv_metric"), n_iter, n_warmup,
                                    max_treedepth, adapt_delta, seed, verbose);
}
