// occ_nested_likelihood.cpp
// Marginalized single-season occupancy LikelihoodSpec for tulpa's nested-Laplace
// spec path (tulpa_nested_laplace(likelihood = )).
//
// tulpa's family enum does not carry an occupancy hook: the occupancy state
// likelihood is owned here. Each site i contributes a Bernoulli on the
// detection indicator D_i = 1{>= 1 detection} with mean mu_i = q_i * sigma(eta_i),
// where sigma(eta_i) is the occupancy probability psi_i and q_i in [0, 1] is the
// per-site probability of detecting at least once given occupancy (read off the
// converged detection estimate). The latent occupancy state is integrated out
// analytically, so the converged Hessian carries the expected (Fisher)
// information q*sigma*(1-sigma)^2/(1-q*sigma) -- the true marginal occupancy
// curvature -- and fitted_eta_var is calibrated with no rescaling. q_i = 0 (an
// unvisited / held-out site) contributes zero score and zero information: it
// drops from the likelihood while keeping its latent value, interpolated by the
// field (the INLA NA-response mechanism). q_i = 1 reduces to a logit Bernoulli.
//
// .tobs_occu_state_marginal_fit() builds this from (D_i, q_i) and passes the
// returned external pointer to tulpa_nested_laplace(likelihood = ).

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include <tulpa/likelihood.h>
#include <tulpa/nested_likelihood.h>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include "tobs_math.h"   // tulpa::math::inv_logit

namespace {

using tulpa::LikelihoodSpec;
using tulpa::ModelData;
using tulpa::NestedLikelihood;
using tulpa::ParamLayout;

// Per-observation response the occupancy spec reads through model_response_data.
struct OccupancyResponse {
    std::vector<double> y;  // [N] detection indicator in {0, 1}
    std::vector<double> q;  // [N] per-site P(>=1 detection | occupied), in [0, 1]
};

// mu = q * sigma(eta) is held inside (eps, 1 - eps) so log(mu) and log(1 - mu)
// stay finite; the (1 - mu) denominator of the score and the information is
// floored separately, since it divides rather than being logged.
constexpr double OCC_MU_EPS      = 1e-15;
constexpr double OCC_DENOM_FLOOR = 1e-12;

// One site of the scaled Bernoulli: log p(D | eta), the eta-space score and the
// expected (Fisher) information, from one inverse link.
//   ll       = D log(mu) + (1 - D) log(1 - mu),       mu = q sigma(eta)
//   grad     = (D - mu) (1 - sigma) / (1 - mu)
//   neg_hess = q sigma (1 - sigma)^2 / (1 - mu)
// q = 0 (no valid visit) contributes nothing on all three.
struct OccSiteTerms {
    double ll;
    double grad;
    double neg_hess;
};

inline OccSiteTerms occ_site_terms(double y, double q, double eta) {
    if (q <= 0.0) return {0.0, 0.0, 0.0};
    const double s     = tulpa::math::inv_logit(eta);
    const double mu    = q * s;
    const double mu_c  = std::max(std::min(mu, 1.0 - OCC_MU_EPS), OCC_MU_EPS);
    const double denom = std::max(1.0 - mu, OCC_DENOM_FLOOR);
    return {
        y ? std::log(mu_c) : std::log(1.0 - mu_c),
        (y - mu) * (1.0 - s) / denom,
        q * s * (1.0 - s) * (1.0 - s) / denom
    };
}

// LikelihoodFn<double>: per-obs log p(D_i | eta_i) of the scaled Bernoulli.
double occ_ll_double(
    int i, const double* eta, const double& /*logit_zi*/,
    const double& /*logit_oi*/, const std::vector<double>& /*params*/,
    const ModelData& /*data*/, const ParamLayout& /*layout*/,
    const void* model_data
) {
    const auto* r = static_cast<const OccupancyResponse*>(model_data);
    return occ_site_terms(r->y[i], r->q[i], eta[0]).ll;
}

// EtaWeightsFn: per-obs eta-space score + expected (Fisher) information.
void occ_eta_weights(
    int i, const double* eta, double /*logit_zi*/, double /*logit_oi*/,
    const std::vector<double>& /*params*/, const ModelData& /*data*/,
    const ParamLayout& /*layout*/, const void* model_data,
    double* grad_eta, double* neg_hess_eta
) {
    const auto* r = static_cast<const OccupancyResponse*>(model_data);
    const OccSiteTerms t = occ_site_terms(r->y[i], r->q[i], eta[0]);
    grad_eta[0]     = t.grad;
    neg_hess_eta[0] = t.neg_hess;
}

// LlEtaWeightsFn: the value and the weights of one evaluation, so the Laplace
// loop's objective evaluation also hands its scatter the working weights.
double occ_ll_eta_weights(
    int i, const double* eta, double /*logit_zi*/, double /*logit_oi*/,
    const std::vector<double>& /*params*/, const ModelData& /*data*/,
    const ParamLayout& /*layout*/, const void* model_data,
    double* grad_eta, double* neg_hess_eta
) {
    const auto* r = static_cast<const OccupancyResponse*>(model_data);
    const OccSiteTerms t = occ_site_terms(r->y[i], r->q[i], eta[0]);
    grad_eta[0]     = t.grad;
    neg_hess_eta[0] = t.neg_hess;
    return t.ll;
}

// Owns the spec object + response so both outlive the XPtr (parked in
// NestedLikelihood::keepalive).
struct OccupancyBundle {
    LikelihoodSpec    spec;
    OccupancyResponse resp;
};

// The spec and its {D, q} response, as every door hands them to the engine.
std::shared_ptr<OccupancyBundle> occ_make_bundle(const Rcpp::NumericVector& y,
                                                 const Rcpp::NumericVector& det_prob,
                                                 const char* caller) {
    if (det_prob.size() != y.size()) {
        Rcpp::stop("%s: det_prob and y must have equal length.", caller);
    }
    auto bundle = std::make_shared<OccupancyBundle>();
    bundle->resp.y.assign(y.begin(), y.end());
    bundle->resp.q.assign(det_prob.begin(), det_prob.end());
    bundle->spec.n_processes       = 1;
    bundle->spec.name              = "occupancy_scaled_bernoulli";
    bundle->spec.ll_double         = &occ_ll_double;
    bundle->spec.eta_weights_fn    = &occ_eta_weights;
    bundle->spec.ll_eta_weights_fn = &occ_ll_eta_weights;
    bundle->spec.n_extra_params    = 0;
    return bundle;
}

} // namespace

// Build the occupancy likelihood for tulpa_nested_laplace(likelihood = ).
// Returns an external pointer to a tulpa::NestedLikelihood owning its spec +
// {D, q} response; the XPtr finalizer frees both at garbage collection.
// [[Rcpp::export]]
SEXP occ_make_nested_likelihood(Rcpp::NumericVector y,
                                Rcpp::NumericVector det_prob) {
    auto bundle = occ_make_bundle(y, det_prob, "occ_make_nested_likelihood");

    auto* lk = new NestedLikelihood;
    lk->spec          = &bundle->spec;
    lk->response_data = &bundle->resp;
    lk->keepalive     = bundle;   // shared_ptr<OccupancyBundle> -> shared_ptr<void>

    return Rcpp::XPtr<NestedLikelihood>(lk, true);
}

// Evaluate the occupancy spec's callbacks at `eta`, one row per site: the
// split pair (ll_double, eta_weights_fn) in columns 1-3 and the fused
// ll_eta_weights_fn in columns 4-6, each as (ll, grad, neg_hess). Reads the
// spec the engine is handed, so it exercises the registered function pointers.
// [[Rcpp::export]]
Rcpp::NumericMatrix occ_nested_likelihood_eval(Rcpp::NumericVector y,
                                               Rcpp::NumericVector det_prob,
                                               Rcpp::NumericVector eta) {
    auto bundle = occ_make_bundle(y, det_prob, "occ_nested_likelihood_eval");
    const int n = y.size();
    if (eta.size() != n) {
        Rcpp::stop("occ_nested_likelihood_eval: eta and y must have equal length.");
    }
    const LikelihoodSpec& spec = bundle->spec;
    const void* resp = &bundle->resp;
    const std::vector<double> params;
    ModelData data;
    ParamLayout layout;
    Rcpp::NumericMatrix out(n, 6);
    Rcpp::colnames(out) = Rcpp::CharacterVector::create(
        "ll_split", "grad_split", "neg_hess_split",
        "ll_fused", "grad_fused", "neg_hess_fused");
    for (int i = 0; i < n; i++) {
        const double e = eta[i];
        double g = 0.0, h = 0.0;
        out(i, 0) = spec.ll_double(i, &e, 0.0, 0.0, params, data, layout, resp);
        spec.eta_weights_fn(i, &e, 0.0, 0.0, params, data, layout, resp, &g, &h);
        out(i, 1) = g;
        out(i, 2) = h;
        g = 0.0; h = 0.0;
        out(i, 3) = spec.ll_eta_weights_fn(i, &e, 0.0, 0.0, params, data, layout,
                                           resp, &g, &h);
        out(i, 4) = g;
        out(i, 5) = h;
    }
    return out;
}
