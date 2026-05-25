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

inline double occ_sigma(double eta) {
    if (eta > 0) return 1.0 / (1.0 + std::exp(-eta));
    double e = std::exp(eta);
    return e / (1.0 + e);
}

// LikelihoodFn<double>: per-obs log p(D_i | eta_i) of the scaled Bernoulli.
double occ_ll_double(
    int i, const double* eta, const double& /*logit_zi*/,
    const double& /*logit_oi*/, const std::vector<double>& /*params*/,
    const ModelData& /*data*/, const ParamLayout& /*layout*/,
    const void* model_data
) {
    const auto* r = static_cast<const OccupancyResponse*>(model_data);
    const double q = r->q[i];
    if (q <= 0.0) return 0.0;
    double mu = q * occ_sigma(eta[0]);
    mu = std::max(std::min(mu, 1.0 - 1e-15), 1e-15);
    return r->y[i] ? std::log(mu) : std::log(1.0 - mu);
}

// EtaWeightsFn: per-obs eta-space score + expected (Fisher) information.
//   grad     = (D - mu) (1 - sigma) / (1 - mu),     mu = q sigma
//   neg_hess = q sigma (1 - sigma)^2 / (1 - mu)      [expected information]
void occ_eta_weights(
    int i, const double* eta, double /*logit_zi*/, double /*logit_oi*/,
    const std::vector<double>& /*params*/, const ModelData& /*data*/,
    const ParamLayout& /*layout*/, const void* model_data,
    double* grad_eta, double* neg_hess_eta
) {
    const auto* r = static_cast<const OccupancyResponse*>(model_data);
    const double q = r->q[i];
    if (q <= 0.0) { grad_eta[0] = 0.0; neg_hess_eta[0] = 0.0; return; }
    const double s     = occ_sigma(eta[0]);
    const double mu    = q * s;
    const double denom = std::max(1.0 - mu, 1e-12);   // 1 - q sigma
    grad_eta[0]     = (r->y[i] - mu) * (1.0 - s) / denom;
    neg_hess_eta[0] = q * s * (1.0 - s) * (1.0 - s) / denom;
}

// Owns the spec object + response so both outlive the XPtr (parked in
// NestedLikelihood::keepalive).
struct OccupancyBundle {
    LikelihoodSpec    spec;
    OccupancyResponse resp;
};

} // namespace

// Build the occupancy likelihood for tulpa_nested_laplace(likelihood = ).
// Returns an external pointer to a tulpa::NestedLikelihood owning its spec +
// {D, q} response; the XPtr finalizer frees both at garbage collection.
// [[Rcpp::export]]
SEXP occ_make_nested_likelihood(Rcpp::NumericVector y,
                                Rcpp::NumericVector det_prob) {
    if (det_prob.size() != y.size()) {
        Rcpp::stop("occ_make_nested_likelihood: det_prob and y must have equal length.");
    }
    auto bundle = std::make_shared<OccupancyBundle>();
    bundle->resp.y.assign(y.begin(), y.end());
    bundle->resp.q.assign(det_prob.begin(), det_prob.end());
    bundle->spec.n_processes    = 1;
    bundle->spec.name           = "occupancy_scaled_bernoulli";
    bundle->spec.ll_double      = &occ_ll_double;
    bundle->spec.eta_weights_fn = &occ_eta_weights;
    bundle->spec.n_extra_params = 0;

    auto* lk = new NestedLikelihood;
    lk->spec          = &bundle->spec;
    lk->response_data = &bundle->resp;
    lk->keepalive     = bundle;   // shared_ptr<OccupancyBundle> -> shared_ptr<void>

    return Rcpp::XPtr<NestedLikelihood>(lk, true);
}
