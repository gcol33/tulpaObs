// nuts_engine.h
// Generic driver around tulpa's NUTS engine for the in-tree FullGradFn targets
// (the count-marginal families in marginal_count_nuts.h and the distance family
// in distance_nuts.cpp). Every such target is a flat real parameter vector with
// a full-gradient closure reaching its model state through
// ModelData.model_response_data, so the engine plumbing -- LikelihoodSpec /
// ModelData / ParamLayout setup, the run_nuts call, and the draws/diagnostics
// marshalling -- is identical and lives here once. Each family supplies only its
// FullGradFn and an opaque model pointer.

#ifndef TULPAOBS_NUTS_ENGINE_H
#define TULPAOBS_NUTS_ENGINE_H

#include <Rcpp.h>
#include <vector>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>

namespace tulpaObs {

// Run tulpa NUTS on a full-gradient target. `grad_fn` is the tulpa FullGradFn
// (it reads its model through ModelData.model_response_data = `model_ptr`);
// `n_params` the parameter dimension; `inv_metric` an optional diagonal inverse
// mass matrix (length n_params). Returns draws + sampler diagnostics.
inline Rcpp::List run_tulpa_nuts(
    decltype(tulpa::LikelihoodSpec::gradient_fn) grad_fn,
    void* model_ptr, int n_params,
    const Rcpp::NumericVector& theta0, double sigma_beta,
    Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
    int n_iter, int n_warmup, int max_treedepth, double adapt_delta,
    int seed, bool verbose
) {
    if ((int) theta0.size() != n_params)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), n_params);

    tulpa::LikelihoodSpec lspec;
    lspec.name = "tulpaobs_marginal";
    lspec.n_processes = 1;
    lspec.gradient_fn = grad_fn;

    tulpa::ModelData data;
    data.N = n_params;
    data.n_processes = 1;
    data.sigma_beta = sigma_beta;
    data.model_response_data = model_ptr;
    data.likelihood_spec = &lspec;
    data.sharing.init(1);
    data.zi_type = tulpa::ZIType::NONE;
    data.p_zi = 0; data.p_oi = 0;

    tulpa::ParamLayout layout;
    layout.total_params = n_params;

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
    run_nuts(&data, &layout, init.data(), n_params, n_iter, n_warmup,
             max_treedepth, adapt_delta, static_cast<unsigned int>(seed),
             verbose ? 1 : 0, im, &result);

    const int n_samples = result.n_sample, np = n_params;
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

}  // namespace tulpaObs

#endif  // TULPAOBS_NUTS_ENGINE_H
