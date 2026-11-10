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
#include <algorithm>
#include <vector>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>

namespace tulpaObs {

// The double-precision log-posterior of a full-gradient target, in the
// per-observation LikelihoodFn<double> form tulpa's generic evaluator sums.
// The layout below carries no process coefficients and no latent blocks, so the
// evaluator adds no prior of its own and passes a zero linear predictor; the
// target's whole log-posterior therefore enters once, at observation 0. The
// value is the one the target's own FullGradFn writes beside its gradient, so
// the engine's runtime gradient check (and a numerical-gradient fallback)
// differentiates exactly the density the hand-coded gradient belongs to.
inline double full_grad_log_post(
    int i, const double* /*eta*/, const double& /*logit_zi*/,
    const double& /*logit_oi*/, const std::vector<double>& params,
    const tulpa::ModelData& data, const tulpa::ParamLayout& layout,
    const void* /*model_data*/
) {
    if (i != 0) return 0.0;
    const auto* spec =
        static_cast<const tulpa::LikelihoodSpec*>(data.likelihood_spec);
    std::vector<double> grad(params.size(), 0.0);
    double lp = 0.0;
    spec->gradient_fn(params, data, layout, grad, &lp);
    return lp;
}

// The chain team of one multi-chain run is min(n_chains, tulpa's OpenMP cap),
// and tulpa's cap lives in tulpa's own OpenMP runtime, so it is set through
// tulpa (OmpCapScope) for the run and restored afterwards.
using ChainTeamScope = tulpa::OmpCapScope;

// Run tulpa NUTS on a full-gradient target. `grad_fn` is the tulpa FullGradFn
// (it reads its model through ModelData.model_response_data = `model_ptr`);
// `n_params` the parameter dimension; `inv_metric` an optional diagonal inverse
// mass matrix, either empty (the engine adapts its own) or of length n_params.
// The engine reads that pointer for n_params entries without checking how many
// it was given, so a short one is rejected here.
//
// The chains run through tulpa's multi-chain runner (get_nuts_chains_fn): every
// chain starts from `theta0` with the same metric, the engine disperses the
// starts and varies the stream by chain index from the one `seed`, and the
// chains share an OpenMP team of at most `n_threads` (each chain is
// deterministic in (seed, chain), so the draws do not depend on the width). A
// target is therefore evaluated concurrently from several chains and must not
// write to its model. Returns one list of draws + diagnostics per chain.
//
// The whole parameter vector is the spec's extra-parameter block behind one
// coefficient-free process, which is the layout tulpa's generic double
// evaluator needs to reach full_grad_log_post; the sampler reads none of those
// block positions.
//
// `spec_name` labels the LikelihoodSpec and `data_N` sets ModelData.N (the
// observation count, which the community targets carry as their site count).
// Neither changes the sampled density: nothing consults the spec name, the
// double evaluator calls full_grad_log_post once per observation and takes the
// target's value at observation 0 only, and the engine otherwise reads
// ModelData.N for a checkpoint fingerprint (gated on a checkpoint path these
// targets never set) and for a latent-factor stride that is multiplied by
// latent_n_factors = 0 here. They are threaded so each caller keeps the values
// it declared rather than relying on that.
inline Rcpp::List run_tulpa_nuts(
    decltype(tulpa::LikelihoodSpec::gradient_fn) grad_fn,
    void* model_ptr, int n_params,
    const Rcpp::NumericVector& theta0, double sigma_beta,
    const std::vector<double>& inv_metric,
    int n_iter, int n_warmup, int max_treedepth, double adapt_delta,
    int seed, bool verbose, int n_chains, int n_threads,
    const char* spec_name = "tulpaobs_marginal",
    int data_N = -1
) {
    if ((int) theta0.size() != n_params)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), n_params);

    tulpa::LikelihoodSpec lspec;
    lspec.name = spec_name;
    lspec.n_processes = 1;
    lspec.gradient_fn = grad_fn;
    lspec.ll_double = &full_grad_log_post;
    lspec.n_extra_params = n_params;

    tulpa::ModelData data;
    data.N = (data_N >= 0) ? data_N : n_params;
    if (data.N < 1)
        Rcpp::stop("run_tulpa_nuts: ModelData.N must be >= 1 (got %d); the "
                   "double log-posterior is read at observation 0.", data.N);
    data.n_processes = 1;
    data.sigma_beta = sigma_beta;
    data.model_response_data = model_ptr;
    data.likelihood_spec = &lspec;
    data.sharing.init(1);
    data.zi_type = tulpa::ZIType::NONE;
    data.p_zi = 0; data.p_oi = 0;
    tulpa::ProcessData proc;
    proc.p = 0;
    data.processes.push_back(proc);

    tulpa::ParamLayout layout = tulpa::compute_layout(data);
    if (layout.total_params != n_params || layout.extra_offset != 0)
        Rcpp::stop("run_tulpa_nuts: the engine laid out %d parameters from "
                   "offset %d for a %d-parameter target.",
                   layout.total_params, layout.extra_offset, n_params);

    tulpa::set_gradient_mode_str("H");

    if (!inv_metric.empty() && (int) inv_metric.size() != n_params)
        Rcpp::stop("inv_metric length %d != expected %d",
                   (int) inv_metric.size(), n_params);

    if (n_chains < 1) Rcpp::stop("run_tulpa_nuts: n_chains must be >= 1 (got %d).", n_chains);
    std::vector<double> init((std::size_t) n_chains * n_params);
    for (int c = 0; c < n_chains; ++c)
        std::copy(theta0.begin(), theta0.end(), init.begin() + (std::size_t) c * n_params);
    std::vector<double> im;
    if (!inv_metric.empty()) {
        im.resize((std::size_t) n_chains * n_params);
        for (int c = 0; c < n_chains; ++c)
            std::copy(inv_metric.begin(), inv_metric.end(),
                      im.begin() + (std::size_t) c * n_params);
    }

    std::vector<tulpa::NUTSResult> results(n_chains);
    {
        ChainTeamScope team(n_threads);
        tulpa::NUTSChainsFn run_chains = tulpa::get_nuts_chains_fn();
        run_chains(&data, &layout, init.data(), n_params, n_chains, n_iter, n_warmup,
                   max_treedepth, adapt_delta, static_cast<unsigned int>(seed),
                   verbose ? 1 : 0, nullptr, im.empty() ? nullptr : im.data(),
                   nullptr, results.data());
    }

    Rcpp::List chains(n_chains);
    const int np = n_params;
    for (int c = 0; c < n_chains; ++c) {
        tulpa::NUTSResult& result = results[c];
        const int n_samples = result.n_sample;
        Rcpp::NumericMatrix draws(n_samples, np);
        Rcpp::NumericVector lp(n_samples), ap(n_samples);
        Rcpp::IntegerVector div(n_samples), td(n_samples);
        for (int s = 0; s < n_samples; ++s) {
            for (int j = 0; j < np; ++j)
                draws(s, j) = result.samples[(std::size_t) s * np + j];
            lp[s] = result.log_prob[s]; ap[s] = result.accept_prob[s];
            div[s] = result.divergent[s]; td[s] = result.treedepth[s];
        }
        const double epsilon = result.epsilon;
        result.free_buffers();
        chains[c] = Rcpp::List::create(
            Rcpp::Named("draws") = draws, Rcpp::Named("log_prob") = lp,
            Rcpp::Named("accept_prob") = ap, Rcpp::Named("divergent") = div,
            Rcpp::Named("treedepth") = td, Rcpp::Named("epsilon") = epsilon,
            Rcpp::Named("n_params") = np);
    }
    return chains;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NUTS_ENGINE_H
