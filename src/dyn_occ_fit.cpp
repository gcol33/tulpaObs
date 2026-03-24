// dyn_occ_fit.cpp
// Rcpp entry point for multi-season dynamic occupancy model
// MacKenzie et al. (2003): colonization-extinction dynamics

#include <Rcpp.h>
#include <vector>
#include <string>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>

#include "dyn_occ_data.h"
#include "dyn_occ_likelihood.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_dyn_occ_fit(
    Rcpp::IntegerVector y_flat_r,       // Detection history [n_sites * n_seasons * max_visits]
    Rcpp::IntegerVector n_visits_r,     // [n_sites * n_seasons]
    Rcpp::LogicalVector any_detected_r, // [n_sites * n_seasons]
    Rcpp::NumericMatrix X_occ_r,        // Initial occupancy covariates [n_sites x p_occ]
    Rcpp::NumericMatrix X_det_r,        // Detection covariates [n_sites x p_det]
    Rcpp::NumericMatrix X_col_r,        // Colonization covariates [n_sites x p_col]
    Rcpp::NumericMatrix X_ext_r,        // Extinction covariates [n_sites x p_ext]
    int n_sites,
    int n_seasons,
    int max_visits,
    double sigma_beta = 10.0,
    int n_iter = 2000,
    int n_warmup = 1000,
    int max_treedepth = 10,
    double adapt_delta = 0.8,
    int seed = 42,
    bool verbose = true
) {
    const int p_occ = X_occ_r.ncol();
    const int p_det = X_det_r.ncol();
    const int p_col = X_col_r.ncol();
    const int p_ext = X_ext_r.ncol();

    // Build DynOccResponseData
    tulpaOcc::DynOccResponseData dyn;
    dyn.n_sites = n_sites;
    dyn.n_seasons = n_seasons;
    dyn.max_visits = max_visits;
    dyn.y = Rcpp::as<std::vector<int>>(y_flat_r);
    dyn.n_visits = Rcpp::as<std::vector<int>>(n_visits_r);
    dyn.any_detected.resize(n_sites * n_seasons);
    for (int i = 0; i < n_sites * n_seasons; i++) {
        dyn.any_detected[i] = any_detected_r[i];
    }

    // Build LikelihoodSpec
    tulpa::LikelihoodSpec spec;
    spec.name = "dynamic_occupancy";
    spec.n_processes = 4;  // psi1, p, gamma, epsilon
    spec.ll_double = tulpaOcc::dyn_occ_log_likelihood<double>;
    spec.ll_arena  = tulpaOcc::dyn_occ_log_likelihood<tulpa::arena::Var>;
    spec.ll_fwd    = tulpaOcc::dyn_occ_log_likelihood<fwd::Dual>;
    spec.n_extra_params = 0;

    // Build ModelData
    tulpa::ModelData data;
    data.N = n_sites;  // One obs per site; likelihood loops over seasons
    data.n_processes = 4;
    data.sigma_beta = sigma_beta;

    // Process 0: initial occupancy
    tulpa::ProcessData proc_occ;
    proc_occ.p = p_occ;
    proc_occ.X_flat.resize(n_sites * p_occ);
    for (int i = 0; i < n_sites; i++)
        for (int j = 0; j < p_occ; j++)
            proc_occ.X_flat[i * p_occ + j] = X_occ_r(i, j);
    data.processes.push_back(proc_occ);

    // Process 1: detection
    tulpa::ProcessData proc_det;
    proc_det.p = p_det;
    proc_det.X_flat.resize(n_sites * p_det);
    for (int i = 0; i < n_sites; i++)
        for (int j = 0; j < p_det; j++)
            proc_det.X_flat[i * p_det + j] = X_det_r(i, j);
    data.processes.push_back(proc_det);

    // Process 2: colonization
    tulpa::ProcessData proc_col;
    proc_col.p = p_col;
    proc_col.X_flat.resize(n_sites * p_col);
    for (int i = 0; i < n_sites; i++)
        for (int j = 0; j < p_col; j++)
            proc_col.X_flat[i * p_col + j] = X_col_r(i, j);
    data.processes.push_back(proc_col);

    // Process 3: extinction
    tulpa::ProcessData proc_ext;
    proc_ext.p = p_ext;
    proc_ext.X_flat.resize(n_sites * p_ext);
    for (int i = 0; i < n_sites; i++)
        for (int j = 0; j < p_ext; j++)
            proc_ext.X_flat[i * p_ext + j] = X_ext_r(i, j);
    data.processes.push_back(proc_ext);

    data.model_response_data = &dyn;
    data.likelihood_spec = &spec;
    data.sharing.init(4);

    // ZI/OI not used
    data.zi_type = tulpa::ZIType::NONE;
    data.p_zi = 0;
    data.p_oi = 0;
    data.zi_prior_sd = 1.0;
    data.oi_prior_sd = 1.0;

    // Compute layout and run NUTS
    tulpa::ParamLayout layout = tulpa::compute_layout(data);
    int n_params = layout.total_params;
    std::vector<double> init(n_params, 0.0);

    tulpa::NUTSFn run_nuts = tulpa::get_nuts_fn();
    tulpa::NUTSResult result = {};

    run_nuts(
        &data, &layout, init.data(), n_params,
        n_iter, n_warmup, max_treedepth, adapt_delta,
        static_cast<unsigned int>(seed),
        verbose ? 1 : 0,
        &result
    );

    int n_samples = result.n_sample;

    // Convert to R
    NumericMatrix draws(n_samples, n_params);
    NumericVector lp(n_samples), ap(n_samples);
    IntegerVector div(n_samples), td(n_samples);

    for (int s = 0; s < n_samples; s++) {
        for (int j = 0; j < n_params; j++)
            draws(s, j) = result.samples[s * n_params + j];
        lp[s] = result.log_prob[s];
        ap[s] = result.accept_prob[s];
        div[s] = result.divergent[s];
        td[s] = result.treedepth[s];
    }

    double epsilon = result.epsilon;
    result.free_buffers();

    // Column names
    CharacterVector col_names(n_params);
    int idx = 0;
    for (int j = 0; j < p_occ; j++)
        col_names[idx++] = "beta_psi1[" + std::to_string(j + 1) + "]";
    for (int j = 0; j < p_det; j++)
        col_names[idx++] = "beta_p[" + std::to_string(j + 1) + "]";
    for (int j = 0; j < p_col; j++)
        col_names[idx++] = "beta_gamma[" + std::to_string(j + 1) + "]";
    for (int j = 0; j < p_ext; j++)
        col_names[idx++] = "beta_epsilon[" + std::to_string(j + 1) + "]";
    for (; idx < n_params; idx++)
        col_names[idx] = "param[" + std::to_string(idx + 1) + "]";
    Rcpp::colnames(draws) = col_names;

    NumericVector means(n_params, 0.0);
    for (int s = 0; s < n_samples; s++)
        for (int j = 0; j < n_params; j++)
            means[j] += draws(s, j) / n_samples;
    means.names() = col_names;

    return List::create(
        Named("draws") = draws,
        Named("means") = means,
        Named("n_samples") = n_samples,
        Named("n_params") = n_params,
        Named("n_sites") = n_sites,
        Named("n_seasons") = n_seasons,
        Named("max_visits") = max_visits,
        Named("p_occ") = p_occ,
        Named("p_det") = p_det,
        Named("p_col") = p_col,
        Named("p_ext") = p_ext,
        Named("col_names") = col_names,
        Named("log_prob") = lp,
        Named("accept_prob") = ap,
        Named("divergent") = div,
        Named("treedepth") = td,
        Named("epsilon") = epsilon
    );
}
