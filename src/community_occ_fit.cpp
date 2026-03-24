// community_occ_fit.cpp
// Rcpp entry point for community (multi-species) occupancy model
// Reuses single-season occupancy likelihood with species-level random effects

#include <Rcpp.h>
#include <vector>
#include <string>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>

#include "occ_data.h"
#include "occ_likelihood.h"

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_community_occ_fit(
    Rcpp::IntegerMatrix y_r,            // Detection history [N x max_visits], N = n_sites * n_species
    Rcpp::NumericMatrix X_occ_r,        // Occupancy covariates [N x p_occ]
    Rcpp::NumericMatrix X_det_r,        // Detection covariates [N x p_det]
    Rcpp::IntegerVector species_group_r, // Species group (1-based) [N]
    int n_species,
    double sigma_beta = 10.0,
    double sigma_re_scale = 1.0,
    int n_iter = 2000,
    int n_warmup = 1000,
    int max_treedepth = 10,
    double adapt_delta = 0.8,
    int seed = 42,
    bool verbose = true
) {
    const int N = y_r.nrow();
    const int max_visits = y_r.ncol();
    const int p_occ = X_occ_r.ncol();
    const int p_det = X_det_r.ncol();

    // Build OccResponseData (same struct as single-season, just N = n_sites * n_species)
    tulpaOcc::OccResponseData occ;
    occ.n_sites = N;  // Each site-species pair is an "observation"
    occ.max_visits = max_visits;
    occ.y.resize(N * max_visits);
    occ.n_visits.resize(N);
    occ.any_detected.resize(N, false);
    occ.n_detections.resize(N, 0);

    for (int i = 0; i < N; i++) {
        int nv = 0;
        for (int j = 0; j < max_visits; j++) {
            int val = y_r(i, j);
            occ.y[i * max_visits + j] = val;
            if (val >= 0) {
                nv++;
                if (val == 1) {
                    occ.any_detected[i] = true;
                    occ.n_detections[i]++;
                }
            }
        }
        occ.n_visits[i] = nv;
    }

    // Build LikelihoodSpec (same as single-season — community structure is in RE)
    tulpa::LikelihoodSpec spec;
    spec.name = "community_occupancy";
    spec.n_processes = 2;
    spec.ll_double = tulpaOcc::occ_log_likelihood<double>;
    spec.ll_arena  = tulpaOcc::occ_log_likelihood<tulpa::arena::Var>;
    spec.ll_fwd    = tulpaOcc::occ_log_likelihood<fwd::Dual>;
    spec.residual_fn = tulpaOcc::occ_residual;
    spec.n_extra_params = 0;

    // Build ModelData
    tulpa::ModelData data;
    data.N = N;
    data.n_processes = 2;
    data.sigma_beta = sigma_beta;
    data.sigma_re_scale = sigma_re_scale;

    // Process 0: occupancy
    tulpa::ProcessData proc_occ;
    proc_occ.p = p_occ;
    proc_occ.X_flat.resize(N * p_occ);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < p_occ; j++)
            proc_occ.X_flat[i * p_occ + j] = X_occ_r(i, j);
    data.processes.push_back(proc_occ);

    // Process 1: detection
    tulpa::ProcessData proc_det;
    proc_det.p = p_det;
    proc_det.X_flat.resize(N * p_det);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < p_det; j++)
            proc_det.X_flat[i * p_det + j] = X_det_r(i, j);
    data.processes.push_back(proc_det);

    data.model_response_data = &occ;
    data.likelihood_spec = &spec;
    data.sharing.init(2);

    // Species random effects — single-term legacy path
    data.re_group.resize(N);
    for (int i = 0; i < N; i++) {
        data.re_group[i] = species_group_r[i];  // 1-based
    }
    data.n_re_groups = n_species;
    data.n_re_terms = 0;  // Legacy single-term
    data.total_re_groups = n_species;
    data.has_re_slopes = false;
    data.has_re_correlated_slopes = false;
    data.total_re_params = n_species;
    data.total_sigma_params = 1;
    data.total_chol_params = 0;
    data.re_parameterization = 1;  // Non-centered

    // RE enters both psi and p
    data.sharing.re[0] = true;
    data.sharing.re[1] = true;

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
        col_names[idx++] = "beta_occ[" + std::to_string(j + 1) + "]";
    for (int j = 0; j < p_det; j++)
        col_names[idx++] = "beta_det[" + std::to_string(j + 1) + "]";
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
        Named("n_species") = n_species,
        Named("p_occ") = p_occ,
        Named("p_det") = p_det,
        Named("col_names") = col_names,
        Named("log_prob") = lp,
        Named("accept_prob") = ap,
        Named("divergent") = div,
        Named("treedepth") = td,
        Named("epsilon") = epsilon
    );
}
