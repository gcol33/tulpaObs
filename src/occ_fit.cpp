// occ_fit.cpp
// Rcpp entry point for single-season occupancy model
// Uses tulpa's NUTS backend via R_GetCCallable cross-package API.

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

// ============================================================================
// Helper: populate spatial fields on ModelData from R list
// ============================================================================
static void populate_spatial(tulpa::ModelData& data, Rcpp::List sp, int n_sites) {
    std::string type = Rcpp::as<std::string>(sp["type"]);

    if (type == "none") return;

    // Sharing configuration
    bool shared_occ = Rcpp::as<bool>(sp["spatial_shared_occ"]);
    bool shared_det = Rcpp::as<bool>(sp["spatial_shared_det"]);
    data.sharing.spatial[0] = shared_occ;
    data.sharing.spatial[1] = shared_det;

    if (type == "icar" || type == "bym2") {
        data.spatial_type = (type == "icar") ? tulpa::SpatialType::ICAR
                                              : tulpa::SpatialType::BYM2;
        data.n_spatial_units = Rcpp::as<int>(sp["n_units"]);
        data.adj_row_ptr = Rcpp::as<std::vector<int>>(sp["adj_row_ptr"]);
        data.adj_col_idx = Rcpp::as<std::vector<int>>(sp["adj_col_idx"]);
        data.n_neighbors = Rcpp::as<std::vector<int>>(sp["n_neighbors"]);

        // 1:1 mapping obs -> spatial unit (occupancy = one obs per site)
        data.spatial_group.resize(n_sites);
        for (int i = 0; i < n_sites; i++) data.spatial_group[i] = i + 1;

        if (type == "bym2") {
            data.bym2_scale_factor = Rcpp::as<double>(sp["scale_factor"]);
        }

    } else if (type == "gp") {
        data.spatial_type = tulpa::SpatialType::GP;
        data.has_gp = true;

        int n_obs = Rcpp::as<int>(sp["n_obs"]);
        int nn = Rcpp::as<int>(sp["nn"]);
        data.gp_data.n_obs = n_obs;
        data.gp_data.nn = nn;
        data.gp_data.coords = Rcpp::as<std::vector<double>>(sp["coords"]);
        data.gp_data.nn_idx = Rcpp::as<std::vector<int>>(sp["nn_idx"]);
        data.gp_data.nn_dist = Rcpp::as<std::vector<double>>(sp["nn_dist"]);
        data.gp_data.nn_neighbor_dist = Rcpp::as<std::vector<double>>(sp["nn_neighbor_dist"]);

        // Convert 1-based R indices to 0-based C++
        std::vector<int> nn_order_r = Rcpp::as<std::vector<int>>(sp["nn_order"]);
        std::vector<int> nn_order_inv_r = Rcpp::as<std::vector<int>>(sp["nn_order_inv"]);
        data.gp_data.nn_order.resize(nn_order_r.size());
        data.gp_data.nn_order_inv.resize(nn_order_inv_r.size());
        for (size_t i = 0; i < nn_order_r.size(); i++) {
            data.gp_data.nn_order[i] = nn_order_r[i] - 1;
            data.gp_data.nn_order_inv[i] = nn_order_inv_r[i] - 1;
        }

        // 1:1 obs-to-location mapping
        data.gp_data.obs_to_loc.resize(n_sites);
        for (int i = 0; i < n_sites; i++) data.gp_data.obs_to_loc[i] = i;

        // Covariance type
        std::string cov_str = Rcpp::as<std::string>(sp["cov_type"]);
        if (cov_str == "matern") data.gp_data.cov_type = tulpa::CovType::MATERN;
        else if (cov_str == "gaussian") data.gp_data.cov_type = tulpa::CovType::GAUSSIAN;
        else if (cov_str == "spherical") data.gp_data.cov_type = tulpa::CovType::SPHERICAL;
        else data.gp_data.cov_type = tulpa::CovType::EXPONENTIAL;

        data.gp_data.nu = Rcpp::as<double>(sp["nu"]);
        data.gp_data.shared = true;
        data.gp_data.solver_config.n_obs = n_obs;

        data.gp_sigma2_prior_U = Rcpp::as<double>(sp["sigma2_prior_U"]);
        data.gp_sigma2_prior_alpha = Rcpp::as<double>(sp["sigma2_prior_alpha"]);
        data.gp_phi_prior_lower = Rcpp::as<double>(sp["phi_prior_lower"]);
        data.gp_phi_prior_upper = Rcpp::as<double>(sp["phi_prior_upper"]);
        data.gp_parameterization = 1;  // Non-centered (default for occupancy)

    } else if (type == "multiscale_gp") {
        data.spatial_type = tulpa::SpatialType::MULTISCALE_GP;
        data.has_multiscale_gp = true;

        int n_obs = Rcpp::as<int>(sp["n_obs"]);
        data.multiscale_gp_data.n_obs = n_obs;
        data.multiscale_gp_data.coords = Rcpp::as<std::vector<double>>(sp["coords"]);

        // 1:1 obs-to-location mapping
        data.multiscale_gp_data.obs_to_loc.resize(n_sites);
        for (int i = 0; i < n_sites; i++) data.multiscale_gp_data.obs_to_loc[i] = i;

        // Local scale
        data.multiscale_gp_data.nn_local = Rcpp::as<int>(sp["nn_local"]);
        data.multiscale_gp_data.nn_idx_local = Rcpp::as<std::vector<int>>(sp["nn_idx_local"]);
        data.multiscale_gp_data.nn_dist_local = Rcpp::as<std::vector<double>>(sp["nn_dist_local"]);
        data.multiscale_gp_data.nn_neighbor_dist_local = Rcpp::as<std::vector<double>>(sp["nn_neighbor_dist_local"]);

        std::vector<int> lo_r = Rcpp::as<std::vector<int>>(sp["nn_order_local"]);
        std::vector<int> loi_r = Rcpp::as<std::vector<int>>(sp["nn_order_inv_local"]);
        data.multiscale_gp_data.nn_order_local.resize(lo_r.size());
        data.multiscale_gp_data.nn_order_inv_local.resize(loi_r.size());
        for (size_t i = 0; i < lo_r.size(); i++) {
            data.multiscale_gp_data.nn_order_local[i] = lo_r[i] - 1;
            data.multiscale_gp_data.nn_order_inv_local[i] = loi_r[i] - 1;
        }

        // Regional scale
        data.multiscale_gp_data.nn_regional = Rcpp::as<int>(sp["nn_regional"]);
        data.multiscale_gp_data.nn_idx_regional = Rcpp::as<std::vector<int>>(sp["nn_idx_regional"]);
        data.multiscale_gp_data.nn_dist_regional = Rcpp::as<std::vector<double>>(sp["nn_dist_regional"]);
        data.multiscale_gp_data.nn_neighbor_dist_regional = Rcpp::as<std::vector<double>>(sp["nn_neighbor_dist_regional"]);

        std::vector<int> ro_r = Rcpp::as<std::vector<int>>(sp["nn_order_regional"]);
        std::vector<int> roi_r = Rcpp::as<std::vector<int>>(sp["nn_order_inv_regional"]);
        data.multiscale_gp_data.nn_order_regional.resize(ro_r.size());
        data.multiscale_gp_data.nn_order_inv_regional.resize(roi_r.size());
        for (size_t i = 0; i < ro_r.size(); i++) {
            data.multiscale_gp_data.nn_order_regional[i] = ro_r[i] - 1;
            data.multiscale_gp_data.nn_order_inv_regional[i] = roi_r[i] - 1;
        }

        // Covariance
        std::string cov_str = Rcpp::as<std::string>(sp["cov_type"]);
        if (cov_str == "matern") data.multiscale_gp_data.cov_type = tulpa::CovType::MATERN;
        else if (cov_str == "gaussian") data.multiscale_gp_data.cov_type = tulpa::CovType::GAUSSIAN;
        else data.multiscale_gp_data.cov_type = tulpa::CovType::EXPONENTIAL;
        data.multiscale_gp_data.nu = Rcpp::as<double>(sp["nu"]);
        data.multiscale_gp_data.shared = true;

        // Range bounds
        data.multiscale_gp_data.range_local_lower = Rcpp::as<double>(sp["range_local_lower"]);
        data.multiscale_gp_data.range_local_upper = Rcpp::as<double>(sp["range_local_upper"]);
        data.multiscale_gp_data.range_regional_lower = Rcpp::as<double>(sp["range_regional_lower"]);
        data.multiscale_gp_data.range_regional_upper = Rcpp::as<double>(sp["range_regional_upper"]);

        // Priors
        data.ms_sigma2_local_prior_U = Rcpp::as<double>(sp["sigma2_local_prior_U"]);
        data.ms_sigma2_local_prior_alpha = Rcpp::as<double>(sp["sigma2_local_prior_alpha"]);
        data.ms_sigma2_regional_prior_U = Rcpp::as<double>(sp["sigma2_regional_prior_U"]);
        data.ms_sigma2_regional_prior_alpha = Rcpp::as<double>(sp["sigma2_regional_prior_alpha"]);

        // Lengthscale priors (geometric mean of range bounds)
        double rl = Rcpp::as<double>(sp["range_local_lower"]);
        double ru = Rcpp::as<double>(sp["range_local_upper"]);
        data.ms_log_ls_local_mean = 0.5 * (std::log(rl) + std::log(ru));
        data.ms_log_ls_local_sd = 0.5;
        rl = Rcpp::as<double>(sp["range_regional_lower"]);
        ru = Rcpp::as<double>(sp["range_regional_upper"]);
        data.ms_log_ls_regional_mean = 0.5 * (std::log(rl) + std::log(ru));
        data.ms_log_ls_regional_sd = 0.5;
    }
}

// ============================================================================
// Helper: build OccResponseData from R matrices
// ============================================================================
static tulpaOcc::OccResponseData build_occ_response(
    Rcpp::IntegerMatrix y_r, int n_sites, int max_visits,
    Rcpp::Nullable<Rcpp::NumericMatrix> X_det_visit_r, int& p_det_visit_out
) {
    tulpaOcc::OccResponseData occ;
    occ.n_sites = n_sites;
    occ.max_visits = max_visits;
    occ.y.resize(n_sites * max_visits);
    occ.n_visits.resize(n_sites);
    occ.any_detected.resize(n_sites, false);
    occ.n_detections.resize(n_sites, 0);

    for (int i = 0; i < n_sites; i++) {
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

    p_det_visit_out = 0;
    if (X_det_visit_r.isNotNull()) {
        NumericMatrix Xv = Rcpp::as<NumericMatrix>(X_det_visit_r);
        p_det_visit_out = Xv.ncol();
        occ.p_det_visit = p_det_visit_out;
        occ.X_det_visit.resize(n_sites * max_visits * p_det_visit_out);
        for (int i = 0; i < n_sites; i++) {
            for (int j = 0; j < max_visits; j++) {
                int row = i * max_visits + j;
                for (int c = 0; c < p_det_visit_out; c++) {
                    occ.X_det_visit[i * max_visits * p_det_visit_out + j * p_det_visit_out + c] =
                        Xv(row, c);
                }
            }
        }
    }

    return occ;
}

// ============================================================================
// Helper: run NUTS and convert results to R list
// ============================================================================
static Rcpp::List run_nuts_and_collect(
    tulpa::ModelData& data, tulpa::ParamLayout& layout,
    int n_iter, int n_warmup, int max_treedepth, double adapt_delta,
    int seed, bool verbose,
    int p_occ, int p_det, int p_det_visit, int n_sites, int max_visits
) {
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

    NumericMatrix draws(n_samples, n_params);
    NumericVector lp(n_samples), ap(n_samples);
    IntegerVector div(n_samples), td(n_samples);

    for (int s = 0; s < n_samples; s++) {
        for (int j = 0; j < n_params; j++) {
            draws(s, j) = result.samples[s * n_params + j];
        }
        lp[s] = result.log_prob[s];
        ap[s] = result.accept_prob[s];
        div[s] = result.divergent[s];
        td[s] = result.treedepth[s];
    }

    double epsilon = result.epsilon;
    result.free_buffers();

    // Column names for fixed effects
    CharacterVector col_names(n_params);
    int idx = 0;
    for (int j = 0; j < p_occ; j++)
        col_names[idx++] = "beta_occ[" + std::to_string(j + 1) + "]";
    for (int j = 0; j < p_det; j++)
        col_names[idx++] = "beta_det[" + std::to_string(j + 1) + "]";
    for (int j = 0; j < p_det_visit; j++)
        col_names[idx++] = "beta_det_visit[" + std::to_string(j + 1) + "]";
    // Remaining params get generic names
    for (; idx < n_params; idx++)
        col_names[idx] = "param[" + std::to_string(idx + 1) + "]";
    Rcpp::colnames(draws) = col_names;

    // Posterior means
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
        Named("max_visits") = max_visits,
        Named("p_occ") = p_occ,
        Named("p_det") = p_det,
        Named("p_det_visit") = p_det_visit,
        Named("col_names") = col_names,
        Named("log_prob") = lp,
        Named("accept_prob") = ap,
        Named("divergent") = div,
        Named("treedepth") = td,
        Named("epsilon") = epsilon
    );
}

// ============================================================================
// Rcpp entry point: single-season occupancy
// ============================================================================

// [[Rcpp::export]]
Rcpp::List cpp_occ_fit(
    Rcpp::IntegerMatrix y_r,
    Rcpp::NumericMatrix X_occ_r,
    Rcpp::NumericMatrix X_det_r,
    Rcpp::Nullable<Rcpp::NumericMatrix> X_det_visit_r = R_NilValue,
    Rcpp::Nullable<Rcpp::List> spatial_params_r = R_NilValue,
    double sigma_beta = 10.0,
    int n_iter = 2000,
    int n_warmup = 1000,
    int max_treedepth = 10,
    double adapt_delta = 0.8,
    int seed = 42,
    bool verbose = true
) {
    const int n_sites = y_r.nrow();
    const int max_visits = y_r.ncol();
    const int p_occ = X_occ_r.ncol();
    const int p_det = X_det_r.ncol();

    // Build OccResponseData
    int p_det_visit = 0;
    tulpaOcc::OccResponseData occ = build_occ_response(
        y_r, n_sites, max_visits, X_det_visit_r, p_det_visit);

    // Build LikelihoodSpec — all three gradient modes
    tulpa::LikelihoodSpec spec;
    spec.name = "occupancy";
    spec.n_processes = 2;
    spec.ll_double = tulpaOcc::occ_log_likelihood<double>;
    spec.ll_arena  = tulpaOcc::occ_log_likelihood<tulpa::arena::Var>;
    spec.ll_fwd    = tulpaOcc::occ_log_likelihood<fwd::Dual>;
    spec.residual_fn = tulpaOcc::occ_residual;
    spec.n_extra_params = p_det_visit;

    // Build ModelData
    tulpa::ModelData data;
    data.N = n_sites;
    data.n_processes = 2;
    data.sigma_beta = sigma_beta;

    // Process 0: occupancy
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

    data.model_response_data = &occ;
    data.likelihood_spec = &spec;
    data.sharing.init(2);

    // ZI/OI not used
    data.zi_type = tulpa::ZIType::NONE;
    data.p_zi = 0;
    data.p_oi = 0;
    data.zi_prior_sd = 1.0;
    data.oi_prior_sd = 1.0;

    // Spatial (optional)
    if (spatial_params_r.isNotNull()) {
        Rcpp::List sp = Rcpp::as<Rcpp::List>(spatial_params_r);
        populate_spatial(data, sp, n_sites);
    }

    // Compute layout and run
    tulpa::ParamLayout layout = tulpa::compute_layout(data);

    return run_nuts_and_collect(
        data, layout, n_iter, n_warmup, max_treedepth, adapt_delta,
        seed, verbose, p_occ, p_det, p_det_visit, n_sites, max_visits);
}
