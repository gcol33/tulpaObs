// populate_helpers.h
// Shared helpers for populating tulpa::ModelData from R lists.
// Used by the unified cpp_occu_fit entry point.

#ifndef TULPAOCC_POPULATE_HELPERS_H
#define TULPAOCC_POPULATE_HELPERS_H

#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <tulpa/model_data.h>

namespace tulpaObs {

// ============================================================================
// Populate spatial fields on ModelData from R list
// ============================================================================
inline void populate_spatial(tulpa::ModelData& data, Rcpp::List sp, int n_units) {
    std::string type = Rcpp::as<std::string>(sp["type"]);

    if (type == "none") return;

    // Sharing configuration
    bool shared_occ = Rcpp::as<bool>(sp["spatial_shared_occ"]);
    bool shared_det = Rcpp::as<bool>(sp["spatial_shared_det"]);
    data.sharing.spatial[0] = shared_occ;
    data.sharing.spatial[1] = shared_det;

    // For models with more than 2 processes (dynamic: 4), propagate sharing
    // to additional processes (default: not shared)
    for (size_t k = 2; k < data.sharing.spatial.size(); k++) {
        data.sharing.spatial[k] = false;
    }

    // Check for extended sharing (dynamic models: gamma, epsilon)
    if (sp.containsElementNamed("spatial_shared_col")) {
        if (data.sharing.spatial.size() > 2)
            data.sharing.spatial[2] = Rcpp::as<bool>(sp["spatial_shared_col"]);
    }
    if (sp.containsElementNamed("spatial_shared_ext")) {
        if (data.sharing.spatial.size() > 3)
            data.sharing.spatial[3] = Rcpp::as<bool>(sp["spatial_shared_ext"]);
    }

    if (type == "icar" || type == "bym2") {
        data.spatial_type = (type == "icar") ? tulpa::SpatialType::ICAR
                                              : tulpa::SpatialType::BYM2;
        data.n_spatial_units = Rcpp::as<int>(sp["n_units"]);
        data.adj_row_ptr = Rcpp::as<std::vector<int>>(sp["adj_row_ptr"]);
        data.adj_col_idx = Rcpp::as<std::vector<int>>(sp["adj_col_idx"]);
        data.n_neighbors = Rcpp::as<std::vector<int>>(sp["n_neighbors"]);

        // 1:1 mapping obs -> spatial unit
        data.spatial_group.resize(n_units);
        for (int i = 0; i < n_units; i++) data.spatial_group[i] = i + 1;

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
        data.gp_data.obs_to_loc.resize(n_units);
        for (int i = 0; i < n_units; i++) data.gp_data.obs_to_loc[i] = i;

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
        data.gp_phi_prior_U = Rcpp::as<double>(sp["phi_prior_U"]);
        data.gp_phi_prior_alpha = Rcpp::as<double>(sp["phi_prior_alpha"]);
        // Centered: compute_gp_spatial_prior evaluates the sampled parameters
        // directly as the field w, and has no z -> w branch, so the stored
        // draws have to be those same values. Setting 1 makes the chain writer
        // forward-transform them as if they were whitened, corrupting every
        // stored field draw while leaving the trajectory untouched.
        data.gp_parameterization = 0;

    } else if (type == "multiscale_gp") {
        data.spatial_type = tulpa::SpatialType::MULTISCALE_GP;
        data.has_multiscale_gp = true;

        int n_obs = Rcpp::as<int>(sp["n_obs"]);
        data.multiscale_gp_data.n_obs = n_obs;
        data.multiscale_gp_data.coords = Rcpp::as<std::vector<double>>(sp["coords"]);

        // 1:1 obs-to-location mapping
        data.multiscale_gp_data.obs_to_loc.resize(n_units);
        for (int i = 0; i < n_units; i++) data.multiscale_gp_data.obs_to_loc[i] = i;

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
// Populate temporal fields on ModelData from R list
// ============================================================================
inline void populate_temporal(tulpa::ModelData& data, Rcpp::List temp_spec) {
    std::string type = Rcpp::as<std::string>(temp_spec["type"]);
    if (type == "none") return;

    if (type == "rw1")      data.temporal_type = tulpa::TemporalType::RW1;
    else if (type == "rw2") data.temporal_type = tulpa::TemporalType::RW2;
    else if (type == "ar1") data.temporal_type = tulpa::TemporalType::AR1;
    else if (type == "iid") data.temporal_type = tulpa::TemporalType::IID;
    else Rcpp::stop("Unknown temporal type: %s", type.c_str());

    data.temporal_time_idx = Rcpp::as<std::vector<int>>(temp_spec["time_idx"]);
    data.n_times = Rcpp::as<int>(temp_spec["n_times"]);
    data.n_temporal_groups = 1;
    data.n_temporal_params = data.n_times * data.n_temporal_groups;

    if (temp_spec.containsElementNamed("group_idx")) {
        data.temporal_group_idx = Rcpp::as<std::vector<int>>(temp_spec["group_idx"]);
        data.n_temporal_groups = Rcpp::as<int>(temp_spec["n_groups"]);
        data.n_temporal_params = data.n_times * data.n_temporal_groups;
    } else {
        // Default: all observations in group 1
        int N = data.temporal_time_idx.size();
        data.temporal_group_idx.resize(N, 1);
    }

    if (temp_spec.containsElementNamed("cyclic"))
        data.temporal_cyclic = Rcpp::as<bool>(temp_spec["cyclic"]);

    if (temp_spec.containsElementNamed("tau_shape"))
        data.tau_temporal_shape = Rcpp::as<double>(temp_spec["tau_shape"]);
    if (temp_spec.containsElementNamed("tau_rate"))
        data.tau_temporal_rate = Rcpp::as<double>(temp_spec["tau_rate"]);

    // Sharing: which processes get temporal effect
    Rcpp::LogicalVector shared = Rcpp::as<Rcpp::LogicalVector>(temp_spec["shared"]);
    for (int k = 0; k < shared.size() && k < (int)data.sharing.temporal.size(); k++) {
        data.sharing.temporal[k] = shared[k];
    }
}

// ============================================================================
// Populate multi-term random effects on ModelData from R list
// ============================================================================
inline void populate_re(tulpa::ModelData& data, Rcpp::List re_spec) {
    int n_terms = Rcpp::as<int>(re_spec["n_terms"]);
    if (n_terms == 0) return;

    data.n_re_terms = n_terms;
    data.re_parameterization = 1;  // Non-centered

    // Per-term group assignments: list of integer vectors
    Rcpp::List group_list = Rcpp::as<Rcpp::List>(re_spec["groups"]);
    Rcpp::IntegerVector n_groups_vec = Rcpp::as<Rcpp::IntegerVector>(re_spec["n_groups"]);

    int N = data.N;
    data.re_group_multi_flat.resize(N * n_terms);
    data.re_n_groups_multi.resize(n_terms);
    data.re_offsets.resize(n_terms);

    int total_groups = 0;
    for (int t = 0; t < n_terms; t++) {
        Rcpp::IntegerVector grp = Rcpp::as<Rcpp::IntegerVector>(group_list[t]);
        data.re_n_groups_multi[t] = n_groups_vec[t];
        data.re_offsets[t] = total_groups;
        for (int i = 0; i < N; i++) {
            data.re_group_multi_flat[i * n_terms + t] = grp[i];
        }
        total_groups += n_groups_vec[t];
    }
    data.total_re_groups = total_groups;
    data.total_re_params = total_groups;  // Updated below if slopes
    data.total_sigma_params = n_terms;
    data.total_chol_params = 0;

    // Also set legacy single-term fields — compute_param_layout uses these
    // when n_re_terms <= 1 (falls into legacy path).
    if (n_terms == 1) {
        data.re_group.resize(N);
        Rcpp::IntegerVector grp0 = Rcpp::as<Rcpp::IntegerVector>(group_list[0]);
        for (int i = 0; i < N; i++) data.re_group[i] = grp0[i];
        data.n_re_groups = n_groups_vec[0];
    }

    // Per-term intercept flag. A block carries the implicit group intercept
    // (coef 0, z = 1) unless it is slope-only (`(0 + x | g)`). Default all-on
    // when the R spec omits it.
    data.re_has_intercept.assign(n_terms, 1);
    if (re_spec.containsElementNamed("re_has_intercept")) {
        Rcpp::IntegerVector hi = Rcpp::as<Rcpp::IntegerVector>(re_spec["re_has_intercept"]);
        for (int t = 0; t < n_terms && t < hi.size(); t++)
            data.re_has_intercept[t] = hi[t];
    }

    // Random slopes (optional)
    data.has_re_slopes = false;
    data.has_re_correlated_slopes = false;
    if (re_spec.containsElementNamed("has_slopes") && Rcpp::as<bool>(re_spec["has_slopes"])) {
        data.has_re_slopes = true;
        Rcpp::IntegerVector n_coefs = Rcpp::as<Rcpp::IntegerVector>(re_spec["n_coefs"]);
        data.re_n_coefs.resize(n_terms);
        data.re_n_slopes.resize(n_terms);
        data.re_slope_matrices.resize(n_terms);

        int total_params = 0;
        int total_sigma = 0;
        for (int t = 0; t < n_terms; t++) {
            data.re_n_coefs[t] = n_coefs[t];
            // Slopes = coefs minus the intercept (0 slopes dropped when none).
            data.re_n_slopes[t] = n_coefs[t] - (data.re_has_intercept[t] ? 1 : 0);
            total_params += data.re_n_groups_multi[t] * n_coefs[t];
            total_sigma += n_coefs[t];

            // Slope design matrix [N x n_slopes]
            if (data.re_n_slopes[t] > 0 && re_spec.containsElementNamed("slope_matrices")) {
                Rcpp::List slope_mats = Rcpp::as<Rcpp::List>(re_spec["slope_matrices"]);
                Rcpp::NumericMatrix Xs = Rcpp::as<Rcpp::NumericMatrix>(slope_mats[t]);
                int ns = Xs.ncol();
                data.re_slope_matrices[t].resize(N * ns);
                for (int i = 0; i < N; i++)
                    for (int j = 0; j < ns; j++)
                        data.re_slope_matrices[t][i * ns + j] = Xs(i, j);
            }
        }
        data.total_re_params = total_params;
        data.total_sigma_params = total_sigma;

        // Correlated slopes via Cholesky, per term. `correlated` is a length-
        // n_terms 0/1 vector from build_re_spec(); a term is correlated only
        // when it asked for it AND has more than one coefficient. Always size
        // re_correlated / re_n_chol to n_terms (the param layout indexes them
        // per term, even for fully-uncorrelated `||` blocks).
        {
            data.re_correlated.assign(n_terms, false);
            data.re_n_chol.assign(n_terms, 0);
            int total_chol = 0;
            bool any_corr = false;
            Rcpp::IntegerVector corr;
            if (re_spec.containsElementNamed("correlated"))
                corr = Rcpp::as<Rcpp::IntegerVector>(re_spec["correlated"]);
            for (int t = 0; t < n_terms; t++) {
                int k = data.re_n_coefs[t];
                bool ct = (t < corr.size()) && (corr[t] != 0) && (k > 1);
                if (ct) {
                    data.re_correlated[t] = true;
                    // Strictly-lower triangle: the tanh-Cholesky prior in tulpa
                    // parameterizes k*(k-1)/2 off-diagonal entries and derives
                    // the diagonal from the unit-norm constraint (see
                    // tulpa_priors_re.h).
                    data.re_n_chol[t] = k * (k - 1) / 2;
                    total_chol += data.re_n_chol[t];
                    any_corr = true;
                }
            }
            data.has_re_correlated_slopes = any_corr;
            data.total_chol_params = total_chol;
        }
    }

    // Sharing: which processes get RE
    Rcpp::LogicalVector shared = Rcpp::as<Rcpp::LogicalVector>(re_spec["shared"]);
    for (int k = 0; k < shared.size() && k < (int)data.sharing.re.size(); k++) {
        data.sharing.re[k] = shared[k];
    }

    if (re_spec.containsElementNamed("sigma_re_scale"))
        data.sigma_re_scale = Rcpp::as<double>(re_spec["sigma_re_scale"]);
}

// ============================================================================
// Populate spatially-varying coefficients on ModelData from R list
// ============================================================================
inline void populate_svc(tulpa::ModelData& data, Rcpp::List svc_spec) {
    data.has_svc = true;

    data.svc_data.n_obs = Rcpp::as<int>(svc_spec["n_obs"]);
    data.svc_data.n_svc = Rcpp::as<int>(svc_spec["n_svc"]);
    data.svc_data.nn = Rcpp::as<int>(svc_spec["nn"]);
    data.svc_data.coords = Rcpp::as<std::vector<double>>(svc_spec["coords"]);
    data.svc_data.svc_indices = Rcpp::as<std::vector<int>>(svc_spec["svc_indices"]);
    data.svc_data.X_svc = Rcpp::as<std::vector<double>>(svc_spec["X_svc"]);
    data.svc_data.nn_idx = Rcpp::as<std::vector<int>>(svc_spec["nn_idx"]);
    data.svc_data.nn_dist = Rcpp::as<std::vector<double>>(svc_spec["nn_dist"]);

    // Convert 1-based R indices to 0-based C++
    std::vector<int> order_r = Rcpp::as<std::vector<int>>(svc_spec["nn_order"]);
    std::vector<int> order_inv_r = Rcpp::as<std::vector<int>>(svc_spec["nn_order_inv"]);
    data.svc_data.nn_order.resize(order_r.size());
    data.svc_data.nn_order_inv.resize(order_inv_r.size());
    for (size_t i = 0; i < order_r.size(); i++) {
        data.svc_data.nn_order[i] = order_r[i] - 1;
        data.svc_data.nn_order_inv[i] = order_inv_r[i] - 1;
    }

    // Covariance type
    std::string cov_str = Rcpp::as<std::string>(svc_spec["cov_type"]);
    if (cov_str == "matern") data.svc_data.cov_type = tulpa::CovType::MATERN;
    else if (cov_str == "gaussian") data.svc_data.cov_type = tulpa::CovType::GAUSSIAN;
    else data.svc_data.cov_type = tulpa::CovType::EXPONENTIAL;

    // Priors
    if (svc_spec.containsElementNamed("sigma2_prior_scale"))
        data.svc_sigma2_prior_scale = Rcpp::as<double>(svc_spec["sigma2_prior_scale"]);
    // Unconditional: svc() requires prior_range, so the spec always carries the
    // anchors. Leaving them at the -1.0 sentinel instead would defer the failure
    // to tulpa's layout gate, which cannot name the term argument.
    data.svc_phi_prior_U = Rcpp::as<double>(svc_spec["phi_prior_U"]);
    data.svc_phi_prior_alpha = Rcpp::as<double>(svc_spec["phi_prior_alpha"]);

    // Sharing: which processes get SVC
    Rcpp::LogicalVector shared = Rcpp::as<Rcpp::LogicalVector>(svc_spec["shared"]);
    for (int k = 0; k < shared.size() && k < (int)data.sharing.svc.size(); k++) {
        data.sharing.svc[k] = shared[k];
    }
}

// ============================================================================
// Populate latent factors on ModelData from R list
// ============================================================================
inline void populate_latent(tulpa::ModelData& data, Rcpp::List latent_spec) {
    data.has_latent = true;
    data.latent_n_factors = Rcpp::as<int>(latent_spec["n_factors"]);

    if (latent_spec.containsElementNamed("shared"))
        data.latent_shared = Rcpp::as<bool>(latent_spec["shared"]);
    if (latent_spec.containsElementNamed("scale"))
        data.latent_scale = Rcpp::as<bool>(latent_spec["scale"]);
    if (latent_spec.containsElementNamed("constraint"))
        data.latent_constraint = Rcpp::as<int>(latent_spec["constraint"]);
    if (latent_spec.containsElementNamed("sigma_prior_rate"))
        data.latent_sigma_prior_rate = Rcpp::as<double>(latent_spec["sigma_prior_rate"]);

    // Sharing
    if (latent_spec.containsElementNamed("shared_vec")) {
        Rcpp::LogicalVector shared = Rcpp::as<Rcpp::LogicalVector>(latent_spec["shared_vec"]);
        for (int k = 0; k < shared.size() && k < (int)data.sharing.latent.size(); k++) {
            data.sharing.latent[k] = shared[k];
        }
    }
}

// ============================================================================
// Build OccResponseData from R matrices
// ============================================================================
// Build OccResponseData from y matrix. Visit-level covariates added separately.
inline tulpaObs::OccResponseData build_occ_response(
    Rcpp::IntegerMatrix y_r, int n_sites, int max_visits
) {
    tulpaObs::OccResponseData occ;
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

    return occ;
}

// Add visit-level detection covariates to an existing OccResponseData
inline void add_visit_covariates(
    tulpaObs::OccResponseData& occ, Rcpp::NumericMatrix Xv
) {
    int n_sites = occ.n_sites;
    int max_visits = occ.max_visits;
    int p = Xv.ncol();
    occ.p_det_visit = p;
    occ.X_det_visit.resize(n_sites * max_visits * p);
    for (int i = 0; i < n_sites; i++) {
        for (int j = 0; j < max_visits; j++) {
            int row = i * max_visits + j;
            for (int c = 0; c < p; c++) {
                occ.X_det_visit[i * max_visits * p + j * p + c] = Xv(row, c);
            }
        }
    }
}

// ============================================================================
// Build DynOccResponseData from R vectors
// ============================================================================
inline tulpaObs::DynOccResponseData build_dyn_occ_response(
    Rcpp::IntegerVector y_flat_r,
    Rcpp::IntegerVector n_visits_r,
    Rcpp::LogicalVector any_detected_r,
    int n_sites, int n_seasons, int max_visits
) {
    tulpaObs::DynOccResponseData dyn;
    dyn.n_sites = n_sites;
    dyn.n_seasons = n_seasons;
    dyn.max_visits = max_visits;
    dyn.y = Rcpp::as<std::vector<int>>(y_flat_r);
    dyn.n_visits = Rcpp::as<std::vector<int>>(n_visits_r);
    dyn.any_detected.resize(n_sites * n_seasons);
    for (int i = 0; i < n_sites * n_seasons; i++) {
        dyn.any_detected[i] = any_detected_r[i];
    }
    return dyn;
}

// ============================================================================
// Add a process to ModelData from an R matrix
// ============================================================================
inline void add_process(tulpa::ModelData& data, Rcpp::NumericMatrix X_r, int N) {
    tulpa::ProcessData proc;
    proc.p = X_r.ncol();
    proc.X_flat.resize(N * proc.p);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < proc.p; j++)
            proc.X_flat[i * proc.p + j] = X_r(i, j);
    data.processes.push_back(proc);
}

// ============================================================================
// Run NUTS and convert results to R list
// ============================================================================
inline Rcpp::List run_nuts_and_collect(
    tulpa::ModelData& data, tulpa::ParamLayout& layout,
    int n_iter, int n_warmup, int max_treedepth, double adapt_delta,
    int seed, bool verbose
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
        nullptr,            // inv_metric_diag: default mass-adaptation
        &result
    );

    int n_samples = result.n_sample;

    Rcpp::NumericMatrix draws(n_samples, n_params);
    Rcpp::NumericVector lp(n_samples), ap(n_samples);
    Rcpp::IntegerVector div(n_samples), td(n_samples);

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

    // Posterior means
    Rcpp::NumericVector means(n_params, 0.0);
    for (int s = 0; s < n_samples; s++)
        for (int j = 0; j < n_params; j++)
            means[j] += draws(s, j) / n_samples;

    return Rcpp::List::create(
        Rcpp::Named("draws") = draws,
        Rcpp::Named("means") = means,
        Rcpp::Named("n_samples") = n_samples,
        Rcpp::Named("n_params") = n_params,
        Rcpp::Named("log_prob") = lp,
        Rcpp::Named("accept_prob") = ap,
        Rcpp::Named("divergent") = div,
        Rcpp::Named("treedepth") = td,
        Rcpp::Named("epsilon") = epsilon
    );
}

} // namespace tulpaObs

#endif // TULPAOCC_POPULATE_HELPERS_H
