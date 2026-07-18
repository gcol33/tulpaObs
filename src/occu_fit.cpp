// occu_fit.cpp
// Unified Rcpp entry point for all tulpaObs model types.
// Selects likelihood by model_type, builds ModelData compositionally
// with orthogonal spatial/temporal/RE/SVC/latent components.

#include <Rcpp.h>
#include <vector>
#include <string>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>

#include "occ_data.h"
#include "dyn_occ_data.h"
#include "integrated_occ_data.h"
#include "occ_likelihood.h"
#include "dyn_occ_likelihood.h"
#include "integrated_occ_likelihood.h"
#include "populate_helpers.h"

using namespace Rcpp;

// ============================================================================
// Unified entry point: all occupancy model types
// ============================================================================

// [[Rcpp::export]]
Rcpp::List cpp_occu_fit(Rcpp::List spec_r) {
    // ---- Unpack model type and sampler settings ----
    std::string model_type = Rcpp::as<std::string>(spec_r["model_type"]);
    double sigma_beta   = Rcpp::as<double>(spec_r["sigma_beta"]);
    int n_iter          = Rcpp::as<int>(spec_r["n_iter"]);
    int n_warmup        = Rcpp::as<int>(spec_r["n_warmup"]);
    int max_treedepth   = Rcpp::as<int>(spec_r["max_treedepth"]);
    double adapt_delta   = Rcpp::as<double>(spec_r["adapt_delta"]);
    int seed            = Rcpp::as<int>(spec_r["seed"]);
    bool verbose        = Rcpp::as<bool>(spec_r["verbose"]);

    // ---- Determine n_processes and build response data ----
    int n_processes;
    int N;
    int p_det_visit = 0;

    // Response data — one of these will be populated
    tulpaObs::OccResponseData occ_response;
    tulpaObs::DynOccResponseData dyn_response;
    tulpaObs::IntegratedOccResponseData int_response;
    void* response_ptr = nullptr;

    // LikelihoodSpec
    tulpa::LikelihoodSpec spec;
    spec.n_extra_params = 0;

    if (model_type == "single") {
        n_processes = 2;
        IntegerMatrix y_r = Rcpp::as<IntegerMatrix>(spec_r["y"]);
        N = y_r.nrow();
        int max_visits = y_r.ncol();

        occ_response = tulpaObs::build_occ_response(y_r, N, max_visits);

        // Add visit-level covariates if present
        if (spec_r.containsElementNamed("X_det_visit")) {
            SEXP xv = spec_r["X_det_visit"];
            if (!Rf_isNull(xv)) {
                NumericMatrix Xv = Rcpp::as<NumericMatrix>(xv);
                tulpaObs::add_visit_covariates(occ_response, Xv);
                p_det_visit = Xv.ncol();
            }
        }

        spec.name = "occupancy";
        spec.ll_double = tulpaObs::occ_log_likelihood<double>;
        spec.ll_arena  = tulpaObs::occ_log_likelihood<tulpa::arena::Var>;
        spec.ll_fwd    = tulpaObs::occ_log_likelihood<fwd::Dual>;
        spec.residual_fn = tulpaObs::occ_residual;
        spec.n_extra_params = p_det_visit;
        response_ptr = &occ_response;

    } else if (model_type == "dynamic") {
        n_processes = 4;
        IntegerVector y_flat_r    = Rcpp::as<IntegerVector>(spec_r["y_flat"]);
        IntegerVector n_visits_r  = Rcpp::as<IntegerVector>(spec_r["n_visits"]);
        LogicalVector any_det_r   = Rcpp::as<LogicalVector>(spec_r["any_detected"]);
        int n_sites   = Rcpp::as<int>(spec_r["n_sites"]);
        int n_seasons = Rcpp::as<int>(spec_r["n_seasons"]);
        int max_visits = Rcpp::as<int>(spec_r["max_visits"]);
        N = n_sites;

        dyn_response = tulpaObs::build_dyn_occ_response(
            y_flat_r, n_visits_r, any_det_r, n_sites, n_seasons, max_visits);

        spec.name = "dynamic_occupancy";
        spec.ll_double = tulpaObs::dyn_occ_log_likelihood<double>;
        spec.ll_arena  = tulpaObs::dyn_occ_log_likelihood<tulpa::arena::Var>;
        spec.ll_fwd    = tulpaObs::dyn_occ_log_likelihood<fwd::Dual>;
        spec.n_extra_params = 0;
        response_ptr = &dyn_response;

    } else if (model_type == "integrated") {
        int n_sources = Rcpp::as<int>(spec_r["n_sources"]);
        n_processes = 1 + n_sources;  // 1 occupancy + S detection processes
        int n_sites_global = Rcpp::as<int>(spec_r["n_sites"]);
        N = n_sites_global;

        int_response.n_sites = n_sites_global;
        int_response.n_sources = n_sources;
        int_response.y.resize(n_sources);
        int_response.max_visits.resize(n_sources);
        int_response.n_sites_per.resize(n_sources);
        int_response.n_visits.resize(n_sources);
        int_response.any_detected.resize(n_sources);
        int_response.site_map.resize(n_sources);

        Rcpp::List y_list = Rcpp::as<Rcpp::List>(spec_r["y_sources"]);
        Rcpp::List site_map_list = Rcpp::as<Rcpp::List>(spec_r["site_maps"]);
        for (int s = 0; s < n_sources; s++) {
            Rcpp::IntegerMatrix ys = Rcpp::as<Rcpp::IntegerMatrix>(y_list[s]);
            int ns = ys.nrow();
            int mv = ys.ncol();
            int_response.n_sites_per[s] = ns;
            int_response.max_visits[s] = mv;

            int_response.y[s].resize(ns * mv);
            int_response.n_visits[s].resize(ns);
            int_response.any_detected[s].resize(ns, false);
            int_response.site_map[s] = Rcpp::as<std::vector<int>>(site_map_list[s]);

            for (int i = 0; i < ns; i++) {
                int nv = 0;
                for (int j = 0; j < mv; j++) {
                    int val = ys(i, j);
                    int_response.y[s][i * mv + j] = val;
                    if (val >= 0) {
                        nv++;
                        if (val == 1) int_response.any_detected[s][i] = true;
                    }
                }
                int_response.n_visits[s][i] = nv;
            }
        }

        spec.name = "integrated_occupancy";
        spec.ll_double = tulpaObs::integrated_occ_log_likelihood<double>;
        spec.ll_arena  = tulpaObs::integrated_occ_log_likelihood<tulpa::arena::Var>;
        spec.ll_fwd    = tulpaObs::integrated_occ_log_likelihood<fwd::Dual>;
        spec.n_extra_params = 0;
        response_ptr = &int_response;

    } else {
        Rcpp::stop("Unknown model_type: %s", model_type.c_str());
    }

    // ---- Build ModelData ----
    tulpa::ModelData data;
    data.N = N;
    data.n_processes = n_processes;
    data.sigma_beta = sigma_beta;
    data.model_response_data = response_ptr;
    data.likelihood_spec = &spec;

    // ---- Add processes from design matrices ----
    List X_list = Rcpp::as<List>(spec_r["X_processes"]);
    for (int k = 0; k < n_processes; k++) {
        NumericMatrix X_r = Rcpp::as<NumericMatrix>(X_list[k]);
        tulpaObs::add_process(data, X_r, N);
    }

    data.sharing.init(n_processes);

    // ---- ZI/OI not used ----
    data.zi_type = tulpa::ZIType::NONE;
    data.p_zi = 0;
    data.p_oi = 0;
    data.zi_prior_sd = 1.0;
    data.oi_prior_sd = 1.0;

    // ---- Spatial (optional) ----
    if (spec_r.containsElementNamed("spatial_params")) {
        SEXP sp_sexp = spec_r["spatial_params"];
        if (!Rf_isNull(sp_sexp)) {
            Rcpp::List sp = Rcpp::as<Rcpp::List>(sp_sexp);
            std::string sp_type = Rcpp::as<std::string>(sp["type"]);
            if (sp_type != "none") {
                tulpaObs::populate_spatial(data, sp, N);
            }
        }
    }

    // ---- Random effects (optional) ----
    // Single-term species RE path (jsdm species random intercept on the one arm)
    if (spec_r.containsElementNamed("re_group")) {
        SEXP re_sexp = spec_r["re_group"];
        if (!Rf_isNull(re_sexp)) {
            data.re_group = Rcpp::as<std::vector<int>>(spec_r["re_group"]);
            data.n_re_groups = Rcpp::as<int>(spec_r["n_re_groups"]);
            data.n_re_terms = 0;
            data.total_re_groups = data.n_re_groups;
            data.has_re_slopes = false;
            data.has_re_correlated_slopes = false;
            data.total_re_params = data.n_re_groups;
            data.total_sigma_params = 1;
            data.total_chol_params = 0;
            data.re_parameterization = 1;

            if (spec_r.containsElementNamed("sigma_re_scale"))
                data.sigma_re_scale = Rcpp::as<double>(spec_r["sigma_re_scale"]);

            bool re_occ = true, re_det = true;
            if (spec_r.containsElementNamed("re_shared_occ"))
                re_occ = Rcpp::as<bool>(spec_r["re_shared_occ"]);
            if (spec_r.containsElementNamed("re_shared_det"))
                re_det = Rcpp::as<bool>(spec_r["re_shared_det"]);
            data.sharing.re[0] = re_occ;
            data.sharing.re[1] = re_det;
        }
    }
    // Multi-term RE path (user-specified via occu_re())
    if (spec_r.containsElementNamed("re_spec")) {
        SEXP re_sexp = spec_r["re_spec"];
        if (!Rf_isNull(re_sexp)) {
            tulpaObs::populate_re(data, Rcpp::as<Rcpp::List>(re_sexp));
        }
    }

    // ---- Temporal (optional) ----
    if (spec_r.containsElementNamed("temporal_spec")) {
        SEXP temp_sexp = spec_r["temporal_spec"];
        if (!Rf_isNull(temp_sexp)) {
            tulpaObs::populate_temporal(data, Rcpp::as<Rcpp::List>(temp_sexp));
        }
    }

    // ---- Spatially-varying coefficients (optional) ----
    if (spec_r.containsElementNamed("svc_spec")) {
        SEXP svc_sexp = spec_r["svc_spec"];
        if (!Rf_isNull(svc_sexp)) {
            tulpaObs::populate_svc(data, Rcpp::as<Rcpp::List>(svc_sexp));
        }
    }

    // ---- Latent factors (optional) ----
    if (spec_r.containsElementNamed("latent_spec")) {
        SEXP lat_sexp = spec_r["latent_spec"];
        if (!Rf_isNull(lat_sexp)) {
            tulpaObs::populate_latent(data, Rcpp::as<Rcpp::List>(lat_sexp));
        }
    }

    // ---- Compute layout and run NUTS ----
    tulpa::ParamLayout layout = tulpa::compute_layout(data);

    Rcpp::List result = tulpaObs::run_nuts_and_collect(
        data, layout, n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose);

    // ---- Add column names for fixed effects ----
    int n_params = layout.total_params;
    CharacterVector col_names(n_params);
    int idx = 0;

    if (spec_r.containsElementNamed("process_names")) {
        List pnames = Rcpp::as<List>(spec_r["process_names"]);
        for (int k = 0; k < n_processes; k++) {
            CharacterVector names_k = Rcpp::as<CharacterVector>(pnames[k]);
            for (int j = 0; j < names_k.size(); j++) {
                if (idx < n_params) col_names[idx++] = names_k[j];
            }
        }
    }

    if (spec_r.containsElementNamed("extra_param_names")) {
        CharacterVector extra = Rcpp::as<CharacterVector>(spec_r["extra_param_names"]);
        for (int j = 0; j < extra.size(); j++) {
            if (idx < n_params) col_names[idx++] = extra[j];
        }
    }

    for (; idx < n_params; idx++)
        col_names[idx] = "param[" + std::to_string(idx + 1) + "]";

    // ---- Name the SVC block (gcol33/tulpaObs#118) ----
    // The engine exports the SVC offsets on ParamLayout, so the block is named
    // here instead of falling through to "param[k]" above. The offsets are
    // absolute positions in the parameter vector, so they are written directly
    // rather than through the running cursor, and after the fallback loop so
    // they are not overwritten by it.
    if (layout.has_svc && data.svc_data.n_svc > 0) {
        const int n_svc = data.svc_data.n_svc;
        const int n_loc = data.svc_data.n_obs;
        for (int j = 0; j < n_svc; j++) {
            const int a = layout.log_sigma2_svc_start + j;
            if (a >= 0 && a < n_params)
                col_names[a] = "log_sigma2_svc[" + std::to_string(j + 1) + "]";
            const int b = layout.log_phi_svc_start + j;
            if (b >= 0 && b < n_params)
                col_names[b] = "log_phi_svc[" + std::to_string(j + 1) + "]";
        }
        // Per-location weights, w_flat[j * n_obs + i] (engine stride: j indexes
        // the SVC term, i the location). Only name per-location when the block
        // really is one weight per location -- the HSGP flavour stores
        // n_svc * m_total basis coefficients instead, and those are not a
        // per-location surface.
        const int w_len = layout.svc_w_end - layout.svc_w_start;
        if (w_len == n_svc * n_loc) {
            for (int j = 0; j < n_svc; j++) {
                for (int i = 0; i < n_loc; i++) {
                    const int c = layout.svc_w_start + j * n_loc + i;
                    if (c >= 0 && c < n_params)
                        col_names[c] = "svc_w[" + std::to_string(i + 1) + "," +
                                       std::to_string(j + 1) + "]";
                }
            }
        }
    }

    NumericMatrix draws = Rcpp::as<NumericMatrix>(result["draws"]);
    Rcpp::colnames(draws) = col_names;

    NumericVector means = Rcpp::as<NumericVector>(result["means"]);
    means.names() = col_names;

    result["draws"] = draws;
    result["means"] = means;
    result["col_names"] = col_names;

    // The SVC block's layout, so R can slice the surface off the draws without
    // re-deriving offsets or parsing column names. 1-based for R.
    if (layout.has_svc && data.svc_data.n_svc > 0) {
        result["svc_layout"] = Rcpp::List::create(
            Rcpp::Named("n_svc")            = data.svc_data.n_svc,
            Rcpp::Named("n_obs")            = data.svc_data.n_obs,
            Rcpp::Named("w_start")          = layout.svc_w_start + 1,
            Rcpp::Named("w_end")            = layout.svc_w_end,
            Rcpp::Named("log_sigma2_start") = layout.log_sigma2_svc_start + 1,
            Rcpp::Named("log_phi_start")    = layout.log_phi_svc_start + 1);
    }

    return result;
}
