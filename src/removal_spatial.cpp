// removal_spatial.cpp
// Nested-Laplace entry points for the areal-spatial removal-sampling abundance
// model: an ICAR or proper-CAR field on the abundance arm of the
// sequential-depletion removal marginal. The removal per-site marginal
// (compute_removal_site) shares the NMixSiteResult moment interface with the
// Royle N-mixture, so it reuses the family-agnostic nested-Laplace driver and
// outer-grid orchestration in nmix_count_spatial_driver.h verbatim -- only the
// per-site kernel differs. `y` is long-form per-pass removals in PASS ORDER per
// site (pass 1 first); `site_idx` maps each pass-row to its site; `X_p` is the
// per-pass detection design. The depletion offset (cumulative prior removals) is
// internal to compute_removal_site, so the abundance-arm field machinery is
// untouched. K_max must clear each site's removal TOTAL (set by the R wrapper).

#include "tobs_shape.h"
#include "nmix_count_spatial_driver.h"
#include "removal_kernel.h"
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

namespace {

// The removal per-site marginal as a callable for the generic spatial driver.
struct RemovalSiteKernel {
    tulpaObs::NMixSiteResult operator()(const int* y, const double* eta_p, int J,
                                        double eta_lambda, int K_max, double r) const {
        return tulpaObs::compute_removal_site(y, eta_p, J, eta_lambda, K_max, r);
    }
};

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_nested_laplace_removal_icar(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::NumericMatrix X_p_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector r_grid,               // NB size grid; c(Inf) for Poisson
    Rcpp::NumericVector beta_lambda_init,
    Rcpp::NumericVector beta_p_init,
    Rcpp::Nullable<Rcpp::NumericVector> z_init = R_NilValue,
    int K_max = -1,
    int max_iter = 100,
    double tol = 1e-6,
    bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = ""
) {
    return tulpaObs::run_count_nested_laplace_icar(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R,
        adj_row_ptr, adj_col_idx, n_neighbors, n_spatial,
        tau_grid, r_grid, beta_lambda_init, beta_p_init,
        tulpaObs::shape::optional_numeric(z_init.get(), "z_init"),
        K_max, max_iter, tol, verbose,
        progress, progress_every, progress_throttle, progress_file,
        RemovalSiteKernel{});
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_laplace_removal_bym2(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::NumericMatrix X_p_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector sigma_grid,
    Rcpp::NumericVector rho_grid,
    Rcpp::NumericVector r_grid,               // NB size grid; c(Inf) for Poisson
    double scale_factor,
    Rcpp::NumericVector beta_lambda_init,
    Rcpp::NumericVector beta_p_init,
    Rcpp::Nullable<Rcpp::NumericVector> v_init = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericVector> w_init = R_NilValue,
    int K_max = -1,
    int max_iter = 100,
    double tol = 1e-6,
    bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = ""
) {
    return tulpaObs::run_count_nested_laplace_bym2(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R,
        adj_row_ptr, adj_col_idx, n_neighbors, n_spatial,
        sigma_grid, rho_grid, r_grid, scale_factor,
        beta_lambda_init, beta_p_init,
        tulpaObs::shape::optional_numeric(v_init.get(), "v_init"),
        tulpaObs::shape::optional_numeric(w_init.get(), "w_init"),
        K_max, max_iter, tol, verbose,
        progress, progress_every, progress_throttle, progress_file,
        RemovalSiteKernel{});
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_laplace_removal_car_proper(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::NumericMatrix X_p_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector rho_grid,
    Rcpp::NumericVector r_grid,               // NB size grid; c(Inf) for Poisson
    Rcpp::NumericVector beta_lambda_init,
    Rcpp::NumericVector beta_p_init,
    Rcpp::Nullable<Rcpp::NumericVector> z_init = R_NilValue,
    int K_max = -1,
    int max_iter = 100,
    double tol = 1e-6,
    bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = ""
) {
    return tulpaObs::run_count_nested_laplace_car_proper(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R,
        adj_row_ptr, adj_col_idx, n_neighbors, n_spatial,
        tau_grid, rho_grid, r_grid, beta_lambda_init, beta_p_init,
        tulpaObs::shape::optional_numeric(z_init.get(), "z_init"),
        K_max, max_iter, tol, verbose,
        progress, progress_every, progress_throttle, progress_file,
        RemovalSiteKernel{});
}
