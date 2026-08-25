// nmix_spatial.cpp
// Nested-Laplace entry points for the spatial Royle (2004) N-mixture model with
// an ICAR or proper-CAR field on the abundance arm. The family-agnostic inner
// Newton / Laplace machinery AND the outer-grid orchestration live in
// nmix_count_spatial_driver.h (shared with removal / distance, #51); these two
// [[Rcpp::export]] functions are thin wrappers that instantiate the shared
// orchestration with the Royle per-site kernel.

#include "tobs_shape.h"
#include "nmix_count_spatial_driver.h"
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List cpp_nested_laplace_nmix_icar(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::IntegerVector map_site_to_unit_R,   // 1-based site -> unit map
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
        tulpaObs::NmixSiteKernel{});
}

// [[Rcpp::export]]
Rcpp::List cpp_nested_laplace_nmix_car_proper(
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
        tulpaObs::NmixSiteKernel{});
}
