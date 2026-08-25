// nmix_spatial_bym2.cpp
// Nested-Laplace entry point for the spatial Royle (2004) N-mixture model with a
// BYM2 (Riebler et al. 2016) field on the abundance arm. The family-agnostic
// inner Newton + outer-grid orchestration live in nmix_count_spatial_driver.h
// (shared with removal, #51); this is a thin wrapper instantiating it with the
// Royle per-site kernel.

#include "tobs_shape.h"
#include "nmix_count_spatial_driver.h"
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Rcpp::List cpp_nested_laplace_nmix_bym2(
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
        tulpaObs::NmixSiteKernel{});
}
