// nmix_count_spatial_driver.h
// Nested-Laplace outer-grid entry points for an areal field (ICAR / BYM2 /
// proper CAR) on the abundance arm of a count-marginal model. The per-site
// marginal is passed as a CountKernelFn, so the Royle N-mixture and removal
// families call the same compiled driver (nmix_count_spatial_driver.cpp); only
// the per-site marginal differs. `K_max < 0` takes the max(y) + 100 default.

#ifndef TULPAOBS_NMIX_COUNT_SPATIAL_DRIVER_H
#define TULPAOBS_NMIX_COUNT_SPATIAL_DRIVER_H

#include "nmix_kernel.h"   // CountKernelFn
#include <Rcpp.h>
#include <string>
#include <vector>

namespace tulpaObs {

Rcpp::List run_count_nested_laplace_icar(
    Rcpp::IntegerVector y, Rcpp::IntegerVector site_idx,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R, Rcpp::NumericMatrix X_p_R,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, int n_spatial,
    Rcpp::NumericVector tau_grid, Rcpp::NumericVector r_grid,
    Rcpp::NumericVector beta_lambda_init, Rcpp::NumericVector beta_p_init,
    const std::vector<double>& z_init,
    int K_max, int max_iter, double tol, bool verbose,
    bool progress, int progress_every, double progress_throttle,
    std::string progress_file, CountKernelFn site_fn);

Rcpp::List run_count_nested_laplace_bym2(
    Rcpp::IntegerVector y, Rcpp::IntegerVector site_idx,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R, Rcpp::NumericMatrix X_p_R,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, int n_spatial,
    Rcpp::NumericVector sigma_grid, Rcpp::NumericVector rho_grid,
    Rcpp::NumericVector r_grid, double scale_factor,
    Rcpp::NumericVector beta_lambda_init, Rcpp::NumericVector beta_p_init,
    const std::vector<double>& v_init,
    const std::vector<double>& w_init,
    int K_max, int max_iter, double tol, bool verbose,
    bool progress, int progress_every, double progress_throttle,
    std::string progress_file, CountKernelFn site_fn);

Rcpp::List run_count_nested_laplace_car_proper(
    Rcpp::IntegerVector y, Rcpp::IntegerVector site_idx,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R, Rcpp::NumericMatrix X_p_R,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, int n_spatial,
    Rcpp::NumericVector tau_grid, Rcpp::NumericVector rho_grid,
    Rcpp::NumericVector r_grid,
    Rcpp::NumericVector beta_lambda_init, Rcpp::NumericVector beta_p_init,
    const std::vector<double>& z_init,
    int K_max, int max_iter, double tol, bool verbose,
    bool progress, int progress_every, double progress_throttle,
    std::string progress_file, CountKernelFn site_fn);

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_COUNT_SPATIAL_DRIVER_H
