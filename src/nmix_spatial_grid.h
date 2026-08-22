#ifndef TULPAOBS_NMIX_SPATIAL_GRID_H
#define TULPAOBS_NMIX_SPATIAL_GRID_H

// Shared scaffolding for the single-species spatial nested-Laplace N-mixture
// entries (cpp_nested_laplace_nmix_{icar,car_proper,bym2}). The per-field-kind
// entries differ only in their spatial Q policy: the hyperparameter axes they
// sweep, the spatial precision builder, and the latent field layout. Everything
// around that -- input validation, the site/obs bucketing, the outer-grid walk
// with its cold restart, the per-cell result -> output-row assembly, progress
// plumbing, and the common return block -- is identical and lives here once.
//
//   prep_nmix_spatial()    validates the shared inputs and builds the obs-by-site
//                          bucketing, the site->unit map, the coef inits, and the
//                          K_max default.
//   run_nmix_spatial_grid()walks a caller-filled theta_grid, calls a per-cell
//                          solver (which owns the field-kind inner Newton), and
//                          assembles the common output list. Each entry appends
//                          its own field-hyperparameter axes to the result.

#include "nmix_progress.h"
#include <RcppEigen.h>
#include <functional>
#include <vector>

namespace tulpaObs {

// Validated, shared inputs common to every single-species spatial entry.
struct NmixSpatialPrep {
    int n_sites = 0, p_lam = 0, n_obs = 0, p_p = 0, K_max = 0;
    std::vector<std::vector<int>> obs_by_site;
    std::vector<int> map_site_to_unit;
    Eigen::VectorXd beta_lam_default, beta_p_default;
};

inline NmixSpatialPrep prep_nmix_spatial(
    const Rcpp::IntegerVector& y, const Rcpp::IntegerVector& site_idx,
    const Rcpp::IntegerVector& map_site_to_unit_R,
    const Rcpp::NumericMatrix& X_lambda_R, const Rcpp::NumericMatrix& X_p_R,
    int n_spatial, const Rcpp::NumericVector& r_grid,
    const Rcpp::NumericVector& beta_lambda_init,
    const Rcpp::NumericVector& beta_p_init,
    int K_max
) {
    NmixSpatialPrep pp;
    pp.n_sites = X_lambda_R.nrow();
    pp.p_lam   = X_lambda_R.ncol();
    pp.n_obs   = X_p_R.nrow();
    pp.p_p     = X_p_R.ncol();
    if ((int)y.size() != pp.n_obs) Rcpp::stop("length(y) must equal nrow(X_p).");
    if ((int)site_idx.size() != pp.n_obs) Rcpp::stop("length(site_idx) must equal nrow(X_p).");
    if ((int)map_site_to_unit_R.size() != pp.n_sites)
        Rcpp::stop("length(map_site_to_unit) must equal nrow(X_lambda).");
    if ((int)beta_lambda_init.size() != pp.p_lam) Rcpp::stop("beta_lambda_init length mismatch.");
    if ((int)beta_p_init.size() != pp.p_p) Rcpp::stop("beta_p_init length mismatch.");
    if (r_grid.size() < 1) Rcpp::stop("r_grid must have length >= 1.");
    if (K_max < 0) {
        int ymax = 0;
        for (int o = 0; o < pp.n_obs; ++o) if (y[o] > ymax) ymax = y[o];
        K_max = ymax + 100;
    }
    pp.K_max = K_max;

    pp.obs_by_site.assign(pp.n_sites, std::vector<int>());
    for (int o = 0; o < pp.n_obs; ++o) {
        int s = site_idx[o] - 1;
        if (s < 0 || s >= pp.n_sites) Rcpp::stop("site_idx out of range at obs %d.", o + 1);
        pp.obs_by_site[s].push_back(o);
    }
    pp.map_site_to_unit.resize(pp.n_sites);
    for (int s = 0; s < pp.n_sites; ++s) {
        int u = map_site_to_unit_R[s] - 1;
        if (u < 0 || u >= n_spatial)
            Rcpp::stop("map_site_to_unit[%d] = %d out of range [1, %d].",
                       s + 1, map_site_to_unit_R[s], n_spatial);
        pp.map_site_to_unit[s] = u;
    }
    pp.beta_lam_default = Eigen::Map<Eigen::VectorXd>(REAL(beta_lambda_init), pp.p_lam);
    pp.beta_p_default   = Eigen::Map<Eigen::VectorXd>(REAL(beta_p_init), pp.p_p);
    return pp;
}

// Normalized per-cell inner-Newton result. `skipped` marks a cell the entry
// declined before the inner solve (CAR non-PD Q(rho), BYM2 out-of-range
// (sigma, rho)): the driver writes the same loud placeholder the hand-written
// loops did and leaves the mode row / covariance untouched.
struct NmixSpatialPoint {
    bool skipped = false;
    double log_marginal = 0.0, grad_norm = 0.0, log_lik = 0.0, boundary_max = 0.0;
    int n_iter = 0;
    bool converged = false;
    Eigen::VectorXd coef;    // [beta_lambda (p_lam); beta_p (p_p)]
    Eigen::VectorXd field;   // z (n_spatial) or [v; w] (2 * n_spatial)
    Eigen::MatrixXd cov_beta;
};

// Walk a caller-filled theta_grid_out, solve each cell, and assemble the common
// return block. `solve_point(k)` owns the cold restart and the field-kind inner
// Newton; the driver owns the buffers, the mode/cov assembly, and progress.
inline Rcpp::List run_nmix_spatial_grid(
    int n_grid, int p_lam, int p_p, int n_spatial, int field_len, int K_max,
    const Rcpp::NumericMatrix& theta_grid_out,
    const std::function<NmixSpatialPoint(int)>& solve_point,
    bool progress, int progress_every, double progress_throttle,
    const std::string& progress_file
) {
    const int n_x = p_lam + p_p + field_len;
    Rcpp::NumericVector log_marginals(n_grid), grad_norms(n_grid),
                        log_liks(n_grid), boundary_maxes(n_grid);
    Rcpp::IntegerVector n_iters(n_grid);
    Rcpp::LogicalVector convergeds(n_grid);
    Rcpp::NumericMatrix modes(n_grid, n_x);
    Rcpp::List cov_blocks(n_grid);   // per-grid marginal coef covariance

    // outer-grid progress
    auto gp = tulpaObs::make_grid_progress("nmix-spatial", n_grid, progress,
                                           progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        NmixSpatialPoint pt = solve_point(k);
        if (pt.skipped) {
            log_marginals[k]  = R_NegInf;
            n_iters[k]        = 0;
            convergeds[k]     = false;
            grad_norms[k]     = R_PosInf;
            log_liks[k]       = R_NegInf;
            boundary_maxes[k] = 0.0;
            // modes row stays 0; cov_blocks[k] stays NULL (matches the
            // hand-written loops' `continue` past the non-PD / invalid cell).
        } else {
            log_marginals[k]  = pt.log_marginal;
            n_iters[k]        = pt.n_iter;
            convergeds[k]     = pt.converged;
            grad_norms[k]     = pt.grad_norm;
            log_liks[k]       = pt.log_lik;
            boundary_maxes[k] = pt.boundary_max;
            for (int j = 0; j < p_lam + p_p; ++j) modes(k, j) = pt.coef(j);
            for (int j = 0; j < field_len; ++j)   modes(k, p_lam + p_p + j) = pt.field(j);
            cov_blocks[k] = Rcpp::wrap(pt.cov_beta);
        }
        if (gp) gp->tick();
    }
    if (gp) gp->finish();

    return Rcpp::List::create(
        Rcpp::Named("theta_grid")   = theta_grid_out,
        Rcpp::Named("log_marginal") = log_marginals,
        Rcpp::Named("modes")        = modes,
        Rcpp::Named("cov_blocks")   = cov_blocks,
        Rcpp::Named("n_iter")       = n_iters,
        Rcpp::Named("converged")    = convergeds,
        Rcpp::Named("grad_norm")    = grad_norms,
        Rcpp::Named("log_lik")      = log_liks,
        Rcpp::Named("boundary_max") = boundary_maxes,
        Rcpp::Named("p_lambda")     = p_lam,
        Rcpp::Named("p_p")          = p_p,
        Rcpp::Named("n_spatial")    = n_spatial,
        Rcpp::Named("n_grid")       = n_grid,
        Rcpp::Named("K_max")        = K_max);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_SPATIAL_GRID_H
