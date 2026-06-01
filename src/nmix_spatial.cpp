// nmix_spatial.cpp
// Nested Laplace for the spatial Royle (2004) N-mixture model with an ICAR
// (intrinsic conditional autoregressive) prior on the abundance-arm spatial
// offset.
//
// Outer grid: 1D over tau_spatial (precision of the ICAR prior).
// Inner Newton at fixed tau_k finds the joint mode of
//   x = (beta_lambda [p_lam], beta_p [p_p], z [n_spatial])
// using complete-data Fisher curvature (PSD, EM-rate convergence) with a
// fallback when the marginal observed-information matrix is non-PSD. At the
// converged mode we form the marginal observed-information Hessian (with the
// Var[N|y] rank-1 correction) and the ICAR prior contribution to compute the
// Laplace log-marginal at that grid point.
//
// Per-grid log marginal (up to a tau-independent additive constant):
//   log p(y | tau_k) ~= log L(mode_k)
//                       + log p(z_mode_k | tau_k)
//                       - 0.5 log|H(mode_k; tau_k)|
// where the prior contribution log p(z | tau) carries the (n_spatial - 1)/2 *
// log(tau) factor reflecting the ICAR rank deficiency.
//
// The Cartesian (2π)^{n_x/2} constant is omitted because it is the same at
// every grid point and only inflates absolute log_marginal values; it has no
// effect on the integration weights or hyperparameter posterior moments.
//
// Reuses the per-site marginal kernel from nmix_kernel.h and the spatial
// scatter/Hessian helpers from nmix_spatial_kernel.h.

#include "nmix_kernel.h"
#include "nmix_spatial_kernel.h"
#include "nmix_linalg.h"
#include "nmix_progress.h"
#include "nmix_spatial_grid.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using tulpaObs::nmix_kernel_sweep_spatial;
using tulpaObs::nmix_kernel_log_lik_only_spatial;
using tulpaObs::compute_eta_lambda_spatial;
using tulpaObs::nmix_assemble_obs_info_spatial;
using tulpaObs::nmix_assemble_complete_fisher_spatial;
using tulpaObs::nmix_add_car_to_spatial_block;
using tulpaObs::nmix_add_car_to_H_only;
using tulpaObs::nmix_add_diagonal_ridge;
using tulpaObs::nmix_icar_log_prior;
using tulpaObs::nmix_car_proper_log_prior;
using tulpaObs::nmix_center_z;

enum class CarPriorKind { ICAR, CAR_PROPER };

inline double nmix_car_log_prior_dispatch(
    CarPriorKind kind,
    int n_spatial, double tau, double rho, double log_det_Q_rho,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const VectorXd& z
) {
    if (kind == CarPriorKind::ICAR) {
        return nmix_icar_log_prior(
            n_spatial, tau, adj_row_ptr, adj_col_idx, n_neighbors, z);
    }
    return nmix_car_proper_log_prior(
        n_spatial, tau, rho, log_det_Q_rho,
        adj_row_ptr, adj_col_idx, n_neighbors, z);
}

// Coefficient covariance for an intrinsic (rank-deficient, sum-to-zero) field
// such as ICAR. The improper prior leaves the (intercept, field-mean) direction
// flat in the joint log-posterior -- the mode is pinned by sum-to-zero
// centering, but the unconstrained Hessian H^{-1} blows up along it and the
// intercept variance is meaningless. We add a large quadratic penalty on
// (sum of the field block)^2, the penalty-method form of the sum-to-zero
// constraint, so the beta-block of the (augmented) inverse is the constrained
// covariance (intercept variance = data precision for the global level). The
// shared `nmix_constrained_top_cov` (nmix_linalg.h) is the single source.
inline MatrixXd nmix_spatial_beta_cov(MatrixXd H, int n_x, int p_beta,
                                      int field_start, int field_len,
                                      bool constrain) {
    return tulpaObs::nmix_constrained_top_cov(std::move(H), n_x, p_beta,
                                              field_start, field_len, constrain);
}

// Per-grid-point inner solve result.
struct SpatialInnerResult {
    VectorXd beta_lambda;
    VectorXd beta_p;
    VectorXd z;
    MatrixXd cov_beta;    // (p_lam+p_p) marginal coefficient covariance at mode
    double log_lik;
    double log_marginal;
    double grad_norm;
    int n_iter;
    bool converged;
    double boundary_max;
};

// Shared inner Newton for ICAR (rho = 1.0, sum-to-zero centering, rank-deficient
// log-prior) and CAR_proper (rho < 1.0, no centering, full-rank log-prior with
// precomputed log_det_Q_rho).
SpatialInnerResult inner_newton_spatial_car(
    CarPriorKind kind,
    int p_lam, int p_p, int n_sites, int n_spatial, int n_obs,
    const Map<MatrixXd>& Xl,
    const Map<MatrixXd>& Xp,
    const Rcpp::IntegerVector& y_R,
    const std::vector<std::vector<int>>& obs_by_site,
    const std::vector<int>& map_site_to_unit,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    double tau, double rho, double log_det_Q_rho,
    double r,                          // NB size; +Inf for Poisson
    int K_max,
    int max_iter,
    double tol,
    const VectorXd& beta_lam_init,
    const VectorXd& beta_p_init,
    const VectorXd& z_init,
    bool verbose
) {
    SpatialInnerResult res;
    res.beta_lambda = beta_lam_init;
    res.beta_p      = beta_p_init;
    res.z           = z_init;
    res.converged   = false;
    res.n_iter      = 0;

    const int n_x = p_lam + p_p + n_spatial;

    VectorXd grad_eta_lam(n_sites), info_eta_lam(n_sites);
    VectorXd mean_N(n_sites), var_N(n_sites), boundary_weight(n_sites);
    VectorXd score_wt_lambda(n_sites);
    VectorXd grad_eta_p(n_obs), info_eta_p(n_obs);
    VectorXd eta_lam(n_sites);
    VectorXd eta_p_long(n_obs);

    double log_lik = R_NegInf;
    double grad_norm = R_PosInf;

    for (int iter = 0; iter < max_iter; ++iter) {
        compute_eta_lambda_spatial(Xl, res.beta_lambda, res.z,
                                   map_site_to_unit, eta_lam);
        eta_p_long.noalias() = Xp * res.beta_p;

        log_lik = nmix_kernel_sweep_spatial(
            obs_by_site, y_R, eta_lam, eta_p_long, K_max, r,
            grad_eta_lam, info_eta_lam, grad_eta_p, info_eta_p,
            mean_N, var_N, boundary_weight, score_wt_lambda
        );

        VectorXd grad = VectorXd::Zero(n_x);
        grad.segment(0, p_lam)         = Xl.transpose() * grad_eta_lam;
        grad.segment(p_lam, p_p)       = Xp.transpose() * grad_eta_p;
        // grad_z[u] = sum_{s: map[s]=u} grad_eta_lam[s]
        for (int s = 0; s < n_sites; ++s) {
            grad(p_lam + p_p + map_site_to_unit[s]) += grad_eta_lam(s);
        }

        // ICAR prior on z: subtract tau * Q z from gradient (assembled inside
        // nmix_add_icar_to_spatial_block below alongside the Hessian add).

        MatrixXd H = MatrixXd::Zero(n_x, n_x);
        nmix_assemble_obs_info_spatial(
            p_lam, p_p, n_spatial,
            Xl, Xp, eta_p_long, obs_by_site, map_site_to_unit,
            info_eta_lam, info_eta_p, var_N, score_wt_lambda, H
        );
        nmix_add_car_to_spatial_block(
            p_lam, p_p, n_spatial, tau, rho,
            adj_row_ptr, adj_col_idx, n_neighbors,
            res.z, grad, H
        );

        grad_norm = grad.norm();
        if (verbose) {
            Rcpp::Rcout << "  iter " << iter
                        << "  log_lik " << log_lik
                        << "  grad_norm " << grad_norm
                        << "  boundary_max " << boundary_weight.maxCoeff() << "\n";
        }
        if (grad_norm < tol) {
            res.converged = true;
            res.n_iter = iter + 1;
            break;
        }

        // Cholesky solve. If H_obs (with var_N correction + ICAR) is not PSD,
        // fall back to complete-data Fisher + ICAR (always PSD when tau > 0).
        // A tiny ridge accommodates the structural (intercept, constant z)
        // null direction that centering pins after each step.
        nmix_add_diagonal_ridge(H);
        VectorXd delta;
        Eigen::LLT<MatrixXd> chol(H);
        if (chol.info() == Eigen::Success) {
            delta = chol.solve(grad);
        } else {
            MatrixXd H_f = MatrixXd::Zero(n_x, n_x);
            nmix_assemble_complete_fisher_spatial(
                p_lam, p_p, n_spatial,
                Xl, Xp, obs_by_site, map_site_to_unit,
                info_eta_lam, info_eta_p, H_f
            );
            nmix_add_car_to_H_only(
                p_lam, p_p, n_spatial, tau, rho,
                adj_row_ptr, adj_col_idx, n_neighbors, H_f
            );
            nmix_add_diagonal_ridge(H_f);
            Eigen::LLT<MatrixXd> chol_f(H_f);
            if (chol_f.info() != Eigen::Success) {
                Rcpp::warning("Cholesky failure (complete-data fallback) at iter %d, tau %.4f, rho %.4f.",
                              iter, tau, rho);
                break;
            }
            delta = chol_f.solve(grad);
            if (verbose) Rcpp::Rcout << "    (Fisher fallback)\n";
        }

        VectorXd delta_lam = delta.segment(0, p_lam);
        VectorXd delta_p   = delta.segment(p_lam, p_p);
        VectorXd delta_z   = delta.segment(p_lam + p_p, n_spatial);

        // Step halving on the joint log-posterior objective.
        double step = 1.0;
        bool stepped = false;
        VectorXd beta_lam_try, beta_p_try, z_try;
        VectorXd eta_lam_try(n_sites), eta_p_try(n_obs);
        for (int h = 0; h < 12; ++h) {
            beta_lam_try = res.beta_lambda + step * delta_lam;
            beta_p_try   = res.beta_p      + step * delta_p;
            z_try        = res.z           + step * delta_z;

            compute_eta_lambda_spatial(Xl, beta_lam_try, z_try,
                                       map_site_to_unit, eta_lam_try);
            eta_p_try.noalias() = Xp * beta_p_try;
            double ll_try = nmix_kernel_log_lik_only_spatial(
                obs_by_site, y_R, eta_lam_try, eta_p_try, K_max, r
            );
            double lp_try = nmix_car_log_prior_dispatch(
                kind, n_spatial, tau, rho, log_det_Q_rho,
                adj_row_ptr, adj_col_idx, n_neighbors, z_try
            );
            double obj_try = ll_try + lp_try;
            double obj_cur = log_lik + nmix_car_log_prior_dispatch(
                kind, n_spatial, tau, rho, log_det_Q_rho,
                adj_row_ptr, adj_col_idx, n_neighbors, res.z
            );
            if (R_finite(obj_try) && obj_try >= obj_cur - 1e-10) {
                res.beta_lambda = beta_lam_try;
                res.beta_p      = beta_p_try;
                res.z           = z_try;
                // Sum-to-zero centering is needed only for the rank-deficient
                // ICAR prior; CAR_proper is full-rank and the (intercept, z)
                // direction is identified by Q(rho) itself.
                if (kind == CarPriorKind::ICAR) {
                    VectorXd x_holder(n_x);
                    x_holder.segment(0, p_lam) = res.beta_lambda;
                    x_holder.segment(p_lam, p_p) = res.beta_p;
                    x_holder.segment(p_lam + p_p, n_spatial) = res.z;
                    nmix_center_z(p_lam, p_p, n_spatial, x_holder);
                    res.z = x_holder.segment(p_lam + p_p, n_spatial);
                }
                stepped = true;
                break;
            }
            step *= 0.5;
        }
        if (!stepped) {
            if (verbose) Rcpp::Rcout << "    (step halving exhausted)\n";
            break;
        }
        res.n_iter = iter + 1;
    }

    // ----- log marginal at the converged mode -----
    compute_eta_lambda_spatial(Xl, res.beta_lambda, res.z,
                               map_site_to_unit, eta_lam);
    eta_p_long.noalias() = Xp * res.beta_p;
    double log_lik_final = nmix_kernel_sweep_spatial(
        obs_by_site, y_R, eta_lam, eta_p_long, K_max, r,
        grad_eta_lam, info_eta_lam, grad_eta_p, info_eta_p,
        mean_N, var_N, boundary_weight, score_wt_lambda
    );
    double log_prior_z_final = nmix_car_log_prior_dispatch(
        kind, n_spatial, tau, rho, log_det_Q_rho,
        adj_row_ptr, adj_col_idx, n_neighbors, res.z
    );

    MatrixXd H_final = MatrixXd::Zero(n_x, n_x);
    nmix_assemble_obs_info_spatial(
        p_lam, p_p, n_spatial,
        Xl, Xp, eta_p_long, obs_by_site, map_site_to_unit,
        info_eta_lam, info_eta_p, var_N, score_wt_lambda, H_final
    );
    // For the log|H| term, use observed info + CAR(rho). If non-PSD, fall
    // back to complete-data Fisher + CAR(rho) (will overstate curvature
    // slightly but keeps the Laplace finite).
    nmix_add_car_to_H_only(
        p_lam, p_p, n_spatial, tau, rho,
        adj_row_ptr, adj_col_idx, n_neighbors, H_final
    );
    nmix_add_diagonal_ridge(H_final);
    const int p_beta = p_lam + p_p;
    Eigen::LLT<MatrixXd> chol(H_final);
    double log_det_H;
    // ICAR is rank-deficient (sum-to-zero); CAR_proper is full-rank and
    // identifies the intercept through Q(rho), so it needs no constraint.
    const bool constrain_field = (kind == CarPriorKind::ICAR);
    if (chol.info() == Eigen::Success) {
        log_det_H = 2.0 * chol.matrixL().toDenseMatrix().diagonal()
                              .array().log().sum();
        res.cov_beta = nmix_spatial_beta_cov(H_final, n_x, p_beta,
                                             p_beta, n_spatial, constrain_field);
    } else {
        MatrixXd H_f = MatrixXd::Zero(n_x, n_x);
        nmix_assemble_complete_fisher_spatial(
            p_lam, p_p, n_spatial,
            Xl, Xp, obs_by_site, map_site_to_unit,
            info_eta_lam, info_eta_p, H_f
        );
        nmix_add_car_to_H_only(
            p_lam, p_p, n_spatial, tau, rho,
            adj_row_ptr, adj_col_idx, n_neighbors, H_f
        );
        nmix_add_diagonal_ridge(H_f);
        Eigen::LLT<MatrixXd> chol_f(H_f);
        if (chol_f.info() != Eigen::Success) {
            res.cov_beta = MatrixXd::Constant(p_beta, p_beta, R_NaN);
            res.log_marginal = R_NegInf;
            res.log_lik = log_lik_final;
            res.grad_norm = grad_norm;
            res.boundary_max = boundary_weight.maxCoeff();
            return res;
        }
        log_det_H = 2.0 * chol_f.matrixL().toDenseMatrix().diagonal()
                                .array().log().sum();
        res.cov_beta = nmix_spatial_beta_cov(H_f, n_x, p_beta,
                                             p_beta, n_spatial, constrain_field);
    }

    res.log_lik = log_lik_final;
    res.log_marginal = log_lik_final + log_prior_z_final - 0.5 * log_det_H;
    res.grad_norm = grad_norm;
    res.boundary_max = boundary_weight.maxCoeff();
    return res;
}

}  // namespace

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
    tulpaObs::NmixSpatialPrep pp = tulpaObs::prep_nmix_spatial(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R, n_spatial,
        r_grid, beta_lambda_init, beta_p_init, K_max);
    const int p_lam = pp.p_lam, p_p = pp.p_p;
    Map<MatrixXd> Xl(REAL(X_lambda_R), pp.n_sites, p_lam);
    Map<MatrixXd> Xp(REAL(X_p_R), pp.n_obs, p_p);

    VectorXd z_default(n_spatial);
    if (z_init.isNotNull()) {
        Rcpp::NumericVector zi(z_init);
        if ((int)zi.size() != n_spatial) Rcpp::stop("length(z_init) must equal n_spatial.");
        for (int s = 0; s < n_spatial; ++s) z_default(s) = zi[s];
    } else {
        z_default.setZero();
    }

    // Grid axes: NB dispersion r (outermost) x ICAR precision tau. For Poisson,
    // r_grid = c(Inf) (single pass, Poisson kernel). Fill the display grid and
    // the per-cell hyperparameters in one nesting; the driver walks them.
    const int n_tau = tau_grid.size(), n_r = r_grid.size();
    const int n_grid = n_tau * n_r;
    Rcpp::NumericMatrix theta_grid_out(n_grid, 2);
    std::vector<double> tau_k(n_grid), r_k(n_grid);
    {
        int k = 0;
        for (int ri = 0; ri < n_r; ++ri)
            for (int t = 0; t < n_tau; ++t, ++k) {
                theta_grid_out(k, 0) = tau_grid[t];
                theta_grid_out(k, 1) = r_grid[ri];
                tau_k[k] = tau_grid[t];
                r_k[k]   = r_grid[ri];
            }
    }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("tau", "r");

    // Cold restart per cell (the (lambda, p) identifiability ridge shifts with
    // the variance z absorbs across tau / r; warm-starting confounds Newton).
    auto solve = [&](int k) -> tulpaObs::NmixSpatialPoint {
        VectorXd beta_lam = pp.beta_lam_default;
        VectorXd beta_p   = pp.beta_p_default;
        VectorXd z        = z_default;
        SpatialInnerResult ir = inner_newton_spatial_car(
            CarPriorKind::ICAR,
            p_lam, p_p, pp.n_sites, n_spatial, pp.n_obs,
            Xl, Xp, y, pp.obs_by_site, pp.map_site_to_unit,
            adj_row_ptr, adj_col_idx, n_neighbors,
            tau_k[k], /*rho=*/1.0, /*log_det_Q_rho=*/0.0, r_k[k],
            pp.K_max, max_iter, tol, beta_lam, beta_p, z, verbose);
        tulpaObs::NmixSpatialPoint pt;
        pt.log_marginal = ir.log_marginal; pt.n_iter = ir.n_iter;
        pt.converged = ir.converged;       pt.grad_norm = ir.grad_norm;
        pt.log_lik = ir.log_lik;           pt.boundary_max = ir.boundary_max;
        pt.coef = VectorXd(p_lam + p_p);
        pt.coef.head(p_lam) = ir.beta_lambda;
        pt.coef.tail(p_p)   = ir.beta_p;
        pt.field = ir.z;
        pt.cov_beta = ir.cov_beta;
        return pt;
    };

    Rcpp::List out = tulpaObs::run_nmix_spatial_grid(
        n_grid, p_lam, p_p, n_spatial, /*field_len=*/n_spatial, pp.K_max,
        theta_grid_out, solve,
        progress, progress_every, progress_throttle, progress_file);
    out["tau_grid"] = tau_grid;
    out["r_grid"]   = r_grid;
    out["n_tau"]    = n_tau;
    out["n_r"]      = n_r;
    return out;
}

namespace {

// Compute log|Q(rho)| = log|D - rho * W| once per rho grid point. Returns
// -INFINITY if Q is not positive definite (callers should treat as a tail
// indicator and skip the grid point).
double log_det_Q_car_proper(
    int n_spatial, double rho,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors
) {
    MatrixXd Q = MatrixXd::Zero(n_spatial, n_spatial);
    for (int s = 0; s < n_spatial; ++s) {
        Q(s, s) = static_cast<double>(n_neighbors[s]);
        for (int kk = adj_row_ptr[s]; kk < adj_row_ptr[s + 1]; ++kk) {
            int t = adj_col_idx[kk];
            Q(s, t) = -rho;
        }
    }
    Eigen::LLT<MatrixXd> chol(Q);
    if (chol.info() != Eigen::Success) return R_NegInf;
    return 2.0 * chol.matrixL().toDenseMatrix().diagonal()
                    .array().log().sum();
}

}  // namespace

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
    tulpaObs::NmixSpatialPrep pp = tulpaObs::prep_nmix_spatial(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R, n_spatial,
        r_grid, beta_lambda_init, beta_p_init, K_max);
    const int p_lam = pp.p_lam, p_p = pp.p_p;
    Map<MatrixXd> Xl(REAL(X_lambda_R), pp.n_sites, p_lam);
    Map<MatrixXd> Xp(REAL(X_p_R), pp.n_obs, p_p);

    VectorXd z_default(n_spatial);
    if (z_init.isNotNull()) {
        Rcpp::NumericVector zi(z_init);
        if ((int)zi.size() != n_spatial) Rcpp::stop("length(z_init) must equal n_spatial.");
        for (int s = 0; s < n_spatial; ++s) z_default(s) = zi[s];
    } else {
        z_default.setZero();
    }

    const int n_tau = tau_grid.size(), n_rho = rho_grid.size(), n_r = r_grid.size();
    const int n_grid = n_tau * n_rho * n_r;

    // Precompute log|Q(rho)| once per rho grid point (the proper-CAR Q policy).
    std::vector<double> log_det_Q_rho(n_rho);
    for (int ir_rho = 0; ir_rho < n_rho; ++ir_rho)
        log_det_Q_rho[ir_rho] = log_det_Q_car_proper(
            n_spatial, rho_grid[ir_rho], adj_row_ptr, adj_col_idx, n_neighbors);

    // Grid axes: r (outermost) x rho x tau. Carry each cell's log|Q(rho)| so the
    // solver can detect the non-PD tail without re-deriving the nesting.
    Rcpp::NumericMatrix theta_grid_out(n_grid, 3);
    std::vector<double> tau_k(n_grid), rho_k(n_grid), r_k(n_grid), logdet_k(n_grid);
    {
        int k = 0;
        for (int ir_disp = 0; ir_disp < n_r; ++ir_disp)
            for (int ir_rho = 0; ir_rho < n_rho; ++ir_rho)
                for (int t = 0; t < n_tau; ++t, ++k) {
                    theta_grid_out(k, 0) = tau_grid[t];
                    theta_grid_out(k, 1) = rho_grid[ir_rho];
                    theta_grid_out(k, 2) = r_grid[ir_disp];
                    tau_k[k]    = tau_grid[t];
                    rho_k[k]    = rho_grid[ir_rho];
                    r_k[k]      = r_grid[ir_disp];
                    logdet_k[k] = log_det_Q_rho[ir_rho];
                }
    }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("tau", "rho", "r");

    auto solve = [&](int k) -> tulpaObs::NmixSpatialPoint {
        tulpaObs::NmixSpatialPoint pt;
        if (!R_finite(logdet_k[k])) { pt.skipped = true; return pt; }  // Q(rho) not PD
        VectorXd beta_lam = pp.beta_lam_default;
        VectorXd beta_p   = pp.beta_p_default;
        VectorXd z        = z_default;
        SpatialInnerResult ir = inner_newton_spatial_car(
            CarPriorKind::CAR_PROPER,
            p_lam, p_p, pp.n_sites, n_spatial, pp.n_obs,
            Xl, Xp, y, pp.obs_by_site, pp.map_site_to_unit,
            adj_row_ptr, adj_col_idx, n_neighbors,
            tau_k[k], rho_k[k], logdet_k[k], r_k[k],
            pp.K_max, max_iter, tol, beta_lam, beta_p, z, verbose);
        pt.log_marginal = ir.log_marginal; pt.n_iter = ir.n_iter;
        pt.converged = ir.converged;       pt.grad_norm = ir.grad_norm;
        pt.log_lik = ir.log_lik;           pt.boundary_max = ir.boundary_max;
        pt.coef = VectorXd(p_lam + p_p);
        pt.coef.head(p_lam) = ir.beta_lambda;
        pt.coef.tail(p_p)   = ir.beta_p;
        pt.field = ir.z;
        pt.cov_beta = ir.cov_beta;
        return pt;
    };

    Rcpp::List out = tulpaObs::run_nmix_spatial_grid(
        n_grid, p_lam, p_p, n_spatial, /*field_len=*/n_spatial, pp.K_max,
        theta_grid_out, solve,
        progress, progress_every, progress_throttle, progress_file);
    out["tau_grid"]      = tau_grid;
    out["rho_grid"]      = rho_grid;
    out["r_grid"]        = r_grid;
    out["log_det_Q_rho"] = Rcpp::wrap(log_det_Q_rho);
    out["n_tau"]         = n_tau;
    out["n_rho"]         = n_rho;
    out["n_r"]           = n_r;
    return out;
}
