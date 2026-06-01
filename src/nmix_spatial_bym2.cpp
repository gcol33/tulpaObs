// nmix_spatial_bym2.cpp
// Nested Laplace for the spatial Royle (2004) N-mixture model with a BYM2
// (Riebler et al. 2016) prior on the abundance-arm spatial offset.
//
// Outer grid: 2D over (sigma, rho). sigma is the joint marginal standard
// deviation of the spatial offset; rho is the spatial fraction of
// variance. Inner Newton at fixed (sigma_k, rho_k) finds the joint mode of
//   x = (beta_lambda [p_lam], beta_p [p_p], v [n_spatial], w [n_spatial])
// using complete-data Fisher curvature with a fallback when the marginal
// observed-information matrix is non-PSD.
//
// Phi (the actual offset) is sigma * (sqrt(rho/scale)*v + sqrt(1-rho)*w);
// the priors on (v, w) are independent of (sigma, rho), so the prior
// contribution to log p(y | sigma, rho) is constant across grid points
// and absorbs into the same tau-independent additive constant as the
// (2 pi)^{n_x/2} Cartesian factor. Only the inner mode and log|H| change
// with (sigma, rho), and so does the data likelihood through the
// coefficients (a, b) of (v, w) in eta_lambda.

#include "nmix_kernel.h"
#include "nmix_spatial_kernel.h"      // nmix_kernel_sweep_spatial, log_lik_only_spatial
#include "nmix_spatial_kernel_bym2.h"
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
using tulpaObs::compute_eta_lambda_bym2;
using tulpaObs::nmix_assemble_obs_info_bym2;
using tulpaObs::nmix_assemble_complete_fisher_bym2;
using tulpaObs::nmix_add_bym2_prior_to_grad_and_H;
using tulpaObs::nmix_add_bym2_prior_to_H_only;
using tulpaObs::nmix_bym2_log_prior;
using tulpaObs::nmix_center_v_bym2;

inline void add_diagonal_ridge_bym2(MatrixXd& H, double rel_ridge = 1e-10) {
    const int n = static_cast<int>(H.rows());
    if (n == 0) return;
    double mean_diag = 0.0;
    for (int i = 0; i < n; ++i) mean_diag += H(i, i);
    mean_diag /= n;
    double ridge = std::max(rel_ridge * mean_diag, 1e-12);
    for (int i = 0; i < n; ++i) H(i, i) += ridge;
}

// Constrained coefficient covariance for BYM2. The structured component v is
// rank-deficient (sum-to-zero), so the (intercept, v-mean) direction is flat in
// the joint posterior and the unconstrained intercept variance is meaningless.
// Pin sum(v)=0 with a large quadratic penalty (the penalty-method form of the
// constraint); the iid component w is proper and is left alone. The shared
// `nmix_constrained_top_cov` (nmix_linalg.h) is the single source.
inline MatrixXd nmix_beta_cov_bym2(MatrixXd H, int n_x, int p_beta,
                                   int v_start, int n_spatial) {
    return tulpaObs::nmix_constrained_top_cov(std::move(H), n_x, p_beta,
                                              v_start, n_spatial, /*constrain=*/true);
}

struct BYM2InnerResult {
    VectorXd beta_lambda;
    VectorXd beta_p;
    VectorXd v;
    VectorXd w;
    MatrixXd cov_beta;    // (p_lam+p_p) marginal coefficient covariance at mode
    double log_lik;
    double log_marginal;
    double grad_norm;
    int n_iter;
    bool converged;
    double boundary_max;
};

BYM2InnerResult inner_newton_bym2(
    int p_lam, int p_p, int n_sites, int n_spatial, int n_obs,
    double a, double b,
    const Map<MatrixXd>& Xl,
    const Map<MatrixXd>& Xp,
    const Rcpp::IntegerVector& y_R,
    const std::vector<std::vector<int>>& obs_by_site,
    const std::vector<int>& map_site_to_unit,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    double r,                          // NB size; +Inf for Poisson
    int K_max,
    int max_iter,
    double tol,
    const VectorXd& beta_lam_init,
    const VectorXd& beta_p_init,
    const VectorXd& v_init,
    const VectorXd& w_init,
    bool verbose
) {
    BYM2InnerResult res;
    res.beta_lambda = beta_lam_init;
    res.beta_p      = beta_p_init;
    res.v           = v_init;
    res.w           = w_init;
    res.converged   = false;
    res.n_iter      = 0;

    const int n_x = p_lam + p_p + 2 * n_spatial;

    VectorXd grad_eta_lam(n_sites), info_eta_lam(n_sites);
    VectorXd mean_N(n_sites), var_N(n_sites), boundary_weight(n_sites);
    VectorXd score_wt_lambda(n_sites);
    VectorXd grad_eta_p(n_obs), info_eta_p(n_obs);
    VectorXd eta_lam(n_sites);
    VectorXd eta_p_long(n_obs);

    double log_lik = R_NegInf;
    double grad_norm = R_PosInf;

    for (int iter = 0; iter < max_iter; ++iter) {
        compute_eta_lambda_bym2(Xl, res.beta_lambda, res.v, res.w, a, b,
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
        // grad_v[u] = a * sum_{s: map[s]=u} grad_eta_lam[s]
        // grad_w[u] = b * sum_{s: map[s]=u} grad_eta_lam[s]
        const int v_start = p_lam + p_p;
        const int w_start = v_start + n_spatial;
        for (int s = 0; s < n_sites; ++s) {
            const int u = map_site_to_unit[s];
            grad(v_start + u) += a * grad_eta_lam(s);
            grad(w_start + u) += b * grad_eta_lam(s);
        }

        MatrixXd H = MatrixXd::Zero(n_x, n_x);
        nmix_assemble_obs_info_bym2(
            p_lam, p_p, n_spatial, a, b,
            Xl, Xp, eta_p_long, obs_by_site, map_site_to_unit,
            info_eta_lam, info_eta_p, var_N, score_wt_lambda, H
        );
        nmix_add_bym2_prior_to_grad_and_H(
            p_lam, p_p, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors,
            res.v, res.w, grad, H
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

        // Cholesky solve. Fall back to complete-data Fisher if non-PSD.
        add_diagonal_ridge_bym2(H);
        VectorXd delta;
        Eigen::LLT<MatrixXd> chol(H);
        if (chol.info() == Eigen::Success) {
            delta = chol.solve(grad);
        } else {
            MatrixXd H_f = MatrixXd::Zero(n_x, n_x);
            nmix_assemble_complete_fisher_bym2(
                p_lam, p_p, n_spatial, a, b,
                Xl, Xp, obs_by_site, map_site_to_unit,
                info_eta_lam, info_eta_p, H_f
            );
            // Re-add the prior to the Fisher fallback (grad already updated
            // above). nmix_add_bym2_prior_to_H_only adds Q + I to the v, w
            // blocks; we only need the Hessian here.
            nmix_add_bym2_prior_to_H_only(
                p_lam, p_p, n_spatial,
                adj_row_ptr, adj_col_idx, n_neighbors, H_f
            );
            add_diagonal_ridge_bym2(H_f);
            Eigen::LLT<MatrixXd> chol_f(H_f);
            if (chol_f.info() != Eigen::Success) {
                Rcpp::warning("Cholesky failure (complete-data fallback) at iter %d, a %.4f, b %.4f.",
                              iter, a, b);
                break;
            }
            delta = chol_f.solve(grad);
            if (verbose) Rcpp::Rcout << "    (Fisher fallback)\n";
        }

        VectorXd delta_lam = delta.segment(0, p_lam);
        VectorXd delta_p   = delta.segment(p_lam, p_p);
        VectorXd delta_v   = delta.segment(v_start, n_spatial);
        VectorXd delta_w   = delta.segment(w_start, n_spatial);

        // Step halving on the joint log-posterior.
        double step = 1.0;
        bool stepped = false;
        VectorXd beta_lam_try, beta_p_try, v_try, w_try;
        VectorXd eta_lam_try(n_sites), eta_p_try(n_obs);
        for (int h = 0; h < 12; ++h) {
            beta_lam_try = res.beta_lambda + step * delta_lam;
            beta_p_try   = res.beta_p      + step * delta_p;
            v_try        = res.v           + step * delta_v;
            w_try        = res.w           + step * delta_w;

            compute_eta_lambda_bym2(Xl, beta_lam_try, v_try, w_try, a, b,
                                    map_site_to_unit, eta_lam_try);
            eta_p_try.noalias() = Xp * beta_p_try;
            double ll_try = nmix_kernel_log_lik_only_spatial(
                obs_by_site, y_R, eta_lam_try, eta_p_try, K_max, r
            );
            double lp_try = nmix_bym2_log_prior(
                n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
                v_try, w_try
            );
            double obj_try = ll_try + lp_try;
            double obj_cur = log_lik + nmix_bym2_log_prior(
                n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
                res.v, res.w
            );
            if (R_finite(obj_try) && obj_try >= obj_cur - 1e-10) {
                res.beta_lambda = beta_lam_try;
                res.beta_p      = beta_p_try;
                res.v           = v_try;
                res.w           = w_try;

                // Sum-to-zero on v (ICAR identifiability). w is identified
                // by its iid prior and is not centred.
                VectorXd x_holder(n_x);
                x_holder.segment(0, p_lam) = res.beta_lambda;
                x_holder.segment(p_lam, p_p) = res.beta_p;
                x_holder.segment(v_start, n_spatial) = res.v;
                x_holder.segment(w_start, n_spatial) = res.w;
                nmix_center_v_bym2(p_lam, p_p, n_spatial, x_holder);
                res.v = x_holder.segment(v_start, n_spatial);

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
    compute_eta_lambda_bym2(Xl, res.beta_lambda, res.v, res.w, a, b,
                            map_site_to_unit, eta_lam);
    eta_p_long.noalias() = Xp * res.beta_p;
    double log_lik_final = nmix_kernel_sweep_spatial(
        obs_by_site, y_R, eta_lam, eta_p_long, K_max, r,
        grad_eta_lam, info_eta_lam, grad_eta_p, info_eta_p,
        mean_N, var_N, boundary_weight, score_wt_lambda
    );
    double log_prior_final = nmix_bym2_log_prior(
        n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        res.v, res.w
    );

    MatrixXd H_final = MatrixXd::Zero(n_x, n_x);
    nmix_assemble_obs_info_bym2(
        p_lam, p_p, n_spatial, a, b,
        Xl, Xp, eta_p_long, obs_by_site, map_site_to_unit,
        info_eta_lam, info_eta_p, var_N, score_wt_lambda, H_final
    );
    nmix_add_bym2_prior_to_H_only(
        p_lam, p_p, n_spatial,
        adj_row_ptr, adj_col_idx, n_neighbors, H_final
    );
    add_diagonal_ridge_bym2(H_final);
    const int p_beta = p_lam + p_p;
    Eigen::LLT<MatrixXd> chol(H_final);
    double log_det_H;
    const int v_start = p_lam + p_p;
    if (chol.info() == Eigen::Success) {
        log_det_H = 2.0 * chol.matrixL().toDenseMatrix().diagonal()
                              .array().log().sum();
        res.cov_beta = nmix_beta_cov_bym2(H_final, n_x, p_beta, v_start, n_spatial);
    } else {
        MatrixXd H_f = MatrixXd::Zero(n_x, n_x);
        nmix_assemble_complete_fisher_bym2(
            p_lam, p_p, n_spatial, a, b,
            Xl, Xp, obs_by_site, map_site_to_unit,
            info_eta_lam, info_eta_p, H_f
        );
        nmix_add_bym2_prior_to_H_only(
            p_lam, p_p, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors, H_f
        );
        add_diagonal_ridge_bym2(H_f);
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
        res.cov_beta = nmix_beta_cov_bym2(H_f, n_x, p_beta, v_start, n_spatial);
    }

    res.log_lik = log_lik_final;
    res.log_marginal = log_lik_final + log_prior_final - 0.5 * log_det_H;
    res.grad_norm = grad_norm;
    res.boundary_max = boundary_weight.maxCoeff();
    return res;
}

}  // namespace

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
    if (scale_factor <= 0) Rcpp::stop("scale_factor must be positive.");
    tulpaObs::NmixSpatialPrep pp = tulpaObs::prep_nmix_spatial(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R, n_spatial,
        r_grid, beta_lambda_init, beta_p_init, K_max);
    const int p_lam = pp.p_lam, p_p = pp.p_p;
    Map<MatrixXd> Xl(REAL(X_lambda_R), pp.n_sites, p_lam);
    Map<MatrixXd> Xp(REAL(X_p_R), pp.n_obs, p_p);

    VectorXd v_default(n_spatial);
    VectorXd w_default(n_spatial);
    if (v_init.isNotNull()) {
        Rcpp::NumericVector vi(v_init);
        if ((int)vi.size() != n_spatial) Rcpp::stop("length(v_init) must equal n_spatial.");
        for (int s = 0; s < n_spatial; ++s) v_default(s) = vi[s];
    } else {
        v_default.setZero();
    }
    if (w_init.isNotNull()) {
        Rcpp::NumericVector wi(w_init);
        if ((int)wi.size() != n_spatial) Rcpp::stop("length(w_init) must equal n_spatial.");
        for (int s = 0; s < n_spatial; ++s) w_default(s) = wi[s];
    } else {
        w_default.setZero();
    }

    // Grid axes: r (outermost) x rho x sigma. BYM2 mixes an ICAR component (v,
    // scaled by a) with an IID component (w, scaled by b); precompute (a, b) and
    // the validity flag per cell so the solver just plugs them into the inner
    // Newton (its field is [v; w], length 2 * n_spatial).
    const int n_sigma = sigma_grid.size(), n_rho = rho_grid.size(), n_r = r_grid.size();
    const int n_grid = n_sigma * n_rho * n_r;
    Rcpp::NumericMatrix theta_grid_out(n_grid, 3);
    std::vector<double> a_k(n_grid), b_k(n_grid), r_k(n_grid);
    std::vector<bool> valid_k(n_grid);
    {
        int k = 0;
        for (int ir_disp = 0; ir_disp < n_r; ++ir_disp)
            for (int ir_rho = 0; ir_rho < n_rho; ++ir_rho)
                for (int sg = 0; sg < n_sigma; ++sg, ++k) {
                    const double sigma = sigma_grid[sg], rho = rho_grid[ir_rho];
                    theta_grid_out(k, 0) = sigma;
                    theta_grid_out(k, 1) = rho;
                    theta_grid_out(k, 2) = r_grid[ir_disp];
                    valid_k[k] = (sigma > 0 && rho >= 0 && rho <= 1);
                    a_k[k] = valid_k[k] ? sigma * std::sqrt(rho / scale_factor) : 0.0;
                    b_k[k] = valid_k[k] ? sigma * std::sqrt(1.0 - rho) : 0.0;
                    r_k[k] = r_grid[ir_disp];
                }
    }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("sigma", "rho", "r");

    auto solve = [&](int k) -> tulpaObs::NmixSpatialPoint {
        tulpaObs::NmixSpatialPoint pt;
        if (!valid_k[k]) { pt.skipped = true; return pt; }  // (sigma, rho) out of range
        VectorXd beta_lam = pp.beta_lam_default;
        VectorXd beta_p   = pp.beta_p_default;
        VectorXd v        = v_default;
        VectorXd w        = w_default;
        BYM2InnerResult ir = inner_newton_bym2(
            p_lam, p_p, pp.n_sites, n_spatial, pp.n_obs,
            a_k[k], b_k[k],
            Xl, Xp, y, pp.obs_by_site, pp.map_site_to_unit,
            adj_row_ptr, adj_col_idx, n_neighbors,
            r_k[k], pp.K_max, max_iter, tol, beta_lam, beta_p, v, w, verbose);
        pt.log_marginal = ir.log_marginal; pt.n_iter = ir.n_iter;
        pt.converged = ir.converged;       pt.grad_norm = ir.grad_norm;
        pt.log_lik = ir.log_lik;           pt.boundary_max = ir.boundary_max;
        pt.coef = VectorXd(p_lam + p_p);
        pt.coef.head(p_lam) = ir.beta_lambda;
        pt.coef.tail(p_p)   = ir.beta_p;
        pt.field = VectorXd(2 * n_spatial);
        pt.field.head(n_spatial) = ir.v;
        pt.field.tail(n_spatial) = ir.w;
        pt.cov_beta = ir.cov_beta;
        return pt;
    };

    Rcpp::List out = tulpaObs::run_nmix_spatial_grid(
        n_grid, p_lam, p_p, n_spatial, /*field_len=*/2 * n_spatial, pp.K_max,
        theta_grid_out, solve,
        progress, progress_every, progress_throttle, progress_file);
    out["sigma_grid"]   = sigma_grid;
    out["rho_grid"]     = rho_grid;
    out["r_grid"]       = r_grid;
    out["scale_factor"] = scale_factor;
    out["n_sigma"]      = n_sigma;
    out["n_rho"]        = n_rho;
    out["n_r"]          = n_r;
    return out;
}
