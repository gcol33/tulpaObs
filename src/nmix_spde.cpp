// nmix_spde.cpp
// Nested Laplace for the spatial Royle (2004) N-mixture model with a continuous
// Matern (SPDE) field on the abundance arm:
//
//   N_i ~ Poisson(lambda_i)              (mixture = "P",  r = +Inf)
//   N_i ~ NegBin(mean lambda_i, size r)  (mixture = "NB", r finite)
//   log lambda_i = X_lambda[i,] . beta_lambda + (A u)_i
//   logit p_ij   = X_p[ij,]   . beta_p
//
// The field u lives at the n_mesh FEM nodes; A (n_sites x n_mesh) projects mesh
// nodes onto site locations. u ~ N(0, Q(range, sigma)^{-1}) with the proper
// Matern FEM precision Q. The outer grid integrates the SPDE hyperparameters
// (range, sigma) -- the continuous-field analogue of the areal (tau[, rho])
// grid in nmix_spatial.cpp. Q (and log|Q|) are built once per grid point on the
// R side (the same spde_precision_Q the occupancy SPDE path uses) and passed in
// as a dense symmetric matrix, so this kernel stays agnostic to how the
// precision is parameterised.
//
// Per-grid log marginal (up to a grid-independent additive constant):
//   log p(y | range_k, sigma_k) ~= log L(mode_k)
//                                  + log p(u_mode_k | range_k, sigma_k)
//                                  - 0.5 log|H(mode_k)|
// with log p(u | .) = 0.5 log|Q| - 0.5 u' Q u (the (2 pi)^{-n_mesh/2} constant
// is grid-independent -- n_mesh is fixed across the grid -- and is dropped, as
// in nmix_spatial.cpp).
//
// Unlike the intrinsic ICAR field, Q is full rank (proper Matern), so the
// (intercept, field-mean) direction is identified by Q itself: no sum-to-zero
// centering and no constrained-covariance projection (the CAR_proper path is
// the areal analogue).

#include "nmix_kernel.h"
#include "nmix_spatial_kernel.h"   // nmix_kernel_sweep_spatial / _log_lik_only_spatial
#include "nmix_linalg.h"
#include "nmix_progress.h"
#include "newton_step.h"           // newton_backtrack / solve_with_fisher_fallback
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
using tulpaObs::compute_nmix_site;
using tulpaObs::NMixSiteResult;
using tulpaObs::nmix_add_diagonal_ridge;

// eta_lambda[s] = X_lambda[s,] . beta_lambda + (A u)_s.
inline void compute_eta_lambda_spde(
    const Map<MatrixXd>& Xl, const VectorXd& beta_lambda,
    const MatrixXd& A, const VectorXd& u,
    VectorXd& eta_lambda /* out, length n_sites */
) {
    eta_lambda.noalias() = Xl * beta_lambda;
    eta_lambda.noalias() += A * u;
}

// Per-site marginal sweep -- fills the per-site / per-visit gradient and Fisher
// vectors plus the Var[N|y] moments. Returns total log-lik.
inline double nmix_sweep_spde(
    const std::vector<std::vector<int>>& obs_by_site,
    const Rcpp::IntegerVector& y_R,
    const VectorXd& eta_lambda, const VectorXd& eta_p_long,
    int K_max, double r,
    VectorXd& grad_eta_lam, VectorXd& info_eta_lam,
    VectorXd& grad_eta_p,   VectorXd& info_eta_p,
    VectorXd& mean_N,       VectorXd& var_N,
    VectorXd& boundary_weight, VectorXd& score_wt_lambda
) {
    return tulpaObs::nmix_kernel_sweep_spatial(
        obs_by_site, y_R, eta_lambda, eta_p_long, K_max, r,
        grad_eta_lam, info_eta_lam, grad_eta_p, info_eta_p,
        mean_N, var_N, boundary_weight, score_wt_lambda);
}

// Assemble the marginal observed-information Hessian over
//   x = (beta_lambda [p_lam], beta_p [p_p], u [n_mesh])
// with the field loading through the dense projection A (a site's eta_lambda
// loads onto every mesh node with a non-zero A row entry -- the SPDE
// generalisation of the areal single-unit loading). Per site s:
//   complete-data Fisher D_s (block (beta_lambda, beta_p) + lambda<->u cross)
//   minus the Var[N|y_s] rank-1 correction u_s u_s' (Louis 1982),
// where the score vector u_s carries -score_wt_lambda along (beta_lambda, A row)
// and +sum_j p_ij x_p on the detection block.
inline void assemble_obs_info_spde(
    int p_lam, int p_p, int n_mesh,
    const Map<MatrixXd>& Xl, const Map<MatrixXd>& Xp,
    const MatrixXd& A,
    const VectorXd& eta_p_long,
    const std::vector<std::vector<int>>& obs_by_site,
    const VectorXd& info_eta_lam, const VectorXd& info_eta_p,
    const VectorXd& var_N, const VectorXd& score_wt_lambda,
    MatrixXd& H /* in/out: zero-initialised [n_x x n_x] */
) {
    const int n_sites = static_cast<int>(obs_by_site.size());
    const int u_start = p_lam + p_p;
    const int n_x = u_start + n_mesh;

    for (int s = 0; s < n_sites; ++s) {
        const auto& idx = obs_by_site[s];
        const int J = static_cast<int>(idx.size());
        if (J == 0) continue;

        const double w_lam = info_eta_lam(s);
        const VectorXd a_s = A.row(s).transpose();   // n_mesh field loading

        // ----- Complete-data Fisher D_s -----
        if (w_lam > 0.0) {
            // beta_lambda x beta_lambda
            H.block(0, 0, p_lam, p_lam)
                .selfadjointView<Eigen::Lower>()
                .rankUpdate(Xl.row(s).transpose(), w_lam);
            // u x u  (field block): w_lam * a_s a_s'
            H.block(u_start, u_start, n_mesh, n_mesh)
                .selfadjointView<Eigen::Lower>()
                .rankUpdate(a_s, w_lam);
            // beta_lambda x u cross into the LOWER-LEFT block (u rows, beta cols)
            H.block(u_start, 0, n_mesh, p_lam).noalias()
                += w_lam * a_s * Xl.row(s);
        }
        for (int j = 0; j < J; ++j) {
            const double w_p = info_eta_p(idx[j]);
            if (w_p > 0.0) {
                H.block(p_lam, p_lam, p_p, p_p)
                    .selfadjointView<Eigen::Lower>()
                    .rankUpdate(Xp.row(idx[j]).transpose(), w_p);
            }
        }

        // ----- Var[N|y_s] rank-1 marginal correction -----
        if (var_N(s) > 0.0) {
            const double swl = score_wt_lambda(s);
            VectorXd u_vec = VectorXd::Zero(n_x);
            u_vec.segment(0, p_lam) = -swl * Xl.row(s).transpose();
            VectorXd ssum = VectorXd::Zero(p_p);
            for (int j = 0; j < J; ++j) {
                const double e = eta_p_long(idx[j]);
                double p_ij;
                if (e > 0.0) p_ij = 1.0 / (1.0 + std::exp(-e));
                else { const double ee = std::exp(e); p_ij = ee / (1.0 + ee); }
                ssum += p_ij * Xp.row(idx[j]).transpose();
            }
            u_vec.segment(p_lam, p_p) = ssum;
            u_vec.segment(u_start, n_mesh) = -swl * a_s;
            H.selfadjointView<Eigen::Lower>().rankUpdate(u_vec, -var_N(s));
        }
    }
    // Every contribution above wrote the lower triangle (the selfadjoint rank
    // updates and the beta_lambda x u cross in the lower-left block). The
    // detection x detection block sits inside the lower triangle too. Mirror
    // the full lower triangle to the upper.
    H = H.selfadjointView<Eigen::Lower>();
}

// Complete-data Fisher Hessian (no var_N correction). Always PSD; the
// Levenberg-Marquardt fallback when the observed-info matrix is not PSD.
inline void assemble_complete_fisher_spde(
    int p_lam, int p_p, int n_mesh,
    const Map<MatrixXd>& Xl, const Map<MatrixXd>& Xp,
    const MatrixXd& A,
    const std::vector<std::vector<int>>& obs_by_site,
    const VectorXd& info_eta_lam, const VectorXd& info_eta_p,
    MatrixXd& H /* in/out: zero-initialised [n_x x n_x] */
) {
    const int n_sites = static_cast<int>(obs_by_site.size());
    const int u_start = p_lam + p_p;

    for (int s = 0; s < n_sites; ++s) {
        const double w_lam = info_eta_lam(s);
        if (w_lam <= 0.0) continue;
        const VectorXd a_s = A.row(s).transpose();
        H.block(0, 0, p_lam, p_lam)
            .selfadjointView<Eigen::Lower>()
            .rankUpdate(Xl.row(s).transpose(), w_lam);
        H.block(u_start, u_start, n_mesh, n_mesh)
            .selfadjointView<Eigen::Lower>()
            .rankUpdate(a_s, w_lam);
        H.block(u_start, 0, n_mesh, p_lam).noalias()
            += w_lam * a_s * Xl.row(s);
    }
    for (int o = 0; o < static_cast<int>(Xp.rows()); ++o) {
        const double w_p = info_eta_p(o);
        if (w_p > 0.0) {
            H.block(p_lam, p_lam, p_p, p_p)
                .selfadjointView<Eigen::Lower>()
                .rankUpdate(Xp.row(o).transpose(), w_p);
        }
    }
    H = H.selfadjointView<Eigen::Lower>();
}

struct SpdeInnerResult {
    VectorXd beta_lambda;
    VectorXd beta_p;
    VectorXd u;
    MatrixXd cov_beta;    // (p_lam + p_p) marginal coefficient covariance at mode
    double log_lik;
    double log_marginal;
    double grad_norm;
    int n_iter;
    bool converged;
    double boundary_max;
};

// Inner Newton at a fixed (range, sigma) -- supplied as the full precision Q
// (n_mesh x n_mesh, already carrying tau_spde^2) and its log-determinant. The
// field is full rank, so no centering / constraint is applied.
SpdeInnerResult inner_newton_spde(
    int p_lam, int p_p, int n_sites, int n_mesh, int n_obs,
    const Map<MatrixXd>& Xl, const Map<MatrixXd>& Xp,
    const MatrixXd& A,
    const Rcpp::IntegerVector& y_R,
    const std::vector<std::vector<int>>& obs_by_site,
    const MatrixXd& Q, double log_det_Q,
    double r, int K_max, int max_iter, double tol,
    const VectorXd& beta_lam_init, const VectorXd& beta_p_init,
    const VectorXd& u_init, bool verbose
) {
    SpdeInnerResult res;
    res.beta_lambda = beta_lam_init;
    res.beta_p      = beta_p_init;
    res.u           = u_init;
    res.converged   = false;
    res.n_iter      = 0;

    const int u_start = p_lam + p_p;
    const int n_x = u_start + n_mesh;

    VectorXd grad_eta_lam(n_sites), info_eta_lam(n_sites);
    VectorXd mean_N(n_sites), var_N(n_sites), boundary_weight(n_sites);
    VectorXd score_wt_lambda(n_sites);
    VectorXd grad_eta_p(n_obs), info_eta_p(n_obs);
    VectorXd eta_lam(n_sites), eta_p_long(n_obs);

    double log_lik = R_NegInf;
    double grad_norm = R_PosInf;

    auto field_log_prior = [&](const VectorXd& u) {
        return 0.5 * log_det_Q - 0.5 * (u.transpose() * (Q * u))(0, 0);
    };

    for (int iter = 0; iter < max_iter; ++iter) {
        compute_eta_lambda_spde(Xl, res.beta_lambda, A, res.u, eta_lam);
        eta_p_long.noalias() = Xp * res.beta_p;

        log_lik = nmix_sweep_spde(
            obs_by_site, y_R, eta_lam, eta_p_long, K_max, r,
            grad_eta_lam, info_eta_lam, grad_eta_p, info_eta_p,
            mean_N, var_N, boundary_weight, score_wt_lambda);

        VectorXd grad = VectorXd::Zero(n_x);
        grad.segment(0, p_lam)   = Xl.transpose() * grad_eta_lam;
        grad.segment(p_lam, p_p) = Xp.transpose() * grad_eta_p;
        grad.segment(u_start, n_mesh) = A.transpose() * grad_eta_lam;
        // Field prior contributes -Q u to the score.
        grad.segment(u_start, n_mesh) -= Q * res.u;

        grad_norm = grad.norm();
        if (verbose) {
            Rcpp::Rcout << "  iter " << iter << "  log_lik " << log_lik
                        << "  grad_norm " << grad_norm
                        << "  boundary_max " << boundary_weight.maxCoeff() << "\n";
        }
        if (grad_norm < tol) {
            res.converged = true;
            res.n_iter = iter + 1;
            break;
        }

        MatrixXd H = MatrixXd::Zero(n_x, n_x);
        assemble_obs_info_spde(p_lam, p_p, n_mesh, Xl, Xp, A, eta_p_long,
                               obs_by_site, info_eta_lam, info_eta_p,
                               var_N, score_wt_lambda, H);
        H.block(u_start, u_start, n_mesh, n_mesh) += Q;
        nmix_add_diagonal_ridge(H);

        VectorXd delta;
        const bool solved = tulpaObs::solve_with_fisher_fallback(
            H, grad,
            [&]() {
                MatrixXd H_f = MatrixXd::Zero(n_x, n_x);
                assemble_complete_fisher_spde(p_lam, p_p, n_mesh, Xl, Xp, A,
                                              obs_by_site, info_eta_lam,
                                              info_eta_p, H_f);
                H_f.block(u_start, u_start, n_mesh, n_mesh) += Q;
                nmix_add_diagonal_ridge(H_f);
                return H_f;
            },
            "the SPDE grid point", iter, verbose, delta);
        if (!solved) break;

        const VectorXd delta_lam = delta.segment(0, p_lam);
        const VectorXd delta_p   = delta.segment(p_lam, p_p);
        const VectorXd delta_u   = delta.segment(u_start, n_mesh);

        VectorXd beta_lam_try, beta_p_try, u_try;
        VectorXd eta_lam_try(n_sites), eta_p_try(n_obs);
        const double obj_cur = log_lik + field_log_prior(res.u);
        const bool stepped = tulpaObs::newton_backtrack(
            obj_cur,
            [&](double step) {
                beta_lam_try = res.beta_lambda + step * delta_lam;
                beta_p_try   = res.beta_p      + step * delta_p;
                u_try        = res.u           + step * delta_u;
                compute_eta_lambda_spde(Xl, beta_lam_try, A, u_try, eta_lam_try);
                eta_p_try.noalias() = Xp * beta_p_try;
                const double ll_try = tulpaObs::nmix_kernel_log_lik_only_spatial(
                    obs_by_site, y_R, eta_lam_try, eta_p_try, K_max, r);
                return ll_try + field_log_prior(u_try);
            },
            [&](double) {
                res.beta_lambda = beta_lam_try;
                res.beta_p      = beta_p_try;
                res.u           = u_try;
            });
        if (!stepped) {
            if (verbose) Rcpp::Rcout << "    (step halving exhausted)\n";
            break;
        }
        res.n_iter = iter + 1;
    }

    // ----- log marginal at the converged mode -----
    compute_eta_lambda_spde(Xl, res.beta_lambda, A, res.u, eta_lam);
    eta_p_long.noalias() = Xp * res.beta_p;
    const double log_lik_final = nmix_sweep_spde(
        obs_by_site, y_R, eta_lam, eta_p_long, K_max, r,
        grad_eta_lam, info_eta_lam, grad_eta_p, info_eta_p,
        mean_N, var_N, boundary_weight, score_wt_lambda);
    const double log_prior_u_final = field_log_prior(res.u);

    MatrixXd H_final = MatrixXd::Zero(n_x, n_x);
    assemble_obs_info_spde(p_lam, p_p, n_mesh, Xl, Xp, A, eta_p_long,
                           obs_by_site, info_eta_lam, info_eta_p,
                           var_N, score_wt_lambda, H_final);
    H_final.block(u_start, u_start, n_mesh, n_mesh) += Q;
    nmix_add_diagonal_ridge(H_final);

    const int p_beta = p_lam + p_p;
    Eigen::LLT<MatrixXd> chol(H_final);
    double log_det_H;
    if (chol.info() == Eigen::Success) {
        log_det_H = 2.0 * chol.matrixL().toDenseMatrix().diagonal()
                              .array().log().sum();
        // Full-rank field: no constraint (CAR_proper analogue).
        res.cov_beta = tulpaObs::nmix_constrained_top_cov(
            H_final, n_x, p_beta, p_beta, n_mesh, /*constrain=*/false);
    } else {
        MatrixXd H_f = MatrixXd::Zero(n_x, n_x);
        assemble_complete_fisher_spde(p_lam, p_p, n_mesh, Xl, Xp, A,
                                      obs_by_site, info_eta_lam,
                                      info_eta_p, H_f);
        H_f.block(u_start, u_start, n_mesh, n_mesh) += Q;
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
        res.cov_beta = tulpaObs::nmix_constrained_top_cov(
            H_f, n_x, p_beta, p_beta, n_mesh, /*constrain=*/false);
    }

    res.log_lik = log_lik_final;
    res.log_marginal = log_lik_final + log_prior_u_final - 0.5 * log_det_H;
    res.grad_norm = grad_norm;
    res.boundary_max = boundary_weight.maxCoeff();
    return res;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_nested_laplace_nmix_spde(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::NumericMatrix X_p_R,
    Rcpp::NumericMatrix A_R,          // dense [n_sites x n_mesh] projection
    Rcpp::List Q_list,                // per-grid-point precision (dense [n_mesh x n_mesh])
    Rcpp::NumericVector log_det_Q,    // per-grid-point log|Q|
    Rcpp::NumericMatrix theta_grid_R, // [n_grid x n_theta], named cols
    Rcpp::NumericVector r_grid,       // NB size per grid point (or +Inf), length n_grid
    Rcpp::NumericVector beta_lambda_init,
    Rcpp::NumericVector beta_p_init,
    Rcpp::Nullable<Rcpp::NumericVector> u_init = R_NilValue,
    int K_max = -1,
    int max_iter = 100,
    double tol = 1e-6,
    bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = ""
) {
    const int n_sites = X_lambda_R.nrow();
    const int p_lam   = X_lambda_R.ncol();
    const int n_obs   = X_p_R.nrow();
    const int p_p     = X_p_R.ncol();
    const int n_mesh  = A_R.ncol();
    if ((int)y.size() != n_obs) Rcpp::stop("length(y) must equal nrow(X_p).");
    if ((int)site_idx.size() != n_obs) Rcpp::stop("length(site_idx) must equal nrow(X_p).");
    if (A_R.nrow() != n_sites) Rcpp::stop("nrow(A) must equal nrow(X_lambda).");
    if ((int)beta_lambda_init.size() != p_lam) Rcpp::stop("beta_lambda_init length mismatch.");
    if ((int)beta_p_init.size() != p_p) Rcpp::stop("beta_p_init length mismatch.");
    const int n_grid = Q_list.size();
    if ((int)log_det_Q.size() != n_grid) Rcpp::stop("length(log_det_Q) must equal length(Q_list).");
    if (theta_grid_R.nrow() != n_grid) Rcpp::stop("nrow(theta_grid) must equal length(Q_list).");
    if ((int)r_grid.size() != n_grid) Rcpp::stop("length(r_grid) must equal length(Q_list).");
    if (K_max < 0) {
        int ymax = 0;
        for (int o = 0; o < n_obs; ++o) if (y[o] > ymax) ymax = y[o];
        K_max = ymax + 100;
    }

    std::vector<std::vector<int>> obs_by_site(n_sites);
    for (int o = 0; o < n_obs; ++o) {
        const int s = site_idx[o] - 1;
        if (s < 0 || s >= n_sites) Rcpp::stop("site_idx out of range at obs %d.", o + 1);
        obs_by_site[s].push_back(o);
    }

    Map<MatrixXd> Xl(REAL(X_lambda_R), n_sites, p_lam);
    Map<MatrixXd> Xp(REAL(X_p_R), n_obs, p_p);
    Map<MatrixXd> A(REAL(A_R), n_sites, n_mesh);

    VectorXd beta_lam_default = Map<VectorXd>(REAL(beta_lambda_init), p_lam);
    VectorXd beta_p_default   = Map<VectorXd>(REAL(beta_p_init), p_p);
    VectorXd u_default(n_mesh);
    if (u_init.isNotNull()) {
        Rcpp::NumericVector ui(u_init);
        if ((int)ui.size() != n_mesh) Rcpp::stop("length(u_init) must equal n_mesh.");
        for (int m = 0; m < n_mesh; ++m) u_default(m) = ui[m];
    } else {
        u_default.setZero();
    }

    const int n_theta = theta_grid_R.ncol();
    const int n_x = p_lam + p_p + n_mesh;

    Rcpp::NumericVector log_marginals(n_grid);
    Rcpp::IntegerVector n_iters(n_grid);
    Rcpp::LogicalVector convergeds(n_grid);
    Rcpp::NumericVector grad_norms(n_grid);
    Rcpp::NumericVector log_liks(n_grid);
    Rcpp::NumericVector boundary_maxes(n_grid);
    Rcpp::NumericMatrix modes(n_grid, n_x);
    Rcpp::List cov_blocks(n_grid);

    // outer-grid progress (tulpa#45)
    auto gp = tulpaObs::make_grid_progress("nmix-spatial", n_grid, progress,
                                           progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        Rcpp::NumericMatrix Qk_R = Q_list[k];
        if (Qk_R.nrow() != n_mesh || Qk_R.ncol() != n_mesh) {
            Rcpp::stop("Q_list[[%d]] must be n_mesh x n_mesh.", k + 1);
        }
        Map<MatrixXd> Qk(REAL(Qk_R), n_mesh, n_mesh);
        const double rr = r_grid[k];

        // Cold restart per grid point (the lambda/p intercept identifiability
        // ridge shifts across the hyperparameters; same rationale as the areal
        // path).
        VectorXd beta_lam = beta_lam_default;
        VectorXd beta_p   = beta_p_default;
        VectorXd u        = u_default;

        SpdeInnerResult ir = inner_newton_spde(
            p_lam, p_p, n_sites, n_mesh, n_obs,
            Xl, Xp, A, y, obs_by_site,
            Qk, log_det_Q[k], rr, K_max, max_iter, tol,
            beta_lam, beta_p, u, verbose);

        log_marginals[k]  = ir.log_marginal;
        n_iters[k]        = ir.n_iter;
        convergeds[k]     = ir.converged;
        grad_norms[k]     = ir.grad_norm;
        log_liks[k]       = ir.log_lik;
        boundary_maxes[k] = ir.boundary_max;
        for (int j = 0; j < p_lam; ++j)    modes(k, j) = ir.beta_lambda(j);
        for (int j = 0; j < p_p; ++j)      modes(k, p_lam + j) = ir.beta_p(j);
        for (int j = 0; j < n_mesh; ++j)   modes(k, p_lam + p_p + j) = ir.u(j);
        cov_blocks[k] = Rcpp::wrap(ir.cov_beta);

        if (verbose) {
            Rcpp::Rcout << "[grid " << k + 1 << "/" << n_grid
                        << "] log_marg=" << ir.log_marginal
                        << " n_iter=" << ir.n_iter
                        << " conv=" << ir.converged << "\n";
        }
        if (gp) gp->tick();
    }
    if (gp) gp->finish();

    return Rcpp::List::create(
        Rcpp::Named("theta_grid")   = theta_grid_R,
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
        Rcpp::Named("n_mesh")       = n_mesh,
        Rcpp::Named("n_grid")       = n_grid,
        Rcpp::Named("n_theta")      = n_theta,
        Rcpp::Named("K_max")        = K_max
    );
}
