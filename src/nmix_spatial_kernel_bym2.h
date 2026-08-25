// nmix_spatial_kernel_bym2.h
// BYM2-specific kernel helpers for the Royle (2004) N-mixture model.
//
// Riebler et al. (2016) reparametrise the BYM model so that the joint
// spatial offset is
//   phi[u] = sigma * ( sqrt(rho / s) * v[u] + sqrt(1 - rho) * w[u] )
// with v ~ ICAR (unscaled, sum-to-zero), w ~ N(0, I) iid, and s the
// geometric mean of the non-zero eigenvalues of Q (the ICAR precision).
// Under this parametrisation sigma is the joint marginal standard
// deviation of phi and rho is the spatial fraction of variance.
//
// State vector layout (n_x = p_lam + p_p + 2 * n_spatial):
//
//   x[0           : p_lam]                            = beta_lambda
//   x[p_lam       : p_lam + p_p]                      = beta_p
//   x[p_lam + p_p : p_lam + p_p + n_spatial]          = v  (ICAR)
//   x[p_lam+p_p+n : p_lam + p_p + 2 * n_spatial]      = w  (iid)
//
// Linear predictor at site s mapped to unit u(s):
//
//   eta_lambda[s] = X_lambda[s,] . beta_lambda
//                 + a * v[u(s)] + b * w[u(s)]
//
// with a = sigma * sqrt(rho / scale_factor) and b = sigma * sqrt(1 - rho).
//
// Because (a, b) only multiply v and w in the linear predictor (not in the
// priors), the priors on v and w are *constant* in (sigma, rho) and drop
// out of the integration up to a tau-independent additive constant; the
// grid integration is over (sigma, rho) only.
//
// The (1, b/a)-style identifiability ridge between (lambda intercept) and
// (constant v, constant w) is pinned by sum-to-zero centering of v after
// each Newton step. The iid component w is identified by its N(0, I)
// prior and is not centred.

#ifndef TULPAOBS_NMIX_SPATIAL_KERNEL_BYM2_H
#define TULPAOBS_NMIX_SPATIAL_KERNEL_BYM2_H

#include "nmix_kernel.h"
#include "nmix_spatial_assemble.h"   // nmix_assemble_obs_info / _complete_fisher
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace tulpaObs {

// eta_lambda[s] = X_lambda[s,] . beta_lambda + a * v[u(s)] + b * w[u(s)]
inline void compute_eta_lambda_bym2(
    const Eigen::Map<Eigen::MatrixXd>& X_lambda,
    const Eigen::VectorXd& beta_lambda,
    const Eigen::VectorXd& v,
    const Eigen::VectorXd& w,
    double a, double b,
    const std::vector<int>& map_site_to_unit,
    Eigen::VectorXd& eta_lambda /* out, length n_sites */
) {
    eta_lambda.noalias() = X_lambda * beta_lambda;
    const int n_sites = static_cast<int>(map_site_to_unit.size());
    for (int s = 0; s < n_sites; ++s) {
        const int u = map_site_to_unit[s];
        eta_lambda(s) += a * v(u) + b * w(u);
    }
}

// Assemble marginal observed-info Hessian over
//   theta = (beta_lambda, beta_p, v, w)
// for BYM2. Per-site complete-data Fisher block:
//   ds_beta_lambda = w_lam * x_lam[s,] x_lam[s,]^T
//   ds_v[u]        = w_lam * a^2          (diagonal)
//   ds_w[u]        = w_lam * b^2          (diagonal)
//   ds_v[u]w[u]    = w_lam * a * b        (cross)
//   ds_beta_lambda,v[u] = w_lam * a * x_lam[s,]
//   ds_beta_lambda,w[u] = w_lam * b * x_lam[s,]
// Plus the per-site Var[N|y_s] rank-1 correction with u_s carrying
//   -x_lam[s,] in the beta_lambda slot,
//   sum_j p_ij x_p[ij,] in the beta_p slot,
//   -a at v[u(s)], -b at w[u(s)].
// BYM2 loading: eta_lambda[s] += a * v[u(s)] + b * w[u(s)].
inline SpatialLoading nmix_loading_bym2(
    int v_start, int w_start, double a, double b,
    const std::vector<int>& map_site_to_unit, int s
) {
    const int u = map_site_to_unit[s];
    SpatialLoading L;
    L.n = 2;
    L.col[0] = v_start + u; L.coef[0] = a;
    L.col[1] = w_start + u; L.coef[1] = b;
    return L;
}

inline void nmix_assemble_obs_info_bym2(
    int p_lam, int p_p, int n_spatial,
    double a, double b,
    const Eigen::Map<Eigen::MatrixXd>& X_lambda,
    const Eigen::Map<Eigen::MatrixXd>& X_p,
    const Eigen::VectorXd& eta_p_long,
    const std::vector<std::vector<int>>& obs_by_site,
    const std::vector<int>& map_site_to_unit,
    const Eigen::VectorXd& info_eta_lam,
    const Eigen::VectorXd& info_eta_p,
    const Eigen::VectorXd& var_N,
    const Eigen::VectorXd& score_wt_lambda,   // N-coeff of s_lambda (1 for Poisson)
    Eigen::MatrixXd& H_obs /* in/out: zero-initialized [n_x x n_x] */
) {
    const int v_start = p_lam + p_p;
    const int w_start = v_start + n_spatial;
    const int n_x = p_lam + p_p + 2 * n_spatial;
    nmix_assemble_obs_info(
        p_lam, p_p, n_x, X_lambda, X_p, eta_p_long, obs_by_site,
        info_eta_lam, info_eta_p, var_N, score_wt_lambda,
        [&](int s) { return nmix_loading_bym2(v_start, w_start, a, b, map_site_to_unit, s); },
        H_obs);
}

// Complete-data Fisher Hessian for BYM2 (no var_N correction); always PSD.
inline void nmix_assemble_complete_fisher_bym2(
    int p_lam, int p_p, int n_spatial,
    double a, double b,
    const Eigen::Map<Eigen::MatrixXd>& X_lambda,
    const Eigen::Map<Eigen::MatrixXd>& X_p,
    const std::vector<std::vector<int>>& obs_by_site,
    const std::vector<int>& map_site_to_unit,
    const Eigen::VectorXd& info_eta_lam,
    const Eigen::VectorXd& info_eta_p,
    Eigen::MatrixXd& H_f /* in/out: zero-initialized [n_x x n_x] */
) {
    const int v_start = p_lam + p_p;
    const int w_start = v_start + n_spatial;
    const int n_x = p_lam + p_p + 2 * n_spatial;
    nmix_assemble_complete_fisher(
        p_lam, p_p, n_x, X_lambda, X_p, obs_by_site,
        info_eta_lam, info_eta_p,
        [&](int s) { return nmix_loading_bym2(v_start, w_start, a, b, map_site_to_unit, s); },
        H_f);
}

// Add BYM2 prior contribution to gradient and Hessian:
//   log p(v) = -0.5 * v' Q v       (ICAR; rank-deficient by 1, constant in sigma/rho)
//   log p(w) = -0.5 * w' w         (iid N(0, I))
// Score: grad_v -= Q v,   grad_w -= w
// Hessian: H_{v,v} += Q,   H_{w,w} += I
inline void nmix_add_bym2_prior_to_grad_and_H(
    int p_lam, int p_p, int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const Eigen::VectorXd& v,
    const Eigen::VectorXd& w,
    Eigen::VectorXd& grad /* in/out [n_x] */,
    Eigen::MatrixXd& H    /* in/out [n_x x n_x] */
) {
    const int v_start = p_lam + p_p;
    const int w_start = v_start + n_spatial;

    // v block: the intrinsic case of the CAR precision, tau = rho = 1. Its
    // z-block starts at p_lam + p_p, which is v_start.
    nmix_add_car_to_spatial_block(p_lam, p_p, n_spatial, /*tau=*/1.0,
                                  /*rho=*/1.0, adj_row_ptr, adj_col_idx,
                                  n_neighbors, v, grad, H);

    // w block: identity (iid N(0, I)).
    for (int s = 0; s < n_spatial; ++s) {
        const int idx_ws = w_start + s;
        H(idx_ws, idx_ws) += 1.0;
        grad(idx_ws) -= w(s);
    }
}

// Hessian-only variant for the final log|H| assembly.
inline void nmix_add_bym2_prior_to_H_only(
    int p_lam, int p_p, int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    Eigen::MatrixXd& H /* in/out [n_x x n_x] */
) {
    const int w_start = p_lam + p_p + n_spatial;

    nmix_add_car_to_H_only(p_lam, p_p, n_spatial, /*tau=*/1.0, /*rho=*/1.0,
                           adj_row_ptr, adj_col_idx, n_neighbors, H);
    for (int s = 0; s < n_spatial; ++s) {
        H(w_start + s, w_start + s) += 1.0;
    }
}

// log p(v) + log p(w) for BYM2, dropping (sigma, rho)-independent additive
// constants. Used inside the line-search objective only -- the absolute
// value cancels in the Laplace marginal because it does not depend on
// (sigma, rho).
inline double nmix_bym2_log_prior(
    int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const Eigen::VectorXd& v,
    const Eigen::VectorXd& w
) {
    const double quad_v = nmix_car_quadratic_form(
        n_spatial, 1.0, adj_row_ptr, adj_col_idx, n_neighbors, v);
    const double quad_w = w.squaredNorm();
    return -0.5 * quad_v - 0.5 * quad_w;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_SPATIAL_KERNEL_BYM2_H
