// nmix_count_spatial_driver.h
// Family-agnostic nested-Laplace inner Newton for an areal (ICAR / proper-CAR)
// field on the abundance arm of a count-marginal model. The per-site marginal
// is supplied as a callable `site_fn` with the compute_nmix_site signature
// (const int* y, const double* eta_p, int J, double eta_lambda, int K_max,
// double r) -> NMixSiteResult, so the Royle N-mixture, removal, and distance
// families share ONE driver -- only the per-site marginal differs. Extracted
// from nmix_spatial.cpp (#51); the N-mixture driver instantiates it with the
// Royle kernel and is byte-identical to before.
//
// State vector x = (beta_lambda [p_lam], beta_p [p_p], z [n_spatial]); the field
// loads onto eta_lambda[s] += z[u(s)]. The gradient, observed-information
// Hessian (with the Var[N|y] rank-1 correction), CAR prior contribution, step-
// halving line search, sum-to-zero centering (ICAR), and Laplace log-marginal
// are all family-agnostic (they read only the per-site moment vectors the sweep
// fills).

#ifndef TULPAOBS_NMIX_COUNT_SPATIAL_DRIVER_H
#define TULPAOBS_NMIX_COUNT_SPATIAL_DRIVER_H

#include "nmix_kernel.h"
#include "nmix_spatial_kernel.h"
#include "nmix_spatial_kernel_bym2.h"  // BYM2 eta / assembler / prior helpers
#include "nmix_linalg.h"
#include "nmix_spatial_grid.h"     // prep_nmix_spatial / run_nmix_spatial_grid / NmixSpatialPoint
#include "newton_step.h"           // newton_backtrack / solve_with_fisher_fallback
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <vector>

namespace tulpaObs {

enum class CarPriorKind { ICAR, CAR_PROPER };

inline double count_car_log_prior_dispatch(
    CarPriorKind kind,
    int n_spatial, double tau, double rho, double log_det_Q_rho,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const Eigen::VectorXd& z
) {
    if (kind == CarPriorKind::ICAR) {
        return nmix_icar_log_prior(
            n_spatial, tau, adj_row_ptr, adj_col_idx, n_neighbors, z);
    }
    return nmix_car_proper_log_prior(
        n_spatial, tau, rho, log_det_Q_rho,
        adj_row_ptr, adj_col_idx, n_neighbors, z);
}

// Which parameterisation the latent field carries. The two enter eta_lambda
// differently and carry different priors; nothing else about the inner solve
// depends on the choice.
//
//   CAR   one block z of n_spatial coordinates, eta += z[unit(s)], prior
//         tau Q(rho) (ICAR is rho = 1, and pins the flat direction by centering).
//   BYM2  two blocks [v | w] of n_spatial each, eta += a v[unit] + b w[unit],
//         prior v ~ ICAR(1), w ~ N(0, I) (Riebler 2016). The structured block v
//         carries the same flat direction, so it centers too.
enum class CountFieldKind { CAR, BYM2 };

struct CountFieldSpec {
    CountFieldKind kind = CountFieldKind::CAR;
    CarPriorKind   car  = CarPriorKind::ICAR;      // CAR only
    double tau = 1.0, rho = 1.0, log_det_Q_rho = 0.0;  // CAR only
    double a = 0.0, b = 0.0;                       // BYM2 loadings

    int  n_blocks() const { return kind == CountFieldKind::BYM2 ? 2 : 1; }
    // Whether the field carries a flat direction the prior does not pin. Both
    // the intrinsic CAR and the BYM2 structured block do.
    bool centers() const {
        return kind == CountFieldKind::BYM2 || car == CarPriorKind::ICAR;
    }
};

// Per-grid-point inner solve result. `field` is the whole latent field: the CAR
// block z, or the BYM2 pair [v | w] laid out back to back.
struct CountSpatialInnerResult {
    Eigen::VectorXd beta_lambda;
    Eigen::VectorXd beta_p;
    Eigen::VectorXd field;
    Eigen::MatrixXd cov_beta;    // (p_lam+p_p) marginal coefficient covariance
    double log_lik;
    double log_marginal;
    double grad_norm;
    int n_iter;
    bool converged;
    double boundary_max;
};

// Inner Newton for one outer grid point: joint mode of (beta_lambda, beta_p,
// field) under the count marginal and the field prior, then the Laplace log
// marginal at that mode.
//
// The loop is one solve for both field kinds. Every step that depends on the
// parameterisation goes through the small lambdas below -- how the field enters
// eta_lambda, how its score scatters, how its prior enters the gradient, the
// observed information and the complete-data Fisher, and what its log prior is
// -- so a third kind is those five lambdas, not a third copy of the loop.
template <class SiteFn>
CountSpatialInnerResult inner_newton_count_field(
    const CountFieldSpec& fs,
    int p_lam, int p_p, int n_sites, int n_spatial, int n_obs,
    const Eigen::Map<Eigen::MatrixXd>& Xl,
    const Eigen::Map<Eigen::MatrixXd>& Xp,
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
    const Eigen::VectorXd& beta_lam_init,
    const Eigen::VectorXd& beta_p_init,
    const Eigen::VectorXd& field_init,
    bool verbose,
    SiteFn site_fn
) {
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const bool bym2    = (fs.kind == CountFieldKind::BYM2);
    const int  n_field = fs.n_blocks() * n_spatial;
    const int  f_start = p_lam + p_p;
    const int  w_start = f_start + n_spatial;      // BYM2 only
    const int  n_x     = f_start + n_field;

    CountSpatialInnerResult res;
    res.beta_lambda = beta_lam_init;
    res.beta_p      = beta_p_init;
    res.field       = field_init;
    res.converged   = false;
    res.n_iter      = 0;

    // --- the five field-kind-dependent steps -------------------------------
    auto eta_lambda_of = [&](const VectorXd& bl, const VectorXd& f,
                             VectorXd& out) {
        if (bym2) {
            compute_eta_lambda_bym2(Xl, bl, f.head(n_spatial), f.tail(n_spatial),
                                    fs.a, fs.b, map_site_to_unit, out);
        } else {
            compute_eta_lambda_spatial(Xl, bl, f, map_site_to_unit, out);
        }
    };
    auto scatter_field_score = [&](const VectorXd& grad_eta_lam, VectorXd& grad) {
        for (int s = 0; s < n_sites; ++s) {
            const int u = map_site_to_unit[s];
            if (bym2) {
                grad(f_start + u) += fs.a * grad_eta_lam(s);
                grad(w_start + u) += fs.b * grad_eta_lam(s);
            } else {
                grad(f_start + u) += grad_eta_lam(s);
            }
        }
    };
    auto add_obs_info = [&](const VectorXd& eta_p_long,
                            const VectorXd& info_eta_lam,
                            const VectorXd& info_eta_p,
                            const VectorXd& var_N,
                            const VectorXd& score_wt_lambda, MatrixXd& H) {
        if (bym2) {
            nmix_assemble_obs_info_bym2(p_lam, p_p, n_spatial, fs.a, fs.b,
                Xl, Xp, eta_p_long, obs_by_site, map_site_to_unit,
                info_eta_lam, info_eta_p, var_N, score_wt_lambda, H);
        } else {
            nmix_assemble_obs_info_spatial(p_lam, p_p, n_spatial,
                Xl, Xp, eta_p_long, obs_by_site, map_site_to_unit,
                info_eta_lam, info_eta_p, var_N, score_wt_lambda, H);
        }
    };
    auto add_complete_fisher = [&](const VectorXd& info_eta_lam,
                                   const VectorXd& info_eta_p, MatrixXd& H_f) {
        if (bym2) {
            nmix_assemble_complete_fisher_bym2(p_lam, p_p, n_spatial, fs.a, fs.b,
                Xl, Xp, obs_by_site, map_site_to_unit,
                info_eta_lam, info_eta_p, H_f);
        } else {
            nmix_assemble_complete_fisher_spatial(p_lam, p_p, n_spatial,
                Xl, Xp, obs_by_site, map_site_to_unit,
                info_eta_lam, info_eta_p, H_f);
        }
    };
    auto add_prior_grad_H = [&](const VectorXd& f, VectorXd& grad, MatrixXd& H) {
        if (bym2) {
            nmix_add_bym2_prior_to_grad_and_H(p_lam, p_p, n_spatial,
                adj_row_ptr, adj_col_idx, n_neighbors,
                f.head(n_spatial), f.tail(n_spatial), grad, H);
        } else {
            nmix_add_car_to_spatial_block(p_lam, p_p, n_spatial, fs.tau, fs.rho,
                adj_row_ptr, adj_col_idx, n_neighbors, f, grad, H);
        }
    };
    auto add_prior_H = [&](MatrixXd& H) {
        if (bym2) {
            nmix_add_bym2_prior_to_H_only(p_lam, p_p, n_spatial,
                adj_row_ptr, adj_col_idx, n_neighbors, H);
        } else {
            nmix_add_car_to_H_only(p_lam, p_p, n_spatial, fs.tau, fs.rho,
                adj_row_ptr, adj_col_idx, n_neighbors, H);
        }
    };
    auto field_log_prior = [&](const VectorXd& f) {
        return bym2
            ? nmix_bym2_log_prior(n_spatial, adj_row_ptr, adj_col_idx,
                                  n_neighbors, f.head(n_spatial), f.tail(n_spatial))
            : count_car_log_prior_dispatch(fs.car, n_spatial, fs.tau, fs.rho,
                                           fs.log_det_Q_rho, adj_row_ptr,
                                           adj_col_idx, n_neighbors, f);
    };
    // Pin the flat field direction by centering the structured block, moving
    // its level into the intercept the coefficient block already carries.
    auto center = [&](const VectorXd& bl, const VectorXd& bp, VectorXd& f) {
        if (!fs.centers()) return;
        VectorXd x_holder(n_x);
        x_holder.segment(0, p_lam)   = bl;
        x_holder.segment(p_lam, p_p) = bp;
        x_holder.segment(f_start, n_field) = f;
        nmix_center_field(p_lam, p_p, n_spatial, x_holder);
        f.head(n_spatial) = x_holder.segment(f_start, n_spatial);
    };

    VectorXd grad_eta_lam(n_sites), info_eta_lam(n_sites);
    VectorXd mean_N(n_sites), var_N(n_sites), boundary_weight(n_sites);
    VectorXd score_wt_lambda(n_sites);
    VectorXd grad_eta_p(n_obs), info_eta_p(n_obs);
    VectorXd eta_lam(n_sites);
    VectorXd eta_p_long(n_obs);

    double log_lik = R_NegInf;
    double grad_norm = R_PosInf;
    const std::string grid_label = bym2
        ? newton_grid_label("a", fs.a, "b", fs.b)
        : newton_grid_label("tau", fs.tau, "rho", fs.rho);

    for (int iter = 0; iter < max_iter; ++iter) {
        eta_lambda_of(res.beta_lambda, res.field, eta_lam);
        eta_p_long.noalias() = Xp * res.beta_p;

        log_lik = count_kernel_sweep_spatial(
            obs_by_site, y_R, eta_lam, eta_p_long, K_max, r,
            grad_eta_lam, info_eta_lam, grad_eta_p, info_eta_p,
            mean_N, var_N, boundary_weight, score_wt_lambda, site_fn
        );

        VectorXd grad = VectorXd::Zero(n_x);
        grad.segment(0, p_lam)   = Xl.transpose() * grad_eta_lam;
        grad.segment(p_lam, p_p) = Xp.transpose() * grad_eta_p;
        scatter_field_score(grad_eta_lam, grad);

        MatrixXd H = MatrixXd::Zero(n_x, n_x);
        add_obs_info(eta_p_long, info_eta_lam, info_eta_p, var_N,
                     score_wt_lambda, H);
        add_prior_grad_H(res.field, grad, H);

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

        nmix_add_diagonal_ridge(H);
        VectorXd delta;
        const bool solved = solve_with_fisher_fallback(
            H, grad,
            [&]() {
                MatrixXd H_f = MatrixXd::Zero(n_x, n_x);
                add_complete_fisher(info_eta_lam, info_eta_p, H_f);
                add_prior_H(H_f);
                nmix_add_diagonal_ridge(H_f);
                return H_f;
            },
            grid_label, iter, verbose, delta);
        if (!solved) break;

        VectorXd delta_lam   = delta.segment(0, p_lam);
        VectorXd delta_p     = delta.segment(p_lam, p_p);
        VectorXd delta_field = delta.segment(f_start, n_field);

        VectorXd beta_lam_try, beta_p_try, field_try;
        VectorXd eta_lam_try(n_sites), eta_p_try(n_obs);
        const double obj_cur = log_lik + field_log_prior(res.field);
        const bool stepped = newton_backtrack(
            obj_cur,
            [&](double step) {
                beta_lam_try = res.beta_lambda + step * delta_lam;
                beta_p_try   = res.beta_p      + step * delta_p;
                field_try    = res.field       + step * delta_field;

                eta_lambda_of(beta_lam_try, field_try, eta_lam_try);
                eta_p_try.noalias() = Xp * beta_p_try;
                double ll_try = count_kernel_log_lik_only_spatial(
                    obs_by_site, y_R, eta_lam_try, eta_p_try, K_max, r, site_fn
                );
                return ll_try + field_log_prior(field_try);
            },
            [&](double) {
                res.beta_lambda = beta_lam_try;
                res.beta_p      = beta_p_try;
                res.field       = field_try;
                center(res.beta_lambda, res.beta_p, res.field);
            });
        if (!stepped) {
            if (verbose) Rcpp::Rcout << "    (step halving exhausted)\n";
            break;
        }
        res.n_iter = iter + 1;
    }

    // ----- log marginal at the converged mode -----
    eta_lambda_of(res.beta_lambda, res.field, eta_lam);
    eta_p_long.noalias() = Xp * res.beta_p;
    double log_lik_final = count_kernel_sweep_spatial(
        obs_by_site, y_R, eta_lam, eta_p_long, K_max, r,
        grad_eta_lam, info_eta_lam, grad_eta_p, info_eta_p,
        mean_N, var_N, boundary_weight, score_wt_lambda, site_fn
    );
    double log_prior_final = field_log_prior(res.field);

    MatrixXd H_final = MatrixXd::Zero(n_x, n_x);
    add_obs_info(eta_p_long, info_eta_lam, info_eta_p, var_N, score_wt_lambda,
                 H_final);
    add_prior_H(H_final);
    nmix_add_diagonal_ridge(H_final);

    const int p_beta = p_lam + p_p;
    // The constrained block is the structured field: the whole CAR block, or
    // BYM2's v alone (w is proper and carries no flat direction).
    const bool constrain_field = fs.centers();
    Eigen::LLT<MatrixXd> chol(H_final);
    double log_det_H;
    if (chol.info() == Eigen::Success) {
        log_det_H = 2.0 * chol.matrixL().toDenseMatrix().diagonal()
                              .array().log().sum();
        res.cov_beta = nmix_constrained_top_cov(H_final, n_x, p_beta, f_start,
                                                n_spatial, constrain_field);
    } else {
        MatrixXd H_f = MatrixXd::Zero(n_x, n_x);
        add_complete_fisher(info_eta_lam, info_eta_p, H_f);
        add_prior_H(H_f);
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
        res.cov_beta = nmix_constrained_top_cov(H_f, n_x, p_beta, f_start,
                                                n_spatial, constrain_field);
    }

    res.log_lik = log_lik_final;
    res.log_marginal = log_lik_final + log_prior_final - 0.5 * log_det_H;
    res.grad_norm = grad_norm;
    res.boundary_max = boundary_weight.maxCoeff();
    return res;
}

// log|Q(rho)| = log|D - rho W| once per rho grid point (proper CAR). Returns
// -INFINITY if Q is not PD (treat as a tail indicator, skip the grid point).
inline double count_log_det_Q_car_proper(
    int n_spatial, double rho,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors
) {
    Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(n_spatial, n_spatial);
    for (int s = 0; s < n_spatial; ++s) {
        Q(s, s) = static_cast<double>(n_neighbors[s]);
        for (int kk = adj_row_ptr[s]; kk < adj_row_ptr[s + 1]; ++kk) {
            int t = adj_col_idx[kk];
            Q(s, t) = -rho;
        }
    }
    Eigen::LLT<Eigen::MatrixXd> chol(Q);
    if (chol.info() != Eigen::Success) return R_NegInf;
    return 2.0 * chol.matrixL().toDenseMatrix().diagonal()
                    .array().log().sum();
}

// ---------------------------------------------------------------------------
// Grid orchestration (shared across families via SiteFn). The N-mixture and
// removal nested-Laplace entry points are thin [[Rcpp::export]] wrappers around
// these; only the per-site marginal `site_fn` and the R-side K_max default
// differ (removal needs the per-site removal total, set by the R wrapper, so the
// orchestration takes K_max as given and only prep applies the max(y)+100 default
// when K_max < 0). Single source of truth for the outer-grid walk.
// ---------------------------------------------------------------------------

template <class SiteFn>
inline Rcpp::List run_count_nested_laplace_icar(
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
    std::string progress_file, SiteFn site_fn
) {
    NmixSpatialPrep pp = prep_nmix_spatial(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R, n_spatial,
        r_grid, beta_lambda_init, beta_p_init, K_max);
    const int p_lam = pp.p_lam, p_p = pp.p_p;
    Eigen::Map<Eigen::MatrixXd> Xl(REAL(X_lambda_R), pp.n_sites, p_lam);
    Eigen::Map<Eigen::MatrixXd> Xp(REAL(X_p_R), pp.n_obs, p_p);

    Eigen::VectorXd z_default(n_spatial);
    if (!z_init.empty()) {
        if ((int)z_init.size() != n_spatial) Rcpp::stop("length(z_init) must equal n_spatial.");
        for (int s = 0; s < n_spatial; ++s) z_default(s) = z_init[s];
    } else {
        z_default.setZero();
    }

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

    auto solve = [&](int k) -> NmixSpatialPoint {
        Eigen::VectorXd beta_lam = pp.beta_lam_default;
        Eigen::VectorXd beta_p   = pp.beta_p_default;
        Eigen::VectorXd z        = z_default;
        CountFieldSpec fs;
        fs.kind = CountFieldKind::CAR;
        fs.car  = CarPriorKind::ICAR;
        fs.tau  = tau_k[k];
        CountSpatialInnerResult ir = inner_newton_count_field(
            fs, p_lam, p_p, pp.n_sites, n_spatial, pp.n_obs,
            Xl, Xp, y, pp.obs_by_site, pp.map_site_to_unit,
            adj_row_ptr, adj_col_idx, n_neighbors, r_k[k],
            pp.K_max, max_iter, tol, beta_lam, beta_p, z, verbose, site_fn);
        NmixSpatialPoint pt;
        pt.log_marginal = ir.log_marginal; pt.n_iter = ir.n_iter;
        pt.converged = ir.converged;       pt.grad_norm = ir.grad_norm;
        pt.log_lik = ir.log_lik;           pt.boundary_max = ir.boundary_max;
        pt.coef = Eigen::VectorXd(p_lam + p_p);
        pt.coef.head(p_lam) = ir.beta_lambda;
        pt.coef.tail(p_p)   = ir.beta_p;
        pt.field = ir.field;
        pt.cov_beta = ir.cov_beta;
        return pt;
    };

    Rcpp::List out = run_nmix_spatial_grid(
        n_grid, p_lam, p_p, n_spatial, /*field_len=*/n_spatial, pp.K_max,
        theta_grid_out, solve,
        progress, progress_every, progress_throttle, progress_file);
    out["tau_grid"] = tau_grid;
    out["r_grid"]   = r_grid;
    out["n_tau"]    = n_tau;
    out["n_r"]      = n_r;
    return out;
}

template <class SiteFn>
inline Rcpp::List run_count_nested_laplace_bym2(
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
    std::string progress_file, SiteFn site_fn
) {
    if (scale_factor <= 0) Rcpp::stop("scale_factor must be positive.");
    NmixSpatialPrep pp = prep_nmix_spatial(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R, n_spatial,
        r_grid, beta_lambda_init, beta_p_init, K_max);
    const int p_lam = pp.p_lam, p_p = pp.p_p;
    Eigen::Map<Eigen::MatrixXd> Xl(REAL(X_lambda_R), pp.n_sites, p_lam);
    Eigen::Map<Eigen::MatrixXd> Xp(REAL(X_p_R), pp.n_obs, p_p);

    Eigen::VectorXd v_default(n_spatial), w_default(n_spatial);
    if (!v_init.empty()) {
        if ((int)v_init.size() != n_spatial) Rcpp::stop("length(v_init) must equal n_spatial.");
        for (int s = 0; s < n_spatial; ++s) v_default(s) = v_init[s];
    } else v_default.setZero();
    if (!w_init.empty()) {
        if ((int)w_init.size() != n_spatial) Rcpp::stop("length(w_init) must equal n_spatial.");
        for (int s = 0; s < n_spatial; ++s) w_default(s) = w_init[s];
    } else w_default.setZero();

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
                    theta_grid_out(k, 0) = sigma; theta_grid_out(k, 1) = rho;
                    theta_grid_out(k, 2) = r_grid[ir_disp];
                    valid_k[k] = (sigma > 0 && rho >= 0 && rho <= 1);
                    a_k[k] = valid_k[k] ? sigma * std::sqrt(rho / scale_factor) : 0.0;
                    b_k[k] = valid_k[k] ? sigma * std::sqrt(1.0 - rho) : 0.0;
                    r_k[k] = r_grid[ir_disp];
                }
    }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("sigma", "rho", "r");

    auto solve = [&](int k) -> NmixSpatialPoint {
        NmixSpatialPoint pt;
        if (!valid_k[k]) { pt.skipped = true; return pt; }
        Eigen::VectorXd beta_lam = pp.beta_lam_default;
        Eigen::VectorXd beta_p   = pp.beta_p_default;
        Eigen::VectorXd field0(2 * n_spatial);
        field0.head(n_spatial) = v_default;
        field0.tail(n_spatial) = w_default;
        CountFieldSpec fs;
        fs.kind = CountFieldKind::BYM2;
        fs.a = a_k[k];
        fs.b = b_k[k];
        CountSpatialInnerResult ir = inner_newton_count_field(
            fs, p_lam, p_p, pp.n_sites, n_spatial, pp.n_obs,
            Xl, Xp, y, pp.obs_by_site, pp.map_site_to_unit,
            adj_row_ptr, adj_col_idx, n_neighbors, r_k[k],
            pp.K_max, max_iter, tol, beta_lam, beta_p, field0, verbose, site_fn);
        pt.log_marginal = ir.log_marginal; pt.n_iter = ir.n_iter;
        pt.converged = ir.converged;       pt.grad_norm = ir.grad_norm;
        pt.log_lik = ir.log_lik;           pt.boundary_max = ir.boundary_max;
        pt.coef = Eigen::VectorXd(p_lam + p_p);
        pt.coef.head(p_lam) = ir.beta_lambda; pt.coef.tail(p_p) = ir.beta_p;
        pt.field = ir.field;
        pt.cov_beta = ir.cov_beta;
        return pt;
    };

    Rcpp::List out = run_nmix_spatial_grid(
        n_grid, p_lam, p_p, n_spatial, /*field_len=*/2 * n_spatial, pp.K_max,
        theta_grid_out, solve,
        progress, progress_every, progress_throttle, progress_file);
    out["sigma_grid"] = sigma_grid; out["rho_grid"] = rho_grid;
    out["r_grid"] = r_grid; out["scale_factor"] = scale_factor;
    out["n_sigma"] = n_sigma; out["n_rho"] = n_rho; out["n_r"] = n_r;
    return out;
}

template <class SiteFn>
inline Rcpp::List run_count_nested_laplace_car_proper(
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
    std::string progress_file, SiteFn site_fn
) {
    NmixSpatialPrep pp = prep_nmix_spatial(
        y, site_idx, map_site_to_unit_R, X_lambda_R, X_p_R, n_spatial,
        r_grid, beta_lambda_init, beta_p_init, K_max);
    const int p_lam = pp.p_lam, p_p = pp.p_p;
    Eigen::Map<Eigen::MatrixXd> Xl(REAL(X_lambda_R), pp.n_sites, p_lam);
    Eigen::Map<Eigen::MatrixXd> Xp(REAL(X_p_R), pp.n_obs, p_p);

    Eigen::VectorXd z_default(n_spatial);
    if (!z_init.empty()) {
        if ((int)z_init.size() != n_spatial) Rcpp::stop("length(z_init) must equal n_spatial.");
        for (int s = 0; s < n_spatial; ++s) z_default(s) = z_init[s];
    } else {
        z_default.setZero();
    }

    const int n_tau = tau_grid.size(), n_rho = rho_grid.size(), n_r = r_grid.size();
    const int n_grid = n_tau * n_rho * n_r;

    std::vector<double> log_det_Q_rho(n_rho);
    for (int ir_rho = 0; ir_rho < n_rho; ++ir_rho)
        log_det_Q_rho[ir_rho] = count_log_det_Q_car_proper(
            n_spatial, rho_grid[ir_rho], adj_row_ptr, adj_col_idx, n_neighbors);

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

    auto solve = [&](int k) -> NmixSpatialPoint {
        NmixSpatialPoint pt;
        if (!R_finite(logdet_k[k])) { pt.skipped = true; return pt; }
        Eigen::VectorXd beta_lam = pp.beta_lam_default;
        Eigen::VectorXd beta_p   = pp.beta_p_default;
        Eigen::VectorXd z        = z_default;
        CountFieldSpec fs;
        fs.kind = CountFieldKind::CAR;
        fs.car  = CarPriorKind::CAR_PROPER;
        fs.tau  = tau_k[k];
        fs.rho  = rho_k[k];
        fs.log_det_Q_rho = logdet_k[k];
        CountSpatialInnerResult ir = inner_newton_count_field(
            fs, p_lam, p_p, pp.n_sites, n_spatial, pp.n_obs,
            Xl, Xp, y, pp.obs_by_site, pp.map_site_to_unit,
            adj_row_ptr, adj_col_idx, n_neighbors, r_k[k],
            pp.K_max, max_iter, tol, beta_lam, beta_p, z, verbose, site_fn);
        pt.log_marginal = ir.log_marginal; pt.n_iter = ir.n_iter;
        pt.converged = ir.converged;       pt.grad_norm = ir.grad_norm;
        pt.log_lik = ir.log_lik;           pt.boundary_max = ir.boundary_max;
        pt.coef = Eigen::VectorXd(p_lam + p_p);
        pt.coef.head(p_lam) = ir.beta_lambda;
        pt.coef.tail(p_p)   = ir.beta_p;
        pt.field = ir.field;
        pt.cov_beta = ir.cov_beta;
        return pt;
    };

    Rcpp::List out = run_nmix_spatial_grid(
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

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_COUNT_SPATIAL_DRIVER_H
