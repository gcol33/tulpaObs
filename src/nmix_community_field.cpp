// nmix_community_field.cpp
// Shared-field mode-find + Laplace marginal for the SPATIAL community
// N-mixture. Given the per-species abundance / detection coefficients (the
// community means + RE deviations, mu + b_s) held fixed, this finds the mode of
// the SHARED areal field z (length n_spatial, one unit per site) on the
// abundance arm and returns the field Laplace log-marginal contribution
//
//   log p(y | z_mode, coefs) - already in the community marginal; HERE we add
//   the field-block correction   log p(z_mode | tau) - 0.5 log|H_zz(z_mode)|
//   plus the part of the data log-lik that depends on z.
//
// The field couples all species through z[u(i)] = z_i: each species contributes
// its per-site abundance score / observed-information to the SAME field
// coordinate, so the field Newton aggregates over species. The per-site
// marginal kernel is the shared nmix_kernel.h; the ICAR / CAR-proper prior
// contributions reuse nmix_spatial_kernel.h. This is the field half of the
// nested-approx + debias split (the community RE block is integrated by the
// engine's tulpa_re_aghq given this z offset; here z and its hyperparameter are
// the outer-integrated latent Gaussian).
//
// Poisson abundance only. One spatial unit per site (map_site_to_unit = identity
// on the caller side; we still take an explicit map for generality).

#include "nmix_kernel.h"
#include "nmix_spatial_kernel.h"
#include "newton_step.h"                // newton_backtrack
#include "tobs_math.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using tulpaObs::clamp_eta;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using tulpaObs::compute_nmix_site;
using tulpaObs::NMixSiteResult;
using tulpaObs::nmix_car_quadratic_form;
using tulpaObs::newton_backtrack;
using tulpaObs::kFieldMaxHalvings;

namespace {

// Per-(species, site) row grouping: for each species, the long-form row indices
// at each site, in input order.
struct SiteVisits {
    int site;                 // 0-based
    std::vector<int> rows;    // long-form indices for this (species, site)
};

}  // namespace

// Inner field-mode Newton + Laplace marginal at a fixed ICAR/CAR precision tau
// (rho = 1 for ICAR; rho < 1 for proper CAR with precomputed log|Q(rho)|),
// holding the per-species coefficients fixed.
//
// Inputs:
//   y, site_idx, species_idx : long form (1-based site / species).
//   X_lambda  : n_sites x p_lam (shared design).
//   X_p       : n_obs   x p_p   (long-form detection design).
//   coef_lambda : n_species x p_lam (mu_lambda + b_lambda_s per species).
//   coef_p      : n_species x p_p   (mu_p     + b_p_s     per species).
//   adj_*       : CSR adjacency on the n_spatial units.
//   tau, rho, log_det_Q_rho : field hyperparameter (rho = 1, log_det = 0 -> ICAR).
//   z_init      : warm start (length n_spatial).
//
// Returns: z_mode, the z-dependent data log-lik, the field log-prior, the field
// Laplace marginal contribution (log_prior - 0.5 log|H_zz|), grad_norm,
// boundary_max, converged.
// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_field_solve(
    Rcpp::IntegerVector y,
    Rcpp::IntegerVector site_idx,
    Rcpp::IntegerVector species_idx,
    Rcpp::NumericMatrix X_lambda,
    Rcpp::NumericMatrix X_p,
    Rcpp::NumericMatrix coef_lambda,
    Rcpp::NumericMatrix coef_p,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    double tau,
    double rho,
    double log_det_Q_rho,
    Rcpp::NumericVector z_init,
    int K_max,
    int max_iter = 100,
    double tol = 1e-6,
    bool verbose = false
) {
    const int n_sites   = X_lambda.nrow();
    const int p_lam     = X_lambda.ncol();
    const int n_obs     = X_p.nrow();
    const int p_p       = X_p.ncol();
    const int n_species = coef_lambda.nrow();
    const bool is_icar  = (rho >= 1.0 - 1e-12);

    if ((int)y.size() != n_obs) Rcpp::stop("length(y) must equal nrow(X_p).");
    if ((int)site_idx.size() != n_obs) Rcpp::stop("length(site_idx) must equal nrow(X_p).");
    if ((int)species_idx.size() != n_obs) Rcpp::stop("length(species_idx) must equal nrow(X_p).");
    if (coef_lambda.ncol() != p_lam) Rcpp::stop("ncol(coef_lambda) must equal p_lam.");
    if (coef_p.ncol() != p_p) Rcpp::stop("ncol(coef_p) must equal p_p.");
    if ((int)map_site_to_unit_R.size() != n_sites) Rcpp::stop("map_site_to_unit length mismatch.");

    std::vector<int> map_site_to_unit(n_sites);
    for (int s = 0; s < n_sites; ++s) {
        int u = map_site_to_unit_R[s] - 1;
        if (u < 0 || u >= n_spatial) Rcpp::stop("map_site_to_unit out of range.");
        map_site_to_unit[s] = u;
    }

    // Group rows by (species, site).
    std::vector<std::vector<SiteVisits>> sp_sites(n_species);
    {
        std::vector<std::vector<std::vector<int>>> rows(
            n_species, std::vector<std::vector<int>>(n_sites));
        for (int r = 0; r < n_obs; ++r)
            rows[species_idx[r] - 1][site_idx[r] - 1].push_back(r);
        for (int s = 0; s < n_species; ++s)
            for (int i = 0; i < n_sites; ++i)
                if (!rows[s][i].empty()) {
                    SiteVisits sv; sv.site = i; sv.rows = rows[s][i];
                    sp_sites[s].push_back(std::move(sv));
                }
    }

    Map<MatrixXd> Xl(REAL(X_lambda), n_sites, p_lam);
    Map<MatrixXd> Xp(REAL(X_p), n_obs, p_p);
    Map<MatrixXd> Cl(REAL(coef_lambda), n_species, p_lam);
    Map<MatrixXd> Cp(REAL(coef_p), n_species, p_p);

    // Site-major baseline abundance predictor eta0_{s,i} = X_lambda_i . coef_lambda_s
    // (NO field), and per-(species,site) detection predictors. Precompute the
    // field-independent parts so the Newton only re-evaluates the marginal as z
    // moves.
    // For each (species, site) entry we keep eta0 and the eta_p vector.
    struct Cell {
        int site;
        std::vector<int> y_site;
        std::vector<double> eta_p;
        double eta0_lambda;
    };
    std::vector<std::vector<Cell>> cells(n_species);
    for (int s = 0; s < n_species; ++s) {
        cells[s].reserve(sp_sites[s].size());
        for (const SiteVisits& sv : sp_sites[s]) {
            Cell c; c.site = sv.site;
            double e0 = 0.0;
            for (int k = 0; k < p_lam; ++k) e0 += Xl(sv.site, k) * Cl(s, k);
            c.eta0_lambda = e0;
            const int J = (int)sv.rows.size();
            c.y_site.resize(J); c.eta_p.resize(J);
            for (int j = 0; j < J; ++j) {
                const int r = sv.rows[j];
                c.y_site[j] = y[r];
                double ep = 0.0;
                for (int k = 0; k < p_p; ++k) ep += Xp(r, k) * Cp(s, k);
                c.eta_p[j] = clamp_eta(ep);
            }
            cells[s].push_back(std::move(c));
        }
    }

    VectorXd z(n_spatial);
    if ((int)z_init.size() == n_spatial)
        for (int u = 0; u < n_spatial; ++u) z(u) = z_init[u];
    else z.setZero();

    auto z_log_prior = [&](const VectorXd& zz) -> double {
        double quad = nmix_car_quadratic_form(
            n_spatial, rho, adj_row_ptr, adj_col_idx, n_neighbors, zz);
        if (is_icar)
            return -0.5 * tau * quad + 0.5 * (n_spatial - 1) * std::log(tau);
        return 0.5 * log_det_Q_rho + 0.5 * n_spatial * std::log(tau)
               - 0.5 * tau * quad;
    };

    // One sweep: at the current z, accumulate the z-dependent data log-lik, the
    // field gradient and field-block observed information (with the Var[N|y]
    // rank-1 z-z correction), and the boundary weight.
    auto sweep = [&](const VectorXd& zz,
                     VectorXd& grad_z, MatrixXd& Hzz,
                     double& boundary_max, bool want_hess) -> double {
        double log_lik = 0.0;
        boundary_max = 0.0;
        grad_z.setZero();
        if (want_hess) Hzz.setZero();
        for (int s = 0; s < n_species; ++s) {
            for (const Cell& c : cells[s]) {
                const int u = map_site_to_unit[c.site];
                const int J = (int)c.y_site.size();
                double eta_lam = c.eta0_lambda + zz(u);
                eta_lam = clamp_eta(eta_lam);
                NMixSiteResult res = compute_nmix_site(
                    c.y_site.data(), c.eta_p.data(), J, eta_lam, K_max);
                log_lik += res.log_lik;
                grad_z(u) += res.grad_eta_lambda;
                if (res.boundary_weight > boundary_max)
                    boundary_max = res.boundary_weight;
                if (want_hess) {
                    // z-z block of the marginal observed info at this site:
                    //   info_eta_lambda - var_N * score_wt_lambda^2
                    // (the loading of z onto eta_lambda has coefficient 1).
                    const double swl = res.score_wt_lambda;
                    Hzz(u, u) += res.info_eta_lambda - res.var_N * swl * swl;
                }
            }
        }
        return log_lik;
    };

    VectorXd grad_z(n_spatial);
    MatrixXd Hzz(n_spatial, n_spatial);
    double boundary_max = 0.0;
    double grad_norm = R_PosInf;
    double log_lik = R_NegInf;
    bool converged = false;
    int n_iter = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        log_lik = sweep(z, grad_z, Hzz, boundary_max, /*want_hess=*/true);

        // Add the prior: grad -= tau Q z; H += tau Q.
        VectorXd grad = grad_z;
        MatrixXd H = Hzz;
        for (int sp = 0; sp < n_spatial; ++sp) {
            const double q_diag = (double)n_neighbors[sp];
            H(sp, sp) += tau * q_diag;
            double nbr = 0.0;
            for (int kk = adj_row_ptr[sp]; kk < adj_row_ptr[sp + 1]; ++kk) {
                int t = adj_col_idx[kk];
                nbr += z(t);
                if (t > sp) {
                    H(sp, t) -= tau * rho;
                    H(t, sp) -= tau * rho;
                }
            }
            grad(sp) -= tau * (q_diag * z(sp) - rho * nbr);
        }

        grad_norm = grad.norm();
        if (verbose)
            Rcpp::Rcout << "  field iter " << iter << " ll " << log_lik
                        << " |g| " << grad_norm << "\n";
        if (grad_norm < tol) { converged = true; n_iter = iter + 1; break; }

        // Ridge for the (intercept, constant-z) null direction under ICAR.
        double md = 0.0;
        for (int sp = 0; sp < n_spatial; ++sp) md += H(sp, sp);
        md /= std::max(1, n_spatial);
        double ridge = std::max(1e-10 * md, 1e-12);
        for (int sp = 0; sp < n_spatial; ++sp) H(sp, sp) += ridge;

        Eigen::LLT<MatrixXd> chol(H);
        if (chol.info() != Eigen::Success) {
            if (verbose) Rcpp::Rcout << "  (field Cholesky failed)\n";
            break;
        }
        VectorXd delta = chol.solve(grad);

        // Step halving on the field objective (z-dependent ll + field prior).
        const double obj_cur = log_lik + z_log_prior(z);
        VectorXd z_try(n_spatial), g_dummy(n_spatial);
        MatrixXd H_dummy(0, 0);
        double bd_dummy;
        const bool stepped = newton_backtrack(
            obj_cur,
            [&](double step) {
                z_try = z + step * delta;
                if (is_icar) {       // sum-to-zero centering
                    double mean = z_try.mean();
                    z_try.array() -= mean;
                }
                const double ll_try = sweep(z_try, g_dummy, H_dummy, bd_dummy,
                                            false);
                return ll_try + z_log_prior(z_try);
            },
            [&](double) { z = z_try; },
            kFieldMaxHalvings);
        if (!stepped) { if (verbose) Rcpp::Rcout << "  (field step halving exhausted)\n"; break; }
        n_iter = iter + 1;
    }

    // Final marginal-correction terms at the mode.
    log_lik = sweep(z, grad_z, Hzz, boundary_max, /*want_hess=*/true);
    MatrixXd H = Hzz;
    for (int sp = 0; sp < n_spatial; ++sp) {
        H(sp, sp) += tau * (double)n_neighbors[sp];
        for (int kk = adj_row_ptr[sp]; kk < adj_row_ptr[sp + 1]; ++kk) {
            int t = adj_col_idx[kk];
            if (t > sp) { H(sp, t) -= tau * rho; H(t, sp) -= tau * rho; }
        }
    }
    // For ICAR, the rank-deficient (constant-z) direction is pinned by the
    // sum-to-zero constraint; add a large quadratic penalty on (sum z)^2 so the
    // log|H_zz| is the constrained determinant (the field Laplace normaliser).
    if (is_icar) {
        double md = 0.0;
        for (int sp = 0; sp < n_spatial; ++sp) md += H(sp, sp);
        md = std::abs(md) / std::max(1, n_spatial);
        if (!(md > 0)) md = 1.0;
        const double kappa = 1e6 * md;
        for (int i = 0; i < n_spatial; ++i)
            for (int j = 0; j < n_spatial; ++j)
                H(i, j) += kappa;
    } else {
        double md = 0.0;
        for (int sp = 0; sp < n_spatial; ++sp) md += H(sp, sp);
        md /= std::max(1, n_spatial);
        double ridge = std::max(1e-10 * md, 1e-12);
        for (int sp = 0; sp < n_spatial; ++sp) H(sp, sp) += ridge;
    }

    double log_prior_z = z_log_prior(z);
    double log_det_Hzz = R_NaN;
    double field_marginal;
    Eigen::LLT<MatrixXd> chol(H);
    if (chol.info() == Eigen::Success) {
        log_det_Hzz = 2.0 * chol.matrixL().toDenseMatrix().diagonal()
                              .array().log().sum();
        // Field Laplace marginal correction: log p(z|tau) - 0.5 log|H_zz|.
        // (The z-dependent data log-lik is already inside the community marginal
        // the R driver computes at this offset; we return it for diagnostics.)
        field_marginal = log_prior_z - 0.5 * log_det_Hzz;
    } else {
        // Singular field Hessian: reject this grid cell (the R driver drops any
        // non-finite field_marginal) rather than propagate a NaN.
        field_marginal = -std::numeric_limits<double>::infinity();
    }

    Rcpp::NumericVector z_out(n_spatial);
    for (int u = 0; u < n_spatial; ++u) z_out[u] = z(u);

    return Rcpp::List::create(
        Rcpp::Named("z")              = z_out,
        Rcpp::Named("log_lik_z")      = log_lik,
        Rcpp::Named("log_prior_z")    = log_prior_z,
        Rcpp::Named("log_det_Hzz")    = log_det_Hzz,
        Rcpp::Named("field_marginal") = field_marginal,
        Rcpp::Named("grad_norm")      = grad_norm,
        Rcpp::Named("boundary_max")   = boundary_max,
        Rcpp::Named("converged")      = converged,
        Rcpp::Named("n_iter")         = n_iter
    );
}
