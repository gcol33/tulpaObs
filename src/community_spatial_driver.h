// community_spatial_driver.h
//
// STATUS: NOT YET WIRED. This is the extraction target for the community
// spatial Laplace-EM written twice, once per response family.
// Nothing includes it yet, so it contributes no code to the build; migrating
// src/nmix_community_spatial.cpp and src/ms_occu_spatial.cpp onto it, and then
// giving the occupancy side the SPDE arm the N-mixture side already has, is the
// remaining work on that issue.
//
// It has been AUDITED against both twins rather than assumed faithful, since it
// was written before any of it could compile:
//
//   * Normalized diff of `community_spatial_em` below against
//     nmix_community_spatial.cpp:322-631 and ms_occu_spatial.cpp:214-435 --
//     comments, blanks and whitespace stripped, the arm-specific names mapped
//     onto each other. 82% of lines identical to the N-mixture twin, 71% to the
//     occupancy one, and every one of the 18 residual hunks is the Family
//     indirection, a rename, or `(double) S` against `(double)S`. No arithmetic
//     differs.
//
//   * The two places the twins genuinely disagree both resolve in favour of
//     what is written here, and neither changes a result:
//       - `add_field_prior` takes the grid `rho` unconditionally on the
//         N-mixture side and forces `rho = 1` for ICAR on the occupancy side.
//         ICAR is intrinsic, so `rho = 1` is the defined value; the form kept
//         here is the occupancy one. It is a no-op today because the N-mixture
//         ICAR plan pins `rho = 1.0` at nmix_community_spatial.cpp:765, and it
//         removes the hazard of a future ICAR grid that varies rho.
//       - `sigma_beta` is threaded through four exported N-mixture entries and
//         then discarded by `(void) sigma_beta;` at
//         nmix_community_spatial.cpp:353, so dropping it from this signature
//         preserves behaviour. It remains an R-visible argument that does
//         nothing, which is a separate cleanup from this extraction.
//
// What the audit does NOT establish is that a migrated build reproduces the
// current one. That needs the migration plus a before/after equivalence run on
// real fits, and until it has been done neither twin should be deleted.
//
// Nested Laplace-EM for a community (multi-species) model with Gaussian
// community hyperpriors on the per-species coefficients AND one field shared
// across species, loading onto the state arm:
//
//   eta_state[s, i] = X_state_i . (mu_state + b_state_s) + f_{u(i)}
//   eta_det  [s, i] = X_det_i   . (mu_det   + b_det_s)
//   b_state_s ~ N(0, Sigma_state),  b_det_s ~ N(0, Sigma_p),  one f shared
//   f ~ ICAR(tau) | BYM2(sigma, rho) | CAR(tau, rho) | Matern SPDE(Q)
//
// Three latent layers, integrated in three ways:
//   * the per-species response -- summed analytically by the family's per-site
//     marginal (a Royle N-mixture site, an occupancy two-state cell, ...);
//   * b_s -- per-species Laplace (Gaussian community prior); the species are
//     conditionally independent given (mu, f), so they fold out by a Schur
//     complement (Louis 1982), exactly as the non-spatial community EM;
//   * (mu, f) -- joint Laplace; f has the field prior, mu is flat. mu and f are
//     NOT conditionally independent (both load eta_state), so they are ONE block.
//
// Per outer grid point (tau[, rho] / sigma, rho / Matern Q [, NB size r]):
//   EM (Sigma fixed within, updated across):
//     mode-find  -- joint Newton over (t, {b_s}) by block elimination on the
//                   complete-data-Fisher (PSD) arrowhead Hessian; step halving
//                   on the joint objective; centre f sum-to-zero (ICAR / BYM2 v).
//     M-step     -- Sigma_k = mean_s [ b_s b_s' + Cov(b_s | y) ],
//                   Cov(b_s | y) = (Fisher_s + Sigma^{-1})^{-1}.
//   final pass   -- observed-info arrowhead; the b-folded (mu, f) info M; the
//                   grid-point Laplace log-marginal and the community-mean
//                   covariance vcov_mu = top-left block of M^{-1} (constrained).
//
// The top block t = (mu_state, mu_det, f) has the SAME layout as the single-
// species spatial state x = (beta_lambda, beta_p, field), so the field prior /
// centering / constrained-covariance helpers in nmix_spatial_kernel.h and
// nmix_spatial_kernel_bym2.h
// apply to t directly (field_start = p_state + p_p).
//
// The RESPONSE FAMILY enters only through the per-species site block. It
// supplies:
//
//   int    n_species() const
//   int    n_records(int s) const        // site records held for species s
//   int    site_of(int s, int i) const   // site index of record i, for the field
//   double site_block(int s, int i, const Eigen::VectorXd& coef,
//                     double field_offset, const double* load, int nfield,
//                     bool want_block, bool want_obs,
//                     Eigen::VectorXd& grad_aug, Eigen::MatrixXd& small,
//                     double* boundary_out) const
//
// `site_block` fills grad_aug (length d + nfield, the d coefficient coordinates
// first and the nfield field loadings last) and, when want_block, the symmetric
// (d + nfield) curvature block; `want_obs` selects the observed information over
// the complete-data Fisher. It returns the record's marginal log-likelihood, and
// writes its boundary diagnostic to `boundary_out` when the family carries one.
// A record with no data contributes nothing -- the field and the community prior
// interpolate it -- which is the single-species held-out convention.
//
// Adding a field structure is a new branch on FieldKind plus a plan builder;
// adding a response family is a new site-block type.
//
// Include AFTER RcppEigen.h: Rcpp::wrap on an Eigen matrix is looked up when
// community_grid_pack.h
// is parsed.

#ifndef TULPAOBS_COMMUNITY_SPATIAL_DRIVER_H
#define TULPAOBS_COMMUNITY_SPATIAL_DRIVER_H

#include "nmix_linalg.h"
#include "nmix_progress.h"
#include "nmix_spatial_kernel.h"        // icar / car_proper prior, Q-adds, centering
#include "nmix_spatial_kernel_bym2.h"   // bym2 prior, Q+I-adds, centering
#include "newton_step.h"                // newton_backtrack
#include "community_grid_pack.h"        // community_pack_grid
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace tulpaObs {

enum class FieldKind { ICAR, CAR_PROPER, BYM2, SPDE };

// SPDE field context: the continuous Matern field lives at n_mesh FEM nodes; the
// dense projection A (n_sites x n_mesh) maps mesh nodes onto sites, and Q
// (n_mesh x n_mesh, carrying tau_spde^2) is the proper Matern precision built
// once per grid point on the R side. The areal kinds carry a single field unit
// per site; SPDE carries the whole mesh row, so the field loading and prior are
// threaded through this struct rather than the (u, a, b) areal triple.
struct SpdeCtx {
    const Eigen::MatrixXd* A = nullptr;   // n_sites x n_mesh
    const Eigen::MatrixXd* Q = nullptr;   // n_mesh x n_mesh (with tau_spde^2)
    double log_det_Q = 0.0;
    int n_mesh = 0;
};

// Field columns one site loads onto eta_state. Areal kinds load a single field
// unit (BYM2 two: v + w); the SPDE field loads the mesh nodes of the site's
// projection row, and a 2D P1 FEM row has at most this many nonzeros -- the
// barycentric weights of the triangle the site falls in.
static const int kMaxFieldLoad = 4;
struct FieldGeom {
    int nfield;
    double load[kMaxFieldLoad];
    int    idx[kMaxFieldLoad];
};

// Areal field geometry for spatial unit u; field_start = d.
inline FieldGeom field_geom(FieldKind kind, int d, int n_spatial, int u,
                            double a, double b) {
    FieldGeom g;
    if (kind == FieldKind::BYM2) {
        g.nfield = 2;
        g.load[0] = a; g.idx[0] = d + u;                 // v
        g.load[1] = b; g.idx[1] = d + n_spatial + u;     // w
    } else {
        g.nfield = 1;
        g.load[0] = 1.0; g.idx[0] = d + u;               // f
    }
    return g;
}

// SPDE field geometry for one site: the nonzeros of A.row(site).
inline FieldGeom field_geom_spde(int d, const SpdeCtx& sp, int site) {
    FieldGeom g;
    g.nfield = 0;
    for (int k = 0; k < sp.n_mesh; ++k) {
        const double w = (*sp.A)(site, k);
        if (w != 0.0) {
            if (g.nfield >= kMaxFieldLoad)
                Rcpp::stop("SPDE projection row has more than %d nonzeros; "
                           "raise kMaxFieldLoad.", kMaxFieldLoad);
            g.load[g.nfield] = w;
            g.idx[g.nfield]  = d + k;
            ++g.nfield;
        }
    }
    return g;
}

inline double field_offset(FieldKind kind, const Eigen::VectorXd& field,
                           int n_spatial, int u, double a, double b) {
    if (kind == FieldKind::BYM2) return a * field(u) + b * field(n_spatial + u);
    return field(u);
}

// SPDE per-site offset: (A u)_site = sum_k A[site, k] u[k].
inline double field_offset_spde(const Eigen::VectorXd& field, const SpdeCtx& sp,
                                int site) {
    double off = 0.0;
    for (int k = 0; k < sp.n_mesh; ++k) {
        const double w = (*sp.A)(site, k);
        if (w != 0.0) off += w * field(k);
    }
    return off;
}

// log p(field | hyper), dropping grid-independent constants. BYM2's (v, w)
// priors are unit-scale (the amplitude is in the loadings), so its value is
// grid-independent -- kept for the line-search objective only.
inline double field_log_prior(FieldKind kind, int n_spatial,
                              double tau, double rho, double log_det_Q_rho,
                              const Rcpp::IntegerVector& adj_row_ptr,
                              const Rcpp::IntegerVector& adj_col_idx,
                              const Rcpp::IntegerVector& n_neighbors,
                              const Eigen::VectorXd& field,
                              const SpdeCtx* spde = nullptr) {
    if (kind == FieldKind::SPDE) {
        // 0.5 log|Q| - 0.5 u' Q u (full rank; the (2 pi)^{-n_mesh/2} constant is
        // grid-independent and dropped, as in the single-species SPDE path).
        return 0.5 * spde->log_det_Q
             - 0.5 * (field.transpose() * ((*spde->Q) * field))(0, 0);
    }
    if (kind == FieldKind::BYM2) {
        const Eigen::VectorXd v = field.head(n_spatial);
        const Eigen::VectorXd w = field.tail(n_spatial);
        return nmix_bym2_log_prior(n_spatial, adj_row_ptr, adj_col_idx,
                                   n_neighbors, v, w);
    }
    if (kind == FieldKind::ICAR) {
        return nmix_icar_log_prior(n_spatial, tau, adj_row_ptr, adj_col_idx,
                                   n_neighbors, field);
    }
    return nmix_car_proper_log_prior(n_spatial, tau, rho, log_det_Q_rho,
                                     adj_row_ptr, adj_col_idx, n_neighbors, field);
}

// Add the field prior to the top-block gradient and Hessian (in-place). The top
// block has the (mu_state, mu_det, field) layout the single-species helpers
// expect (field_start = d = p_state + p_p). ICAR is the intrinsic form, so its
// CAR precision is read at rho = 1 whatever the grid point carries; proper-CAR
// takes the grid rho; BYM2 splits the (v, w) blocks.
inline void add_field_prior(FieldKind kind, int p_state, int p_p, int n_spatial,
                            double tau, double rho,
                            const Rcpp::IntegerVector& adj_row_ptr,
                            const Rcpp::IntegerVector& adj_col_idx,
                            const Rcpp::IntegerVector& n_neighbors,
                            const Eigen::VectorXd& field,
                            Eigen::VectorXd& grad_top, Eigen::MatrixXd& T,
                            const SpdeCtx* spde = nullptr) {
    if (kind == FieldKind::SPDE) {
        // Field prior -Q u to the score, +Q to the negative-Hessian field block.
        const int d = p_state + p_p;
        grad_top.segment(d, spde->n_mesh).noalias() -= (*spde->Q) * field;
        T.block(d, d, spde->n_mesh, spde->n_mesh).noalias() += (*spde->Q);
        return;
    }
    if (kind == FieldKind::BYM2) {
        const Eigen::VectorXd v = field.head(n_spatial);
        const Eigen::VectorXd w = field.tail(n_spatial);
        nmix_add_bym2_prior_to_grad_and_H(p_state, p_p, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors, v, w, grad_top, T);
        return;
    }
    const double rho_use = (kind == FieldKind::CAR_PROPER) ? rho : 1.0;
    nmix_add_car_to_spatial_block(p_state, p_p, n_spatial, tau, rho_use,
        adj_row_ptr, adj_col_idx, n_neighbors, field, grad_top, T);
}

// Sum-to-zero centering of the (intrinsic) field component against the global
// state-arm intercept. ICAR centres f; BYM2 centres v only (w is proper).
// CAR_proper and SPDE are full rank: the (intercept, field-mean) direction is
// identified by Q itself, so they need none.
inline void center_field(FieldKind kind, int p_state, int p_p, int n_spatial,
                         Eigen::VectorXd& field) {
    if (kind == FieldKind::CAR_PROPER || kind == FieldKind::SPDE ||
        n_spatial <= 0) return;
    const int len = (kind == FieldKind::BYM2) ? 2 * n_spatial : n_spatial;
    Eigen::VectorXd holder(p_state + p_p + len);
    holder.setZero();
    holder.segment(p_state + p_p, len) = field;
    nmix_center_field(p_state, p_p, n_spatial, holder);
    field = holder.segment(p_state + p_p, len);
}

struct CommSpatialResult {
    Eigen::VectorXd mu;      // d community means
    Eigen::VectorXd field;   // nfield * n_spatial (ICAR/CAR/SPDE: f; BYM2: [v; w])
    Eigen::MatrixXd Sigma_state, Sigma_p;
    Eigen::MatrixXd blup_state, blup_p;   // S x p_state, S x p_p
    Eigen::MatrixXd vcov_mu;              // d x d, b- and field-folded
    double log_marginal = R_NegInf;
    double log_lik = R_NegInf;
    bool   converged = false;
    int    n_iter = 0;
    double boundary_max = 0.0;
};

// One outer-grid-point Laplace-EM fit. `field` is laid out [f] (ICAR / CAR /
// SPDE) or [v; w] (BYM2); `a`, `b` are the BYM2 loadings.
template <class Family>
CommSpatialResult community_spatial_em(
    FieldKind kind,
    const Family& fam,
    int p_state, int p_p,
    const std::vector<int>& map_site_to_unit,
    int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    double tau, double rho, double log_det_Q_rho, double a, double b,
    const Eigen::VectorXd& mu_init, const Eigen::MatrixXd& Sigma_state_init,
    const Eigen::MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose,
    const SpdeCtx* spde = nullptr) {

    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int d      = p_state + p_p;
    const int S      = fam.n_species();
    const int nfield = (kind == FieldKind::BYM2) ? 2 : 1;
    const int field_len = nfield * n_spatial;
    const int m = d + field_len;
    // No ridge on the community means (mu). The state-arm intercept and the
    // field constant mode form an exactly flat (data + intrinsic-prior)
    // direction; any ridge on mu is the ONLY curvature along it and drives the
    // intercept to its prior mean 0, with the field absorbing the level (then
    // deleted by the sum-to-zero centering). The single-species spatial path
    // keeps the fixed effects flat for the same reason -- mu is identified by
    // the data once the field is anchored sum-to-zero.

    CommSpatialResult out;
    VectorXd mu = mu_init;
    VectorXd field = VectorXd::Zero(field_len);
    std::vector<VectorXd> bvec(S, VectorXd::Zero(d));
    MatrixXd Sig_state = Sigma_state_init, Sig_p = Sigma_p_init;

    auto block_prec = [&](const MatrixXd& Ss, const MatrixXd& Sp) {
        MatrixXd P = MatrixXd::Zero(d, d);
        P.topLeftCorner(p_state, p_state) = nmix_safe_inverse(Ss);
        P.bottomRightCorner(p_p, p_p)     = nmix_safe_inverse(Sp);
        return P;
    };

    // Per-site (geometry, offset) -- branches on the field kind once. SPDE reads
    // the site's projection row; the areal kinds use the single-unit map.
    auto site_geom = [&](int site, const VectorXd& field_) {
        if (kind == FieldKind::SPDE)
            return std::make_pair(field_geom_spde(d, *spde, site),
                                  field_offset_spde(field_, *spde, site));
        const int u = map_site_to_unit[site];
        return std::make_pair(field_geom(kind, d, n_spatial, u, a, b),
                              field_offset(kind, field_, n_spatial, u, a, b));
    };

    // Data log-lik at the current state (line-search objective helper).
    auto data_loglik = [&](const VectorXd& mu_, const VectorXd& field_,
                           const std::vector<VectorXd>& b_) {
        double ll = 0.0;
        VectorXd grad_aug; MatrixXd small;
        for (int s = 0; s < S; ++s) {
            VectorXd coef = mu_ + b_[s];
            const int n_rec = fam.n_records(s);
            for (int i = 0; i < n_rec; ++i) {
                auto go = site_geom(fam.site_of(s, i), field_);
                const FieldGeom& g = go.first;
                double l = fam.site_block(s, i, coef, go.second,
                                          g.load, g.nfield, false, false,
                                          grad_aug, small, nullptr);
                if (!R_finite(l)) return R_NegInf;
                ll += l;
            }
        }
        return ll;
    };

    // Joint objective = data loglik + field prior + b prior (mu is flat).
    auto objective = [&](const VectorXd& mu_, const VectorXd& field_,
                         const std::vector<VectorXd>& b_, const MatrixXd& P) {
        double obj = data_loglik(mu_, field_, b_);
        if (!R_finite(obj)) return R_NegInf;
        obj += field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                               adj_row_ptr, adj_col_idx, n_neighbors, field_, spde);
        for (int s = 0; s < S; ++s) obj -= 0.5 * b_[s].dot(P * b_[s]);
        return obj;
    };

    // Assemble the arrowhead system at the current state. Fills grad_top (m),
    // T (m x m top block), and per-species grad_b / D / C. `want_obs` selects
    // the observed-info (true) vs complete-data Fisher (false) curvature.
    auto assemble = [&](const VectorXd& mu_, const VectorXd& field_,
                        const std::vector<VectorXd>& b_, const MatrixXd& P,
                        bool want_obs,
                        VectorXd& grad_top, MatrixXd& T,
                        std::vector<VectorXd>& grad_b,
                        std::vector<MatrixXd>& D,
                        std::vector<MatrixXd>& C,
                        double* log_lik_out, double* boundary_out) {
        grad_top.setZero(m);
        T.setZero(m, m);
        double ll = 0.0, bmax = 0.0;
        VectorXd grad_aug; MatrixXd small;
        for (int s = 0; s < S; ++s) {
            grad_b[s].setZero(d);
            D[s].setZero(d, d);
            C[s].setZero(m, d);
            VectorXd coef = mu_ + b_[s];
            const int n_rec = fam.n_records(s);
            for (int i = 0; i < n_rec; ++i) {
                auto go = site_geom(fam.site_of(s, i), field_);
                const FieldGeom& g = go.first;
                double bw = 0.0;
                ll += fam.site_block(s, i, coef, go.second,
                                     g.load, g.nfield, true, want_obs,
                                     grad_aug, small, &bw);
                if (bw > bmax) bmax = bw;
                // d-block: per-species coefficient curvature + community mean.
                grad_b[s].head(d) += grad_aug.head(d);
                grad_top.head(d)  += grad_aug.head(d);
                D[s].noalias()                  += small.topLeftCorner(d, d);
                T.topLeftCorner(d, d).noalias() += small.topLeftCorner(d, d);
                C[s].topLeftCorner(d, d).noalias() += small.topLeftCorner(d, d);
                // field rows/cols and crosses.
                for (int f = 0; f < g.nfield; ++f) {
                    const int gi = g.idx[f];
                    grad_top(gi) += grad_aug(d + f);
                    for (int c = 0; c < d; ++c) {
                        const double cross = small(d + f, c);
                        T(c, gi) += cross;
                        T(gi, c) += cross;
                        C[s](gi, c) += cross;
                    }
                    for (int f2 = 0; f2 < g.nfield; ++f2)
                        T(gi, g.idx[f2]) += small(d + f, d + f2);
                }
            }
            // Per-species RE prior.
            D[s].noalias() += P;
            grad_b[s].noalias() -= P * b_[s];
        }
        // Field prior on the top block (mu is flat -- see the no-ridge note).
        add_field_prior(kind, p_state, p_p, n_spatial, tau, rho,
                        adj_row_ptr, adj_col_idx, n_neighbors, field_,
                        grad_top, T, spde);
        if (log_lik_out) *log_lik_out = ll;
        if (boundary_out) *boundary_out = bmax;
    };

    VectorXd grad_top(m);
    MatrixXd T(m, m);
    std::vector<VectorXd> grad_b(S, VectorXd::Zero(d));
    std::vector<MatrixXd> D(S, MatrixXd::Zero(d, d)), C(S, MatrixXd::Zero(m, d));
    std::vector<MatrixXd> Dinv(S, MatrixXd::Zero(d, d));

    for (int em_it = 0; em_it < max_iter_em; ++em_it) {
        out.n_iter = em_it + 1;
        const MatrixXd P = block_prec(Sig_state, Sig_p);

        // ---- mode-find: joint Newton over (mu, field, {b_s}) ----
        for (int nit = 0; nit < inner_max; ++nit) {
            double ll_cur = 0.0;
            assemble(mu, field, bvec, P, /*want_obs=*/false,
                     grad_top, T, grad_b, D, C, &ll_cur, nullptr);
            // Block elimination: M = T - sum_s C_s D_s^{-1} C_s'.
            MatrixXd M = T;
            VectorXd rhs = grad_top;
            for (int s = 0; s < S; ++s) {
                Dinv[s] = nmix_safe_inverse(D[s]);
                const MatrixXd CD = C[s] * Dinv[s];          // m x d
                M.noalias()   -= CD * C[s].transpose();
                rhs.noalias() -= CD * grad_b[s];
            }
            // M is a Schur complement of a complete-data-Fisher arrowhead, so it
            // is PSD by construction and the ridged safe inverse is the whole
            // fallback -- there is no separate Fisher matrix to retry on.
            nmix_add_diagonal_ridge(M);
            const VectorXd delta_top = nmix_safe_inverse(M) * rhs;
            std::vector<VectorXd> delta_b(S, VectorXd::Zero(d));
            for (int s = 0; s < S; ++s)
                delta_b[s] = Dinv[s] * (grad_b[s] - C[s].transpose() * delta_top);

            const VectorXd dmu = delta_top.head(d);
            const VectorXd dfield = delta_top.segment(d, field_len);

            // Step halving on the joint objective. obj_cur reuses the data
            // log-lik just assembled (no redundant likelihood pass).
            double obj_cur = ll_cur
                + field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                  adj_row_ptr, adj_col_idx, n_neighbors, field, spde);
            for (int s = 0; s < S; ++s) obj_cur -= 0.5 * bvec[s].dot(P * bvec[s]);
            VectorXd mu_try, field_try;
            std::vector<VectorXd> b_try(S);
            double max_step = 0.0;
            const bool stepped = newton_backtrack(
                obj_cur,
                [&](double step) {
                    mu_try = mu + step * dmu;
                    field_try = field + step * dfield;
                    for (int s = 0; s < S; ++s) b_try[s] = bvec[s] + step * delta_b[s];
                    return objective(mu_try, field_try, b_try, P);
                },
                [&](double step) {
                    mu = mu_try; field = field_try; bvec = b_try;
                    center_field(kind, p_state, p_p, n_spatial, field);
                    double dmax = std::max(dmu.cwiseAbs().maxCoeff(),
                                           dfield.cwiseAbs().maxCoeff());
                    for (int s = 0; s < S; ++s)
                        dmax = std::max(dmax, delta_b[s].cwiseAbs().maxCoeff());
                    max_step = step * dmax;
                });
            if (!stepped) break;
            if (max_step < inner_tol) break;
        }

        // ---- recompute Fisher blocks at the converged mode for the M-step ----
        assemble(mu, field, bvec, P, /*want_obs=*/false,
                 grad_top, T, grad_b, D, C, nullptr, nullptr);
        for (int s = 0; s < S; ++s) Dinv[s] = nmix_safe_inverse(D[s]);

        // ---- M-step: EM update of the community covariances ----
        MatrixXd Sstate_new = MatrixXd::Zero(p_state, p_state);
        MatrixXd Sp_new     = MatrixXd::Zero(p_p, p_p);
        for (int s = 0; s < S; ++s) {
            const VectorXd bs = bvec[s].head(p_state);
            const VectorXd bp = bvec[s].tail(p_p);
            Sstate_new += bs * bs.transpose() + Dinv[s].topLeftCorner(p_state, p_state);
            Sp_new     += bp * bp.transpose() + Dinv[s].bottomRightCorner(p_p, p_p);
        }
        Sstate_new /= (double) S;
        Sp_new     /= (double) S;
        const double dSig = std::max((Sstate_new - Sig_state).cwiseAbs().maxCoeff(),
                                     (Sp_new     - Sig_p).cwiseAbs().maxCoeff());
        Sig_state = Sstate_new;
        Sig_p     = Sp_new;
        if (verbose) Rcpp::Rcout << "  em " << out.n_iter << " dSigma=" << dSig << "\n";
        if (dSig < tol_em) { out.converged = true; break; }
    }

    // ---- final pass: observed-info arrowhead, log-marginal, vcov_mu ----
    const MatrixXd P = block_prec(Sig_state, Sig_p);
    const double logdetP = nmix_logdet_spd(P);

    // Assemble observed-info blocks; fall back to complete-data Fisher if the
    // observed-info determinant is not finite (indefinite away from the mode).
    auto final_assemble = [&](bool want_obs, double& loglik_marg,
                              MatrixXd& vcov_mu, double& boundary_max) -> bool {
        VectorXd gtop(m); MatrixXd Tt(m, m);
        std::vector<VectorXd> gb(S, VectorXd::Zero(d));
        std::vector<MatrixXd> Dd(S, MatrixXd::Zero(d, d)), Cc(S, MatrixXd::Zero(m, d));
        double ll = 0.0;
        assemble(mu, field, bvec, P, want_obs, gtop, Tt, gb, Dd, Cc, &ll, &boundary_max);

        double sum_logdet_D = 0.0;
        MatrixXd M = Tt;
        for (int s = 0; s < S; ++s) {
            const double ldD = nmix_logdet_spd(Dd[s]);
            if (!R_finite(ldD)) return false;
            sum_logdet_D += ldD;
            const MatrixXd Dinv_s = nmix_safe_inverse(Dd[s]);
            M.noalias() -= Cc[s] * Dinv_s * Cc[s].transpose();
        }
        // log|M| (rank-deficient field null pinned by a tiny ridge, as in the
        // single-species path).
        MatrixXd M_det = M;
        nmix_add_diagonal_ridge(M_det);
        const double ldM = nmix_logdet_spd(M_det);
        if (!R_finite(ldM)) return false;

        // Per-species b-prior normaliser + the b-Schur determinant.
        double lp = field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                    adj_row_ptr, adj_col_idx, n_neighbors, field, spde);
        double bquad = 0.0;
        for (int s = 0; s < S; ++s) bquad += bvec[s].dot(P * bvec[s]);
        loglik_marg = ll + lp
                      + 0.5 * (double) S * logdetP - 0.5 * bquad
                      - 0.5 * sum_logdet_D - 0.5 * ldM;

        // vcov_mu = top-left d-block of M^{-1}, constrained (sum field = 0) for
        // the rank-deficient intrinsic fields (ICAR / BYM2 v); CAR_proper and
        // SPDE are full rank.
        const bool constrain = (kind != FieldKind::CAR_PROPER &&
                                kind != FieldKind::SPDE);
        MatrixXd cov_top = nmix_constrained_top_cov(
            M, m, d, /*field_start=*/d, /*field_len=*/field_len, constrain);
        if (!cov_top.allFinite()) return false;
        vcov_mu = cov_top;
        return true;
    };

    double loglik_marg = R_NegInf, boundary_max = 0.0;
    MatrixXd vcov_mu;
    bool ok = final_assemble(/*want_obs=*/true, loglik_marg, vcov_mu, boundary_max);
    if (!ok)
        ok = final_assemble(/*want_obs=*/false, loglik_marg, vcov_mu, boundary_max);

    out.mu = mu;
    out.field = field;
    out.Sigma_state = Sig_state;
    out.Sigma_p = Sig_p;
    out.blup_state = MatrixXd(S, p_state);
    out.blup_p     = MatrixXd(S, p_p);
    for (int s = 0; s < S; ++s) {
        out.blup_state.row(s) = bvec[s].head(p_state).transpose();
        out.blup_p.row(s)     = bvec[s].tail(p_p).transpose();
    }
    out.vcov_mu = ok ? vcov_mu : MatrixXd::Constant(d, d, R_NaN);
    out.log_marginal = ok ? loglik_marg : R_NegInf;
    out.log_lik = data_loglik(mu, field, bvec);
    out.boundary_max = boundary_max;
    return out;
}

// One outer-grid point: the field hyperparameters (tau / rho / log|Q(rho)|, the
// BYM2 mixing loadings a, b, or the Matern precision the SPDE field carries)
// that the inner Laplace-EM conditions on. Each entry builds a plan of these in
// its own grid-nesting order; the driver below just walks the plan.
struct CommGridPoint {
    double tau = 1.0, rho = 1.0, log_det_Q_rho = 0.0, a = 0.0, b = 0.0;
    SpdeCtx spde;
};

// Outer-grid driver for the community shared-field nested-Laplace-EM. Owns the
// init unpacking, the grid walk, progress plumbing, and the community_pack_grid
// output assembly. Each entry supplies its FieldKind, its field length, the
// per-point plan + theta_grid_out it filled, and `family_at(point)` returning
// the response family to fit that point with (the N-mixture family carries the
// point's NB size, so the family is built per point).
template <class MakeFamily>
Rcpp::List run_community_spatial_grid(
    FieldKind kind, int field_len,
    MakeFamily family_at,
    int p_state, int p_p,
    const std::vector<int>& map_site_to_unit,
    int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const std::vector<CommGridPoint>& plan,
    const Rcpp::NumericMatrix& theta_grid_out,
    const Rcpp::NumericVector& mu_init,
    const Rcpp::NumericMatrix& Sigma_state_init,
    const Rcpp::NumericMatrix& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose,
    const char* sigma_state_name, const char* blup_state_name,
    const char* p_state_name, bool report_boundary,
    bool progress, int progress_every, double progress_throttle,
    const std::string& progress_file) {

    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;

    const int d = p_state + p_p;
    if ((int) mu_init.size() != d)
        Rcpp::stop("mu_init length must equal the state-arm plus detection-arm width.");
    if (Sigma_state_init.nrow() != p_state || Sigma_state_init.ncol() != p_state)
        Rcpp::stop("The state-arm community covariance init must be %d x %d.",
                   p_state, p_state);
    if (Sigma_p_init.nrow() != p_p || Sigma_p_init.ncol() != p_p)
        Rcpp::stop("The detection-arm community covariance init must be %d x %d.",
                   p_p, p_p);
    VectorXd mu0     = Map<VectorXd>(REAL(mu_init), d);
    MatrixXd Sstate0 = Map<MatrixXd>(REAL(Sigma_state_init), p_state, p_state);
    MatrixXd Sp0     = Map<MatrixXd>(REAL(Sigma_p_init), p_p, p_p);

    const int n_grid = (int) plan.size();
    std::vector<CommSpatialResult> results(n_grid);
    // outer-grid progress
    auto gp = make_grid_progress("community-spatial", n_grid, progress,
                                 progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        const CommGridPoint& g = plan[k];
        results[k] = community_spatial_em(
            kind, family_at(g), p_state, p_p, map_site_to_unit, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors,
            g.tau, g.rho, g.log_det_Q_rho, g.a, g.b,
            mu0, Sstate0, Sp0,
            max_iter_em, tol_em, inner_max, inner_tol, verbose,
            kind == FieldKind::SPDE ? &g.spde : nullptr);
        if (gp) gp->tick();
    }
    if (gp) gp->finish();
    return community_pack_grid(
        d, field_len, n_grid, results, theta_grid_out,
        &CommSpatialResult::Sigma_state, &CommSpatialResult::blup_state,
        sigma_state_name, blup_state_name, p_state_name, p_state, p_p, n_spatial,
        report_boundary ? &CommSpatialResult::boundary_max : nullptr);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_COMMUNITY_SPATIAL_DRIVER_H
