// community_spatial_em.h
// Shared field-structure layer + Laplace-EM driver for a per-species marginal
// with Gaussian community hyperpriors on the per-species coefficients AND one
// shared ICAR / BYM2 / proper-CAR / SPDE field on the STATE arm (the sfMsNMix /
// svcMsPGOcc pattern). Two response families instantiate this: the community
// Royle N-mixture (nmix_community_spatial.cpp) and the community two-state
// occupancy marginal (ms_occu_spatial.cpp). The algorithm is identical between
// them -- a Laplace-EM over species random effects with a shared areal field,
// solved by a Schur-complement Newton step with step-halving, then an EM update
// of the community covariances, wrapped in an outer hyperparameter grid; only
// the per-species marginal (Royle N-mixture vs occupancy two-state) differs.
//
// The family-specific piece is ONE callable, `SiteBlockFn`:
//
//   double site_block_fn(int s, int site, const VectorXd& coef,
//                        double field_offset, const double* load, int nfield,
//                        bool want_block, bool want_obs,
//                        VectorXd& grad_aug, MatrixXd& small,
//                        double* boundary_out)
//
// evaluating one (species, site) cell: `coef` is that species' full
// (state, detection) coefficient vector, `field_offset` the shared field's
// current contribution to eta_state at that site, `load`/`nfield` the field
// loading geometry (see FieldGeom below). It fills `grad_aug` (length
// d + nfield: coefficient columns then field columns) and, when `want_block`,
// the symmetric curvature `small` ((d + nfield) x (d + nfield)) via the
// design-sandwiched per-site eta-space block; returns the site's marginal
// log-lik. Each family builds this from its own per-site kernel
// (nmix_kernel.h's compute_nmix_site_cached, ms_occu_kernel.h's
// ms_occu_site_cell) and adapts it to this signature with a small lambda that
// captures its own data container + design; that adapter is the ONLY
// family-specific code left outside this header.
//
// Everything else -- field geometry / offset / log-prior / prior-gradient /
// sum-to-zero centering (areal + SPDE), the arrowhead block-elimination Newton,
// the EM covariance update, the observed-info final pass + community-mean
// covariance, and the outer hyperparameter-grid driver (areal and SPDE) -- is
// family-agnostic and lives here once.

#ifndef TULPAOBS_COMMUNITY_SPATIAL_EM_H
#define TULPAOBS_COMMUNITY_SPATIAL_EM_H

#include "nmix_spatial_kernel.h"        // icar/car_proper prior, Q-adds, centering
#include "nmix_spatial_kernel_bym2.h"   // bym2 prior, Q+I-adds, centering
#include "nmix_linalg.h"
#include "nmix_progress.h"
#include "community_grid_pack.h"
#include "newton_step.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace tulpaObs {

using Eigen::MatrixXd;
using Eigen::VectorXd;

enum class CommFieldKind { ICAR, CAR_PROPER, BYM2, SPDE };

// SPDE field context: the continuous Matern field lives at n_mesh FEM nodes;
// the dense projection A (n_sites x n_mesh) maps mesh nodes onto sites, and Q
// (n_mesh x n_mesh, carrying tau_spde^2) is the proper Matern precision built
// once per grid point on the R side. The areal kinds carry a single field unit
// per site; SPDE carries the whole mesh row, so the field loading and prior are
// threaded through this struct rather than the (u, a, b) areal triple.
struct CommSpdeCtx {
    const MatrixXd* A = nullptr;     // n_sites x n_mesh
    const MatrixXd* Q = nullptr;     // n_mesh x n_mesh (with tau_spde^2)
    double log_det_Q = 0.0;
    int n_mesh = 0;
};

// Field geometry for one site: the field-column loadings on eta_state and their
// global indices in the top block (field_start = d). Areal kinds load a single
// field unit (BYM2 two: v + w); the SPDE field loads the mesh nodes of the
// site's projection row (a 2D P1 FEM row has at most MAX_FIELD_LOAD nonzeros --
// the barycentric weights of the triangle the site falls in).
constexpr int kCommMaxFieldLoad = 4;
struct CommFieldGeom {
    int nfield;
    double load[kCommMaxFieldLoad];
    int    idx[kCommMaxFieldLoad];
};

// Areal field geometry for spatial unit u.
inline CommFieldGeom comm_field_geom(CommFieldKind kind, int d, int n_spatial,
                                     int u, double a, double b) {
    CommFieldGeom g;
    if (kind == CommFieldKind::BYM2) {
        g.nfield = 2;
        g.load[0] = a; g.idx[0] = d + u;                 // v
        g.load[1] = b; g.idx[1] = d + n_spatial + u;     // w
    } else {
        g.nfield = 1;
        g.load[0] = 1.0; g.idx[0] = d + u;               // f
    }
    return g;
}

// SPDE field geometry for one site: read the nonzeros of A.row(site).
inline CommFieldGeom comm_field_geom_spde(int d, const CommSpdeCtx& sp, int site) {
    CommFieldGeom g;
    g.nfield = 0;
    for (int k = 0; k < sp.n_mesh; ++k) {
        const double w = (*sp.A)(site, k);
        if (w != 0.0) {
            if (g.nfield >= kCommMaxFieldLoad)
                Rcpp::stop("SPDE projection row has more than %d nonzeros; "
                           "raise kCommMaxFieldLoad.", kCommMaxFieldLoad);
            g.load[g.nfield] = w;
            g.idx[g.nfield]  = d + k;
            ++g.nfield;
        }
    }
    return g;
}

inline double comm_field_offset(CommFieldKind kind, const VectorXd& field,
                                int n_spatial, int u, double a, double b) {
    if (kind == CommFieldKind::BYM2) return a * field(u) + b * field(n_spatial + u);
    return field(u);
}

// SPDE per-site offset: (A u)_site = sum_k A[site, k] u[k].
inline double comm_field_offset_spde(const VectorXd& field, const CommSpdeCtx& sp,
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
inline double comm_field_log_prior(CommFieldKind kind, int n_spatial,
                                   double tau, double rho, double log_det_Q_rho,
                                   const Rcpp::IntegerVector& adj_row_ptr,
                                   const Rcpp::IntegerVector& adj_col_idx,
                                   const Rcpp::IntegerVector& n_neighbors,
                                   const VectorXd& field,
                                   const CommSpdeCtx* spde = nullptr) {
    if (kind == CommFieldKind::SPDE) {
        // 0.5 log|Q| - 0.5 u' Q u (full rank; the (2 pi)^{-n_mesh/2} constant is
        // grid-independent and dropped, as in the single-species SPDE path).
        return 0.5 * spde->log_det_Q
             - 0.5 * (field.transpose() * ((*spde->Q) * field))(0, 0);
    }
    if (kind == CommFieldKind::BYM2) {
        const VectorXd v = field.head(n_spatial);
        const VectorXd w = field.tail(n_spatial);
        return tulpaObs::nmix_bym2_log_prior(n_spatial, adj_row_ptr, adj_col_idx,
                                             n_neighbors, v, w);
    }
    if (kind == CommFieldKind::ICAR) {
        return tulpaObs::nmix_icar_log_prior(n_spatial, tau, adj_row_ptr,
                                             adj_col_idx, n_neighbors, field);
    }
    return tulpaObs::nmix_car_proper_log_prior(n_spatial, tau, rho, log_det_Q_rho,
                                               adj_row_ptr, adj_col_idx,
                                               n_neighbors, field);
}

// Add the field prior to the top-block gradient and Hessian (in-place). The top
// block has the (beta_state, beta_p, field) layout the single-species helpers
// expect (field_start = d = p_state + p_p).
inline void comm_add_field_prior(CommFieldKind kind, int p_state, int p_p,
                                 int n_spatial, double tau, double rho,
                                 const Rcpp::IntegerVector& adj_row_ptr,
                                 const Rcpp::IntegerVector& adj_col_idx,
                                 const Rcpp::IntegerVector& n_neighbors,
                                 const VectorXd& field,
                                 VectorXd& grad_top, MatrixXd& T,
                                 const CommSpdeCtx* spde = nullptr) {
    if (kind == CommFieldKind::SPDE) {
        // Field prior -Q u to the score, +Q to the negative-Hessian field block.
        const int d = p_state + p_p;
        grad_top.segment(d, spde->n_mesh).noalias() -= (*spde->Q) * field;
        T.block(d, d, spde->n_mesh, spde->n_mesh).noalias() += (*spde->Q);
        return;
    }
    if (kind == CommFieldKind::BYM2) {
        const VectorXd v = field.head(n_spatial);
        const VectorXd w = field.tail(n_spatial);
        tulpaObs::nmix_add_bym2_prior_to_grad_and_H(p_state, p_p, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors, v, w, grad_top, T);
        return;
    }
    const double rho_use = (kind == CommFieldKind::CAR_PROPER) ? rho : 1.0;
    tulpaObs::nmix_add_car_to_spatial_block(p_state, p_p, n_spatial, tau, rho_use,
        adj_row_ptr, adj_col_idx, n_neighbors, field, grad_top, T);
}

// Sum-to-zero centering of the (intrinsic) field component against the global
// state-arm intercept. ICAR centres f; BYM2 centres v only (w is proper);
// CAR_proper and SPDE are full-rank and need none.
inline void comm_center_field(CommFieldKind kind, int p_state, int p_p,
                              int n_spatial, VectorXd& field) {
    if (kind == CommFieldKind::CAR_PROPER || kind == CommFieldKind::SPDE ||
        n_spatial <= 0) return;
    const int len = field.size();
    VectorXd holder(p_state + p_p + len);
    holder.setZero();
    holder.segment(p_state + p_p, len) = field;
    tulpaObs::nmix_center_field(p_state, p_p, n_spatial, holder);
    field = holder.segment(p_state + p_p, len);
}

struct CommSpatialResult {
    VectorXd mu;             // d community means
    VectorXd field;          // n_field_total = nfield * n_spatial (ICAR/CAR/SPDE: f; BYM2: [v; w])
    MatrixXd Sigma_state, Sigma_p;
    MatrixXd blup_state, blup_p;   // S x p_state, S x p_p
    MatrixXd vcov_mu;              // d x d community-mean covariance (b- and field-folded)
    double log_marginal = R_NegInf;
    double log_lik = R_NegInf;
    bool   converged = false;
    int    n_iter = 0;
    double boundary_max = 0.0;
};

// One outer-grid-point Laplace-EM fit, over the marginal `site_block_fn`
// evaluates. `field` is laid out [f] (ICAR/CAR/SPDE) or [v; w] (BYM2); `a`, `b`
// are the BYM2 loadings (unused otherwise); `map_site_to_unit` is unused (may be
// empty) for SPDE, which reads its geometry from `spde` instead.
template <class SiteBlockFn>
CommSpatialResult community_spatial_em(
    CommFieldKind kind,
    int S, int n_sites, int p_state, int p_p,
    const std::vector<int>& map_site_to_unit,
    int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    double tau, double rho, double log_det_Q_rho, double a, double b,
    const VectorXd& mu_init, const MatrixXd& Sigma_state_init,
    const MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose, const char* verbose_label,
    SiteBlockFn site_block_fn,
    const CommSpdeCtx* spde = nullptr) {

    const int d = p_state + p_p;
    const int nfield = (kind == CommFieldKind::BYM2) ? 2 : 1;
    const int field_len = nfield * n_spatial;
    const int m = d + field_len;
    // No ridge on the community means (mu). The state-arm intercept and the
    // field constant mode form an exactly flat (data + ICAR-prior) direction;
    // any ridge on mu is the ONLY curvature along it and drives the intercept to
    // its prior mean 0, with the field absorbing the level (then deleted by the
    // sum-to-zero centering). The single-species spatial path keeps the fixed
    // effects flat for the same reason -- mu is identified by the data once the
    // field is anchored sum-to-zero.

    CommSpatialResult out;
    VectorXd mu = mu_init;
    VectorXd field = VectorXd::Zero(field_len);
    std::vector<VectorXd> bvec(S, VectorXd::Zero(d));
    MatrixXd Sig_state = Sigma_state_init, Sig_p = Sigma_p_init;

    auto block_prec = [&](const MatrixXd& Sl, const MatrixXd& Sp) {
        MatrixXd P = MatrixXd::Zero(d, d);
        P.topLeftCorner(p_state, p_state) = nmix_safe_inverse(Sl);
        P.bottomRightCorner(p_p, p_p) = nmix_safe_inverse(Sp);
        return P;
    };

    // Per-site (geometry, offset) -- branches on the field kind once. SPDE reads
    // the site's projection row; the areal kinds use the single-unit map.
    auto site_geom = [&](int site, const VectorXd& field_) {
        if (kind == CommFieldKind::SPDE)
            return std::make_pair(comm_field_geom_spde(d, *spde, site),
                                  comm_field_offset_spde(field_, *spde, site));
        const int u = map_site_to_unit[site];
        return std::make_pair(comm_field_geom(kind, d, n_spatial, u, a, b),
                              comm_field_offset(kind, field_, n_spatial, u, a, b));
    };

    // Data log-lik at the current state (line-search objective helper).
    auto data_loglik = [&](const VectorXd& mu_, const VectorXd& field_,
                           const std::vector<VectorXd>& b_) {
        double ll = 0.0;
        VectorXd grad_aug; MatrixXd small;
        for (int s = 0; s < S; ++s) {
            VectorXd coef = mu_ + b_[s];
            for (int site = 0; site < n_sites; ++site) {
                auto go = site_geom(site, field_);
                const CommFieldGeom& g = go.first;
                double l = site_block_fn(s, site, coef, go.second, g.load, g.nfield,
                                         false, false, grad_aug, small, nullptr);
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
        obj += comm_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
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
            for (int site = 0; site < n_sites; ++site) {
                auto go = site_geom(site, field_);
                const CommFieldGeom& g = go.first;
                double bw = 0.0;
                ll += site_block_fn(s, site, coef, go.second, g.load, g.nfield,
                                    true, want_obs, grad_aug, small, &bw);
                if (bw > bmax) bmax = bw;
                // d-block: per-species coefficient curvature + community mean.
                grad_b[s].head(d) += grad_aug.head(d);
                grad_top.head(d)  += grad_aug.head(d);
                D[s].noalias()                += small.topLeftCorner(d, d);
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
        comm_add_field_prior(kind, p_state, p_p, n_spatial, tau, rho,
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
            tulpaObs::nmix_add_diagonal_ridge(M);
            const VectorXd delta_top = nmix_safe_inverse(M) * rhs;
            std::vector<VectorXd> delta_b(S, VectorXd::Zero(d));
            for (int s = 0; s < S; ++s)
                delta_b[s] = Dinv[s] * (grad_b[s] - C[s].transpose() * delta_top);

            const VectorXd dmu = delta_top.head(d);
            const VectorXd dfield = delta_top.segment(d, field_len);

            // Step halving on the joint objective. obj_cur reuses the data
            // log-lik just assembled (no redundant likelihood pass).
            double obj_cur = ll_cur
                + comm_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                       adj_row_ptr, adj_col_idx, n_neighbors, field, spde);
            for (int s = 0; s < S; ++s) obj_cur -= 0.5 * bvec[s].dot(P * bvec[s]);
            VectorXd mu_try, field_try;
            std::vector<VectorXd> b_try(S);
            double max_step = 0.0;
            const bool stepped = tulpaObs::newton_backtrack(
                obj_cur,
                [&](double step) {
                    mu_try = mu + step * dmu;
                    field_try = field + step * dfield;
                    for (int s = 0; s < S; ++s) b_try[s] = bvec[s] + step * delta_b[s];
                    return objective(mu_try, field_try, b_try, P);
                },
                [&](double step) {
                    mu = mu_try; field = field_try; bvec = b_try;
                    comm_center_field(kind, p_state, p_p, n_spatial, field);
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
        MatrixXd Sp_new = MatrixXd::Zero(p_p, p_p);
        for (int s = 0; s < S; ++s) {
            const VectorXd bst = bvec[s].head(p_state);
            const VectorXd bp = bvec[s].tail(p_p);
            Sstate_new += bst * bst.transpose() + Dinv[s].topLeftCorner(p_state, p_state);
            Sp_new     += bp * bp.transpose()   + Dinv[s].bottomRightCorner(p_p, p_p);
        }
        Sstate_new /= (double)S;
        Sp_new /= (double)S;
        const double dSig = std::max((Sstate_new - Sig_state).cwiseAbs().maxCoeff(),
                                     (Sp_new - Sig_p).cwiseAbs().maxCoeff());
        Sig_state = Sstate_new;
        Sig_p = Sp_new;
        if (verbose)
            Rcpp::Rcout << "  " << verbose_label << " " << out.n_iter
                       << " dSigma=" << dSig << "\n";
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
        tulpaObs::nmix_add_diagonal_ridge(M_det);
        const double ldM = nmix_logdet_spd(M_det);
        if (!R_finite(ldM)) return false;

        // Per-species b-prior normaliser + the b-Schur determinant.
        double lp = comm_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                         adj_row_ptr, adj_col_idx, n_neighbors, field, spde);
        double bquad = 0.0;
        for (int s = 0; s < S; ++s) bquad += bvec[s].dot(P * bvec[s]);
        loglik_marg = ll + lp
                      + 0.5 * (double)S * logdetP - 0.5 * bquad
                      - 0.5 * sum_logdet_D - 0.5 * ldM;

        // vcov_mu = top-left d-block of M^{-1}, constrained (sum field = 0) for
        // the rank-deficient intrinsic fields (ICAR / BYM2 v); CAR_proper and
        // SPDE are full rank.
        const bool constrain = (kind != CommFieldKind::CAR_PROPER &&
                                kind != CommFieldKind::SPDE);
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
    out.blup_p = MatrixXd(S, p_p);
    for (int s = 0; s < S; ++s) {
        out.blup_state.row(s) = bvec[s].head(p_state).transpose();
        out.blup_p.row(s) = bvec[s].tail(p_p).transpose();
    }
    out.vcov_mu = ok ? vcov_mu : MatrixXd::Constant(d, d, R_NaN);
    out.log_marginal = ok ? loglik_marg : R_NegInf;
    out.log_lik = data_loglik(mu, field, bvec);
    out.boundary_max = boundary_max;
    return out;
}

// One outer-grid point: the field hyperparameters (tau / rho / log|Q(rho)| and
// the BYM2 mixing loadings a, b) plus, for count families, a per-point nuisance
// (NB size r) the inner Laplace-EM conditions on. `r` is ignored by families
// with no such nuisance (occupancy passes a placeholder). Each entry builds a
// plan of these in its own grid-nesting order; the driver below just walks it.
struct CommGridPoint { double tau, rho, log_det_Q_rho, a, b, r; };

// Shared outer-grid driver for an AREAL community spatial nested-Laplace-EM
// (icar / car_proper / bym2). `make_site_block_fn(r)` builds the family-
// specific per-site evaluator closure for one grid point's nuisance `r`
// (ignored by families with none). Progress reporting is opt-in per call
// (occupancy has none today; count families report under "nmix-spatial").
template <class SiteBlockFactory>
Rcpp::List run_community_spatial_grid(
    CommFieldKind kind, int field_len,
    int S, int n_sites, int p_state, int p_p,
    const std::vector<int>& map_site_to_unit, int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const std::vector<CommGridPoint>& plan,
    const Rcpp::NumericMatrix& theta_grid_out,
    const VectorXd& mu_init, const MatrixXd& Sigma_state_init,
    const MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose, const char* verbose_label,
    const char* sigma_state_name, const char* blup_state_name,
    const char* p_state_name,
    SiteBlockFactory make_site_block_fn,
    bool progress, int progress_every, double progress_throttle,
    const std::string& progress_file, const char* progress_grid_label,
    bool report_boundary) {

    const int n_grid = (int)plan.size();
    std::vector<CommSpatialResult> results(n_grid);
    auto gp = tulpaObs::make_grid_progress(progress_grid_label, n_grid, progress,
                                           progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        const CommGridPoint& g = plan[k];
        results[k] = community_spatial_em(
            kind, S, n_sites, p_state, p_p, map_site_to_unit, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors,
            g.tau, g.rho, g.log_det_Q_rho, g.a, g.b,
            mu_init, Sigma_state_init, Sigma_p_init,
            max_iter_em, tol_em, inner_max, inner_tol, verbose, verbose_label,
            make_site_block_fn(g.r));
        if (gp) gp->tick();
    }
    if (gp) gp->finish();
    return tulpaObs::community_pack_grid(
        p_state + p_p, field_len, n_grid, results, theta_grid_out,
        &CommSpatialResult::Sigma_state, &CommSpatialResult::blup_state,
        sigma_state_name, blup_state_name, p_state_name, p_state, p_p, n_spatial,
        report_boundary ? &CommSpatialResult::boundary_max : nullptr);
}

// Shared outer-grid driver for the continuous Matern (SPDE) shared field. The
// field lives at n_mesh FEM nodes; the dense projection A (n_sites x n_mesh)
// maps mesh nodes onto sites, shared across species exactly as the areal field
// is. Per grid point the proper Matern precision Q(range, sigma) (carrying
// tau_spde^2) and its log|Q| are built once on the R side and passed in via
// `Q_list` / `log_det_Q`; the outer grid axes are (range, sigma) [, a per-point
// nuisance `r_grid`].
template <class SiteBlockFactory>
Rcpp::List run_community_spatial_grid_spde(
    int S, const MatrixXd& A, int n_mesh, int p_state, int p_p,
    const Rcpp::List& Q_list, const Rcpp::NumericVector& log_det_Q,
    const Rcpp::NumericMatrix& theta_grid_R, const Rcpp::NumericVector& r_grid,
    const VectorXd& mu_init, const MatrixXd& Sigma_state_init,
    const MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose, const char* verbose_label,
    const char* sigma_state_name, const char* blup_state_name,
    const char* p_state_name,
    SiteBlockFactory make_site_block_fn,
    bool progress, int progress_every, double progress_throttle,
    const std::string& progress_file, const char* progress_grid_label,
    bool report_boundary) {

    const int n_grid = Q_list.size();
    if ((int)log_det_Q.size() != n_grid)
        Rcpp::stop("length(log_det_Q) must equal length(Q_list).");
    if (theta_grid_R.nrow() != n_grid)
        Rcpp::stop("nrow(theta_grid) must equal length(Q_list).");
    if ((int)r_grid.size() != n_grid)
        Rcpp::stop("length(r_grid) must equal length(Q_list).");

    // Empty CSR adjacency (unused on the SPDE path).
    Rcpp::IntegerVector empty_int(0);
    const int n_sites = (int)A.rows();

    std::vector<CommSpatialResult> results(n_grid);
    std::vector<MatrixXd> Qmats(n_grid);   // keep alive for the CommSpdeCtx pointers
    auto gp = tulpaObs::make_grid_progress(progress_grid_label, n_grid, progress,
                                           progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        Rcpp::NumericMatrix Qk_R = Q_list[k];
        if (Qk_R.nrow() != n_mesh || Qk_R.ncol() != n_mesh)
            Rcpp::stop("Q_list[[%d]] must be n_mesh x n_mesh.", k + 1);
        Qmats[k] = Eigen::Map<MatrixXd>(REAL(Qk_R), n_mesh, n_mesh);
        CommSpdeCtx sp;
        sp.A = &A; sp.Q = &Qmats[k]; sp.log_det_Q = log_det_Q[k]; sp.n_mesh = n_mesh;
        results[k] = community_spatial_em(
            CommFieldKind::SPDE, S, n_sites, p_state, p_p,
            /*map=*/std::vector<int>(), /*n_spatial=*/n_mesh,
            empty_int, empty_int, empty_int,
            /*tau=*/1.0, /*rho=*/1.0, /*log_det_Q_rho=*/0.0, /*a=*/0.0, /*b=*/0.0,
            mu_init, Sigma_state_init, Sigma_p_init,
            max_iter_em, tol_em, inner_max, inner_tol, verbose, verbose_label,
            make_site_block_fn(r_grid[k]), &sp);
        if (gp) gp->tick();
    }
    if (gp) gp->finish();
    return tulpaObs::community_pack_grid(
        p_state + p_p, n_mesh, n_grid, results, theta_grid_R,
        &CommSpatialResult::Sigma_state, &CommSpatialResult::blup_state,
        sigma_state_name, blup_state_name, p_state_name, p_state, p_p, n_mesh,
        report_boundary ? &CommSpatialResult::boundary_max : nullptr);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_COMMUNITY_SPATIAL_EM_H
