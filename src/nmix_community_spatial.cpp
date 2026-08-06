// nmix_community_spatial.cpp
// Nested Laplace-EM for the SPATIAL community / multispecies N-mixture
// (spAbundance sfMsNMix analogue): a per-species Royle (2004) N-mixture with
// Gaussian community hyperpriors on the per-species coefficients AND one shared
// areal field f on the abundance arm.
//
//   N_{s,i}        ~ Poisson(lambda_{s,i})            (or NB, global size r)
//   y_{s,i,j} | N  ~ Binomial(N_{s,i}, p_{s,i,j})
//   log lambda_{s,i} = X_lambda_i . (mu_lambda + b_lambda_s) + f_{u(i)}
//   logit p_{s,i,j}  = X_p_{ij}   . (mu_p      + b_p_s)
//   b_lambda_s ~ N(0, Sigma_lambda),  b_p_s ~ N(0, Sigma_p)
//   f ~ ICAR(tau) | BYM2(sigma,rho) | CAR(tau,rho)
//
// Three latent layers, integrated in three ways:
//   * N        -- summed analytically per site by the nmix_kernel.h marginal.
//   * b_s      -- per-species Laplace (Gaussian community prior); the species
//                 are conditionally independent given (mu, f), so they fold out
//                 by a Schur complement (Louis 1982), exactly as the non-spatial
//                 community EM (nmix_community_em.cpp).
//   * (mu, f)  -- joint Laplace; f has the field prior, mu a weak ridge. mu and
//                 f are NOT conditionally independent (both load eta_lambda), so
//                 they are solved as ONE block.
//
// The "top" block t = (mu_lambda, mu_p, f) has the SAME layout as the single-
// species spatial state x = (beta_lambda, beta_p, field), so the field prior /
// centering / constrained-covariance helpers in nmix_spatial_kernel.h and
// nmix_spatial_kernel_bym2.h apply to t directly (field_start = p_lam + p_p).
//
// Per outer grid point (tau[, rho] / sigma, rho [, NB size r]):
//   EM (Sigma fixed within, updated across):
//     mode-find  -- joint Newton over (t, {b_s}) by block elimination on the
//                   complete-data-Fisher (PSD) arrowhead Hessian; step-halving
//                   on the joint objective; centre f sum-to-zero (ICAR / BYM2 v).
//     M-step     -- Sigma_k = mean_s [ b_s b_s' + Cov(b_s | y) ],
//                   Cov(b_s | y) = (Fisher_s + Sigma^{-1})^{-1}.
//   final pass   -- observed-info arrowhead; the b-folded (mu, f) info M; the
//                   grid-point Laplace log-marginal and the community-mean
//                   covariance vcov_mu = top-left block of M^{-1} (constrained).
//
// The per-site marginal math is nmix_kernel.h (the single source); the oracle
// (NMixCommunityOracle) is reused purely as the pre-grouped, lgamma-cached data
// container its constructor builds.

#include "nmix_kernel.h"
#include "nmix_community_oracle.h"
#include "nmix_linalg.h"
#include "community_grid_pack.h"
#include "nmix_progress.h"
#include "nmix_spatial_kernel.h"        // nmix_icar/car_proper prior, Q-adds, centering
#include "nmix_spatial_kernel_bym2.h"   // bym2 prior, Q+I-adds, centering
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

namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using tulpaObs::NMixCommunityOracle;
using tulpaObs::NMixSiteResult;
using tulpaObs::compute_nmix_site_cached;
using tulpaObs::nmix_safe_inverse;
using tulpaObs::nmix_logdet_spd;
using tulpaObs::nmix_constrained_top_cov;

enum class FieldKind { ICAR, CAR_PROPER, BYM2, SPDE };

// SPDE field context: the continuous Matern field lives at n_mesh FEM nodes;
// the dense projection A (n_sites x n_mesh) maps mesh nodes onto sites, and Q
// (n_mesh x n_mesh, carrying tau_spde^2) is the proper Matern precision built
// once per grid point on the R side. The areal kinds carry a single field unit
// per site; SPDE carries the whole mesh row, so the field loading and prior are
// threaded through this struct rather than the (u, a, b) areal triple.
struct SpdeCtx {
    const MatrixXd* A = nullptr;     // n_sites x n_mesh
    const MatrixXd* Q = nullptr;     // n_mesh x n_mesh (with tau_spde^2)
    double log_det_Q = 0.0;
    int n_mesh = 0;
};


// Per-site augmented assembly. na = d + nfield, with the d coefficient coords
// first and the nfield field loadings last (eta coord 0 = lambda, 1..J =
// visits). Fills grad_aug (na) and, when want_block, the symmetric small block
// (na x na): small = Z_aug' B Z_aug with B the per-site eta-space curvature
// (complete-data Fisher diagonal, plus -Var[N|y] rank-1 when want_obs). The
// d x d top-left of small IS the per-species coefficient curvature (the
// nmix_community_em block); the field rows/cols and crosses extend it to the
// shared field. Returns the per-site marginal log-lik. A no-visit (J = 0)
// species-site contributes nothing (interpolated by the field / prior), matching
// the single-species held-out convention.
double site_blocks(const NMixCommunityOracle::SiteRec& rec,
                   const Map<MatrixXd>& Xlam, int p_lam, int p_p,
                   const VectorXd& coef, double field_offset,
                   const double* load, int nfield, double r,
                   bool want_block, bool want_obs,
                   VectorXd& grad_aug, MatrixXd& small,
                   double* boundary_out = nullptr) {
    const int d  = p_lam + p_p;
    const int na = d + nfield;
    const int J  = rec.cache.n_visits;
    grad_aug.setZero(na);
    if (want_block) small.setZero(na, na);
    if (J == 0) return 0.0;

    double eta_lam = field_offset;
    for (int c = 0; c < p_lam; ++c) eta_lam += Xlam(rec.site, c) * coef(c);
    eta_lam = clamp_eta(eta_lam);
    std::vector<double> eta_p(J);
    for (int j = 0; j < J; ++j) {
        double v = 0.0;
        for (int c = 0; c < p_p; ++c) v += rec.Xp(j, c) * coef(p_lam + c);
        eta_p[j] = clamp_eta(v);
    }
    const NMixSiteResult res = compute_nmix_site_cached(rec.cache, eta_p.data(),
                                                        eta_lam, r);
    if (boundary_out) *boundary_out = res.boundary_weight;
    // Gradient: eta coord 0 -> mu_lambda cols (Xlam) and field cols (loadings);
    // eta coord 1+j -> mu_p cols (Xp).
    for (int c = 0; c < p_lam; ++c)
        grad_aug(c) += Xlam(rec.site, c) * res.grad_eta_lambda;
    for (int f = 0; f < nfield; ++f)
        grad_aug(d + f) += load[f] * res.grad_eta_lambda;
    for (int j = 0; j < J; ++j)
        for (int c = 0; c < p_p; ++c)
            grad_aug(p_lam + c) += rec.Xp(j, c) * res.grad_eta_p[j];

    if (!want_block) return res.log_lik;

    // Augmented design Z_aug ((1+J) x na): row 0 = lambda, rows 1+j = visits.
    MatrixXd Z = MatrixXd::Zero(1 + J, na);
    for (int c = 0; c < p_lam; ++c) Z(0, c) = Xlam(rec.site, c);
    for (int f = 0; f < nfield; ++f) Z(0, d + f) = load[f];
    for (int j = 0; j < J; ++j)
        for (int c = 0; c < p_p; ++c) Z(1 + j, p_lam + c) = rec.Xp(j, c);

    // Per-site eta-space block B ((1+J)^2): complete-data Fisher diagonal, plus
    // the -Var[N|y] rank-1 coupling (Louis 1982) when want_obs.
    MatrixXd B = MatrixXd::Zero(1 + J, 1 + J);
    B(0, 0) = res.info_eta_lambda;
    for (int j = 0; j < J; ++j) B(1 + j, 1 + j) = res.info_eta_p[j];
    if (want_obs && J > 0) {
        VectorXd vv(1 + J);
        vv(0) = -res.score_wt_lambda;
        for (int j = 0; j < J; ++j) {
            const double e = eta_p[j];
            vv(1 + j) = (e > 0.0) ? 1.0 / (1.0 + std::exp(-e))
                                  : std::exp(e) / (1.0 + std::exp(e));
        }
        B.noalias() -= res.var_N * (vv * vv.transpose());
    }
    small.noalias() = Z.transpose() * B * Z;
    return res.log_lik;
}

// Field geometry for one site: the field-column loadings on eta_lambda and
// their global indices in the top block (field_start = d). Areal kinds load a
// single field unit (BYM2 two: v + w); the SPDE field loads the mesh nodes of
// the site's projection row (a 2D P1 FEM row has at most MAX_FIELD_LOAD
// nonzeros -- the barycentric weights of the triangle the site falls in).
static const int MAX_FIELD_LOAD = 4;
struct FieldGeom {
    int nfield;
    double load[MAX_FIELD_LOAD];
    int    idx[MAX_FIELD_LOAD];
};

// Areal field geometry for spatial unit u.
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

// SPDE field geometry for one site: read the nonzeros of A.row(site).
inline FieldGeom field_geom_spde(int d, const SpdeCtx& sp, int site) {
    FieldGeom g;
    g.nfield = 0;
    for (int k = 0; k < sp.n_mesh; ++k) {
        const double w = (*sp.A)(site, k);
        if (w != 0.0) {
            if (g.nfield >= MAX_FIELD_LOAD)
                Rcpp::stop("SPDE projection row has more than %d nonzeros; "
                           "raise MAX_FIELD_LOAD.", MAX_FIELD_LOAD);
            g.load[g.nfield] = w;
            g.idx[g.nfield]  = d + k;
            ++g.nfield;
        }
    }
    return g;
}

inline double field_offset(FieldKind kind, const VectorXd& field, int n_spatial,
                           int u, double a, double b) {
    if (kind == FieldKind::BYM2) return a * field(u) + b * field(n_spatial + u);
    return field(u);
}

// SPDE per-site offset: (A u)_site = sum_k A[site, k] u[k].
inline double field_offset_spde(const VectorXd& field, const SpdeCtx& sp,
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
                              const VectorXd& field,
                              const SpdeCtx* spde = nullptr) {
    if (kind == FieldKind::SPDE) {
        // 0.5 log|Q| - 0.5 u' Q u (full rank; the (2 pi)^{-n_mesh/2} constant is
        // grid-independent and dropped, as in the single-species SPDE path).
        return 0.5 * spde->log_det_Q
             - 0.5 * (field.transpose() * ((*spde->Q) * field))(0, 0);
    }
    if (kind == FieldKind::BYM2) {
        const VectorXd v = field.head(n_spatial);
        const VectorXd w = field.tail(n_spatial);
        return tulpaObs::nmix_bym2_log_prior(n_spatial, adj_row_ptr, adj_col_idx,
                                             n_neighbors, v, w);
    }
    if (kind == FieldKind::ICAR) {
        return tulpaObs::nmix_icar_log_prior(n_spatial, tau, adj_row_ptr,
                                             adj_col_idx, n_neighbors, field);
    }
    return tulpaObs::nmix_car_proper_log_prior(n_spatial, tau, rho, log_det_Q_rho,
                                               adj_row_ptr, adj_col_idx,
                                               n_neighbors, field);
}

// Add the field prior to the top-block gradient and Hessian (in-place). The top
// block has the (beta_lambda, beta_p, field) layout the single-species helpers
// expect (field_start = d = p_lam + p_p).
inline void add_field_prior(FieldKind kind, int p_lam, int p_p, int n_spatial,
                            double tau, double rho,
                            const Rcpp::IntegerVector& adj_row_ptr,
                            const Rcpp::IntegerVector& adj_col_idx,
                            const Rcpp::IntegerVector& n_neighbors,
                            const VectorXd& field,
                            VectorXd& grad_top, MatrixXd& T,
                            const SpdeCtx* spde = nullptr) {
    if (kind == FieldKind::SPDE) {
        // Field prior -Q u to the score, +Q to the negative-Hessian field block.
        const int d = p_lam + p_p;
        grad_top.segment(d, spde->n_mesh).noalias() -= (*spde->Q) * field;
        T.block(d, d, spde->n_mesh, spde->n_mesh).noalias() += (*spde->Q);
        return;
    }
    if (kind == FieldKind::BYM2) {
        const VectorXd v = field.head(n_spatial);
        const VectorXd w = field.tail(n_spatial);
        tulpaObs::nmix_add_bym2_prior_to_grad_and_H(p_lam, p_p, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors, v, w, grad_top, T);
    } else {
        tulpaObs::nmix_add_car_to_spatial_block(p_lam, p_p, n_spatial, tau, rho,
            adj_row_ptr, adj_col_idx, n_neighbors, field, grad_top, T);
    }
}

// Sum-to-zero centering of the (intrinsic) field component against the global
// abundance intercept. ICAR centres f; BYM2 centres v only (w is proper);
// CAR_proper is full-rank and needs none.
inline void center_field(FieldKind kind, int p_lam, int p_p, int n_spatial,
                         VectorXd& field) {
    // CAR_proper and SPDE are full rank: the (intercept, field-mean) direction
    // is identified by Q itself, so no sum-to-zero centering.
    if (kind == FieldKind::CAR_PROPER || kind == FieldKind::SPDE ||
        n_spatial <= 0) return;
    const int len = (kind == FieldKind::BYM2) ? 2 * n_spatial : n_spatial;
    VectorXd holder(p_lam + p_p + len);
    holder.setZero();
    holder.segment(p_lam + p_p, len) = field;
    tulpaObs::nmix_center_field(p_lam, p_p, n_spatial, holder);
    field = holder.segment(p_lam + p_p, len);
}

struct CommSpatialResult {
    VectorXd mu;          // d community means
    VectorXd field;       // n_field_total = nfield * n_spatial (ICAR/CAR: f; BYM2: [v; w])
    MatrixXd Sigma_l, Sigma_p;
    MatrixXd blup_l, blup_p;   // S x p_lam, S x p_p
    MatrixXd vcov_mu;          // d x d community-mean covariance (b- and field-folded)
    double log_marginal = R_NegInf;
    double log_lik = R_NegInf;
    bool   converged = false;
    int    n_iter = 0;
    double boundary_max = 0.0;
};

// One outer-grid-point Laplace-EM fit. `field` is laid out [f] (ICAR/CAR) or
// [v; w] (BYM2); `a`, `b` are the BYM2 loadings (unused for ICAR/CAR).
CommSpatialResult community_spatial_em(
    FieldKind kind,
    const NMixCommunityOracle& orc,
    const Map<MatrixXd>& Xlam,
    const std::vector<int>& map_site_to_unit,
    int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    double tau, double rho, double log_det_Q_rho, double a, double b,
    double r,
    const VectorXd& mu_init, const MatrixXd& Sigma_l_init,
    const MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    double sigma_beta, bool verbose,
    const SpdeCtx* spde = nullptr) {

    const int p_lam = orc.p_lam;
    const int p_p   = orc.p_p;
    const int d     = p_lam + p_p;
    const int S     = orc.n_groups;
    const int nfield = (kind == FieldKind::BYM2) ? 2 : 1;
    const int field_len = nfield * n_spatial;
    const int m = d + field_len;
    // No ridge on the community means (mu). The abundance intercept and the
    // field constant mode form an exactly flat (data + ICAR-prior) direction;
    // any ridge on mu is the ONLY curvature along it and drives the intercept to
    // its prior mean 0, with the field absorbing the level (then deleted by the
    // sum-to-zero centering). The single-species spatial path keeps the fixed
    // effects flat for the same reason -- mu is identified by the data once the
    // field is anchored sum-to-zero.
    (void) sigma_beta;

    CommSpatialResult out;
    VectorXd mu = mu_init;
    VectorXd field = VectorXd::Zero(field_len);
    std::vector<VectorXd> bvec(S, VectorXd::Zero(d));
    MatrixXd Sig_l = Sigma_l_init, Sig_p = Sigma_p_init;

    auto block_prec = [&](const MatrixXd& Sl, const MatrixXd& Sp) {
        MatrixXd P = MatrixXd::Zero(d, d);
        P.topLeftCorner(p_lam, p_lam) = nmix_safe_inverse(Sl);
        P.bottomRightCorner(p_p, p_p) = nmix_safe_inverse(Sp);
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
            for (const auto& rec : orc.sp_sites[s]) {
                auto go = site_geom(rec.site, field_);
                const FieldGeom& g = go.first;
                double l = site_blocks(rec, Xlam, p_lam, p_p, coef, go.second,
                                       g.load, g.nfield, r, false, false,
                                       grad_aug, small);
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
            for (const auto& rec : orc.sp_sites[s]) {
                auto go = site_geom(rec.site, field_);
                const FieldGeom& g = go.first;
                double bw = 0.0;
                ll += site_blocks(rec, Xlam, p_lam, p_p, coef, go.second,
                                  g.load, g.nfield, r, true, want_obs,
                                  grad_aug, small, &bw);
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
        add_field_prior(kind, p_lam, p_p, n_spatial, tau, rho,
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
        const MatrixXd P = block_prec(Sig_l, Sig_p);

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
            double step = 1.0;
            bool stepped = false;
            double obj_cur = ll_cur
                + field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                  adj_row_ptr, adj_col_idx, n_neighbors, field, spde);
            for (int s = 0; s < S; ++s) obj_cur -= 0.5 * bvec[s].dot(P * bvec[s]);
            double max_step = 0.0;
            for (int h = 0; h < 12; ++h) {
                VectorXd mu_try = mu + step * dmu;
                VectorXd field_try = field + step * dfield;
                std::vector<VectorXd> b_try(S);
                for (int s = 0; s < S; ++s) b_try[s] = bvec[s] + step * delta_b[s];
                const double obj_try = objective(mu_try, field_try, b_try, P);
                if (R_finite(obj_try) && obj_try >= obj_cur - 1e-10) {
                    mu = mu_try; field = field_try; bvec = b_try;
                    center_field(kind, p_lam, p_p, n_spatial, field);
                    double dmax = std::max(dmu.cwiseAbs().maxCoeff(),
                                           dfield.cwiseAbs().maxCoeff());
                    for (int s = 0; s < S; ++s)
                        dmax = std::max(dmax, delta_b[s].cwiseAbs().maxCoeff());
                    max_step = step * dmax;
                    stepped = true;
                    break;
                }
                step *= 0.5;
            }
            if (!stepped) break;
            if (max_step < inner_tol) break;
        }

        // ---- recompute Fisher blocks at the converged mode for the M-step ----
        assemble(mu, field, bvec, P, /*want_obs=*/false,
                 grad_top, T, grad_b, D, C, nullptr, nullptr);
        for (int s = 0; s < S; ++s) Dinv[s] = nmix_safe_inverse(D[s]);

        // ---- M-step: EM update of the community covariances ----
        MatrixXd Sl_new = MatrixXd::Zero(p_lam, p_lam);
        MatrixXd Sp_new = MatrixXd::Zero(p_p, p_p);
        for (int s = 0; s < S; ++s) {
            const VectorXd bl = bvec[s].head(p_lam);
            const VectorXd bp = bvec[s].tail(p_p);
            Sl_new += bl * bl.transpose() + Dinv[s].topLeftCorner(p_lam, p_lam);
            Sp_new += bp * bp.transpose() + Dinv[s].bottomRightCorner(p_p, p_p);
        }
        Sl_new /= (double)S;
        Sp_new /= (double)S;
        const double dSig = std::max((Sl_new - Sig_l).cwiseAbs().maxCoeff(),
                                     (Sp_new - Sig_p).cwiseAbs().maxCoeff());
        Sig_l = Sl_new;
        Sig_p = Sp_new;
        if (verbose) Rcpp::Rcout << "  em " << out.n_iter << " dSigma=" << dSig << "\n";
        if (dSig < tol_em) { out.converged = true; break; }
    }

    // ---- final pass: observed-info arrowhead, log-marginal, vcov_mu ----
    const MatrixXd P = block_prec(Sig_l, Sig_p);
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
        double lp = field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                    adj_row_ptr, adj_col_idx, n_neighbors, field, spde);
        double bquad = 0.0;
        for (int s = 0; s < S; ++s) bquad += bvec[s].dot(P * bvec[s]);
        loglik_marg = ll + lp
                      + 0.5 * (double)S * logdetP - 0.5 * bquad
                      - 0.5 * sum_logdet_D - 0.5 * ldM;

        // vcov_mu = top-left d-block of M^{-1}, constrained (sum field = 0) for
        // the rank-deficient intrinsic fields (ICAR / BYM2 v); CAR_proper is
        // full rank.
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
    out.Sigma_l = Sig_l;
    out.Sigma_p = Sig_p;
    out.blup_l = MatrixXd(S, p_lam);
    out.blup_p = MatrixXd(S, p_p);
    for (int s = 0; s < S; ++s) {
        out.blup_l.row(s) = bvec[s].head(p_lam).transpose();
        out.blup_p.row(s) = bvec[s].tail(p_p).transpose();
    }
    out.vcov_mu = ok ? vcov_mu : MatrixXd::Constant(d, d, R_NaN);
    out.log_marginal = ok ? loglik_marg : R_NegInf;
    out.log_lik = data_loglik(mu, field, bvec);
    out.boundary_max = boundary_max;
    return out;
}

// Resolve the oracle XPtr to a concrete NMixCommunityOracle (the pre-grouped,
// lgamma-cached data container its constructor built).
const NMixCommunityOracle& as_community_oracle(SEXP oracle) {
    Rcpp::XPtr<tulpa::REGroupOracle> base(oracle);
    NMixCommunityOracle* orcp = dynamic_cast<NMixCommunityOracle*>(base.get());
    if (orcp == nullptr)
        Rcpp::stop("oracle is not an NMixCommunityOracle.");
    return *orcp;
}

std::vector<int> resolve_map(const Rcpp::IntegerVector& map_R, int n_sites,
                             int n_spatial) {
    if ((int)map_R.size() != n_sites)
        Rcpp::stop("length(map_site_to_unit) must equal n_sites.");
    std::vector<int> map(n_sites);
    for (int s = 0; s < n_sites; ++s) {
        const int u = map_R[s] - 1;
        if (u < 0 || u >= n_spatial)
            Rcpp::stop("map_site_to_unit out of range [1, n_spatial].");
        map[s] = u;
    }
    return map;
}

// One outer-grid point: the field hyperparameters (tau / rho / log|Q(rho)| and
// the BYM2 mixing loadings a, b) plus the NB size r that the inner Laplace-EM
// conditions on. Each entry builds a plan of these in its own grid-nesting
// order; the driver below just walks the plan.
struct CommGridPoint { double tau, rho, log_det_Q_rho, a, b, r; };

// Shared outer-grid driver for the community spatial nested-Laplace-EM. Owns the
// oracle / map / init unpacking, the grid walk, progress plumbing, and the
// community_pack_grid output assembly; each field kind supplies only its
// FieldKind, its field length, and the per-point plan + theta_grid_out it
// filled. community_
// spatial_em already dispatches the field prior on FieldKind, so adding a field
// kind is a new plan-builder, not a new copy of this loop.
Rcpp::List run_community_spatial_grid(
    FieldKind kind, int field_len,
    SEXP oracle,
    const Rcpp::IntegerVector& map_site_to_unit_R,
    const Rcpp::NumericMatrix& X_lambda_R,
    int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const std::vector<CommGridPoint>& plan,
    const Rcpp::NumericMatrix& theta_grid_out,
    const Rcpp::NumericVector& mu_init,
    const Rcpp::NumericMatrix& Sigma_lambda_init,
    const Rcpp::NumericMatrix& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    double sigma_beta, bool verbose,
    bool progress, int progress_every, double progress_throttle,
    const std::string& progress_file) {

    const NMixCommunityOracle& orc = as_community_oracle(oracle);
    const int n_sites = X_lambda_R.nrow();
    const int p_lam = orc.p_lam, p_p = orc.p_p, d = p_lam + p_p;
    Map<MatrixXd> Xl(REAL(X_lambda_R), n_sites, p_lam);
    std::vector<int> map = resolve_map(map_site_to_unit_R, n_sites, n_spatial);
    if ((int)mu_init.size() != d)
        Rcpp::stop("mu_init length must equal p_lam + p_p.");
    if (Sigma_lambda_init.nrow() != p_lam || Sigma_lambda_init.ncol() != p_lam)
        Rcpp::stop("Sigma_lambda_init must be p_lam x p_lam.");
    if (Sigma_p_init.nrow() != p_p || Sigma_p_init.ncol() != p_p)
        Rcpp::stop("Sigma_p_init must be p_p x p_p.");
    VectorXd mu0 = Map<VectorXd>(REAL(mu_init), d);
    MatrixXd Sl0 = Map<MatrixXd>(REAL(Sigma_lambda_init), p_lam, p_lam);
    MatrixXd Sp0 = Map<MatrixXd>(REAL(Sigma_p_init), p_p, p_p);

    const int n_grid = (int)plan.size();
    std::vector<CommSpatialResult> results(n_grid);
    // outer-grid progress (tulpa#45)
    auto gp = tulpaObs::make_grid_progress("nmix-spatial", n_grid, progress,
                                           progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        const CommGridPoint& g = plan[k];
        results[k] = community_spatial_em(
            kind, orc, Xl, map, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors,
            g.tau, g.rho, g.log_det_Q_rho, g.a, g.b,
            g.r, mu0, Sl0, Sp0,
            max_iter_em, tol_em, inner_max, inner_tol, sigma_beta, verbose);
        if (gp) gp->tick();
    }
    if (gp) gp->finish();
    return tulpaObs::community_pack_grid(
        d, field_len, n_grid, results, theta_grid_out,
        &CommSpatialResult::Sigma_l, &CommSpatialResult::blup_l,
        "Sigma_lambda", "b_lambda", "p_lambda", p_lam, p_p, n_spatial,
        &CommSpatialResult::boundary_max);
}

}  // namespace

// ---------------------------------------------------------------------------
// Entry points (one per field kind). `oracle` is the XPtr from
// cpp_nmix_community_oracle(); the outer grid integrates the field hyperparameter
// (tau / sigma, rho) and, under NB, the size r (outermost axis). Each grid point
// is a cold-restart Laplace-EM (the (lambda, p) identifiability ridge shifts as
// the field absorbs different variance, so warm-starting confounds the inner
// Newton -- matching the single-species spatial path).
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_spatial_icar(
    SEXP oracle,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector r_grid,
    Rcpp::NumericVector mu_init,
    Rcpp::NumericMatrix Sigma_lambda_init,
    Rcpp::NumericMatrix Sigma_p_init,
    int max_iter_em = 100, double tol_em = 1e-4,
    int inner_max = 50, double inner_tol = 1e-6,
    double sigma_beta = 100.0, bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = "") {

    const int n_tau = tau_grid.size(), n_r = r_grid.size();
    const int n_grid = n_tau * n_r;
    Rcpp::NumericMatrix theta_grid_out(n_grid, 2);
    std::vector<CommGridPoint> plan(n_grid);
    int k = 0;
    for (int ri = 0; ri < n_r; ++ri)
        for (int t = 0; t < n_tau; ++t, ++k) {
            theta_grid_out(k, 0) = tau_grid[t];
            theta_grid_out(k, 1) = r_grid[ri];
            plan[k] = CommGridPoint{ tau_grid[t], /*rho=*/1.0, /*log_det_Q_rho=*/0.0,
                                     /*a=*/0.0, /*b=*/0.0, r_grid[ri] };
        }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("tau", "r");
    return run_community_spatial_grid(
        FieldKind::ICAR, /*field_len=*/n_spatial, oracle, map_site_to_unit_R,
        X_lambda_R, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, theta_grid_out, mu_init, Sigma_lambda_init, Sigma_p_init,
        max_iter_em, tol_em, inner_max, inner_tol, sigma_beta, verbose,
        progress, progress_every, progress_throttle, progress_file);
}

// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_spatial_car_proper(
    SEXP oracle,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector rho_grid,
    Rcpp::NumericVector log_det_Q_rho,   // one per rho grid point
    Rcpp::NumericVector r_grid,
    Rcpp::NumericVector mu_init,
    Rcpp::NumericMatrix Sigma_lambda_init,
    Rcpp::NumericMatrix Sigma_p_init,
    int max_iter_em = 100, double tol_em = 1e-4,
    int inner_max = 50, double inner_tol = 1e-6,
    double sigma_beta = 100.0, bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = "") {

    const int n_tau = tau_grid.size(), n_rho = rho_grid.size(), n_r = r_grid.size();
    if ((int)log_det_Q_rho.size() != n_rho)
        Rcpp::stop("length(log_det_Q_rho) must equal length(rho_grid).");
    const int n_grid = n_tau * n_rho * n_r;
    Rcpp::NumericMatrix theta_grid_out(n_grid, 3);
    std::vector<CommGridPoint> plan(n_grid);
    int k = 0;
    for (int ri = 0; ri < n_r; ++ri)
        for (int rh = 0; rh < n_rho; ++rh)
            for (int t = 0; t < n_tau; ++t, ++k) {
                theta_grid_out(k, 0) = tau_grid[t];
                theta_grid_out(k, 1) = rho_grid[rh];
                theta_grid_out(k, 2) = r_grid[ri];
                plan[k] = CommGridPoint{ tau_grid[t], rho_grid[rh], log_det_Q_rho[rh],
                                         /*a=*/0.0, /*b=*/0.0, r_grid[ri] };
            }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("tau", "rho", "r");
    return run_community_spatial_grid(
        FieldKind::CAR_PROPER, /*field_len=*/n_spatial, oracle, map_site_to_unit_R,
        X_lambda_R, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, theta_grid_out, mu_init, Sigma_lambda_init, Sigma_p_init,
        max_iter_em, tol_em, inner_max, inner_tol, sigma_beta, verbose,
        progress, progress_every, progress_throttle, progress_file);
}

// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_spatial_bym2(
    SEXP oracle,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector sigma_grid,
    Rcpp::NumericVector rho_grid,
    double scale_factor,
    Rcpp::NumericVector r_grid,
    Rcpp::NumericVector mu_init,
    Rcpp::NumericMatrix Sigma_lambda_init,
    Rcpp::NumericMatrix Sigma_p_init,
    int max_iter_em = 100, double tol_em = 1e-4,
    int inner_max = 50, double inner_tol = 1e-6,
    double sigma_beta = 100.0, bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = "") {

    const int n_sig = sigma_grid.size(), n_rho = rho_grid.size(), n_r = r_grid.size();
    const int n_grid = n_sig * n_rho * n_r;
    Rcpp::NumericMatrix theta_grid_out(n_grid, 3);
    std::vector<CommGridPoint> plan(n_grid);
    int k = 0;
    for (int ri = 0; ri < n_r; ++ri)
        for (int rh = 0; rh < n_rho; ++rh)
            for (int sg = 0; sg < n_sig; ++sg, ++k) {
                const double sigma = sigma_grid[sg], rho = rho_grid[rh];
                const double a = sigma * std::sqrt(rho / scale_factor);
                const double b = sigma * std::sqrt(1.0 - rho);
                theta_grid_out(k, 0) = sigma;
                theta_grid_out(k, 1) = rho;
                theta_grid_out(k, 2) = r_grid[ri];
                plan[k] = CommGridPoint{ /*tau=*/1.0, rho, /*log_det_Q_rho=*/0.0,
                                         a, b, r_grid[ri] };
            }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("sigma", "rho", "r");
    return run_community_spatial_grid(
        FieldKind::BYM2, /*field_len=*/2 * n_spatial, oracle, map_site_to_unit_R,
        X_lambda_R, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, theta_grid_out, mu_init, Sigma_lambda_init, Sigma_p_init,
        max_iter_em, tol_em, inner_max, inner_tol, sigma_beta, verbose,
        progress, progress_every, progress_throttle, progress_file);
}

// Continuous Matern (SPDE) shared field on the abundance arm. The field lives
// at n_mesh FEM nodes; the dense projection A (n_sites x n_mesh) maps mesh nodes
// onto sites, shared across species exactly as the areal field is. Per grid
// point the proper Matern precision Q(range, sigma) (carrying tau_spde^2) and
// its log|Q| are built once on the R side and passed in; the outer grid axes are
// (range, sigma) [, NB size r]. n_spatial is the mesh node count here.
// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_spatial_spde(
    SEXP oracle,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::NumericMatrix A_R,            // dense [n_sites x n_mesh]
    Rcpp::List Q_list,                  // per-grid-point precision [n_mesh x n_mesh]
    Rcpp::NumericVector log_det_Q,      // per-grid-point log|Q|
    Rcpp::NumericMatrix theta_grid_R,   // [n_grid x n_theta] (range, sigma, r)
    Rcpp::NumericVector r_grid,         // NB size per grid point (or +Inf)
    Rcpp::NumericVector mu_init,
    Rcpp::NumericMatrix Sigma_lambda_init,
    Rcpp::NumericMatrix Sigma_p_init,
    int max_iter_em = 100, double tol_em = 1e-4,
    int inner_max = 50, double inner_tol = 1e-6,
    double sigma_beta = 100.0, bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = "") {

    const NMixCommunityOracle& orc = as_community_oracle(oracle);
    const int n_sites = X_lambda_R.nrow();
    const int n_mesh  = A_R.ncol();
    const int p_lam = orc.p_lam, p_p = orc.p_p, d = p_lam + p_p;
    if (A_R.nrow() != n_sites)
        Rcpp::stop("nrow(A) must equal nrow(X_lambda).");
    Map<MatrixXd> Xl(REAL(X_lambda_R), n_sites, p_lam);
    const MatrixXd A = Map<MatrixXd>(REAL(A_R), n_sites, n_mesh);
    if ((int)mu_init.size() != d)
        Rcpp::stop("mu_init length must equal p_lam + p_p.");
    if (Sigma_lambda_init.nrow() != p_lam || Sigma_lambda_init.ncol() != p_lam)
        Rcpp::stop("Sigma_lambda_init must be p_lam x p_lam.");
    if (Sigma_p_init.nrow() != p_p || Sigma_p_init.ncol() != p_p)
        Rcpp::stop("Sigma_p_init must be p_p x p_p.");
    VectorXd mu0 = Map<VectorXd>(REAL(mu_init), d);
    MatrixXd Sl0 = Map<MatrixXd>(REAL(Sigma_lambda_init), p_lam, p_lam);
    MatrixXd Sp0 = Map<MatrixXd>(REAL(Sigma_p_init), p_p, p_p);

    const int n_grid = Q_list.size();
    if ((int)log_det_Q.size() != n_grid)
        Rcpp::stop("length(log_det_Q) must equal length(Q_list).");
    if (theta_grid_R.nrow() != n_grid)
        Rcpp::stop("nrow(theta_grid) must equal length(Q_list).");
    if ((int)r_grid.size() != n_grid)
        Rcpp::stop("length(r_grid) must equal length(Q_list).");

    // Empty CSR adjacency (unused on the SPDE path).
    Rcpp::IntegerVector empty_int(0);

    std::vector<CommSpatialResult> results(n_grid);
    std::vector<MatrixXd> Qmats(n_grid);   // keep alive for the SpdeCtx pointers
    // outer-grid progress (tulpa#45)
    auto gp = tulpaObs::make_grid_progress("nmix-spatial", n_grid, progress,
                                           progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        Rcpp::NumericMatrix Qk_R = Q_list[k];
        if (Qk_R.nrow() != n_mesh || Qk_R.ncol() != n_mesh)
            Rcpp::stop("Q_list[[%d]] must be n_mesh x n_mesh.", k + 1);
        Qmats[k] = Map<MatrixXd>(REAL(Qk_R), n_mesh, n_mesh);
        SpdeCtx sp;
        sp.A = &A; sp.Q = &Qmats[k]; sp.log_det_Q = log_det_Q[k]; sp.n_mesh = n_mesh;
        results[k] = community_spatial_em(
            FieldKind::SPDE, orc, Xl, /*map=*/std::vector<int>(), /*n_spatial=*/n_mesh,
            empty_int, empty_int, empty_int,
            /*tau=*/1.0, /*rho=*/1.0, /*log_det_Q_rho=*/0.0, /*a=*/0.0, /*b=*/0.0,
            r_grid[k], mu0, Sl0, Sp0,
            max_iter_em, tol_em, inner_max, inner_tol, sigma_beta, verbose, &sp);
        if (gp) gp->tick();
    }
    if (gp) gp->finish();
    return tulpaObs::community_pack_grid(
        d, n_mesh, n_grid, results, theta_grid_R,
        &CommSpatialResult::Sigma_l, &CommSpatialResult::blup_l,
        "Sigma_lambda", "b_lambda", "p_lambda", p_lam, p_p, n_mesh,
        &CommSpatialResult::boundary_max);
}
