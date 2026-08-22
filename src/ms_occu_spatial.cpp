// ms_occu_spatial.cpp
// Areal-spatial community single-season occupancy (ms_occu() + shared field;
// the occupancy analogue of sfMsNMix, tulpaObs#75). A per-species two-state
// occupancy model with Gaussian community hyperpriors on the per-species
// coefficients AND one shared ICAR / BYM2 / proper-CAR field on the OCCUPANCY
// arm:
//
//   logit psi_{s,i} = X_psi_i . (mu_psi + b_psi_s) + f_{u(i)}
//   logit p_{s,i}   = X_p_i   . (mu_p   + b_p_s)
//   b_psi_s ~ N(0, Sigma_psi),  b_p_s ~ N(0, Sigma_p),  one f shared
//
// The latent z integrates out per species-site in closed form (the occupancy
// two-state marginal, ms_occu_kernel.h, site-level detection). Per outer grid
// point a Laplace-EM iterates the joint (mu, f, {b_s}) mode-find (block-elim
// Newton, b_s Schur-folded) + a closed-form Sigma M-step, and R grid-integrates
// the community means / covariance / field over the field-hyperparameter
// posterior. The top block carries the (mu_psi, mu_p, f) layout the single-
// species areal helpers expect (field_start = d = p_psi + p_p), so the shared
// field-prior / centering / constrained-covariance helpers in
// nmix_spatial_kernel.h / nmix_spatial_kernel_bym2.h / nmix_linalg.h apply
// directly -- the ONLY occupancy-specific piece is the per-site cell (score +
// curvature), which comes from ms_occu_site_cell. This mirrors the count
// community-spatial driver (nmix_community_spatial.cpp) but over the occupancy
// marginal; community means stay flat (the intercept + field constant mode is
// an exactly flat direction the sum-to-zero centering of the intrinsic field
// resolves).

#include "tobs_math.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "ms_occu_kernel.h"
#include "nmix_spatial_kernel.h"
#include "nmix_spatial_kernel_bym2.h"
#include "nmix_linalg.h"
#include "community_grid_pack.h"
#include "newton_step.h"

using namespace Rcpp;
using tulpaObs::clamp_eta;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

namespace tulpaObs {
namespace {

enum class OccFieldKind { ICAR, CAR_PROPER, BYM2 };


// Per-(species) site records: the site index, occupancy / detection design rows
// are shared (site-level) across species, so only the summary (n_valid, n_det,
// any_det) is per species. Pre-grouped by species for the EM site loop.
struct OccCommData {
    int n_sites = 0, n_species = 0, p_psi = 0, p_p = 0;
    Map<MatrixXd> X_psi;     // n_sites x p_psi
    Map<MatrixXd> X_p;       // n_sites x p_p
    // per species: n_valid / n_det / any_det over sites (length n_sites)
    std::vector<std::vector<int>> n_valid, n_det;
    std::vector<std::vector<char>> any_det;
    std::vector<int> map_site_to_unit;   // length n_sites (0-based field unit)

    OccCommData(const NumericMatrix& Xpsi_R, const NumericMatrix& Xp_R)
        : X_psi(REAL(Xpsi_R), Xpsi_R.nrow(), Xpsi_R.ncol()),
          X_p(REAL(Xp_R), Xp_R.nrow(), Xp_R.ncol()) {}
};

// Field geometry for one site: loadings on eta_psi + global field indices.
static const int OCC_MAX_FIELD_LOAD = 4;
struct OccFieldGeom { int nfield; double load[OCC_MAX_FIELD_LOAD]; int idx[OCC_MAX_FIELD_LOAD]; };

inline OccFieldGeom occ_field_geom(OccFieldKind kind, int d, int n_spatial, int u,
                                   double a, double b) {
    OccFieldGeom g;
    if (kind == OccFieldKind::BYM2) {
        g.nfield = 2;
        g.load[0] = a; g.idx[0] = d + u;
        g.load[1] = b; g.idx[1] = d + n_spatial + u;
    } else {
        g.nfield = 1;
        g.load[0] = 1.0; g.idx[0] = d + u;
    }
    return g;
}

inline double occ_field_offset(OccFieldKind kind, const VectorXd& field,
                               int n_spatial, int u, double a, double b) {
    if (kind == OccFieldKind::BYM2) return a * field(u) + b * field(n_spatial + u);
    return field(u);
}

// One site's contribution: writes grad_aug (length na = d + nfield) and, when
// want_block, the curvature `small` (na x na) via the Z' B Z sandwich. The
// occupancy arm (eta_psi) carries the field loadings + X_psi; the detection arm
// (eta_p) carries X_p. Mirrors site_blocks in nmix_community_spatial.cpp but
// over the occupancy cell. Returns the per-site marginal log-lik.
inline double occ_site_blocks(const OccCommData& d_, int s, int site,
                              const VectorXd& coef, double field_offset,
                              const double* load, int nfield,
                              bool want_block, bool want_obs,
                              VectorXd& grad_aug, MatrixXd& small) {
    const int p_psi = d_.p_psi, p_p = d_.p_p, d = p_psi + p_p;
    const int na = d + nfield;
    grad_aug.setZero(na);
    if (want_block) small.setZero(na, na);
    const int nv = d_.n_valid[s][site];
    if (nv == 0) return 0.0;

    double eta_psi = field_offset;
    for (int c = 0; c < p_psi; ++c) eta_psi += d_.X_psi(site, c) * coef(c);
    eta_psi = clamp_eta(eta_psi);
    double eta_p = 0.0;
    for (int c = 0; c < p_p; ++c) eta_p += d_.X_p(site, c) * coef(p_psi + c);
    eta_p = clamp_eta(eta_p);

    const MsOccuSiteCell cell = ms_occu_site_cell(
        eta_psi, eta_p, nv, d_.n_det[s][site], d_.any_det[s][site] != 0, want_obs);

    // gradient: eta_psi -> mu_psi cols (X_psi) + field cols (loadings);
    //           eta_p   -> mu_p cols (X_p).
    for (int c = 0; c < p_psi; ++c) grad_aug(c) += d_.X_psi(site, c) * cell.g_psi;
    for (int f = 0; f < nfield; ++f) grad_aug(d + f) += load[f] * cell.g_psi;
    for (int c = 0; c < p_p; ++c) grad_aug(p_psi + c) += d_.X_p(site, c) * cell.g_p;

    if (!want_block) return cell.log_lik;

    // Augmented design Z (2 x na): row 0 = psi (mu_psi + field), row 1 = p.
    MatrixXd Z = MatrixXd::Zero(2, na);
    for (int c = 0; c < p_psi; ++c) Z(0, c) = d_.X_psi(site, c);
    for (int f = 0; f < nfield; ++f) Z(0, d + f) = load[f];
    for (int c = 0; c < p_p; ++c) Z(1, p_psi + c) = d_.X_p(site, c);

    // Per-site eta-space 2x2 negative-Hessian B.
    MatrixXd B(2, 2);
    B(0, 0) = cell.B_pp_psi;
    B(1, 1) = cell.B_p_p;
    B(0, 1) = cell.B_cross;
    B(1, 0) = cell.B_cross;
    small.noalias() = Z.transpose() * B * Z;
    return cell.log_lik;
}

struct OccCommResult {
    VectorXd mu;
    VectorXd field;
    MatrixXd Sigma_psi, Sigma_p;
    MatrixXd blup_psi, blup_p;
    MatrixXd vcov_mu;
    double log_marginal = R_NegInf;
    double log_lik = R_NegInf;
    bool converged = false;
    int n_iter = 0;
};

inline double occ_field_log_prior(OccFieldKind kind, int n_spatial, double tau,
                                  double rho, double log_det_Q_rho,
                                  const IntegerVector& adj_row_ptr,
                                  const IntegerVector& adj_col_idx,
                                  const IntegerVector& n_neighbors,
                                  const VectorXd& field) {
    if (kind == OccFieldKind::BYM2) {
        const VectorXd v = field.head(n_spatial);
        const VectorXd w = field.tail(n_spatial);
        return nmix_bym2_log_prior(n_spatial, adj_row_ptr, adj_col_idx,
                                   n_neighbors, v, w);
    }
    if (kind == OccFieldKind::ICAR)
        return nmix_icar_log_prior(n_spatial, tau, adj_row_ptr, adj_col_idx,
                                   n_neighbors, field);
    return nmix_car_proper_log_prior(n_spatial, tau, rho, log_det_Q_rho,
                                     adj_row_ptr, adj_col_idx, n_neighbors, field);
}

// Add the field prior to the top-block gradient + Hessian (field_start = d).
// ICAR uses rho = 1 (intrinsic Q); proper-CAR uses the grid rho; BYM2 splits the
// (v, w) blocks. Mirrors add_field_prior in nmix_community_spatial.cpp.
inline void occ_add_field_prior(OccFieldKind kind, int p_psi, int p_p,
                                int n_spatial, double tau, double rho,
                                const IntegerVector& adj_row_ptr,
                                const IntegerVector& adj_col_idx,
                                const IntegerVector& n_neighbors,
                                const VectorXd& field, VectorXd& grad_top,
                                MatrixXd& T) {
    if (kind == OccFieldKind::BYM2) {
        const VectorXd v = field.head(n_spatial);
        const VectorXd w = field.tail(n_spatial);
        nmix_add_bym2_prior_to_grad_and_H(p_psi, p_p, n_spatial, adj_row_ptr,
                                          adj_col_idx, n_neighbors, v, w,
                                          grad_top, T);
        return;
    }
    const double rho_use = (kind == OccFieldKind::CAR_PROPER) ? rho : 1.0;
    nmix_add_car_to_spatial_block(p_psi, p_p, n_spatial, tau, rho_use,
                                  adj_row_ptr, adj_col_idx, n_neighbors, field,
                                  grad_top, T);
}

// Sum-to-zero centering of the intrinsic field against the global occupancy
// intercept. ICAR centres f; BYM2 centres v only (w proper); CAR_proper is full
// rank and needs none.
inline void occ_center_field(OccFieldKind kind, int p_psi, int p_p, int n_spatial,
                             VectorXd& field) {
    if (kind == OccFieldKind::CAR_PROPER || n_spatial <= 0) return;
    const int len = field.size();
    VectorXd holder(p_psi + p_p + len);
    holder.setZero();
    holder.segment(p_psi + p_p, len) = field;
    nmix_center_field(p_psi, p_p, n_spatial, holder);
    field = holder.segment(p_psi + p_p, len);
}

// One outer-grid-point Laplace-EM fit over the occupancy marginal.
OccCommResult occ_community_spatial_em(
    OccFieldKind kind, const OccCommData& d_, int n_spatial,
    const IntegerVector& adj_row_ptr, const IntegerVector& adj_col_idx,
    const IntegerVector& n_neighbors,
    double tau, double rho, double log_det_Q_rho, double a, double b,
    const VectorXd& mu_init, const MatrixXd& Sigma_psi_init,
    const MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose) {

    const int p_psi = d_.p_psi, p_p = d_.p_p, d = p_psi + p_p, S = d_.n_species;
    const int nfield = (kind == OccFieldKind::BYM2) ? 2 : 1;
    const int field_len = nfield * n_spatial;
    const int m = d + field_len;

    OccCommResult out;
    VectorXd mu = mu_init;
    VectorXd field = VectorXd::Zero(field_len);
    std::vector<VectorXd> bvec(S, VectorXd::Zero(d));
    MatrixXd Sig_psi = Sigma_psi_init, Sig_p = Sigma_p_init;

    auto block_prec = [&](const MatrixXd& Sl, const MatrixXd& Sp) {
        MatrixXd P = MatrixXd::Zero(d, d);
        P.topLeftCorner(p_psi, p_psi) = nmix_safe_inverse(Sl);
        P.bottomRightCorner(p_p, p_p) = nmix_safe_inverse(Sp);
        return P;
    };
    auto site_geom = [&](int site, const VectorXd& field_) {
        const int u = d_.map_site_to_unit[site];
        return std::make_pair(occ_field_geom(kind, d, n_spatial, u, a, b),
                              occ_field_offset(kind, field_, n_spatial, u, a, b));
    };
    auto data_loglik = [&](const VectorXd& mu_, const VectorXd& field_,
                           const std::vector<VectorXd>& b_) {
        double ll = 0.0; VectorXd grad_aug; MatrixXd small;
        for (int s = 0; s < S; ++s) {
            VectorXd coef = mu_ + b_[s];
            for (int site = 0; site < d_.n_sites; ++site) {
                auto go = site_geom(site, field_);
                double l = occ_site_blocks(d_, s, site, coef, go.second,
                                           go.first.load, go.first.nfield,
                                           false, false, grad_aug, small);
                if (!R_finite(l)) return R_NegInf;
                ll += l;
            }
        }
        return ll;
    };
    auto objective = [&](const VectorXd& mu_, const VectorXd& field_,
                         const std::vector<VectorXd>& b_, const MatrixXd& P) {
        double obj = data_loglik(mu_, field_, b_);
        if (!R_finite(obj)) return R_NegInf;
        obj += occ_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                   adj_row_ptr, adj_col_idx, n_neighbors, field_);
        for (int s = 0; s < S; ++s) obj -= 0.5 * b_[s].dot(P * b_[s]);
        return obj;
    };
    auto assemble = [&](const VectorXd& mu_, const VectorXd& field_,
                        const std::vector<VectorXd>& b_, const MatrixXd& P,
                        bool want_obs, VectorXd& grad_top, MatrixXd& T,
                        std::vector<VectorXd>& grad_b, std::vector<MatrixXd>& D,
                        std::vector<MatrixXd>& C, double* log_lik_out) {
        grad_top.setZero(m); T.setZero(m, m);
        double ll = 0.0; VectorXd grad_aug; MatrixXd small;
        for (int s = 0; s < S; ++s) {
            grad_b[s].setZero(d); D[s].setZero(d, d); C[s].setZero(m, d);
            VectorXd coef = mu_ + b_[s];
            for (int site = 0; site < d_.n_sites; ++site) {
                auto go = site_geom(site, field_);
                const OccFieldGeom& g = go.first;
                ll += occ_site_blocks(d_, s, site, coef, go.second,
                                      g.load, g.nfield, true, want_obs,
                                      grad_aug, small);
                grad_b[s].head(d) += grad_aug.head(d);
                grad_top.head(d)  += grad_aug.head(d);
                D[s].noalias() += small.topLeftCorner(d, d);
                T.topLeftCorner(d, d).noalias() += small.topLeftCorner(d, d);
                C[s].topLeftCorner(d, d).noalias() += small.topLeftCorner(d, d);
                for (int f = 0; f < g.nfield; ++f) {
                    const int gi = g.idx[f];
                    grad_top(gi) += grad_aug(d + f);
                    for (int c = 0; c < d; ++c) {
                        const double cross = small(d + f, c);
                        T(c, gi) += cross; T(gi, c) += cross;
                        C[s](gi, c) += cross;
                    }
                    for (int f2 = 0; f2 < g.nfield; ++f2)
                        T(gi, g.idx[f2]) += small(d + f, d + f2);
                }
            }
            D[s].noalias() += P;
            grad_b[s].noalias() -= P * b_[s];
        }
        occ_add_field_prior(kind, p_psi, p_p, n_spatial, tau, rho,
                            adj_row_ptr, adj_col_idx, n_neighbors, field_,
                            grad_top, T);
        if (log_lik_out) *log_lik_out = ll;
    };

    VectorXd grad_top(m); MatrixXd T(m, m);
    std::vector<VectorXd> grad_b(S, VectorXd::Zero(d));
    std::vector<MatrixXd> D(S, MatrixXd::Zero(d, d)), C(S, MatrixXd::Zero(m, d));
    std::vector<MatrixXd> Dinv(S, MatrixXd::Zero(d, d));

    for (int em_it = 0; em_it < max_iter_em; ++em_it) {
        out.n_iter = em_it + 1;
        const MatrixXd P = block_prec(Sig_psi, Sig_p);

        for (int nit = 0; nit < inner_max; ++nit) {
            double ll_cur = 0.0;
            assemble(mu, field, bvec, P, false, grad_top, T, grad_b, D, C, &ll_cur);
            MatrixXd M = T; VectorXd rhs = grad_top;
            for (int s = 0; s < S; ++s) {
                Dinv[s] = nmix_safe_inverse(D[s]);
                const MatrixXd CD = C[s] * Dinv[s];
                M.noalias()   -= CD * C[s].transpose();
                rhs.noalias() -= CD * grad_b[s];
            }
            nmix_add_diagonal_ridge(M);
            const VectorXd delta_top = nmix_safe_inverse(M) * rhs;
            std::vector<VectorXd> delta_b(S, VectorXd::Zero(d));
            for (int s = 0; s < S; ++s)
                delta_b[s] = Dinv[s] * (grad_b[s] - C[s].transpose() * delta_top);
            const VectorXd dmu = delta_top.head(d);
            const VectorXd dfield = delta_top.segment(d, field_len);

            double obj_cur = ll_cur
                + occ_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                      adj_row_ptr, adj_col_idx, n_neighbors, field);
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
                    occ_center_field(kind, p_psi, p_p, n_spatial, field);
                    double dmax = std::max(dmu.cwiseAbs().maxCoeff(),
                                           dfield.cwiseAbs().maxCoeff());
                    for (int s = 0; s < S; ++s)
                        dmax = std::max(dmax, delta_b[s].cwiseAbs().maxCoeff());
                    max_step = step * dmax;
                });
            if (!stepped) break;
            if (max_step < inner_tol) break;
        }

        assemble(mu, field, bvec, P, false, grad_top, T, grad_b, D, C, nullptr);
        for (int s = 0; s < S; ++s) Dinv[s] = nmix_safe_inverse(D[s]);

        MatrixXd Spsi_new = MatrixXd::Zero(p_psi, p_psi);
        MatrixXd Sp_new   = MatrixXd::Zero(p_p, p_p);
        for (int s = 0; s < S; ++s) {
            const VectorXd bpsi = bvec[s].head(p_psi);
            const VectorXd bp   = bvec[s].tail(p_p);
            Spsi_new += bpsi * bpsi.transpose() + Dinv[s].topLeftCorner(p_psi, p_psi);
            Sp_new   += bp   * bp.transpose()   + Dinv[s].bottomRightCorner(p_p, p_p);
        }
        Spsi_new /= (double) S; Sp_new /= (double) S;
        const double dSig = std::max((Spsi_new - Sig_psi).cwiseAbs().maxCoeff(),
                                     (Sp_new   - Sig_p).cwiseAbs().maxCoeff());
        Sig_psi = Spsi_new; Sig_p = Sp_new;
        if (verbose) Rcpp::Rcout << "  occ-em " << out.n_iter << " dSigma=" << dSig << "\n";
        if (dSig < tol_em) { out.converged = true; break; }
    }

    const MatrixXd P = block_prec(Sig_psi, Sig_p);
    const double logdetP = nmix_logdet_spd(P);

    auto final_assemble = [&](bool want_obs, double& loglik_marg,
                              MatrixXd& vcov_mu) -> bool {
        VectorXd gtop(m); MatrixXd Tt(m, m);
        std::vector<VectorXd> gb(S, VectorXd::Zero(d));
        std::vector<MatrixXd> Dd(S, MatrixXd::Zero(d, d)), Cc(S, MatrixXd::Zero(m, d));
        double ll = 0.0;
        assemble(mu, field, bvec, P, want_obs, gtop, Tt, gb, Dd, Cc, &ll);
        double sum_logdet_D = 0.0; MatrixXd M = Tt;
        for (int s = 0; s < S; ++s) {
            const double ldD = nmix_logdet_spd(Dd[s]);
            if (!R_finite(ldD)) return false;
            sum_logdet_D += ldD;
            M.noalias() -= Cc[s] * nmix_safe_inverse(Dd[s]) * Cc[s].transpose();
        }
        MatrixXd M_det = M; nmix_add_diagonal_ridge(M_det);
        const double ldM = nmix_logdet_spd(M_det);
        if (!R_finite(ldM)) return false;
        double lp = occ_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                        adj_row_ptr, adj_col_idx, n_neighbors, field);
        double bquad = 0.0;
        for (int s = 0; s < S; ++s) bquad += bvec[s].dot(P * bvec[s]);
        loglik_marg = ll + lp + 0.5 * (double) S * logdetP - 0.5 * bquad
                    - 0.5 * sum_logdet_D - 0.5 * ldM;
        const bool constrain = (kind != OccFieldKind::CAR_PROPER);
        MatrixXd cov_top = nmix_constrained_top_cov(M, m, d, d, field_len, constrain);
        if (!cov_top.allFinite()) return false;
        vcov_mu = cov_top;
        return true;
    };

    double loglik_marg = R_NegInf; MatrixXd vcov_mu;
    bool ok = final_assemble(true, loglik_marg, vcov_mu);
    if (!ok) ok = final_assemble(false, loglik_marg, vcov_mu);

    out.mu = mu; out.field = field;
    out.Sigma_psi = Sig_psi; out.Sigma_p = Sig_p;
    out.blup_psi = MatrixXd(S, p_psi); out.blup_p = MatrixXd(S, p_p);
    for (int s = 0; s < S; ++s) {
        out.blup_psi.row(s) = bvec[s].head(p_psi).transpose();
        out.blup_p.row(s)   = bvec[s].tail(p_p).transpose();
    }
    out.vcov_mu = ok ? vcov_mu : MatrixXd::Constant(d, d, R_NaN);
    out.log_marginal = ok ? loglik_marg : R_NegInf;
    out.log_lik = data_loglik(mu, field, bvec);
    return out;
}

// Build the occupancy community data container from R inputs.
OccCommData build_occ_data(const NumericMatrix& X_psi_R, const NumericMatrix& X_p_R,
                           const IntegerMatrix& n_valid_R, const IntegerMatrix& n_det_R,
                           const IntegerVector& map_R, int n_spatial) {
    OccCommData d_(X_psi_R, X_p_R);
    d_.n_sites   = X_psi_R.nrow();
    d_.p_psi     = X_psi_R.ncol();
    d_.p_p       = X_p_R.ncol();
    d_.n_species = n_valid_R.ncol();
    if (n_valid_R.nrow() != d_.n_sites || n_det_R.nrow() != d_.n_sites ||
        n_det_R.ncol() != d_.n_species)
        Rcpp::stop("n_valid / n_det must be n_sites x n_species.");
    d_.n_valid.assign(d_.n_species, std::vector<int>(d_.n_sites, 0));
    d_.n_det.assign(d_.n_species, std::vector<int>(d_.n_sites, 0));
    d_.any_det.assign(d_.n_species, std::vector<char>(d_.n_sites, 0));
    for (int s = 0; s < d_.n_species; ++s)
        for (int i = 0; i < d_.n_sites; ++i) {
            d_.n_valid[s][i] = n_valid_R(i, s);
            d_.n_det[s][i]   = n_det_R(i, s);
            d_.any_det[s][i] = (n_det_R(i, s) > 0) ? 1 : 0;
        }
    if ((int)map_R.size() != d_.n_sites)
        Rcpp::stop("length(map_site_to_unit) must equal n_sites.");
    d_.map_site_to_unit.assign(d_.n_sites, 0);
    for (int i = 0; i < d_.n_sites; ++i) {
        const int u = map_R[i] - 1;
        if (u < 0 || u >= n_spatial)
            Rcpp::stop("map_site_to_unit out of range [1, n_spatial].");
        d_.map_site_to_unit[i] = u;
    }
    return d_;
}

struct OccGridPoint { double tau, rho, log_det_Q_rho, a, b; };

Rcpp::List run_occ_spatial_grid(
    OccFieldKind kind, int field_len,
    const NumericMatrix& X_psi_R, const NumericMatrix& X_p_R,
    const IntegerMatrix& n_valid_R, const IntegerMatrix& n_det_R,
    const IntegerVector& map_R, int n_spatial,
    const IntegerVector& adj_row_ptr, const IntegerVector& adj_col_idx,
    const IntegerVector& n_neighbors,
    const std::vector<OccGridPoint>& plan, const NumericMatrix& theta_grid_out,
    const NumericVector& mu_init, const NumericMatrix& Sigma_psi_init,
    const NumericMatrix& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose) {

    OccCommData d_ = build_occ_data(X_psi_R, X_p_R, n_valid_R, n_det_R, map_R, n_spatial);
    const int p_psi = d_.p_psi, p_p = d_.p_p, d = p_psi + p_p;
    if ((int)mu_init.size() != d)
        Rcpp::stop("mu_init length must equal p_psi + p_p.");
    if (Sigma_psi_init.nrow() != p_psi || Sigma_psi_init.ncol() != p_psi)
        Rcpp::stop("Sigma_psi_init must be p_psi x p_psi.");
    if (Sigma_p_init.nrow() != p_p || Sigma_p_init.ncol() != p_p)
        Rcpp::stop("Sigma_p_init must be p_p x p_p.");
    VectorXd mu0 = Map<VectorXd>(REAL(mu_init), d);
    MatrixXd Spsi0 = Map<MatrixXd>(REAL(Sigma_psi_init), p_psi, p_psi);
    MatrixXd Sp0   = Map<MatrixXd>(REAL(Sigma_p_init), p_p, p_p);

    const int n_grid = (int) plan.size();
    std::vector<OccCommResult> results(n_grid);
    for (int k = 0; k < n_grid; ++k) {
        const OccGridPoint& g = plan[k];
        results[k] = occ_community_spatial_em(
            kind, d_, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
            g.tau, g.rho, g.log_det_Q_rho, g.a, g.b,
            mu0, Spsi0, Sp0, max_iter_em, tol_em, inner_max, inner_tol, verbose);
    }
    return tulpaObs::community_pack_grid(
        d, field_len, n_grid, results, theta_grid_out,
        &OccCommResult::Sigma_psi, &OccCommResult::blup_psi,
        "Sigma_psi", "b_psi", "p_psi", p_psi, p_p, n_spatial);
}

}  // namespace
}  // namespace tulpaObs


// ---------------------------------------------------------------------------
// Cross-check entry: per-site occupancy cell score + curvature vs the R oracle
// (FD-validated). Returns log-lik, score (g_psi, g_p), and the 2x2 -Hessian.
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_site_cell(double eta_psi, double eta_p, int n_valid,
                                 int n_det, bool observed) {
    const tulpaObs::MsOccuSiteCell c = tulpaObs::ms_occu_site_cell(
        eta_psi, eta_p, n_valid, n_det, n_det > 0, observed);
    Rcpp::NumericMatrix B(2, 2);
    B(0, 0) = c.B_pp_psi; B(1, 1) = c.B_p_p;
    B(0, 1) = c.B_cross;  B(1, 0) = c.B_cross;
    return Rcpp::List::create(
        Rcpp::Named("log_lik") = c.log_lik,
        Rcpp::Named("grad") = Rcpp::NumericVector::create(c.g_psi, c.g_p),
        Rcpp::Named("neg_hess") = B);
}

// ICAR
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_spatial_icar(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p,
    Rcpp::IntegerMatrix n_valid, Rcpp::IntegerMatrix n_det,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector mu_init, Rcpp::NumericMatrix Sigma_psi_init,
    Rcpp::NumericMatrix Sigma_p_init, int max_iter_em, bool verbose) {
    std::vector<tulpaObs::OccGridPoint> plan;
    Rcpp::NumericMatrix tg(tau_grid.size(), 1);
    for (int i = 0; i < tau_grid.size(); ++i) {
        plan.push_back({tau_grid[i], 1.0, 0.0, 0.0, 0.0});
        tg(i, 0) = tau_grid[i];
    }
    colnames(tg) = Rcpp::CharacterVector::create("tau");
    return tulpaObs::run_occ_spatial_grid(
        tulpaObs::OccFieldKind::ICAR, n_spatial, X_psi, X_p, n_valid, n_det,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, mu_init, Sigma_psi_init, Sigma_p_init,
        max_iter_em, 1e-4, 50, 1e-7, verbose);
}

// Proper CAR
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_spatial_car_proper(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p,
    Rcpp::IntegerMatrix n_valid, Rcpp::IntegerMatrix n_det,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector rho_grid, Rcpp::NumericVector log_det_Q_rho,
    Rcpp::NumericVector mu_init, Rcpp::NumericMatrix Sigma_psi_init,
    Rcpp::NumericMatrix Sigma_p_init, int max_iter_em, bool verbose) {
    std::vector<tulpaObs::OccGridPoint> plan;
    std::vector<double> tg_tau, tg_rho;
    for (int ri = 0; ri < rho_grid.size(); ++ri)
        for (int ti = 0; ti < tau_grid.size(); ++ti) {
            plan.push_back({tau_grid[ti], rho_grid[ri], log_det_Q_rho[ri], 0.0, 0.0});
            tg_tau.push_back(tau_grid[ti]); tg_rho.push_back(rho_grid[ri]);
        }
    Rcpp::NumericMatrix tg(plan.size(), 2);
    for (std::size_t i = 0; i < plan.size(); ++i) { tg(i, 0) = tg_tau[i]; tg(i, 1) = tg_rho[i]; }
    colnames(tg) = Rcpp::CharacterVector::create("tau", "rho");
    return tulpaObs::run_occ_spatial_grid(
        tulpaObs::OccFieldKind::CAR_PROPER, n_spatial, X_psi, X_p, n_valid, n_det,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, mu_init, Sigma_psi_init, Sigma_p_init,
        max_iter_em, 1e-4, 50, 1e-7, verbose);
}

// BYM2
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_spatial_bym2(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p,
    Rcpp::IntegerMatrix n_valid, Rcpp::IntegerMatrix n_det,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector sigma_grid,
    Rcpp::NumericVector rho_grid, double scale_factor,
    Rcpp::NumericVector mu_init, Rcpp::NumericMatrix Sigma_psi_init,
    Rcpp::NumericMatrix Sigma_p_init, int max_iter_em, bool verbose) {
    std::vector<tulpaObs::OccGridPoint> plan;
    std::vector<double> tg_sig, tg_rho;
    for (int ri = 0; ri < rho_grid.size(); ++ri)
        for (int si = 0; si < sigma_grid.size(); ++si) {
            const double sigma = sigma_grid[si], rho = rho_grid[ri];
            const double a = sigma * std::sqrt(rho / scale_factor);
            const double b = sigma * std::sqrt(1.0 - rho);
            plan.push_back({1.0, rho, 0.0, a, b});
            tg_sig.push_back(sigma); tg_rho.push_back(rho);
        }
    Rcpp::NumericMatrix tg(plan.size(), 2);
    for (std::size_t i = 0; i < plan.size(); ++i) { tg(i, 0) = tg_sig[i]; tg(i, 1) = tg_rho[i]; }
    colnames(tg) = Rcpp::CharacterVector::create("sigma", "rho");
    return tulpaObs::run_occ_spatial_grid(
        tulpaObs::OccFieldKind::BYM2, 2 * n_spatial, X_psi, X_p, n_valid, n_det,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, mu_init, Sigma_psi_init, Sigma_p_init,
        max_iter_em, 1e-4, 50, 1e-7, verbose);
}
