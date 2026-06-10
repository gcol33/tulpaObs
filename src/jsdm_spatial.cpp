// jsdm_spatial.cpp
// Areal-spatial joint species distribution model (jsdm() + shared field;
// tulpaObs#76). The JSDM observes presence/absence directly (no detection
// process), with shared fixed-effect coefficients, a scalar per-species random
// intercept, and one shared ICAR / BYM2 / proper-CAR areal field on the latent
// occupancy:
//
//   logit psi_{s,i} = X_i . beta + b_s + f_{u(i)}
//   b_s ~ N(0, sigma_re^2),   one f shared across species
//
// Per outer grid point a Laplace-EM iterates the joint (beta, f, {b_s}) mode-find
// (block-elim Newton, the scalar b_s Schur-folded) + a closed-form sigma_re^2
// M-step, and R grid-integrates the fixed effects / their covariance / the field
// over the field-hyperparameter posterior (law of total covariance). One spatial
// unit per site. The top block carries the (beta, f) layout the single-species
// areal helpers expect (field_start = p_occ), so the shared field-prior /
// centering / constrained-covariance helpers in nmix_spatial_kernel.h /
// nmix_spatial_kernel_bym2.h / nmix_linalg.h apply directly -- the ONLY
// JSDM-specific piece is the per-(species, site) Bernoulli cell (jsdm_site_cell).
// This mirrors the community-occupancy areal driver (ms_occu_spatial.cpp) but
// over a single occupancy arm with a scalar species random intercept.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "jsdm_spatial_kernel.h"
#include "nmix_spatial_kernel.h"
#include "nmix_spatial_kernel_bym2.h"
#include "nmix_linalg.h"

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

namespace tulpaObs {
namespace {

enum class JsdmFieldKind { ICAR, CAR_PROPER, BYM2 };

inline double clamp30(double x) { return x < -30.0 ? -30.0 : (x > 30.0 ? 30.0 : x); }

// JSDM community data: a shared site-level occupancy design X (n_sites x p_occ),
// the per-species presence/absence matrix y (n_sites x n_species), and the
// site -> spatial-unit map.
struct JsdmCommData {
    int n_sites = 0, n_species = 0, p_occ = 0;
    Map<MatrixXd> X;                       // n_sites x p_occ
    std::vector<std::vector<int>> y;       // [species][site] presence/absence
    std::vector<int> map_site_to_unit;     // length n_sites (0-based field unit)

    explicit JsdmCommData(const NumericMatrix& X_R)
        : X(REAL(X_R), X_R.nrow(), X_R.ncol()) {}
};

// Field geometry for one site: loadings on eta + global field indices. The top
// layout is [beta (p_occ), field], so the field starts at column p_occ.
static const int JSDM_MAX_FIELD_LOAD = 2;
struct JsdmFieldGeom { int nfield; double load[JSDM_MAX_FIELD_LOAD]; int idx[JSDM_MAX_FIELD_LOAD]; };

inline JsdmFieldGeom jsdm_field_geom(JsdmFieldKind kind, int p_occ, int n_spatial,
                                     int u, double a, double b) {
    JsdmFieldGeom g;
    if (kind == JsdmFieldKind::BYM2) {
        g.nfield = 2;
        g.load[0] = a; g.idx[0] = p_occ + u;
        g.load[1] = b; g.idx[1] = p_occ + n_spatial + u;
    } else {
        g.nfield = 1;
        g.load[0] = 1.0; g.idx[0] = p_occ + u;
    }
    return g;
}

inline double jsdm_field_offset(JsdmFieldKind kind, const VectorXd& field,
                                int n_spatial, int u, double a, double b) {
    if (kind == JsdmFieldKind::BYM2) return a * field(u) + b * field(n_spatial + u);
    return field(u);
}

// One (species, site) contribution. The scalar species random intercept b_s
// loads onto eta with coefficient 1; the field loads with `load`; the fixed
// effects load with X_i. Writes grad_top (length m = p_occ + field_len), the top
// curvature `Ttop` (m x m, when want_block), the scalar species-block gradient
// `grad_b` and curvature `D_b`, and the species-top cross column `C_col`
// (length m). Returns the per-cell marginal log-lik. The Bernoulli cell's
// observed negative Hessian equals the Fisher info, so there is no
// observed/complete-data branch.
inline double jsdm_site_blocks(const JsdmCommData& d_, int s, int site,
                               const VectorXd& beta, double b_s,
                               double field_offset, const double* load,
                               const int* fidx, int nfield,
                               bool want_block,
                               VectorXd& grad_top, MatrixXd& Ttop,
                               double& grad_b, double& D_b, VectorXd& C_col) {
    const int p_occ = d_.p_occ;
    double eta = field_offset + b_s;
    for (int c = 0; c < p_occ; ++c) eta += d_.X(site, c) * beta(c);
    eta = clamp30(eta);

    const JsdmSiteCell cell = jsdm_site_cell(eta, d_.y[s][site]);

    // Score: eta -> beta cols (X_i), field cols (loadings), and the scalar b_s.
    for (int c = 0; c < p_occ; ++c) grad_top(c) += d_.X(site, c) * cell.g;
    for (int f = 0; f < nfield; ++f) grad_top(fidx[f]) += load[f] * cell.g;
    grad_b += cell.g;

    if (!want_block) return cell.log_lik;

    // Augmented top-design row z (length m): X_i on the beta cols, loadings on
    // the field cols. The curvature scatter is B * z z' (top) and B for b_s, with
    // B * z the species-top cross column.
    const double B = cell.B;
    for (int c = 0; c < p_occ; ++c) {
        const double xc = d_.X(site, c);
        for (int c2 = 0; c2 < p_occ; ++c2) Ttop(c, c2) += B * xc * d_.X(site, c2);
        for (int f = 0; f < nfield; ++f) {
            const double cross = B * xc * load[f];
            Ttop(c, fidx[f]) += cross; Ttop(fidx[f], c) += cross;
        }
        C_col(c) += B * xc;   // cross between beta_c and b_s
    }
    for (int f = 0; f < nfield; ++f) {
        for (int f2 = 0; f2 < nfield; ++f2)
            Ttop(fidx[f], fidx[f2]) += B * load[f] * load[f2];
        C_col(fidx[f]) += B * load[f];   // cross between field and b_s
    }
    D_b += B;                            // -d^2 / d b_s^2
    return cell.log_lik;
}

struct JsdmCommResult {
    VectorXd beta;
    VectorXd field;
    double sigma_re2 = 0.0;
    VectorXd blup;            // length n_species (b_s)
    MatrixXd vcov_beta;       // p_occ x p_occ
    double log_marginal = R_NegInf;
    double log_lik = R_NegInf;
    bool converged = false;
    int n_iter = 0;
};

inline double jsdm_field_log_prior(JsdmFieldKind kind, int n_spatial, double tau,
                                   double rho, double log_det_Q_rho,
                                   const IntegerVector& adj_row_ptr,
                                   const IntegerVector& adj_col_idx,
                                   const IntegerVector& n_neighbors,
                                   const VectorXd& field) {
    if (kind == JsdmFieldKind::BYM2) {
        const VectorXd v = field.head(n_spatial);
        const VectorXd w = field.tail(n_spatial);
        return nmix_bym2_log_prior(n_spatial, adj_row_ptr, adj_col_idx,
                                   n_neighbors, v, w);
    }
    if (kind == JsdmFieldKind::ICAR)
        return nmix_icar_log_prior(n_spatial, tau, adj_row_ptr, adj_col_idx,
                                   n_neighbors, field);
    return nmix_car_proper_log_prior(n_spatial, tau, rho, log_det_Q_rho,
                                     adj_row_ptr, adj_col_idx, n_neighbors, field);
}

// Add the field prior to the top-block gradient + Hessian (field_start = p_occ,
// no detection arm so p_p = 0). ICAR uses rho = 1; proper-CAR the grid rho; BYM2
// splits the (v, w) blocks.
inline void jsdm_add_field_prior(JsdmFieldKind kind, int p_occ, int n_spatial,
                                 double tau, double rho,
                                 const IntegerVector& adj_row_ptr,
                                 const IntegerVector& adj_col_idx,
                                 const IntegerVector& n_neighbors,
                                 const VectorXd& field, VectorXd& grad_top,
                                 MatrixXd& T) {
    if (kind == JsdmFieldKind::BYM2) {
        const VectorXd v = field.head(n_spatial);
        const VectorXd w = field.tail(n_spatial);
        nmix_add_bym2_prior_to_grad_and_H(p_occ, 0, n_spatial, adj_row_ptr,
                                          adj_col_idx, n_neighbors, v, w,
                                          grad_top, T);
        return;
    }
    const double rho_use = (kind == JsdmFieldKind::CAR_PROPER) ? rho : 1.0;
    nmix_add_car_to_spatial_block(p_occ, 0, n_spatial, tau, rho_use,
                                  adj_row_ptr, adj_col_idx, n_neighbors, field,
                                  grad_top, T);
}

// Sum-to-zero centering of the intrinsic field against the global intercept.
// ICAR centres f; BYM2 centres v only (w proper); CAR_proper is full rank.
inline void jsdm_center_field(JsdmFieldKind kind, int p_occ, int n_spatial,
                              VectorXd& field) {
    if (kind == JsdmFieldKind::CAR_PROPER || n_spatial <= 0) return;
    const int len = field.size();
    VectorXd holder(p_occ + len);
    holder.setZero();
    holder.segment(p_occ, len) = field;
    if (kind == JsdmFieldKind::BYM2)
        nmix_center_v_bym2(p_occ, 0, n_spatial, holder);
    else
        nmix_center_z(p_occ, 0, n_spatial, holder);
    field = holder.segment(p_occ, len);
}

// One outer-grid-point Laplace-EM fit over the JSDM Bernoulli marginal.
JsdmCommResult jsdm_spatial_em(
    JsdmFieldKind kind, const JsdmCommData& d_, int n_spatial,
    const IntegerVector& adj_row_ptr, const IntegerVector& adj_col_idx,
    const IntegerVector& n_neighbors,
    double tau, double rho, double log_det_Q_rho, double a, double b,
    const VectorXd& beta_init, double sigma_re2_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose) {

    const int p_occ = d_.p_occ, S = d_.n_species;
    const int nfield = (kind == JsdmFieldKind::BYM2) ? 2 : 1;
    const int field_len = nfield * n_spatial;
    const int m = p_occ + field_len;

    JsdmCommResult out;
    VectorXd beta = beta_init;
    VectorXd field = VectorXd::Zero(field_len);
    std::vector<double> bvec(S, 0.0);
    double sig2 = sigma_re2_init;

    auto site_geom = [&](int site, const VectorXd& field_) {
        const int u = d_.map_site_to_unit[site];
        return std::make_pair(jsdm_field_geom(kind, p_occ, n_spatial, u, a, b),
                              jsdm_field_offset(kind, field_, n_spatial, u, a, b));
    };
    auto data_loglik = [&](const VectorXd& beta_, const VectorXd& field_,
                           const std::vector<double>& b_) {
        double ll = 0.0;
        VectorXd gtop = VectorXd::Zero(m), ccol = VectorXd::Zero(m);
        MatrixXd Tt; double gb = 0.0, db = 0.0;
        for (int s = 0; s < S; ++s)
            for (int site = 0; site < d_.n_sites; ++site) {
                auto go = site_geom(site, field_);
                gtop.setZero(); gb = 0.0;
                double l = jsdm_site_blocks(d_, s, site, beta_, b_[s], go.second,
                                            go.first.load, go.first.idx,
                                            go.first.nfield, false,
                                            gtop, Tt, gb, db, ccol);
                if (!R_finite(l)) return R_NegInf;
                ll += l;
            }
        return ll;
    };
    auto objective = [&](const VectorXd& beta_, const VectorXd& field_,
                         const std::vector<double>& b_, double prec) {
        double obj = data_loglik(beta_, field_, b_);
        if (!R_finite(obj)) return R_NegInf;
        obj += jsdm_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                    adj_row_ptr, adj_col_idx, n_neighbors, field_);
        for (int s = 0; s < S; ++s) obj -= 0.5 * prec * b_[s] * b_[s];
        return obj;
    };
    // Assemble the top gradient/Hessian, the per-species scalar gradient/curvature
    // and the species-top cross column. The species prior precision `prec` =
    // 1/sig2 enters the scalar block.
    auto assemble = [&](const VectorXd& beta_, const VectorXd& field_,
                        const std::vector<double>& b_, double prec,
                        VectorXd& grad_top, MatrixXd& T,
                        std::vector<double>& grad_b, std::vector<double>& D,
                        std::vector<VectorXd>& C, double* log_lik_out) {
        grad_top.setZero(m); T.setZero(m, m);
        double ll = 0.0;
        VectorXd gtop_s(m), ccol(m); MatrixXd Tt(m, m);
        for (int s = 0; s < S; ++s) {
            grad_b[s] = 0.0; D[s] = 0.0; C[s].setZero(m);
            for (int site = 0; site < d_.n_sites; ++site) {
                auto go = site_geom(site, field_);
                ll += jsdm_site_blocks(d_, s, site, beta_, b_[s], go.second,
                                       go.first.load, go.first.idx,
                                       go.first.nfield, true,
                                       grad_top, T, grad_b[s], D[s], C[s]);
            }
            D[s] += prec;
            grad_b[s] -= prec * b_[s];
        }
        jsdm_add_field_prior(kind, p_occ, n_spatial, tau, rho,
                             adj_row_ptr, adj_col_idx, n_neighbors, field_,
                             grad_top, T);
        if (log_lik_out) *log_lik_out = ll;
    };

    VectorXd grad_top(m); MatrixXd T(m, m);
    std::vector<double> grad_b(S, 0.0), D(S, 0.0), Dinv(S, 0.0);
    std::vector<VectorXd> C(S, VectorXd::Zero(m));

    for (int em_it = 0; em_it < max_iter_em; ++em_it) {
        out.n_iter = em_it + 1;
        const double prec = (sig2 > 1e-12) ? 1.0 / sig2 : 1e12;

        for (int nit = 0; nit < inner_max; ++nit) {
            double ll_cur = 0.0;
            assemble(beta, field, bvec, prec, grad_top, T, grad_b, D, C, &ll_cur);
            MatrixXd M = T; VectorXd rhs = grad_top;
            for (int s = 0; s < S; ++s) {
                Dinv[s] = (D[s] > 0.0) ? 1.0 / D[s] : 0.0;
                M.noalias()   -= Dinv[s] * (C[s] * C[s].transpose());
                rhs.noalias() -= Dinv[s] * grad_b[s] * C[s];
            }
            nmix_add_diagonal_ridge(M);
            const VectorXd delta_top = nmix_safe_inverse(M) * rhs;
            std::vector<double> delta_b(S, 0.0);
            for (int s = 0; s < S; ++s)
                delta_b[s] = Dinv[s] * (grad_b[s] - C[s].dot(delta_top));
            const VectorXd dbeta  = delta_top.head(p_occ);
            const VectorXd dfield = delta_top.segment(p_occ, field_len);

            double step = 1.0; bool stepped = false;
            double obj_cur = ll_cur
                + jsdm_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                       adj_row_ptr, adj_col_idx, n_neighbors, field);
            for (int s = 0; s < S; ++s) obj_cur -= 0.5 * prec * bvec[s] * bvec[s];
            double max_step = 0.0;
            for (int h = 0; h < 12; ++h) {
                VectorXd beta_try = beta + step * dbeta;
                VectorXd field_try = field + step * dfield;
                std::vector<double> b_try(S);
                for (int s = 0; s < S; ++s) b_try[s] = bvec[s] + step * delta_b[s];
                const double obj_try = objective(beta_try, field_try, b_try, prec);
                if (R_finite(obj_try) && obj_try >= obj_cur - 1e-10) {
                    beta = beta_try; field = field_try; bvec = b_try;
                    jsdm_center_field(kind, p_occ, n_spatial, field);
                    double dmax = std::max(dbeta.cwiseAbs().maxCoeff(),
                                           dfield.cwiseAbs().maxCoeff());
                    for (int s = 0; s < S; ++s)
                        dmax = std::max(dmax, std::abs(delta_b[s]));
                    max_step = step * dmax; stepped = true; break;
                }
                step *= 0.5;
            }
            if (!stepped) break;
            if (max_step < inner_tol) break;
        }

        assemble(beta, field, bvec, prec, grad_top, T, grad_b, D, C, nullptr);
        for (int s = 0; s < S; ++s) Dinv[s] = (D[s] > 0.0) ? 1.0 / D[s] : 0.0;

        double sig2_new = 0.0;
        for (int s = 0; s < S; ++s) sig2_new += bvec[s] * bvec[s] + Dinv[s];
        sig2_new /= (double) S;
        const double dSig = std::abs(sig2_new - sig2);
        sig2 = sig2_new;
        if (verbose) Rcpp::Rcout << "  jsdm-em " << out.n_iter << " dsigma2=" << dSig << "\n";
        if (dSig < tol_em) { out.converged = true; break; }
    }

    const double prec = (sig2 > 1e-12) ? 1.0 / sig2 : 1e12;
    const double logdetP = std::log(prec);

    auto final_assemble = [&](double& loglik_marg, MatrixXd& vcov_beta) -> bool {
        VectorXd gtop(m); MatrixXd Tt(m, m);
        std::vector<double> gb(S, 0.0), Dd(S, 0.0);
        std::vector<VectorXd> Cc(S, VectorXd::Zero(m));
        double ll = 0.0;
        assemble(beta, field, bvec, prec, gtop, Tt, gb, Dd, Cc, &ll);
        double sum_logdet_D = 0.0; MatrixXd M = Tt;
        for (int s = 0; s < S; ++s) {
            if (!(Dd[s] > 0.0)) return false;
            sum_logdet_D += std::log(Dd[s]);
            M.noalias() -= (1.0 / Dd[s]) * (Cc[s] * Cc[s].transpose());
        }
        MatrixXd M_det = M; nmix_add_diagonal_ridge(M_det);
        const double ldM = nmix_logdet_spd(M_det);
        if (!R_finite(ldM)) return false;
        double lp = jsdm_field_log_prior(kind, n_spatial, tau, rho, log_det_Q_rho,
                                         adj_row_ptr, adj_col_idx, n_neighbors, field);
        double bquad = 0.0;
        for (int s = 0; s < S; ++s) bquad += prec * bvec[s] * bvec[s];
        loglik_marg = ll + lp + 0.5 * (double) S * logdetP - 0.5 * bquad
                    - 0.5 * sum_logdet_D - 0.5 * ldM;
        const bool constrain = (kind != JsdmFieldKind::CAR_PROPER);
        MatrixXd cov_top = nmix_constrained_top_cov(M, m, p_occ, p_occ, field_len,
                                                    constrain);
        if (!cov_top.allFinite()) return false;
        vcov_beta = cov_top;
        return true;
    };

    double loglik_marg = R_NegInf; MatrixXd vcov_beta;
    bool ok = final_assemble(loglik_marg, vcov_beta);

    out.beta = beta; out.field = field; out.sigma_re2 = sig2;
    out.blup = VectorXd(S);
    for (int s = 0; s < S; ++s) out.blup(s) = bvec[s];
    out.vcov_beta = ok ? vcov_beta : MatrixXd::Constant(p_occ, p_occ, R_NaN);
    out.log_marginal = ok ? loglik_marg : R_NegInf;
    out.log_lik = data_loglik(beta, field, bvec);
    return out;
}

JsdmCommData build_jsdm_data(const NumericMatrix& X_R, const IntegerMatrix& y_R,
                             const IntegerVector& map_R, int n_spatial) {
    JsdmCommData d_(X_R);
    d_.n_sites   = X_R.nrow();
    d_.p_occ     = X_R.ncol();
    d_.n_species = y_R.ncol();
    if (y_R.nrow() != d_.n_sites)
        Rcpp::stop("y must be n_sites x n_species.");
    d_.y.assign(d_.n_species, std::vector<int>(d_.n_sites, 0));
    for (int s = 0; s < d_.n_species; ++s)
        for (int i = 0; i < d_.n_sites; ++i)
            d_.y[s][i] = y_R(i, s);
    d_.map_site_to_unit.assign(d_.n_sites, 0);
    for (int i = 0; i < d_.n_sites; ++i) {
        const int u = map_R[i] - 1;
        if (u < 0 || u >= n_spatial)
            Rcpp::stop("map_site_to_unit out of range [1, n_spatial].");
        d_.map_site_to_unit[i] = u;
    }
    return d_;
}

Rcpp::List pack_jsdm_grid(int p_occ, int field_len, int n_grid,
                          const std::vector<JsdmCommResult>& results,
                          const NumericMatrix& theta_grid_out, int n_spatial) {
    NumericVector log_marginals(n_grid), log_liks(n_grid), sigma_re2(n_grid);
    IntegerVector n_iters(n_grid);
    LogicalVector convergeds(n_grid);
    NumericMatrix modes(n_grid, p_occ + field_len);
    Rcpp::List vcov_beta(n_grid), blup(n_grid);
    for (int k = 0; k < n_grid; ++k) {
        const JsdmCommResult& rr = results[k];
        log_marginals[k] = rr.log_marginal;
        log_liks[k]      = rr.log_lik;
        sigma_re2[k]     = rr.sigma_re2;
        n_iters[k]       = rr.n_iter;
        convergeds[k]    = rr.converged;
        for (int j = 0; j < p_occ; ++j)     modes(k, j) = rr.beta(j);
        for (int j = 0; j < field_len; ++j) modes(k, p_occ + j) = rr.field(j);
        vcov_beta[k] = Rcpp::wrap(rr.vcov_beta);
        blup[k]      = Rcpp::wrap(rr.blup);
    }
    return Rcpp::List::create(
        Rcpp::Named("theta_grid")   = theta_grid_out,
        Rcpp::Named("log_marginal") = log_marginals,
        Rcpp::Named("log_lik")      = log_liks,
        Rcpp::Named("sigma_re2")    = sigma_re2,
        Rcpp::Named("modes")        = modes,
        Rcpp::Named("vcov_beta")    = vcov_beta,
        Rcpp::Named("blup")         = blup,
        Rcpp::Named("n_iter")       = n_iters,
        Rcpp::Named("converged")    = convergeds,
        Rcpp::Named("p_occ")        = p_occ,
        Rcpp::Named("n_spatial")    = n_spatial,
        Rcpp::Named("n_grid")       = n_grid);
}

struct JsdmGridPoint { double tau, rho, log_det_Q_rho, a, b; };

Rcpp::List run_jsdm_spatial_grid(
    JsdmFieldKind kind, int field_len,
    const NumericMatrix& X_R, const IntegerMatrix& y_R,
    const IntegerVector& map_R, int n_spatial,
    const IntegerVector& adj_row_ptr, const IntegerVector& adj_col_idx,
    const IntegerVector& n_neighbors,
    const std::vector<JsdmGridPoint>& plan, const NumericMatrix& theta_grid_out,
    const NumericVector& beta_init, double sigma_re2_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose) {

    JsdmCommData d_ = build_jsdm_data(X_R, y_R, map_R, n_spatial);
    const int p_occ = d_.p_occ;
    VectorXd beta0(p_occ);
    for (int j = 0; j < p_occ; ++j) beta0(j) = beta_init[j];

    const int n_grid = (int) plan.size();
    std::vector<JsdmCommResult> results(n_grid);
    for (int k = 0; k < n_grid; ++k) {
        const JsdmGridPoint& g = plan[k];
        results[k] = jsdm_spatial_em(
            kind, d_, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
            g.tau, g.rho, g.log_det_Q_rho, g.a, g.b,
            beta0, sigma_re2_init, max_iter_em, tol_em, inner_max, inner_tol,
            verbose);
    }
    return pack_jsdm_grid(p_occ, field_len, n_grid, results, theta_grid_out,
                          n_spatial);
}

}  // namespace
}  // namespace tulpaObs


// ---------------------------------------------------------------------------
// Cross-check entry: per-(species, site) Bernoulli cell log-lik, score, and
// negative Hessian vs an R oracle (FD-validated).
// [[Rcpp::export]]
Rcpp::List cpp_jsdm_site_cell(double eta, int y) {
    const tulpaObs::JsdmSiteCell c = tulpaObs::jsdm_site_cell(eta, y);
    return Rcpp::List::create(
        Rcpp::Named("log_lik") = c.log_lik,
        Rcpp::Named("grad") = c.g,
        Rcpp::Named("neg_hess") = c.B);
}

// ICAR
// [[Rcpp::export]]
Rcpp::List cpp_jsdm_spatial_icar(
    Rcpp::NumericMatrix X, Rcpp::IntegerMatrix y,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector beta_init, double sigma_re2_init,
    int max_iter_em, bool verbose) {
    std::vector<tulpaObs::JsdmGridPoint> plan;
    Rcpp::NumericMatrix tg(tau_grid.size(), 1);
    for (int i = 0; i < tau_grid.size(); ++i) {
        plan.push_back({tau_grid[i], 1.0, 0.0, 0.0, 0.0});
        tg(i, 0) = tau_grid[i];
    }
    colnames(tg) = Rcpp::CharacterVector::create("tau");
    return tulpaObs::run_jsdm_spatial_grid(
        tulpaObs::JsdmFieldKind::ICAR, n_spatial, X, y,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, beta_init, sigma_re2_init, max_iter_em, 1e-4, 50, 1e-7, verbose);
}

// Proper CAR
// [[Rcpp::export]]
Rcpp::List cpp_jsdm_spatial_car_proper(
    Rcpp::NumericMatrix X, Rcpp::IntegerMatrix y,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector rho_grid, Rcpp::NumericVector log_det_Q_rho,
    Rcpp::NumericVector beta_init, double sigma_re2_init,
    int max_iter_em, bool verbose) {
    std::vector<tulpaObs::JsdmGridPoint> plan;
    std::vector<double> tg_tau, tg_rho;
    for (int ri = 0; ri < rho_grid.size(); ++ri)
        for (int ti = 0; ti < tau_grid.size(); ++ti) {
            plan.push_back({tau_grid[ti], rho_grid[ri], log_det_Q_rho[ri], 0.0, 0.0});
            tg_tau.push_back(tau_grid[ti]); tg_rho.push_back(rho_grid[ri]);
        }
    Rcpp::NumericMatrix tg(plan.size(), 2);
    for (std::size_t i = 0; i < plan.size(); ++i) { tg(i, 0) = tg_tau[i]; tg(i, 1) = tg_rho[i]; }
    colnames(tg) = Rcpp::CharacterVector::create("tau", "rho");
    return tulpaObs::run_jsdm_spatial_grid(
        tulpaObs::JsdmFieldKind::CAR_PROPER, n_spatial, X, y,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, beta_init, sigma_re2_init, max_iter_em, 1e-4, 50, 1e-7, verbose);
}

// BYM2
// [[Rcpp::export]]
Rcpp::List cpp_jsdm_spatial_bym2(
    Rcpp::NumericMatrix X, Rcpp::IntegerMatrix y,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector sigma_grid,
    Rcpp::NumericVector rho_grid, double scale_factor,
    Rcpp::NumericVector beta_init, double sigma_re2_init,
    int max_iter_em, bool verbose) {
    std::vector<tulpaObs::JsdmGridPoint> plan;
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
    return tulpaObs::run_jsdm_spatial_grid(
        tulpaObs::JsdmFieldKind::BYM2, 2 * n_spatial, X, y,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, beta_init, sigma_re2_init, max_iter_em, 1e-4, 50, 1e-7, verbose);
}
