// cell_coupling_occu_multiscale_cover.h
// Stateful `CellCouplingSpec` implementing the per-cell log-density of the
// three-level occupancy + cover hurdle: a cell-level occupancy gate (psi), a
// plot-level availability gate (theta), a visit-level detection (p) and the
// cover hurdle (pos) on detected visits.
//
//   z_c        ~ Bernoulli(psi_c)                  # cell occupancy
//   a_cj | z=1 ~ Bernoulli(theta_cj)               # plot availability
//   y_cjv|a=1  ~ Bernoulli(p_cjv)                   # detection
//   cover|y=1  ~ f_pos(.; eta_pos, phi)             # cover hurdle
//
// Arm layout (kk indexes arm_ids() = {0, 1, 2, 3}):
//   kk = 0 psi   : 1 row per cell,             eta = logit psi_c
//   kk = 1 theta : M_c rows per cell (plots),  eta = logit theta_cj
//   kk = 2 p     : sum_j J_cj rows per cell,   eta = logit p_cjv (plot-major)
//   kk = 3 pos   : same rows as arm 2,         positive observation at detected
//                  visits (0 elsewhere); phi(3) is the pos dispersion.
//
// The spec is STATEFUL: it carries `cell_plot_sizes_[c]` (visits per plot for
// cell c) captured at construction so it can partition arm 2/3's flat per-cell
// visit rows into plots aligned with arm 1's rows. The visit rows of a cell
// are laid plot-major (plot 0's visits, then plot 1's, ...), matching the row
// order `build_cell_rows_from_arms` produces from the input ordering.
//
// Per-cell density:
//   Q_j    = theta_j G_j + (1 - theta_j) H_j         (marginalize a_cj | z=1)
//   L_cell = psi prod_j Q_j + (1 - psi) 1[all y = 0]  (marginalize z_c)
// splits on whether ANY visit in the cell detected:
//   * any_det : log L = log psi + sum_j log Q_j  (plots independent given z=1;
//               a detected plot factorizes, a non-detected plot is a within-
//               plot occupancy mixture -> shared nodet_mixture_block(w=theta_j))
//   * no_det  : L = psi M + (1 - psi), M = prod_j m_j, the nested mixture whose
//               plots couple through psi.
//
// All closed-form first + second + cross derivatives below are FD-checked in
// tulpaObs/tests/testthat/test-occu-multiscale-cover-coupling.R against
// numerical derivatives of the cell density. Derivation:
// dev_notes/occu_multiscale_cover_derivation.md.

#ifndef TULPAOBS_CELL_COUPLING_OCCU_MULTISCALE_COVER_H
#define TULPAOBS_CELL_COUPLING_OCCU_MULTISCALE_COVER_H

#include <tulpa/cell_coupling.h>
#include "occu_coupling_shared.h"  // sigmoid_ / PosPolicy / nodet_mixture_block
#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace tulpaObs {

template <class PosPolicy>
class OccuMultiscaleCoverCoupling final : public tulpa::CellCouplingSpec {
public:
    // `cell_plot_sizes[c]` = visit counts of cell c's plots, in plot order.
    // Captured by value so the registry's shared_ptr keeps it alive for the
    // fit (the CellCouplingSpec lifetime contract).
    explicit OccuMultiscaleCoverCoupling(std::vector<std::vector<int>> cell_plot_sizes)
        : cell_plot_sizes_(std::move(cell_plot_sizes)) {}

    std::vector<int> arm_ids() const override { return {0, 1, 2, 3}; }

    // The blocks this spec actually writes: (psi, theta), (psi, p),
    // (theta, theta), (theta, p) and the dense (p, p). The cover arm (3)
    // enters only the factorized detected-plot terms, which carry no cross, so
    // every (., pos) block stays zero -- stated at the write site below and
    // enforced here. The engine's default allocates every kk <= ll pair, which
    // for four arms means four slabs this spec never touches, two of them
    // Jc x Jc. There is no rank-1 self-cross path here: (p, p) is written
    // densely either way, so the list does not depend on the flag.
    std::vector<std::pair<int, int>> dense_cross_pairs(
            int /*n_coupled*/, bool /*rank1_self_supported*/) const override {
        return {{0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}};
    }

    double evaluate_cell(int                        cell_idx,
                         const tulpa::CellEtas&     etas,
                         const tulpa::CellResponse& y_cell,
                         tulpa::CellDerivs&         out) const override {
        const int M_c     = etas.n_rows_in_arm(1);   // plots in this cell
        const int Jc      = etas.n_rows_in_arm(2);   // total visits in this cell
        const double psi  = sigmoid_(etas.eta(0, 0));
        const double phi  = y_cell.phi(3);
        const bool want_hess = !out.grad_only;
        const bool expected  = (out.curvature == tulpa::CurvatureMode::Expected);

        const std::vector<int>& sizes = cell_plot_sizes_[cell_idx];

        // Per-plot visit offsets into the flat arm-2/3 row list, plus the
        // per-plot any-detection flag and the cell-wide flag.
        std::vector<int>  off(M_c);
        std::vector<char> plot_det(M_c, 0);
        bool any_det_cell = false;
        {
            int running = 0;
            for (int j = 0; j < M_c; j++) {
                off[j] = running;
                const int sz = sizes[j];
                for (int v = 0; v < sz; v++) {
                    if (y_cell.y(2, running + v) > 0.5) { plot_det[j] = 1; any_det_cell = true; }
                }
                running += sz;
            }
        }

        if (any_det_cell) {
            return eval_any_det_(M_c, Jc, sizes, off, plot_det, psi, phi,
                                 want_hess, expected, etas, y_cell, out);
        }
        return eval_no_det_(M_c, Jc, sizes, off, psi,
                            want_hess, expected, etas, out);
    }

    std::string name() const override {
        return std::string(PosPolicy::multiscale_spec_name());
    }

    bool thread_safe() const override { return true; }

private:
    std::vector<std::vector<int>> cell_plot_sizes_;

    // ---- Branch A: any detection in the cell (z = 1 forced) --------------
    // log L = log psi + sum_j log Q_j. log psi separates; each plot is
    // independent. A detected plot factorizes (standard GLM terms, no cross);
    // a non-detected plot is a within-plot occupancy mixture handled by the
    // shared nodet_mixture_block with w = theta_j.
    double eval_any_det_(int M_c, int Jc,
                         const std::vector<int>& sizes,
                         const std::vector<int>& off,
                         const std::vector<char>& plot_det,
                         double psi, double phi,
                         bool want_hess, bool expected,
                         const tulpa::CellEtas&     etas,
                         const tulpa::CellResponse& y_cell,
                         tulpa::CellDerivs&         out) const {
        double cell_ll = log_safe(psi);
        out.arm_grad[0][0] = 1.0 - psi;
        if (want_hess) out.arm_neg_hess_diag[0][0] = psi * (1.0 - psi);

        for (int j = 0; j < M_c; j++) {
            const int sz = sizes[j];
            const double eta_theta = etas.eta(1, j);
            const double theta = sigmoid_(eta_theta);

            if (plot_det[j]) {
                // Detected plot: log theta_j + sum_v [detected: log p + log f_pos;
                // undetected: log(1-p)]. Fully factorized -> no cross-Hessians.
                cell_ll += log_safe(theta);
                out.arm_grad[1][j] = 1.0 - theta;
                if (want_hess) out.arm_neg_hess_diag[1][j] = theta * (1.0 - theta);

                for (int v = 0; v < sz; v++) {
                    const int row = off[j] + v;
                    const double p_v = sigmoid_(etas.eta(2, row));
                    if (y_cell.y(2, row) > 0.5) {
                        cell_ll += log_safe(p_v);
                        out.arm_grad[2][row] = 1.0 - p_v;
                        if (want_hess) out.arm_neg_hess_diag[2][row] = p_v * (1.0 - p_v);

                        const double y_pos   = y_cell.y(3, row);
                        const double eta_pos = etas.eta(3, row);
                        cell_ll += PosPolicy::log_density(y_pos, eta_pos, phi);
                        double g_pos = 0.0, h_pos = 0.0;
                        PosPolicy::grad_hess_eta(y_pos, eta_pos, phi,
                                                 want_hess, g_pos, h_pos);
                        out.arm_grad[3][row] = g_pos;
                        if (want_hess) out.arm_neg_hess_diag[3][row] = h_pos;
                    } else {
                        cell_ll += log_safe(1.0 - p_v);
                        out.arm_grad[2][row] = -p_v;
                        if (want_hess) out.arm_neg_hess_diag[2][row] = p_v * (1.0 - p_v);
                    }
                }
            } else {
                // Non-detected plot inside an occupied cell: within-plot
                // occupancy mixture m_j = theta_j P0_j + (1 - theta_j). Reuse
                // the shared nodet block with w = theta_j over this plot's
                // visits; place its compact outputs into the cell buffers.
                std::vector<double> eta_p_buf(sz);
                for (int v = 0; v < sz; v++) eta_p_buf[v] = etas.eta(2, off[j] + v);

                // arm_grad[2] / arm_neg_hess_diag[2] are contiguous over the
                // plot's visit rows, so write through directly at offset.
                double* g_p  = out.arm_grad[2] + off[j];
                double* nh_p = out.arm_neg_hess_diag[2] + off[j];
                // (theta_j, p_jv): arm_cross_hess[1][2] row j, cols off[j]..
                // are contiguous -> write through. (p_jv, p_jw): not contiguous
                // (row stride Jc, block stride sz) -> compact local + scatter.
                double* cross_t_p = nullptr;
                std::vector<double> cpp_loc;
                if (!expected && want_hess && out.arm_cross_hess) {
                    if (out.arm_cross_hess[1] && out.arm_cross_hess[1][2])
                        cross_t_p = out.arm_cross_hess[1][2] + (std::size_t)j * Jc + off[j];
                    if (out.arm_cross_hess[2] && out.arm_cross_hess[2][2])
                        cpp_loc.assign((std::size_t)sz * sz, 0.0);
                }

                double g_th = 0.0, nh_th = 0.0;
                cell_ll += nodet_mixture_block(
                    theta, eta_p_buf.data(), sz, want_hess, expected,
                    g_th, nh_th, g_p, nh_p,
                    cross_t_p, cpp_loc.empty() ? nullptr : cpp_loc.data());
                out.arm_grad[1][j] = g_th;
                if (want_hess) out.arm_neg_hess_diag[1][j] = nh_th;

                if (!cpp_loc.empty()) {
                    double* cpp = out.arm_cross_hess[2][2];
                    for (int v = 0; v < sz; v++) {
                        for (int w = v + 1; w < sz; w++) {
                            const double val = cpp_loc[(std::size_t)v * sz + w];
                            const int rv = off[j] + v, rw = off[j] + w;
                            cpp[(std::size_t)rv * Jc + rw] = val;
                            cpp[(std::size_t)rw * Jc + rv] = val;
                        }
                    }
                }
            }
        }
        // (.,pos) cross-Hessians: zero (pos enters only the factorized det
        // terms, which carry no cross). Buffers come zeroed.
        return cell_ll;
    }

    // ---- Branch B: no detection anywhere in the cell ---------------------
    // L = psi M + (1 - psi), M = prod_j m_j, m_j = theta_j P0_j + (1 - theta_j),
    // P0_j = prod_v (1 - p_jv). The nested mixture; plots couple through psi.
    double eval_no_det_(int M_c, int Jc,
                        const std::vector<int>& sizes,
                        const std::vector<int>& off,
                        double psi,
                        bool want_hess, bool expected,
                        const tulpa::CellEtas& etas,
                        tulpa::CellDerivs&     out) const {
        // The cell mixture, its intermediates and its scores are the same
        // derivation the NUTS target evaluates over its own layout, so it is
        // computed once in occu_coupling_shared.h and read here through the
        // CellEtas accessors. Everything below is this layout's own business:
        // writing the scores into the arm buffers, and the curvature.
        const MscaleNoDetCell cell = mscale_nodet_cell(
            psi, M_c,
            [&](int j)   { return etas.eta(1, j); },
            off.data(), sizes.data(),
            [&](int row) { return etas.eta(2, row); });

        const double  M      = cell.M;
        const double  logM   = cell.logM;
        const double  L      = cell.L;
        const double  inv_L  = cell.inv_L;
        const double  psi_d  = cell.psi_d;
        const double  s_psi  = cell.s_psi;
        const std::vector<double>& theta   = cell.theta;
        const std::vector<double>& theta_d = cell.theta_d;
        const std::vector<double>& P0      = cell.P0;
        const std::vector<double>& m       = cell.m;
        const std::vector<double>& logm    = cell.logm;
        const std::vector<double>& A       = cell.A;
        const std::vector<double>& Mmj     = cell.Mmj;
        const std::vector<double>& p_val   = cell.p_val;
        const std::vector<double>& B       = cell.B;
        const std::vector<double>& s_theta = cell.s_theta;
        const std::vector<double>& s_p     = cell.s_p;
        const std::vector<int>&    plot_of = cell.plot_of;

        out.arm_grad[0][0] = s_psi;
        for (int j = 0; j < M_c; j++) out.arm_grad[1][j] = s_theta[j];
        for (int row = 0; row < Jc; row++) out.arm_grad[2][row] = s_p[row];

        if (!want_hess) return log_safe(L);

        if (expected) {
            // Complete-data Fisher: block-diagonal, responsibility-weighted.
            const double gamma_c = (L > 0.0) ? (psi * M * inv_L) : 0.0;
            out.arm_neg_hess_diag[0][0] = psi_d;
            for (int j = 0; j < M_c; j++) {
                out.arm_neg_hess_diag[1][j] = gamma_c * theta_d[j];
                const double gamma_j = (m[j] > 0.0) ? (theta[j] * P0[j] / m[j]) : 0.0;
                const int sz = sizes[j];
                for (int v = 0; v < sz; v++) {
                    const int row = off[j] + v;
                    out.arm_neg_hess_diag[2][row] =
                        gamma_c * gamma_j * p_val[row] * (1.0 - p_val[row]);
                }
            }
            return log_safe(L);
        }

        // Observed (true mixture) Hessian. Diagonal: -d2 logL = -(d2L)/L + s^2.
        out.arm_neg_hess_diag[0][0] =
            -(psi_d * (1.0 - 2.0 * psi) * (M - 1.0)) * inv_L + s_psi * s_psi;
        for (int j = 0; j < M_c; j++) {
            const double d2L_th =
                psi * Mmj[j] * (-(1.0 - P0[j]) * theta_d[j] * (1.0 - 2.0 * theta[j]));
            out.arm_neg_hess_diag[1][j] = -d2L_th * inv_L + s_theta[j] * s_theta[j];
            const int sz = sizes[j];
            for (int v = 0; v < sz; v++) {
                const int row = off[j] + v;
                const double d2L_pp =
                    psi * Mmj[j] * (-theta[j] * P0[j] * p_val[row] * (1.0 - 2.0 * p_val[row]));
                out.arm_neg_hess_diag[2][row] = -d2L_pp * inv_L + s_p[row] * s_p[row];
            }
        }

        // Cross-Hessians: -d2 logL / d eta_x d eta_y = -(d2L)/L + s_x s_y.
        double* const* const* CH = out.arm_cross_hess;
        if (!CH) return log_safe(L);
        double* ch_psi_theta = (CH[0] && CH[0][1]) ? CH[0][1] : nullptr;  // 1 x M_c
        double* ch_psi_p     = (CH[0] && CH[0][2]) ? CH[0][2] : nullptr;  // 1 x Jc
        double* ch_theta_theta = (CH[1] && CH[1][1]) ? CH[1][1] : nullptr; // M_c x M_c
        double* ch_theta_p   = (CH[1] && CH[1][2]) ? CH[1][2] : nullptr;  // M_c x Jc
        double* ch_p_p       = (CH[2] && CH[2][2]) ? CH[2][2] : nullptr;  // Jc x Jc

        // (psi, theta_j) and (psi, p_jv).
        for (int j = 0; j < M_c; j++) {
            if (ch_psi_theta) {
                const double d2L = psi_d * Mmj[j] * A[j];
                ch_psi_theta[j] = -d2L * inv_L + s_psi * s_theta[j];
            }
        }
        if (ch_psi_p) {
            for (int j = 0; j < M_c; j++) {
                const int sz = sizes[j];
                for (int v = 0; v < sz; v++) {
                    const int row = off[j] + v;
                    const double d2L = psi_d * Mmj[j] * B[row];
                    ch_psi_p[row] = -d2L * inv_L + s_psi * s_p[row];
                }
            }
        }

        // (theta_j, theta_k), j < k (cross-plot only; M_{-j,-k}).
        if (ch_theta_theta) {
            for (int j = 0; j < M_c; j++) {
                for (int k = j + 1; k < M_c; k++) {
                    const double Mmjk = std::exp(logM - logm[j] - logm[k]);
                    const double d2L  = psi * Mmjk * A[j] * A[k];
                    const double val  = -d2L * inv_L + s_theta[j] * s_theta[k];
                    ch_theta_theta[(std::size_t)j * M_c + k] = val;
                    ch_theta_theta[(std::size_t)k * M_c + j] = val;
                }
            }
        }

        // (theta_j, p_kw) for all j, all visits w of plot k (same or cross plot).
        if (ch_theta_p) {
            for (int j = 0; j < M_c; j++) {
                for (int k = 0; k < M_c; k++) {
                    const int sz = sizes[k];
                    const double Mmjk = (j == k) ? 0.0
                                       : std::exp(logM - logm[j] - logm[k]);
                    for (int w = 0; w < sz; w++) {
                        const int row = off[k] + w;
                        double d2L;
                        if (j == k) {
                            // psi M_{-j} d2 m_j / d theta_j d p_jw
                            d2L = psi * Mmj[j] * (-theta_d[j] * P0[j] * p_val[row]);
                        } else {
                            d2L = psi * Mmjk * A[j] * B[row];
                        }
                        ch_theta_p[(std::size_t)j * Jc + row] =
                            -d2L * inv_L + s_theta[j] * s_p[row];
                    }
                }
            }
        }

        // (p_jv, p_kw), upper triangle in flat row index (same plot v!=w or
        // cross plot). Write both triangles for parity.
        if (ch_p_p) {
            for (int rv = 0; rv < Jc; rv++) {
                for (int rw = rv + 1; rw < Jc; rw++) {
                    const int j = plot_of[rv], k = plot_of[rw];
                    double d2L;
                    if (j == k) {
                        // psi M_{-j} d2 m_j / d p_jv d p_jw = psi M_{-j} theta_j P0_j p_v p_w
                        d2L = psi * Mmj[j] * (theta[j] * P0[j] * p_val[rv] * p_val[rw]);
                    } else {
                        const double Mmjk = std::exp(logM - logm[j] - logm[k]);
                        d2L = psi * Mmjk * B[rv] * B[rw];
                    }
                    const double val = -d2L * inv_L + s_p[rv] * s_p[rw];
                    ch_p_p[(std::size_t)rv * Jc + rw] = val;
                    ch_p_p[(std::size_t)rw * Jc + rv] = val;
                }
            }
        }
        return log_safe(L);
    }
};

typedef OccuMultiscaleCoverCoupling<LognormalPositive> OccuMultiscaleCoverLognormalCoupling;
typedef OccuMultiscaleCoverCoupling<BetaPositive>      OccuMultiscaleCoverBetaCoupling;
typedef OccuMultiscaleCoverCoupling<GaussianPositive>  OccuMultiscaleCoverGaussianCoupling;

} // namespace tulpaObs

#endif // TULPAOBS_CELL_COUPLING_OCCU_MULTISCALE_COVER_H
