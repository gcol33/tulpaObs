// cell_coupling_occu_cover.h
// Stateless `CellCouplingSpec` implementing the per-cell log-density of
// the occu-cover hurdle (psi/p/pos arms) for the joint nested-Laplace
// path in tulpa (gcol33/tulpa#32 Layer B.2 consumer).
//
// Arm layout (kk indexes the spec's arm_ids() return):
//   kk = 0 -> psi arm: 1 row per cell, carries no y data, eta = logit psi_c
//   kk = 1 -> p   arm: J_c rows per cell, y(1, v) in {0, 1} detection,
//                       eta = logit p_cv
//   kk = 2 -> pos arm: J_c rows per cell, y(2, v) = positive observation at
//                       detected visits (0 elsewhere). Interpretation +
//                       sigma_pos / phi role are policy-defined; the
//                       lognormal policy reads y_cell.phi(2) as the SD on
//                       the log scale, the beta policy reads it as the
//                       precision.
//
// Per-cell density splits on `any_det = any(y_det_cv == 1)`:
//
//   det case (any_det):
//     log p_cell = log psi + sum_v [
//                    (y_det = 1) (log p_v + log f_pos(y_pos_v; eta_pos_v, phi))
//                  + (y_det = 0)  log (1 - p_v) ]
//
//   nodet case (all y_det = 0):
//     log p_cell = log L,   L = psi P0 + (1 - psi),   P0 = prod_v (1 - p_v)
//
// Only the det branch's pos contribution changes between positive families
// — the psi / p derivatives and the nodet branch are family-independent.
// The templated `OccuCoverCoupling<PosPolicy, Aggregated>` factors that out;
// each concrete spec (`occu_cover_lognormal`, `occu_cover_beta`) is a typedef
// over a positive-arm policy supplying log_density / grad_eta / neg_hess_eta.
//
// `Aggregated` (tulpaObs#33) selects the cover arm's granularity:
//   * Aggregated = false (per-visit): the pos arm carries one row per visit
//     aligned with the detection arm; the det branch adds log f_pos at every
//     detected visit (the J_det-factor cover likelihood).
//   * Aggregated = true  (cell-aggregated): the pos arm carries ONE row per
//     occupancy unit holding the mean / median cover over that unit's detected
//     visits; the det branch adds a single log f_pos(ybar; eta_pos_cell) when
//     the cell has any detection. The cover arm then contributes one
//     observation per cell rather than one per detected visit, so a per-cell
//     cover signal no longer outweighs the occupancy arm on a shared field.
// Both share the psi / p derivatives and the nodet branch; only the pos
// evaluation count differs, and the pos term factorises from psi / p either
// way (cross-Hessians stay zero in the det branch).
//
// All closed-form first + second + cross derivatives below are FD-checked
// in tulpaObs/tests/testthat/test-occu-cover-coupling.R against numerical
// derivatives of the cell density above.

#ifndef TULPAOBS_CELL_COUPLING_OCCU_COVER_H
#define TULPAOBS_CELL_COUPLING_OCCU_COVER_H

#include <tulpa/cell_coupling.h>
#include "occu_coupling_shared.h"  // sigmoid_ / log_safe_ / PosPolicy / nodet_mixture_block
#include <cmath>
#include <string>
#include <vector>

namespace tulpaObs {

template <class PosPolicy, bool Aggregated = false>
class OccuCoverCoupling final : public tulpa::CellCouplingSpec {
public:
    std::vector<int> arm_ids() const override { return {0, 1, 2}; }

    double evaluate_cell(int                       /*cell_idx*/,
                         const tulpa::CellEtas&     etas,
                         const tulpa::CellResponse& y_cell,
                         tulpa::CellDerivs&         out) const override {
        const int Jc = etas.n_rows_in_arm(1);
        const double eta_psi = etas.eta(0, 0);
        const double psi     = sigmoid_(eta_psi);
        const double phi_pos = y_cell.phi(2);

        // Grad-only request (cached-factor reuse step): the kernel discards the
        // Hessian this call would produce, so skip every negative-Hessian and
        // cross-Hessian term -- including the positive arm's trigamma -- and
        // leave the pre-zeroed buffers untouched. The gradient and the returned
        // log-density are still exact.
        const bool want_hess = !out.grad_only;

        bool any_det = false;
        for (int v = 0; v < Jc; v++) {
            if (y_cell.y(1, v) > 0.5) { any_det = true; break; }
        }

        if (any_det) {
            double cell_ll = log_safe_(psi);

            out.arm_grad[0][0]          = 1.0 - psi;
            if (want_hess) out.arm_neg_hess_diag[0][0] = psi * (1.0 - psi);

            for (int v = 0; v < Jc; v++) {
                const double eta_p = etas.eta(1, v);
                const double p_v   = sigmoid_(eta_p);
                const double y_det = y_cell.y(1, v);

                if (y_det > 0.5) {
                    cell_ll += log_safe_(p_v);
                    out.arm_grad[1][v]          = 1.0 - p_v;
                    if (want_hess) out.arm_neg_hess_diag[1][v] = p_v * (1.0 - p_v);

                    if (!Aggregated) {
                        // Per-visit cover: one log f_pos per detected visit, the
                        // pos arm row aligned with this detection visit.
                        const double y_pos   = y_cell.y(2, v);
                        const double eta_pos = etas.eta(2, v);
                        cell_ll += PosPolicy::log_density(y_pos, eta_pos, phi_pos);
                        double g_pos = 0.0, h_pos = 0.0;
                        PosPolicy::grad_hess_eta(y_pos, eta_pos, phi_pos,
                                                 want_hess, g_pos, h_pos);
                        out.arm_grad[2][v] = g_pos;
                        if (want_hess) out.arm_neg_hess_diag[2][v] = h_pos;
                    }
                } else {
                    cell_ll += log_safe_(1.0 - p_v);
                    out.arm_grad[1][v]          = -p_v;
                    if (want_hess) out.arm_neg_hess_diag[1][v] = p_v * (1.0 - p_v);
                }
            }
            if (Aggregated) {
                // Cell-aggregated cover: a single log f_pos at the occupancy
                // unit's one pos row (the mean / median cover over its detected
                // visits), evaluated once because any_det holds here. The pos
                // arm carries exactly one row per detected occupancy unit, so
                // eta(2, 0) / y(2, 0) is the cell's aggregated cover predictor /
                // observation.
                const double y_pos   = y_cell.y(2, 0);
                const double eta_pos = etas.eta(2, 0);
                cell_ll += PosPolicy::log_density(y_pos, eta_pos, phi_pos);
                double g_pos = 0.0, h_pos = 0.0;
                PosPolicy::grad_hess_eta(y_pos, eta_pos, phi_pos,
                                         want_hess, g_pos, h_pos);
                out.arm_grad[2][0] = g_pos;
                if (want_hess) out.arm_neg_hess_diag[2][0] = h_pos;
            }
            // Cross-Hessians all zero in the det case (psi, p, pos arms
            // factorise: log psi + sum_v log h_v + (pos term), each depending
            // on disjoint etas). Buffers come zeroed.
            return cell_ll;
        }

        // ---- nodet case (family-independent: pos arm doesn't contribute) ----
        // The whole cell collapses to the single all-undetected occupancy
        // mixture L = psi P0 + (1 - psi). Score / curvature (Observed mixture
        // Hessian or complete-data Fisher) come from the shared
        // nodet_mixture_block with w = psi and the cell's Jc visits; psi -> arm
        // 0, the detection visits -> arm 1, and the (psi, p) / (p, p) cross
        // blocks land in arm_cross_hess[0][1] / [1][1]. The pos arm contributes
        // nothing (no detected visits) -- its buffers stay zeroed.
        const bool expected = (out.curvature == tulpa::CurvatureMode::Expected);
        std::vector<double> eta_p_buf(Jc);
        for (int v = 0; v < Jc; v++) eta_p_buf[v] = etas.eta(1, v);

        double* cross_w_p = (!expected && want_hess && out.arm_cross_hess
                             && out.arm_cross_hess[0] && out.arm_cross_hess[0][1])
                            ? out.arm_cross_hess[0][1] : nullptr;
        double* cross_p_p = (!expected && want_hess && out.arm_cross_hess
                             && out.arm_cross_hess[1] && out.arm_cross_hess[1][1])
                            ? out.arm_cross_hess[1][1] : nullptr;

        double g_psi = 0.0, nh_psi = 0.0;
        const double cell_ll = nodet_mixture_block(
            psi, eta_p_buf.data(), Jc, want_hess, expected,
            g_psi, nh_psi, out.arm_grad[1], out.arm_neg_hess_diag[1],
            cross_w_p, cross_p_p);
        out.arm_grad[0][0] = g_psi;
        if (want_hess) out.arm_neg_hess_diag[0][0] = nh_psi;
        // (psi, pos), (p, pos), (pos, pos) cross-Hessians: zero in nodet.
        return cell_ll;
    }

    std::string name() const override {
        return std::string(PosPolicy::spec_name()) + (Aggregated ? "_agg" : "");
    }

    bool thread_safe() const override { return true; }
};

typedef OccuCoverCoupling<LognormalPositive, false> OccuCoverLognormalCoupling;
typedef OccuCoverCoupling<BetaPositive,      false> OccuCoverBetaCoupling;
typedef OccuCoverCoupling<LognormalPositive, true>  OccuCoverLognormalAggCoupling;
typedef OccuCoverCoupling<BetaPositive,      true>  OccuCoverBetaAggCoupling;

} // namespace tulpaObs

#endif // TULPAOBS_CELL_COUPLING_OCCU_COVER_H
