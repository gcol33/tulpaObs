// cell_coupling_occu_only.h
// Stateless `CellCouplingSpec` implementing the per-cell log-density of the
// single-season occupancy mixture (psi/p arms only, NO cover arm) for the joint
// nested-Laplace path in tulpa ( consumer).
//
// This is `OccuCoverCoupling` with the positive (cover) arm removed: the same
// site-level occupancy + detection mixture, fitted single-arm through the joint
// direct-grid engine so the shared areal field's hyperparameters integrate on
// the outer grid rather than through the EM fixed-point iteration that
// oscillates at EVA scale.
//
// Arm layout (kk indexes the spec's arm_ids() return):
//   kk = 0 -> psi arm: 1 row per cell (occupancy unit), no y data,
//                       eta = logit psi_c
//   kk = 1 -> p   arm: J_c rows per cell, y(1, v) in {0, 1} detection,
//                       eta = logit p_cv
//
// Per-cell density splits on `any_det = any(y_det_cv == 1)`:
//
//   det case (any_det):
//     log p_cell = log psi + sum_v [ (y_det = 1) log p_v
//                                  + (y_det = 0) log (1 - p_v) ]
//
//   nodet case (all y_det = 0):
//     log p_cell = log L,   L = psi P0 + (1 - psi),   P0 = prod_v (1 - p_v)
//
// Both branches reuse `occu_det_psi_p_block` / `occu_nodet_block` from
// occu_coupling_shared.h -- the SAME family-independent occupancy / detection
// derivatives the occu_cover specs use (single source of truth). The occupancy
// and detection arms factorise in the det branch, so cross-Hessians stay zero
// there; the nodet branch carries the (psi, p) and (p, p) cross-Hessian blocks
// under the Observed curvature exactly as occu_nodet_block writes them.
//
// All closed-form first + second + cross derivatives are FD-checked in
// tulpaObs/tests/testthat/test-occu-only-coupling.R against numerical
// derivatives of the cell density above (the cover arm dropped from the
// occu_cover coupling FD harness).

#ifndef TULPAOBS_CELL_COUPLING_OCCU_ONLY_H
#define TULPAOBS_CELL_COUPLING_OCCU_ONLY_H

#include <tulpa/cell_coupling.h>
#include "occu_coupling_shared.h"  // occu_det_psi_p_block / occu_nodet_block
#include <string>
#include <utility>
#include <vector>

namespace tulpaObs {

class OccuOnlyCoupling final : public tulpa::CellCouplingSpec {
public:
    std::vector<int> arm_ids() const override { return {0, 1}; }

    // Dense cross-Hessian slabs this spec writes (arms: psi = 0, p = 1). Same
    // set as OccuCoverCoupling with the cover arm dropped, because both run the
    // SAME occu_det_psi_p_block / occu_nodet_block helpers, and those write only
    // arm_cross_hess[0][1] and [1][1]. The det branch factorises (every cross is
    // zero); the nodet branch writes the (psi, p) cross (1 x J) plus the (p, p)
    // block, which rides the rank-1 self-cross when the engine supplies it and
    // needs a dense slab otherwise. Without this override the default allocates
    // every kk <= ll pair, so an all-undetected cell allocates a J x J (p, p)
    // slab it never writes -- O(J^2) per cell instead of O(J).
    std::vector<std::pair<int, int>> dense_cross_pairs(
            int /*n_coupled*/, bool rank1_self_supported) const override {
        if (rank1_self_supported) return {{0, 1}};
        return {{0, 1}, {1, 1}};
    }

    double evaluate_cell(int                       /*cell_idx*/,
                         const tulpa::CellEtas&     etas,
                         const tulpa::CellResponse& y_cell,
                         tulpa::CellDerivs&         out) const override {
        const int Jc = etas.n_rows_in_arm(1);
        const int B  = etas.n_batch();   // 1 (single species) or B (batched)

        const bool want_hess = !out.grad_only;

        double total_ll = 0.0;
        for (int s = 0; s < B; s++) {
            const double psi = sigmoid_(etas.eta(0, 0, s));

            bool any_det = false;
            for (int v = 0; v < Jc; v++) {
                if (y_cell.y(1, v, s) > 0.5) { any_det = true; break; }
            }

            if (!any_det) {
                total_ll += occu_nodet_block(psi, etas, y_cell, Jc, want_hess, out, s);
            } else {
                total_ll += occu_det_psi_p_block(psi, etas, y_cell, Jc,
                                                 want_hess, out, s);
            }
        }
        return total_ll;
    }

    std::string name() const override { return "occu_only"; }

    bool thread_safe() const override { return true; }
};

} // namespace tulpaObs

#endif // TULPAOBS_CELL_COUPLING_OCCU_ONLY_H
