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
#include <utility>
#include <vector>

namespace tulpaObs {

template <class PosPolicy, bool Aggregated = false>
class OccuCoverCoupling final : public tulpa::CellCouplingSpec {
public:
    std::vector<int> arm_ids() const override { return {0, 1, 2}; }

    // Dense cross-Hessian slabs this spec writes (arms: psi = 0, p = 1, pos = 2).
    // The det branch factorises -- every cross is zero -- and the nodet branch
    // writes only the (psi, p) cross (1 x J) plus the (p, p) block. The (p, p)
    // block rides the rank-1 self-cross when the engine supplies it (single
    // response); batched has no rank-1 path and needs it dense. The (p, pos),
    // (pos, pos) and (psi, pos) blocks are always zero. Omitting them is what
    // bounds a cell with J visits to O(J) instead of O(J^2) dense slabs -- the
    // all-undetected cells in a large grid would otherwise allocate a J x J
    // (p, p) block per cell.
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

        // Grad-only request (cached-factor reuse step): the kernel discards the
        // Hessian this call would produce, so skip every negative-Hessian and
        // cross-Hessian term -- including the positive arm's trigamma -- and
        // leave the pre-zeroed buffers untouched. The gradient and the returned
        // log-density are still exact.
        const bool want_hess = !out.grad_only;

        // Batched multi-response (gcol33/tulpa#66): the visit-design row is
        // species-invariant, so loop species INNER -- the kernel scatters each
        // species' derivatives into its own block. With B = 1 the loop runs
        // once at s = 0 and every (.., s) accessor / species-major write
        // [s * rc + j] reduces to the pre-batch path, byte-identical.
        const int base2_stride = out.n_rows_in_arm(2);
        double total_ll = 0.0;
        for (int s = 0; s < B; s++) {
            const double psi     = sigmoid_(etas.eta(0, 0, s));
            const double phi_pos = y_cell.phi(2, s);

            bool any_det = false;
            for (int v = 0; v < Jc; v++) {
                if (y_cell.y(1, v, s) > 0.5) { any_det = true; break; }
            }

            if (!any_det) {
                // nodet case (family-independent: pos arm doesn't contribute).
                total_ll += occu_nodet_block(psi, etas, y_cell, Jc, want_hess, out, s);
                continue;
            }

            // Occupancy + detection arms (shared with the latent spec); cross-
            // Hessians stay zero (the psi / p / pos arms factorise in the det
            // branch).
            double cell_ll = occu_det_psi_p_block(psi, etas, y_cell, Jc,
                                                  want_hess, out, s);
            const int base2 = s * base2_stride;

            if (!Aggregated) {
                // Per-visit cover: one log f_pos per detected visit, the pos
                // arm row aligned with the detection visit.
                for (int v = 0; v < Jc; v++) {
                    if (y_cell.y(1, v, s) <= 0.5) continue;
                    const double y_pos   = y_cell.y(2, v, s);
                    const double eta_pos = etas.eta(2, v, s);
                    cell_ll += PosPolicy::log_density(y_pos, eta_pos, phi_pos);
                    double g_pos = 0.0, h_pos = 0.0;
                    PosPolicy::grad_hess_eta(y_pos, eta_pos, phi_pos,
                                             want_hess, g_pos, h_pos);
                    out.arm_grad[2][base2 + v] = g_pos;
                    if (want_hess) out.arm_neg_hess_diag[2][base2 + v] = h_pos;
                }
            } else {
                // Cell-aggregated cover: a single log f_pos at the occupancy
                // unit's one pos row (the mean / median cover over its detected
                // visits), evaluated once because any_det holds here.
                const double y_pos   = y_cell.y(2, 0, s);
                const double eta_pos = etas.eta(2, 0, s);
                cell_ll += PosPolicy::log_density(y_pos, eta_pos, phi_pos);
                double g_pos = 0.0, h_pos = 0.0;
                PosPolicy::grad_hess_eta(y_pos, eta_pos, phi_pos,
                                         want_hess, g_pos, h_pos);
                out.arm_grad[2][base2] = g_pos;
                if (want_hess) out.arm_neg_hess_diag[2][base2] = h_pos;
            }
            total_ll += cell_ll;
        }
        return total_ll;
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
