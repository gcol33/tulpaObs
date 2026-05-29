// cell_coupling_occu_cover.h
// Stateless `CellCouplingSpec` implementing the per-cell log-density of
// the lognormal cover-hurdle (psi/p/pos arms) for the joint nested-Laplace
// path in tulpa (gcol33/tulpa#32 Layer B.2 consumer).
//
// Arm layout (kk indexes the spec's arm_ids() return):
//   kk = 0 -> psi arm: 1 row per cell, carries no y data, eta = logit psi_c
//   kk = 1 -> p   arm: J_c rows per cell, y(1, v) in {0, 1} detection,
//                       eta = logit p_cv
//   kk = 2 -> pos arm: J_c rows per cell, y(2, v) = raw y_pos > 0 at
//                       detected visits (0 elsewhere), eta = mu of log y,
//                       y_cell.phi(2) = lognormal SD on the log scale
//
// Per-cell density splits on `any_det = any(y_det_cv == 1)`:
//
//   det case (any_det):
//     log p_cell = log psi + sum_v [
//                    (y_det = 1) (log p_v + log f_pos(y_pos_v; eta_pos_v, sigma_pos))
//                  + (y_det = 0)  log (1 - p_v) ]
//
//   nodet case (all y_det = 0):
//     log p_cell = log L,   L = psi P0 + (1 - psi),   P0 = prod_v (1 - p_v)
//
// All closed-form first + second + cross derivatives below are FD-checked
// in tulpaObs/tests/testthat/test-occu-cover-coupling.R against numerical
// derivatives of the cell density above.

#ifndef TULPAOBS_CELL_COUPLING_OCCU_COVER_H
#define TULPAOBS_CELL_COUPLING_OCCU_COVER_H

#include <tulpa/cell_coupling.h>
#include <cmath>
#include <string>
#include <vector>

namespace tulpaObs {

inline double sigmoid_(double eta) {
    return 1.0 / (1.0 + std::exp(-eta));
}

inline double log_safe_(double x) {
    return (x > 0.0) ? std::log(x) : -1e300;
}

class OccuCoverLognormalCoupling final : public tulpa::CellCouplingSpec {
public:
    std::vector<int> arm_ids() const override { return {0, 1, 2}; }

    double evaluate_cell(int                       /*cell_idx*/,
                         const tulpa::CellEtas&     etas,
                         const tulpa::CellResponse& y_cell,
                         tulpa::CellDerivs&         out) const override {
        const int Jc = etas.n_rows_in_arm(1);
        const double eta_psi = etas.eta(0, 0);
        const double psi     = sigmoid_(eta_psi);
        const double sigma_pos = y_cell.phi(2);
        const double sigma2    = sigma_pos * sigma_pos;
        const double inv_sigma2 = (sigma2 > 0.0) ? 1.0 / sigma2 : 0.0;
        const double log_sigma  = log_safe_(sigma_pos);
        const double log_2pi    = std::log(2.0 * M_PI);

        bool any_det = false;
        for (int v = 0; v < Jc; v++) {
            if (y_cell.y(1, v) > 0.5) { any_det = true; break; }
        }

        if (any_det) {
            double cell_ll = log_safe_(psi);

            out.arm_grad[0][0]          = 1.0 - psi;
            out.arm_neg_hess_diag[0][0] = psi * (1.0 - psi);

            for (int v = 0; v < Jc; v++) {
                const double eta_p = etas.eta(1, v);
                const double p_v   = sigmoid_(eta_p);
                const double y_det = y_cell.y(1, v);

                if (y_det > 0.5) {
                    cell_ll += log_safe_(p_v);
                    out.arm_grad[1][v]          = 1.0 - p_v;
                    out.arm_neg_hess_diag[1][v] = p_v * (1.0 - p_v);

                    const double y_pos     = y_cell.y(2, v);
                    const double log_y_pos = log_safe_(y_pos);
                    const double eta_pos   = etas.eta(2, v);
                    const double r         = log_y_pos - eta_pos;
                    cell_ll += -log_y_pos - log_sigma - 0.5 * log_2pi
                               - 0.5 * r * r * inv_sigma2;
                    out.arm_grad[2][v]          = r * inv_sigma2;
                    out.arm_neg_hess_diag[2][v] = inv_sigma2;
                } else {
                    cell_ll += log_safe_(1.0 - p_v);
                    out.arm_grad[1][v]          = -p_v;
                    out.arm_neg_hess_diag[1][v] = p_v * (1.0 - p_v);
                }
            }
            // Cross-Hessians all zero in the det case (psi, p, pos arms
            // factorise: log psi + sum_v log h_v, with each log h_v
            // depending on disjoint etas). Buffers come zeroed.
            return cell_ll;
        }

        // ---- nodet case ----
        double log_P0 = 0.0;
        std::vector<double> p_cache(Jc);
        for (int v = 0; v < Jc; v++) {
            p_cache[v] = sigmoid_(etas.eta(1, v));
            log_P0    += log_safe_(1.0 - p_cache[v]);
        }
        const double P0       = std::exp(log_P0);
        const double L        = psi * P0 + (1.0 - psi);
        const double inv_L    = (L > 0.0) ? (1.0 / L) : 0.0;
        const double inv_L2   = inv_L * inv_L;
        const double one_m_P0 = 1.0 - P0;
        const double psi_1mp  = psi * (1.0 - psi);

        // psi: grad and neg-hess diagonal
        const double g_psi = -psi_1mp * one_m_P0 * inv_L;
        out.arm_grad[0][0] = g_psi;
        // -d^2/d eta_psi^2 = psi(1-psi)(1-2psi)(1-P0)/L + g_psi^2
        out.arm_neg_hess_diag[0][0] =
            psi_1mp * (1.0 - 2.0 * psi) * one_m_P0 * inv_L
            + g_psi * g_psi;

        // p arm: per-visit grad + neg-hess diagonal
        for (int v = 0; v < Jc; v++) {
            const double p_v        = p_cache[v];
            const double psi_P0_p   = psi * P0 * p_v;
            const double g_p_v      = -psi_P0_p * inv_L;
            out.arm_grad[1][v] = g_p_v;
            // -d^2/d eta_p_v^2 = psi*P0*p_v*(1 - 2 p_v)/L + g_p_v^2
            out.arm_neg_hess_diag[1][v] =
                psi_P0_p * (1.0 - 2.0 * p_v) * inv_L
                + g_p_v * g_p_v;
        }

        // pos arm contributes nothing in the nodet case (no detected
        // visits). Gradient and neg-hess buffers already zeroed.

        // Cross-Hessian (psi, p_v): -d^2/d eta_psi d eta_p_v = P0 p_v psi(1-psi)/L^2
        if (out.arm_cross_hess && out.arm_cross_hess[0]
            && out.arm_cross_hess[0][1]) {
            for (int v = 0; v < Jc; v++) {
                out.arm_cross_hess[0][1][v] =
                    P0 * p_cache[v] * psi_1mp * inv_L2;
            }
        }
        // Cross-Hessian (p_v, p_w) for v != w in nodet:
        //   -d^2/d eta_p_v d eta_p_w = -psi P0 p_v p_w (1-psi) / L^2
        if (out.arm_cross_hess && out.arm_cross_hess[1]
            && out.arm_cross_hess[1][1]) {
            const double a = -psi * P0 * (1.0 - psi) * inv_L2;
            for (int v = 0; v < Jc; v++) {
                for (int w = v + 1; w < Jc; w++) {
                    const double val = a * p_cache[v] * p_cache[w];
                    out.arm_cross_hess[1][1][(std::size_t)v * Jc + w] = val;
                    out.arm_cross_hess[1][1][(std::size_t)w * Jc + v] = val;
                }
            }
        }
        // (psi, pos), (p, pos), (pos, pos) cross-Hessians: zero in nodet.

        return log_safe_(L);
    }

    std::string name() const override {
        return std::string("occu_cover_lognormal");
    }

    bool thread_safe() const override { return true; }
};

} // namespace tulpaObs

#endif // TULPAOBS_CELL_COUPLING_OCCU_COVER_H
