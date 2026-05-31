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
// The templated `OccuCoverCoupling<PosPolicy>` factors that out; each
// concrete spec (`occu_cover_lognormal`, `occu_cover_beta`) is a typedef
// over a positive-arm policy supplying log_density / grad_eta /
// neg_hess_eta.
//
// All closed-form first + second + cross derivatives below are FD-checked
// in tulpaObs/tests/testthat/test-occu-cover-coupling.R against numerical
// derivatives of the cell density above.

#ifndef TULPAOBS_CELL_COUPLING_OCCU_COVER_H
#define TULPAOBS_CELL_COUPLING_OCCU_COVER_H

#include <tulpa/cell_coupling.h>
#include <tulpa/portable_math.h>   // portable_digamma / portable_trigamma (R-free, inlinable)
#include <Rcpp.h>                  // M_PI and R headers (no R::digamma in the hot path now)
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


// ---------------------------------------------------------------------------
// Positive-arm policies: log-density and first/second derivatives wrt the
// arm's linear predictor `eta_pos`. The cell-coupling spec is templated on
// one of these so the det-branch wiring stays single-source.
// ---------------------------------------------------------------------------

struct LognormalPositive {
    static constexpr const char* spec_name() { return "occu_cover_lognormal"; }

    // log f(y; mu = eta_pos, sigma = phi) on the natural scale:
    //   -log y - log sigma - 0.5 log(2 pi) - 0.5 ((log y - eta) / sigma)^2
    static double log_density(double y_pos, double eta_pos, double phi) {
        const double log_y     = log_safe_(y_pos);
        const double sigma     = phi;
        const double log_sigma = log_safe_(sigma);
        const double inv_s2    = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
        const double r         = log_y - eta_pos;
        return -log_y - log_sigma - 0.5 * std::log(2.0 * M_PI)
               - 0.5 * r * r * inv_s2;
    }
    // d log f / d eta_pos = (log y - eta) / sigma^2
    static double grad_eta(double y_pos, double eta_pos, double phi) {
        double g = 0.0, h = 0.0;
        grad_hess_eta(y_pos, eta_pos, phi, false, g, h);
        return g;
    }
    // -d^2 log f / d eta_pos^2 = 1 / sigma^2
    static double neg_hess_eta(double y_pos, double eta_pos, double phi) {
        double g = 0.0, h = 0.0;
        grad_hess_eta(y_pos, eta_pos, phi, true, g, h);
        return h;
    }
    // Combined eta-gradient and (constant) negative Hessian in one pass.
    // `want_hess` gates the curvature; a grad-only inner-Newton step
    // (cached-factor reuse) passes false. Single source for grad_eta /
    // neg_hess_eta and the cell-coupling hot path.
    static void grad_hess_eta(double y_pos, double eta_pos, double phi,
                              bool want_hess, double& grad, double& neg_hess) {
        const double sigma  = phi;
        const double inv_s2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
        grad = (log_safe_(y_pos) - eta_pos) * inv_s2;
        if (want_hess) neg_hess = inv_s2;
    }
};


struct BetaPositive {
    static constexpr const char* spec_name() { return "occu_cover_beta"; }

    // log f(y; mu = sigmoid(eta_pos), phi = precision) on (0, 1):
    //   lgamma(phi) - lgamma(mu phi) - lgamma((1-mu) phi)
    //   + (mu phi - 1) log y + ((1 - mu) phi - 1) log(1 - y)
    static double log_density(double y_pos, double eta_pos, double phi) {
        const double mu      = sigmoid_(eta_pos);
        const double log_y   = log_safe_(y_pos);
        const double log_1my = log_safe_(1.0 - y_pos);
        return std::lgamma(phi) - std::lgamma(mu * phi)
             - std::lgamma((1.0 - mu) * phi)
             + (mu * phi - 1.0) * log_y
             + ((1.0 - mu) * phi - 1.0) * log_1my;
    }
    // d log f / d eta = phi * mu (1 - mu) * g
    //   where g = -digamma(mu phi) + digamma((1-mu) phi) + log y - log(1-y).
    static double grad_eta(double y_pos, double eta_pos, double phi) {
        double g = 0.0, h = 0.0;
        grad_hess_eta(y_pos, eta_pos, phi, false, g, h);
        return g;
    }
    // -d^2 log f / d eta^2 = phi^2 mu^2 (1-mu)^2 (trigamma(mu phi)
    //                          + trigamma((1-mu) phi))
    //                       - phi mu (1-mu) g (1 - 2 mu)
    static double neg_hess_eta(double y_pos, double eta_pos, double phi) {
        double g = 0.0, h = 0.0;
        grad_hess_eta(y_pos, eta_pos, phi, true, g, h);
        return h;
    }
    // Combined eta-gradient + negative Hessian in one pass. The digamma terms
    // are shared between the score and the curvature (computed once); the
    // trigamma terms are computed ONLY when `want_hess` is true, so a grad-only
    // inner-Newton step (cached-factor reuse) skips them entirely. Uses tulpa's
    // portable, inlinable, OpenMP-safe digamma/trigamma rather than the R math
    // library calls. Single source for grad_eta / neg_hess_eta and the
    // cell-coupling hot path.
    static void grad_hess_eta(double y_pos, double eta_pos, double phi,
                              bool want_hess, double& grad, double& neg_hess) {
        const double mu      = sigmoid_(eta_pos);
        const double log_y   = log_safe_(y_pos);
        const double log_1my = log_safe_(1.0 - y_pos);
        const double a       = mu * phi;
        const double b       = (1.0 - mu) * phi;
        const double g       = -tulpa::math::portable_digamma(a)
                              +  tulpa::math::portable_digamma(b)
                              + log_y - log_1my;
        const double m1m     = mu * (1.0 - mu);
        grad = phi * m1m * g;
        if (want_hess) {
            const double trig = tulpa::math::portable_trigamma(a)
                              + tulpa::math::portable_trigamma(b);
            neg_hess = phi * phi * m1m * m1m * trig
                     - phi * m1m * g * (1.0 - 2.0 * mu);
        }
    }
};


template <class PosPolicy>
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

                    const double y_pos   = y_cell.y(2, v);
                    const double eta_pos = etas.eta(2, v);
                    cell_ll += PosPolicy::log_density(y_pos, eta_pos, phi_pos);
                    double g_pos = 0.0, h_pos = 0.0;
                    PosPolicy::grad_hess_eta(y_pos, eta_pos, phi_pos,
                                             want_hess, g_pos, h_pos);
                    out.arm_grad[2][v] = g_pos;
                    if (want_hess) out.arm_neg_hess_diag[2][v] = h_pos;
                } else {
                    cell_ll += log_safe_(1.0 - p_v);
                    out.arm_grad[1][v]          = -p_v;
                    if (want_hess) out.arm_neg_hess_diag[1][v] = p_v * (1.0 - p_v);
                }
            }
            // Cross-Hessians all zero in the det case (psi, p, pos arms
            // factorise: log psi + sum_v log h_v, with each log h_v
            // depending on disjoint etas). Buffers come zeroed.
            return cell_ll;
        }

        // ---- nodet case (family-independent: pos arm doesn't contribute) ----
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

        // Gradients are identical in both curvature modes (Fisher scoring
        // changes the curvature, not the score). Only the neg-Hessian and the
        // cross-Hessian differ: Observed writes the true mixture Hessian (which
        // is indefinite away from the mode); Expected writes the complete-data
        // expected information (Fisher), which is block-diagonal and PSD by
        // construction. gamma_c = P(z_c = 1 | all undetected) = psi*P0/L is the
        // occupancy responsibility that weights the detection arm's information.
        const bool expected = (out.curvature == tulpa::CurvatureMode::Expected);
        const double gamma_c = (L > 0.0) ? (psi * P0 * inv_L) : 0.0;

        // psi: grad and neg-hess diagonal
        const double g_psi = -psi_1mp * one_m_P0 * inv_L;
        out.arm_grad[0][0] = g_psi;
        if (want_hess) {
            out.arm_neg_hess_diag[0][0] = expected
                ? psi_1mp                                        // complete-data Fisher: psi(1-psi)
                // -d^2/d eta_psi^2 = psi(1-psi)(1-2psi)(1-P0)/L + g_psi^2
                : psi_1mp * (1.0 - 2.0 * psi) * one_m_P0 * inv_L + g_psi * g_psi;
        }

        // p arm: per-visit grad + neg-hess diagonal
        for (int v = 0; v < Jc; v++) {
            const double p_v        = p_cache[v];
            const double psi_P0_p   = psi * P0 * p_v;
            const double g_p_v      = -psi_P0_p * inv_L;
            out.arm_grad[1][v] = g_p_v;
            if (want_hess) {
                out.arm_neg_hess_diag[1][v] = expected
                    ? gamma_c * p_v * (1.0 - p_v)                // complete-data Fisher
                    // -d^2/d eta_p_v^2 = psi*P0*p_v*(1 - 2 p_v)/L + g_p_v^2
                    : psi_P0_p * (1.0 - 2.0 * p_v) * inv_L + g_p_v * g_p_v;
            }
        }

        // pos arm contributes nothing in the nodet case (no detected
        // visits). Gradient and neg-hess buffers already zeroed.

        // Cross-Hessians: the complete-data Fisher (Expected) is block-diagonal
        // -- given z the arms are independent GLMs -- so it writes none. The
        // observed (Expected = false) cross terms are the missing-information
        // contribution that makes the mixture Hessian indefinite. A grad-only
        // step writes none either (the kernel discards them).
        if (!expected && want_hess) {
            // (psi, p_v): -d^2/d eta_psi d eta_p_v = P0 p_v psi(1-psi)/L^2
            if (out.arm_cross_hess && out.arm_cross_hess[0]
                && out.arm_cross_hess[0][1]) {
                for (int v = 0; v < Jc; v++) {
                    out.arm_cross_hess[0][1][v] =
                        P0 * p_cache[v] * psi_1mp * inv_L2;
                }
            }
            // (p_v, p_w) for v != w:
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
        }
        // (psi, pos), (p, pos), (pos, pos) cross-Hessians: zero in nodet.

        return log_safe_(L);
    }

    std::string name() const override {
        return std::string(PosPolicy::spec_name());
    }

    bool thread_safe() const override { return true; }
};

typedef OccuCoverCoupling<LognormalPositive> OccuCoverLognormalCoupling;
typedef OccuCoverCoupling<BetaPositive>      OccuCoverBetaCoupling;

} // namespace tulpaObs

#endif // TULPAOBS_CELL_COUPLING_OCCU_COVER_H
