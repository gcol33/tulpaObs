// occu_coupling_shared.h
// Shared building blocks for the occu-cover cell-coupling specs (the 2-level
// `OccuCoverCoupling` in cell_coupling_occu_cover.h and the 3-level
// `OccuMultiscaleCoverCoupling` in cell_coupling_occu_multiscale_cover.h):
//
//   * sigmoid_ / log_safe_                     -- scalar link helpers
//   * LognormalPositive / BetaPositive         -- positive-arm policies
//                                                 (log-density + eta grad/hess)
//   * nodet_mixture_block                       -- closed-form score + curvature
//                                                 of one all-undetected
//                                                 occupancy/availability mixture
//
// The no-detection mixture below is the per-cell density of the 2-level
// occu_cover when no visit detects, and the per-PLOT density of a 3-level
// multiscale cell when a plot has no detection but the cell does. Both reduce
// to L = w * P0 + (1 - w) with w = sigmoid(eta_w) the occupancy/availability
// probability and P0 = prod_v (1 - p_v) the all-undetected probability given
// w = 1. Writing it once here keeps the two specs on a single source of truth;
// the derivation is in dev_notes/occu_multiscale_cover_derivation.md and every
// Observed derivative is FD-checked in the two coupling tests.

#ifndef TULPAOBS_OCCU_COUPLING_SHARED_H
#define TULPAOBS_OCCU_COUPLING_SHARED_H

#include <tulpa/portable_math.h>   // portable_digamma / portable_trigamma
#include <tulpa/cell_coupling.h>   // CellEtas / CellResponse / CellDerivs / CurvatureMode
#include <Rcpp.h>                  // M_PI
#include <cmath>
#include <cstddef>
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
// arm's linear predictor `eta_pos`. The cell-coupling specs are templated on
// one of these so the detected-visit cover wiring stays single-source.
// ---------------------------------------------------------------------------

struct LognormalPositive {
    static constexpr const char* spec_name()            { return "occu_cover_lognormal"; }
    static constexpr const char* multiscale_spec_name() { return "occu_multiscale_cover_lognormal"; }

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
    static double grad_eta(double y_pos, double eta_pos, double phi) {
        double g = 0.0, h = 0.0;
        grad_hess_eta(y_pos, eta_pos, phi, false, g, h);
        return g;
    }
    static double neg_hess_eta(double y_pos, double eta_pos, double phi) {
        double g = 0.0, h = 0.0;
        grad_hess_eta(y_pos, eta_pos, phi, true, g, h);
        return h;
    }
    // d log f / d eta_pos = (log y - eta) / sigma^2; -d^2 = 1 / sigma^2.
    static void grad_hess_eta(double y_pos, double eta_pos, double phi,
                              bool want_hess, double& grad, double& neg_hess) {
        const double sigma  = phi;
        const double inv_s2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
        grad = (log_safe_(y_pos) - eta_pos) * inv_s2;
        if (want_hess) neg_hess = inv_s2;
    }
};


struct BetaPositive {
    static constexpr const char* spec_name()            { return "occu_cover_beta"; }
    static constexpr const char* multiscale_spec_name() { return "occu_multiscale_cover_beta"; }

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
    static double grad_eta(double y_pos, double eta_pos, double phi) {
        double g = 0.0, h = 0.0;
        grad_hess_eta(y_pos, eta_pos, phi, false, g, h);
        return g;
    }
    static double neg_hess_eta(double y_pos, double eta_pos, double phi) {
        double g = 0.0, h = 0.0;
        grad_hess_eta(y_pos, eta_pos, phi, true, g, h);
        return h;
    }
    // d log f / d eta = phi mu(1-mu) g,
    //   g = -digamma(mu phi) + digamma((1-mu) phi) + log y - log(1-y).
    // -d^2 log f / d eta^2 = phi^2 mu^2 (1-mu)^2 (trigamma(mu phi)
    //                          + trigamma((1-mu) phi)) - phi mu(1-mu) g (1-2mu).
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


// ---------------------------------------------------------------------------
// nodet_mixture_block -- closed-form score + curvature of one all-undetected
// occupancy/availability mixture over `nv` visits:
//
//     L = w * P0 + (1 - w),   P0 = prod_v (1 - p_v),
//
// with w = sigmoid(eta_w) and p_v = sigmoid(eta_p[v]). Returns log L and
// writes (into compact length-nv local buffers the caller then places into
// its per-cell derivative buffers):
//
//   g_w, nh_w               score / neg-Hessian of eta_w
//   g_p[v], nh_p[v]         score / neg-Hessian of each eta_p_v
//   cross_w_p[v]            -d^2 logL / d eta_w d eta_p_v        (observed only)
//   cross_p_p[v*nv + w]     -d^2 logL / d eta_p_v d eta_p_w      (observed only,
//                                                                 v != w, both
//                                                                 triangles)
//
// Curvature flags:
//   want_hess == false -> only g_w / g_p written (grad-only inner step).
//   expected   == true -> complete-data Fisher (block diagonal):
//                         nh_w = w(1-w), nh_p[v] = gamma p_v(1-p_v) with
//                         gamma = w P0 / L; cross buffers left untouched.
//   expected   == false -> observed (true mixture) Hessian incl. cross terms.
//
// `cross_w_p` / `cross_p_p` may be nullptr when the caller does not need them
// (grad-only or Expected). The Observed compact forms use the identity
// L + w(1 - P0) = 1 to collapse cross(w, p) -> w(1-w) P0 p_v / L^2 and
// cross(p_v, p_w) -> -w(1-w) P0 p_v p_w / L^2 (see derivation note).
// ---------------------------------------------------------------------------
inline double nodet_mixture_block(
    double        w,
    const double* eta_p,
    int           nv,
    bool          want_hess,
    bool          expected,
    double&       g_w,
    double&       nh_w,
    double*       g_p,
    double*       nh_p,
    double*       cross_w_p,
    double*       cross_p_p,
    double*       p_out = nullptr)
{
    double log_P0 = 0.0;
    // Reuse the caller's p_out as the p cache when supplied; otherwise a small
    // stack buffer (nv is the number of visits in one plot/cell -- tiny).
    double  p_stack[64];
    double* p = (p_out != nullptr) ? p_out
              : (nv <= 64 ? p_stack : nullptr);
    std::vector<double> p_heap;
    if (p == nullptr) { p_heap.assign(nv, 0.0); p = p_heap.data(); }

    for (int v = 0; v < nv; v++) {
        p[v]    = sigmoid_(eta_p[v]);
        log_P0 += log_safe_(1.0 - p[v]);
    }
    const double P0     = std::exp(log_P0);
    const double L      = w * P0 + (1.0 - w);
    const double inv_L  = (L > 0.0) ? (1.0 / L) : 0.0;
    const double inv_L2 = inv_L * inv_L;
    const double one_m_P0 = 1.0 - P0;
    const double w_1mw  = w * (1.0 - w);

    // Scores (identical in both curvature modes).
    g_w = -w_1mw * one_m_P0 * inv_L;
    for (int v = 0; v < nv; v++) {
        g_p[v] = -w * P0 * p[v] * inv_L;
    }
    if (!want_hess) return log_safe_(L);

    if (expected) {
        // Complete-data Fisher: block-diagonal, responsibility-weighted.
        const double gamma = (L > 0.0) ? (w * P0 * inv_L) : 0.0;
        nh_w = w_1mw;
        for (int v = 0; v < nv; v++) {
            nh_p[v] = gamma * p[v] * (1.0 - p[v]);
        }
        return log_safe_(L);
    }

    // Observed (true mixture) Hessian.
    nh_w = w_1mw * (1.0 - 2.0 * w) * one_m_P0 * inv_L + g_w * g_w;
    for (int v = 0; v < nv; v++) {
        const double w_P0_p = w * P0 * p[v];
        nh_p[v] = w_P0_p * (1.0 - 2.0 * p[v]) * inv_L + g_p[v] * g_p[v];
    }
    if (cross_w_p) {
        for (int v = 0; v < nv; v++) {
            cross_w_p[v] = P0 * p[v] * w_1mw * inv_L2;
        }
    }
    if (cross_p_p) {
        const double a = -w * P0 * (1.0 - w) * inv_L2;
        for (int v = 0; v < nv; v++) {
            for (int wv = v + 1; wv < nv; wv++) {
                const double val = a * p[v] * p[wv];
                cross_p_p[(std::size_t)v * nv + wv] = val;
                cross_p_p[(std::size_t)wv * nv + v] = val;
            }
        }
    }
    return log_safe_(L);
}


// ---------------------------------------------------------------------------
// occu_det_psi_p_block -- shared det-branch occupancy + detection accumulation
// for the occu_cover specs (per-visit, aggregated, and latent). At a cell with
// >= 1 detection the occupancy / detection arms factorise from the cover arm,
// so this writes arm 0 (psi) and arm 1 (p) score / observed information and
// returns  log psi + sum_v [ y_det log p_v + (1 - y_det) log(1 - p_v) ]
// WITHOUT any cover (pos) contribution -- the caller adds its positive-arm
// term (one log f_pos per detected visit, one aggregated log f_pos, or one
// latent marginal log M). Cross-Hessians are zero in the det branch.
inline double occu_det_psi_p_block(double                     psi,
                                   const tulpa::CellEtas&     etas,
                                   const tulpa::CellResponse& y_cell,
                                   int                        Jc,
                                   bool                       want_hess,
                                   tulpa::CellDerivs&         out) {
    double cell_ll = log_safe_(psi);
    out.arm_grad[0][0] = 1.0 - psi;
    if (want_hess) out.arm_neg_hess_diag[0][0] = psi * (1.0 - psi);
    for (int v = 0; v < Jc; ++v) {
        const double p_v   = sigmoid_(etas.eta(1, v));
        const double y_det = y_cell.y(1, v);
        if (y_det > 0.5) {
            cell_ll += log_safe_(p_v);
            out.arm_grad[1][v] = 1.0 - p_v;
        } else {
            cell_ll += log_safe_(1.0 - p_v);
            out.arm_grad[1][v] = -p_v;
        }
        if (want_hess) out.arm_neg_hess_diag[1][v] = p_v * (1.0 - p_v);
    }
    return cell_ll;
}


// ---------------------------------------------------------------------------
// occu_nodet_block -- shared no-detection branch for the occu_cover specs. The
// whole cell collapses to the all-undetected occupancy mixture
// L = psi P0 + (1 - psi); the cover arm contributes nothing. Drives the shared
// nodet_mixture_block with w = psi over the cell's Jc visits, places the
// occupancy score / curvature into arm 0 and the detection terms into arm 1
// (incl. the (psi, p) and (p, p) cross-Hessian blocks under the Observed
// curvature), and returns log L.
inline double occu_nodet_block(double                     psi,
                               const tulpa::CellEtas&     etas,
                               int                        Jc,
                               bool                       want_hess,
                               tulpa::CellDerivs&         out) {
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
    return cell_ll;
}

} // namespace tulpaObs

#endif // TULPAOBS_OCCU_COUPLING_SHARED_H
