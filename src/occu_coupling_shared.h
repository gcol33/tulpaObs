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

// Per-row detection multiplicity (n_trials) for the detection arm (index 1),
// with a safe fallback of 1 when the arm carries no trial-count buffer (e.g. the
// batched fused path, whose response list may leave arm_n_trials null). Weight 1
// is the uncompressed path and is byte-identical to the pre-compression kernel.
inline double occu_det_weight_(const tulpa::CellResponse& y_cell, int j) {
    return (y_cell.arm_n_trials && y_cell.arm_n_trials[1])
           ? static_cast<double>(y_cell.n_trials(1, j)) : 1.0;
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
    // d log f / d log_sigma = (log y - eta)^2 / sigma^2 - 1 = rstd^2 - 1. The
    // analytic dispersion score the NUTS targets need (the Laplace / coupling
    // paths grid-integrate sigma and never differentiate it). Single source for
    // the per-obs lognormal log-dispersion gradient (occu_cover_nuts.cpp and the
    // spatial-factor community sampler ms_occu_cover_spatial_nuts.cpp).
    static double grad_logdisp(double y_pos, double eta_pos, double phi) {
        const double sigma = phi;
        const double r     = (sigma > 0.0) ? (log_safe_(y_pos) - eta_pos) / sigma : 0.0;
        return r * r - 1.0;
    }
};


// Identity-link Gaussian positive arm (gcol33/tulpaObs#112): the delta-normal
// hurdle's magnitude part, for a positive response already on a real,
// unbounded scale. It is LognormalPositive with the response taken as-is (no
// log transform) and no change-of-variable Jacobian: residual is (y - eta),
// mean on the response scale is mu = eta. Dispersion is the SD sigma (= phi),
// shared handling with the lognormal arm everywhere outside the density.
struct GaussianPositive {
    static constexpr const char* spec_name()            { return "occu_cover_gaussian"; }
    static constexpr const char* multiscale_spec_name() { return "occu_multiscale_cover_gaussian"; }

    // log f(y; mu = eta_pos, sigma = phi) on the natural scale:
    //   -log sigma - 0.5 log(2 pi) - 0.5 ((y - eta) / sigma)^2
    static double log_density(double y_pos, double eta_pos, double phi) {
        const double sigma     = phi;
        const double log_sigma = log_safe_(sigma);
        const double inv_s2    = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
        const double r         = y_pos - eta_pos;
        return -log_sigma - 0.5 * std::log(2.0 * M_PI)
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
    // d log f / d eta_pos = (y - eta) / sigma^2; -d^2 = 1 / sigma^2.
    static void grad_hess_eta(double y_pos, double eta_pos, double phi,
                              bool want_hess, double& grad, double& neg_hess) {
        const double sigma  = phi;
        const double inv_s2 = (sigma > 0.0) ? 1.0 / (sigma * sigma) : 0.0;
        grad = (y_pos - eta_pos) * inv_s2;
        if (want_hess) neg_hess = inv_s2;
    }
    // d log f / d log_sigma = ((y - eta) / sigma)^2 - 1 = rstd^2 - 1.
    static double grad_logdisp(double y_pos, double eta_pos, double phi) {
        const double sigma = phi;
        const double r     = (sigma > 0.0) ? (y_pos - eta_pos) / sigma : 0.0;
        return r * r - 1.0;
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
    // The first term is the (always-positive) expected information; the second
    // is the data-dependent part that can drive the observed Hessian indefinite.
    // When `fisher` is non-null it receives that first term -- the per-obs Fisher
    // information -- so the latent marginal can build a PSD curvature without
    // recomputing the digamma / trigamma evaluation (single source of truth).
    static void grad_hess_eta(double y_pos, double eta_pos, double phi,
                              bool want_hess, double& grad, double& neg_hess,
                              double* fisher = nullptr) {
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
        if (want_hess || fisher) {
            const double trig = tulpa::math::portable_trigamma(a)
                              + tulpa::math::portable_trigamma(b);
            const double fish = phi * phi * m1m * m1m * trig;
            if (fisher) *fisher = fish;
            if (want_hess) neg_hess = fish - phi * m1m * g * (1.0 - 2.0 * mu);
        }
    }
    // d log f / d log_phi = phi [ digamma(phi) - mu digamma(mu phi)
    //   - (1-mu) digamma((1-mu) phi) + mu log y + (1-mu) log(1-y) ]. The analytic
    // dispersion (precision) score the NUTS targets need; single source for the
    // per-obs beta log-precision gradient (occu_cover_nuts.cpp and
    // ms_occu_cover_spatial_nuts.cpp).
    static double grad_logdisp(double y_pos, double eta_pos, double phi) {
        const double mu    = sigmoid_(eta_pos);
        const double a     = mu * phi;
        const double b     = (1.0 - mu) * phi;
        const double ly    = log_safe_(y_pos);
        const double l1my  = log_safe_(1.0 - y_pos);
        return phi * (tulpa::math::portable_digamma(phi)
                      - mu * tulpa::math::portable_digamma(a)
                      - (1.0 - mu) * tulpa::math::portable_digamma(b)
                      + mu * ly + (1.0 - mu) * l1my);
    }
};


// ---------------------------------------------------------------------------
// Positive-arm code dispatch. The cover positive arm is selected at runtime by
// an integer code shared with the coupling / ploglik layers:
//   0 = lognormal, 3 = beta, 4 = gaussian (gcol33/tulpaObs#112).
// Single source for the NUTS targets (cover_nuts, occu_cover_nuts, the
// community spatial-factor sampler, multiscale) so no positive-density branch
// is copy-pasted across the samplers -- each calls these and the compiler
// inlines the switch.
inline double pos_log_density(int code, double y, double eta, double phi) {
    switch (code) {
        case 3:  return BetaPositive::log_density(y, eta, phi);
        case 4:  return GaussianPositive::log_density(y, eta, phi);
        default: return LognormalPositive::log_density(y, eta, phi);
    }
}
inline double pos_grad_eta(int code, double y, double eta, double phi) {
    switch (code) {
        case 3:  return BetaPositive::grad_eta(y, eta, phi);
        case 4:  return GaussianPositive::grad_eta(y, eta, phi);
        default: return LognormalPositive::grad_eta(y, eta, phi);
    }
}
inline double pos_grad_logdisp(int code, double y, double eta, double phi) {
    switch (code) {
        case 3:  return BetaPositive::grad_logdisp(y, eta, phi);
        case 4:  return GaussianPositive::grad_logdisp(y, eta, phi);
        default: return LognormalPositive::grad_logdisp(y, eta, phi);
    }
}


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
//
// Rank-1 (p, p) emission (gcol33/tulpaObs#94): the observed (p, p) off-diagonal
// cross is exactly the rank-1 a * p p^T, a = -w(1-w) P0 / L^2. When
// `rank1_coef_out` and `rank1_p_out` are both non-null (and the curvature is
// Observed) the block writes (a, p) there for the engine's rank-1 self-cross
// path INSTEAD of the dense `cross_p_p`, and folds the rank-1's own diagonal
// a p_v^2 into `nh_p` (storing the true diagonal minus a p_v^2) so the engine
// adds the full a p p^T. The two paths are mutually exclusive; pass
// `cross_p_p = nullptr` when requesting rank-1. With both rank-1 pointers null
// the dense `cross_p_p` path is byte-identical to before.
// ---------------------------------------------------------------------------
// `wt` (gcol33/tulpaObs detection-pattern compression) is an optional per-visit
// integer multiplicity: row v stands for wt[v] exchangeable non-detected visits
// that share this detection row (identical eta_p[v]). Passing nullptr (or all
// ones) is the uncompressed path and is byte-identical to before. The reduction
// is exact: P0 = prod_v (1 - p_v)^{wt_v}, and the mixture's score / observed
// Hessian in the compressed (unique-row) basis is S^T H S with S the row->pattern
// selection, which collapses to per-row weight factors on g_p / cross_w_p and to
// the rank-1 vector wt_v * p_v (see the derivation below and
// dev_notes/occu_multiscale_cover_derivation.md). Only all-undetected rows are
// ever compressed (the caller keeps detected visits individual for their cover),
// so wt applies uniformly across this block's rows.
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
    double*       p_out = nullptr,
    double*       rank1_coef_out = nullptr,
    double*       rank1_p_out = nullptr,
    const double* wt = nullptr)
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
        const double wv = wt ? wt[v] : 1.0;
        log_P0 += wv * log_safe_(1.0 - p[v]);
    }
    const double P0     = std::exp(log_P0);
    const double L      = w * P0 + (1.0 - w);
    const double inv_L  = (L > 0.0) ? (1.0 / L) : 0.0;
    const double inv_L2 = inv_L * inv_L;
    const double one_m_P0 = 1.0 - P0;
    const double w_1mw  = w * (1.0 - w);

    // Scores. Per-visit score scales linearly with the row multiplicity wt_v:
    // d log L / d eta_v aggregates the wt_v identical visits' contributions.
    g_w = -w_1mw * one_m_P0 * inv_L;
    for (int v = 0; v < nv; v++) {
        const double wv = wt ? wt[v] : 1.0;
        g_p[v] = wv * (-w * P0 * p[v] * inv_L);
    }
    if (!want_hess) return log_safe_(L);

    if (expected) {
        // Complete-data Fisher: block-diagonal, responsibility-weighted; the
        // wt_v identical visits contribute wt_v times the single-visit Fisher.
        const double gamma = (L > 0.0) ? (w * P0 * inv_L) : 0.0;
        nh_w = w_1mw;
        for (int v = 0; v < nv; v++) {
            const double wv = wt ? wt[v] : 1.0;
            nh_p[v] = wv * gamma * p[v] * (1.0 - p[v]);
        }
        return log_safe_(L);
    }

    // Observed (true mixture) Hessian. The reduced block is S^T H S with H the
    // full expanded V x V observed Hessian = diag(Hs_v) + a p p^T. For a unique
    // row v of multiplicity wt_v: the diagonal own-term aggregates to wt_v * Hs_v
    // and the rank-1 vector entry aggregates to wt_v * p_v (so the engine's
    // a (wt.p)(wt.p)^T reproduces both the wt_v(wt_v-1) intra-row and the
    // inter-row cross pairs). `Hs_v` is the single-visit observed diagonal and
    // `s` the single-visit score (NOT the wt-scaled g_p[v]).
    nh_w = w_1mw * (1.0 - 2.0 * w) * one_m_P0 * inv_L + g_w * g_w;
    const double a = -w * P0 * (1.0 - w) * inv_L2;
    const bool use_rank1 = (rank1_coef_out && rank1_p_out);
    for (int v = 0; v < nv; v++) {
        const double wv     = wt ? wt[v] : 1.0;
        const double s      = -w * P0 * p[v] * inv_L;
        const double Hs     = w * P0 * p[v] * (1.0 - 2.0 * p[v]) * inv_L + s * s;
        if (use_rank1) {
            // Rank-1 path: engine adds a * (wt.p)(wt.p)^T; fold the row's own
            // rank-1 diagonal out of nh_p (wv * (Hs - a p^2), so nh_p + a(wt p)^2
            // reconstructs the full reduced (v, v) entry wv*Hs + wv(wv-1) a p^2).
            nh_p[v] = wv * (Hs - a * p[v] * p[v]);
        } else {
            // Dense path: full reduced diagonal here, off-diagonal in cross_p_p.
            nh_p[v] = wv * Hs + wv * (wv - 1.0) * a * p[v] * p[v];
        }
    }
    if (cross_w_p) {
        for (int v = 0; v < nv; v++) {
            const double wv = wt ? wt[v] : 1.0;
            cross_w_p[v] = wv * P0 * p[v] * w_1mw * inv_L2;
        }
    }
    if (use_rank1) {
        // Rank-1 emission: hand the engine (a, wt.p); the fold above already
        // removed each row's own rank-1 diagonal from nh_p.
        *rank1_coef_out = a;
        for (int v = 0; v < nv; v++) {
            const double wv = wt ? wt[v] : 1.0;
            rank1_p_out[v] = wv * p[v];
        }
    } else if (cross_p_p) {
        for (int v = 0; v < nv; v++) {
            const double wv = wt ? wt[v] : 1.0;
            for (int wvv = v + 1; wvv < nv; wvv++) {
                const double ww  = wt ? wt[wvv] : 1.0;
                const double val = a * (wv * p[v]) * (ww * p[wvv]);
                cross_p_p[(std::size_t)v * nv + wvv] = val;
                cross_p_p[(std::size_t)wvv * nv + v] = val;
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
// `s` is the batch (species) index (gcol33/tulpa#66). It offsets the reads
// (eta / y species column) and the writes (the rc * n_batch species-major
// buffer slot [s * rc + j]). s = 0 with a B=1 CellEtas/Derivs is byte-identical
// to the pre-batch path, so the existing per-visit / latent / multiscale
// callers (which omit s) are unchanged.
inline double occu_det_psi_p_block(double                     psi,
                                   const tulpa::CellEtas&     etas,
                                   const tulpa::CellResponse& y_cell,
                                   int                        Jc,
                                   bool                       want_hess,
                                   tulpa::CellDerivs&         out,
                                   int                        s = 0) {
    const int base0 = s * out.n_rows_in_arm(0);
    const int base1 = s * out.n_rows_in_arm(1);
    double cell_ll = log_safe_(psi);
    out.arm_grad[0][base0] = 1.0 - psi;
    if (want_hess) out.arm_neg_hess_diag[0][base0] = psi * (1.0 - psi);
    for (int v = 0; v < Jc; ++v) {
        const double p_v   = sigmoid_(etas.eta(1, v, s));
        const double y_det = y_cell.y(1, v, s);
        // Detection-pattern compression: a non-detected row can stand for wv
        // exchangeable visits sharing this detection row (n_trials = wv). Detected
        // rows are always individual (wv = 1) so their per-visit cover stays
        // aligned. wv = 1 (the uncompressed path) is byte-identical to before.
        const double wv    = occu_det_weight_(y_cell, v);
        if (y_det > 0.5) {
            cell_ll += log_safe_(p_v);
            out.arm_grad[1][base1 + v] = 1.0 - p_v;
        } else {
            cell_ll += wv * log_safe_(1.0 - p_v);
            out.arm_grad[1][base1 + v] = -wv * p_v;
        }
        if (want_hess) out.arm_neg_hess_diag[1][base1 + v] = wv * p_v * (1.0 - p_v);
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
// `s` (gcol33/tulpa#66) offsets the species column of the reads and the
// species-major slot of the writes; s = 0 / B=1 is byte-identical, so the
// pre-batch callers (which omit s) are unchanged. The compact (psi, p) and
// (p, p) cross buffers are sliced to species s by the same s * (rc_k * rc_l)
// offset the kernel uses to lay them out.
inline double occu_nodet_block(double                     psi,
                               const tulpa::CellEtas&     etas,
                               const tulpa::CellResponse& y_cell,
                               int                        Jc,
                               bool                       want_hess,
                               tulpa::CellDerivs&         out,
                               int                        s = 0) {
    const bool expected = (out.curvature == tulpa::CurvatureMode::Expected);
    const int rc0   = out.n_rows_in_arm(0);
    const int rc1   = out.n_rows_in_arm(1);
    const int base0 = s * rc0;
    const int base1 = s * rc1;
    // eta and the per-row detection multiplicity (n_trials) for this cell's
    // rows. wt_v = 1 (the uncompressed path) leaves the mixture byte-identical.
    std::vector<double> eta_p_buf(Jc);
    std::vector<double> wt_buf(Jc);
    for (int v = 0; v < Jc; v++) {
        eta_p_buf[v] = etas.eta(1, v, s);
        wt_buf[v]    = occu_det_weight_(y_cell, v);
    }

    double* cross_w_p = (!expected && want_hess && out.arm_cross_hess
                         && out.arm_cross_hess[0] && out.arm_cross_hess[0][1])
                        ? out.arm_cross_hess[0][1] + (std::size_t) s * rc0 * rc1
                        : nullptr;

    // Detection (p, p) cross-Hessian: prefer the engine's rank-1 self-cross
    // path (gcol33/tulpaObs#94) on the single-response path, where the dense
    // V x V block is exactly rank-1; fall back to the dense cross_p_p buffer
    // otherwise (batched, or an engine that did not supply the descriptor).
    double* rank1_coef = nullptr;
    double* rank1_p    = nullptr;
    double* cross_p_p  = nullptr;
    if (!expected && want_hess && out.n_batch() == 1
        && out.arm_cross_rank1_coef && out.arm_cross_rank1_vec
        && out.arm_cross_rank1_vec[1]) {
        rank1_coef = &out.arm_cross_rank1_coef[1];
        rank1_p    = out.arm_cross_rank1_vec[1];
    } else {
        cross_p_p = (!expected && want_hess && out.arm_cross_hess
                     && out.arm_cross_hess[1] && out.arm_cross_hess[1][1])
                    ? out.arm_cross_hess[1][1] + (std::size_t) s * rc1 * rc1
                    : nullptr;
    }

    double g_psi = 0.0, nh_psi = 0.0;
    const double cell_ll = nodet_mixture_block(
        psi, eta_p_buf.data(), Jc, want_hess, expected,
        g_psi, nh_psi, out.arm_grad[1] + base1, out.arm_neg_hess_diag[1] + base1,
        cross_w_p, cross_p_p, /*p_out=*/nullptr,
        rank1_coef, rank1_p, wt_buf.data());
    out.arm_grad[0][base0] = g_psi;
    if (want_hess) out.arm_neg_hess_diag[0][base0] = nh_psi;
    return cell_ll;
}

} // namespace tulpaObs

#endif // TULPAOBS_OCCU_COUPLING_SHARED_H
