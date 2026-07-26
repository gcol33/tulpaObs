// occu_cover_latent.h
// Per-site latent-cover marginal integrator for the latent-cover-per-unit
// occu_cover model (tulpaObs latent granularity).
//
// Generative structure (per occupancy unit i, conditional on z_i = 1):
//
//   u_i ~ N(0, sigma_u^2)                       one cover latent per unit
//   y_pos_ij | (detected j) ~ f_pos(eta_pos_i + u_i, disp2)
//
// where eta_pos_i is the UNIT-level cover linear predictor (constant across the
// unit's detected visits), `disp2` is the second (within-unit) dispersion --
// the lognormal residual SD sigma_eps or the beta precision phi -- pre-fit and
// fixed, and sigma_u is the integrated cover-latent SD (the outer-grid axis,
// carried on the pos arm's phi slot).
//
// The cover term contributed by unit i is the marginal over u_i:
//
//   log M_i = log integral_u [ prod_{j in det} f_pos(y_ij | eta + u, disp2) ]
//                              * N(u; 0, sigma_u^2) du
//
// Because eta is unit-level, log M_i is a function of the SINGLE scalar eta, so
// the cell-coupling kernel needs only a scalar score d logM/d eta and scalar
// observed information -d^2 logM/d eta^2 per unit -- no within-arm off-diagonal.
//
// Two policies share this interface:
//   * LognormalLatent -- closed form. log y_ij = eta + u + eps, eps ~ N(0,a),
//     u ~ N(0,b) (a = sigma_eps^2, b = sigma_u^2). The detected log-covers are
//     jointly N(eta, a I + b 11') (compound symmetry), whose log-density,
//     score and information have elementary closed forms in the per-unit
//     sufficient statistics (m, T1 = sum log y, T2 = sum (log y)^2). EXACT.
//   * BetaLatent -- no conjugacy; adaptive Gauss-Hermite over the scalar u_i,
//     reusing BetaPositive's per-observation density / eta-derivatives at the
//     shifted predictor eta + u. The marginal eta-derivatives use the standard
//     identities  d logM/d eta = E_post[s(u)],
//     -d^2 logM/d eta^2 = E_post[sum -d^2 log f/d pred^2] - Var_post(s(u)),
//     with s(u) = sum d log f/d pred and the posterior
//     post(u) ∝ exp(ell(u)) N(u; 0, b).
//
// All three quantities (log_m, score, neg_hess) are FD-checked against
// numerical derivatives of log M_i in test-occu-cover-latent.R.

#ifndef TULPAOBS_OCCU_COVER_LATENT_H
#define TULPAOBS_OCCU_COVER_LATENT_H

#include "occu_coupling_shared.h"   // LognormalPositive / BetaPositive / sigmoid_
#include <tulpa/gauss_hermite.h>    // tulpa::GaussHermite / gauss_hermite_prob
#include <Rcpp.h>                   // M_PI
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>

namespace tulpaObs {

// Result of integrating one unit's cover latent: the marginal log-density and
// its first two derivatives with respect to the unit-level cover predictor eta.
struct LatentMarginal {
    double log_m    = 0.0;   // log M_i
    double score    = 0.0;   //  d   log M_i / d eta
    double neg_hess = 0.0;   // marginal information: observed -d^2 log M_i/d eta^2,
                             // or the PSD Fisher form under CurvatureMode::Expected
};

// Per-unit lognormal sufficient statistics over the detected visits.
struct LnSuffStat {
    int    m  = 0;     // number of detected visits for the unit
    double t1 = 0.0;   // sum_j log y_ij
    double t2 = 0.0;   // sum_j (log y_ij)^2
};

inline LnSuffStat ln_suff_stat(const double* y, int m) {
    LnSuffStat s; s.m = m;
    for (int j = 0; j < m; ++j) {
        const double ly = log_safe(y[j]);   // consistent boundary with the other cover paths
        s.t1 += ly;
        s.t2 += ly * ly;
    }
    return s;
}


// ---------------------------------------------------------------------------
// LognormalLatent -- closed-form compound-symmetry marginal.
//
//   Sigma = a I_m + b 11',  a = sigma_eps^2,  b = sigma_u^2
//   r_j   = log y_ij - eta
//   log|Sigma|        = (m-1) log a + log(a + m b)
//   r' Sigma^{-1} r   = (1/a) [ S2c - (b/(a + m b)) S1^2 ]
//     with S1 = sum r_j = T1 - m eta,  S2c = sum r_j^2 = T2 - 2 eta T1 + m eta^2
//   log M = -m/2 log(2 pi) - 1/2 log|Sigma| - 1/2 r'Sigma^{-1}r - T1
//     (the trailing -T1 = -sum log y_ij is the lognormal change-of-variables
//      Jacobian)
//   d  logM / d eta   =  S1 / (a + m b)
//  -d^2 logM / d eta^2 =  m  / (a + m b)
// ---------------------------------------------------------------------------
struct LognormalLatent {
    static constexpr const char* spec_name() { return "occu_cover_lognormal_latent"; }
    typedef LnSuffStat SiteData;

    static SiteData make_site(const double* y, int m) { return ln_suff_stat(y, m); }

    // `disp2` is the within-unit residual SD sigma_eps; `gh` is unused. The
    // marginal is exactly quadratic, so observed == expected information and the
    // curvature mode is ignored (accepted only for interface parity with the
    // beta policy).
    static LatentMarginal marginal(const SiteData& s, double eta,
                                   double disp2, double sigma_u,
                                   const tulpa::GaussHermite& /*gh*/,
                                   tulpa::CurvatureMode /*curv*/
                                       = tulpa::CurvatureMode::Observed) {
        const int    m = s.m;
        const double a = disp2 * disp2;             // sigma_eps^2
        const double b = sigma_u * sigma_u;         // sigma_u^2
        const double denom = a + m * b;             // a + m b

        const double s1   = s.t1 - m * eta;                       // sum r_j
        const double s2c  = s.t2 - 2.0 * eta * s.t1 + m * eta * eta;  // sum r_j^2
        const double quad = (a > 0.0)
            ? (s2c - (b / denom) * s1 * s1) / a
            : 0.0;
        const double logdet = (m - 1) * std::log(a) + std::log(denom);

        LatentMarginal out;
        out.log_m    = -0.5 * m * std::log(2.0 * M_PI)
                       - 0.5 * logdet - 0.5 * quad - s.t1;
        out.score    = s1 / denom;
        out.neg_hess = (double) m / denom;
        return out;
    }
};


// ---------------------------------------------------------------------------
// BetaLatent -- adaptive Gauss-Hermite over the scalar cover latent.
//
// ell(u)   = sum_j log Beta(y_ij | sigmoid(eta + u), phi)
// s(u)     = sum_j d log Beta / d pred   (= ell'(u), since pred = eta + u)
// sneg(u)  = sum_j (-d^2 log Beta / d pred^2)   (= -ell''(u))
// q(u)     = ell(u) + log N(u; 0, b)   posterior log-kernel
//
// Mode-find q'(u) = s(u) - u/b = 0 by Newton (q'' = -sneg - 1/b < 0), then
// place adaptive GH nodes u_k = uhat + tau z_k, tau^2 = 1/(sneg(uhat) + 1/b).
// Posterior weights pi_k ∝ w_k exp(z_k^2/2) exp(ell(u_k)) N(u_k; 0, b).
//
//   log M           = log tau + 1/2 log(2 pi) + logsumexp_k(log_term_k)
//   d  logM / d eta = E_pi[s(u)]
//  -d^2 logM/d eta^2 = E_pi[sneg(u)] - Var_pi(s(u))     (Observed)
//
// The observed marginal information E_pi[sneg] - Var_pi(s) can go indefinite at
// extreme outer-grid cells (large sigma_u driving the beta mean toward 0/1),
// where it stalls the inner Newton and the cell returns log_m = -Inf. Under
// CurvatureMode::Expected the inner step instead uses the always-positive
// marginal Fisher information E_pi[sum_j fisher_beta(eta + u)] -- the per-obs
// Fisher term integrated over the latent posterior, dropping both the
// data-dependent part of the per-obs Hessian and the -Var_pi(s) term. This is
// PSD by construction, so those cells converge. The reported log_m / score are
// unchanged; only neg_hess differs, and the engine's final mode-pass always
// re-evaluates with Observed curvature, so log_marginal / SEs are unaffected.
// ---------------------------------------------------------------------------
struct BetaLatent {
    static constexpr const char* spec_name() { return "occu_cover_beta_latent"; }
    typedef std::vector<double> SiteData;

    static SiteData make_site(const double* y, int m) {
        return SiteData(y, y + m);
    }

    // ell(u), s(u) = ell'(u), sneg(u) = -ell''(u) at predictor eta + u. `sfish`
    // accumulates the per-obs Fisher information sum_j fisher_beta(eta + u): the
    // always-positive curvature used to build the Expected marginal Hessian.
    static void ell_and_derivs(const SiteData& y, double pred, double phi,
                               double& ell, double& s, double& sneg,
                               double& sfish) {
        ell = 0.0; s = 0.0; sneg = 0.0; sfish = 0.0;
        for (double yj : y) {
            ell += BetaPositive::log_density(yj, pred, phi);
            double g = 0.0, h = 0.0, f = 0.0;
            BetaPositive::grad_hess_eta(yj, pred, phi, true, g, h, &f);
            s     += g;   // d log Beta / d pred
            sneg  += h;   // -d^2 log Beta / d pred^2  (observed)
            sfish += f;   // Fisher information         (expected, PSD)
        }
    }

    // `disp2` is the beta precision phi; `gh` carries the probabilist GH rule.
    // `curv` selects the marginal curvature returned in neg_hess: Observed (the
    // exact -d^2 logM/d eta^2) or Expected (the PSD marginal Fisher information).
    static LatentMarginal marginal(const SiteData& y, double eta,
                                   double disp2, double sigma_u,
                                   const tulpa::GaussHermite& gh,
                                   tulpa::CurvatureMode curv
                                       = tulpa::CurvatureMode::Observed) {
        const double phi = disp2;
        const double b   = sigma_u * sigma_u;
        const double inv_b = 1.0 / b;

        const double NEG_INF = -std::numeric_limits<double>::infinity();

        // Newton mode-find on q(u) = ell(eta + u) + logN(u; 0, b). At extreme
        // outer-grid cells eta + u can drive the beta mean to 0/1, where
        // lgamma/digamma diverge; guard every evaluation so a degenerate cell
        // resolves to log_m = -Inf (zero outer-grid weight) instead of a NaN
        // that would poison the whole grid's weight vector.
        double u = 0.0, ell, s, sneg, sfish;
        for (int it = 0; it < 50; ++it) {
            ell_and_derivs(y, eta + u, phi, ell, s, sneg, sfish);
            if (!std::isfinite(s) || !std::isfinite(sneg)) break;
            const double qpp = -sneg - inv_b;          // q''(u) < 0
            if (qpp >= 0.0 || !std::isfinite(qpp)) break;
            const double step = (s - u * inv_b) / qpp;
            if (!std::isfinite(step)) break;
            u -= step;
            if (std::abs(step) < 1e-10) break;
        }
        double uhat = u;
        ell_and_derivs(y, eta + uhat, phi, ell, s, sneg, sfish);
        double tau2 = 1.0 / (sneg + inv_b);            // -1 / q''(uhat)
        if (!std::isfinite(tau2) || tau2 <= 0.0) { uhat = 0.0; tau2 = b; }
        const double tau = std::sqrt(tau2);

        const int   n = (int) gh.nodes.size();
        const double half_log_2pi = 0.5 * std::log(2.0 * M_PI);
        const double log_norm_const = 0.5 * std::log(2.0 * M_PI * b);

        // First pass: log weights for logsumexp. A node whose conditional
        // log-likelihood is non-finite (eta + u_k in the divergent tail) is
        // dropped to -Inf with zeroed score/info so 0 * Inf never produces NaN.
        std::vector<double> log_term(n), sk(n), snegk(n), sfishk(n);
        double max_lt = NEG_INF;
        for (int k = 0; k < n; ++k) {
            const double zk = gh.nodes[k];
            const double u_k = uhat + tau * zk;
            double ell_k, s_k, sneg_k, sfish_k;
            ell_and_derivs(y, eta + u_k, phi, ell_k, s_k, sneg_k, sfish_k);
            if (!std::isfinite(ell_k)) {
                log_term[k] = NEG_INF; sk[k] = 0.0; snegk[k] = 0.0; sfishk[k] = 0.0;
                continue;
            }
            const double log_h = ell_k - log_norm_const - 0.5 * u_k * u_k * inv_b;
            log_term[k] = std::log(gh.weights[k]) + 0.5 * zk * zk + log_h;
            sk[k] = s_k; snegk[k] = sneg_k; sfishk[k] = sfish_k;
            if (log_term[k] > max_lt) max_lt = log_term[k];
        }

        LatentMarginal out;
        if (!std::isfinite(max_lt)) {                  // every node degenerate
            out.log_m = NEG_INF; out.score = 0.0; out.neg_hess = 1e-8;
            return out;
        }

        double sum_w = 0.0;
        for (int k = 0; k < n; ++k) sum_w += std::exp(log_term[k] - max_lt);
        const double log_sum = max_lt + std::log(sum_w);

        // Posterior moments under pi_k = exp(log_term_k) / sum.
        double e_s = 0.0, e_s2 = 0.0, e_sneg = 0.0, e_sfish = 0.0;
        for (int k = 0; k < n; ++k) {
            const double pk = std::exp(log_term[k] - log_sum);
            e_s     += pk * sk[k];
            e_s2    += pk * sk[k] * sk[k];
            e_sneg  += pk * snegk[k];
            e_sfish += pk * sfishk[k];
        }
        const double var_s = e_s2 - e_s * e_s;

        out.log_m    = std::log(tau) + half_log_2pi + log_sum;
        out.score    = e_s;
        // Observed: E_pi[sneg] - Var_pi(s), the exact marginal information (can
        // be indefinite at extreme cells). Expected: E_pi[sum fisher], PSD by
        // construction, so the inner Newton step is well-conditioned there.
        out.neg_hess = (curv == tulpa::CurvatureMode::Expected)
                           ? e_sfish
                           : e_sneg - var_s;
        if (!std::isfinite(out.log_m)) {
            out.log_m = NEG_INF; out.score = 0.0; out.neg_hess = 1e-8;
        }
        return out;
    }
};

} // namespace tulpaObs

#endif // TULPAOBS_OCCU_COVER_LATENT_H
