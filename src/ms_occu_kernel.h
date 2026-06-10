// ms_occu_kernel.h
// Per-species occupancy two-state marginal for the community single-season
// occupancy NUTS target (ms_occu(), method = "nuts"; tulpaObs#69). Site-level
// detection: one detection probability per (species, site), constant across that
// site's visits, summarised by (n_valid, n_det) per site. The latent presence z
// is integrated out exactly per site:
//
//   psi_i = sigmoid(eta_psi_i),  p_i = sigmoid(eta_p_i)
//   any-detection site i (n_det_i > 0):
//     log L_i = log psi_i + n_det_i log p_i + (n_valid_i - n_det_i) log(1 - p_i)
//   no-detection site i:
//     L_i = psi_i (1 - p_i)^{n_valid_i} + (1 - psi_i),   log L_i = log L_i
//
// The marginal log-likelihood is sum_i log L_i; the per-site scores wrt eta_psi
// and eta_p are returned for the design-sandwiched coefficient gradient. This is
// the single-source detection case of the integrated-community marginal and
// mirrors the R oracle .ms_int_occu_sp_ll / .ms_int_occu_sp_grad byte-for-byte.

#ifndef TULPAOBS_MS_OCCU_KERNEL_H
#define TULPAOBS_MS_OCCU_KERNEL_H

#include <vector>
#include <cmath>
#include <cstddef>

namespace tulpaObs {

// Clamp a probability away from the {0, 1} boundary, matching the R kernel's
// pmin(pmax(x, 1e-12), 1 - 1e-12).
inline double msocc_clamp_(double x) {
    const double lo = 1e-12, hi = 1.0 - 1e-12;
    return x < lo ? lo : (x > hi ? hi : x);
}

inline double msocc_sigmoid_(double eta) { return 1.0 / (1.0 + std::exp(-eta)); }

// Per-(species) site-level summaries: n_valid_i visits, n_det_i detections.
struct MsOccuSiteSummary {
    std::vector<int> n_valid;   // length n_sites
    std::vector<int> n_det;     // length n_sites
    std::vector<char> any_det;  // length n_sites (n_det_i > 0)
};

// Result of one species' occupancy marginal sweep: the marginal log-likelihood
// and the per-site eta scores (length n_sites each).
struct MsOccuSiteResult {
    double log_lik = 0.0;
    std::vector<double> grad_eta_psi;
    std::vector<double> grad_eta_p;
};

// One site's occupancy two-state marginal contribution in eta-space: the log
// likelihood, the scores g = (d ll / d eta_psi, d ll / d eta_p), and the 2x2
// negative-Hessian B = -d^2 ll / d eta d eta'. The areal nested-Laplace driver
// sandwiches B through the per-site design (psi loads the field, p the detection
// design). `observed` selects the true observed Hessian (final marginal pass);
// when false the responsibility-weighted complete-data Fisher block-diagonal
// curvature is returned (the EM mode-find, always PSD).
struct MsOccuSiteCell {
    double log_lik = 0.0;
    double g_psi = 0.0, g_p = 0.0;
    double B_pp_psi = 0.0;   // -d^2 / d eta_psi^2
    double B_p_p   = 0.0;    // -d^2 / d eta_p^2
    double B_cross = 0.0;    // -d^2 / d eta_psi d eta_p
};

// Evaluate one site's occupancy two-state marginal cell (log-lik, score,
// curvature) at (eta_psi, eta_p) given the site summary (n_valid, n_det). A
// no-visit site (n_valid == 0) contributes nothing.
inline MsOccuSiteCell ms_occu_site_cell(double eta_psi, double eta_p,
                                        int n_valid, int n_det, bool any_det,
                                        bool observed) {
    MsOccuSiteCell out;
    if (n_valid <= 0) return out;
    const double psi = msocc_clamp_(msocc_sigmoid_(eta_psi));
    const double p   = msocc_clamp_(msocc_sigmoid_(eta_p));
    const double w1mw = psi * (1.0 - psi);
    const double p1mp = p * (1.0 - p);
    if (any_det) {
        // z = 1 forced; psi and p factorise.
        out.log_lik = std::log(psi) + n_det * std::log(p)
                    + (n_valid - n_det) * std::log1p(-p);
        out.g_psi = 1.0 - psi;
        out.g_p   = (double) n_det - (double) n_valid * p;
        // -d^2: psi(1-psi) on the occupancy arm; n_valid p(1-p) on detection.
        out.B_pp_psi = w1mw;
        out.B_p_p    = (double) n_valid * p1mp;
        out.B_cross  = 0.0;
        return out;
    }
    // No-detection: L = psi q + (1 - psi), q = (1 - p)^n_valid.
    const double q = std::exp((double) n_valid * std::log1p(-p));
    const double A = psi * q;            // occupied-undetected mass
    const double L = A + (1.0 - psi);
    const double invL = (L > 0.0) ? 1.0 / L : 0.0;
    out.log_lik = std::log(L);
    out.g_psi = w1mw * (q - 1.0) * invL;
    // dq/d eta_p = -n_valid p q.
    const double dq = -(double) n_valid * p * q;
    out.g_p = psi * dq * invL;
    if (!observed) {
        // Complete-data (responsibility-weighted) Fisher, block-diagonal PSD:
        //   gamma = P(z = 1 | y, undetected) = A / L,
        //   occupancy info = psi(1-psi); detection info = gamma n_valid p(1-p).
        const double gamma = A * invL;
        out.B_pp_psi = w1mw;
        out.B_p_p    = gamma * (double) n_valid * p1mp;
        out.B_cross  = 0.0;
        return out;
    }
    // Observed negative Hessian of log L. Let f_psi = dL/d eta_psi = w1mw (q - 1),
    // f_p = dL/d eta_p = psi dq.
    const double f_psi = w1mw * (q - 1.0);
    const double f_p   = psi * dq;
    // d^2 L / d eta_psi^2 = w1mw (1 - 2 psi) (q - 1).
    const double L_psipsi = w1mw * (1.0 - 2.0 * psi) * (q - 1.0);
    // d^2 q / d eta_p^2 = n_valid p q [ n_valid p - (1 - p) ] / ... derive:
    //   dq = -n_valid p q;  d(dq)/d eta_p = -n_valid q [ p1mp + p * (dq/q) ]
    //      = -n_valid q [ p(1-p) - n_valid p^2 ].
    const double d2q = -(double) n_valid * q * (p1mp - (double) n_valid * p * p);
    const double L_pp = psi * d2q;
    // d^2 L / d eta_psi d eta_p = w1mw dq.
    const double L_psip = w1mw * dq;
    // -d^2 logL = -[ L'' / L - (L'/L)(L'/L) ].
    out.B_pp_psi = -(L_psipsi * invL - f_psi * f_psi * invL * invL);
    out.B_p_p    = -(L_pp   * invL - f_p   * f_p   * invL * invL);
    out.B_cross  = -(L_psip * invL - f_psi * f_p   * invL * invL);
    return out;
}

// Evaluate one species' occupancy two-state marginal over all sites. `eta_psi`
// and `eta_p` are the per-site occupancy / detection linear predictors (length
// n_sites). Writes the per-site scores into `res` (resized as needed).
inline void compute_ms_occu_site(const double* eta_psi, const double* eta_p,
                                 const MsOccuSiteSummary& summ, int n_sites,
                                 bool want_grad, MsOccuSiteResult& res) {
    res.log_lik = 0.0;
    if (want_grad) {
        res.grad_eta_psi.assign((std::size_t) n_sites, 0.0);
        res.grad_eta_p.assign((std::size_t) n_sites, 0.0);
    }
    for (int i = 0; i < n_sites; ++i) {
        const double psi = msocc_clamp_(msocc_sigmoid_(eta_psi[i]));
        const double p   = msocc_clamp_(msocc_sigmoid_(eta_p[i]));
        const int nv = summ.n_valid[i];
        const int nd = summ.n_det[i];
        if (summ.any_det[i]) {
            // z = 1 forced.
            res.log_lik += std::log(psi) + nd * std::log(p)
                         + (nv - nd) * std::log1p(-p);
            if (want_grad) {
                res.grad_eta_psi[i] = 1.0 - psi;                 // d/d eta_psi
                res.grad_eta_p[i]   = (double) nd - (double) nv * p;
            }
        } else {
            // L = psi (1 - p)^nv + (1 - psi).
            const double prodterm = std::exp((double) nv * std::log1p(-p));
            const double A = psi * prodterm;       // occupied-undetected mass
            const double B = 1.0 - psi;            // unoccupied mass
            const double L = A + B;
            res.log_lik += std::log(L);
            if (want_grad) {
                res.grad_eta_psi[i] = psi * (1.0 - psi) * (prodterm - 1.0) / L;
                res.grad_eta_p[i]   = -(A / L) * (double) nv * p;
            }
        }
    }
}

}  // namespace tulpaObs

#endif  // TULPAOBS_MS_OCCU_KERNEL_H
