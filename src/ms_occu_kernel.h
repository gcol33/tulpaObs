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
