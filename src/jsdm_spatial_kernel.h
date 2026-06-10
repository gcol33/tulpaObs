// jsdm_spatial_kernel.h
// Per-(species, site) Bernoulli cell for the areal-spatial joint species
// distribution model (jsdm() + shared field; tulpaObs#76). The JSDM has no
// detection process: the observed presence/absence enters the latent occupancy
// directly, so the per-(species, site) contribution is a plain Bernoulli on
//
//   logit psi_{s,i} = X_i . beta + b_s + f_{u(i)}
//
// with beta the shared fixed-effect coefficients, b_s ~ N(0, sigma_re^2) a
// scalar per-species random intercept, and f one shared ICAR / BYM2 / proper-CAR
// areal field. There is no latent state to integrate out (y is observed), so the
// cell is convex and the observed negative Hessian equals the Fisher information
// psi(1 - psi); no observed/complete-data branch is needed.

#ifndef TULPAOBS_JSDM_SPATIAL_KERNEL_H
#define TULPAOBS_JSDM_SPATIAL_KERNEL_H

#include <cmath>

namespace tulpaObs {

// Clamp a probability away from {0, 1}, matching the occupancy kernels.
inline double jsdm_clamp_(double x) {
    const double lo = 1e-12, hi = 1.0 - 1e-12;
    return x < lo ? lo : (x > hi ? hi : x);
}

inline double jsdm_sigmoid_(double eta) { return 1.0 / (1.0 + std::exp(-eta)); }

// One (species, site) Bernoulli cell in eta-space: log-likelihood, score
// g = d ll / d eta = y - psi, and the negative Hessian B = -d^2 ll / d eta^2 =
// psi(1 - psi) (exact: the Bernoulli log-likelihood is concave, observed info =
// Fisher info).
struct JsdmSiteCell {
    double log_lik = 0.0;
    double g = 0.0;     // d ll / d eta
    double B = 0.0;     // -d^2 ll / d eta^2
};

inline JsdmSiteCell jsdm_site_cell(double eta, int y) {
    JsdmSiteCell out;
    const double psi = jsdm_clamp_(jsdm_sigmoid_(eta));
    out.log_lik = (y == 1) ? std::log(psi) : std::log1p(-psi);
    out.g = (double) y - psi;
    out.B = psi * (1.0 - psi);
    return out;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_JSDM_SPATIAL_KERNEL_H
