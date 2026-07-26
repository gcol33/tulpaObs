// nuts_re_block.h
// Shared single-grouping intercept random-effect block for the observation-family
// NUTS targets (gcol33/tulpaObs#51, #82). The effect is non-centered: one whitened
// coordinate z_g ~ N(0, 1) per group plus one hyperparameter log_sigma_re, and the
// per-site offset
//   eta_arm[site] += sigma_re * z[group[site]],   sigma_re = exp(log_sigma_re)
// loads additively onto ONE arm's linear predictor.
//
// This is the arithmetic the N-mixture / removal (marginal_count_nuts.h), open
// N-mixture (dyn_abun_nuts.cpp), distance (distance_nuts.cpp), and false-positive
// occupancy (fp_occu_nuts.cpp) samplers each carried inline. It is the companion
// of nuts_field_block.h and follows the same build / forward / backward shape:
//   ReBlock rb = re_block_build(spec, base, n_sites);   // base = first free coord
//   const double s_re = re_block_sigma(rb, theta);
//   ... per site: eta_arm += re_block_offset(rb, s_re, theta, site);
//                 re_block_accumulate(rb, s_re, d log L / d eta_arm, site,
//                                     theta, grad, grad_logsig);
//   lp += re_block_backward(rb, theta, grad_logsig, grad);

#ifndef TULPAOBS_NUTS_RE_BLOCK_H
#define TULPAOBS_NUTS_RE_BLOCK_H

#include <Rcpp.h>
#include <cmath>
#include <vector>

namespace tulpaObs {

// Marshalled random-effect block. arm < 0 => no RE (every helper is then a no-op
// and the flat vector carries no z / log_sigma coordinates).
struct ReBlock {
    int arm = -1;               // -1 none; 0 = state arm, 1 = detection arm
    int n_groups = 0;           // grouping-factor levels
    int o_z = 0, o_logsig = 0;  // offsets of the z block / log_sigma_re in theta
    double sigma_lsd = 1.5;     // prior SD on log_sigma_re
    std::vector<int> group;     // 0-based group per site (length n_sites)

    bool active() const { return arm >= 0; }
};

// Parse the optional RE block. `base` is the first free flat coordinate; the
// n_groups whitened coordinates and the trailing log_sigma_re follow it.
// `max_arm` is the highest arm index the caller supports (0 = state arm only);
// an out-of-range arm is rejected rather than silently ignored.
inline ReBlock re_block_build(const Rcpp::List& spec, int base, int n_sites,
                              int max_arm = 0) {
    ReBlock rb;
    if (!spec.containsElementNamed("re_arm")) return rb;
    const int arm = Rcpp::as<int>(spec["re_arm"]);
    if (arm < 0) return rb;
    if (arm > max_arm)
        Rcpp::stop("re_arm = %d is not supported by this family (max %d)",
                   arm, max_arm);
    rb.arm = arm;
    Rcpp::IntegerVector rg = spec["re_group"];          // 1-based, length n_sites
    if ((int) rg.size() != n_sites)
        Rcpp::stop("re_group must have length n_sites");
    rb.group.resize(n_sites);
    for (int i = 0; i < n_sites; ++i) rb.group[i] = rg[i] - 1;
    rb.n_groups = Rcpp::as<int>(spec["n_re_groups"]);
    if (spec.containsElementNamed("sigma_re_lsd"))
        rb.sigma_lsd = Rcpp::as<double>(spec["sigma_re_lsd"]);
    rb.o_z = base;
    rb.o_logsig = base + rb.n_groups;
    return rb;
}

// Flat coordinates the block occupies (0 when inactive).
inline int re_block_size(const ReBlock& rb) {
    return rb.active() ? rb.n_groups + 1 : 0;
}

// sigma_re = exp(log_sigma_re); 0 when inactive, which makes every offset 0.
inline double re_block_sigma(const ReBlock& rb, const double* theta) {
    return rb.active() ? std::exp(theta[rb.o_logsig]) : 0.0;
}

// Per-site offset added to the loaded arm's linear predictor.
inline double re_block_offset(const ReBlock& rb, double sigma,
                              const double* theta, int site) {
    return rb.active() ? sigma * theta[rb.o_z + rb.group[site]] : 0.0;
}

// Chain d log L / d eta_arm[site] into the z gradient, and accumulate the
// log_sigma score (d eta / d log_sigma = sigma_re * z = the offset itself).
inline void re_block_accumulate(const ReBlock& rb, double sigma, double grad_eta,
                                int site, const double* theta, double* grad,
                                double& grad_logsig) {
    if (!rb.active()) return;
    const int g = rb.o_z + rb.group[site];
    grad[g] += sigma * grad_eta;
    grad_logsig += grad_eta * (sigma * theta[g]);
}

// Whitened prior z ~ N(0, I) + log_sigma_re ~ N(0, sigma_lsd^2). Adds the
// accumulated data score for log_sigma_re. Returns the log-prior contribution.
inline double re_block_backward(const ReBlock& rb, const double* theta,
                                double grad_logsig, double* grad) {
    if (!rb.active()) return 0.0;
    double lp = 0.0;
    for (int g = 0; g < rb.n_groups; ++g) {
        const double zg = theta[rb.o_z + g];
        lp -= 0.5 * zg * zg;
        grad[rb.o_z + g] -= zg;
    }
    const double ls = theta[rb.o_logsig];
    const double ils2 = 1.0 / (rb.sigma_lsd * rb.sigma_lsd);
    lp -= 0.5 * ils2 * ls * ls;
    grad[rb.o_logsig] += grad_logsig - ils2 * ls;
    return lp;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NUTS_RE_BLOCK_H
