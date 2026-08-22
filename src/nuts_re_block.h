// nuts_re_block.h
// Shared single-grouping intercept random-effect block for the
// observation-family NUTS targets. The effect is non-centered: one whitened
// coordinate z_g ~ N(0, 1) per group plus one hyperparameter log_sigma_re, and
// the per-row offset
//   eta_arm[row] += sigma_re * z[group[row]],   sigma_re = exp(log_sigma_re)
// loads additively onto ONE arm's linear predictor. The GROUP SD IS SAMPLED: the
// block owns the coordinate log_sigma_re under a N(0, sigma_lsd^2) prior, so the
// posterior integrates the variance component rather than conditioning on a point
// estimate of it.
//
// A "row" is whatever the caller's arm indexes: the count / occupancy families
// carry one code per SITE, while the joint occupancy + cover target
// (occu_cover_nuts.cpp) carries one per (site, visit) observation row. A code of
// 0 marks a row the effect does not load on (a padded visit, or an unseen level),
// matching the deterministic engine's 0 scatter sentinel; it contributes no
// offset and no score.
//
// This is the arithmetic the N-mixture / removal (marginal_count_nuts.h), open
// N-mixture (dyn_abun_nuts.cpp), distance (distance_nuts.cpp), and false-positive
// occupancy (fp_occu_nuts.cpp) samplers each carried inline. It is the companion
// of nuts_field_block.h and follows the same build / forward / backward shape:
//   ReBlock rb = re_block_build(spec, base, n_rows);   // base = first free coord
//   const double s_re = re_block_sigma(rb, theta);
//   ... per row:  eta_arm += re_block_offset(rb, s_re, theta, row);
//                 re_block_accumulate(rb, s_re, d log L / d eta_arm, row,
//                                     theta, grad, grad_logsig);
//   lp += re_block_backward(rb, theta, grad_logsig, grad);
//
// re_block_build_list() marshals SEVERAL such blocks from one list entry (crossed
// / nested groupings, or one block per arm), laying them out back to back.

#ifndef TULPAOBS_NUTS_RE_BLOCK_H
#define TULPAOBS_NUTS_RE_BLOCK_H

#include <Rcpp.h>
#include <cmath>
#include <cstddef>
#include <vector>

namespace tulpaObs {

// Marshalled random-effect block. arm < 0 => no RE (every helper is then a no-op
// and the flat vector carries no z / log_sigma coordinates).
struct ReBlock {
    int arm = -1;               // -1 none; the caller's own arm index
    int n_groups = 0;           // grouping-factor levels
    int o_z = 0, o_logsig = 0;  // offsets of the z block / log_sigma_re in theta
    double sigma_lsd = 1.5;     // prior SD on log_sigma_re
    std::vector<int> group;     // 0-based group per design row; -1 = no effect

    bool active() const { return arm >= 0; }
};

// Parse the optional RE block. `base` is the first free flat coordinate; the
// n_groups whitened coordinates and the trailing log_sigma_re follow it.
// `n_rows` is the length of the arm's design (one group code per row).
// `max_arm` is the highest arm index the caller supports (0 = state arm only);
// an out-of-range arm is rejected rather than silently ignored.
inline ReBlock re_block_build(const Rcpp::List& spec, int base, int n_rows,
                              int max_arm = 0) {
    ReBlock rb;
    if (!spec.containsElementNamed("re_arm")) return rb;
    const int arm = Rcpp::as<int>(spec["re_arm"]);
    if (arm < 0) return rb;
    if (arm > max_arm)
        Rcpp::stop("re_arm = %d is not supported by this family (max %d)",
                   arm, max_arm);
    rb.arm = arm;
    rb.n_groups = Rcpp::as<int>(spec["n_re_groups"]);
    Rcpp::IntegerVector rg = spec["re_group"];      // 1-based; 0 = no effect
    if ((int) rg.size() != n_rows)
        Rcpp::stop("re_group must carry one code per design row (%d), got %d",
                   n_rows, (int) rg.size());
    rb.group.resize(n_rows);
    for (int i = 0; i < n_rows; ++i) {
        const int g = rg[i];
        if (g < 0 || g > rb.n_groups)
            Rcpp::stop("re_group codes must lie in [0, n_re_groups = %d]",
                       rb.n_groups);
        rb.group[i] = g - 1;
    }
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

// Marshal a LIST of RE block descriptions (crossed / nested groupings, or one
// block per arm) laid out back to back from `base`. Each element is read by
// re_block_build, so the blocks share one parameterisation and one chain rule.
// An absent / empty entry yields no blocks and no coordinates.
inline std::vector<ReBlock> re_block_build_list(const Rcpp::List& spec,
                                                const char* key, int base,
                                                int n_rows, int max_arm) {
    std::vector<ReBlock> out;
    if (!spec.containsElementNamed(key)) return out;
    Rcpp::List blocks = Rcpp::as<Rcpp::List>(spec[key]);
    out.reserve(blocks.size());
    for (int b = 0; b < blocks.size(); ++b) {
        ReBlock rb = re_block_build(Rcpp::as<Rcpp::List>(blocks[b]), base,
                                    n_rows, max_arm);
        base += re_block_size(rb);
        out.push_back(rb);
    }
    return out;
}

// Total flat coordinates a list of blocks occupies.
inline int re_block_list_size(const std::vector<ReBlock>& blocks) {
    int n = 0;
    for (std::size_t b = 0; b < blocks.size(); ++b) n += re_block_size(blocks[b]);
    return n;
}

// sigma_re = exp(log_sigma_re); 0 when inactive, which makes every offset 0.
inline double re_block_sigma(const ReBlock& rb, const double* theta) {
    return rb.active() ? std::exp(theta[rb.o_logsig]) : 0.0;
}

// Per-row offset added to the loaded arm's linear predictor. A row whose group
// code was 0 (padded visit / unseen level) carries no effect.
inline double re_block_offset(const ReBlock& rb, double sigma,
                              const double* theta, int row) {
    if (!rb.active()) return 0.0;
    const int g = rb.group[row];
    return g < 0 ? 0.0 : sigma * theta[rb.o_z + g];
}

// Chain d log L / d eta_arm[row] into the z gradient, and accumulate the
// log_sigma score (d eta / d log_sigma = sigma_re * z = the offset itself).
inline void re_block_accumulate(const ReBlock& rb, double sigma, double grad_eta,
                                int row, const double* theta, double* grad,
                                double& grad_logsig) {
    if (!rb.active()) return;
    const int g0 = rb.group[row];
    if (g0 < 0) return;
    const int g = rb.o_z + g0;
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
