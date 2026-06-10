// nuts_field_block.h
// Shared fixed-hyper non-centered areal field block for the count / occupancy
// observation-family NUTS targets (gcol33/tulpaObs#72). The field is a whitened
// Gaussian z = Linv %*% raw, raw ~ N(0, I), with Linv the inverse Cholesky of the
// FIXED field precision tau Q(rho) (a small ridge proper-ises an intrinsic ICAR).
// The field covariance is fixed at the nested-Laplace posterior estimate; NUTS
// samples only the whitened raw (and the family coefficients). The field loads
// additively onto one arm's per-site linear predictor,
//   eta_arm[site] += z[field_map[site]].
//
// This is the same forward / backward arithmetic the N-mixture spatial NUTS
// (marginal_count_nuts.h) carries inline; factoring it here lets the removal,
// distance, false-positive occupancy, and open N-mixture samplers reuse one
// implementation rather than each re-deriving the chain rule. A family's eval:
//   FieldBlock fb = field_block_build(spec, base);       // base = first free coord
//   std::vector<double> z, gz;
//   field_block_forward(fb, theta, z);                   // z[unit]; add z[map[s]]
//   ... accumulate gz[map[s]] += d log L / d eta_arm[s] ...
//   lp += field_block_backward(fb, theta, gz, grad);     // grad over raw + prior

#ifndef TULPAOBS_NUTS_FIELD_BLOCK_H
#define TULPAOBS_NUTS_FIELD_BLOCK_H

#include <Rcpp.h>
#include <vector>
#include <cstddef>

namespace tulpaObs {

// Marshalled fixed-hyper areal field. n_field_units == 0 => no field (every
// helper is then a no-op and the flat vector carries no raw coordinates).
struct FieldBlock {
    int n_field_units = 0;            // field nodes (== spatial units)
    int o_raw = 0;                    // first whitened-field coordinate in theta
    std::vector<int> field_map;       // 0-based unit per site (length n_sites)
    std::vector<double> Linv;         // row-major n_field_units x n_field_units

    bool active() const { return n_field_units > 0; }
};

// Parse the optional field block from a NUTS spec. `base` is the first free flat
// coordinate (after the coefficient + dispersion + RE blocks); the n_field_units
// whitened-field coordinates follow it. `n_sites` is the per-site count the
// field_map must cover. Returns an inactive block when the spec carries no field.
inline FieldBlock field_block_build(const Rcpp::List& spec, int base, int n_sites) {
    FieldBlock fb;
    if (!spec.containsElementNamed("n_field_units")) return fb;
    fb.n_field_units = Rcpp::as<int>(spec["n_field_units"]);
    if (fb.n_field_units <= 0) { fb.n_field_units = 0; return fb; }
    Rcpp::IntegerVector fm = spec["field_map"];        // 1-based site -> unit
    if ((int) fm.size() != n_sites)
        Rcpp::stop("field_map must have length n_sites");
    fb.field_map.resize(n_sites);
    for (int i = 0; i < n_sites; ++i) {
        const int u = fm[i] - 1;
        if (u < 0 || u >= fb.n_field_units)
            Rcpp::stop("field_map values must lie in [1, n_field_units]");
        fb.field_map[i] = u;
    }
    Rcpp::NumericMatrix Li = spec["field_Linv"];
    if (Li.nrow() != fb.n_field_units || Li.ncol() != fb.n_field_units)
        Rcpp::stop("field_Linv must be n_field_units x n_field_units");
    fb.Linv.resize((std::size_t) fb.n_field_units * fb.n_field_units);
    for (int u = 0; u < fb.n_field_units; ++u)
        for (int v = 0; v < fb.n_field_units; ++v)
            fb.Linv[(std::size_t) u * fb.n_field_units + v] = Li(u, v);
    fb.o_raw = base;
    return fb;
}

// Number of flat coordinates the field block contributes (0 when inactive).
inline int field_block_size(const FieldBlock& fb) { return fb.n_field_units; }

// Forward map z = Linv %*% raw into `z` (resized to n_field_units). The caller
// adds z[field_map[site]] to the chosen arm's eta. No-op (clears z) when inactive.
inline void field_block_forward(const FieldBlock& fb, const double* theta,
                                std::vector<double>& z) {
    if (!fb.active()) { z.clear(); return; }
    const int n = fb.n_field_units;
    z.assign(n, 0.0);
    for (int u = 0; u < n; ++u) {
        double zz = 0.0;
        const double* Lu = &fb.Linv[(std::size_t) u * n];
        for (int v = 0; v < n; ++v) zz += Lu[v] * theta[fb.o_raw + v];
        z[u] = zz;
    }
}

// Resize and zero a grad_z accumulator to n_field_units (no-op when inactive).
inline void field_block_init_grad(const FieldBlock& fb, std::vector<double>& grad_z) {
    if (!fb.active()) { grad_z.clear(); return; }
    grad_z.assign(fb.n_field_units, 0.0);
}

// Backward pass: with grad_z[u] = sum over sites mapped to u of d log L / d eta,
// add the chain-rule whitened-field gradient grad_raw = Linv^T grad_z and the
// raw ~ N(0, I) prior gradient into `grad`, and return the prior log-density
// contribution -0.5 ||raw||^2. No-op returning 0 when inactive.
inline double field_block_backward(const FieldBlock& fb, const double* theta,
                                   const std::vector<double>& grad_z, double* grad) {
    if (!fb.active()) return 0.0;
    const int n = fb.n_field_units;
    double lp = 0.0;
    for (int v = 0; v < n; ++v) {
        double gr = 0.0;
        for (int u = 0; u < n; ++u)
            gr += fb.Linv[(std::size_t) u * n + v] * grad_z[u];
        const double rv = theta[fb.o_raw + v];
        lp -= 0.5 * rv * rv;
        grad[fb.o_raw + v] += gr - rv;
    }
    return lp;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NUTS_FIELD_BLOCK_H
