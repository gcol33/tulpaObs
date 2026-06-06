// nmix_oracle_emit.h
// Flat copy-out of a per-group Eval into the REGroupOracle output pointers.
// Every native N-mixture oracle (NMixCommunityOracle, NMixGroupedOracle,
// SpatialNMixCommunityOracle) assembles its per-group quantities in one site
// loop returning an Eval struct exposing .logL, .grad (length d), .negH and
// .fisher (d x d). These templates are the single source for the grad_hess /
// newton_hess copy-out; templated on Eval so each call inlines at the call site
// (no virtual dispatch, no per-oracle duplicate of the column-major fill).

#ifndef TULPAOBS_NMIX_ORACLE_EMIT_H
#define TULPAOBS_NMIX_ORACLE_EMIT_H

#include <cstddef>

namespace tulpaObs {

// REGroupOracle::grad_hess output: value, b-space score, marginal observed info.
template <class Eval>
inline void emit_grad_hess(const Eval& e, int d, double& logL,
                           double* grad, double* negH) {
    logL = e.logL;
    for (int i = 0; i < d; ++i) grad[i] = e.grad(i);
    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j) negH[(std::size_t)i * d + j] = e.negH(i, j);
}

// REGroupOracle::newton_hess output: the PSD complete-data Fisher curvature.
template <class Eval>
inline void emit_fisher(const Eval& e, int d, double* H) {
    for (int i = 0; i < d; ++i)
        for (int j = 0; j < d; ++j) H[(std::size_t)i * d + j] = e.fisher(i, j);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_ORACLE_EMIT_H
