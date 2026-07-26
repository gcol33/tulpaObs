// tobs_math.h
// Scalar numerical primitives shared across the tulpaObs kernels.
//
// Every one of these was independently reinvented in a handful of translation
// units before it lived here; they are small enough to inline and generic
// enough that no family owns them. Nothing in this header knows about a model
// type, an arm layout, or a likelihood.

#ifndef TULPAOBS_TOBS_MATH_H
#define TULPAOBS_TOBS_MATH_H

#include <cmath>
#include <cstddef>

namespace tulpaObs {

// Linear-predictor clamp bound. A logit clamped to [-30, 30] gives detection /
// occupancy probabilities within ~9.4e-14 of {0, 1}, far inside double
// precision, while keeping exp() away from overflow. The same bound guards
// log-mean predictors on the count arms before exp(). The R twin is
// `.TOBS_ETA_BOUND` in R/utils_numeric.R.
constexpr double kEtaClampBound = 30.0;

inline double clamp_eta(double e, double bound = kEtaClampBound) {
    return e > bound ? bound : (e < -bound ? -bound : e);
}

// Logistic function, branch-split on the sign of `x` so the exponential is
// always taken of a non-positive argument. This is the form R's plogis() uses,
// so a diagnostic / pointwise-likelihood kernel written against an R oracle
// reproduces it to the bit. The fit kernels take the plain 1/(1+exp(-eta))
// instead (occu_coupling_shared.h::sigmoid_): the two agree to within one ULP
// and neither overflows, but they do not round identically for negative eta, so
// the R-matching path keeps this one.
inline double stable_plogis(double x) {
    if (x >= 0.0) { const double z = std::exp(-x); return 1.0 / (1.0 + z); }
    const double z = std::exp(x); return z / (1.0 + z);
}

// log(x) guarded at x <= 0, returning -1e300 rather than -Inf / NaN, so a
// density evaluated at the boundary of its support stays finite and comparable.
// The R twin is `.tobs_log_safe()`.
inline double log_safe(double x) {
    return (x > 0.0) ? std::log(x) : -1e300;
}

// log(exp(a) + exp(b)), max-shifted. The shifted term is <= 0, so the other
// contributes exp(0) = 1 exactly and the sum is log1p(exp(min - max)) with no
// catastrophic cancellation. Two infinite equal arguments return that infinity
// rather than NaN.
inline double logsumexp2(double a, double b) {
    const double m = a > b ? a : b;
    if (std::isinf(m) && a == b) return m;
    return m + std::log1p(std::exp((a > b ? b : a) - m));
}

// Posterior-predictive-check discrepancy between an observed count `o` and its
// expectation `e`: the Freeman-Tukey statistic (sqrt(o) - sqrt(e))^2, or the
// chi-squared (o - e)^2 / e with a floor on the denominator so an expectation
// of zero does not divide by zero.
inline double ppc_stat(double o, double e, bool freeman) {
    if (freeman) { const double t = std::sqrt(o) - std::sqrt(e); return t * t; }
    const double t = o - e; return t * t / (e + 1e-10);
}

// Row `i` of a column-major [nrow x p] design dotted with the coefficient slice
// `dr[idx, off + 0 .. p-1]` of a column-major [ndr x .] draw matrix. The design
// and the draws are both R matrices read through their raw pointers, so this is
// the one place the column-major index arithmetic is written out.
inline double row_draw_dot(const double* X, int nrow, int i, const double* dr,
                           int ndr, int idx, int off, int p) {
    double a = 0.0;
    for (int k = 0; k < p; ++k)
        a += X[(std::size_t) k * nrow + i] *
             dr[(std::size_t) (off + k) * ndr + idx];
    return a;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_TOBS_MATH_H
