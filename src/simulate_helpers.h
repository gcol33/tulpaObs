// simulate_helpers.h
// Pieces shared by the family generators (simulate_count.cpp,
// simulate_distance.cpp, simulate_single.cpp, simulate_ms.cpp). Each generator
// draws ONE replicate at a given parameter vector theta; which posterior draw
// theta is, is chosen on the R side.

#ifndef TULPAOBS_SIMULATE_HELPERS_H
#define TULPAOBS_SIMULATE_HELPERS_H

#include <Rcpp.h>
#include "tobs_math.h"

namespace tulpaObs {

// Latent abundance: negative binomial in the (size, mu) parameterisation when
// `r_size` is finite, Poisson otherwise. R's rnbinom(size, mu) IS
// rpois(rgamma(size, mu / size)) internally, so drawing it through the exposed
// R:: samplers consumes the same two-draw RNG stream and is byte-identical to
// the R-side simulator.
inline int draw_latent_N(double lambda, double r_size) {
    if (R_finite(r_size)) {
        if (lambda <= 0.0) return 0;
        return (int) R::rpois(R::rgamma(r_size, lambda / r_size));
    }
    return (int) R::rpois(lambda);
}

// Linear predictor of row i of an arm design X (nrow rows) at the coefficient
// block theta[off .. off + p - 1], plus that row's realised offset (a field or
// random effect) when the arm carries one, i.e. when `offs` is non-empty.
inline double arm_eta(const double* X, int nrow, int i, const double* theta,
                      int off, int p, const Rcpp::NumericVector& offs) {
    double a = row_draw_dot(X, nrow, i, theta, 1, 0, off, p);
    return offs.size() ? a + offs[i] : a;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_SIMULATE_HELPERS_H
