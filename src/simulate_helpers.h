// simulate_helpers.h
// Latent-state draws shared by the posterior-predictive simulators
// (simulate_count.cpp, simulate_distance.cpp, simulate_ms.cpp).

#ifndef TULPAOBS_SIMULATE_HELPERS_H
#define TULPAOBS_SIMULATE_HELPERS_H

#include <Rcpp.h>

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

}  // namespace tulpaObs

#endif  // TULPAOBS_SIMULATE_HELPERS_H
