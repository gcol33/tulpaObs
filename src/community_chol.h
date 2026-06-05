// community_chol.h
// Dense linear algebra on a community-covariance Cholesky factor, shared by the
// community NUTS targets (the spatial-factor occu_cover model,
// ms_occu_cover_spatial_nuts.cpp, and the community N-mixture, ms_abun_nuts.cpp).
//
// A P x P community covariance Sigma = C C' is carried by its lower-triangular
// Cholesky factor C, packed column-major over the lower triangle with the
// diagonal on the log scale (the unconstrained map onto a PD matrix; P is the
// number of coefficients on one arm, so it is tiny). These helpers unpack the
// packed vector, invert / square the factor, and map the gradient of
//   T = -0.5 tr(Sigma^{-1} M) - 0.5 S log|Sigma|
// (plus the coordinate hyperprior) back to the packed coordinates. They mirror
// the R helpers .ms_ocs_chol_unpack / .ms_ocs_chol_block_grad and are the single
// C++ source for both community NUTS targets.

#ifndef TULPAOBS_COMMUNITY_CHOL_H
#define TULPAOBS_COMMUNITY_CHOL_H

#include <vector>
#include <cmath>
#include <cstddef>

namespace tulpaObs {

// Packed column-major lower-triangle (diagonal on the log scale) -> lower
// Cholesky factor C (row-major), matching .ms_ocs_chol_unpack.
inline void chol_unpack_cpp(const double* vec, int P, std::vector<double>& C) {
    C.assign((std::size_t) P * P, 0.0);
    int pos = 0;
    for (int j = 0; j < P; ++j) {
        C[(std::size_t) j * P + j] = std::exp(vec[pos++]);
        for (int i = j + 1; i < P; ++i) C[(std::size_t) i * P + j] = vec[pos++];
    }
}

// Inverse of a lower-triangular matrix (row-major), by forward substitution.
inline void lower_tri_inv(const std::vector<double>& C, int P,
                          std::vector<double>& M) {
    M.assign((std::size_t) P * P, 0.0);
    for (int j = 0; j < P; ++j) {
        M[(std::size_t) j * P + j] = 1.0 / C[(std::size_t) j * P + j];
        for (int i = j + 1; i < P; ++i) {
            double s = 0.0;
            for (int k = j; k < i; ++k)
                s += C[(std::size_t) i * P + k] * M[(std::size_t) k * P + j];
            M[(std::size_t) i * P + j] = -s / C[(std::size_t) i * P + i];
        }
    }
}

// Sigma^{-1} = Cinv' Cinv (Cinv lower-tri, row-major) -> row-major P*P.
inline void sinv_from_cinv(const std::vector<double>& Mi, int P,
                           std::vector<double>& Si) {
    Si.assign((std::size_t) P * P, 0.0);
    for (int a = 0; a < P; ++a)
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int k = 0; k < P; ++k)
                s += Mi[(std::size_t) k * P + a] * Mi[(std::size_t) k * P + b];
            Si[(std::size_t) a * P + b] = s;
        }
}

// Cholesky-coordinate gradient of T = -0.5 tr(Sigma^{-1} M) - 0.5 S log|Sigma|
// plus the coordinate hyperprior, written into `out` (packed column-major lower
// triangle). Mirrors .ms_ocs_chol_block_grad: G = 0.5 Si M Si - 0.5 S Si,
// dC = 2 G C, with the log-link chain on the diagonal.
inline void chol_block_grad_cpp(const std::vector<double>& C,
                                const std::vector<double>& Si,
                                const std::vector<double>& M, int P, double S,
                                const double* vec, double logdiag_mean,
                                double logdiag_sd, double offdiag_sd,
                                double* out) {
    std::vector<double> SM((std::size_t) P * P, 0.0), G((std::size_t) P * P, 0.0),
                        dC((std::size_t) P * P, 0.0);
    for (int a = 0; a < P; ++a)
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int k = 0; k < P; ++k) s += Si[a * P + k] * M[k * P + b];
            SM[a * P + b] = s;
        }
    for (int a = 0; a < P; ++a)
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int k = 0; k < P; ++k) s += SM[a * P + k] * Si[k * P + b];
            G[a * P + b] = 0.5 * s - 0.5 * S * Si[a * P + b];
        }
    for (int a = 0; a < P; ++a)
        for (int b = 0; b < P; ++b) {
            double s = 0.0;
            for (int k = 0; k < P; ++k) s += G[a * P + k] * C[k * P + b];
            dC[a * P + b] = 2.0 * s;
        }
    int pos = 0;
    for (int j = 0; j < P; ++j) {
        const double cjj = C[(std::size_t) j * P + j];
        out[pos] = dC[(std::size_t) j * P + j] * cjj
                 - (vec[pos] - logdiag_mean) / (logdiag_sd * logdiag_sd);
        ++pos;
        for (int i = j + 1; i < P; ++i) {
            out[pos] = dC[(std::size_t) i * P + j]
                     - vec[pos] / (offdiag_sd * offdiag_sd);
            ++pos;
        }
    }
}

}  // namespace tulpaObs

#endif  // TULPAOBS_COMMUNITY_CHOL_H
