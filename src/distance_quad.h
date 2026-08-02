// distance_quad.h
// Distance-sampling detection key functions and the per-bin Gauss-Legendre
// quadrature the distance kernel (distance_kernel.h) integrates them against.
//
// For a covered region of half-width W (line transect) or radius W (point
// transect), an individual's perpendicular / radial distance x has density
//   f(x) = 1 / W            (line:  uniform in x)
//   f(x) = 2 x / W^2        (point: area element 2 pi x, normalised on [0, W])
// and is detected with probability g(x; sigma[, b]) given by a key function. The
// probability an individual falls in distance bin b = [c_{b-1}, c_b] AND is
// detected is the integral
//   pi_b = integral_{c_{b-1}}^{c_b} g(x; theta) f(x) dx,
// and p = sum_b pi_b is the average detection probability over the region. These
// integrals (and their first and second derivatives in the log-scale detection
// parameters) have no common closed form across keys (hazard-rate has none), so
// they are evaluated by a fixed high-order Gauss-Legendre rule per bin: the
// detection functions are smooth, so the rule is effectively exact. The x nodes
// and the f(x)-folded base weights depend only on the bins and the transect
// geometry, so they are built once per fit; only g(x; theta) is re-evaluated as
// the detection parameters change.
//
// Key functions (theta on the log scale: eta_sigma = log sigma, eta_b = log b):
//   half-normal:  g(x) = exp(-x^2 / (2 sigma^2))
//   hazard-rate:  g(x) = 1 - exp(-(x / sigma)^(-b)),   b > 0 shape
// Derivatives are taken w.r.t. eta_sigma (and eta_b for hazard-rate), so the
// detection arm enters the marginal exactly like a log-linear predictor; see the
// derivation in distance_kernel.h.
//
// References:
//   Buckland, Anderson, Burnham, Laake, Borchers, Thomas (2001) Introduction to
//     Distance Sampling. Oxford University Press.
//   Royle, Dawson & Bates (2004) Ecology 85: 1591-1597 (distance + N-mixture).

#ifndef TULPAOBS_DISTANCE_QUAD_H
#define TULPAOBS_DISTANCE_QUAD_H

#include <Rcpp.h>
#include <cmath>
#include <vector>

namespace tulpaObs {

enum DistKey { DIST_HALFNORMAL = 0, DIST_HAZARD = 1 };
enum DistTransect { DIST_LINE = 0, DIST_POINT = 1 };

// Gauss-Legendre nodes/weights on [-1, 1] via Newton-Raphson on the Legendre
// polynomial P_n (standard three-term recurrence for P_n and P_n'). Generated at
// runtime so the rule order is free and there are no hand-transcribed tables.
inline void gauss_legendre(int n, std::vector<double>& nodes,
                           std::vector<double>& weights) {
    nodes.assign(n, 0.0);
    weights.assign(n, 0.0);
    const double PI = 3.14159265358979323846;
    const int m = (n + 1) / 2;                 // roots are symmetric about 0
    for (int i = 0; i < m; ++i) {
        double x = std::cos(PI * (i + 0.75) / (n + 0.5));   // initial guess
        double dp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = x;            // P_0, P_1
            for (int k = 2; k <= n; ++k) {
                double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
                p0 = p1; p1 = p2;
            }
            dp = n * (x * p1 - p0) / (x * x - 1.0);          // P_n'(x)
            double dx = p1 / dp;
            x -= dx;
            if (std::fabs(dx) < 1e-15) break;
        }
        const double w = 2.0 / ((1.0 - x * x) * dp * dp);
        nodes[i]         = -x;
        nodes[n - 1 - i] =  x;
        weights[i]         = w;
        weights[n - 1 - i] = w;
    }
}

// Per-fit quadrature: for each of n_bins distance bins, the mapped x nodes and
// the base weights (Gauss weight x bin half-width x f(x)) so that
//   pi_b(theta) = sum_q base_w[b][q] * g(x[b][q]; theta).
struct DistQuad {
    int n_bins = 0, order = 0, transect = DIST_LINE;
    double W = 0.0;
    std::vector<double> cut;                 // length n_bins + 1
    std::vector<std::vector<double>> x;      // [n_bins][order] node positions
    std::vector<std::vector<double>> base_w; // [n_bins][order] f-folded weights
};

inline DistQuad dist_build_quad(const std::vector<double>& cut, int transect,
                                int order = 64) {
    DistQuad q;
    q.n_bins   = (int)cut.size() - 1;
    q.order    = order;
    q.transect = transect;
    q.cut      = cut;
    q.W        = cut.back();
    std::vector<double> gnode, gw;
    gauss_legendre(order, gnode, gw);
    q.x.assign(q.n_bins, std::vector<double>(order));
    q.base_w.assign(q.n_bins, std::vector<double>(order));
    const double W = q.W;
    for (int b = 0; b < q.n_bins; ++b) {
        const double lo = cut[b], hi = cut[b + 1];
        const double half = 0.5 * (hi - lo), mid = 0.5 * (hi + lo);
        for (int qi = 0; qi < order; ++qi) {
            const double xx = half * gnode[qi] + mid;
            q.x[b][qi] = xx;
            const double fx = (transect == DIST_POINT)
                ? (2.0 * xx / (W * W))         // point: 2x / W^2
                : (1.0 / W);                   // line:  1 / W
            q.base_w[b][qi] = gw[qi] * half * fx;
        }
    }
    return q;
}

// Detection-function value and its eta-derivatives at one node. `g` is g(x);
// `g_e` / `g_ee` are d/d eta_sigma and d^2/d eta_sigma^2; for the hazard key the
// shape derivatives `g_b`, `g_bb`, `g_eb` (d/d eta_b, d^2/d eta_b^2, cross) are
// also filled. Shape outputs are 0 for the half-normal key (it has no b).
struct KeyDeriv {
    double g, g_e, g_ee;          // value, d/d eta_sigma, d^2/d eta_sigma^2
    double g_b, g_bb, g_eb;       // d/d eta_b, d^2/d eta_b^2, d^2/d eta_sigma d eta_b
};

// Detection-function value alone, with none of KeyDeriv's derivative terms --
// for the hazard key this skips the log() that only feeds g_b / g_eb / g_bb.
inline double dist_key_value(double x, int key, double sigma, double b) {
    if (key == DIST_HALFNORMAL) {
        const double u = (x * x) / (sigma * sigma);
        return std::exp(-0.5 * u);
    }
    // hazard-rate: g = 1 - exp(-z), z = (x / sigma)^(-b) = (sigma / x)^b
    if (x <= 0.0) return 1.0;                                  // g -> 1
    const double z = std::pow(sigma / x, b);
    if (!std::isfinite(z) || z > 700.0) return 1.0;           // exp(-z) underflows
    return 1.0 - std::exp(-z);
}

inline KeyDeriv dist_key_deriv(double x, int key, double sigma, double b) {
    KeyDeriv k;
    k.g_b = k.g_bb = k.g_eb = 0.0;
    if (key == DIST_HALFNORMAL) {
        const double u = (x * x) / (sigma * sigma);          // x^2 / sigma^2
        const double g = std::exp(-0.5 * u);
        k.g    = g;
        k.g_e  = g * u;                                       // d/d eta_sigma
        k.g_ee = g * (u * u - 2.0 * u);                      // d^2/d eta_sigma^2
        return k;
    }
    // hazard-rate: g = 1 - exp(-z),  z = (x / sigma)^(-b) = (sigma / x)^b
    if (x <= 0.0) {                                           // g -> 1, derivs -> 0
        k.g = 1.0; k.g_e = k.g_ee = 0.0;
        return k;
    }
    const double ratio = sigma / x;
    const double L = std::log(ratio);                        // ln(sigma / x)
    const double z = std::pow(ratio, b);
    if (!std::isfinite(z) || z > 700.0) {                    // exp(-z) underflows; g -> 1, all derivs -> 0
        k.g = 1.0; k.g_e = k.g_ee = 0.0;
        return k;
    }
    const double w = std::exp(-z);                           // exp(-z)
    const double bz = b * z;
    k.g    = 1.0 - w;
    k.g_e  = w * bz;                                          // d/d eta_sigma
    k.g_ee = b * b * z * w * (1.0 - z);                     // d^2/d eta_sigma^2
    k.g_b  = w * bz * L;                                      // d/d eta_b
    const double common = w * bz * (1.0 + b * L - bz * L);   // d^2/d eta_sigma d eta_b
    k.g_eb = common;
    k.g_bb = L * common;                                     // d^2/d eta_b^2
    return k;
}

// External-pointer wrapping so a fit builds the quadrature ONCE (the R driver
// calls cpp_distance_build_quad() a single time per fit) and every repeated
// .Call() into the compiled kernel reuses it, instead of Newton-Raphson
// root-finding the Gauss-Legendre nodes fresh on every call (gcol33/tulpaObs#165).
inline SEXP dist_quad_wrap(const DistQuad& q) {
    return Rcpp::XPtr<DistQuad>(new DistQuad(q), true);
}

inline Rcpp::XPtr<DistQuad> dist_quad_from_xptr(SEXP quad_xptr) {
    return Rcpp::XPtr<DistQuad>(quad_xptr);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_DISTANCE_QUAD_H
