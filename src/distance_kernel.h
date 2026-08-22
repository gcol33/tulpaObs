// distance_kernel.h
// Per-site marginal log-likelihood + gradients + observed-information pieces for
// binned distance-sampling abundance models with a Poisson OR negative-binomial
// abundance mixing distribution.
//
// Per site i with B distance bins and bin counts y_{i1}, ..., y_{iB}:
//   N_i ~ Poisson(lambda_i)                       (mixture = "P",  r = +Inf)
//   N_i ~ NegBin(mean = lambda_i, size = r)       (mixture = "NB", r finite)
//   lambda_i = exp(eta_lambda_i)                  (expected total in the region)
//   (y_{i1}, ..., y_{iB}, N_i - R_i) ~ Multinomial(N_i, (pi_{i1}, ..., pi_{iB}, 1 - p_i))
//   pi_{ib} = integral_bin g(x; sigma_i[, b]) f(x) dx,   p_i = sum_b pi_{ib}
//   sigma_i = exp(eta_sigma_i),  shape b = exp(eta_b)   (hazard-rate only)
//
// R_i = sum_b y_{ib} is the detected total; an individual lands in bin b and is
// detected with prob pi_{ib}, or is undetected with prob 1 - p_i. The latent N is
// summed out in closed form, and -- as for the N-mixture and removal families --
// the per-N weight a_N is an N-linear slope plus an eta-independent term plus a
// per-N combinatorial term, so the SAME accumulate_count_moments /
// fill_nb_dispersion compute the abundance posterior moments and the NB
// dispersion row/col (nmix_kernel.h). The marginal sums N from R_i to K_max:
//
//   log L_i = LSE_{N=R_i..K_max} { N*slope + base_const + comb_N }
//   slope      = log lambda_i + log(1 - p_i)        (Poisson; NB: eta_lambda - log(r+lambda) + log(1-p))
//   base_const = -lambda_i + sum_b y_{ib} log pi_{ib} - R_i log(1-p_i) - sum_b lgamma(y_{ib}+1)
//   comb_N     = -lgamma(N - R_i + 1)               (NB adds lgamma(N + r); see nmix_kernel.h)
//
// Detection arm (parameters eta_d = (eta_sigma[, eta_b]); pi_{b,k} = d pi_b / d eta_dk
// from quadrature, p_k = sum_b pi_{b,k}). The detection score is affine in N
// (coefficient -B_k = -p_k/(1-p)), so the marginal observed information is the
// Louis (1982) E[I_c|y] minus the rank-1 Var[N|y] v v', with v = (score_wt_lambda,
// -B_sigma[, -B_b]) -- the same structure the (lambda, logit-p) families use, with
// the detection N-coefficient now -p_k/(1-p) instead of -p_ij:
//   d log L_i / d eta_dk = sum_b y_b (pi_{b,k}/pi_b) - (E[N|y]-R) p_k/(1-p)
//   E[I_c]_{jk}          = -sum_b y_b [pi_{b,jk}/pi_b - (pi_{b,j}/pi_b)(pi_{b,k}/pi_b)]
//                          + (E[N|y]-R) [p_{jk}/(1-p) + (p_j p_k)/(1-p)^2]
// A positive-definite Fisher-scoring fallback (info_eta_d_fs) replaces E[I_c] with
// the multinomial expected information E[N|y] [sum_b pi_{b,j}pi_{b,k}/pi_b +
// p_j p_k/(1-p)] for the inner Newton when the observed block is not PSD.
//
// References:
//   Buckland et al. (2001) Introduction to Distance Sampling. OUP.
//   Royle, Dawson & Bates (2004) Ecology 85: 1591-1597.
//   Louis (1982) JRSS-B 44: 226-233.

#ifndef TULPAOBS_DISTANCE_KERNEL_H
#define TULPAOBS_DISTANCE_KERNEL_H

#include "nmix_kernel.h"     // NMixSiteResult, NMixMoments, accumulate_count_moments, fill_nb_dispersion
#include "distance_quad.h"   // DistQuad, DistKey, dist_key_deriv
#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

namespace tulpaObs {

// Combinatorial term of the latent-N marginal sum: the sum runs N = R_i ..
// K_max with K_lo = R_i exactly (the site's own detected total), so the per-N
// term `-lgamma((N - R_i) + 1)` depends only on the offset k = N - R_i, NEVER
// on R_i itself -- it is the SAME table for every site. Building it once per
// fit (size K_max + 1, covering the widest range any site can need) and
// indexing into it replaces up to n_sites * K_max redundant R::lgammafn()
// calls per sweep with array lookups; a caller that sweeps the same K_max
// repeatedly (a Newton iteration, a NUTS leapfrog step, an AGHQ node) builds
// this ONCE outside its loop and passes it to every compute_distance_site()
// call. Byte-identical to the inline computation -- same R::lgammafn(), just
// memoized.
inline std::vector<double> dist_build_comb_table(int K_max) {
    std::vector<double> t((std::size_t) K_max + 1);
    for (int k = 0; k <= K_max; ++k) t[k] = -R::lgammafn((double)k + 1.0);
    return t;
}

// Reusable scratch for compute_distance_site()'s per-call heap allocations
// (the per-bin detection vectors and the per-N marginal-weight vector). A
// caller that evaluates many sites -- or the same site many times at changing
// eta -- builds ONE DistScratch outside its loop(s) and passes it by pointer;
// the vectors are resized (not reallocated) to each call's n_bins / K_grid, so
// after the first call warms them to the loop's maximum sizes no further heap
// traffic occurs. In OpenMP code each thread must own its own DistScratch.
struct DistScratch {
    std::vector<double> pi, d1_0, d2_00, d1_1, d2_11, d2_01, a;
};

// Per-site distance marginal result. The abundance / NB-dispersion fields mirror
// NMixSiteResult; the detection arm carries up to two parameters (sigma, and the
// hazard-rate shape b), so its gradient / Fisher are small fixed-size blocks.
struct DistSiteResult {
    double log_lik;
    double grad_eta_lambda, info_eta_lambda, score_wt_lambda;
    double mean_N, var_N, boundary_weight;
    int n_dparam;                 // 1 (half-normal) or 2 (hazard-rate)
    double grad_eta_d[2];         // detection-parameter gradients (eta scale)
    double info_eta_d[2][2];      // Louis E[I_c|y] detection block (may be indefinite)
    double info_eta_d_fs[2][2];   // PSD Fisher-scoring detection block (fallback)
    double vN_d[2];               // detection score N-coefficients (-p_k/(1-p))
    double p_det;                 // overall detection prob p_i (diagnostics / fitted)
    double grad_theta, info_theta, info_lambda_theta, cov_N_stheta, var_stheta;
};

// Per-site binned distance marginal. `y_bins` are the B bin counts (bin 1 =
// nearest); `eta_sigma` the log-scale detection predictor; `eta_b` the log-shape
// (hazard-rate only, ignored for half-normal); `key` the detection key; `quad`
// the per-fit bin quadrature; `K_max` the marginal-sum truncation (>= sum(y));
// `r` the NB size (+Inf -> Poisson). `value_only`: skip every derivative
// computation (the per-bin quadrature second-derivative accumulation, its five
// per-bin derivative vectors, and the whole detection-arm gradient/Fisher block)
// and fill only `log_lik` (plus the cheap `mean_N` / `var_N` / `boundary_weight`
// / `p_det` moments the log-likelihood sum already produces). Every
// gradient/info field of the returned struct is left at its zero-initialized
// default in this mode -- the caller must not read them. For a caller that only
// needs `ll_cell` (the per-site marginal value at a trial point, e.g. the
// mode-adaptation backtracking line search in R/community_latent.R), this drops
// the 5 extra per-bin heap-allocated vectors and the O(n_dparam^2 * n_bins)
// detection block that `working()`'s score/curvature need but a value lookup
// does not. `headroom`: caps this site's own ceiling at `K_lo + headroom` rather
// than the shared `K_max`, mirroring nmix_precompute_site()'s per-site cap
// (nmix_kernel.h). K_lo == R (the site's detected total) exactly for every
// distance site, so unlike the N-mixture the comb_table stays valid unchanged --
// it is already indexed by the offset `k = N - K_lo`, never by K_lo itself
// (dist_build_comb_table(), #167), so capping here needs only a smaller K_grid,
// no separate per-site cache. A negative headroom (the default) disables the
// cap: every site still sums to the shared K_max, the historical behaviour.
inline DistSiteResult compute_distance_site(
    const int* y_bins, int n_bins,
    double eta_lambda, double eta_sigma, double eta_b,
    int key, const DistQuad& quad, int K_max,
    double r = std::numeric_limits<double>::infinity(),
    bool value_only = false,
    const std::vector<double>* comb_table = nullptr,
    DistScratch* scratch = nullptr,
    int headroom = -1
) {
    const bool is_nb = std::isfinite(r);
    const int nd = (key == DIST_HAZARD) ? 2 : 1;
    DistSiteResult res;
    res.n_dparam = nd;
    for (int j = 0; j < 2; ++j) {
        res.grad_eta_d[j] = res.vN_d[j] = 0.0;
        for (int k = 0; k < 2; ++k) res.info_eta_d[j][k] = res.info_eta_d_fs[j][k] = 0.0;
    }
    res.grad_theta = res.info_theta = res.info_lambda_theta = 0.0;
    res.cov_N_stheta = res.var_stheta = 0.0; res.score_wt_lambda = 1.0;

    int R = 0;
    for (int b = 0; b < n_bins; ++b) R += y_bins[b];
    const int K_lo = R;
    if (K_max < K_lo) {
        res.log_lik = -std::numeric_limits<double>::infinity();
        res.grad_eta_lambda = 0.0; res.info_eta_lambda = 0.0;
        res.mean_N = 0.0; res.var_N = 0.0; res.boundary_weight = 0.0; res.p_det = 0.0;
        return res;
    }

    const double sigma = std::exp(eta_sigma);
    const double b_shape = (key == DIST_HAZARD) ? std::exp(eta_b) : 0.0;

    // Per-bin detection integrals pi_b and, unless `value_only`, their
    // eta-derivatives (summed over the bin's quadrature nodes) and the region
    // totals p, p_k, p_jk. The five per-bin derivative vectors and the second-
    // derivative quadrature accumulation are skipped when only `pi_b` / `p`
    // (and hence `log_lik`) are needed.
    std::vector<double> pi_local, d1_0_local, d2_00_local, d1_1_local, d2_11_local, d2_01_local;
    std::vector<double>& pi    = scratch ? scratch->pi    : pi_local;
    std::vector<double>& d1_0  = scratch ? scratch->d1_0  : d1_0_local;
    std::vector<double>& d2_00 = scratch ? scratch->d2_00 : d2_00_local;
    std::vector<double>& d1_1  = scratch ? scratch->d1_1  : d1_1_local;
    std::vector<double>& d2_11 = scratch ? scratch->d2_11 : d2_11_local;
    std::vector<double>& d2_01 = scratch ? scratch->d2_01 : d2_01_local;
    pi.resize(n_bins);
    if (!value_only) {
        d1_0.assign(n_bins, 0.0); d2_00.assign(n_bins, 0.0);
        d1_1.assign(n_bins, 0.0); d2_11.assign(n_bins, 0.0); d2_01.assign(n_bins, 0.0);
    }
    double p = 0.0, p0 = 0.0, p1 = 0.0;          // p, dp/deta_sigma, dp/deta_b
    double p00 = 0.0, p11 = 0.0, p01 = 0.0;      // second derivatives of p
    for (int b = 0; b < n_bins; ++b) {
        double s = 0.0, s0 = 0.0, s00 = 0.0, s1 = 0.0, s11 = 0.0, s01 = 0.0;
        const std::vector<double>& xb = quad.x[b];
        const std::vector<double>& wb = quad.base_w[b];
        if (value_only) {
            for (int q = 0; q < quad.order; ++q)
                s += wb[q] * dist_key_value(xb[q], key, sigma, b_shape);
        } else {
            for (int q = 0; q < quad.order; ++q) {
                const KeyDeriv k = dist_key_deriv(xb[q], key, sigma, b_shape);
                const double w = wb[q];
                s   += w * k.g;
                s0  += w * k.g_e;
                s00 += w * k.g_ee;
                if (nd == 2) { s1 += w * k.g_b; s11 += w * k.g_bb; s01 += w * k.g_eb; }
            }
        }
        pi[b] = s;
        p += s;
        if (!value_only) {
            d1_0[b] = s0;  d2_00[b] = s00;
            if (nd == 2) { d1_1[b] = s1; d2_11[b] = s11; d2_01[b] = s01; }
            p0  += s0;  p00 += s00;
            if (nd == 2) { p1 += s1; p11 += s11; p01 += s01; }
        }
    }

    // Gauss-Legendre quadrature can round the total detection probability to
    // 1+eps; clamp so 1-p stays positive on a detection-ceiling site.
    if (p > 1.0 - 1e-12) p = 1.0 - 1e-12;
    const double omp = 1.0 - p;                  // 1 - p (undetected prob)
    // A detected bin with vanishing cell probability makes the site impossible.
    bool impossible = (omp <= 0.0);
    for (int b = 0; b < n_bins && !impossible; ++b)
        if (y_bins[b] > 0 && pi[b] <= 0.0) impossible = true;
    if (impossible) {
        res.log_lik = -std::numeric_limits<double>::infinity();
        res.grad_eta_lambda = 0.0; res.info_eta_lambda = 0.0;
        res.mean_N = std::exp(eta_lambda); res.var_N = res.mean_N;
        res.boundary_weight = 0.0; res.p_det = p;
        return res;
    }
    const double log1mp = std::log(omp);
    res.p_det = p;

    double det_const = 0.0, sum_lgam_yfact = 0.0;
    for (int b = 0; b < n_bins; ++b) {
        if (y_bins[b] > 0) det_const += (double)y_bins[b] * std::log(pi[b]);
        sum_lgam_yfact += R::lgammafn((double)y_bins[b] + 1.0);
    }

    const double lambda = std::exp(eta_lambda);
    double slope, base_const;
    if (is_nb) {
        const double log_rpl = std::log(r + lambda);
        slope      = (eta_lambda - log_rpl) + log1mp;
        base_const = -std::lgamma(r) + r * std::log(r) - r * log_rpl
                     + det_const - (double)R * log1mp - sum_lgam_yfact;
    } else {
        slope      = eta_lambda + log1mp;
        base_const = -lambda + det_const - (double)R * log1mp - sum_lgam_yfact;
    }

    // Per-site ceiling K_hi: headroom >= 0 caps this site's sum at K_lo +
    // headroom rather than the shared K_max, whenever that cap is tighter -- a
    // site whose own total already sits within `headroom` of K_max keeps the
    // shared ceiling, never grows past it.
    const int K_hi = (headroom >= 0 && K_lo <= K_max - headroom)
        ? K_lo + headroom : K_max;
    const int K_grid = K_hi - K_lo + 1;
    std::vector<double> a_local;
    std::vector<double>& a = scratch ? scratch->a : a_local;
    a.resize(K_grid);
    // K_lo == R exactly (see above), so N - R == k always: the combinatorial
    // term depends only on the offset k, never on R itself -- comb_table[k],
    // when supplied, is a memoized -lgamma(k+1) built once per fit rather than
    // recomputed at every (site, call) pair. Indexing k up to K_hi - K_lo <=
    // K_max - K_lo stays in the table's bounds regardless of the cap.
    double max_a = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < K_grid; ++k) {
        const int N = K_lo + k;
        const double comb = comb_table ? (*comb_table)[k] : -R::lgammafn((double)k + 1.0);
        a[k] = (double)N * slope + base_const + comb;
        if (is_nb) a[k] += std::lgamma((double)N + r);
        if (a[k] > max_a) max_a = a[k];
    }
    const NMixMoments m =
        accumulate_count_moments(a.data(), K_lo, K_grid, max_a, r, is_nb);
    res.log_lik = m.log_lik; res.mean_N = m.mean_N; res.var_N = m.var_N;
    res.boundary_weight = m.boundary_weight;

    if (value_only) {
        // `log_lik` (and the mean_N/var_N/boundary_weight already filled above)
        // is everything a value-only caller reads; every gradient/info field
        // stays at its zero-initialized default.
        return res;
    }

    // Detection arm. d1[k] / d2[j][k] are the region first/second derivatives of
    // the per-bin / total detection probabilities in eta space.
    const double Nmr = m.mean_N - (double)R;     // E[N|y] - R
    double d1[2]   = { p0, p1 };
    double pjk[2][2] = { { p00, p01 }, { p01, p11 } };
    double B[2]    = { p0 / omp, p1 / omp };     // p_k / (1 - p)
    for (int kk = 0; kk < nd; ++kk) {
        // gradient
        double g = -Nmr * B[kk];
        for (int b = 0; b < n_bins; ++b)
            if (y_bins[b] > 0) {
                const double d1b = (kk == 0) ? d1_0[b] : d1_1[b];
                g += (double)y_bins[b] * (d1b / pi[b]);
            }
        res.grad_eta_d[kk] = g;
        res.vN_d[kk] = -B[kk];
    }
    for (int j = 0; j < nd; ++j) {
        for (int k = 0; k < nd; ++k) {
            double louis = Nmr * (pjk[j][k] / omp + (d1[j] * d1[k]) / (omp * omp));
            double fs    = m.mean_N * (d1[j] * d1[k]) / omp;
            for (int b = 0; b < n_bins; ++b) {
                const double d1bj = (j == 0) ? d1_0[b] : d1_1[b];
                const double d1bk = (k == 0) ? d1_0[b] : d1_1[b];
                const double d2bjk = (j == 0 && k == 0) ? d2_00[b]
                                   : (j == 1 && k == 1) ? d2_11[b] : d2_01[b];
                if (y_bins[b] > 0) {
                    const double rj = d1bj / pi[b], rk = d1bk / pi[b];
                    louis += -(double)y_bins[b] * (d2bjk / pi[b] - rj * rk);
                }
                if (pi[b] > 0.0) fs += m.mean_N * (d1bj * d1bk / pi[b]);
            }
            res.info_eta_d[j][k]    = louis;
            res.info_eta_d_fs[j][k] = fs;
        }
    }

    if (!is_nb) {
        res.grad_eta_lambda = m.mean_N - lambda;
        res.info_eta_lambda = lambda;
        res.score_wt_lambda = 1.0;
        return res;
    }
    NMixSiteResult ab;
    fill_nb_dispersion(ab, lambda, r, m);
    res.grad_eta_lambda   = ab.grad_eta_lambda;
    res.info_eta_lambda   = ab.info_eta_lambda;
    res.score_wt_lambda   = ab.score_wt_lambda;
    res.grad_theta        = ab.grad_theta;
    res.info_theta        = ab.info_theta;
    res.info_lambda_theta = ab.info_lambda_theta;
    res.cov_N_stheta      = ab.cov_N_stheta;
    res.var_stheta        = ab.var_stheta;
    return res;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_DISTANCE_KERNEL_H
