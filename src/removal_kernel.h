// removal_kernel.h
// Per-site marginal log-likelihood + gradients + Fisher info for the
// removal-sampling (sequential-depletion) abundance model with a Poisson OR
// negative-binomial abundance mixing distribution.
//
// Per site i with K ordered removal passes and removals y_{i1}, ..., y_{iK}:
//   N_i ~ Poisson(lambda_i)                         (mixture = "P",  r = +Inf)
//   N_i ~ NegBin(mean = lambda_i, size = r)         (mixture = "NB", r finite)
//   lambda_i = exp(eta_lambda_i)
//   y_{ik} | (N_i, y_{i,<k}) ~ Binomial(A_{ik}, p_{ik}),
//       A_{ik} = N_i - sum_{l<k} y_{il}   (individuals still present at pass k)
//   p_{ik} = plogis(eta_p_{ik})
//
// Each pass removes the detected individuals, so the trials available at pass k
// deplete by the cumulative prior removals. The product of the depleting
// binomials equals the multinomial-removal likelihood (Royle 2004,
// unmarked::multinomPois with removalPiFun): an individual is first removed at
// pass k with probability pi_k = p_k prod_{l<k}(1-p_l), or never with
// pi_0 = prod_l (1-p_l). The marginal sums the (unobserved) total abundance out:
//
//   log L_i = LSE_{N = R_i .. K_max} { log f(N | lambda_i, r) + log mult(y_i | N) }
//
// with R_i = sum_k y_{ik} the total removed (the smallest admissible N) and
// f = Poisson or NB. The per-N weight is laid out exactly like the N-mixture
// (nmix_kernel.h): an N-linear slope, an eta-independent combinatorial term, and
// a base constant, so the SAME accumulate_count_moments / fill_nb_dispersion the
// N-mixture uses compute the abundance posterior moments and the NB dispersion
// row/col. Only the per-N combinatorics and the detection arm differ, because
// pass k sees the depleted count A_{ik} rather than the full N (Poisson:
// E[A_{ik}] = E[N|y] - C_{i,<k}, C_{i,<k} = sum_{l<k} y_{il}):
//
//   d log L_i / d eta_lambda  = E[N|y_i] - lambda_i           (Poisson; NB as nmix)
//   d log L_i / d eta_p_{ik}  = y_{ik} - p_{ik} (E[N|y_i] - C_{i,<k})
//   I_{p_ik, p_ik}            = (E[N|y_i] - C_{i,<k}) p_{ik}(1-p_{ik})
//
// The marginal observed information is the complete-data Fisher diagonal minus
// the rank-1 shared-latent score covariance Var[N|y] v v', with the same v as
// the N-mixture (v_lambda = score weight, v_p_ik = p_ik): the detection score's
// N-coefficient is -p_ik regardless of the constant depletion offset (Louis
// 1982). So the Laplace driver (marginal_count_laplace.h) is shared verbatim.
//
// References:
//   Royle (2004) Biometrics 60: 108-115 (N-mixture; multinomPois removal).
//   Dorazio, Jelks & Jordan (2005) Biometrics 61: 1093-1101 (removal N-mixture).
//   Louis (1982) JRSS-B 44: 226-233.

#ifndef TULPAOBS_REMOVAL_KERNEL_H
#define TULPAOBS_REMOVAL_KERNEL_H

#include "nmix_kernel.h"   // NMixSiteResult, logit_log_probs, accumulate_count_moments, fill_nb_dispersion
#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

namespace tulpaObs {

// Per-site removal marginal. `y` are the per-pass removals in pass order (pass 1
// first); `eta_p` the per-pass detection logit predictor in the same order;
// `n_pass` the number of passes K; `eta_lambda` the log abundance predictor;
// `K_max` the marginal-sum truncation (must be >= sum(y)); `r` the NB size
// (+Inf -> Poisson). The pass order matters: depletion accumulates the removals
// of earlier passes, so callers must pass the passes in order with no gaps.
inline NMixSiteResult compute_removal_site(
    const int* y,
    const double* eta_p,
    int n_pass,
    double eta_lambda,
    int K_max,
    double r = std::numeric_limits<double>::infinity()
) {
    const bool is_nb = std::isfinite(r);
    NMixSiteResult res;
    res.grad_eta_p.assign(n_pass, 0.0);
    res.info_eta_p.assign(n_pass, 0.0);
    res.grad_theta = 0.0; res.info_theta = 0.0; res.info_lambda_theta = 0.0;
    res.cov_N_stheta = 0.0; res.var_stheta = 0.0; res.score_wt_lambda = 1.0;

    // Cumulative prior removals C_{<k} (offset) and total removed R = K_lo.
    std::vector<int> offset(n_pass, 0);
    int R = 0;
    for (int k = 0; k < n_pass; ++k) { offset[k] = R; R += y[k]; }
    const int K_lo = R;
    const bool admissible = (K_max >= K_lo);
    if (!admissible) {
        res.log_lik = -std::numeric_limits<double>::infinity();
        res.grad_eta_lambda = 0.0; res.info_eta_lambda = 0.0;
        res.mean_N = 0.0; res.var_N = 0.0; res.boundary_weight = 0.0;
        return res;
    }

    const double lambda = std::exp(eta_lambda);
    std::vector<double> p_vec(n_pass);
    double sum_log_1mp = 0.0;        // sum_k log(1 - p_k)  (the N-linear det term)
    double det_const   = 0.0;        // sum_k [y_k log p_k - C_{<=k} log(1-p_k)]
    double const_log_yfact = 0.0;    // -sum_k lgamma(y_k + 1)
    for (int k = 0; k < n_pass; ++k) {
        double lp, l1mp;
        logit_log_probs(eta_p[k], lp, l1mp);
        sum_log_1mp += l1mp;
        const int C_le = offset[k] + y[k];          // C_{<=k} = sum_{l<=k} y_l
        det_const += (double)y[k] * lp - (double)C_le * l1mp;
        const_log_yfact -= R::lgammafn((double)y[k] + 1.0);
        if (eta_p[k] > 0.0) p_vec[k] = 1.0 / (1.0 + std::exp(-eta_p[k]));
        else { double e = std::exp(eta_p[k]); p_vec[k] = e / (1.0 + e); }
    }

    // Abundance prior: Poisson vs NB changes the N-slope and the base constant;
    // NB also carries the per-N lgamma(N + r) (added in the grid loop). The
    // detection contribution (det_const, sum_log_1mp) is identical for both.
    double slope, base_const;
    if (is_nb) {
        const double log_rpl = std::log(r + lambda);
        slope      = (eta_lambda - log_rpl) + sum_log_1mp;
        base_const = -std::lgamma(r) + r * std::log(r) - r * log_rpl
                     + det_const + const_log_yfact;
    } else {
        slope      = eta_lambda + sum_log_1mp;
        base_const = -lambda + det_const + const_log_yfact;
    }

    const int K_grid = K_max - K_lo + 1;
    std::vector<double> a(K_grid);
    double max_a = -std::numeric_limits<double>::infinity();
    for (int k = 0; k < K_grid; ++k) {
        const int N = K_lo + k;
        // Combinatorial term: sum_pass lgamma(N - off + 1) - lgamma(N+1)
        //                     - sum_pass lgamma(N - off - y + 1)   (eta-independent)
        double term = -R::lgammafn((double)N + 1.0);
        for (int j = 0; j < n_pass; ++j) {
            term += R::lgammafn((double)(N - offset[j]) + 1.0)
                  - R::lgammafn((double)(N - offset[j] - y[j]) + 1.0);
        }
        a[k] = (double)N * slope + base_const + term;
        if (is_nb) a[k] += std::lgamma((double)N + r);
        if (a[k] > max_a) max_a = a[k];
    }
    const NMixMoments m =
        accumulate_count_moments(a.data(), K_lo, K_grid, max_a, r, is_nb);
    res.log_lik = m.log_lik; res.mean_N = m.mean_N; res.var_N = m.var_N;
    res.boundary_weight = m.boundary_weight;

    // Detection arm: pass k sees the depleted available count A_k = N - C_{<k},
    // so the binomial score / info use E[A_k] = E[N|y] - offset_k.
    for (int k = 0; k < n_pass; ++k) {
        const double avail = std::max(m.mean_N - (double)offset[k], 0.0);
        res.grad_eta_p[k] = (double)y[k] - avail * p_vec[k];
        res.info_eta_p[k] = avail * p_vec[k] * (1.0 - p_vec[k]);
    }

    if (!is_nb) {
        res.grad_eta_lambda = m.mean_N - lambda;
        res.info_eta_lambda = lambda;
        res.score_wt_lambda = 1.0;
        return res;
    }
    fill_nb_dispersion(res, lambda, r, m);
    return res;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_REMOVAL_KERNEL_H
