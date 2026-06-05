// dyn_abun_kernel.h
// Per-site marginal log-likelihood + analytic gradients for the Dail-Madsen
// (2011) open-population N-mixture, by an exact HMM forward recursion over the
// latent abundance states 0..K with forward-mode differentiation.
//
// Per site i, primary seasons t = 1..T, secondary visits j:
//   N_{i,1} ~ Poisson(lambda_i),          log lambda_i = X_lambda beta_lambda
//   N_{i,t} = S_{i,t} + G_{i,t}  (t >= 2):
//     S_{i,t} ~ Binomial(N_{i,t-1}, omega_i)   apparent survival, logit omega
//     G_{i,t} ~ Poisson(gamma_i)               recruitment (constant), log gamma
//   y_{i,t,j} | N_{i,t} ~ Binomial(N_{i,t}, p_i),  logit p_i = X_p beta_p
//
// The latent abundance sequence is summed out by the forward algorithm with a
// per-season transition matrix Tr(n -> n') = sum_s Binom(s|n,omega)
// Poisson(n'-s|gamma) (the survival-binomial convolved with the recruitment-
// Poisson). The marginal is not closed form (unlike the static N-mixture), so the
// gradient is obtained by propagating d alpha_t(n) / d eta_k alongside the scaled
// forward state for each of the four arms (eta_lambda, eta_p, eta_omega,
// eta_gamma). Only the transition depends on (eta_omega, eta_gamma); only the
// observation on eta_p; only the initial on eta_lambda. Per-season scaling keeps
// the recursion numerically stable, and log L = sum_t log c_t with
// d log L / d eta_k = sum_t (d c_t / d eta_k) / c_t.
//
// All four arms are site-level here (one logit/log predictor per arm per site;
// constant across a site's seasons), so every gradient is a scalar per site.
// Poisson initial abundance and constant recruitment (the pcountOpen
// "constant" dynamics); a missing visit (y < 0) is dropped from its season.
//
// References:
//   Dail, D., Madsen, L. (2011) Biometrics 67: 577-587.
//   Royle, J. A. (2004) Biometrics 60: 108-115 (the static observation arm).

#ifndef TULPAOBS_DYN_ABUN_KERNEL_H
#define TULPAOBS_DYN_ABUN_KERNEL_H

#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>

namespace tulpaObs {

struct DynAbunSiteResult {
    double log_lik;
    double grad_eta_lambda, grad_eta_p, grad_eta_omega, grad_eta_gamma;
    double mean_N1;            // E[N_1 | y] (diagnostic / fitted)
};

inline double da_inv_logit(double e) {
    if (e > 0.0) return 1.0 / (1.0 + std::exp(-e));
    double ee = std::exp(e);
    return ee / (1.0 + ee);
}

// Per-site Dail-Madsen forward marginal. `y` is laid out season-major
// (y[t*J + j]) with -1 marking a missing visit; T seasons, J secondary visits,
// K the abundance truncation (states 0..K).
inline DynAbunSiteResult compute_dyn_abun_site(
    const int* y, int T, int J, int K,
    double eta_lambda, double eta_p, double eta_omega, double eta_gamma
) {
    const int S = K + 1;                       // number of abundance states
    const double lambda = std::exp(eta_lambda);
    const double p      = da_inv_logit(eta_p);
    const double omega  = da_inv_logit(eta_omega);
    const double gamma  = std::exp(eta_gamma);
    const double logp = std::log(p), log1mp = std::log1p(-p);
    const double logom = std::log(omega), log1mom = std::log1p(-omega);
    const double loggam = std::log(gamma);

    DynAbunSiteResult res;
    res.grad_eta_lambda = res.grad_eta_p = res.grad_eta_omega = res.grad_eta_gamma = 0.0;

    // Recruitment pmf pois(g) = Poisson(g | gamma) and d/d eta_gamma = pois*(g-gamma).
    std::vector<double> pois(S), dpois_g(S);
    for (int g = 0; g < S; ++g) {
        const double lp = -gamma + (double)g * loggam - R::lgammafn((double)g + 1.0);
        pois[g] = std::exp(lp);
        dpois_g[g] = pois[g] * ((double)g - gamma);
    }

    // Observation obs_t(n) = prod_j Binom(y_tj | n, p) and d/d eta_p =
    // obs * sum_j (y_tj - n p). Precompute per season: max y, sum y, count.
    // obs(n) = 0 for n < max_j y_tj.
    auto obs_season = [&](int t, std::vector<double>& obs, std::vector<double>& dobs) {
        int ymax = 0, ysum = 0, nv = 0;
        for (int j = 0; j < J; ++j) {
            const int yy = y[t * J + j];
            if (yy < 0) continue;
            nv++; ysum += yy; if (yy > ymax) ymax = yy;
        }
        for (int n = 0; n < S; ++n) {
            if (n < ymax) { obs[n] = 0.0; dobs[n] = 0.0; continue; }
            double lo = 0.0;
            for (int j = 0; j < J; ++j) {
                const int yy = y[t * J + j];
                if (yy < 0) continue;
                lo += R::lgammafn((double)n + 1.0) - R::lgammafn((double)yy + 1.0)
                    - R::lgammafn((double)(n - yy) + 1.0)
                    + (double)yy * logp + (double)(n - yy) * log1mp;
            }
            obs[n] = std::exp(lo);
            dobs[n] = obs[n] * ((double)ysum - (double)nv * (double)n * p);
        }
    };

    // Forward state and its 4 derivatives (alpha normalised after each season).
    std::vector<double> a(S), da_l(S), da_p(S), da_o(S), da_g(S);
    std::vector<double> obs(S), dobs(S);

    // --- Season 1: initial Poisson(lambda) x observation. ---
    obs_season(0, obs, dobs);
    double c1 = 0.0, dc_l = 0.0, dc_p = 0.0;
    for (int n = 0; n < S; ++n) {
        const double lpn = -lambda + (double)n * eta_lambda - R::lgammafn((double)n + 1.0);
        const double pi_n = std::exp(lpn);
        const double dpi_l = pi_n * ((double)n - lambda);     // d/d eta_lambda
        a[n]   = pi_n * obs[n];
        da_l[n] = dpi_l * obs[n];
        da_p[n] = pi_n * dobs[n];
        da_o[n] = 0.0; da_g[n] = 0.0;
        c1 += a[n]; dc_l += da_l[n]; dc_p += da_p[n];
    }
    if (!(c1 > 0.0)) {           // impossible history (e.g. counts above K)
        res.log_lik = -std::numeric_limits<double>::infinity();
        res.mean_N1 = 0.0; return res;
    }
    double log_lik = std::log(c1);
    res.grad_eta_lambda += dc_l / c1;
    res.grad_eta_p      += dc_p / c1;
    // mean_N1 from the season-1 posterior (before normalisation effects cancel).
    double mN1 = 0.0;
    for (int n = 0; n < S; ++n) mN1 += (double)n * (a[n] / c1);
    res.mean_N1 = mN1;
    // Normalise alpha and its derivatives.
    for (int n = 0; n < S; ++n) {
        a[n]   /= c1;
        da_l[n] = (da_l[n] - a[n] * dc_l) / c1;
        da_p[n] = (da_p[n] - a[n] * dc_p) / c1;
        da_o[n] = 0.0; da_g[n] = 0.0;       // no omega/gamma dependence yet
    }

    // --- Seasons 2..T: transition then observation. ---
    std::vector<double> pre(S), dpre_l(S), dpre_p(S), dpre_o(S), dpre_g(S);
    std::vector<double> binom(S), dbinom(S);   // survivor pmf for the current row n
    for (int t = 1; t < T; ++t) {
        for (int n2 = 0; n2 < S; ++n2) {
            pre[n2] = dpre_l[n2] = dpre_p[n2] = dpre_o[n2] = dpre_g[n2] = 0.0;
        }
        for (int n = 0; n < S; ++n) {
            if (a[n] == 0.0 && da_l[n] == 0.0 && da_p[n] == 0.0 &&
                da_o[n] == 0.0 && da_g[n] == 0.0) continue;
            // Survivor pmf Binom(s | n, omega), s = 0..n, and d/d eta_omega =
            // binom * (s - n omega).
            for (int s = 0; s <= n; ++s) {
                const double lb = R::lgammafn((double)n + 1.0) - R::lgammafn((double)s + 1.0)
                    - R::lgammafn((double)(n - s) + 1.0)
                    + (double)s * logom + (double)(n - s) * log1mom;
                binom[s] = std::exp(lb);
                dbinom[s] = binom[s] * ((double)s - (double)n * omega);
            }
            const double an = a[n], dl = da_l[n], dp = da_p[n], dom = da_o[n], dg = da_g[n];
            // Tr(n -> n') = sum_s binom[s] pois[n'-s]; scatter into pre.
            for (int s = 0; s <= n; ++s) {
                const double bs = binom[s], dbs = dbinom[s];
                for (int gn = 0; gn + s < S; ++gn) {
                    const int n2 = s + gn;
                    const double tr = bs * pois[gn];
                    pre[n2]   += an * tr;
                    // chain rule: d(alpha_{t-1}(n) Tr)/d eta_k
                    dpre_l[n2] += dl * tr;
                    dpre_p[n2] += dp * tr;
                    dpre_o[n2] += dom * tr + an * (dbs * pois[gn]);
                    dpre_g[n2] += dg * tr + an * (bs * dpois_g[gn]);
                }
            }
        }
        obs_season(t, obs, dobs);
        double ct = 0.0, dct_l = 0.0, dct_p = 0.0, dct_o = 0.0, dct_g = 0.0;
        for (int n2 = 0; n2 < S; ++n2) {
            a[n2]   = pre[n2] * obs[n2];
            da_l[n2] = dpre_l[n2] * obs[n2];
            da_p[n2] = dpre_p[n2] * obs[n2] + pre[n2] * dobs[n2];
            da_o[n2] = dpre_o[n2] * obs[n2];
            da_g[n2] = dpre_g[n2] * obs[n2];
            ct += a[n2]; dct_l += da_l[n2]; dct_p += da_p[n2];
            dct_o += da_o[n2]; dct_g += da_g[n2];
        }
        if (!(ct > 0.0)) {
            res.log_lik = -std::numeric_limits<double>::infinity();
            return res;
        }
        log_lik += std::log(ct);
        res.grad_eta_lambda += dct_l / ct;
        res.grad_eta_p      += dct_p / ct;
        res.grad_eta_omega  += dct_o / ct;
        res.grad_eta_gamma  += dct_g / ct;
        for (int n2 = 0; n2 < S; ++n2) {
            a[n2]   /= ct;
            da_l[n2] = (da_l[n2] - a[n2] * dct_l) / ct;
            da_p[n2] = (da_p[n2] - a[n2] * dct_p) / ct;
            da_o[n2] = (da_o[n2] - a[n2] * dct_o) / ct;
            da_g[n2] = (da_g[n2] - a[n2] * dct_g) / ct;
        }
    }

    res.log_lik = log_lik;
    return res;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_DYN_ABUN_KERNEL_H
