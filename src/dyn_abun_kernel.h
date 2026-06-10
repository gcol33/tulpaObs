// dyn_abun_kernel.h
// Per-site marginal log-likelihood + analytic gradients for the Dail-Madsen
// (2011) open-population N-mixture, by an exact HMM forward recursion over the
// latent abundance states 0..K with forward-mode differentiation.
//
// Per site i, primary seasons t = 1..T, secondary visits j:
//   N_{i,1} ~ Poisson(lambda_i),          log lambda_i = X_lambda beta_lambda
//   N_{i,t} = S_{i,t} + G_{i,t}  (t >= 2):
//     S_{i,t} ~ Binomial(N_{i,t-1}, omega_{i,t-1})  apparent survival, logit omega
//     G_{i,t} ~ Poisson(gamma_{i,t-1})              recruitment, log gamma
//   y_{i,t,j} | N_{i,t} ~ Binomial(N_{i,t}, p_i),  logit p_i = X_p beta_p
//
// The transition from season t-1 to t uses the survival and recruitment rates
// of the INTERVAL (t-1): there are T-1 interval-specific eta_omega / eta_gamma
// per site, so a season covariate on omega / gamma drives the dynamics. The
// latent abundance sequence is summed out by the forward algorithm with the
// per-interval transition matrix
// Tr_{t-1}(n -> n') = sum_s Binom(s | n, omega_{t-1}) Poisson(n'-s | gamma_{t-1})
// (the survival-binomial convolved with the recruitment-Poisson). The marginal
// is not closed form (unlike the static N-mixture), so the gradient is obtained
// by propagating d alpha_t(n) / d eta_k alongside the scaled forward state.
// Only the transition depends on (eta_omega, eta_gamma); only the observation on
// eta_p; only the initial on eta_lambda. Per-season scaling keeps the recursion
// numerically stable, and log L = sum_t log c_t with
// d log L / d eta_k = sum_t (d c_t / d eta_k) / c_t.
//
// The eta_lambda / eta_p / eta_logr arms are site-level (one predictor per site),
// so their gradients are scalars per site. The eta_omega / eta_gamma arms are
// interval-level: a length-(T-1) gradient vector each, since the t-th interval's
// rate enters only its own transition. The forward-mode pass therefore carries
// T-1 derivative directions for omega and for gamma -- direction iv is "born" at
// the transition that uses interval iv (from an * d Tr_iv) and then propagates
// through every later transition / observation unchanged (linear operators),
// exactly as the initial-only lambda direction does. Poisson initial abundance
// and a missing visit (y < 0) is dropped from its season.
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
    // Per-interval survival / recruitment scores, length T-1. grad_eta_omega /
    // grad_eta_gamma are their sums (the scalar constant-rate score).
    std::vector<double> grad_eta_omega_vec, grad_eta_gamma_vec;
    double grad_eta_logr;      // d log L / d log r (negbin initial; 0 under Poisson)
    double mean_N1;            // E[N_1 | y] (diagnostic / fitted)
};

inline double da_inv_logit(double e) {
    if (e > 0.0) return 1.0 / (1.0 + std::exp(-e));
    double ee = std::exp(e);
    return ee / (1.0 + ee);
}

// Shared per-season observation pmf obs[n] = prod_j Binom(y_tj | n, p), with
// obs[n] = 0 for n < max_j y_tj. Season t is laid out y[t*J + j], -1 = missing.
// Single source of truth for the forward gradient kernel and the backward
// conditional-likelihood pass below.
inline void da_obs_season_pmf(const int* y, int t, int J, int S,
                              double logp, double log1mp, std::vector<double>& obs) {
    int ymax = 0;
    for (int j = 0; j < J; ++j) {
        const int yy = y[t * J + j];
        if (yy >= 0 && yy > ymax) ymax = yy;
    }
    for (int n = 0; n < S; ++n) {
        if (n < ymax) { obs[n] = 0.0; continue; }
        double lo = 0.0;
        for (int j = 0; j < J; ++j) {
            const int yy = y[t * J + j];
            if (yy < 0) continue;
            lo += R::lgammafn((double)n + 1.0) - R::lgammafn((double)yy + 1.0)
                - R::lgammafn((double)(n - yy) + 1.0)
                + (double)yy * logp + (double)(n - yy) * log1mp;
        }
        obs[n] = std::exp(lo);
    }
}

// Shared recruitment pmf pois[g] = Poisson(g | gamma), g = 0..S-1.
inline void da_recruit_pmf(int S, double gamma, double loggam,
                           std::vector<double>& pois) {
    for (int g = 0; g < S; ++g) {
        const double lp = -gamma + (double)g * loggam - R::lgammafn((double)g + 1.0);
        pois[g] = std::exp(lp);
    }
}

// Shared survivor pmf binom[s] = Binom(s | n, omega), s = 0..n.
inline void da_binom_pmf_row(int n, double logom, double log1mom,
                             std::vector<double>& binom) {
    for (int s = 0; s <= n; ++s) {
        const double lb = R::lgammafn((double)n + 1.0) - R::lgammafn((double)s + 1.0)
            - R::lgammafn((double)(n - s) + 1.0)
            + (double)s * logom + (double)(n - s) * log1mom;
        binom[s] = std::exp(lb);
    }
}

// Core forward recursion with interval-indexed survival / recruitment rates:
// `eta_omega` / `eta_gamma` point to length-(T-1) arrays, eta[iv] the rate for
// the transition from season iv to iv+1. `use_nb` switches the season-1 initial
// abundance from Poisson(lambda) to negative-binomial NB(mean = lambda, size =
// r), r = exp(eta_logr) -- the pcountOpen "NB" mixture. Only the initial
// distribution changes; survival, recruitment and detection are identical. The
// dispersion enters the marginal solely through season 1, so its derivative
// `da_r` propagates exactly like `da_l` (lambda), which is also initial-only.
//
// The forward-mode derivative state for omega / gamma is per-interval: da_o[iv]
// and da_g[iv] are the length-S derivative of the (normalised) forward state in
// the direction eta_omega[iv] / eta_gamma[iv]. Direction iv is zero until the
// transition that uses interval iv, where the survivor / recruitment pmf carries
// the only eta[iv] dependence; from then on it propagates through later
// transitions and observations like any other forward derivative. res grads:
// grad_eta_lambda / grad_eta_p / grad_eta_logr are scalars; grad_eta_omega_vec /
// grad_eta_gamma_vec are length T-1 (and grad_eta_omega / grad_eta_gamma their
// sums, the constant-rate scalar score). With T-1 == 1 the constant-rate path is
// recovered (the per-interval vectors hold one element).
inline DynAbunSiteResult compute_dyn_abun_site(
    const int* y, int T, int J, int K,
    double eta_lambda, double eta_p,
    const double* eta_omega, const double* eta_gamma,
    bool use_nb = false, double eta_logr = 0.0
) {
    const int S = K + 1;                       // number of abundance states
    const int nIv = T - 1;                      // number of transition intervals
    const double lambda = std::exp(eta_lambda);
    const double p      = da_inv_logit(eta_p);
    const double logp = std::log(p), log1mp = std::log1p(-p);
    const double rr = use_nb ? std::exp(eta_logr) : 0.0;  // NB size

    // Per-interval survival omega_iv, recruitment gamma_iv and their pmfs.
    std::vector<double> omega(nIv), gamma(nIv), logom(nIv), log1mom(nIv), loggam(nIv);
    std::vector<std::vector<double> > pois(nIv), dpois_g(nIv);
    for (int iv = 0; iv < nIv; ++iv) {
        omega[iv]  = da_inv_logit(eta_omega[iv]);
        gamma[iv]  = std::exp(eta_gamma[iv]);
        logom[iv]  = std::log(omega[iv]); log1mom[iv] = std::log1p(-omega[iv]);
        loggam[iv] = std::log(gamma[iv]);
        pois[iv].resize(S); dpois_g[iv].resize(S);
        da_recruit_pmf(S, gamma[iv], loggam[iv], pois[iv]);
        for (int g = 0; g < S; ++g) dpois_g[iv][g] = pois[iv][g] * ((double)g - gamma[iv]);
    }

    DynAbunSiteResult res;
    res.grad_eta_lambda = res.grad_eta_p = res.grad_eta_omega = res.grad_eta_gamma = 0.0;
    res.grad_eta_logr = 0.0;
    res.grad_eta_omega_vec.assign(nIv, 0.0);
    res.grad_eta_gamma_vec.assign(nIv, 0.0);

    // Observation obs_t(n) = prod_j Binom(y_tj | n, p) and d/d eta_p =
    // obs * sum_j (y_tj - n p). Precompute per season: max y, sum y, count.
    auto obs_season = [&](int t, std::vector<double>& obs, std::vector<double>& dobs) {
        int ysum = 0, nv = 0;
        for (int j = 0; j < J; ++j) {
            const int yy = y[t * J + j];
            if (yy < 0) continue;
            nv++; ysum += yy;
        }
        da_obs_season_pmf(y, t, J, S, logp, log1mp, obs);
        for (int n = 0; n < S; ++n)
            dobs[n] = obs[n] * ((double)ysum - (double)nv * (double)n * p);
    };

    // Forward state and its derivatives (alpha normalised after each season). da_r
    // is the dispersion derivative, carried alongside da_l (both initial-only). The
    // omega / gamma derivatives are per-interval: da_o[iv], da_g[iv] each length S.
    std::vector<double> a(S), da_l(S), da_p(S), da_r(S);
    std::vector<std::vector<double> > da_o(nIv), da_g(nIv);
    for (int iv = 0; iv < nIv; ++iv) { da_o[iv].assign(S, 0.0); da_g[iv].assign(S, 0.0); }
    std::vector<double> obs(S), dobs(S);

    // --- Season 1: initial Poisson / NB (lambda[, r]) x observation. ---
    obs_season(0, obs, dobs);
    double c1 = 0.0, dc_l = 0.0, dc_p = 0.0, dc_r = 0.0;
    for (int n = 0; n < S; ++n) {
        double pi_n, dpi_l, dpi_r;
        if (use_nb) {
            // NB(mean = lambda, size = rr): log pi = lgamma(n+r) - lgamma(r)
            //   - lgamma(n+1) + r log(r/(r+mu)) + n log(mu/(r+mu)).
            const double rpm = rr + lambda;
            const double lpn = R::lgammafn((double)n + rr) - R::lgammafn(rr)
                - R::lgammafn((double)n + 1.0)
                + rr * std::log(rr / rpm) + (double)n * std::log(lambda / rpm);
            pi_n = std::exp(lpn);
            // d log pi / d eta_lambda = n - mu (n+r)/(r+mu).
            dpi_l = pi_n * ((double)n - lambda * ((double)n + rr) / rpm);
            // d log pi / d log r = r [psi(n+r) - psi(r) + log(r/(r+mu)) + 1
            //   - (r+n)/(r+mu)].
            const double dlog_dlogr = rr * (R::digamma((double)n + rr) - R::digamma(rr)
                + std::log(rr / rpm) + 1.0 - (rr + (double)n) / rpm);
            dpi_r = pi_n * dlog_dlogr;
        } else {
            const double lpn = -lambda + (double)n * eta_lambda - R::lgammafn((double)n + 1.0);
            pi_n = std::exp(lpn);
            dpi_l = pi_n * ((double)n - lambda);     // d/d eta_lambda
            dpi_r = 0.0;
        }
        a[n]   = pi_n * obs[n];
        da_l[n] = dpi_l * obs[n];
        da_p[n] = pi_n * dobs[n];
        da_r[n] = dpi_r * obs[n];
        c1 += a[n]; dc_l += da_l[n]; dc_p += da_p[n]; dc_r += da_r[n];
    }
    if (!(c1 > 0.0)) {           // impossible history (e.g. counts above K)
        res.log_lik = -std::numeric_limits<double>::infinity();
        res.mean_N1 = 0.0; return res;
    }
    double log_lik = std::log(c1);
    res.grad_eta_lambda += dc_l / c1;
    res.grad_eta_p      += dc_p / c1;
    res.grad_eta_logr   += dc_r / c1;
    // mean_N1 from the season-1 posterior (before normalisation effects cancel).
    double mN1 = 0.0;
    for (int n = 0; n < S; ++n) mN1 += (double)n * (a[n] / c1);
    res.mean_N1 = mN1;
    // Normalise alpha and its derivatives (omega/gamma directions still zero).
    for (int n = 0; n < S; ++n) {
        a[n]   /= c1;
        da_l[n] = (da_l[n] - a[n] * dc_l) / c1;
        da_p[n] = (da_p[n] - a[n] * dc_p) / c1;
        da_r[n] = (da_r[n] - a[n] * dc_r) / c1;
    }

    // --- Seasons 2..T: transition (interval iv = t-1) then observation. ---
    std::vector<double> pre(S), dpre_l(S), dpre_p(S), dpre_r(S);
    std::vector<std::vector<double> > dpre_o(nIv), dpre_g(nIv);
    for (int iv = 0; iv < nIv; ++iv) { dpre_o[iv].resize(S); dpre_g[iv].resize(S); }
    std::vector<double> binom(S), dbinom(S);   // survivor pmf for the current row n
    for (int t = 1; t < T; ++t) {
        const int iv = t - 1;                       // active interval for this step
        const double om = omega[iv];
        const std::vector<double>& pois_iv    = pois[iv];
        const std::vector<double>& dpois_g_iv = dpois_g[iv];
        for (int n2 = 0; n2 < S; ++n2) {
            pre[n2] = dpre_l[n2] = dpre_p[n2] = dpre_r[n2] = 0.0;
            for (int jv = 0; jv < nIv; ++jv) { dpre_o[jv][n2] = 0.0; dpre_g[jv][n2] = 0.0; }
        }
        for (int n = 0; n < S; ++n) {
            // Survivor pmf Binom(s | n, omega_iv) and d/d eta_omega = binom*(s - n omega).
            da_binom_pmf_row(n, logom[iv], log1mom[iv], binom);
            for (int s = 0; s <= n; ++s) dbinom[s] = binom[s] * ((double)s - (double)n * om);
            const double an = a[n], dl = da_l[n], dp = da_p[n], dr = da_r[n];
            // Tr(n -> n') = sum_s binom[s] pois_iv[n'-s]; scatter into pre.
            for (int s = 0; s <= n; ++s) {
                const double bs = binom[s], dbs = dbinom[s];
                for (int gn = 0; gn + s < S; ++gn) {
                    const int n2 = s + gn;
                    const double tr = bs * pois_iv[gn];
                    pre[n2]   += an * tr;
                    // chain rule: d(alpha_{t-1}(n) Tr_iv)/d eta_k. The transition
                    // carries no lambda/r dependence, so those propagate the
                    // upstream derivative alone (like lambda).
                    dpre_l[n2] += dl * tr;
                    dpre_p[n2] += dp * tr;
                    dpre_r[n2] += dr * tr;
                    // omega / gamma interval iv is born here (an * d Tr_iv); every
                    // other interval propagates its upstream derivative through Tr_iv.
                    dpre_o[iv][n2] += an * (dbs * pois_iv[gn]);
                    dpre_g[iv][n2] += an * (bs * dpois_g_iv[gn]);
                    for (int jv = 0; jv < nIv; ++jv) {
                        if (da_o[jv][n] != 0.0) dpre_o[jv][n2] += da_o[jv][n] * tr;
                        if (da_g[jv][n] != 0.0) dpre_g[jv][n2] += da_g[jv][n] * tr;
                    }
                }
            }
        }
        obs_season(t, obs, dobs);
        double ct = 0.0, dct_l = 0.0, dct_p = 0.0, dct_r = 0.0;
        std::vector<double> dct_o(nIv, 0.0), dct_g(nIv, 0.0);
        for (int n2 = 0; n2 < S; ++n2) {
            a[n2]   = pre[n2] * obs[n2];
            da_l[n2] = dpre_l[n2] * obs[n2];
            da_p[n2] = dpre_p[n2] * obs[n2] + pre[n2] * dobs[n2];
            da_r[n2] = dpre_r[n2] * obs[n2];
            ct += a[n2]; dct_l += da_l[n2]; dct_p += da_p[n2]; dct_r += da_r[n2];
            for (int jv = 0; jv < nIv; ++jv) {
                da_o[jv][n2] = dpre_o[jv][n2] * obs[n2];
                da_g[jv][n2] = dpre_g[jv][n2] * obs[n2];
                dct_o[jv] += da_o[jv][n2]; dct_g[jv] += da_g[jv][n2];
            }
        }
        if (!(ct > 0.0)) {
            res.log_lik = -std::numeric_limits<double>::infinity();
            return res;
        }
        log_lik += std::log(ct);
        res.grad_eta_lambda += dct_l / ct;
        res.grad_eta_p      += dct_p / ct;
        res.grad_eta_logr   += dct_r / ct;
        for (int jv = 0; jv < nIv; ++jv) {
            res.grad_eta_omega_vec[jv] += dct_o[jv] / ct;
            res.grad_eta_gamma_vec[jv] += dct_g[jv] / ct;
        }
        for (int n2 = 0; n2 < S; ++n2) {
            a[n2]   /= ct;
            da_l[n2] = (da_l[n2] - a[n2] * dct_l) / ct;
            da_p[n2] = (da_p[n2] - a[n2] * dct_p) / ct;
            da_r[n2] = (da_r[n2] - a[n2] * dct_r) / ct;
            for (int jv = 0; jv < nIv; ++jv) {
                da_o[jv][n2] = (da_o[jv][n2] - a[n2] * dct_o[jv]) / ct;
                da_g[jv][n2] = (da_g[jv][n2] - a[n2] * dct_g[jv]) / ct;
            }
        }
    }

    for (int iv = 0; iv < nIv; ++iv) {
        res.grad_eta_omega += res.grad_eta_omega_vec[iv];
        res.grad_eta_gamma += res.grad_eta_gamma_vec[iv];
    }
    res.log_lik = log_lik;
    return res;
}

// Constant-rate overload: one survival omega and recruitment gamma shared across
// a site's T-1 transitions. Broadcasts the scalar to every interval and calls the
// interval-indexed core; the scalar grad_eta_omega / grad_eta_gamma (the sum over
// intervals, which is exactly d log L / d eta for a shared rate) is what callers
// read, so this is byte-identical to the original constant-rate kernel.
inline DynAbunSiteResult compute_dyn_abun_site(
    const int* y, int T, int J, int K,
    double eta_lambda, double eta_p, double eta_omega, double eta_gamma,
    bool use_nb = false, double eta_logr = 0.0
) {
    const int nIv = T - 1;
    std::vector<double> eo(nIv, eta_omega), eg(nIv, eta_gamma);
    return compute_dyn_abun_site(y, T, J, K, eta_lambda, eta_p,
                                 eo.data(), eg.data(), use_nb, eta_logr);
}

struct DynAbunPCurv {
    double log_lik;
    double d1;   // d log L / d eta_p
    double d2;   // d^2 log L / d eta_p^2
};

// Per-site log marginal L(eta_p) and its first / second derivatives in the
// site-level detection offset eta_p, by a SECOND-ORDER forward-mode pass through
// the same exact HMM forward recursion as compute_dyn_abun_site. eta_p is a single
// scalar per site (the detection intercept shifts eta_p uniformly across the
// site's visits), so a grouped detection random effect adds a scalar offset to it
// -- exactly the per-row separability the AGHQ make_site contract needs. Unlike
// the initial-abundance arm (where eta_lambda enters ONLY the season-1 initial
// distribution, so the data-conditional weights c(n1) are eta-independent and
// precomputed once), eta_p enters the observation pmf at EVERY season, so the full
// O(K^2 T) forward marginal is re-evaluated per call (the per-node-cost obstacle of
// tulpaObs#82). Detection enters only the observation operator obs_t(n); the
// survival / recruitment transition and the initial distribution are p-independent,
// so the forward-mode first and second derivatives propagate through the linear
// transition unchanged and the only source terms are d obs / d eta_p and
// d^2 obs / d eta_p^2, both closed form:
//   obs_t(n) = prod_j Binom(y_tj | n, p),  p = invlogit(eta_p)
//   d  log obs_t(n) / d eta_p = sum_j (y_tj - n p) =: g_t(n)
//   d2 log obs_t(n) / d eta_p^2 = - (#visits) n p(1-p)
//   d  obs = obs g,  d2 obs = obs (g^2 + d2 log obs).
// Per season log L += log c_t, with c_t the forward normaliser; the first
// derivative is sum_t c_t' / c_t (= grad_eta_p of compute_dyn_abun_site) and the
// second is sum_t [c_t'' / c_t - (c_t' / c_t)^2]. The forward state and its first
// and second eta_p derivatives are renormalised each season exactly as the
// likelihood is, keeping the recursion stable. With want_deriv = false only log L
// is formed (the AGHQ lmat path needs the marginal alone). Poisson or NB initial
// abundance; eta_omega / eta_gamma are interval-indexed (length T-1).
inline DynAbunPCurv compute_dyn_abun_p_curv(
    const int* y, int T, int J, int K,
    double eta_lambda, double eta_p,
    const double* eta_omega, const double* eta_gamma,
    bool use_nb = false, double eta_logr = 0.0, bool want_deriv = true
) {
    const int S = K + 1;
    const int nIv = T - 1;
    const double lambda = std::exp(eta_lambda);
    const double p      = da_inv_logit(eta_p);
    const double logp = std::log(p), log1mp = std::log1p(-p);
    const double pq = p * (1.0 - p);
    const double rr = use_nb ? std::exp(eta_logr) : 0.0;

    // Per-interval survival / recruitment pmfs (p-independent).
    std::vector<double> logom(nIv), log1mom(nIv);
    std::vector<std::vector<double> > pois(nIv);
    for (int iv = 0; iv < nIv; ++iv) {
        const double omega = da_inv_logit(eta_omega[iv]);
        const double gamma = std::exp(eta_gamma[iv]);
        logom[iv] = std::log(omega); log1mom[iv] = std::log1p(-omega);
        pois[iv].resize(S);
        da_recruit_pmf(S, gamma, std::log(gamma), pois[iv]);
    }

    DynAbunPCurv res; res.log_lik = 0.0; res.d1 = 0.0; res.d2 = 0.0;

    // obs_t(n) and (want_deriv) its 1st / 2nd eta_p derivatives.
    auto obs_season = [&](int t, std::vector<double>& obs,
                          std::vector<double>& dobs, std::vector<double>& d2obs) {
        int ysum = 0, nv = 0;
        for (int j = 0; j < J; ++j) {
            const int yy = y[t * J + j];
            if (yy < 0) continue;
            nv++; ysum += yy;
        }
        da_obs_season_pmf(y, t, J, S, logp, log1mp, obs);
        if (want_deriv) {
            for (int n = 0; n < S; ++n) {
                const double g = (double)ysum - (double)nv * (double)n * p;
                dobs[n]  = obs[n] * g;
                d2obs[n] = obs[n] * (g * g - (double)nv * (double)n * pq);
            }
        }
    };

    std::vector<double> a(S), da(S), d2a(S), obs(S), dobs(S), d2obs(S);

    // --- Season 1: initial (p-independent) x observation. ---
    obs_season(0, obs, dobs, d2obs);
    double c = 0.0, dc = 0.0, d2c = 0.0;
    for (int n = 0; n < S; ++n) {
        double pi_n;
        if (use_nb) {
            const double rpm = rr + lambda;
            const double lpn = R::lgammafn((double)n + rr) - R::lgammafn(rr)
                - R::lgammafn((double)n + 1.0)
                + rr * std::log(rr / rpm) + (double)n * std::log(lambda / rpm);
            pi_n = std::exp(lpn);
        } else {
            const double lpn = -lambda + (double)n * eta_lambda - R::lgammafn((double)n + 1.0);
            pi_n = std::exp(lpn);
        }
        a[n] = pi_n * obs[n]; c += a[n];
        if (want_deriv) {
            da[n]  = pi_n * dobs[n];  d2a[n] = pi_n * d2obs[n];
            dc += da[n]; d2c += d2a[n];
        }
    }
    if (!(c > 0.0)) {
        res.log_lik = -std::numeric_limits<double>::infinity();
        return res;
    }
    res.log_lik = std::log(c);
    if (want_deriv) {
        const double r1 = dc / c, r2 = d2c / c;
        res.d1 += r1; res.d2 += r2 - r1 * r1;
        // Renormalise alpha and its first / second derivatives:
        //   ahat = a/c,  ahat' = a'/c - ahat (c'/c),
        //   ahat'' = a''/c - 2 ahat' (c'/c) - ahat (c''/c).
        for (int n = 0; n < S; ++n) {
            const double ah  = a[n] / c;
            const double dah = da[n] / c - ah * r1;
            const double d2ah = d2a[n] / c - 2.0 * dah * r1 - ah * r2;
            a[n] = ah; da[n] = dah; d2a[n] = d2ah;
        }
    } else {
        for (int n = 0; n < S; ++n) a[n] /= c;
    }

    // --- Seasons 2..T: transition (p-independent) then observation. ---
    std::vector<double> pre(S), dpre(S), d2pre(S), binom(S);
    for (int t = 1; t < T; ++t) {
        const int iv = t - 1;
        const std::vector<double>& pois_iv = pois[iv];
        for (int n2 = 0; n2 < S; ++n2) {
            pre[n2] = 0.0;
            if (want_deriv) { dpre[n2] = 0.0; d2pre[n2] = 0.0; }
        }
        for (int n = 0; n < S; ++n) {
            da_binom_pmf_row(n, logom[iv], log1mom[iv], binom);
            const double an = a[n];
            const double dan  = want_deriv ? da[n]  : 0.0;
            const double d2an = want_deriv ? d2a[n] : 0.0;
            for (int s = 0; s <= n; ++s) {
                const double bs = binom[s];
                for (int gn = 0; gn + s < S; ++gn) {
                    const int n2 = s + gn;
                    const double tr = bs * pois_iv[gn];
                    pre[n2] += an * tr;
                    if (want_deriv) { dpre[n2] += dan * tr; d2pre[n2] += d2an * tr; }
                }
            }
        }
        obs_season(t, obs, dobs, d2obs);
        double ct = 0.0, dct = 0.0, d2ct = 0.0;
        for (int n2 = 0; n2 < S; ++n2) {
            a[n2] = pre[n2] * obs[n2]; ct += a[n2];
            if (want_deriv) {
                da[n2]  = dpre[n2] * obs[n2] + pre[n2] * dobs[n2];
                d2a[n2] = d2pre[n2] * obs[n2] + 2.0 * dpre[n2] * dobs[n2]
                          + pre[n2] * d2obs[n2];
                dct += da[n2]; d2ct += d2a[n2];
            }
        }
        if (!(ct > 0.0)) {
            res.log_lik = -std::numeric_limits<double>::infinity();
            res.d1 = 0.0; res.d2 = 0.0; return res;
        }
        res.log_lik += std::log(ct);
        if (want_deriv) {
            const double r1 = dct / ct, r2 = d2ct / ct;
            res.d1 += r1; res.d2 += r2 - r1 * r1;
            for (int n2 = 0; n2 < S; ++n2) {
                const double ah  = a[n2] / ct;
                const double dah = da[n2] / ct - ah * r1;
                const double d2ah = d2a[n2] / ct - 2.0 * dah * r1 - ah * r2;
                a[n2] = ah; da[n2] = dah; d2a[n2] = d2ah;
            }
        } else {
            for (int n2 = 0; n2 < S; ++n2) a[n2] /= ct;
        }
    }
    return res;
}

// Constant-rate overload: one survival / recruitment shared across the site's
// T-1 transitions (broadcasts to the interval-indexed core above).
inline DynAbunPCurv compute_dyn_abun_p_curv(
    const int* y, int T, int J, int K,
    double eta_lambda, double eta_p, double eta_omega, double eta_gamma,
    bool use_nb = false, double eta_logr = 0.0, bool want_deriv = true
) {
    const int nIv = T - 1;
    std::vector<double> eo(nIv, eta_omega), eg(nIv, eta_gamma);
    return compute_dyn_abun_p_curv(y, T, J, K, eta_lambda, eta_p,
                                   eo.data(), eg.data(), use_nb, eta_logr, want_deriv);
}

// Per-site conditional likelihood c(n1) = P(y_1, ..., y_T | N_1 = n1), the data
// likelihood given the season-1 abundance, INDEPENDENT of the initial-abundance
// predictor eta_lambda (it conditions on N_1). This is the workhorse of the
// grouped random-effect AGHQ path on the initial-abundance arm (tulpaObs#51): a
// site-level RE shifts only eta_lambda, which enters solely the season-1 initial
// distribution pi_{n1}(eta_lambda), so the per-site marginal is
//   L(eta_lambda) = sum_{n1} pi_{n1}(eta_lambda) c(n1),
// and its first / second eta_lambda derivatives are sum_{n1} pi'_{n1} c(n1) and
// sum_{n1} pi''_{n1} c(n1) -- O(K) dot products once c is known. The expensive
// O(K^2 T) work (the survival / recruitment transition and the per-season
// observation) lives entirely in c, which the engine precomputes ONCE per
// make_site call (the detection / survival / recruitment predictors are held
// fixed during the RE integration); each quadrature-node / per-group-Newton
// evaluation is then a cheap dot product. c is obtained by the HMM backward
// recursion b_t(n) = sum_{n'} Tr_t(n -> n') obs_{t+1}(n') b_{t+1}(n'), with
// c(n1) = obs_1(n1) b_1(n1). Tr_t uses the interval-t survival / recruitment, so
// eta_omega / eta_gamma point to length-(T-1) arrays. Writes c_out[0..K].
inline void compute_dyn_abun_init_weights(
    const int* y, int T, int J, int K,
    double eta_p, const double* eta_omega, const double* eta_gamma, double* c_out
) {
    const int S = K + 1;
    const int nIv = T - 1;
    const double p     = da_inv_logit(eta_p);
    const double logp = std::log(p), log1mp = std::log1p(-p);

    std::vector<std::vector<double> > pois(nIv);
    std::vector<double> logom(nIv), log1mom(nIv);
    for (int iv = 0; iv < nIv; ++iv) {
        const double omega = da_inv_logit(eta_omega[iv]);
        const double gamma = std::exp(eta_gamma[iv]);
        logom[iv] = std::log(omega); log1mom[iv] = std::log1p(-omega);
        pois[iv].resize(S);
        da_recruit_pmf(S, gamma, std::log(gamma), pois[iv]);
    }

    // backward: b_{T-1}(n) = 1 (no future observations); for t = T-2 .. 0,
    // b_t(n) = sum_{n'} Tr_t(n -> n') obs_{t+1}(n') b_{t+1}(n'). The transition
    // from season t to t+1 uses interval t.
    std::vector<double> b(S, 1.0), bprev(S), obs(S), w(S), binom(S);
    for (int t = T - 2; t >= 0; --t) {
        const std::vector<double>& pois_t = pois[t];
        da_obs_season_pmf(y, t + 1, J, S, logp, log1mp, obs);
        for (int np = 0; np < S; ++np) w[np] = obs[np] * b[np];
        for (int n = 0; n < S; ++n) {
            da_binom_pmf_row(n, logom[t], log1mom[t], binom);
            double acc = 0.0;
            for (int s = 0; s <= n; ++s) {
                const double bs = binom[s];
                for (int gn = 0; gn + s < S; ++gn) acc += bs * pois_t[gn] * w[s + gn];
            }
            bprev[n] = acc;
        }
        b.swap(bprev);
    }
    da_obs_season_pmf(y, 0, J, S, logp, log1mp, obs);
    for (int n = 0; n < S; ++n) c_out[n] = obs[n] * b[n];
}

// Constant-rate overload: one survival / recruitment shared across the site's
// T-1 transitions (broadcasts to the interval-indexed core).
inline void compute_dyn_abun_init_weights(
    const int* y, int T, int J, int K,
    double eta_p, double eta_omega, double eta_gamma, double* c_out
) {
    const int nIv = T - 1;
    std::vector<double> eo(nIv, eta_omega), eg(nIv, eta_gamma);
    compute_dyn_abun_init_weights(y, T, J, K, eta_p, eo.data(), eg.data(), c_out);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_DYN_ABUN_KERNEL_H
