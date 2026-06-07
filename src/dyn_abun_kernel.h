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

// `use_nb` switches the season-1 initial abundance from Poisson(lambda) to
// negative-binomial NB(mean = lambda, size = r), r = exp(eta_logr) -- the
// pcountOpen "NB" mixture. Only the initial distribution changes; survival,
// recruitment and detection are identical. The dispersion enters the marginal
// solely through season 1, so its derivative `da_r` propagates exactly like
// `da_l` (lambda), which is also initial-only. With `use_nb = false` the
// Poisson path is byte-identical to before (da_r stays 0, grad_eta_logr = 0).
inline DynAbunSiteResult compute_dyn_abun_site(
    const int* y, int T, int J, int K,
    double eta_lambda, double eta_p, double eta_omega, double eta_gamma,
    bool use_nb = false, double eta_logr = 0.0
) {
    const int S = K + 1;                       // number of abundance states
    const double lambda = std::exp(eta_lambda);
    const double p      = da_inv_logit(eta_p);
    const double omega  = da_inv_logit(eta_omega);
    const double gamma  = std::exp(eta_gamma);
    const double logp = std::log(p), log1mp = std::log1p(-p);
    const double logom = std::log(omega), log1mom = std::log1p(-omega);
    const double loggam = std::log(gamma);
    const double rr = use_nb ? std::exp(eta_logr) : 0.0;  // NB size

    DynAbunSiteResult res;
    res.grad_eta_lambda = res.grad_eta_p = res.grad_eta_omega = res.grad_eta_gamma = 0.0;
    res.grad_eta_logr = 0.0;

    // Recruitment pmf and its eta_gamma derivative dpois_g = pois * (g - gamma).
    std::vector<double> pois(S), dpois_g(S);
    da_recruit_pmf(S, gamma, loggam, pois);
    for (int g = 0; g < S; ++g) dpois_g[g] = pois[g] * ((double)g - gamma);

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
    // is the dispersion derivative, carried alongside da_l (both initial-only).
    std::vector<double> a(S), da_l(S), da_p(S), da_o(S), da_g(S), da_r(S);
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
        da_o[n] = 0.0; da_g[n] = 0.0;
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
    // Normalise alpha and its derivatives.
    for (int n = 0; n < S; ++n) {
        a[n]   /= c1;
        da_l[n] = (da_l[n] - a[n] * dc_l) / c1;
        da_p[n] = (da_p[n] - a[n] * dc_p) / c1;
        da_r[n] = (da_r[n] - a[n] * dc_r) / c1;
        da_o[n] = 0.0; da_g[n] = 0.0;       // no omega/gamma dependence yet
    }

    // --- Seasons 2..T: transition then observation. ---
    std::vector<double> pre(S), dpre_l(S), dpre_p(S), dpre_o(S), dpre_g(S), dpre_r(S);
    std::vector<double> binom(S), dbinom(S);   // survivor pmf for the current row n
    for (int t = 1; t < T; ++t) {
        for (int n2 = 0; n2 < S; ++n2) {
            pre[n2] = dpre_l[n2] = dpre_p[n2] = dpre_o[n2] = dpre_g[n2] = dpre_r[n2] = 0.0;
        }
        for (int n = 0; n < S; ++n) {
            if (a[n] == 0.0 && da_l[n] == 0.0 && da_p[n] == 0.0 &&
                da_o[n] == 0.0 && da_g[n] == 0.0 && da_r[n] == 0.0) continue;
            // Survivor pmf Binom(s | n, omega) and d/d eta_omega = binom*(s - n omega).
            da_binom_pmf_row(n, logom, log1mom, binom);
            for (int s = 0; s <= n; ++s) dbinom[s] = binom[s] * ((double)s - (double)n * omega);
            const double an = a[n], dl = da_l[n], dp = da_p[n], dom = da_o[n], dg = da_g[n];
            const double dr = da_r[n];
            // Tr(n -> n') = sum_s binom[s] pois[n'-s]; scatter into pre.
            for (int s = 0; s <= n; ++s) {
                const double bs = binom[s], dbs = dbinom[s];
                for (int gn = 0; gn + s < S; ++gn) {
                    const int n2 = s + gn;
                    const double tr = bs * pois[gn];
                    pre[n2]   += an * tr;
                    // chain rule: d(alpha_{t-1}(n) Tr)/d eta_k. The transition
                    // carries no lambda/r dependence, so those propagate the
                    // upstream derivative alone (like lambda).
                    dpre_l[n2] += dl * tr;
                    dpre_p[n2] += dp * tr;
                    dpre_o[n2] += dom * tr + an * (dbs * pois[gn]);
                    dpre_g[n2] += dg * tr + an * (bs * dpois_g[gn]);
                    dpre_r[n2] += dr * tr;
                }
            }
        }
        obs_season(t, obs, dobs);
        double ct = 0.0, dct_l = 0.0, dct_p = 0.0, dct_o = 0.0, dct_g = 0.0, dct_r = 0.0;
        for (int n2 = 0; n2 < S; ++n2) {
            a[n2]   = pre[n2] * obs[n2];
            da_l[n2] = dpre_l[n2] * obs[n2];
            da_p[n2] = dpre_p[n2] * obs[n2] + pre[n2] * dobs[n2];
            da_o[n2] = dpre_o[n2] * obs[n2];
            da_g[n2] = dpre_g[n2] * obs[n2];
            da_r[n2] = dpre_r[n2] * obs[n2];
            ct += a[n2]; dct_l += da_l[n2]; dct_p += da_p[n2];
            dct_o += da_o[n2]; dct_g += da_g[n2]; dct_r += da_r[n2];
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
        res.grad_eta_logr   += dct_r / ct;
        for (int n2 = 0; n2 < S; ++n2) {
            a[n2]   /= ct;
            da_l[n2] = (da_l[n2] - a[n2] * dct_l) / ct;
            da_p[n2] = (da_p[n2] - a[n2] * dct_p) / ct;
            da_o[n2] = (da_o[n2] - a[n2] * dct_o) / ct;
            da_g[n2] = (da_g[n2] - a[n2] * dct_g) / ct;
            da_r[n2] = (da_r[n2] - a[n2] * dct_r) / ct;
        }
    }

    res.log_lik = log_lik;
    return res;
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
// recursion b_t(n) = sum_{n'} Tr(n -> n') obs_{t+1}(n') b_{t+1}(n'), with
// c(n1) = obs_1(n1) b_1(n1). Writes c_out[0..K].
inline void compute_dyn_abun_init_weights(
    const int* y, int T, int J, int K,
    double eta_p, double eta_omega, double eta_gamma, double* c_out
) {
    const int S = K + 1;
    const double p     = da_inv_logit(eta_p);
    const double omega = da_inv_logit(eta_omega);
    const double gamma = std::exp(eta_gamma);
    const double logp = std::log(p), log1mp = std::log1p(-p);
    const double logom = std::log(omega), log1mom = std::log1p(-omega);
    const double loggam = std::log(gamma);

    std::vector<double> pois(S);
    da_recruit_pmf(S, gamma, loggam, pois);

    // backward: b_{T-1}(n) = 1 (no future observations); for t = T-2 .. 0,
    // b_t(n) = sum_{n'} Tr(n -> n') obs_{t+1}(n') b_{t+1}(n').
    std::vector<double> b(S, 1.0), bprev(S), obs(S), w(S), binom(S);
    for (int t = T - 2; t >= 0; --t) {
        da_obs_season_pmf(y, t + 1, J, S, logp, log1mp, obs);
        for (int np = 0; np < S; ++np) w[np] = obs[np] * b[np];
        for (int n = 0; n < S; ++n) {
            da_binom_pmf_row(n, logom, log1mom, binom);
            double acc = 0.0;
            for (int s = 0; s <= n; ++s) {
                const double bs = binom[s];
                for (int gn = 0; gn + s < S; ++gn) acc += bs * pois[gn] * w[s + gn];
            }
            bprev[n] = acc;
        }
        b.swap(bprev);
    }
    da_obs_season_pmf(y, 0, J, S, logp, log1mp, obs);
    for (int n = 0; n < S; ++n) c_out[n] = obs[n] * b[n];
}

}  // namespace tulpaObs

#endif  // TULPAOBS_DYN_ABUN_KERNEL_H
