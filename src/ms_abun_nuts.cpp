// ms_abun_nuts.cpp
// C++ joint log-posterior + gradient for the community / multispecies N-mixture
// (ms_abun()) NUTS target. The R reference .tobs_ms_abun_nuts_logpost
// (R/ms_abun_nuts.R) is the oracle; this port mirrors it and is cross-checked
// against it before driving tulpa's NUTS engine.
//
// The non-spatial community N-mixture is per-species Royle (2004) with Gaussian
// community hyperpriors on the per-species coefficients:
//   N_{s,i}       ~ Poisson(lambda_{s,i})            (or NegBin(lambda, r_s))
//   y_{s,i,j}|N   ~ Binomial(N_{s,i}, p_{s,i,j})
//   log lambda_{s,i} = X_lambda_i  . (mu_lambda + b_lambda_s)
//   logit p_{s,i,j}  = X_p_{ij}    . (mu_p      + b_p_s)
//   b_lambda_s ~ N(0, Sigma_lambda),  b_p_s ~ N(0, Sigma_p)
// Under NB the dispersion is a per-species random effect log_r_s ~
// N(mu_log_r, sigma_log_r^2), with r_s = exp(mu_log_r + b_logr_s) constant across
// that species' sites. N_{s,i} integrates out per species-site in closed form;
// the Laplace-EM (nmix_laplace_re) profiles the per-species deviations and
// community covariances out, returning a Gaussian community-mean posterior. NUTS
// instead samples EVERYTHING jointly -- the community means, the per-species
// deviations, AND the community covariances -- from the exact joint posterior,
// which gives calibrated (non-Gaussian) community-mean / covariance intervals and
// the per-(species, site) pointwise likelihood WAIC / LOO need.
//
// NON-CENTERED parameterisation: the per-species block carries the whitened
// standard-normal z_s, and the deviation is reconstructed per arm as
// b_{s,arm} = C_arm z_{s,arm} (C_arm the log-Cholesky factor of Sigma_arm). The
// community covariance then leaves the b-prior entirely and enters ONLY the data
// term through b = C z. This breaks the centered b<->Sigma funnel that otherwise
// saturated the NUTS treedepth. The target factorises as
//   log p = sum_{s,i} log m_{s,i}(theta)                 # per-species-site marginal
//         - 0.5 ||mu_coef||^2 / sigma.beta^2             # community-mean priors
//         [ - 0.5 mu_log_r^2 / sigma.logr^2 ]            # (NB)
//         - 0.5 sum_s ||z_s||^2                          # whitened RE prior (N(0,I))
//         + log p(Sigma coords)                          # log-Cholesky hyperpriors
// where m_{s,i} is the Royle marginal exposed by compute_nmix_site()
// (nmix_kernel.h, the same kernel the Laplace fit and the AGHQ RE path use); it
// returns grad_eta_lambda, grad_eta_p, and (NB) grad_theta = d log m / d log_r, so
// the coefficient gradient is the design-sandwiched eta-gradient -- no new
// likelihood math. The chol gradient flows from the data term via b = C z
// (chol_data_grad_noncentered) plus the coordinate hyperprior; the arm Cholesky
// factors use the shared helpers in community_chol.h. The R reference
// .tobs_ms_abun_nuts_logpost mirrors this and is the cross-check oracle.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <tulpa/nuts_api.h>
#include "nmix_kernel.h"
#include "community_chol.h"

#include "nuts_engine.h"
using namespace Rcpp;

namespace tulpaObs {

// Marshalled per-fit data for the community N-mixture NUTS target. The long form
// is what .tobs_ms_nmix_longform() produces: y / site_idx / species_idx in the
// stacked per-(species, site, visit) order, X_lambda site-level (n_sites rows,
// shared across species), X_p long form (n_obs rows). Per-(species, site) visit
// row lists are built once at construction so the eval loop gathers each
// species-site's counts / detection etas contiguously for compute_nmix_site().
struct MsNmixNutsData {
    int n_sites = 0, n_obs = 0, n_species = 0, p_lam = 0, p_p = 0, K_max = 0;

    // Threads for the per-species gradient loop below. 0 leaves the count to
    // OpenMP, which is what an omitted `n_threads` means; a positive value is
    // the ceiling the caller asked for. The reduction after the loop is
    // serial and order-fixed, so the gradient is the same at any count.
    int n_threads = 0;
    bool is_nb = false;
    std::vector<int> y;                                 // length n_obs (long form)
    NumericMatrix X_lambda;                             // n_sites x p_lam
    NumericMatrix X_p;                                  // n_obs   x p_p
    std::vector<std::vector<std::vector<int>>> obs;     // [s][i] -> rows into y / X_p
    // Per-(species, site) lgamma cache. The combinatorial lgamma terms of the
    // Royle marginal are eta-independent, so they are built once here (NOT on
    // every leapfrog step, which dominated the runtime); the eval loop calls
    // compute_nmix_site_cached with the precomputed cache.
    std::vector<std::vector<NMixSiteCache>> cache;      // [s][i]
    // Optional per-(species, site) latent-N ceiling, row-major [s * n_sites + i]
    // (empty = every site truncates at the shared K_max). The latent sum is
    // re-evaluated on every leapfrog step, so its width is the per-step cost;
    // keying each species-site's ceiling to that cell's own fitted abundance
    // stops one heavy cell from widening all of them.
    std::vector<int> K_site;

    // Optional fixed-hyper areal field SHARED across species on the abundance
    // arm: the non-centered Gaussian field f = Linv %*% raw, raw ~ N(0, I), with
    // Linv the inverse Cholesky of the FIXED field precision tau Q(rho) (the
    // nested-Laplace #12 estimate). eta_lambda_{s,i} += f[u(i)] for every
    // species. The flat vector grows by `raw` (length n_field_units) appended
    // after the chol blocks. The field is f = L %*% raw, raw ~ N(0, I_{n_raw}),
    // with L a (possibly non- square) n_field_units x n_raw loading: the square
    // inverse Cholesky for a full-rank proper-CAR field (n_raw == n_field_units),
    // or the sum-to-zero eigen-loading for an intrinsic icar (n_raw =
    // n_field_units - 1) / bym2 (n_raw = 2 n_field_units - 1) field.
    int n_field_units = 0, n_raw = 0, o_raw = 0;
    std::vector<int> field_map;       // 0-based unit per site (length n_sites)
    std::vector<double> Linv;         // row-major n_field_units x n_raw loading

    // Packed-coordinate layout (mirrors .tobs_ms_abun_nuts_layout).
    int P_tot = 0, q_lam = 0, q_p = 0, q_logr = 0, total = 0;
    int mu_off = 0, b_off = 0, chol_lam_off = 0, chol_p_off = 0, chol_logr_off = 0;
};

inline void ms_abun_nuts_layout(MsNmixNutsData& d) {
    d.P_tot   = d.p_lam + d.p_p + (d.is_nb ? 1 : 0);
    d.q_lam   = d.p_lam * (d.p_lam + 1) / 2;
    d.q_p     = d.p_p   * (d.p_p   + 1) / 2;
    d.q_logr  = d.is_nb ? 1 : 0;
    d.mu_off  = 0;
    d.b_off   = d.P_tot;
    d.chol_lam_off  = d.P_tot + d.n_species * d.P_tot;
    d.chol_p_off    = d.chol_lam_off + d.q_lam;
    d.chol_logr_off = d.chol_p_off + d.q_p;
    d.o_raw   = d.chol_logr_off + d.q_logr;
    d.total   = d.o_raw + d.n_raw;
}

inline MsNmixNutsData ms_abun_nuts_build_data(const Rcpp::List& spec) {
    MsNmixNutsData d;
    d.n_threads = spec.containsElementNamed("n_threads")
                  ? Rcpp::as<int>(spec["n_threads"]) : 0;
    IntegerVector y        = spec["y"];
    IntegerVector site_idx = spec["site_idx"];        // 1-based, into X_lambda rows
    IntegerVector sp_idx   = spec["species_idx"];     // 1-based
    d.X_lambda  = Rcpp::as<NumericMatrix>(spec["X_lambda"]);
    d.X_p       = Rcpp::as<NumericMatrix>(spec["X_p"]);
    d.n_sites   = Rcpp::as<int>(spec["n_sites"]);
    d.n_species = Rcpp::as<int>(spec["n_species"]);
    d.K_max     = Rcpp::as<int>(spec["K_max"]);
    d.is_nb     = Rcpp::as<bool>(spec["is_nb"]);
    d.n_obs     = y.size();
    d.p_lam     = d.X_lambda.ncol();
    d.p_p       = d.X_p.ncol();
    d.y.assign(y.begin(), y.end());
    if (spec.containsElementNamed("K_site")) {
        IntegerMatrix Ks = Rcpp::as<IntegerMatrix>(spec["K_site"]);
        if (Ks.nrow() != d.n_species || Ks.ncol() != d.n_sites) {
            Rcpp::stop("K_site must be n_species x n_sites (%d x %d, got %d x %d).",
                       d.n_species, d.n_sites, Ks.nrow(), Ks.ncol());
        }
        d.K_site.assign((size_t)d.n_species * d.n_sites, 0);
        for (int s = 0; s < d.n_species; ++s)
            for (int i = 0; i < d.n_sites; ++i)
                d.K_site[(size_t)s * d.n_sites + i] = Ks(s, i);
    }
    // Optional shared areal field block.
    if (spec.containsElementNamed("n_field_units")) {
        d.n_field_units = Rcpp::as<int>(spec["n_field_units"]);
        if (d.n_field_units > 0) {
            IntegerVector fmap = spec["field_map"];     // 1-based unit per site
            if ((int) fmap.size() != d.n_sites)
                Rcpp::stop("field_map must have one entry per site.");
            d.field_map.assign(d.n_sites, 0);
            for (int i = 0; i < d.n_sites; ++i) {
                const int u = fmap[i] - 1;
                if (u < 0 || u >= d.n_field_units)
                    Rcpp::stop("field_map out of range [1, n_field_units].");
                d.field_map[i] = u;
            }
            // Accept a general n_field_units x n_raw loading (field_load, #113);
            // the legacy square inverse Cholesky (field_Linv) is n_raw == NF.
            NumericMatrix Li = spec.containsElementNamed("field_load")
                ? Rcpp::as<NumericMatrix>(spec["field_load"])
                : Rcpp::as<NumericMatrix>(spec["field_Linv"]);
            if (Li.nrow() != d.n_field_units)
                Rcpp::stop("field loading must have n_field_units rows.");
            d.n_raw = Li.ncol();
            d.Linv.assign((std::size_t) d.n_field_units * d.n_raw, 0.0);
            for (int a = 0; a < d.n_field_units; ++a)
                for (int b = 0; b < d.n_raw; ++b)
                    d.Linv[(std::size_t) a * d.n_raw + b] = Li(a, b);
        }
    }
    d.obs.assign(d.n_species,
                 std::vector<std::vector<int>>(d.n_sites, std::vector<int>()));
    for (int o = 0; o < d.n_obs; ++o) {
        const int s = sp_idx[o] - 1, i = site_idx[o] - 1;
        if (s < 0 || s >= d.n_species)
            Rcpp::stop("species_idx out of range in ms_abun_nuts_build_data");
        if (i < 0 || i >= d.n_sites)
            Rcpp::stop("site_idx out of range in ms_abun_nuts_build_data");
        d.obs[s][i].push_back(o);
    }
    // Build the eta-independent lgamma cache for every observed (species, site).
    d.cache.assign(d.n_species, std::vector<NMixSiteCache>(d.n_sites));
    std::vector<int> y_site;
    for (int s = 0; s < d.n_species; ++s) {
        for (int i = 0; i < d.n_sites; ++i) {
            const std::vector<int>& obs = d.obs[s][i];
            const int J = (int) obs.size();
            if (J == 0) continue;
            y_site.resize(J);
            for (int jj = 0; jj < J; ++jj) y_site[jj] = d.y[obs[jj]];
            const int K_si = d.K_site.empty()
                ? 0 : d.K_site[(size_t)s * d.n_sites + i];
            if (K_si > 0) {
                int y_max = 0;
                for (int jj = 0; jj < J; ++jj)
                    if (y_site[jj] > y_max) y_max = y_site[jj];
                if (K_si < y_max) {
                    Rcpp::stop("K_site[%d, %d] = %d is below that cell's max(y) = %d.",
                               s + 1, i + 1, K_si, y_max);
                }
            }
            d.cache[s][i] = nmix_precompute_site(y_site.data(), J, d.K_max, -1, K_si);
        }
    }
    ms_abun_nuts_layout(d);
    return d;
}

// Joint log-posterior + gradient over the packed vector
//   theta = (mu, {b_s} species-major, chol_lambda, chol_p [, chol_logr]).
// NUTS maximises, so this returns the log-posterior (no negation) and writes the
// gradient into `g` (length d.total). Mirrors .tobs_ms_abun_nuts_logpost.
inline double ms_abun_nuts_eval(const MsNmixNutsData& d, const double* th,
                                double sigma_beta, double sigma_logr,
                                const CommunityCholPri& pr, double* g) {
    const int P = d.P_tot, S = d.n_species, p_lam = d.p_lam, p_p = d.p_p;
    const bool nb = d.is_nb;
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;

    const double* mu = th + d.mu_off;
    const double* z  = th + d.b_off;       // species-major z (whitened), length P each
    double* g_mu = g + d.mu_off;
    double* g_z  = g + d.b_off;

    // Cholesky factors per arm (Sigma_arm = C_arm C_arm', row-major); the logr
    // arm is the 1x1 scalar SD = exp(chol_logr). Under the non-centered map the
    // per-species deviation is b_{s,arm} = C_arm z_{s,arm}, so the covariance
    // enters only the data term (z carries a standard-normal prior).
    std::vector<double> C_lam, C_p;
    chol_unpack_cpp(th + d.chol_lam_off, p_lam, C_lam);
    chol_unpack_cpp(th + d.chol_p_off,   p_p,   C_p);
    const double C_lr = nb ? chol_diag_exp(th[d.chol_logr_off]) : 0.0;

    // chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
    std::vector<double> A_lam((std::size_t) p_lam * p_lam, 0.0),
                        A_p((std::size_t) p_p * p_p, 0.0);
    double A_lr = 0.0;

    // Shared areal field: reconstruct f = Linv %*% raw once. Every species'
    // eta_lambda gains the per-site offset f[field_map[site]]; the field score
    // accumulates grad_eta_lambda over species (and over the sites mapping to
    // each unit), then maps back to raw via Linv'.
    const int NF = d.n_field_units;
    const int NR = d.n_raw;
    std::vector<double> f_field(NF, 0.0);
    if (NF > 0) {
        const double* raw = th + d.o_raw;
        for (int a = 0; a < NF; ++a) {
            double v = 0.0;
            for (int b = 0; b < NR; ++b)
                v += d.Linv[(std::size_t) a * NR + b] * raw[b];
            f_field[a] = v;
        }
    }
    std::vector<double> gf_s((std::size_t) (NF > 0 ? S : 0) * NF, 0.0);

    // ---- data log-lik + inner gradient (per species, per site) ----
    //
    // The per-species work is independent: each species reconstructs b = C z,
    // scores its sites, writes its own (disjoint) g_z block, and produces partial
    // (g_mu, A_arm, log-lik) contributions. The species loop is parallelised over
    // cores; the partials land in per-species slots and are reduced SERIALLY in
    // species order afterwards, so the reduction is deterministic (thread-count
    // independent). The per-species partial-sum reordering differs from the fully
    // interleaved serial sum only at the floating-point reduction-order level
    // (~1e-13), well inside the 1e-9 cross-check against the R oracle.
    std::vector<double> gmu_s((std::size_t) S * P, 0.0), lp_s(S, 0.0);
    std::vector<double> Alam_s((std::size_t) S * p_lam * p_lam, 0.0);
    std::vector<double> Ap_s((std::size_t) S * p_p * p_p, 0.0);
    std::vector<double> Alr_s(S, 0.0);

    // Scratch sized once per thread and reused across species. b_* are fully
    // overwritten each pass and eta_p_site is resized then written per site;
    // gbl/gbp accumulate and are re-zeroed explicitly.
    #ifdef _OPENMP
    const int omp_n = d.n_threads > 0 ? d.n_threads : omp_get_max_threads();
    #pragma omp parallel num_threads(omp_n)
    #endif
    {
    std::vector<double> b_lam(p_lam), b_p(p_p), gbl(p_lam), gbp(p_p), eta_p_site;

    #pragma omp for schedule(static)
    for (int s = 0; s < S; ++s) {
        const double* z_s = z + s * P;
        const double* zl  = z_s;             // length p_lam
        const double* zp  = z_s + p_lam;     // length p_p
        const double  zr  = nb ? z_s[p_lam + p_p] : 0.0;
        std::fill(gbl.begin(), gbl.end(), 0.0);
        std::fill(gbp.begin(), gbp.end(), 0.0);
        // reconstruct b = C z (lower-triangular C)
        for (int i = 0; i < p_lam; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j) v += C_lam[(std::size_t) i * p_lam + j] * zl[j];
            b_lam[i] = v;
        }
        for (int i = 0; i < p_p; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j) v += C_p[(std::size_t) i * p_p + j] * zp[j];
            b_p[i] = v;
        }
        const double b_lr = nb ? C_lr * zr : 0.0;
        const double r = nb ? std::exp(mu[p_lam + p_p] + b_lr)
                            : std::numeric_limits<double>::infinity();
        double* gmu_loc = &gmu_s[(std::size_t) s * P];
        double lp_loc = 0.0, gblr = 0.0;
        for (int i = 0; i < d.n_sites; ++i) {
            const std::vector<int>& obs = d.obs[s][i];
            const int J = (int) obs.size();
            if (J == 0) continue;
            double eta_lambda = 0.0;
            for (int k = 0; k < p_lam; ++k)
                eta_lambda += d.X_lambda(i, k) * (mu[k] + b_lam[k]);
            if (NF > 0) eta_lambda += f_field[d.field_map[i]];
            eta_p_site.resize(J);
            for (int jj = 0; jj < J; ++jj) {
                const int o = obs[jj];
                double e = 0.0;
                for (int k = 0; k < p_p; ++k)
                    e += d.X_p(o, k) * (mu[p_lam + k] + b_p[k]);
                eta_p_site[jj] = e;
            }
            const NMixSiteResult res = compute_nmix_site_cached(
                d.cache[s][i], eta_p_site.data(), eta_lambda, r);
            lp_loc += res.log_lik;
            // abundance arm: mu_lambda and b_lambda_s share grad_eta_lambda
            for (int k = 0; k < p_lam; ++k) {
                const double gx = res.grad_eta_lambda * d.X_lambda(i, k);
                gmu_loc[k] += gx;
                gbl[k]     += gx;
            }
            // shared field: this site's grad_eta_lambda accumulates onto its unit.
            if (NF > 0)
                gf_s[(std::size_t) s * NF + d.field_map[i]] += res.grad_eta_lambda;
            // detection arm: mu_p and b_p_s share the per-visit grad_eta_p
            for (int jj = 0; jj < J; ++jj) {
                const int o = obs[jj];
                const double ge = res.grad_eta_p[jj];
                for (int k = 0; k < p_p; ++k) {
                    const double gx = ge * d.X_p(o, k);
                    gmu_loc[p_lam + k] += gx;
                    gbp[k]             += gx;
                }
            }
            // dispersion arm: mu_log_r and b_logr_s share grad_theta (theta = log r)
            if (nb) {
                gmu_loc[p_lam + p_p] += res.grad_theta;
                gblr                 += res.grad_theta;
            }
        }
        lp_s[s] = lp_loc;
        // z gradient (data part) = C' grad_b -- disjoint per-species write.
        double* gz_s = g_z + s * P;
        for (int v = 0; v < p_lam; ++v) {
            double sg = 0.0;
            for (int i = v; i < p_lam; ++i) sg += C_lam[(std::size_t) i * p_lam + v] * gbl[i];
            gz_s[v] += sg;
        }
        for (int v = 0; v < p_p; ++v) {
            double sg = 0.0;
            for (int i = v; i < p_p; ++i) sg += C_p[(std::size_t) i * p_p + v] * gbp[i];
            gz_s[p_lam + v] += sg;
        }
        if (nb) gz_s[p_lam + p_p] += C_lr * gblr;
        // per-species chol accumulators A_arm = grad_b z' (reduced in order below).
        double* Al = &Alam_s[(std::size_t) s * p_lam * p_lam];
        for (int i = 0; i < p_lam; ++i)
            for (int j = 0; j <= i; ++j) Al[(std::size_t) i * p_lam + j] = gbl[i] * zl[j];
        double* Ap = &Ap_s[(std::size_t) s * p_p * p_p];
        for (int i = 0; i < p_p; ++i)
            for (int j = 0; j <= i; ++j) Ap[(std::size_t) i * p_p + j] = gbp[i] * zp[j];
        if (nb) Alr_s[s] = gblr * zr;
    }
    }  // omp parallel

    // serial reduction in species order -> byte-identical to the serial path.
    double lp = 0.0;
    for (int s = 0; s < S; ++s) {
        const double* gmu_loc = &gmu_s[(std::size_t) s * P];
        for (int k = 0; k < P; ++k) g_mu[k] += gmu_loc[k];
        lp += lp_s[s];
        const double* Al = &Alam_s[(std::size_t) s * p_lam * p_lam];
        for (int t = 0; t < p_lam * p_lam; ++t) A_lam[t] += Al[t];
        const double* Ap = &Ap_s[(std::size_t) s * p_p * p_p];
        for (int t = 0; t < p_p * p_p; ++t) A_p[t] += Ap[t];
        if (nb) A_lr += Alr_s[s];
    }

    // ---- z prior: standard normal over the whole per-species block ----
    for (int j = 0; j < S * P; ++j) {
        const double zz = z[j];
        g_z[j] -= zz;
        lp     += -0.5 * zz * zz;
    }

    // ---- shared field: reduce per-species field scores (unit order), map to raw
    //      via Linv', add the whitened-raw N(0, I) prior. ----
    if (NF > 0) {
        std::vector<double> g_f(NF, 0.0);
        for (int s = 0; s < S; ++s)
            for (int u = 0; u < NF; ++u)
                g_f[u] += gf_s[(std::size_t) s * NF + u];
        const double* raw = th + d.o_raw;
        double* g_raw = g + d.o_raw;
        // d lp / d raw_b = sum_a L[a, b] g_f[a] - raw_b  (raw ~ N(0, I_{n_raw})).
        for (int b = 0; b < NR; ++b) {
            double v = 0.0;
            for (int a = 0; a < NF; ++a)
                v += d.Linv[(std::size_t) a * NR + b] * g_f[a];
            g_raw[b] += v - raw[b];
            lp += -0.5 * raw[b] * raw[b];
        }
    }

    // ---- chol coords: data gradient (via b = C z) + coordinate hyperprior ----
    lp += chol_data_grad_noncentered(A_lam, C_lam, p_lam, th + d.chol_lam_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_lam_off);
    lp += chol_data_grad_noncentered(A_p, C_p, p_p, th + d.chol_p_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_p_off);
    if (nb) {
        const std::vector<double> A1(1, A_lr), C1(1, C_lr);
        lp += chol_data_grad_noncentered(A1, C1, 1, th + d.chol_logr_off,
                                         pr.logdiag_mean, pr.logdiag_sd,
                                         pr.offdiag_sd, g + d.chol_logr_off);
    }

    // ---- community-mean priors: N(0, sigma.beta^2) on the coefficient means,
    //      N(0, sigma.logr^2) on the log-dispersion community mean (NB). ----
    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < p_lam + p_p; ++k) {
        g_mu[k] -= inv_sb2 * mu[k];
        lp      += -0.5 * inv_sb2 * mu[k] * mu[k];
    }
    if (nb) {
        const double inv_sl2 = 1.0 / (sigma_logr * sigma_logr);
        const double m = mu[p_lam + p_p];
        g_mu[p_lam + p_p] -= inv_sl2 * m;
        lp                += -0.5 * inv_sl2 * m * m;
    }
    return lp;
}

// NUTS model carrying the marshalled data + prior scales; the FullGradFn reaches
// it through ModelData.model_response_data.
struct MsNmixNutsModel {
    MsNmixNutsData d;
    double sigma_beta = 10.0, sigma_logr = 1.5;
    CommunityCholPri pr;
};

// FullGradFn: log-posterior + gradient over the entire parameter vector.
inline void ms_abun_nuts_full_grad(const std::vector<double>& params,
                                   const tulpa::ModelData& data,
                                   const tulpa::ParamLayout& /*layout*/,
                                   std::vector<double>& grad, double* log_post_out) {
    const MsNmixNutsModel* m =
        static_cast<const MsNmixNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = ms_abun_nuts_eval(m->d, params.data(), m->sigma_beta,
                                        m->sigma_logr, m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

} // namespace tulpaObs

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .tobs_ms_abun_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_ms_abun_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                          Rcpp::List pri, double sigma_beta,
                                          double sigma_logr) {
    tulpaObs::MsNmixNutsData d = tulpaObs::ms_abun_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    tulpaObs::CommunityCholPri pr = tulpaObs::community_chol_pri_from_list(pri);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::ms_abun_nuts_eval(
        d, theta.begin(), sigma_beta, sigma_logr, pr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the community N-mixture target via tulpa's engine and the
// FullGradFn (gradient mode "H"). `theta0` is the warm-start (the Laplace-EM
// mode); `inv_metric` an optional length-n_params inverse-mass diagonal (the
// Laplace curvature). Returns draws + diagnostics.
// [[Rcpp::export]]
Rcpp::List cpp_ms_abun_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                            Rcpp::List pri, double sigma_beta, double sigma_logr,
                            Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                            int n_iter, int n_warmup, int max_treedepth,
                            double adapt_delta, int seed, bool verbose) {
    tulpaObs::MsNmixNutsModel m;
    m.d = tulpaObs::ms_abun_nuts_build_data(spec);
    m.sigma_beta = sigma_beta; m.sigma_logr = sigma_logr;
    m.pr = tulpaObs::community_chol_pri_from_list(pri);
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    return tulpaObs::run_tulpa_nuts(
        &tulpaObs::ms_abun_nuts_full_grad, &m, m.d.total,
        theta0, sigma_beta, inv_metric,
        n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose,
        "ms_abun_nmix", m.d.n_sites);
}
