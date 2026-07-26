// ms_dyn_occu_nuts.cpp
// C++ joint log-posterior + gradient for the community / multispecies DYNAMIC
// (multi-season, HMM) occupancy (ms_dyn_occu()) NUTS target. The R reference
// .tobs_ms_dyn_occu_nuts_logpost (R/ms_dyn_occu_nuts.R) is the oracle; this port
// mirrors it and is cross-checked against it before driving tulpa's NUTS engine.
//
// The community dynamic occupancy is per-species HMM occupancy with Gaussian
// community hyperpriors on the per-species first-season occupancy / detection
// coefficients and SHARED community colonisation / extinction transition
// coefficients:
//   z_{s,i,1}       ~ Bernoulli(psi1_{s,i})
//   z_{s,i,t}|z,..  ~ transition(gamma_i, eps_i)              (t = 2..T)
//   y_{s,i,t,j}|z=1 ~ Bernoulli(p_{s,i})
//   logit psi1_{s,i} = X_psi1_i . (mu_psi1 + b_psi1_s)
//   logit p_{s,i}    = X_p_i    . (mu_p    + b_p_s)
//   logit gamma_i / eps_i = X_gamma_i . beta_gamma / X_eps_i . beta_eps  (shared)
//   b_psi1_s ~ N(0, Sigma_psi1),  b_p_s ~ N(0, Sigma_p)
// The latent occupancy path integrates out per species by a scaled HMM forward
// filter; a forward-backward smoothing pass yields the Fisher-identity gradient
// of the marginal wrt the four site-level linear predictors (season-1 smoothed
// posterior w1, detection score, and colonisation / extinction sufficient
// statistics), identical to the analytic gradient in the stMsPGOcc field fitter
// (R/ms_dyn_occu_spatial.R, .ms_dyn_occu_fb_vec).
//
// NON-CENTERED parameterisation (mirrors ms_occu_nuts.cpp): the per-species block
// carries the whitened standard-normal z_s, the deviation is reconstructed per
// arm as b_{s,arm} = C_arm z_{s,arm}, so the community covariance leaves the
// b-prior entirely and enters ONLY the data term. The two SHARED transition arms
// (gamma, eps) carry no per-species random effect; they accumulate their gradient
// across the full species loop. The target factorises as
//   log p = sum_s log L_s(theta)                     # per-species HMM forward
//         - 0.5 ||mu||^2 / sigma.beta^2              # community-mean priors
//         - 0.5 ||global||^2 / sigma.beta^2          # shared-transition priors
//         - 0.5 sum_s ||z_s||^2                      # whitened RE prior (N(0,I))
//         + log p(Sigma coords)                      # log-Cholesky hyperpriors
// The chol gradient flows from the data term via b = C z
// (chol_data_grad_noncentered) plus the coordinate hyperprior; the arm Cholesky
// factors use the shared helpers in community_chol.h.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
#include "community_chol.h"

#include "nuts_engine.h"
using namespace Rcpp;

namespace tulpaObs {

static const double MDO_CLAMP_LO = 1e-12;
static const double MDO_CLAMP_HI = 1.0 - 1e-12;
inline double mdo_clamp(double x) {
    return x < MDO_CLAMP_LO ? MDO_CLAMP_LO : (x > MDO_CLAMP_HI ? MDO_CLAMP_HI : x);
}
inline double mdo_plogis(double e) { return 1.0 / (1.0 + std::exp(-e)); }

// Log-Cholesky coordinate hyperprior scalars for the community covariances.
struct MsDynOccuPri {
    double logdiag_mean = 0.0, logdiag_sd = 1.5, offdiag_sd = 1.0;
};

// Marshalled per-fit data. X_* are site-level designs (n_sites rows, shared
// across species). Per species the (n_valid, n_det) [n_sites x n_seasons]
// detection sufficient statistics drive the HMM emissions; stored flattened
// column-major (index i + t * n_sites) to match the R matrices.
struct MsDynOccuNutsData {
    int n_sites = 0, n_seasons = 0, n_species = 0;
    int p_psi1 = 0, p_p = 0, p_gam = 0, p_eps = 0;
    NumericMatrix X_psi1, X_p, X_gamma, X_eps;
    std::vector<std::vector<int>> nvalid;         // [s][i + t*n_sites]
    std::vector<std::vector<int>> ndet;           // [s][i + t*n_sites]

    // Packed-coordinate layout (mirrors .tobs_ms_dyn_occu_nuts_layout).
    int P_tot = 0, q_psi1 = 0, q_p = 0, G = 0, total = 0;
    int mu_off = 0, b_off = 0, chol_psi1_off = 0, chol_p_off = 0, global_off = 0;
};

inline void ms_dyn_occu_nuts_layout(MsDynOccuNutsData& d) {
    d.P_tot         = d.p_psi1 + d.p_p;
    d.q_psi1        = d.p_psi1 * (d.p_psi1 + 1) / 2;
    d.q_p           = d.p_p    * (d.p_p    + 1) / 2;
    d.G             = d.p_gam + d.p_eps;
    d.mu_off        = 0;
    d.b_off         = d.P_tot;
    d.chol_psi1_off = d.P_tot + d.n_species * d.P_tot;
    d.chol_p_off    = d.chol_psi1_off + d.q_psi1;
    d.global_off    = d.chol_p_off + d.q_p;
    d.total         = d.global_off + d.G;
}

inline MsDynOccuNutsData ms_dyn_occu_nuts_build_data(const Rcpp::List& spec) {
    MsDynOccuNutsData d;
    d.X_psi1    = Rcpp::as<NumericMatrix>(spec["X_psi1"]);
    d.X_p       = Rcpp::as<NumericMatrix>(spec["X_p"]);
    d.X_gamma   = Rcpp::as<NumericMatrix>(spec["X_gamma"]);
    d.X_eps     = Rcpp::as<NumericMatrix>(spec["X_eps"]);
    d.n_sites   = Rcpp::as<int>(spec["n_sites"]);
    d.n_seasons = Rcpp::as<int>(spec["n_seasons"]);
    d.n_species = Rcpp::as<int>(spec["n_species"]);
    d.p_psi1    = d.X_psi1.ncol();
    d.p_p       = d.X_p.ncol();
    d.p_gam     = d.X_gamma.ncol();
    d.p_eps     = d.X_eps.ncol();

    // Per-species (n_valid, n_det) [n_sites x n_seasons] integer matrices,
    // supplied as R lists of matrices (length n_species).
    Rcpp::List nv = Rcpp::as<Rcpp::List>(spec["n_valid"]);
    Rcpp::List nd = Rcpp::as<Rcpp::List>(spec["n_det"]);
    if ((int) nv.size() != d.n_species || (int) nd.size() != d.n_species)
        Rcpp::stop("n_valid / n_det lists must have length n_species");
    const std::size_t sz = (std::size_t) d.n_sites * d.n_seasons;
    d.nvalid.assign(d.n_species, std::vector<int>());
    d.ndet.assign(d.n_species, std::vector<int>());
    for (int s = 0; s < d.n_species; ++s) {
        IntegerMatrix vs = Rcpp::as<IntegerMatrix>(nv[s]);
        IntegerMatrix ds = Rcpp::as<IntegerMatrix>(nd[s]);
        if (vs.nrow() != d.n_sites || vs.ncol() != d.n_seasons ||
            ds.nrow() != d.n_sites || ds.ncol() != d.n_seasons)
            Rcpp::stop("n_valid / n_det matrices must be n_sites x n_seasons");
        d.nvalid[s].assign(sz, 0);
        d.ndet[s].assign(sz, 0);
        for (int t = 0; t < d.n_seasons; ++t)
            for (int i = 0; i < d.n_sites; ++i) {
                d.nvalid[s][(std::size_t) i + (std::size_t) t * d.n_sites] = vs(i, t);
                d.ndet[s][(std::size_t) i + (std::size_t) t * d.n_sites]   = ds(i, t);
            }
    }
    ms_dyn_occu_nuts_layout(d);
    return d;
}

// Per-species forward-backward smoothing (mirrors .ms_dyn_occu_fb_vec, per-site).
// Inputs: per-site psi1 / p (unclamped plogis values) and the SHARED per-site
// gamma / eps; the species' per-(site,season) (nvalid, ndet). Accumulates:
//   ll                          # marginal log-lik summed over finite sites
//   score_psi1[i] = w1[i] - psi1[i]                 (season-1 smoothed posterior)
//   score_p[i]    = sum_t w_it (ndet - nvalid p)    (p CLAMPED, as in R)
//   g_gam_add[i]  = col_y[i] - gamma[i] col_n[i]    (gamma UNCLAMPED)
//   g_eps_add[i]  = ext_y[i] - eps[i]   ext_n[i]    (eps   UNCLAMPED)
// The score subtractions use the UNCLAMPED outer values (psi1, gamma, eps) while
// the HMM recursion and p_score use internally clamped copies -- exactly the R
// oracle's split. Returns the per-species ll; writes the four score vectors.
inline double ms_dyn_occu_fb(const MsDynOccuNutsData& d, int s,
                             const double* psi1, const double* p,
                             const double* gamma, const double* eps,
                             double* score_psi1, double* score_p,
                             double* g_gam_add, double* g_eps_add) {
    const int Ns = d.n_sites, T = d.n_seasons;
    const std::vector<int>& nv = d.nvalid[s];
    const std::vector<int>& nd = d.ndet[s];
    std::vector<double> cs(T), A0(T), A1(T), emocc(T), emunocc(T), wsm(T);
    double ll = 0.0;
    for (int i = 0; i < Ns; ++i) {
        const double ps1c = mdo_clamp(psi1[i]);
        const double pc   = mdo_clamp(p[i]);
        const double gc   = mdo_clamp(gamma[i]);
        const double ec   = mdo_clamp(eps[i]);
        const double lpc  = std::log(pc), l1pc = std::log(1.0 - pc);
        // emissions per season
        for (int t = 0; t < T; ++t) {
            const std::size_t o = (std::size_t) i + (std::size_t) t * Ns;
            const int nvt = nv[o], ndt = nd[o];
            emocc[t]   = std::exp(ndt * lpc + (nvt - ndt) * l1pc);
            emunocc[t] = (ndt > 0) ? 0.0 : 1.0;
        }
        // forward (scaled)
        bool finite = true;
        double v0 = (1.0 - ps1c) * emunocc[0], v1 = ps1c * emocc[0];
        double ct = v0 + v1;
        cs[0] = ct;
        if (ct > 0.0) { A0[0] = v0 / ct; A1[0] = v1 / ct; }
        else { A0[0] = 0.0; A1[0] = 0.0; finite = false; }
        double site_ll = (ct > 0.0) ? std::log(ct) : -INFINITY;
        for (int t = 1; t < T; ++t) {
            const double pr1 = A1[t - 1] * (1.0 - ec) + A0[t - 1] * gc;
            const double pr0 = A1[t - 1] * ec         + A0[t - 1] * (1.0 - gc);
            v0 = pr0 * emunocc[t]; v1 = pr1 * emocc[t];
            ct = v0 + v1; cs[t] = ct;
            if (ct > 0.0) { A0[t] = v0 / ct; A1[t] = v1 / ct; }
            else { A0[t] = 0.0; A1[t] = 0.0; finite = false; }
            site_ll += (ct > 0.0) ? std::log(ct) : -INFINITY;
        }
        if (finite && std::isfinite(site_ll)) ll += site_ll;
        // backward (scaled) + smoothed marginals / pairwise joints
        wsm[T - 1] = A1[T - 1];
        double bw0 = 1.0, bw1 = 1.0;
        double col_y = 0.0, col_n = 0.0, ext_y = 0.0, ext_n = 0.0;
        for (int t = T - 2; t >= 0; --t) {
            const double bb0 = emunocc[t + 1] * bw0;
            const double bb1 = emocc[t + 1]   * bw1;
            const double invc = (cs[t + 1] > 0.0) ? 1.0 / cs[t + 1] : 0.0;
            const double xi01 = A0[t] * gc         * bb1 * invc;   // colonisation
            const double xi00 = A0[t] * (1.0 - gc) * bb0 * invc;
            const double xi10 = A1[t] * ec         * bb0 * invc;   // extinction
            const double xi11 = A1[t] * (1.0 - ec) * bb1 * invc;
            col_y += xi01; col_n += (xi00 + xi01);
            ext_y += xi10; ext_n += (xi10 + xi11);
            const double nb0 = ((1.0 - gc) * bb0 + gc         * bb1) * invc;
            const double nb1 = (ec         * bb0 + (1.0 - ec) * bb1) * invc;
            bw0 = nb0; bw1 = nb1;
            wsm[t] = A1[t] * bw1;
        }
        // detection score: sum_t w_it (ndet - nvalid * p_clamped)
        double ps = 0.0;
        for (int t = 0; t < T; ++t) {
            const std::size_t o = (std::size_t) i + (std::size_t) t * Ns;
            ps += wsm[t] * (nd[o] - nv[o] * pc);
        }
        score_psi1[i] = wsm[0] - psi1[i];        // unclamped psi1
        score_p[i]    = ps;
        g_gam_add[i]  = col_y - gamma[i] * col_n;  // unclamped gamma
        g_eps_add[i]  = ext_y - eps[i]   * ext_n;  // unclamped eps
    }
    return ll;
}

// Joint log-posterior + gradient over the packed vector
//   theta = (mu, {z_s} species-major, chol_psi1, chol_p, global(gam,eps)).
// NUTS maximises, so this returns the log-posterior and writes the gradient into
// `g` (length d.total). Mirrors .tobs_ms_dyn_occu_nuts_logpost.
inline double ms_dyn_occu_nuts_eval(const MsDynOccuNutsData& d, const double* th,
                                    double sigma_beta, const MsDynOccuPri& pr,
                                    double* g) {
    const int P = d.P_tot, S = d.n_species;
    const int p_psi1 = d.p_psi1, p_p = d.p_p, p_gam = d.p_gam, p_eps = d.p_eps;
    const int n_sites = d.n_sites;
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;

    const double* mu     = th + d.mu_off;
    const double* z      = th + d.b_off;
    const double* global = th + d.global_off;
    double* g_mu     = g + d.mu_off;
    double* g_z      = g + d.b_off;
    double* g_global = g + d.global_off;

    // Shared per-site transition probabilities (gamma, eps).
    std::vector<double> gamma(n_sites), eps(n_sites);
    for (int i = 0; i < n_sites; ++i) {
        double eg = 0.0, ee = 0.0;
        for (int k = 0; k < p_gam; ++k) eg += d.X_gamma(i, k) * global[k];
        for (int k = 0; k < p_eps; ++k) ee += d.X_eps(i, k)   * global[p_gam + k];
        gamma[i] = mdo_plogis(eg);
        eps[i]   = mdo_plogis(ee);
    }

    std::vector<double> C_psi1, C_p;
    chol_unpack_cpp(th + d.chol_psi1_off, p_psi1, C_psi1);
    chol_unpack_cpp(th + d.chol_p_off,    p_p,    C_p);

    std::vector<double> A_psi1((std::size_t) p_psi1 * p_psi1, 0.0),
                        A_p((std::size_t) p_p * p_p, 0.0);

    // Per-species partials (reduced serially in species order -> deterministic /
    // byte-exact against the R oracle, and safe under OpenMP).
    std::vector<double> gmu_s((std::size_t) S * P, 0.0), lp_s(S, 0.0);
    std::vector<double> Apsi_s((std::size_t) S * p_psi1 * p_psi1, 0.0);
    std::vector<double> Ap_s((std::size_t) S * p_p * p_p, 0.0);
    std::vector<double> ggam_s((std::size_t) S * n_sites, 0.0);
    std::vector<double> geps_s((std::size_t) S * n_sites, 0.0);

    #pragma omp parallel for schedule(static)
    for (int s = 0; s < S; ++s) {
        const double* z_s  = z + s * P;
        const double* zpsi1 = z_s;
        const double* zp    = z_s + p_psi1;
        std::vector<double> b_psi1(p_psi1), b_p(p_p);
        for (int i = 0; i < p_psi1; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j)
                v += C_psi1[(std::size_t) i * p_psi1 + j] * zpsi1[j];
            b_psi1[i] = v;
        }
        for (int i = 0; i < p_p; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j) v += C_p[(std::size_t) i * p_p + j] * zp[j];
            b_p[i] = v;
        }
        std::vector<double> psi1(n_sites), pp(n_sites);
        for (int i = 0; i < n_sites; ++i) {
            double e_psi1 = 0.0, e_p = 0.0;
            for (int k = 0; k < p_psi1; ++k)
                e_psi1 += d.X_psi1(i, k) * (mu[k] + b_psi1[k]);
            for (int k = 0; k < p_p; ++k)
                e_p += d.X_p(i, k) * (mu[p_psi1 + k] + b_p[k]);
            psi1[i] = mdo_plogis(e_psi1);
            pp[i]   = mdo_plogis(e_p);
        }
        std::vector<double> s_psi1(n_sites), s_p(n_sites);
        double* ggam = &ggam_s[(std::size_t) s * n_sites];
        double* geps = &geps_s[(std::size_t) s * n_sites];
        lp_s[s] = ms_dyn_occu_fb(d, s, psi1.data(), pp.data(),
                                 gamma.data(), eps.data(),
                                 s_psi1.data(), s_p.data(), ggam, geps);

        // design-sandwiched eta-gradient on the psi1 / p arms (grad_b).
        std::vector<double> gbpsi1(p_psi1, 0.0), gbp(p_p, 0.0);
        double* gmu_loc = &gmu_s[(std::size_t) s * P];
        for (int i = 0; i < n_sites; ++i) {
            const double gpsi1 = s_psi1[i];
            for (int k = 0; k < p_psi1; ++k) {
                const double gx = gpsi1 * d.X_psi1(i, k);
                gmu_loc[k] += gx;
                gbpsi1[k]  += gx;
            }
            const double gp = s_p[i];
            for (int k = 0; k < p_p; ++k) {
                const double gx = gp * d.X_p(i, k);
                gmu_loc[p_psi1 + k] += gx;
                gbp[k]              += gx;
            }
        }
        // z gradient (data part) = C' grad_b (disjoint per-species write).
        double* gz_s = g_z + s * P;
        for (int v = 0; v < p_psi1; ++v) {
            double sg = 0.0;
            for (int i = v; i < p_psi1; ++i)
                sg += C_psi1[(std::size_t) i * p_psi1 + v] * gbpsi1[i];
            gz_s[v] += sg;
        }
        for (int v = 0; v < p_p; ++v) {
            double sg = 0.0;
            for (int i = v; i < p_p; ++i)
                sg += C_p[(std::size_t) i * p_p + v] * gbp[i];
            gz_s[p_psi1 + v] += sg;
        }
        // per-species chol accumulators A_arm = grad_b z'.
        double* Apsi = &Apsi_s[(std::size_t) s * p_psi1 * p_psi1];
        for (int i = 0; i < p_psi1; ++i)
            for (int j = 0; j <= i; ++j)
                Apsi[(std::size_t) i * p_psi1 + j] = gbpsi1[i] * zpsi1[j];
        double* Ap = &Ap_s[(std::size_t) s * p_p * p_p];
        for (int i = 0; i < p_p; ++i)
            for (int j = 0; j <= i; ++j)
                Ap[(std::size_t) i * p_p + j] = gbp[i] * zp[j];
    }

    // serial reduction in species order -> byte-identical to the serial path.
    double lp = 0.0;
    std::vector<double> g_gam(n_sites, 0.0), g_eps(n_sites, 0.0);
    for (int s = 0; s < S; ++s) {
        const double* gmu_loc = &gmu_s[(std::size_t) s * P];
        for (int k = 0; k < P; ++k) g_mu[k] += gmu_loc[k];
        lp += lp_s[s];
        const double* Apsi = &Apsi_s[(std::size_t) s * p_psi1 * p_psi1];
        for (int t = 0; t < p_psi1 * p_psi1; ++t) A_psi1[t] += Apsi[t];
        const double* Ap = &Ap_s[(std::size_t) s * p_p * p_p];
        for (int t = 0; t < p_p * p_p; ++t) A_p[t] += Ap[t];
        const double* ggam = &ggam_s[(std::size_t) s * n_sites];
        const double* geps = &geps_s[(std::size_t) s * n_sites];
        for (int i = 0; i < n_sites; ++i) { g_gam[i] += ggam[i]; g_eps[i] += geps[i]; }
    }

    // ---- z prior: standard normal over the whole per-species block ----
    for (int j = 0; j < S * P; ++j) {
        const double zz = z[j];
        g_z[j] -= zz;
        lp     += -0.5 * zz * zz;
    }

    // ---- chol coords: data gradient (via b = C z) + coordinate hyperprior ----
    lp += chol_data_grad_noncentered(A_psi1, C_psi1, p_psi1, th + d.chol_psi1_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_psi1_off);
    lp += chol_data_grad_noncentered(A_p, C_p, p_p, th + d.chol_p_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_p_off);

    // ---- community-mean priors: N(0, sigma.beta^2) over all mu ----
    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < P; ++k) {
        g_mu[k] -= inv_sb2 * mu[k];
        lp      += -0.5 * inv_sb2 * mu[k] * mu[k];
    }

    // ---- shared-transition (global) priors + design-sandwiched gradient ----
    for (int k = 0; k < d.G; ++k) lp += -0.5 * inv_sb2 * global[k] * global[k];
    for (int k = 0; k < p_gam; ++k) {
        double sg = 0.0;
        for (int i = 0; i < n_sites; ++i) sg += d.X_gamma(i, k) * g_gam[i];
        g_global[k] = sg - inv_sb2 * global[k];
    }
    for (int k = 0; k < p_eps; ++k) {
        double sg = 0.0;
        for (int i = 0; i < n_sites; ++i) sg += d.X_eps(i, k) * g_eps[i];
        g_global[p_gam + k] = sg - inv_sb2 * global[p_gam + k];
    }
    return lp;
}

struct MsDynOccuNutsModel {
    MsDynOccuNutsData d;
    double sigma_beta = 5.0;
    MsDynOccuPri pr;
};

inline void ms_dyn_occu_nuts_full_grad(const std::vector<double>& params,
                                       const tulpa::ModelData& data,
                                       const tulpa::ParamLayout& /*layout*/,
                                       std::vector<double>& grad,
                                       double* log_post_out) {
    const MsDynOccuNutsModel* m =
        static_cast<const MsDynOccuNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = ms_dyn_occu_nuts_eval(m->d, params.data(), m->sigma_beta,
                                            m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

inline MsDynOccuPri ms_dyn_occu_pri_from_list(const Rcpp::List& pri) {
    MsDynOccuPri pr;
    pr.logdiag_mean = Rcpp::as<double>(pri["chol_logdiag_mean"]);
    pr.logdiag_sd   = Rcpp::as<double>(pri["chol_logdiag_sd"]);
    pr.offdiag_sd   = Rcpp::as<double>(pri["chol_offdiag_sd"]);
    return pr;
}

} // namespace tulpaObs

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .tobs_ms_dyn_occu_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_ms_dyn_occu_nuts_joint_logpost(Rcpp::List spec,
                                              Rcpp::NumericVector theta,
                                              Rcpp::List pri, double sigma_beta) {
    tulpaObs::MsDynOccuNutsData d = tulpaObs::ms_dyn_occu_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    tulpaObs::MsDynOccuPri pr = tulpaObs::ms_dyn_occu_pri_from_list(pri);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::ms_dyn_occu_nuts_eval(
        d, theta.begin(), sigma_beta, pr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the community dynamic occupancy target via tulpa's engine and the
// FullGradFn (gradient mode "H"). `theta0` is the warm start (the Laplace-EM
// mode); `inv_metric` an optional length-n_params inverse-mass diagonal.
// [[Rcpp::export]]
Rcpp::List cpp_ms_dyn_occu_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                                Rcpp::List pri, double sigma_beta,
                                Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                                int n_iter, int n_warmup, int max_treedepth,
                                double adapt_delta, int seed, bool verbose) {
    tulpaObs::MsDynOccuNutsModel m;
    m.d = tulpaObs::ms_dyn_occu_nuts_build_data(spec);
    m.sigma_beta = sigma_beta;
    m.pr = tulpaObs::ms_dyn_occu_pri_from_list(pri);
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    return tulpaObs::run_tulpa_nuts(
        &tulpaObs::ms_dyn_occu_nuts_full_grad, &m, m.d.total,
        theta0, sigma_beta, inv_metric,
        n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose,
        "ms_dyn_occu", m.d.n_sites);
}
