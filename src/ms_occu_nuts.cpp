// ms_occu_nuts.cpp
// C++ joint log-posterior + gradient for the community / multispecies
// single-season occupancy (ms_occu()) NUTS target. The R reference
// .tobs_ms_occu_nuts_logpost (R/ms_occu_nuts.R) is the oracle; this port mirrors
// it and is cross-checked against it before driving tulpa's NUTS engine.
//
// The community single-season occupancy is per-species two-state occupancy with
// Gaussian community hyperpriors on the per-species coefficients:
//   z_{s,i}        ~ Bernoulli(psi_{s,i})
//   y_{s,i,j}|z=1  ~ Bernoulli(p_{s,i})
//   logit psi_{s,i} = X_psi_i . (mu_psi + b_psi_s)
//   logit p_{s,i}   = X_p_i   . (mu_p   + b_p_s)
//   b_psi_s ~ N(0, Sigma_psi),  b_p_s ~ N(0, Sigma_p)
// The latent z integrates out per species-site in closed form (the two-state
// mixture, ms_occu_kernel.h); the Laplace-EM (.tobs_community_em) profiles the
// per-species deviations and per-arm community covariances out, returning a
// Gaussian community-mean posterior. NUTS instead samples EVERYTHING jointly --
// the community means, the per-species deviations, AND the two INDEPENDENT
// per-arm community covariances Sigma_psi / Sigma_p -- from the exact joint
// posterior, which gives calibrated (non-Gaussian) community-mean / covariance
// intervals and removes the Laplace small-cluster attenuation of the variance
// components.
//
// NON-CENTERED parameterisation (mirrors ms_abun_nuts.cpp): the per-species
// block carries the whitened standard-normal z_s, the deviation is reconstructed
// per arm as b_{s,arm} = C_arm z_{s,arm} (C_arm the log-Cholesky factor of
// Sigma_arm), so the community covariance leaves the b-prior entirely and enters
// ONLY the data term through b = C z. This breaks the centered b<->Sigma funnel.
// The target factorises as
//   log p = sum_{s,i} log L_{s,i}(theta)                # per-species-site marginal
//         - 0.5 ||mu_coef||^2 / sigma.beta^2            # community-mean priors
//         - 0.5 sum_s ||z_s||^2                         # whitened RE prior (N(0,I))
//         + log p(Sigma coords)                         # log-Cholesky hyperpriors
// where L_{s,i} is the occupancy two-state marginal exposed by
// compute_ms_occu_site() (ms_occu_kernel.h, matching the R oracle); it returns
// grad_eta_psi and grad_eta_p, so the coefficient gradient is the
// design-sandwiched eta-gradient -- no new likelihood math. The chol gradient
// flows from the data term via b = C z (chol_data_grad_noncentered) plus the
// coordinate hyperprior; the arm Cholesky factors use the shared helpers in
// community_chol.h.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
#include "ms_occu_kernel.h"
#include "community_chol.h"

#include "nuts_engine.h"
using namespace Rcpp;

namespace tulpaObs {

// Log-Cholesky coordinate hyperprior scalars for the community covariances.
struct MsOccuPri {
    double logdiag_mean = 0.0, logdiag_sd = 1.5, offdiag_sd = 1.0;
};

// Marshalled per-fit data for the community occupancy NUTS target. X_psi and
// X_p are site-level designs (n_sites rows, shared across species); the
// per-species site summaries (n_valid, n_det) drive the occupancy marginal.
struct MsOccuNutsData {
    int n_sites = 0, n_species = 0, p_psi = 0, p_p = 0;
    NumericMatrix X_psi;                          // n_sites x p_psi
    NumericMatrix X_p;                            // n_sites x p_p
    std::vector<MsOccuSiteSummary> summ;          // [s]

    // Packed-coordinate layout (mirrors .tobs_ms_occu_nuts_layout).
    int P_tot = 0, q_psi = 0, q_p = 0, total = 0;
    int mu_off = 0, b_off = 0, chol_psi_off = 0, chol_p_off = 0;
};

inline void ms_occu_nuts_layout(MsOccuNutsData& d) {
    d.P_tot      = d.p_psi + d.p_p;
    d.q_psi      = d.p_psi * (d.p_psi + 1) / 2;
    d.q_p        = d.p_p   * (d.p_p   + 1) / 2;
    d.mu_off     = 0;
    d.b_off      = d.P_tot;
    d.chol_psi_off = d.P_tot + d.n_species * d.P_tot;
    d.chol_p_off   = d.chol_psi_off + d.q_psi;
    d.total      = d.chol_p_off + d.q_p;
}

inline MsOccuNutsData ms_occu_nuts_build_data(const Rcpp::List& spec) {
    MsOccuNutsData d;
    d.X_psi     = Rcpp::as<NumericMatrix>(spec["X_psi"]);
    d.X_p       = Rcpp::as<NumericMatrix>(spec["X_p"]);
    d.n_sites   = Rcpp::as<int>(spec["n_sites"]);
    d.n_species = Rcpp::as<int>(spec["n_species"]);
    d.p_psi     = d.X_psi.ncol();
    d.p_p       = d.X_p.ncol();

    // Per-species (n_valid, n_det) site summaries, packed as integer matrices
    // [n_sites x n_species] (column-major species, matching the R packer).
    IntegerMatrix nv = Rcpp::as<IntegerMatrix>(spec["n_valid"]);
    IntegerMatrix nd = Rcpp::as<IntegerMatrix>(spec["n_det"]);
    if (nv.nrow() != d.n_sites || nv.ncol() != d.n_species)
        Rcpp::stop("n_valid must be n_sites x n_species");
    if (nd.nrow() != d.n_sites || nd.ncol() != d.n_species)
        Rcpp::stop("n_det must be n_sites x n_species");
    d.summ.assign(d.n_species, MsOccuSiteSummary());
    for (int s = 0; s < d.n_species; ++s) {
        MsOccuSiteSummary& su = d.summ[s];
        su.n_valid.assign((std::size_t) d.n_sites, 0);
        su.n_det.assign((std::size_t) d.n_sites, 0);
        su.any_det.assign((std::size_t) d.n_sites, 0);
        for (int i = 0; i < d.n_sites; ++i) {
            su.n_valid[i] = nv(i, s);
            su.n_det[i]   = nd(i, s);
            su.any_det[i] = (nd(i, s) > 0) ? 1 : 0;
        }
    }
    ms_occu_nuts_layout(d);
    return d;
}

// Joint log-posterior + gradient over the packed vector
//   theta = (mu, {b_s} species-major, chol_psi, chol_p).
// NUTS maximises, so this returns the log-posterior (no negation) and writes the
// gradient into `g` (length d.total). Mirrors .tobs_ms_occu_nuts_logpost.
inline double ms_occu_nuts_eval(const MsOccuNutsData& d, const double* th,
                                double sigma_beta, const MsOccuPri& pr,
                                double* g) {
    const int P = d.P_tot, S = d.n_species, p_psi = d.p_psi, p_p = d.p_p;
    const int n_sites = d.n_sites;
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;

    const double* mu = th + d.mu_off;
    const double* z  = th + d.b_off;       // species-major z (whitened)
    double* g_mu = g + d.mu_off;
    double* g_z  = g + d.b_off;

    // Cholesky factors per arm (Sigma_arm = C_arm C_arm', row-major). Under the
    // non-centered map the per-species deviation is b_{s,arm} = C_arm z_{s,arm},
    // so the covariance enters only the data term.
    std::vector<double> C_psi, C_p;
    chol_unpack_cpp(th + d.chol_psi_off, p_psi, C_psi);
    chol_unpack_cpp(th + d.chol_p_off,   p_p,   C_p);

    // chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
    std::vector<double> A_psi((std::size_t) p_psi * p_psi, 0.0),
                        A_p((std::size_t) p_p * p_p, 0.0);

    // Per-species partial sums (reduced serially in species order afterwards, so
    // the reduction is deterministic / thread-count independent and byte-exact
    // against the R oracle).
    std::vector<double> gmu_s((std::size_t) S * P, 0.0), lp_s(S, 0.0);
    std::vector<double> Apsi_s((std::size_t) S * p_psi * p_psi, 0.0);
    std::vector<double> Ap_s((std::size_t) S * p_p * p_p, 0.0);

    // Scratch is sized once per thread and reused across species: the widths are
    // known before the loop, and the per-species work below is O(n_sites), so a
    // fresh allocation per iteration buys nothing. b_* and eta_* are fully
    // overwritten each pass; gb* accumulate and are re-zeroed explicitly.
    #pragma omp parallel
    {
    std::vector<double> b_psi(p_psi), b_p(p_p), gbpsi(p_psi), gbp(p_p);
    std::vector<double> eta_psi(n_sites), eta_p(n_sites);
    MsOccuSiteResult res;

    #pragma omp for schedule(static)
    for (int s = 0; s < S; ++s) {
        const double* z_s = z + s * P;
        const double* zpsi = z_s;             // length p_psi
        const double* zp   = z_s + p_psi;     // length p_p
        std::fill(gbpsi.begin(), gbpsi.end(), 0.0);
        std::fill(gbp.begin(), gbp.end(), 0.0);
        // reconstruct b = C z (lower-triangular C)
        for (int i = 0; i < p_psi; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j)
                v += C_psi[(std::size_t) i * p_psi + j] * zpsi[j];
            b_psi[i] = v;
        }
        for (int i = 0; i < p_p; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j)
                v += C_p[(std::size_t) i * p_p + j] * zp[j];
            b_p[i] = v;
        }
        // per-site eta over the shared site-level designs
        for (int i = 0; i < n_sites; ++i) {
            double e_psi = 0.0, e_p = 0.0;
            for (int k = 0; k < p_psi; ++k)
                e_psi += d.X_psi(i, k) * (mu[k] + b_psi[k]);
            for (int k = 0; k < p_p; ++k)
                e_p += d.X_p(i, k) * (mu[p_psi + k] + b_p[k]);
            eta_psi[i] = e_psi;
            eta_p[i]   = e_p;
        }
        compute_ms_occu_site(eta_psi.data(), eta_p.data(), d.summ[s], n_sites,
                             true, res);
        lp_s[s] = res.log_lik;

        // design-sandwiched eta-gradient: shared by mu and b on each arm.
        double* gmu_loc = &gmu_s[(std::size_t) s * P];
        for (int i = 0; i < n_sites; ++i) {
            const double gpsi = res.grad_eta_psi[i];
            for (int k = 0; k < p_psi; ++k) {
                const double gx = gpsi * d.X_psi(i, k);
                gmu_loc[k] += gx;
                gbpsi[k]   += gx;
            }
            const double gp = res.grad_eta_p[i];
            for (int k = 0; k < p_p; ++k) {
                const double gx = gp * d.X_p(i, k);
                gmu_loc[p_psi + k] += gx;
                gbp[k]             += gx;
            }
        }
        // z gradient (data part) = C' grad_b -- disjoint per-species write.
        double* gz_s = g_z + s * P;
        for (int v = 0; v < p_psi; ++v) {
            double sg = 0.0;
            for (int i = v; i < p_psi; ++i)
                sg += C_psi[(std::size_t) i * p_psi + v] * gbpsi[i];
            gz_s[v] += sg;
        }
        for (int v = 0; v < p_p; ++v) {
            double sg = 0.0;
            for (int i = v; i < p_p; ++i)
                sg += C_p[(std::size_t) i * p_p + v] * gbp[i];
            gz_s[p_psi + v] += sg;
        }
        // per-species chol accumulators A_arm = grad_b z' (reduced below).
        double* Apsi = &Apsi_s[(std::size_t) s * p_psi * p_psi];
        for (int i = 0; i < p_psi; ++i)
            for (int j = 0; j <= i; ++j)
                Apsi[(std::size_t) i * p_psi + j] = gbpsi[i] * zpsi[j];
        double* Ap = &Ap_s[(std::size_t) s * p_p * p_p];
        for (int i = 0; i < p_p; ++i)
            for (int j = 0; j <= i; ++j)
                Ap[(std::size_t) i * p_p + j] = gbp[i] * zp[j];
    }
    }  // omp parallel

    // serial reduction in species order -> byte-identical to the serial path.
    double lp = 0.0;
    for (int s = 0; s < S; ++s) {
        const double* gmu_loc = &gmu_s[(std::size_t) s * P];
        for (int k = 0; k < P; ++k) g_mu[k] += gmu_loc[k];
        lp += lp_s[s];
        const double* Apsi = &Apsi_s[(std::size_t) s * p_psi * p_psi];
        for (int t = 0; t < p_psi * p_psi; ++t) A_psi[t] += Apsi[t];
        const double* Ap = &Ap_s[(std::size_t) s * p_p * p_p];
        for (int t = 0; t < p_p * p_p; ++t) A_p[t] += Ap[t];
    }

    // ---- z prior: standard normal over the whole per-species block ----
    for (int j = 0; j < S * P; ++j) {
        const double zz = z[j];
        g_z[j] -= zz;
        lp     += -0.5 * zz * zz;
    }

    // ---- chol coords: data gradient (via b = C z) + coordinate hyperprior ----
    lp += chol_data_grad_noncentered(A_psi, C_psi, p_psi, th + d.chol_psi_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_psi_off);
    lp += chol_data_grad_noncentered(A_p, C_p, p_p, th + d.chol_p_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_p_off);

    // ---- community-mean priors: N(0, sigma.beta^2) on the coefficient means ----
    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < p_psi + p_p; ++k) {
        g_mu[k] -= inv_sb2 * mu[k];
        lp      += -0.5 * inv_sb2 * mu[k] * mu[k];
    }
    return lp;
}

// NUTS model carrying the marshalled data + prior scales; the FullGradFn reaches
// it through ModelData.model_response_data.
struct MsOccuNutsModel {
    MsOccuNutsData d;
    double sigma_beta = 5.0;
    MsOccuPri pr;
};

// FullGradFn: log-posterior + gradient over the entire parameter vector.
inline void ms_occu_nuts_full_grad(const std::vector<double>& params,
                                   const tulpa::ModelData& data,
                                   const tulpa::ParamLayout& /*layout*/,
                                   std::vector<double>& grad, double* log_post_out) {
    const MsOccuNutsModel* m =
        static_cast<const MsOccuNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = ms_occu_nuts_eval(m->d, params.data(), m->sigma_beta,
                                        m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

inline MsOccuPri ms_occu_pri_from_list(const Rcpp::List& pri) {
    MsOccuPri pr;
    pr.logdiag_mean = Rcpp::as<double>(pri["chol_logdiag_mean"]);
    pr.logdiag_sd   = Rcpp::as<double>(pri["chol_logdiag_sd"]);
    pr.offdiag_sd   = Rcpp::as<double>(pri["chol_offdiag_sd"]);
    return pr;
}

} // namespace tulpaObs

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .tobs_ms_occu_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                          Rcpp::List pri, double sigma_beta) {
    tulpaObs::MsOccuNutsData d = tulpaObs::ms_occu_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    tulpaObs::MsOccuPri pr = tulpaObs::ms_occu_pri_from_list(pri);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::ms_occu_nuts_eval(
        d, theta.begin(), sigma_beta, pr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the community occupancy target via tulpa's engine and the
// FullGradFn (gradient mode "H"). `theta0` is the warm-start (the Laplace-EM
// mode); `inv_metric` an optional length-n_params inverse-mass diagonal (the
// Laplace curvature). Returns draws + diagnostics.
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                            Rcpp::List pri, double sigma_beta,
                            Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                            int n_iter, int n_warmup, int max_treedepth,
                            double adapt_delta, int seed, bool verbose) {
    tulpaObs::MsOccuNutsModel m;
    m.d = tulpaObs::ms_occu_nuts_build_data(spec);
    m.sigma_beta = sigma_beta;
    m.pr = tulpaObs::ms_occu_pri_from_list(pri);
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    return tulpaObs::run_tulpa_nuts(
        &tulpaObs::ms_occu_nuts_full_grad, &m, m.d.total,
        theta0, sigma_beta, inv_metric,
        n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose,
        "ms_occu", m.d.n_sites);
}
