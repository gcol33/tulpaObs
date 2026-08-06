// Community relative-abundance (count) NUTS: the non-spatial community GLMM
// sampled over the exact joint posterior (community means, per-species
// coefficient deviations, and the community covariance), via tulpa's in-tree
// FullGradFn engine. This is the reduced counterpart of ms_abun_nuts.cpp: the
// community count has NO detection arm and NO latent-N marginalisation, so the
// per-(species, site) contribution is a plain GLMM log-likelihood over the
// observed response.
//
// Three response families share one target:
//   Poisson   log p = sum_{s,i} [ y eta - exp(eta) - lgamma(y+1) ]
//   NegBin    log p = sum_{s,i} dnbinom(y; r_s, mu = exp(eta)); r_s a per-species
//             dispersion RE log_r_s ~ N(mu_log_r, sigma_log_r^2) (a second
//             community arm, matching ms_abun / the ms_count Laplace-EM)
//   Gaussian  log p = sum_{s,i} dnorm(y; eta, sqrt(phi_s)); phi_s a per-species
//             FREE residual variance (matching the ms_count Laplace outer loop,
//             which estimates each phi_s with no community prior)
// with eta_{s,i} = X_i . (mu_beta + b_beta_s) (identity for Gaussian, log else).
// NON-CENTERED: b_{s,arm} = C_arm z_{s,arm}, z ~ N(0, I), so the community
// covariance (log-Cholesky C) enters only the data term. The arm Cholesky helpers
// are the shared community_chol.h. An NA entry in the response matrix is a missing
// site x species observation: it is masked out of the (species, site) data sum
// (matching the Laplace-EM per-species valid subsets), so a species contributes
// only its observed sites. Byte-exact vs the R oracle .tobs_ms_count_nuts_logpost
// (gcol33/tulpaObs#117).

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
#include "community_chol.h"
#include "tobs_math.h"

#include "nuts_engine.h"
using namespace Rcpp;

namespace tulpaObs {

// Response family codes.
// MSC_BERN is the jsdm() response: presence/absence observed directly, logit
// link. It carries no dispersion parameter, so it shares the Poisson parameter
// layout (community means + per-species deviations only).
enum MsCountFamily { MSC_POIS = 0, MSC_NB = 1, MSC_GAUSS = 2, MSC_BERN = 3 };

struct MsCountNutsData {
    int n_sites = 0, n_species = 0, p_beta = 0;
    int family = MSC_POIS;
    std::vector<double> y;          // n_sites x n_species, column-major (species-major)
    std::vector<double> lgy;        // lgamma(y + 1) (Poisson only), same layout
    std::vector<unsigned char> valid; // 1 = observed, 0 = missing (NA), same layout
    NumericMatrix X;                // n_sites x p_beta
    // Gaussian: weakly-informative N(logphi_mean, logphi_sd^2) on each free log_phi_s.
    double logphi_mean = 0.0, logphi_sd = 2.0;
    // layout: mu (P), z species-major (S*P), chol_beta (q_beta),
    //         [chol_logr (1) if NB], [log_phi (S) if Gaussian].
    int P = 0, q_beta = 0, total = 0;
    int mu_off = 0, b_off = 0, chol_beta_off = 0, chol_logr_off = 0, logphi_off = 0;
};

inline void ms_count_nuts_layout(MsCountNutsData& d) {
    const bool nb = (d.family == MSC_NB), gauss = (d.family == MSC_GAUSS);
    d.P             = d.p_beta + (nb ? 1 : 0);
    d.q_beta        = d.p_beta * (d.p_beta + 1) / 2;
    d.mu_off        = 0;
    d.b_off         = d.P;
    d.chol_beta_off = d.P + d.n_species * d.P;
    d.chol_logr_off = d.chol_beta_off + d.q_beta;   // valid only under NB
    d.logphi_off    = d.chol_beta_off + d.q_beta;   // valid only under Gaussian
    d.total         = d.chol_beta_off + d.q_beta + (nb ? 1 : 0)
                    + (gauss ? d.n_species : 0);
}

inline MsCountNutsData ms_count_nuts_build_data(const Rcpp::List& spec) {
    MsCountNutsData d;
    d.X = Rcpp::as<NumericMatrix>(spec["X"]);
    NumericMatrix ym = Rcpp::as<NumericMatrix>(spec["y"]);   // n_sites x n_species
    d.n_sites = d.X.nrow(); d.p_beta = d.X.ncol(); d.n_species = ym.ncol();
    if (ym.nrow() != d.n_sites) Rcpp::stop("y must have n_sites rows.");
    std::string fam = spec.containsElementNamed("family")
                    ? Rcpp::as<std::string>(spec["family"]) : "poisson";
    d.family = (fam == "negbin")    ? MSC_NB
             : (fam == "gaussian")  ? MSC_GAUSS
             : (fam == "bernoulli") ? MSC_BERN
                                    : MSC_POIS;
    if (d.family == MSC_GAUSS) {
        if (spec.containsElementNamed("logphi_mean"))
            d.logphi_mean = Rcpp::as<double>(spec["logphi_mean"]);
        if (spec.containsElementNamed("logphi_sd"))
            d.logphi_sd = Rcpp::as<double>(spec["logphi_sd"]);
    }
    d.y.assign((std::size_t) d.n_sites * d.n_species, 0.0);
    d.valid.assign((std::size_t) d.n_sites * d.n_species, 1);
    if (d.family == MSC_POIS) d.lgy.assign((std::size_t) d.n_sites * d.n_species, 0.0);
    for (int s = 0; s < d.n_species; ++s)
        for (int i = 0; i < d.n_sites; ++i) {
            const double yy = ym(i, s);
            const std::size_t idx = (std::size_t) s * d.n_sites + i;
            if (ISNAN(yy)) { d.valid[idx] = 0; continue; }   // missing site x species
            d.y[idx] = yy;
            if (d.family == MSC_POIS) d.lgy[idx] = std::lgamma(yy + 1.0);
        }
    ms_count_nuts_layout(d);
    return d;
}

// Joint log-posterior + gradient over theta = (mu, {z_s}, chol_beta,
// [chol_logr], [log_phi]). NUTS maximises. Mirrors .tobs_ms_count_nuts_logpost.
inline double ms_count_nuts_eval(const MsCountNutsData& d, const double* th,
                                 double sigma_beta, double sigma_logr,
                                 const CommunityCholPri& pr, double* g) {
    const int pb = d.p_beta, S = d.n_species, N = d.n_sites, P = d.P;
    const bool nb = (d.family == MSC_NB), gauss = (d.family == MSC_GAUSS);
    const bool bern = (d.family == MSC_BERN);
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;
    const double* mu     = th + d.mu_off;
    const double* z      = th + d.b_off;
    const double* logphi = gauss ? th + d.logphi_off : nullptr;
    double* g_mu     = g + d.mu_off;
    double* g_z      = g + d.b_off;
    double* g_logphi = gauss ? g + d.logphi_off : nullptr;

    std::vector<double> C_beta;
    chol_unpack_cpp(th + d.chol_beta_off, pb, C_beta);
    const double C_lr = nb ? chol_diag_exp(th[d.chol_logr_off]) : 0.0;

    std::vector<double> A_beta((std::size_t) pb * pb, 0.0);   // A[i,j] = sum_s gb_i z_j
    double A_lr = 0.0;

    // per-species partials for a deterministic (thread-count independent) reduction
    std::vector<double> gmu_s((std::size_t) S * P, 0.0), lp_s(S, 0.0);
    std::vector<double> Abeta_s((std::size_t) S * pb * pb, 0.0);
    std::vector<double> Alr_s(S, 0.0), glogphi_s(gauss ? S : 0, 0.0);

    // Scratch sized once per thread and reused across species (widths are known
    // before the loop). b_beta is fully overwritten each pass; gbb accumulates
    // and is re-zeroed explicitly.
    #pragma omp parallel
    {
    std::vector<double> b_beta(pb), gbb(pb);

    #pragma omp for schedule(static)
    for (int s = 0; s < S; ++s) {
        const double* z_s = z + s * P;
        const double* zb  = z_s;               // length pb
        const double  zr  = nb ? z_s[pb] : 0.0;
        std::fill(gbb.begin(), gbb.end(), 0.0);
        for (int i = 0; i < pb; ++i) {         // b_beta = C_beta z (lower-tri C)
            double v = 0.0;
            for (int j = 0; j <= i; ++j) v += C_beta[(std::size_t) i * pb + j] * zb[j];
            b_beta[i] = v;
        }
        const double log_r_s = nb ? mu[pb] + C_lr * zr : 0.0;
        const double r   = nb    ? std::exp(log_r_s < 30.0 ? log_r_s : 30.0) : 0.0;
        const double phi = gauss ? std::exp(logphi[s] < 30.0 ? logphi[s] : 30.0) : 0.0;
        const double* ys = &d.y[(std::size_t) s * N];
        const double* lg = (d.family == MSC_POIS) ? &d.lgy[(std::size_t) s * N] : nullptr;
        const unsigned char* vs = &d.valid[(std::size_t) s * N];
        double* gmu_loc = &gmu_s[(std::size_t) s * P];
        double lp_loc = 0.0, gblr = 0.0, glp = 0.0;
        for (int i = 0; i < N; ++i) {
            if (!vs[i]) continue;              // skip missing site x species
            double eta = 0.0;
            for (int k = 0; k < pb; ++k) eta += d.X(i, k) * (mu[k] + b_beta[k]);
            double ge;                         // d log p / d eta
            if (nb) {
                double m = std::exp(eta < kExpArgBound ? eta : kExpArgBound);
                if (m < kMinCountMean) m = kMinCountMean;
                lp_loc += R::dnbinom_mu(ys[i], r, m, 1);
                ge = r * (ys[i] - m) / (r + m);
                const double dLL_dr = R::digamma(ys[i] + r) - R::digamma(r)
                                    + std::log(r / (r + m)) + 1.0 - (ys[i] + r) / (r + m);
                gblr += r * dLL_dr;            // chain rule log_r -> r
            } else if (gauss) {
                const double resid = ys[i] - eta;
                lp_loc += R::dnorm(ys[i], eta, std::sqrt(phi), 1);
                ge  = resid / phi;
                glp += -0.5 + resid * resid / (2.0 * phi);   // d log p / d log_phi
            } else if (bern) {                 // Bernoulli (logit), jsdm()
                // log p = y eta - log(1 + exp(eta)), the log-sum-exp computed in
                // the numerically stable branch; d log p / d eta = y - plogis(eta).
                const double lse = (eta > 0.0) ? eta + std::log1p(std::exp(-eta))
                                               : std::log1p(std::exp(eta));
                lp_loc += ys[i] * eta - lse;
                ge = ys[i] - 1.0 / (1.0 + std::exp(-eta));
            } else {                           // Poisson
                const double lam = std::exp(eta < kExpArgBound ? eta : kExpArgBound);
                lp_loc += ys[i] * eta - lam - lg[i];
                ge = ys[i] - lam;
            }
            for (int k = 0; k < pb; ++k) {
                const double gx = ge * d.X(i, k);
                gmu_loc[k] += gx; gbb[k] += gx;
            }
        }
        lp_s[s] = lp_loc;
        double* gz_s = g_z + s * P;            // z grad (beta, data) = C_beta' gbb
        for (int v = 0; v < pb; ++v) {
            double sg = 0.0;
            for (int i = v; i < pb; ++i) sg += C_beta[(std::size_t) i * pb + v] * gbb[i];
            gz_s[v] += sg;
        }
        double* Ab = &Abeta_s[(std::size_t) s * pb * pb];   // A_beta = gbb z'
        for (int i = 0; i < pb; ++i)
            for (int j = 0; j <= i; ++j) Ab[(std::size_t) i * pb + j] = gbb[i] * zb[j];
        if (nb) {
            gmu_loc[pb] += gblr;               // d/d mu_log_r
            gz_s[pb]    += C_lr * gblr;         // z grad (logr, data) = C_lr gblr
            Alr_s[s]     = gblr * zr;
        }
        if (gauss) glogphi_s[s] = glp;
    }
    }  // omp parallel

    double lp = 0.0;                           // serial reduction (species order)
    for (int s = 0; s < S; ++s) {
        const double* gmu_loc = &gmu_s[(std::size_t) s * P];
        for (int k = 0; k < P; ++k) g_mu[k] += gmu_loc[k];
        lp += lp_s[s];
        const double* Ab = &Abeta_s[(std::size_t) s * pb * pb];
        for (int t = 0; t < pb * pb; ++t) A_beta[t] += Ab[t];
        if (nb) A_lr += Alr_s[s];
        if (gauss) g_logphi[s] += glogphi_s[s];
    }

    for (int j = 0; j < S * P; ++j) {          // z prior N(0, I)
        const double zz = z[j]; g_z[j] -= zz; lp += -0.5 * zz * zz;
    }

    // chol_beta coords: data gradient (b = C z) + log-Cholesky hyperprior.
    lp += chol_data_grad_noncentered(A_beta, C_beta, pb, th + d.chol_beta_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_beta_off);
    if (nb) {                                  // scalar log-dispersion covariance
        const std::vector<double> A1(1, A_lr), C1(1, C_lr);
        lp += chol_data_grad_noncentered(A1, C1, 1, th + d.chol_logr_off,
                                         pr.logdiag_mean, pr.logdiag_sd,
                                         pr.offdiag_sd, g + d.chol_logr_off);
    }
    if (gauss) {                               // free log_phi prior N(mean, sd^2)
        const double inv = 1.0 / (d.logphi_sd * d.logphi_sd);
        for (int s = 0; s < S; ++s) {
            const double dd = logphi[s] - d.logphi_mean;
            g_logphi[s] -= inv * dd;
            lp += -0.5 * inv * dd * dd;
        }
    }

    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < pb; ++k) {             // community-mean prior
        g_mu[k] -= inv_sb2 * mu[k]; lp += -0.5 * inv_sb2 * mu[k] * mu[k];
    }
    if (nb) {                                  // log-dispersion community-mean prior
        const double inv_sl2 = 1.0 / (sigma_logr * sigma_logr);
        const double m = mu[pb];
        g_mu[pb] -= inv_sl2 * m; lp += -0.5 * inv_sl2 * m * m;
    }
    return lp;
}

struct MsCountNutsModel {
    MsCountNutsData d;
    double sigma_beta = 10.0, sigma_logr = 1.5;
    CommunityCholPri pr;
};

inline void ms_count_nuts_full_grad(const std::vector<double>& params,
                                    const tulpa::ModelData& data,
                                    const tulpa::ParamLayout& /*layout*/,
                                    std::vector<double>& grad, double* log_post_out) {
    const MsCountNutsModel* m =
        static_cast<const MsCountNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = ms_count_nuts_eval(m->d, params.data(), m->sigma_beta,
                                         m->sigma_logr, m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

} // namespace tulpaObs

// [[Rcpp::export]]
Rcpp::List cpp_ms_count_nuts_joint_logpost(Rcpp::List spec, Rcpp::NumericVector theta,
                                           Rcpp::List pri, double sigma_beta,
                                           double sigma_logr) {
    tulpaObs::MsCountNutsData d = tulpaObs::ms_count_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    tulpaObs::CommunityCholPri pr = tulpaObs::community_chol_pri_from_list(pri);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::ms_count_nuts_eval(
        d, theta.begin(), sigma_beta, sigma_logr, pr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// [[Rcpp::export]]
Rcpp::List cpp_ms_count_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                             Rcpp::List pri, double sigma_beta, double sigma_logr,
                             Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                             int n_iter, int n_warmup, int max_treedepth,
                             double adapt_delta, int seed, bool verbose) {
    tulpaObs::MsCountNutsModel m;
    m.d = tulpaObs::ms_count_nuts_build_data(spec);
    m.sigma_beta = sigma_beta; m.sigma_logr = sigma_logr;
    m.pr = tulpaObs::community_chol_pri_from_list(pri);
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    return tulpaObs::run_tulpa_nuts(
        &tulpaObs::ms_count_nuts_full_grad, &m, m.d.total,
        theta0, sigma_beta, inv_metric,
        n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose,
        "ms_count", m.d.n_sites);
}
