// ms_int_occu_nuts.cpp
// C++ joint log-posterior + gradient for the community / multispecies INTEGRATED
// occupancy (ms_int_occu()) NUTS target. The R reference
// .tobs_ms_int_occu_nuts_logpost (R/ms_int_occu_nuts.R) is the oracle; this port
// mirrors it and is cross-checked against it before driving tulpa's NUTS engine.
//
// The multi-source generalisation of the community single-season occupancy target
// (ms_occu_nuts.cpp): per-species two-state occupancy with D detection SOURCES,
// each source a plain per-species random-effect detection arm with its own
// Gaussian community covariance. There are NO shared globals (unlike the dynamic
// family's gamma / eps). The latent z integrates out per (species, site) in closed
// form -- a detected site (any source detects) forces z = 1; a no-detection site
// pools all sources' occupied-undetected mass:
//   L_i = psi_i * prod_d (1 - p_{d,i})^{n_valid_{d,i}} + (1 - psi_i).
//
// NON-CENTERED parameterisation (mirrors ms_occu_nuts.cpp): the per-species block
// carries the whitened z_s, the deviation is b_{s,arm} = C_arm z_{s,arm} per arm
// (psi + D detection), so each community covariance leaves the b-prior and enters
// ONLY the data term. The target factorises as
//   log p = sum_{s,i} log L_{s,i}(theta)
//         - 0.5 ||mu_coef||^2 / sigma.beta^2
//         - 0.5 sum_s ||z_s||^2
//         + sum_arm log p(Sigma_arm coords).
// The multi-source marginal + its eta-gradients are computed inline here (matching
// .ms_int_occu_sp_ll / _grad); the chol gradient flows from the data term via
// b = C z (chol_data_grad_noncentered) plus the coordinate hyperprior; the arm
// Cholesky factors use the shared helpers in community_chol.h.

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

static inline double mio_sigmoid_(double x) { return 1.0 / (1.0 + std::exp(-x)); }
static inline double mio_clamp_(double p) {
    return std::min(std::max(p, 1e-12), 1.0 - 1e-12);
}

struct MsIntOccuPri {
    double logdiag_mean = 0.0, logdiag_sd = 1.5, offdiag_sd = 1.0;
};

// Per-(species, source) site summary: n_valid_{d,i} visits, n_det_{d,i} detections
// (each length n_sites). `any_det` is per-site pooled over sources.
struct MioSpecSumm {
    std::vector<std::vector<int> > n_valid;   // [D][n_sites]
    std::vector<std::vector<int> > n_det;     // [D][n_sites]
    std::vector<char> any_det;                // [n_sites]
};

struct MsIntOccuNutsData {
    int n_sites = 0, n_species = 0, D = 0, p_psi = 0;
    NumericMatrix X_psi;                       // n_sites x p_psi
    std::vector<NumericMatrix> X_p;            // [D] each n_sites x p_p_d
    std::vector<int> p_p;                      // [D]
    std::vector<MioSpecSumm> summ;             // [s]

    // Packed-coordinate layout (mirrors .tobs_ms_int_occu_nuts_layout). Arms are
    // psi + the D detection sources; z blocks are species-major, width P.
    int P_tot = 0, q_psi = 0, total = 0;
    int mu_off = 0, b_off = 0, chol_psi_off = 0;
    std::vector<int> p_off;                    // [D] within-P offset of source d
    std::vector<int> q_p;                      // [D] chol length of source d
    std::vector<int> chol_p_off;               // [D] packed offset of source d chol
};

inline void ms_int_occu_nuts_layout(MsIntOccuNutsData& d) {
    d.P_tot = d.p_psi;
    d.p_off.assign(d.D, 0);
    for (int s = 0; s < d.D; ++s) { d.p_off[s] = d.P_tot; d.P_tot += d.p_p[s]; }
    d.q_psi = d.p_psi * (d.p_psi + 1) / 2;
    d.mu_off = 0;
    d.b_off  = d.P_tot;
    int coff = d.P_tot + d.n_species * d.P_tot;
    d.chol_psi_off = coff; coff += d.q_psi;
    d.q_p.assign(d.D, 0);
    d.chol_p_off.assign(d.D, 0);
    for (int s = 0; s < d.D; ++s) {
        d.q_p[s] = d.p_p[s] * (d.p_p[s] + 1) / 2;
        d.chol_p_off[s] = coff; coff += d.q_p[s];
    }
    d.total = coff;
}

inline MsIntOccuNutsData ms_int_occu_nuts_build_data(const Rcpp::List& spec) {
    MsIntOccuNutsData d;
    d.X_psi     = Rcpp::as<NumericMatrix>(spec["X_psi"]);
    d.n_sites   = Rcpp::as<int>(spec["n_sites"]);
    d.n_species = Rcpp::as<int>(spec["n_species"]);
    d.D         = Rcpp::as<int>(spec["D"]);
    d.p_psi     = d.X_psi.ncol();

    Rcpp::List Xp = Rcpp::as<Rcpp::List>(spec["X_p"]);   // list of D designs
    Rcpp::List nv = Rcpp::as<Rcpp::List>(spec["n_valid"]);  // list of D int matrices
    Rcpp::List nd = Rcpp::as<Rcpp::List>(spec["n_det"]);
    d.X_p.reserve(d.D); d.p_p.assign(d.D, 0);
    std::vector<IntegerMatrix> NV, ND; NV.reserve(d.D); ND.reserve(d.D);
    for (int s = 0; s < d.D; ++s) {
        NumericMatrix Xd = Rcpp::as<NumericMatrix>(Xp[s]);
        d.X_p.push_back(Xd);
        d.p_p[s] = Xd.ncol();
        IntegerMatrix vd = Rcpp::as<IntegerMatrix>(nv[s]);
        IntegerMatrix dd = Rcpp::as<IntegerMatrix>(nd[s]);
        if (vd.nrow() != d.n_sites || vd.ncol() != d.n_species)
            Rcpp::stop("n_valid[[d]] must be n_sites x n_species");
        if (dd.nrow() != d.n_sites || dd.ncol() != d.n_species)
            Rcpp::stop("n_det[[d]] must be n_sites x n_species");
        NV.push_back(vd); ND.push_back(dd);
    }

    d.summ.assign(d.n_species, MioSpecSumm());
    for (int s = 0; s < d.n_species; ++s) {
        MioSpecSumm& su = d.summ[s];
        su.n_valid.assign(d.D, std::vector<int>((std::size_t) d.n_sites, 0));
        su.n_det.assign(d.D,   std::vector<int>((std::size_t) d.n_sites, 0));
        su.any_det.assign((std::size_t) d.n_sites, 0);
        for (int i = 0; i < d.n_sites; ++i) {
            int tot_det = 0;
            for (int dd = 0; dd < d.D; ++dd) {
                su.n_valid[dd][i] = NV[dd](i, s);
                su.n_det[dd][i]   = ND[dd](i, s);
                tot_det += ND[dd](i, s);
            }
            su.any_det[i] = (tot_det > 0) ? 1 : 0;
        }
    }
    ms_int_occu_nuts_layout(d);
    return d;
}

// Joint log-posterior + gradient over theta = (mu, {z_s}, chol_psi, chol_p1..pD).
// NUTS maximises; returns the log-posterior (no negation) and writes grad into g.
inline double ms_int_occu_nuts_eval(const MsIntOccuNutsData& d, const double* th,
                                    double sigma_beta, const MsIntOccuPri& pr,
                                    double* g) {
    const int P = d.P_tot, S = d.n_species, p_psi = d.p_psi, D = d.D;
    const int n_sites = d.n_sites;
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;

    const double* mu = th + d.mu_off;
    const double* z  = th + d.b_off;
    double* g_mu = g + d.mu_off;
    double* g_z  = g + d.b_off;

    // Cholesky factors per arm (row-major, lower-triangular).
    std::vector<double> C_psi;
    chol_unpack_cpp(th + d.chol_psi_off, p_psi, C_psi);
    std::vector<std::vector<double> > C_p(D);
    for (int dd = 0; dd < D; ++dd)
        chol_unpack_cpp(th + d.chol_p_off[dd], d.p_p[dd], C_p[dd]);

    // chol data-gradient accumulators A_arm[i,j] = sum_s grad_b_{s,i} z_{s,j}.
    std::vector<double> A_psi((std::size_t) p_psi * p_psi, 0.0);
    std::vector<std::vector<double> > A_p(D);
    for (int dd = 0; dd < D; ++dd)
        A_p[dd].assign((std::size_t) d.p_p[dd] * d.p_p[dd], 0.0);

    double lp = 0.0;
    // Serial species loop (byte-exact vs the R oracle; sizes are modest).
    for (int s = 0; s < S; ++s) {
        const double* z_s = z + s * P;
        const MioSpecSumm& su = d.summ[s];

        // reconstruct b_psi = C_psi z_psi and eta_psi.
        std::vector<double> b_psi(p_psi, 0.0);
        for (int i = 0; i < p_psi; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j)
                v += C_psi[(std::size_t) i * p_psi + j] * z_s[j];
            b_psi[i] = v;
        }
        std::vector<double> eta_psi(n_sites, 0.0);
        for (int i = 0; i < n_sites; ++i) {
            double e = 0.0;
            for (int k = 0; k < p_psi; ++k)
                e += d.X_psi(i, k) * (mu[k] + b_psi[k]);
            eta_psi[i] = e;
        }
        // per-source b_pd + eta_pd.
        std::vector<std::vector<double> > eta_p(D), b_p(D);
        for (int dd = 0; dd < D; ++dd) {
            const int ppd = d.p_p[dd];
            const double* z_pd = z_s + d.p_off[dd];
            b_p[dd].assign(ppd, 0.0);
            for (int i = 0; i < ppd; ++i) {
                double v = 0.0;
                for (int j = 0; j <= i; ++j)
                    v += C_p[dd][(std::size_t) i * ppd + j] * z_pd[j];
                b_p[dd][i] = v;
            }
            eta_p[dd].assign(n_sites, 0.0);
            for (int i = 0; i < n_sites; ++i) {
                double e = 0.0;
                for (int k = 0; k < ppd; ++k)
                    e += d.X_p[dd](i, k) * (mu[d.p_off[dd] + k] + b_p[dd][k]);
                eta_p[dd][i] = e;
            }
        }

        // ---- multi-source two-state marginal + eta-gradients (per site) ----
        std::vector<double> g_eta_psi(n_sites, 0.0);
        std::vector<std::vector<double> > g_eta_pd(D, std::vector<double>(n_sites, 0.0));
        for (int i = 0; i < n_sites; ++i) {
            const double psi = mio_clamp_(mio_sigmoid_(eta_psi[i]));
            // per-source p_d + log(1-p_d) accumulation.
            double log_undet = 0.0, log_det_term = 0.0;
            std::vector<double> pd(D);
            for (int dd = 0; dd < D; ++dd) {
                const double p = mio_clamp_(mio_sigmoid_(eta_p[dd][i]));
                pd[dd] = p;
                const int nvd = su.n_valid[dd][i];
                const int ndd = su.n_det[dd][i];
                log_undet    += (double) nvd * std::log1p(-p);
                log_det_term += (double) ndd * std::log(p)
                              + (double) (nvd - ndd) * std::log1p(-p);
            }
            if (su.any_det[i]) {
                lp += std::log(psi) + log_det_term;
                g_eta_psi[i] = 1.0 - psi;
                for (int dd = 0; dd < D; ++dd)
                    g_eta_pd[dd][i] = (double) su.n_det[dd][i]
                                    - (double) su.n_valid[dd][i] * pd[dd];
            } else {
                const double prodterm = std::exp(log_undet);
                const double A = psi * prodterm;
                const double L = A + (1.0 - psi);
                lp += std::log(L);
                g_eta_psi[i] = psi * (1.0 - psi) * (prodterm - 1.0) / L;
                for (int dd = 0; dd < D; ++dd)
                    g_eta_pd[dd][i] = -(A / L) * (double) su.n_valid[dd][i] * pd[dd];
            }
        }

        // ---- design-sandwich the eta-gradients into mu / b, per arm ----
        double* gmu = g_mu;
        double* gz_s = g_z + s * P;
        // psi arm.
        std::vector<double> gbpsi(p_psi, 0.0);
        for (int i = 0; i < n_sites; ++i) {
            const double gp = g_eta_psi[i];
            for (int k = 0; k < p_psi; ++k) {
                const double gx = gp * d.X_psi(i, k);
                gmu[k]    += gx;
                gbpsi[k]  += gx;
            }
        }
        for (int v = 0; v < p_psi; ++v) {
            double sg = 0.0;
            for (int i = v; i < p_psi; ++i)
                sg += C_psi[(std::size_t) i * p_psi + v] * gbpsi[i];
            gz_s[v] += sg;
        }
        for (int i = 0; i < p_psi; ++i)
            for (int j = 0; j <= i; ++j)
                A_psi[(std::size_t) i * p_psi + j] += gbpsi[i] * z_s[j];
        // detection arms.
        for (int dd = 0; dd < D; ++dd) {
            const int ppd = d.p_p[dd];
            const int poff = d.p_off[dd];
            const double* z_pd = z_s + poff;
            std::vector<double> gbpd(ppd, 0.0);
            for (int i = 0; i < n_sites; ++i) {
                const double gp = g_eta_pd[dd][i];
                for (int k = 0; k < ppd; ++k) {
                    const double gx = gp * d.X_p[dd](i, k);
                    gmu[poff + k] += gx;
                    gbpd[k]       += gx;
                }
            }
            double* gz_pd = gz_s + poff;
            for (int v = 0; v < ppd; ++v) {
                double sg = 0.0;
                for (int i = v; i < ppd; ++i)
                    sg += C_p[dd][(std::size_t) i * ppd + v] * gbpd[i];
                gz_pd[v] += sg;
            }
            for (int i = 0; i < ppd; ++i)
                for (int j = 0; j <= i; ++j)
                    A_p[dd][(std::size_t) i * ppd + j] += gbpd[i] * z_pd[j];
        }
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
    for (int dd = 0; dd < D; ++dd)
        lp += chol_data_grad_noncentered(A_p[dd], C_p[dd], d.p_p[dd],
                                         th + d.chol_p_off[dd], pr.logdiag_mean,
                                         pr.logdiag_sd, pr.offdiag_sd,
                                         g + d.chol_p_off[dd]);

    // ---- community-mean priors: N(0, sigma.beta^2) on every coefficient mean ----
    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < P; ++k) {
        g_mu[k] -= inv_sb2 * mu[k];
        lp      += -0.5 * inv_sb2 * mu[k] * mu[k];
    }
    return lp;
}

struct MsIntOccuNutsModel {
    MsIntOccuNutsData d;
    double sigma_beta = 5.0;
    MsIntOccuPri pr;
};

inline void ms_int_occu_nuts_full_grad(const std::vector<double>& params,
                                       const tulpa::ModelData& data,
                                       const tulpa::ParamLayout& /*layout*/,
                                       std::vector<double>& grad,
                                       double* log_post_out) {
    const MsIntOccuNutsModel* m =
        static_cast<const MsIntOccuNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = ms_int_occu_nuts_eval(m->d, params.data(), m->sigma_beta,
                                            m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

inline MsIntOccuPri ms_int_occu_pri_from_list(const Rcpp::List& pri) {
    MsIntOccuPri pr;
    pr.logdiag_mean = Rcpp::as<double>(pri["chol_logdiag_mean"]);
    pr.logdiag_sd   = Rcpp::as<double>(pri["chol_logdiag_sd"]);
    pr.offdiag_sd   = Rcpp::as<double>(pri["chol_offdiag_sd"]);
    return pr;
}

} // namespace tulpaObs

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .tobs_ms_int_occu_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_ms_int_occu_nuts_joint_logpost(Rcpp::List spec,
                                              Rcpp::NumericVector theta,
                                              Rcpp::List pri, double sigma_beta) {
    tulpaObs::MsIntOccuNutsData d = tulpaObs::ms_int_occu_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    tulpaObs::MsIntOccuPri pr = tulpaObs::ms_int_occu_pri_from_list(pri);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::ms_int_occu_nuts_eval(
        d, theta.begin(), sigma_beta, pr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the community integrated occupancy target via tulpa's engine and
// the FullGradFn. Returns draws + diagnostics.
// [[Rcpp::export]]
Rcpp::List cpp_ms_int_occu_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                                Rcpp::List pri, double sigma_beta,
                                Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                                int n_iter, int n_warmup, int max_treedepth,
                                double adapt_delta, int seed, bool verbose) {
    tulpaObs::MsIntOccuNutsModel m;
    m.d = tulpaObs::ms_int_occu_nuts_build_data(spec);
    m.sigma_beta = sigma_beta;
    m.pr = tulpaObs::ms_int_occu_pri_from_list(pri);
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    return tulpaObs::run_tulpa_nuts(
        &tulpaObs::ms_int_occu_nuts_full_grad, &m, m.d.total,
        theta0, sigma_beta, inv_metric,
        n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose,
        "ms_int_occu", m.d.n_sites);
}
