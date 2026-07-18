// ms_occu_cover_nuts.cpp
// C++ joint log-posterior + gradient for the community / multispecies joint
// occupancy-detection + cover (ms_occu_cover()) NUTS target (the NON-spatial
// per-arm-covariance community model; #115 part B7). The R reference
// .tobs_ms_occu_cover_nuts_logpost (R/ms_occu_cover_nuts.R) is the oracle; this
// port mirrors it and is cross-checked against it before driving tulpa's NUTS.
//
// The joint-cover analogue of the community occupancy targets (ms_occu_nuts.cpp /
// ms_int_occu_nuts.cpp): three non-centered per-species arms (occ + p + pos), each
// b_{s,arm} = C_arm z_{s,arm} with its own log-Cholesky community covariance, plus
// ONE shared community log-dispersion scalar (the positive-arm beta precision /
// lognormal-or-gaussian residual SD, on the log scale) that carries no per-species
// random effect. The per-(species, cell) two-state occu_cover marginal + its
// eta / log-dispersion gradients are the SAME as the single-model occu_cover NUTS
// (occu_cover_nuts.cpp): the occupancy/detection mixture via nodet_mixture_block
// and the cover arm via pos_log_density / pos_grad_eta / pos_grad_logdisp
// (occu_coupling_shared.h). This community port wraps that per-cell body in the
// non-centered per-species loop and adds the community-covariance gradients.
//
// The joint log-posterior is
//   log p = sum_s log L_s(theta)                      # per-species occu_cover marginal
//         - 0.5 ||mu_coef||^2 / sigma.beta^2          # community-mean priors
//         - 0.5 sum_s ||z_s||^2                       # whitened RE prior (N(0,I))
//         + sum_arm log p(Sigma_arm coords)           # log-Cholesky hyperpriors
//         + log p(log_disp)                           # weakly-informative dispersion prior
// under b_{s,arm} = C_arm z_{s,arm}. Detection / cover designs are per-visit (a
// site-level block broadcast across visits + an optional visit-level block,
// site-major), matching .occu_cover_eta_from_par.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <tulpa/model_data.h>
#include <tulpa/param_layout.h>
#include <tulpa/likelihood.h>
#include <tulpa/nuts_api.h>
#include "occu_coupling_shared.h"
#include "community_chol.h"

using namespace Rcpp;

namespace tulpaObs {

struct MsOccuCoverPri {
    double logdiag_mean = 0.0, logdiag_sd = 1.5, offdiag_sd = 1.0;
    double log_disp_mean = 0.0, log_disp_sd = 2.0;
};

// Per-species response matrices [n_sites x max_visits].
struct MsocSpec {
    IntegerMatrix y;      // detection 0/1 (NA -> 0)
    NumericMatrix y_pos;  // cover (non-finite where absent)
    IntegerMatrix valid;  // visit observed 0/1
};

struct MsOccuCoverNutsData {
    int n_sites = 0, max_visits = 0, n_species = 0, pos_code = 0;
    // Shared designs across species.
    NumericMatrix X_occ, X_det_site, X_det_visit, X_pos_site, X_pos_visit;
    std::vector<MsocSpec> sp;                 // [n_species]

    int p_occ = 0, p_det_site = 0, p_det_visit = 0, p_pos_site = 0, p_pos_visit = 0;
    // Arm coefficient widths (site + visit blocks combined per arm).
    int P_occ = 0, P_p = 0, P_pos = 0, P_tot = 0;
    // Within-P arm offsets (occ at 0, p at P_occ, pos at P_occ+P_p).
    int off_occ = 0, off_p = 0, off_pos = 0;
    // Packed layout.
    int b_off = 0, chol_occ_off = 0, chol_p_off = 0, chol_pos_off = 0, log_disp_off = 0;
    int q_occ = 0, q_p = 0, q_pos = 0, total = 0;
};

inline void ms_occu_cover_nuts_layout(MsOccuCoverNutsData& d) {
    d.P_occ = d.p_occ;
    d.P_p   = d.p_det_site + d.p_det_visit;
    d.P_pos = d.p_pos_site + d.p_pos_visit;
    d.P_tot = d.P_occ + d.P_p + d.P_pos;
    d.off_occ = 0; d.off_p = d.P_occ; d.off_pos = d.P_occ + d.P_p;
    d.q_occ = d.P_occ * (d.P_occ + 1) / 2;
    d.q_p   = d.P_p   * (d.P_p   + 1) / 2;
    d.q_pos = d.P_pos * (d.P_pos + 1) / 2;
    d.b_off = d.P_tot;
    int coff = d.P_tot + d.n_species * d.P_tot;
    d.chol_occ_off = coff; coff += d.q_occ;
    d.chol_p_off   = coff; coff += d.q_p;
    d.chol_pos_off = coff; coff += d.q_pos;
    d.log_disp_off = coff; coff += 1;
    d.total = coff;
}

inline MsOccuCoverNutsData ms_occu_cover_nuts_build_data(const Rcpp::List& spec) {
    MsOccuCoverNutsData d;
    d.n_sites    = Rcpp::as<int>(spec["n_sites"]);
    d.max_visits = Rcpp::as<int>(spec["max_visits"]);
    d.n_species  = Rcpp::as<int>(spec["n_species"]);
    d.pos_code   = Rcpp::as<int>(spec["pos_code"]);
    d.X_occ       = Rcpp::as<NumericMatrix>(spec["X_occ"]);
    d.X_det_site  = Rcpp::as<NumericMatrix>(spec["X_det_site"]);
    d.X_det_visit = Rcpp::as<NumericMatrix>(spec["X_det_visit"]);
    d.X_pos_site  = Rcpp::as<NumericMatrix>(spec["X_pos_site"]);
    d.X_pos_visit = Rcpp::as<NumericMatrix>(spec["X_pos_visit"]);
    d.p_occ       = d.X_occ.ncol();
    d.p_det_site  = d.X_det_site.ncol();
    d.p_det_visit = d.X_det_visit.ncol();
    d.p_pos_site  = d.X_pos_site.ncol();
    d.p_pos_visit = d.X_pos_visit.ncol();

    Rcpp::List Y  = Rcpp::as<Rcpp::List>(spec["y"]);      // list of n_species
    Rcpp::List YP = Rcpp::as<Rcpp::List>(spec["y_pos"]);
    Rcpp::List V  = Rcpp::as<Rcpp::List>(spec["valid"]);
    d.sp.reserve(d.n_species);
    for (int s = 0; s < d.n_species; ++s) {
        MsocSpec ms;
        ms.y     = Rcpp::as<IntegerMatrix>(Y[s]);
        ms.y_pos = Rcpp::as<NumericMatrix>(YP[s]);
        ms.valid = Rcpp::as<IntegerMatrix>(V[s]);
        d.sp.push_back(ms);
    }
    ms_occu_cover_nuts_layout(d);
    return d;
}

// Joint log-posterior + gradient over theta = (mu, {z_s}, chol_occ, chol_p,
// chol_pos, log_disp). NUTS maximises; returns the log-posterior (no negation)
// and writes the gradient into g. Mirrors .tobs_ms_occu_cover_nuts_logpost.
inline double ms_occu_cover_nuts_eval(const MsOccuCoverNutsData& d, const double* th,
                                      double sigma_beta, const MsOccuCoverPri& pr,
                                      double* g) {
    const int N = d.n_sites, J = d.max_visits, S = d.n_species, P = d.P_tot;
    const int P_occ = d.P_occ, P_p = d.P_p, P_pos = d.P_pos;
    const int pds = d.p_det_site, pdv = d.p_det_visit;
    const int pps = d.p_pos_site, ppv = d.p_pos_visit;
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;

    const double* mu = th;                        // length P
    const double* z  = th + d.b_off;              // species-major z
    const double  log_disp = th[d.log_disp_off];
    const double  disp = std::exp(log_disp);
    double* g_mu = g;
    double* g_z  = g + d.b_off;

    std::vector<double> C_occ, C_p, C_pos;
    chol_unpack_cpp(th + d.chol_occ_off, P_occ, C_occ);
    chol_unpack_cpp(th + d.chol_p_off,   P_p,   C_p);
    chol_unpack_cpp(th + d.chol_pos_off, P_pos, C_pos);

    std::vector<double> A_occ((std::size_t) P_occ * P_occ, 0.0),
                        A_p((std::size_t) P_p * P_p, 0.0),
                        A_pos((std::size_t) P_pos * P_pos, 0.0);

    double lp = 0.0, g_logdisp = 0.0;
    std::vector<double> eta_p(J), g_eta_p(J), g_eta_pos(J), eta_pc(J), g_pc(J);

    for (int s = 0; s < S; ++s) {
        const double* z_s   = z + s * P;
        const double* z_occ = z_s + d.off_occ;
        const double* z_p   = z_s + d.off_p;
        const double* z_pos = z_s + d.off_pos;
        const MsocSpec& ms  = d.sp[s];

        // reconstruct b = mu + C z per arm.
        std::vector<double> b_occ(P_occ, 0.0), b_p(P_p, 0.0), b_pos(P_pos, 0.0);
        for (int i = 0; i < P_occ; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j) v += C_occ[(std::size_t) i * P_occ + j] * z_occ[j];
            b_occ[i] = mu[d.off_occ + i] + v;
        }
        for (int i = 0; i < P_p; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j) v += C_p[(std::size_t) i * P_p + j] * z_p[j];
            b_p[i] = mu[d.off_p + i] + v;
        }
        for (int i = 0; i < P_pos; ++i) {
            double v = 0.0;
            for (int j = 0; j <= i; ++j) v += C_pos[(std::size_t) i * P_pos + j] * z_pos[j];
            b_pos[i] = mu[d.off_pos + i] + v;
        }
        const double* bp_site   = b_p.data();
        const double* bp_visit  = b_p.data() + pds;
        const double* bpos_site = b_pos.data();
        const double* bpos_visit= b_pos.data() + pps;

        // per-arm b-space coefficient gradients (design-sandwiched).
        std::vector<double> gbo(P_occ, 0.0), gbp(P_p, 0.0), gbpos(P_pos, 0.0);

        for (int i = 0; i < N; ++i) {
            double eta_psi = 0.0;
            for (int k = 0; k < P_occ; ++k) eta_psi += d.X_occ(i, k) * b_occ[k];
            const double psi = sigmoid_(eta_psi);

            double eta_p_site = 0.0;
            for (int k = 0; k < pds; ++k) eta_p_site += d.X_det_site(i, k) * bp_site[k];
            bool any_det = false;
            for (int v = 0; v < J; ++v) {
                g_eta_p[v] = 0.0; g_eta_pos[v] = 0.0;
                if (ms.valid(i, v) == 0) { eta_p[v] = 0.0; continue; }
                double e = eta_p_site;
                if (pdv > 0) {
                    const int row = i * J + v;
                    for (int k = 0; k < pdv; ++k) e += d.X_det_visit(row, k) * bp_visit[k];
                }
                eta_p[v] = e;
                if (ms.y(i, v) == 1) any_det = true;
            }

            double g_eta_psi = 0.0;
            if (any_det) {
                lp += log_safe_(psi);
                g_eta_psi = 1.0 - psi;
                for (int v = 0; v < J; ++v) {
                    if (ms.valid(i, v) == 0) continue;
                    const double pv = sigmoid_(eta_p[v]);
                    if (ms.y(i, v) == 1) { lp += log_safe_(pv);       g_eta_p[v] = 1.0 - pv; }
                    else                 { lp += log_safe_(1.0 - pv); g_eta_p[v] = -pv; }
                }
                for (int v = 0; v < J; ++v) {
                    if (ms.valid(i, v) == 0 || ms.y(i, v) != 1) continue;
                    const double yp = ms.y_pos(i, v);
                    if (!std::isfinite(yp)) continue;
                    double eta_pos = 0.0;
                    for (int k = 0; k < pps; ++k) eta_pos += d.X_pos_site(i, k) * bpos_site[k];
                    if (ppv > 0) {
                        const int row = i * J + v;
                        for (int k = 0; k < ppv; ++k) eta_pos += d.X_pos_visit(row, k) * bpos_visit[k];
                    }
                    lp += pos_log_density(d.pos_code, yp, eta_pos, disp);
                    g_eta_pos[v] = pos_grad_eta(d.pos_code, yp, eta_pos, disp);
                    g_logdisp   += pos_grad_logdisp(d.pos_code, yp, eta_pos, disp);
                }
            } else {
                int nv = 0;
                for (int v = 0; v < J; ++v) {
                    if (ms.valid(i, v) == 0) continue;
                    eta_pc[nv] = eta_p[v]; ++nv;
                }
                double g_w = 0.0, nh_w = 0.0;
                const double cell_ll = nodet_mixture_block(
                    psi, eta_pc.data(), nv, false, false,
                    g_w, nh_w, g_pc.data(), nullptr, nullptr, nullptr);
                lp += cell_ll;
                g_eta_psi = g_w;
                int j = 0;
                for (int v = 0; v < J; ++v) {
                    if (ms.valid(i, v) == 0) continue;
                    g_eta_p[v] = g_pc[j]; ++j;
                }
            }

            // design-sandwich onto the per-arm b-space coefficient gradients.
            for (int k = 0; k < P_occ; ++k) gbo[k] += g_eta_psi * d.X_occ(i, k);
            double g_eta_p_sum = 0.0, g_eta_pos_sum = 0.0;
            for (int v = 0; v < J; ++v) { g_eta_p_sum += g_eta_p[v]; g_eta_pos_sum += g_eta_pos[v]; }
            for (int k = 0; k < pds; ++k) gbp[k]           += g_eta_p_sum   * d.X_det_site(i, k);
            for (int k = 0; k < pps; ++k) gbpos[k]         += g_eta_pos_sum * d.X_pos_site(i, k);
            if (pdv > 0) {
                for (int v = 0; v < J; ++v) {
                    if (g_eta_p[v] == 0.0) continue;
                    const int row = i * J + v;
                    for (int k = 0; k < pdv; ++k) gbp[pds + k] += g_eta_p[v] * d.X_det_visit(row, k);
                }
            }
            if (ppv > 0) {
                for (int v = 0; v < J; ++v) {
                    if (g_eta_pos[v] == 0.0) continue;
                    const int row = i * J + v;
                    for (int k = 0; k < ppv; ++k) gbpos[pps + k] += g_eta_pos[v] * d.X_pos_visit(row, k);
                }
            }
        }

        // ---- non-centered: g_mu += g_b; g_z = C' g_b; A_arm += g_b z' ----
        double* gz_s = g_z + s * P;
        // occ arm
        for (int k = 0; k < P_occ; ++k) g_mu[d.off_occ + k] += gbo[k];
        for (int vc = 0; vc < P_occ; ++vc) {
            double sg = 0.0;
            for (int i = vc; i < P_occ; ++i) sg += C_occ[(std::size_t) i * P_occ + vc] * gbo[i];
            gz_s[d.off_occ + vc] += sg;
        }
        for (int i = 0; i < P_occ; ++i)
            for (int j = 0; j <= i; ++j) A_occ[(std::size_t) i * P_occ + j] += gbo[i] * z_occ[j];
        // p arm
        for (int k = 0; k < P_p; ++k) g_mu[d.off_p + k] += gbp[k];
        for (int vc = 0; vc < P_p; ++vc) {
            double sg = 0.0;
            for (int i = vc; i < P_p; ++i) sg += C_p[(std::size_t) i * P_p + vc] * gbp[i];
            gz_s[d.off_p + vc] += sg;
        }
        for (int i = 0; i < P_p; ++i)
            for (int j = 0; j <= i; ++j) A_p[(std::size_t) i * P_p + j] += gbp[i] * z_p[j];
        // pos arm
        for (int k = 0; k < P_pos; ++k) g_mu[d.off_pos + k] += gbpos[k];
        for (int vc = 0; vc < P_pos; ++vc) {
            double sg = 0.0;
            for (int i = vc; i < P_pos; ++i) sg += C_pos[(std::size_t) i * P_pos + vc] * gbpos[i];
            gz_s[d.off_pos + vc] += sg;
        }
        for (int i = 0; i < P_pos; ++i)
            for (int j = 0; j <= i; ++j) A_pos[(std::size_t) i * P_pos + j] += gbpos[i] * z_pos[j];
    }

    // ---- z prior: standard normal over the whole per-species block ----
    for (int j = 0; j < S * P; ++j) {
        const double zz = z[j];
        g_z[j] -= zz;
        lp     += -0.5 * zz * zz;
    }

    // ---- chol coords: data gradient (via b = C z) + coordinate hyperprior ----
    lp += chol_data_grad_noncentered(A_occ, C_occ, P_occ, th + d.chol_occ_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_occ_off);
    lp += chol_data_grad_noncentered(A_p, C_p, P_p, th + d.chol_p_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_p_off);
    lp += chol_data_grad_noncentered(A_pos, C_pos, P_pos, th + d.chol_pos_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_pos_off);

    // ---- community-mean priors: N(0, sigma.beta^2) on every coefficient mean ----
    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < P; ++k) {
        g_mu[k] -= inv_sb2 * mu[k];
        lp      += -0.5 * inv_sb2 * mu[k] * mu[k];
    }

    // ---- shared log-dispersion prior ----
    const double ld = log_disp - pr.log_disp_mean;
    const double ild2 = 1.0 / (pr.log_disp_sd * pr.log_disp_sd);
    lp += -0.5 * ild2 * ld * ld;
    g[d.log_disp_off] = g_logdisp - ild2 * ld;

    return lp;
}

struct MsOccuCoverNutsModel {
    MsOccuCoverNutsData d;
    double sigma_beta = 5.0;
    MsOccuCoverPri pr;
};

inline void ms_occu_cover_nuts_full_grad(const std::vector<double>& params,
                                         const tulpa::ModelData& data,
                                         const tulpa::ParamLayout& /*layout*/,
                                         std::vector<double>& grad,
                                         double* log_post_out) {
    const MsOccuCoverNutsModel* m =
        static_cast<const MsOccuCoverNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->d.total, 0.0);
    const double lp = ms_occu_cover_nuts_eval(m->d, params.data(), m->sigma_beta,
                                              m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

inline MsOccuCoverPri ms_occu_cover_pri_from_list(const Rcpp::List& pri) {
    MsOccuCoverPri pr;
    pr.logdiag_mean  = Rcpp::as<double>(pri["chol_logdiag_mean"]);
    pr.logdiag_sd    = Rcpp::as<double>(pri["chol_logdiag_sd"]);
    pr.offdiag_sd    = Rcpp::as<double>(pri["chol_offdiag_sd"]);
    if (pri.containsElementNamed("log_disp_mean"))
        pr.log_disp_mean = Rcpp::as<double>(pri["log_disp_mean"]);
    if (pri.containsElementNamed("log_disp_sd"))
        pr.log_disp_sd = Rcpp::as<double>(pri["log_disp_sd"]);
    return pr;
}

} // namespace tulpaObs

// Full-vector joint log-posterior + gradient, the Rcpp cross-check entry for the
// R oracle .tobs_ms_occu_cover_nuts_logpost.
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_cover_nuts_joint_logpost(Rcpp::List spec,
                                                Rcpp::NumericVector theta,
                                                Rcpp::List pri, double sigma_beta) {
    tulpaObs::MsOccuCoverNutsData d = tulpaObs::ms_occu_cover_nuts_build_data(spec);
    if ((int) theta.size() != d.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), d.total);
    tulpaObs::MsOccuCoverPri pr = tulpaObs::ms_occu_cover_pri_from_list(pri);
    Rcpp::NumericVector grad(d.total);
    const double lp = tulpaObs::ms_occu_cover_nuts_eval(
        d, theta.begin(), sigma_beta, pr, grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// Run NUTS on the community joint occu+cover target via tulpa's engine and the
// in-tree FullGradFn. Returns draws + diagnostics.
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_cover_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                                  Rcpp::List pri, double sigma_beta,
                                  Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                                  int n_iter, int n_warmup, int max_treedepth,
                                  double adapt_delta, int seed, bool verbose) {
    tulpaObs::MsOccuCoverNutsModel m;
    m.d = tulpaObs::ms_occu_cover_nuts_build_data(spec);
    m.sigma_beta = sigma_beta;
    m.pr = tulpaObs::ms_occu_cover_pri_from_list(pri);
    if ((int) theta0.size() != m.d.total)
        Rcpp::stop("theta0 length %d != expected %d", (int) theta0.size(), m.d.total);

    tulpa::LikelihoodSpec lspec;
    lspec.name = "ms_occu_cover";
    lspec.n_processes = 1;
    lspec.gradient_fn = &tulpaObs::ms_occu_cover_nuts_full_grad;

    tulpa::ModelData data;
    data.N = m.d.n_sites;
    data.n_processes = 1;
    data.sigma_beta = sigma_beta;
    data.model_response_data = &m;
    data.likelihood_spec = &lspec;
    data.sharing.init(1);
    data.zi_type = tulpa::ZIType::NONE;
    data.p_zi = 0; data.p_oi = 0;

    tulpa::ParamLayout layout;
    layout.total_params = m.d.total;

    tulpa::set_gradient_mode_str("H");

    std::vector<double> init(theta0.begin(), theta0.end());
    std::vector<double> imv;
    const double* im = nullptr;
    if (inv_metric.isNotNull()) {
        Rcpp::NumericVector v(inv_metric);
        imv.assign(v.begin(), v.end()); im = imv.data();
    }

    tulpa::NUTSFn run_nuts = tulpa::get_nuts_fn();
    tulpa::NUTSResult result = {};
    run_nuts(&data, &layout, init.data(), m.d.total, n_iter, n_warmup,
             max_treedepth, adapt_delta, static_cast<unsigned int>(seed),
             verbose ? 1 : 0, im, &result);

    const int n_samples = result.n_sample, np = m.d.total;
    Rcpp::NumericMatrix draws(n_samples, np);
    Rcpp::NumericVector lp(n_samples), ap(n_samples);
    Rcpp::IntegerVector div(n_samples), td(n_samples);
    for (int s = 0; s < n_samples; ++s) {
        for (int j = 0; j < np; ++j) draws(s, j) = result.samples[s * np + j];
        lp[s] = result.log_prob[s]; ap[s] = result.accept_prob[s];
        div[s] = result.divergent[s]; td[s] = result.treedepth[s];
    }
    const double epsilon = result.epsilon;
    result.free_buffers();
    return Rcpp::List::create(
        Rcpp::Named("draws") = draws, Rcpp::Named("log_prob") = lp,
        Rcpp::Named("accept_prob") = ap, Rcpp::Named("divergent") = div,
        Rcpp::Named("treedepth") = td, Rcpp::Named("epsilon") = epsilon,
        Rcpp::Named("n_params") = np);
}
