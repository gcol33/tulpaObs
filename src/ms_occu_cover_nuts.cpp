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
#include "tobs_shape.h"
#include "occu_coupling_shared.h"
#include "community_chol.h"

#include "nuts_engine.h"
using namespace Rcpp;

namespace tulpaObs {

// The three community-covariance hyperprior scalars plus the two this family
// adds for the shared cover log-dispersion.
struct MsOccuCoverPri : CommunityCholPri {
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

    // Dispersion-RE variant (#115 B7 follow-up): a fourth 1-D community arm.
    // `Pz` is the per-species z-block stride (P_coef + 1 for the RE, else P_tot);
    // `mu_ld_off` the community-mean log-dispersion coord; `chol_ld_off` the
    // 1x1 log-Cholesky (= log sigma_ld).
    bool re_disp = false;
    int Pz = 0, mu_ld_off = 0, chol_ld_off = 0;
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
    if (d.re_disp) {
        // theta = (mu_coef[P_tot], mu_ld[1], {z_s}[S*(P_tot+1)], chol_occ, chol_p,
        //          chol_pos, chol_ld[1]); z_s = (z_occ, z_p, z_pos, z_ld).
        d.Pz = d.P_tot + 1;
        d.mu_ld_off = d.P_tot;
        d.b_off = d.P_tot + 1;
        int coff = (d.P_tot + 1) + d.n_species * d.Pz;
        d.chol_occ_off = coff; coff += d.q_occ;
        d.chol_p_off   = coff; coff += d.q_p;
        d.chol_pos_off = coff; coff += d.q_pos;
        d.chol_ld_off  = coff; coff += 1;
        d.total = coff;
    } else {
        d.Pz = d.P_tot;
        d.b_off = d.P_tot;
        int coff = d.P_tot + d.n_species * d.P_tot;
        d.chol_occ_off = coff; coff += d.q_occ;
        d.chol_p_off   = coff; coff += d.q_p;
        d.chol_pos_off = coff; coff += d.q_pos;
        d.log_disp_off = coff; coff += 1;
        d.total = coff;
    }
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

    namespace sh = tulpaObs::shape;
    sh::check_dim_arg(d.n_sites, "n_sites");
    sh::check_dim_arg(d.max_visits, "max_visits");
    sh::check_dim_arg(d.n_species, "n_species");
    const R_xlen_t n_rows = (R_xlen_t) d.n_sites * d.max_visits;
    sh::check_nrow(d.X_occ, d.n_sites, "X_occ");
    sh::check_nrow(d.X_det_site, d.n_sites, "X_det_site");
    sh::check_nrow(d.X_pos_site, d.n_sites, "X_pos_site");
    if (d.p_det_visit > 0) sh::check_nrow(d.X_det_visit, n_rows, "X_det_visit");
    if (d.p_pos_visit > 0) sh::check_nrow(d.X_pos_visit, n_rows, "X_pos_visit");

    Rcpp::List Y  = Rcpp::as<Rcpp::List>(spec["y"]);      // list of n_species
    Rcpp::List YP = Rcpp::as<Rcpp::List>(spec["y_pos"]);
    Rcpp::List V  = Rcpp::as<Rcpp::List>(spec["valid"]);
    sh::check_len(Y, d.n_species, "y");
    sh::check_len(YP, d.n_species, "y_pos");
    sh::check_len(V, d.n_species, "valid");
    d.sp.reserve(d.n_species);
    for (int s = 0; s < d.n_species; ++s) {
        MsocSpec ms;
        ms.y     = Rcpp::as<IntegerMatrix>(Y[s]);
        ms.y_pos = Rcpp::as<NumericMatrix>(YP[s]);
        ms.valid = Rcpp::as<IntegerMatrix>(V[s]);
        sh::check_dim(ms.y, d.n_sites, d.max_visits, "y[[s]]");
        sh::check_dim(ms.y_pos, d.n_sites, d.max_visits, "y_pos[[s]]");
        sh::check_dim(ms.valid, d.n_sites, d.max_visits, "valid[[s]]");
        d.sp.push_back(ms);
    }
    if (spec.containsElementNamed("re_disp"))
        d.re_disp = Rcpp::as<bool>(spec["re_disp"]);
    ms_occu_cover_nuts_layout(d);
    return d;
}

// One species' data contribution: given the reconstructed per-arm coefficients
// b_occ / b_p / b_pos and the (per-species) dispersion `disp`, sweep every cell's
// two-state occu_cover marginal + the cover arm, writing the design-sandwiched
// per-arm b-space gradients into gbo / gbp / gbpos (zeroed here) and the species
// dispersion score into g_ld_s. Returns the species' data log-likelihood. The
// per-cell math is shared by the shared-dispersion and dispersion-RE targets;
// they differ only in how `disp` is formed and how g_ld_s is used. The eta_*
// buffers are caller-owned scratch (length max_visits).
inline double msoc_cell_sweep(const MsOccuCoverNutsData& d, const MsocSpec& ms,
                              const double* b_occ, const double* b_p,
                              const double* b_pos, double disp,
                              double* gbo, double* gbp, double* gbpos,
                              double& g_ld_s,
                              std::vector<double>& eta_p,
                              std::vector<double>& g_eta_p,
                              std::vector<double>& g_eta_pos,
                              std::vector<double>& eta_pc,
                              std::vector<double>& g_pc) {
    const int N = d.n_sites, J = d.max_visits;
    const int P_occ = d.P_occ, P_p = d.P_p, P_pos = d.P_pos;
    const int pds = d.p_det_site, pdv = d.p_det_visit;
    const int pps = d.p_pos_site, ppv = d.p_pos_visit;
    const double* bp_site   = b_p;
    const double* bp_visit  = b_p + pds;
    const double* bpos_site = b_pos;
    const double* bpos_visit= b_pos + pps;
    for (int k = 0; k < P_occ; ++k) gbo[k]   = 0.0;
    for (int k = 0; k < P_p;   ++k) gbp[k]   = 0.0;
    for (int k = 0; k < P_pos; ++k) gbpos[k] = 0.0;
    double lp = 0.0;
    g_ld_s = 0.0;
    const PosCodeAccess pos(d.pos_code);

    for (int i = 0; i < N; ++i) {
        double eta_psi = 0.0;
        for (int k = 0; k < P_occ; ++k) eta_psi += d.X_occ(i, k) * b_occ[k];
        const double psi = sigmoid_(eta_psi);

        double eta_p_site = 0.0;
        for (int k = 0; k < pds; ++k) eta_p_site += d.X_det_site(i, k) * bp_site[k];
        auto eta_p_of = [&](int v) {
            double e = eta_p_site;
            if (pdv > 0) {
                const int row = i * J + v;
                for (int k = 0; k < pdv; ++k) e += d.X_det_visit(row, k) * bp_visit[k];
            }
            return e;
        };
        auto eta_pos_of = [&](int v) {
            double e = 0.0;
            for (int k = 0; k < pps; ++k) e += d.X_pos_site(i, k) * bpos_site[k];
            if (ppv > 0) {
                const int row = i * J + v;
                for (int k = 0; k < ppv; ++k) e += d.X_pos_visit(row, k) * bpos_visit[k];
            }
            return e;
        };

        double g_eta_psi = 0.0;
        occu_cover_cell_sweep(
            pos, psi, J,
            [&](int v) { return ms.valid(i, v) != 0; },
            [&](int v) { return ms.y(i, v) == 1; },
            [&](int v) { return ms.y_pos(i, v); },
            eta_p_of, eta_pos_of, disp,
            eta_p.data(), eta_pc.data(), g_pc.data(),
            lp, g_eta_psi, g_eta_p.data(), g_eta_pos.data(), g_ld_s);

        occu_cover_cell_sandwich(
            i, J, g_eta_psi, g_eta_p.data(), g_eta_pos.data(),
            d.X_occ, d.X_det_site, d.X_det_visit, d.X_pos_site, d.X_pos_visit,
            gbo, gbp, gbp + pds, gbpos, gbpos + pps);
    }
    return lp;
}

// Joint log-posterior + gradient over theta = (mu, {z_s}, chol_occ, chol_p,
// chol_pos, log_disp). NUTS maximises; returns the log-posterior (no negation)
// and writes the gradient into g. Mirrors .tobs_ms_occu_cover_nuts_logpost.
inline double ms_occu_cover_nuts_eval(const MsOccuCoverNutsData& d, const double* th,
                                      double sigma_beta, const MsOccuCoverPri& pr,
                                      double* g) {
    const int J = d.max_visits, S = d.n_species, P = d.P_tot;
    const int P_occ = d.P_occ, P_p = d.P_p, P_pos = d.P_pos;
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

        std::vector<double> b_occ(P_occ), b_p(P_p), b_pos(P_pos);
        chol_noncentered_recon(C_occ.data(), z_occ, P_occ, b_occ.data(), mu + d.off_occ);
        chol_noncentered_recon(C_p.data(), z_p, P_p, b_p.data(), mu + d.off_p);
        chol_noncentered_recon(C_pos.data(), z_pos, P_pos, b_pos.data(), mu + d.off_pos);

        std::vector<double> gbo(P_occ), gbp(P_p), gbpos(P_pos);
        double g_ld_s = 0.0;
        lp += msoc_cell_sweep(d, d.sp[s], b_occ.data(), b_p.data(), b_pos.data(),
                              disp, gbo.data(), gbp.data(), gbpos.data(), g_ld_s,
                              eta_p, g_eta_p, g_eta_pos, eta_pc, g_pc);
        g_logdisp += g_ld_s;

        double* gz_s = g_z + s * P;
        chol_noncentered_push(C_occ.data(), gbo.data(), z_occ, P_occ,
                              gz_s + d.off_occ, A_occ.data(), g_mu + d.off_occ);
        chol_noncentered_push(C_p.data(), gbp.data(), z_p, P_p,
                              gz_s + d.off_p, A_p.data(), g_mu + d.off_p);
        chol_noncentered_push(C_pos.data(), gbpos.data(), z_pos, P_pos,
                              gz_s + d.off_pos, A_pos.data(), g_mu + d.off_pos);
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

// Dispersion-RE variant: the shared log-dispersion becomes a fourth 1-D community
// arm log_disp_s = mu_ld + sigma_ld * z_ld_s. Reuses msoc_cell_sweep (per-species
// disp), chol_noncentered_recon / chol_noncentered_push; adds the 1-D ld arm
// push-back + priors.
// Mirrors .tobs_ms_occu_cover_re_disp_logpost.
inline double ms_occu_cover_re_disp_nuts_eval(const MsOccuCoverNutsData& d,
                                              const double* th, double sigma_beta,
                                              const MsOccuCoverPri& pr, double* g) {
    const int J = d.max_visits, S = d.n_species;
    const int P_occ = d.P_occ, P_p = d.P_p, P_pos = d.P_pos, P_coef = d.P_tot;
    const int Pz = d.Pz;
    for (int j = 0; j < d.total; ++j) g[j] = 0.0;

    const double* mu = th;                    // mu_coef [P_coef]
    const double  mu_ld = th[d.mu_ld_off];
    const double* z  = th + d.b_off;          // species-major z, stride Pz
    double* g_mu = g;
    double* g_z  = g + d.b_off;
    const double sigma_ld = std::exp(th[d.chol_ld_off]);

    std::vector<double> C_occ, C_p, C_pos;
    chol_unpack_cpp(th + d.chol_occ_off, P_occ, C_occ);
    chol_unpack_cpp(th + d.chol_p_off,   P_p,   C_p);
    chol_unpack_cpp(th + d.chol_pos_off, P_pos, C_pos);

    std::vector<double> A_occ((std::size_t) P_occ * P_occ, 0.0),
                        A_p((std::size_t) P_p * P_p, 0.0),
                        A_pos((std::size_t) P_pos * P_pos, 0.0);
    double A_ld = 0.0, g_mld = 0.0, lp = 0.0;
    std::vector<double> eta_p(J), g_eta_p(J), g_eta_pos(J), eta_pc(J), g_pc(J);

    for (int s = 0; s < S; ++s) {
        const double* z_s   = z + s * Pz;
        const double* z_occ = z_s + d.off_occ;
        const double* z_p   = z_s + d.off_p;
        const double* z_pos = z_s + d.off_pos;
        const double  z_ld  = z_s[P_coef];         // z_ld at the end of the block
        std::vector<double> b_occ(P_occ), b_p(P_p), b_pos(P_pos);
        chol_noncentered_recon(C_occ.data(), z_occ, P_occ, b_occ.data(), mu + d.off_occ);
        chol_noncentered_recon(C_p.data(), z_p, P_p, b_p.data(), mu + d.off_p);
        chol_noncentered_recon(C_pos.data(), z_pos, P_pos, b_pos.data(), mu + d.off_pos);
        const double disp_s = std::exp(mu_ld + sigma_ld * z_ld);

        std::vector<double> gbo(P_occ), gbp(P_p), gbpos(P_pos);
        double g_ld_s = 0.0;
        lp += msoc_cell_sweep(d, d.sp[s], b_occ.data(), b_p.data(), b_pos.data(),
                              disp_s, gbo.data(), gbp.data(), gbpos.data(), g_ld_s,
                              eta_p, g_eta_p, g_eta_pos, eta_pc, g_pc);

        double* gz_s = g_z + s * Pz;
        chol_noncentered_push(C_occ.data(), gbo.data(), z_occ, P_occ,
                              gz_s + d.off_occ, A_occ.data(), g_mu + d.off_occ);
        chol_noncentered_push(C_p.data(), gbp.data(), z_p, P_p,
                              gz_s + d.off_p, A_p.data(), g_mu + d.off_p);
        chol_noncentered_push(C_pos.data(), gbpos.data(), z_pos, P_pos,
                              gz_s + d.off_pos, A_pos.data(), g_mu + d.off_pos);
        // 1-D dispersion arm: log_disp_s = mu_ld + sigma_ld * z_ld.
        g_mld        += g_ld_s;
        gz_s[P_coef] += sigma_ld * g_ld_s;
        A_ld         += g_ld_s * z_ld;
    }

    for (int j = 0; j < S * Pz; ++j) {
        const double zz = z[j]; g_z[j] -= zz; lp += -0.5 * zz * zz;
    }

    lp += chol_data_grad_noncentered(A_occ, C_occ, P_occ, th + d.chol_occ_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_occ_off);
    lp += chol_data_grad_noncentered(A_p, C_p, P_p, th + d.chol_p_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_p_off);
    lp += chol_data_grad_noncentered(A_pos, C_pos, P_pos, th + d.chol_pos_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_pos_off);
    std::vector<double> A_ld_v(1, A_ld), C_ld_v(1, sigma_ld);
    lp += chol_data_grad_noncentered(A_ld_v, C_ld_v, 1, th + d.chol_ld_off,
                                     pr.logdiag_mean, pr.logdiag_sd, pr.offdiag_sd,
                                     g + d.chol_ld_off);

    const double inv_sb2 = 1.0 / (sigma_beta * sigma_beta);
    for (int k = 0; k < P_coef; ++k) {
        g_mu[k] -= inv_sb2 * mu[k];
        lp      += -0.5 * inv_sb2 * mu[k] * mu[k];
    }
    const double ld = mu_ld - pr.log_disp_mean;
    const double ild2 = 1.0 / (pr.log_disp_sd * pr.log_disp_sd);
    lp += -0.5 * ild2 * ld * ld;
    g[d.mu_ld_off] = g_mld - ild2 * ld;

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
    const double lp = m->d.re_disp
        ? ms_occu_cover_re_disp_nuts_eval(m->d, params.data(), m->sigma_beta,
                                          m->pr, grad.data())
        : ms_occu_cover_nuts_eval(m->d, params.data(), m->sigma_beta,
                                  m->pr, grad.data());
    if (log_post_out) *log_post_out = lp;
}

inline MsOccuCoverPri ms_occu_cover_pri_from_list(const Rcpp::List& pri) {
    MsOccuCoverPri pr;
    community_chol_pri_read(pri, pr);
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
    const double lp = d.re_disp
        ? tulpaObs::ms_occu_cover_re_disp_nuts_eval(d, theta.begin(), sigma_beta, pr, grad.begin())
        : tulpaObs::ms_occu_cover_nuts_eval(d, theta.begin(), sigma_beta, pr, grad.begin());
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

    return tulpaObs::run_tulpa_nuts(
        &tulpaObs::ms_occu_cover_nuts_full_grad, &m, m.d.total,
        theta0, sigma_beta, tulpaObs::shape::optional_numeric(inv_metric.get(), "inv_metric"),
        n_iter, n_warmup, max_treedepth, adapt_delta, seed, verbose,
        "ms_occu_cover", m.d.n_sites);
}
