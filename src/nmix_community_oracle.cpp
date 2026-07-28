// nmix_community_oracle.cpp
// Out-of-line bodies for NMixCommunityOracle (declared in nmix_community_oracle.h).
// Kept out of the header so a translation unit that only calls the oracle (the EM
// driver nmix_community_em.cpp, and the engine shim aghq_re.cpp) does not
// re-instantiate the per-site Eigen assembly inline -- which overflows MinGW g++
// under -O2. The per-site marginal math is nmix_kernel.h (the single source).

#include "nmix_community_oracle.h"
#include "nmix_oracle_emit.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <cstddef>
#include <vector>

namespace tulpaObs {

NMixCommunityOracle::NMixCommunityOracle(const Rcpp::IntegerVector& y,
                                         const Rcpp::IntegerVector& site_idx,
                                         const Rcpp::IntegerVector& species_idx,
                                         const Rcpp::NumericMatrix& X_lambda,
                                         const Rcpp::NumericMatrix& X_p,
                                         int n_sites, int n_species, int K_max,
                                         bool nb, bool zi, int headroom) {
    this->n_sites = n_sites;
    p_lam    = X_lambda.ncol();
    p_p      = X_p.ncol();
    is_nb    = nb;
    is_zi    = zi;
    // Under NB the per-species RE vector carries a trailing log_r_s coordinate
    // (dispersion is a per-species random effect); Poisson omits it. Under ZI a
    // further trailing logit_omega_s coordinate carries the per-species
    // structural-zero probability. Coordinate order: [lambda | p | log_r? |
    // omega?]; n_theta == d in every case (the community means, plus mu_log_r
    // under NB and mu_omega under ZI).
    idx_logr  = nb ? (p_lam + p_p) : -1;
    idx_omega = zi ? (p_lam + p_p + (nb ? 1 : 0)) : -1;
    d        = p_lam + p_p + (nb ? 1 : 0) + (zi ? 1 : 0);
    n_theta  = d;
    n_groups = n_species;
    mu       = Eigen::VectorXd::Zero(d);

    Xlam.resize(n_sites, p_lam);
    for (int i = 0; i < n_sites; ++i)
        for (int c = 0; c < p_lam; ++c) Xlam(i, c) = X_lambda(i, c);

    // Group the long-form rows by (species, site), preserving input order.
    const int n_obs = y.size();
    std::vector<std::vector<std::vector<int>>> rows(
        n_species, std::vector<std::vector<int>>(n_sites));
    for (int r = 0; r < n_obs; ++r)
        rows[species_idx[r] - 1][site_idx[r] - 1].push_back(r);

    sp_sites.assign(n_species, std::vector<SiteRec>());
    for (int s = 0; s < n_species; ++s) {
        sp_sites[s].reserve(n_sites);
        for (int i = 0; i < n_sites; ++i) {
            SiteRec rec;
            rec.site = i;
            const std::vector<int>& rr = rows[s][i];
            const int J = (int)rr.size();
            std::vector<int> yv(J);
            rec.Xp.resize(J, p_p);
            for (int j = 0; j < J; ++j) {
                yv[j] = y[rr[j]];
                for (int c = 0; c < p_p; ++c) rec.Xp(j, c) = X_p(rr[j], c);
            }
            rec.cache = nmix_precompute_site(yv.data(), J, K_max, headroom);
            sp_sites[s].push_back(std::move(rec));
        }
    }
}

NMixCommunityOracle::SpeciesEval
NMixCommunityOracle::eval_species(int g, const double* b,
                                  bool want_negH, bool want_fisher) const {
    SpeciesEval e;
    e.grad = Eigen::VectorXd::Zero(d);
    if (want_negH)   e.negH   = Eigen::MatrixXd::Zero(d, d);
    if (want_fisher) e.fisher = Eigen::MatrixXd::Zero(d, d);

    Eigen::VectorXd coef(d);
    for (int i = 0; i < d; ++i) coef(i) = mu(i) + b[i];

    // Per-species NB size r_s = exp(mu_log_r + b_logr_s) = exp(coef[idx_logr]);
    // Poisson is the r = +Inf limit (idx_logr < 0). The log_r coordinate enters
    // every site directly (design = identity), so its row of the score / curvature
    // is the per-site dispersion outputs summed over the species' sites, with the
    // lambda<->log_r cross sandwiched only through the abundance design.
    const double r_s = is_nb ? std::exp(coef(idx_logr))
                             : std::numeric_limits<double>::infinity();

    // Per-species structural-zero probability omega_s = plogis(logit_omega_s)
    // (ZI only); the omega coordinate enters every all-zero site's marginal
    // through the mixture wrap, design = identity.
    double om = 0.0, log_om = 0.0, log1m_om = 0.0;
    if (is_zi) {
        logit_log_probs(coef(idx_omega), log_om, log1m_om);
        om = std::exp(log_om);
    }

    std::vector<double> eta_p;
    for (const SiteRec& rec : sp_sites[g]) {
        const int J = rec.cache.n_visits;
        double eta_lam = site_offset(rec.site);
        for (int c = 0; c < p_lam; ++c) eta_lam += Xlam(rec.site, c) * coef(c);
        eta_lam = tulpaObs::clamp_eta(eta_lam);

        eta_p.assign(J, 0.0);
        for (int j = 0; j < J; ++j) {
            double v = 0.0;
            for (int c = 0; c < p_p; ++c) v += rec.Xp(j, c) * coef(p_lam + c);
            eta_p[j] = tulpaObs::clamp_eta(v);
        }

        const NMixSiteResult res =
            compute_nmix_site_cached(rec.cache, eta_p.data(), eta_lam, r_s);

        // Structural-zero mixture wrap (ZI). pi = posterior structural-zero weight
        // exp(c0 - m), c0 = log(omega) (all-zero site only), c1 = log(1-omega) +
        // L_i; pi == 0 at any detection site (a detection rules out N = 0). w =
        // 1 - pi scales the abundance / detection / dispersion score. Without ZI,
        // pi = 0 / w = 1 and this is byte-identical to the plain Royle path.
        double pi = 0.0, w = 1.0;
        if (is_zi) {
            const double c1 = log1m_om + res.log_lik;
            if (rec.cache.K_lo == 0) {                 // max(y_i) == 0: all-zero site
                const double c0 = log_om;
                const double mx = c0 > c1 ? c0 : c1;
                const double m  = mx + std::log(std::exp(c0 - mx) + std::exp(c1 - mx));
                pi = std::exp(c0 - m);
                w  = 1.0 - pi;
                e.logL += m;
            } else {
                e.logL += c1;                          // pi = 0
            }
        } else {
            e.logL += res.log_lik;
        }

        // Score: inner (lambda / p / log_r) scaled by w; the omega coordinate is
        // d m_i / d logit_omega = pi - omega_s (identity design).
        for (int c = 0; c < p_lam; ++c)
            e.grad(c) += w * Xlam(rec.site, c) * res.grad_eta_lambda;
        for (int j = 0; j < J; ++j)
            for (int c = 0; c < p_p; ++c)
                e.grad(p_lam + c) += w * rec.Xp(j, c) * res.grad_eta_p[j];
        if (is_nb) e.grad(idx_logr)  += w * res.grad_theta;
        if (is_zi) e.grad(idx_omega) += (pi - om);

        if (!want_negH && !want_fisher) continue;

        // Per-site inner block. Coords: 0 = lambda, 1..J = visits, (NB) log_r at
        // tt = 1 + J. The ZI wrap extends by an omega coord at oo = dd (du coords).
        const int dd = 1 + J + (is_nb ? 1 : 0);
        const int tt = 1 + J;                       // log_r inner-coord index (NB only)
        const int du = dd + (is_zi ? 1 : 0);
        const int oo = dd;                          // omega inner-coord index (ZI only)

        // Plain marginal inner block Bm (observed info) and complete-data Fisher
        // Bf, exactly as the non-ZI path; the ZI wrap rescales / extends them.
        Eigen::MatrixXd Bm, Bf;
        if (want_negH)   { Bm = Eigen::MatrixXd::Zero(dd, dd); Bm(0, 0) = res.info_eta_lambda; }
        if (want_fisher) { Bf = Eigen::MatrixXd::Zero(dd, dd); Bf(0, 0) = res.info_eta_lambda; }
        for (int j = 0; j < J; ++j) {
            if (want_negH)   Bm(1 + j, 1 + j) = res.info_eta_p[j];
            if (want_fisher) Bf(1 + j, 1 + j) = res.info_eta_p[j];
        }
        if (is_nb) {
            // Complete-data Fisher: log_r diagonal info_theta and the lambda<->log_r
            // cross info_lambda_theta (Poisson-neutral / zero under r = +Inf).
            if (want_fisher) {
                Bf(tt, tt) = res.info_theta;
                Bf(0, tt)  = Bf(tt, 0) = res.info_lambda_theta;
            }
            if (want_negH) {
                Bm(tt, tt) = res.info_theta;
                Bm(0, tt)  = Bm(tt, 0) = res.info_lambda_theta;
            }
        }
        if (want_negH && (J > 0 || is_nb)) {
            // Marginal observed info = complete-data Fisher - Cov(score | y). The eta
            // scores are affine in N, contributing the rank-1 var_N * vv vv'
            // (vv = N-coefficients, sign squares out); the log_r score is non-affine
            // (psi(N+r)), so its row uses the kernel's cov_N_stheta / var_stheta.
            Eigen::VectorXd vv(dd);
            vv(0) = -res.score_wt_lambda;
            for (int j = 0; j < J; ++j) {
                const double pj = (eta_p[j] > 0.0)
                    ? 1.0 / (1.0 + std::exp(-eta_p[j]))
                    : std::exp(eta_p[j]) / (1.0 + std::exp(eta_p[j]));
                vv(1 + j) = pj;
            }
            if (is_nb) {
                vv(tt) = 0.0;
                Bm.noalias() -= res.var_N * (vv * vv.transpose());
                for (int c = 0; c < tt; ++c) {
                    const double cv = -vv(c) * res.cov_N_stheta;
                    Bm(c, tt) -= cv;
                    Bm(tt, c) -= cv;
                }
                Bm(tt, tt) -= res.var_stheta;
            } else if (J > 0) {
                Bm.noalias() -= res.var_N * (vv * vv.transpose());
            }
        }

        // Inner plain score g_i, needed for the ZI rank-1 / omega-cross terms.
        Eigen::VectorXd ginner;
        if (is_zi) {
            ginner = Eigen::VectorXd::Zero(dd);
            ginner(0) = res.grad_eta_lambda;
            for (int j = 0; j < J; ++j) ginner(1 + j) = res.grad_eta_p[j];
            if (is_nb) ginner(tt) = res.grad_theta;
        }

        // Assemble the du x du observed-info / Fisher blocks. Marginal observed
        // info: inner -> (1-pi) Bm - pi(1-pi) g g', omega/omega -> om(1-om) -
        // pi(1-pi), omega/inner -> pi(1-pi) g. Complete-data Fisher (PSD Newton
        // curvature, latent z_i observed): block-diagonal, inner -> (1-pi) Bf,
        // omega -> om(1-om). Without ZI, du == dd, w == 1 and these are the plain
        // blocks unchanged.
        Eigen::MatrixXd Bobs, Bfis;
        if (want_negH) {
            Bobs = Eigen::MatrixXd::Zero(du, du);
            Bobs.topLeftCorner(dd, dd) = w * Bm;
            if (is_zi) {
                Bobs.topLeftCorner(dd, dd).noalias() -=
                    (pi * w) * (ginner * ginner.transpose());
                Bobs(oo, oo) = om * (1.0 - om) - pi * w;
                for (int c = 0; c < dd; ++c) {
                    const double cv = (pi * w) * ginner(c);
                    Bobs(oo, c) = cv; Bobs(c, oo) = cv;
                }
            }
        }
        if (want_fisher) {
            Bfis = Eigen::MatrixXd::Zero(du, du);
            Bfis.topLeftCorner(dd, dd) = w * Bf;
            if (is_zi) Bfis(oo, oo) = om * (1.0 - om);
        }

        // Design map Z_i (du x d): eta coord 0 -> Xlam_i over the lambda coefs,
        // eta coord 1+j -> Xp row j over the p coefs, (NB) log_r -> identity, and
        // (ZI) the omega coord -> identity on the omega RE coordinate.
        Eigen::MatrixXd Zi = Eigen::MatrixXd::Zero(du, d);
        for (int c = 0; c < p_lam; ++c) Zi(0, c) = Xlam(rec.site, c);
        for (int j = 0; j < J; ++j)
            for (int c = 0; c < p_p; ++c) Zi(1 + j, p_lam + c) = rec.Xp(j, c);
        if (is_nb) Zi(tt, idx_logr)  = 1.0;
        if (is_zi) Zi(oo, idx_omega) = 1.0;

        const Eigen::MatrixXd Zt = Zi.transpose();
        if (want_negH) {
            const Eigen::MatrixXd BZ = Bobs * Zi;
            e.negH.noalias() += Zt * BZ;
        }
        if (want_fisher) {
            const Eigen::MatrixXd BZ = Bfis * Zi;
            e.fisher.noalias() += Zt * BZ;
        }
    }
    return e;
}

void NMixCommunityOracle::grad_hess(int g, const double* b, double& logL,
                                    double* grad, double* negH) const {
    emit_grad_hess(eval_species(g, b, /*want_negH=*/true, /*want_fisher=*/false),
                   d, logL, grad, negH);
}

void NMixCommunityOracle::node_ll(int g, const double* B, int n_nodes,
                                  double* out) const {
    Eigen::VectorXd coef(d);
    std::vector<double> eta_p;
    std::vector<double> a_scratch;   // per-call log-weight buffer (node_ll only
                                     // needs log L_i, so it uses the fast path
                                     // that skips the moment / dispersion pass);
                                     // local, so thread-safe across group sweeps.
    for (int k = 0; k < n_nodes; ++k) {
        const double* bk = B + (std::size_t)k * d;
        for (int i = 0; i < d; ++i) coef(i) = mu(i) + bk[i];
        const double r_s = is_nb ? std::exp(coef(idx_logr))
                                 : std::numeric_limits<double>::infinity();
        double log_om = 0.0, log1m_om = 0.0;
        if (is_zi) logit_log_probs(coef(idx_omega), log_om, log1m_om);
        double ll = 0.0;
        for (const SiteRec& rec : sp_sites[g]) {
            const int J = rec.cache.n_visits;
            double eta_lam = site_offset(rec.site);
            for (int c = 0; c < p_lam; ++c) eta_lam += Xlam(rec.site, c) * coef(c);
            eta_lam = tulpaObs::clamp_eta(eta_lam);
            eta_p.assign(J, 0.0);
            for (int j = 0; j < J; ++j) {
                double v = 0.0;
                for (int c = 0; c < p_p; ++c) v += rec.Xp(j, c) * coef(p_lam + c);
                eta_p[j] = tulpaObs::clamp_eta(v);
            }
            const double llr =
                nmix_loglik_cached(rec.cache, eta_p.data(), eta_lam, r_s, a_scratch);
            if (is_zi) {
                const double c1 = log1m_om + llr;
                if (rec.cache.K_lo == 0) {             // all-zero site: mix in the 0
                    const double c0 = log_om;
                    const double mx = c0 > c1 ? c0 : c1;
                    ll += mx + std::log(std::exp(c0 - mx) + std::exp(c1 - mx));
                } else {
                    ll += c1;
                }
            } else {
                ll += llr;
            }
        }
        out[k] = ll;
    }
}

void NMixCommunityOracle::theta_score(int g, const double* b,
                                      double* dl_dtheta) const {
    // n_theta == d. The community means' score is the eta-coordinate RE score; the
    // trailing mu_log_r score equals the b_logr_s score (chain rule derivative 1,
    // log_r_s = mu_log_r + b_logr_s), already carried in e.grad[idx_logr].
    const SpeciesEval e = eval_species(g, b, /*want_negH=*/false, /*want_fisher=*/false);
    for (int i = 0; i < d; ++i) dl_dtheta[i] = e.grad(i);
}

bool NMixCommunityOracle::newton_hess(int g, const double* b, double* H) const {
    emit_fisher(eval_species(g, b, /*want_negH=*/false, /*want_fisher=*/true), d, H);
    return true;
}


// Attach the shared per-site abundance offset (= sigma * f). The outer
// nested-Laplace driver conditions on one field value per grid point and calls
// this before each inner AGHQ solve.
void NMixCommunityOracle::set_offset(const Rcpp::NumericVector& z) {
    if ((int) z.size() != n_sites)
        Rcpp::stop("set_offset: length(z) must equal n_sites.");
    offset.resize(n_sites);
    for (int i = 0; i < n_sites; ++i) offset(i) = z[i];
}

}  // namespace tulpaObs

// Rcpp factory: build the native community / multispecies N-mixture oracle and
// return it as an XPtr<tulpa::REGroupOracle>. tulpa::tulpa_re_aghq() consumes
// the XPtr through the engine's REGroupOracle interface; the per-species
// marginal / score / observed-info assembly runs entirely in tulpaObs, with no
// per-group / per-node round trip into R.
// [[Rcpp::export]]
SEXP cpp_nmix_community_oracle(Rcpp::IntegerVector y, Rcpp::IntegerVector site_idx,
                               Rcpp::IntegerVector species_idx,
                               Rcpp::NumericMatrix X_lambda,
                               Rcpp::NumericMatrix X_p,
                               int n_sites, int n_species, int K_max,
                               bool nb = false, bool zi = false,
                               int headroom = -1) {
    return Rcpp::XPtr<tulpa::REGroupOracle>(
        new tulpaObs::NMixCommunityOracle(y, site_idx, species_idx,
                                          X_lambda, X_p,
                                          n_sites, n_species, K_max, nb, zi,
                                          headroom),
        true);
}
