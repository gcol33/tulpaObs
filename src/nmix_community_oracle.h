// nmix_community_oracle.h
// Native compiled REGroupOracle for the community / multispecies N-mixture
// (spAbundance msNMix): per-species Royle (2004) N-mixture with Gaussian
// community hyperpriors on the per-species abundance / detection coefficients.
// The grouping factor is the species; the per-group RE vector is
// b_s = (b_lambda_s, b_p_s) (dimension d = p_lambda + p_p), entering the
// abundance and detection linear predictors as coef = mu + b_s where mu is the
// fixed community mean (the engine's theta).
//
// This is the compiled counterpart of the R make_group oracle .nmix_re_oracle()
// in R/nmix_laplace_re.R: it computes, per species, the marginal value / score
// / observed-information by summing the per-site N-mixture marginal
// (nmix_kernel.h) over the species' sites and sandwiching the per-site eta-space
// blocks with the abundance / detection designs. Routing tulpa_nmix_laplace_re
// through this native oracle removes the per-group / per-node round trip into R
// the RClosureOracle bridge incurs, so the integration math stays the single
// shared engine (aghq_re_core) but the per-species likelihood is compiled.
//
// The per-site marginal factorizes over sites given the coefficients, so the
// per-site accumulation reproduces the whole-species marginal exactly; the
// curvature is the design-sandwiched per-site observed-information block
//   B_i = diag(I^lambda_i, I^p_ij) - Var(N_i|y_i) v_i v_i',
//   v_i = (-score_wt_lambda_i, p_i1, ..., p_iJ)
// (Louis 1982; the abundance/detection coupling). The complete-data Fisher (the
// PSD Newton curvature for the mode-find, where the observed info is indefinite
// away from the mode) is the same sandwich with the diagonal Fisher only.
//
// The substantive method bodies live in nmix_community_oracle.cpp (declared
// here): they instantiate a large amount of Eigen template code, and keeping
// them out-of-line means a TU that only *calls* the oracle (e.g. the EM driver
// nmix_community_em.cpp) does not re-instantiate the per-site assembly inline,
// which overflows the MinGW g++ compiler under -O2.
//
// Poisson OR negative-binomial abundance. Under NB the dispersion is a
// PER-SPECIES random effect log_r_s ~ N(mu_log_r, sigma_log_r), so the
// per-species RE vector widens to b_s = (b_lambda_s, b_p_s, b_logr_s) of
// dimension d = p_lambda + p_p + 1 and the community mean log_r enters theta as
// the trailing fixed effect mu_log_r (n_theta = d). The per-species size is
// r_s = exp(mu_log_r + b_logr_s); the engine integrates b_logr_s as a third
// (scalar, diagonal) covariance block. Poisson keeps d = p_lambda + p_p (no
// dispersion coordinate) and is the r = +Inf limit. The per-site NB marginal /
// score / dispersion machinery is nmix_kernel.h (the single source); the log_r
// coordinate's design is the identity (it enters every site directly), so the
// log_r row of the score / observed-info block is assembled from the kernel's
// dispersion outputs (grad_theta / info_theta / info_lambda_theta /
// cov_N_stheta / var_stheta).
//
// Zero-inflation (ZIP / ZINB, spAbundance-style structural absence) is a
// PER-SPECIES structural-zero random effect logit_omega_s ~ N(mu_omega,
// sigma_omega), mirroring the log_r_s design: the per-species RE vector gains a
// trailing logit_omega_s coordinate (identity design, index idx_omega), the
// community mean mu_omega joins theta, and a scalar covariance block integrates
// b_omega_s. The per-site marginal wraps the plain Royle marginal L_i in a
// structural-zero mixture -- for an all-zero site,
//   m_i = log( omega_s + (1 - omega_s) exp(L_i) ) = LSE(log omega_s,
//                                                       log(1 - omega_s) + L_i),
// a detection site rules out N = 0 so m_i = log(1 - omega_s) + L_i. Writing this
// as LSE(c0, c1) with c0 = log(omega_s) (all-zero only, else -Inf) and c1 =
// log(1 - omega_s) + L_i and pi_i = the posterior structural-zero weight
// exp(c0 - m_i) (== 0 at a detection site), the composition is closed form on top
// of the kernel outputs: the abundance/detection/log_r score scales by (1 - pi_i);
//   d m_i / d logit_omega = pi_i - omega_s;
// the marginal observed-info inner block becomes
//   (1 - pi_i) B_i - pi_i (1 - pi_i) g_i g_i'   (g_i the plain inner score),
//   info(omega, omega) = omega_s(1 - omega_s) - pi_i(1 - pi_i),
//   info(omega, inner) = pi_i(1 - pi_i) g_i;
// the complete-data Fisher (the PSD Newton curvature, latent zero-indicator z_i
// observed) is block-diagonal: (1 - pi_i) F_i on the inner block, omega_s(1 -
// omega_s) on the omega diagonal, zero cross. Poisson / NB with no ZI is the
// idx_omega < 0 limit (pi_i = 0, no omega coordinate) -- byte-identical to the
// plain path. ZI has no closed-form EM (like NB), so it is fit by the joint AGHQ
// optimizer only.

#ifndef TULPAOBS_NMIX_COMMUNITY_ORACLE_H
#define TULPAOBS_NMIX_COMMUNITY_ORACLE_H

#include "tulpa/aghq_oracle.h"
#include "nmix_kernel.h"
#include "tobs_math.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <limits>
#include <vector>

namespace tulpaObs {

struct NMixCommunityOracle : tulpa::REGroupOracle {
    int p_lam = 0, p_p = 0;
    int idx_logr = -1;                    // RE/coef index of log_r_s (NB only; -1 = Poisson)
    int idx_omega = -1;                   // RE/coef index of logit_omega_s (ZI only; -1 = no ZI)
    Eigen::MatrixXd Xlam;                 // n_sites x p_lambda (shared across species)
    Eigen::VectorXd mu;                   // active community means (theta), length d
    bool   is_nb = false;                 // negative-binomial abundance (per-species log_r RE)
    bool   is_zi = false;                 // zero-inflated abundance (per-species logit_omega RE)
    int    n_sites = 0;

    // Optional SHARED per-site abundance offset (spAbundance sfMsNMix): the
    // value sigma * f of the outer-integrated areal field, identical across
    // species, added to eta_lambda before the per-site marginal. This is the
    // ONLY thing that distinguishes the spatial community N-mixture from the
    // non-spatial one, so both run on this class. Empty => no field, and
    // eta_lambda starts from 0.0 exactly as on the non-spatial path.
    Eigen::VectorXd offset;

    // 0.0 when no field is attached; the field value at `site` otherwise.
    double site_offset(int site) const {
        return offset.size() ? offset(site) : 0.0;
    }
    // Attach the field values the outer nested-Laplace grid point conditions on.
    void set_offset(const Rcpp::NumericVector& z);

    // Per species, per site: the cached Poisson marginal (lgamma precompute) and
    // the detection design rows for that site's visits, in input order.
    struct SiteRec {
        int site = 0;                     // 0-based row into Xlam
        NMixSiteCache cache;              // eta-independent lgamma terms
        Eigen::MatrixXd Xp;               // J_i x p_p
    };
    std::vector<std::vector<SiteRec>> sp_sites;   // [n_species][n_sites]

    // Per-species value / score, and optionally the marginal observed info
    // (negH) and / or the PSD complete-data Fisher; one site loop is the single
    // source for grad_hess / newton_hess / theta_score.
    struct SpeciesEval {
        double logL = 0.0;
        Eigen::VectorXd grad;             // d ell_g / db  (== d ell_g / d mu); under NB the
                                          // trailing entry is d ell_g / d log_r_s (chain
                                          // rule derivative 1: log_r_s = mu_log_r + b_logr_s)
        Eigen::MatrixXd negH;             // -d^2 ell_g / db^2 (marginal observed info)
        Eigen::MatrixXd fisher;           // complete-data Fisher (PSD)
    };

    NMixCommunityOracle(const Rcpp::IntegerVector& y,
                        const Rcpp::IntegerVector& site_idx,
                        const Rcpp::IntegerVector& species_idx,
                        const Rcpp::NumericMatrix& X_lambda,
                        const Rcpp::NumericMatrix& X_p,
                        int n_sites, int n_species, int K_max,
                        bool nb = false, bool zi = false);

    SpeciesEval eval_species(int g, const double* b,
                             bool want_negH = true, bool want_fisher = true) const;

    void rebind(const double* theta) override {
        // n_theta == d: the community means (mu_lambda, mu_p) and, under NB, the
        // trailing community log-dispersion mu_log_r. The per-species size is
        // r_s = exp(mu_log_r + b_logr_s), assembled per group in eval_species.
        for (int i = 0; i < d; ++i) mu(i) = theta[i];
    }
    void grad_hess(int g, const double* b, double& logL,
                   double* grad, double* negH) const override;
    void node_ll(int g, const double* B, int n_nodes, double* out) const override;
    void theta_score(int g, const double* b, double* dl_dtheta) const override;
    bool newton_hess(int g, const double* b, double* H) const override;
    bool thread_safe() const override { return true; }
};

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_COMMUNITY_ORACLE_H
