// count_grouped_oracle.h
// Shared base for single-species count-model grouped-RE oracles (the native
// REGroupOracle that drives tulpa::tulpa_re_aghq for abun / removal / ... with
// a site-level grouped random effect on EITHER the abundance arm or the
// detection arm).
//
// Every such oracle has the SAME shape: a per-group RE vector b (dim d) enters
// ONE arm through a site-level design row Z_i, so the per-site eta-space block
// collapses to a SCALAR per site (the abundance arm's (1,1) entry, or the
// detection arm's contracted ones'-vector form), the marginal observed info
// carries the rank-1 Var(N|y) coupling (Louis 1982), and each (g, b) call
// iterates only over the group's sites. The ONLY family-specific pieces are the
// per-site marginal kernel and the per-site data it needs. This base owns the
// Z-sandwich, the group loop, rebind, and every REGroupOracle method; a derived
// class supplies the per-site data in its constructor (via build_common) and the
// per-site marginal through eval_site().
//
// The per-site marginal is whatever returns an NMixSiteResult: the N-mixture
// (nmix_kernel.h) and the removal-sampling kernel (removal_kernel.h) both do,
// with identical detection-arm score/info semantics, so they share this base.

#ifndef TULPAOBS_COUNT_GROUPED_ORACLE_H
#define TULPAOBS_COUNT_GROUPED_ORACLE_H

#include "tulpa/aghq_oracle.h"
#include "nmix_kernel.h"      // NMixSiteResult
#include "tobs_math.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <limits>
#include <vector>

namespace tulpaObs {

struct CountGroupedOracle : tulpa::REGroupOracle {
    // Arm selector: 0 = lambda (RE shifts the abundance log-predictor),
    //               1 = p      (RE is a per-site scalar offset applied
    //                           uniformly to every visit / pass at that site).
    int arm = 0;
    int p_lam = 0, p_p = 0;
    bool   is_nb = false;
    double r     = std::numeric_limits<double>::infinity();  // active NB size
    int n_sites = 0;

    // Abundance design (one row per site) and per-site detection design.
    Eigen::MatrixXd Xlam;                          // n_sites x p_lam
    std::vector<Eigen::MatrixXd> Xp_site;          // [n_sites] of J_i x p_p
    Eigen::MatrixXd Zsite;                          // n_sites x d (RE design rows)

    // Group membership: per-group list of site indices (0-based).
    std::vector<std::vector<int>> sites_by_group;

    // Active fixed-effect predictors (rebuilt in rebind), clamped to [-30, 30].
    Eigen::VectorXd eta_lambda_base;               // length n_sites
    std::vector<std::vector<double>> eta_p_site;   // per-site J_i entries

    // The single family-specific hook: the per-site marginal at this site's
    // current eta. `eta_p_ptr` is the per-visit/per-pass detection logit vector
    // (length Xp_site[i].rows(), or null when the site has no detection rows);
    // `eta_lam` the site log-abundance predictor; `r` is the member dispersion.
    virtual NMixSiteResult eval_site(int i, const double* eta_p_ptr,
                                     double eta_lam) const = 0;

    // The per-site predictors at RE value b. The abundance arm shifts the site
    // log-abundance; the detection arm shifts every visit's logit by the same
    // scalar, which needs a shifted copy of the site's detection vector --
    // `scratch` holds it, and the returned `eta_p_ptr` aliases that buffer, so
    // it stays valid until the next call on the same scratch. Callers reuse one
    // scratch across a group's sites to keep the loop allocation-free.
    void site_eta(int i, const double* b, std::vector<double>& scratch,
                  double& eta_lam, const double*& eta_p_ptr) const {
        double shift = 0.0;
        for (int c = 0; c < d; ++c) shift += Zsite(i, c) * b[c];
        if (arm == 0) {
            eta_lam   = clamp_eta(eta_lambda_base(i) + shift);
            eta_p_ptr = eta_p_site[i].empty() ? nullptr : eta_p_site[i].data();
        } else {
            eta_lam = eta_lambda_base(i);
            const int J = (int)Xp_site[i].rows();
            const double shift_cl = clamp_eta(shift);
            scratch.assign(J, 0.0);
            for (int j = 0; j < J; ++j) scratch[j] = eta_p_site[i][j] + shift_cl;
            eta_p_ptr = scratch.empty() ? nullptr : scratch.data();
        }
    }

    // Per-site evaluation request flags (see nmix_re_oracle.h history): the
    // Newton mode-find asks for negH (marginal observed info) and the PSD Fisher
    // at the SAME (g, b), so one site loop fills both on demand.
    struct GroupEval {
        double logL = 0.0;
        Eigen::VectorXd grad;          // d ell_g / db, length d
        Eigen::MatrixXd negH;          // -d^2 ell_g / db^2, marginal observed info
        Eigen::MatrixXd fisher;        // complete-data Fisher (PSD), Newton curvature
        Eigen::VectorXd theta_grad;    // d ell_g / d theta, length n_theta
    };
    GroupEval eval_group(int g, const double* b,
                         bool want_negH, bool want_fisher,
                         bool want_theta_grad) const;

    void rebind(const double* theta) override;
    void grad_hess(int g, const double* b, double& logL,
                   double* grad, double* negH) const override;
    void node_ll(int g, const double* B, int n_nodes, double* out) const override;
    void theta_score(int g, const double* b, double* dl_dtheta) const override;
    bool newton_hess(int g, const double* b, double* H) const override;
    bool thread_safe() const override { return true; }

protected:
    // Fill the shared members from the common inputs and hand back the per-site
    // count vectors (input row order preserved) so the derived constructor can
    // build its own per-site kernel data. Called from derived constructors
    // (a base constructor cannot dispatch to a derived eval_site).
    void build_common(int arm_,
                      const Rcpp::IntegerVector& y,
                      const Rcpp::IntegerVector& site_idx,
                      const Rcpp::NumericMatrix& X_lambda,
                      const Rcpp::NumericMatrix& X_p,
                      const Rcpp::NumericMatrix& Z_site,
                      const Rcpp::IntegerVector& site_group,
                      int n_sites_, int n_groups_, bool nb,
                      std::vector<std::vector<int>>& y_by_site_out);
};

}  // namespace tulpaObs

#endif  // TULPAOBS_COUNT_GROUPED_ORACLE_H
