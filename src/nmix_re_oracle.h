// nmix_re_oracle.h
// Native compiled REGroupOracle for the single-species N-mixture with
// site-level grouped random effects on EITHER the abundance (lambda) or
// detection (p) arm (gcol33/tulpaObs#13). The grouping factor is something
// other than species (station, observer-per-site, site cluster); the per-group
// RE vector b has dimension d = sum over RE-blocks of n_coefs (1..3), entering
// ONE arm via the per-site design row Z_i (length d).
//
// This is the single-species analogue of NMixCommunityOracle: same per-site
// kernel (nmix_kernel.h), same Z-sandwich pattern, but theta is the plain
// fixed-effect vector (no community-mean concept) and the RE enters one arm
// only (so the per-site eta-space block collapses to a SCALAR per site, the
// abundance arm's (1,1) entry or the p arm's contracted ones'-vector form;
// the design map Z_i is rank-1 in eta). Each (g, b) call iterates only over
// the group's sites; the integration math stays in the shared engine
// (aghq_re_core) and the per-group / per-node loop never crosses into R.
//
// Poisson OR negative-binomial abundance. Under NB the global dispersion size
// r is the (d+1)-th theta entry log_r (carried as theta[n_theta - 1] when
// is_nb); Poisson is the r = +Inf limit, threaded through nmix_kernel.h.

#ifndef TULPAOBS_NMIX_RE_ORACLE_H
#define TULPAOBS_NMIX_RE_ORACLE_H

#include "tulpa/aghq_oracle.h"
#include "nmix_kernel.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <limits>
#include <vector>

namespace tulpaObs {

struct NMixGroupedOracle : tulpa::REGroupOracle {
    // Arm selector: 0 = lambda (RE shifts the abundance log-predictor),
    //               1 = p      (RE is a per-site scalar offset applied
    //                           uniformly to every visit at that site).
    int arm = 0;
    int p_lam = 0, p_p = 0;
    bool   is_nb = false;
    double r     = std::numeric_limits<double>::infinity();  // active NB size

    // Designs. Xlam: n_sites x p_lam (one row per site). Xp_site[i]: J_i x p_p
    // for site i, rows in input visit order. y_site[i]: J_i counts.
    Eigen::MatrixXd Xlam;
    std::vector<Eigen::MatrixXd> Xp_site;
    std::vector<NMixSiteCache>    site_cache;  // eta-independent lgamma terms
    int n_sites = 0;

    // Per-site RE design row (length d), site-level. Materialised from cbind of
    // per-term Z's in R (intercept columns expanded to 1s). Row i is Z_i.
    Eigen::MatrixXd Zsite;

    // Group membership: per-group list of site indices (0-based).
    std::vector<std::vector<int>> sites_by_group;

    // Active fixed-effect predictors (rebuilt in rebind). Both are clamped to
    // [-30, 30] in eta-space, matching the R-closure path.
    Eigen::VectorXd eta_lambda_base;                 // length n_sites
    std::vector<std::vector<double>> eta_p_site;     // per-site J_i entries

    NMixGroupedOracle(int arm_,
                      const Rcpp::IntegerVector& y,
                      const Rcpp::IntegerVector& site_idx,
                      const Rcpp::NumericMatrix& X_lambda,
                      const Rcpp::NumericMatrix& X_p,
                      const Rcpp::NumericMatrix& Z_site,
                      const Rcpp::IntegerVector& site_group,
                      int n_sites_, int n_groups_, int K_max,
                      bool nb);

    static inline double clamp30(double e) {
        return e < -30.0 ? -30.0 : (e > 30.0 ? 30.0 : e);
    }

    // Per-site evaluation request flags. The Newton mode-find calls grad_hess
    // (negH = marginal observed info) THEN newton_hess (PSD Fisher) at the
    // SAME (g, b); we expose both via one site loop with selector flags so the
    // hot path runs the kernel once even when the engine asks for both shapes.
    struct GroupEval {
        double logL = 0.0;
        Eigen::VectorXd grad;          // d ell_g / db, length d
        Eigen::MatrixXd negH;          // -d^2 ell_g / db^2, marginal observed info
        Eigen::MatrixXd fisher;        // complete-data Fisher (PSD), Newton curvature
        Eigen::VectorXd theta_grad;    // d ell_g / d theta, length n_theta (data score)
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
};

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_RE_ORACLE_H
