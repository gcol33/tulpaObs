// distance_re_oracle.h
// Native compiled REGroupOracle for the binned distance-sampling abundance
// model with a site-level grouped random effect on the ABUNDANCE (lambda) arm
// (gcol33/tulpaObs#51). Shares the entire Z-sandwich + group-loop +
// REGroupOracle plumbing with the N-mixture / removal grouped oracles via
// CountGroupedOracle; the only distance-specific pieces are the per-site
// binned-multinomial marginal (distance_kernel.h, integrated by the per-fit
// Gauss-Legendre quadrature) and the per-site bin counts.
//
// The detection arm is the per-site log-scale predictor `eta_sigma` -- ONE row
// per site -- so the half-normal key's theta vector is exactly the count-family
// layout [beta_lambda | beta_sigma | log_r?]: the lone detection row at each
// site carries the sigma gradient as grad_eta_p[0], and the abundance-arm RE
// never touches the detection-curvature path. The hazard-rate key adds a global
// scalar shape coordinate that is not a per-site design column, so it is not
// expressible in this base; the R wrapper gates it (half-normal key only). The
// areal-spatial path DOES carry the hazard shape -- it folds the scalar log-shape
// into the areal-BFGS fixed block (tulpaObs#79) rather than this per-site
// grouped-RE theta layout.
//
// Detection-arm RE is likewise out of scope here: the latent N couples a site's
// bins, so a detection RE would not factorize into the per-site scalar offset
// the base assumes. The R wrapper routes only an abundance-arm RE through this
// oracle. Poisson OR negative-binomial abundance (the r = +Inf limit).

#ifndef TULPAOBS_DISTANCE_RE_ORACLE_H
#define TULPAOBS_DISTANCE_RE_ORACLE_H

#include "count_grouped_oracle.h"
#include "distance_kernel.h"
#include "distance_quad.h"
#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

struct DistanceGroupedOracle : CountGroupedOracle {
    std::vector<std::vector<int>> y_bins_site;  // per-site bin counts
    int n_bins = 0;
    int K_max  = 0;
    DistQuad quad;                              // per-fit detection quadrature

    DistanceGroupedOracle(int arm_,
                          const Rcpp::IntegerMatrix& y_bins,
                          const Rcpp::NumericMatrix& X_lambda,
                          const Rcpp::NumericMatrix& X_sigma,
                          const Rcpp::NumericMatrix& Z_site,
                          const Rcpp::IntegerVector& site_group,
                          int n_sites_, int n_groups_,
                          const Rcpp::NumericVector& cutpoints,
                          int transect, int quad_order, int K_max_,
                          bool nb);

    // The RE arm is the abundance arm; eta_p_ptr[0] is the site's log-sigma.
    // Pack the distance marginal's abundance + NB-dispersion fields into the
    // NMixSiteResult the base consumes, with the lone detection row carrying the
    // sigma gradient (grad_eta_p[0]) and its PSD Fisher block (info_eta_p[0]).
    NMixSiteResult eval_site(int i, const double* eta_p_ptr,
                             double eta_lam) const override {
        const double eta_sigma = (eta_p_ptr != nullptr) ? eta_p_ptr[0] : 0.0;
        const DistSiteResult d = compute_distance_site(
            y_bins_site[i].data(), n_bins, eta_lam, eta_sigma, /*eta_b=*/0.0,
            DIST_HALFNORMAL, quad, K_max, r);
        NMixSiteResult res;
        res.log_lik         = d.log_lik;
        res.grad_eta_lambda = d.grad_eta_lambda;
        res.info_eta_lambda = d.info_eta_lambda;
        res.score_wt_lambda = d.score_wt_lambda;
        res.mean_N          = d.mean_N;
        res.var_N           = d.var_N;
        res.boundary_weight = d.boundary_weight;
        res.grad_eta_p.assign(1, d.grad_eta_d[0]);
        res.info_eta_p.assign(1, d.info_eta_d_fs[0][0]);
        res.grad_theta        = d.grad_theta;
        res.info_theta        = d.info_theta;
        res.info_lambda_theta = d.info_lambda_theta;
        res.cov_N_stheta      = d.cov_N_stheta;
        res.var_stheta        = d.var_stheta;
        return res;
    }
};

}  // namespace tulpaObs

#endif  // TULPAOBS_DISTANCE_RE_ORACLE_H
