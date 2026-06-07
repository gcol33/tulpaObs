// removal_re_oracle.h
// Native compiled REGroupOracle for the single-species removal-sampling
// (sequential-depletion) abundance model with a site-level grouped random
// effect on EITHER the abundance (lambda) or detection (p) arm
// (gcol33/tulpaObs#51). Shares the entire Z-sandwich + group-loop + REGroupOracle
// plumbing with the N-mixture grouped oracle via CountGroupedOracle; the only
// removal-specific pieces are the per-site marginal (removal_kernel.h, which
// depletes the available count per pass) and the per-site pass data.
//
// The detection design rows MUST be in pass order per site (depletion
// accumulates earlier passes); the long-form row order is preserved by
// build_common, so the caller orders the passes. Poisson OR negative-binomial
// abundance, the same r = +Inf limit as the N-mixture.

#ifndef TULPAOBS_REMOVAL_RE_ORACLE_H
#define TULPAOBS_REMOVAL_RE_ORACLE_H

#include "count_grouped_oracle.h"
#include "removal_kernel.h"
#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

struct RemovalGroupedOracle : CountGroupedOracle {
    std::vector<std::vector<int>> y_site;  // per-site removals, pass order
    int K_max = 0;

    RemovalGroupedOracle(int arm_,
                         const Rcpp::IntegerVector& y,
                         const Rcpp::IntegerVector& site_idx,
                         const Rcpp::NumericMatrix& X_lambda,
                         const Rcpp::NumericMatrix& X_p,
                         const Rcpp::NumericMatrix& Z_site,
                         const Rcpp::IntegerVector& site_group,
                         int n_sites_, int n_groups_, int K_max_,
                         bool nb);

    NMixSiteResult eval_site(int i, const double* eta_p_ptr,
                             double eta_lam) const override {
        const int n_pass = (int)y_site[i].size();
        return compute_removal_site(y_site[i].data(), eta_p_ptr, n_pass,
                                    eta_lam, K_max, r);
    }
};

}  // namespace tulpaObs

#endif  // TULPAOBS_REMOVAL_RE_ORACLE_H
