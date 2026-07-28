// nmix_re_oracle.h
// Native compiled REGroupOracle for the single-species N-mixture with
// site-level grouped random effects on EITHER the abundance (lambda) or
// detection (p) arm (gcol33/tulpaObs#13). The grouping factor is something
// other than species (station, observer-per-site, site cluster); the per-group
// RE vector b has dimension d = sum over RE-blocks of n_coefs (1..3), entering
// ONE arm via the per-site design row Z_i (length d).
//
// All of the Z-sandwich + group-loop + REGroupOracle plumbing is shared with
// the other count-model grouped-RE oracles in CountGroupedOracle
// (count_grouped_oracle.h); this class only supplies the per-site N-mixture
// marginal (nmix_kernel.h, cached lgamma) and the per-site cache.
//
// Poisson OR negative-binomial abundance. Under NB the global dispersion size
// r is the (d+1)-th theta entry log_r (carried as theta[n_theta - 1] when
// is_nb); Poisson is the r = +Inf limit, threaded through nmix_kernel.h.

#ifndef TULPAOBS_NMIX_RE_ORACLE_H
#define TULPAOBS_NMIX_RE_ORACLE_H

#include "count_grouped_oracle.h"
#include "nmix_kernel.h"
#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

struct NMixGroupedOracle : CountGroupedOracle {
    std::vector<NMixSiteCache> site_cache;  // eta-independent lgamma terms

    NMixGroupedOracle(int arm_,
                      const Rcpp::IntegerVector& y,
                      const Rcpp::IntegerVector& site_idx,
                      const Rcpp::NumericMatrix& X_lambda,
                      const Rcpp::NumericMatrix& X_p,
                      const Rcpp::NumericMatrix& Z_site,
                      const Rcpp::IntegerVector& site_group,
                      int n_sites_, int n_groups_, int K_max,
                      bool nb, int headroom = -1);

    NMixSiteResult eval_site(int i, const double* eta_p_ptr,
                             double eta_lam) const override {
        return compute_nmix_site_cached(site_cache[i], eta_p_ptr, eta_lam, r);
    }
};

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_RE_ORACLE_H
