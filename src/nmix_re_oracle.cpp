// nmix_re_oracle.cpp
// Constructor + Rcpp factory for NMixGroupedOracle (declared in
// nmix_re_oracle.h). The shared Z-sandwich / group-loop bodies live in
// count_grouped_oracle.cpp; this file only builds the per-site N-mixture cache
// and exposes the factory.

#include "nmix_re_oracle.h"
#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

NMixGroupedOracle::NMixGroupedOracle(int arm_,
                                     const Rcpp::IntegerVector& y,
                                     const Rcpp::IntegerVector& site_idx,
                                     const Rcpp::NumericMatrix& X_lambda,
                                     const Rcpp::NumericMatrix& X_p,
                                     const Rcpp::NumericMatrix& Z_site,
                                     const Rcpp::IntegerVector& site_group,
                                     int n_sites_, int n_groups_, int K_max,
                                     bool nb) {
    std::vector<std::vector<int>> y_by_site;
    build_common(arm_, y, site_idx, X_lambda, X_p, Z_site, site_group,
                 n_sites_, n_groups_, nb, y_by_site);

    site_cache.assign(n_sites, NMixSiteCache());
    for (int i = 0; i < n_sites; ++i) {
        const int J = (int)y_by_site[i].size();
        site_cache[i] = nmix_precompute_site(y_by_site[i].data(), J, K_max);
    }
}

}  // namespace tulpaObs

// Rcpp factory: build the native single-species grouped-RE N-mixture oracle and
// return it as an XPtr<tulpa::REGroupOracle>. tulpa::tulpa_re_aghq() consumes
// the XPtr through the engine's REGroupOracle interface; the per-group
// marginal / score / observed-info assembly runs entirely in tulpaObs.
// [[Rcpp::export]]
SEXP cpp_nmix_grouped_oracle(int arm,
                             Rcpp::IntegerVector y,
                             Rcpp::IntegerVector site_idx,
                             Rcpp::NumericMatrix X_lambda,
                             Rcpp::NumericMatrix X_p,
                             Rcpp::NumericMatrix Z_site,
                             Rcpp::IntegerVector site_group,
                             int n_sites, int n_groups, int K_max,
                             bool nb = false) {
    return Rcpp::XPtr<tulpa::REGroupOracle>(
        new tulpaObs::NMixGroupedOracle(arm, y, site_idx, X_lambda, X_p,
                                        Z_site, site_group,
                                        n_sites, n_groups, K_max, nb),
        true);
}
