// removal_re_oracle.cpp
// Constructor + Rcpp factory for RemovalGroupedOracle (declared in
// removal_re_oracle.h). The shared Z-sandwich / group-loop bodies live in
// count_grouped_oracle.cpp; this file only stores the per-site pass removals
// and exposes the factory.

#include "removal_re_oracle.h"
#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

RemovalGroupedOracle::RemovalGroupedOracle(int arm_,
                                           const Rcpp::IntegerVector& y,
                                           const Rcpp::IntegerVector& site_idx,
                                           const Rcpp::NumericMatrix& X_lambda,
                                           const Rcpp::NumericMatrix& X_p,
                                           const Rcpp::NumericMatrix& Z_site,
                                           const Rcpp::IntegerVector& site_group,
                                           int n_sites_, int n_groups_,
                                           int K_max_, bool nb) {
    K_max = K_max_;
    build_common(arm_, y, site_idx, X_lambda, X_p, Z_site, site_group,
                 n_sites_, n_groups_, nb, y_site);
}

}  // namespace tulpaObs

// Rcpp factory: build the native single-species grouped-RE removal oracle and
// return it as an XPtr<tulpa::REGroupOracle>, consumed by tulpa::tulpa_re_aghq()
// exactly like the N-mixture grouped oracle.
// [[Rcpp::export]]
SEXP cpp_removal_grouped_oracle(int arm,
                                Rcpp::IntegerVector y,
                                Rcpp::IntegerVector site_idx,
                                Rcpp::NumericMatrix X_lambda,
                                Rcpp::NumericMatrix X_p,
                                Rcpp::NumericMatrix Z_site,
                                Rcpp::IntegerVector site_group,
                                int n_sites, int n_groups, int K_max,
                                bool nb = false) {
    return Rcpp::XPtr<tulpa::REGroupOracle>(
        new tulpaObs::RemovalGroupedOracle(arm, y, site_idx, X_lambda, X_p,
                                           Z_site, site_group,
                                           n_sites, n_groups, K_max, nb),
        true);
}
