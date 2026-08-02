// distance_re_oracle.cpp
// Constructor + Rcpp factory for DistanceGroupedOracle (declared in
// distance_re_oracle.h). The shared Z-sandwich / group-loop bodies live in
// count_grouped_oracle.cpp; this file builds the per-fit detection quadrature,
// stores the per-site bin counts, and exposes the factory.

#include "distance_re_oracle.h"
#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

DistanceGroupedOracle::DistanceGroupedOracle(
        int arm_,
        const Rcpp::IntegerMatrix& y_bins,
        const Rcpp::NumericMatrix& X_lambda,
        const Rcpp::NumericMatrix& X_sigma,
        const Rcpp::NumericMatrix& Z_site,
        const Rcpp::IntegerVector& site_group,
        int n_sites_, int n_groups_,
        SEXP quad_xptr, int K_max_, bool nb,
        int key_code_, double eta_b_) {
    K_max    = K_max_;
    key_code = key_code_;
    eta_b    = eta_b_;
    n_bins   = y_bins.ncol();

    // The detection design is per-site log-sigma: one "row" per site. Feed a
    // one-row-per-site long form (a dummy count vector + site_idx = 1..n_sites)
    // to the shared builder so Xp_site[i] is the i-th X_sigma row and the
    // per-site detection eta is the lone scalar eta_sigma_i.
    Rcpp::IntegerVector dummy_y(n_sites_, 0);
    Rcpp::IntegerVector site_idx(n_sites_);
    for (int i = 0; i < n_sites_; ++i) site_idx[i] = i + 1;
    std::vector<std::vector<int>> y_by_site_ignored;
    build_common(arm_, dummy_y, site_idx, X_lambda, X_sigma, Z_site,
                 site_group, n_sites_, n_groups_, nb, y_by_site_ignored);

    quad = *dist_quad_from_xptr(quad_xptr);

    y_bins_site.assign(n_sites_, std::vector<int>(n_bins, 0));
    for (int i = 0; i < n_sites_; ++i)
        for (int b = 0; b < n_bins; ++b) y_bins_site[i][b] = y_bins(i, b);
}

}  // namespace tulpaObs

// Rcpp factory: build the native single-species grouped-RE distance oracle and
// return it as an XPtr<tulpa::REGroupOracle>, consumed by tulpa::tulpa_re_aghq()
// exactly like the N-mixture / removal grouped oracles. Abundance-arm RE only
// (arm = 0). key = DIST_HALFNORMAL (0) or DIST_HAZARD (1); under the hazard key
// the shape eta_b is FIXED (the R wrapper profiles it over the outer log-marginal).
// [[Rcpp::export]]
SEXP cpp_distance_grouped_oracle(int arm,
                                 Rcpp::IntegerMatrix y_bins,
                                 Rcpp::NumericMatrix X_lambda,
                                 Rcpp::NumericMatrix X_sigma,
                                 Rcpp::NumericMatrix Z_site,
                                 Rcpp::IntegerVector site_group,
                                 int n_sites, int n_groups,
                                 SEXP quad_xptr, int K_max,
                                 bool nb = false, int key = 0, double eta_b = 0.0) {
    return Rcpp::XPtr<tulpa::REGroupOracle>(
        new tulpaObs::DistanceGroupedOracle(arm, y_bins, X_lambda, X_sigma,
                                            Z_site, site_group,
                                            n_sites, n_groups, quad_xptr,
                                            K_max, nb, key, eta_b),
        true);
}
