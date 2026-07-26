// nmix_community_spatial_oracle.cpp
// R entry points for the SPATIAL community / multispecies N-mixture oracle
// (spAbundance sfMsNMix). The oracle itself IS NMixCommunityOracle: the spatial
// model differs from the non-spatial one only in seeding eta_lambda with the
// shared per-site field offset, which that class carries as an optional
// `offset` (nmix_community_oracle.h). Poisson only -- the global-NB extension
// belongs to the non-spatial community path.

#include "nmix_community_oracle.h"
#include <Rcpp.h>

// [[Rcpp::export]]
SEXP cpp_nmix_spatial_community_oracle(Rcpp::IntegerVector y,
                                       Rcpp::IntegerVector site_idx,
                                       Rcpp::IntegerVector species_idx,
                                       Rcpp::NumericMatrix X_lambda,
                                       Rcpp::NumericMatrix X_p,
                                       int n_sites, int n_species, int K_max) {
    return Rcpp::XPtr<tulpa::REGroupOracle>(
        new tulpaObs::NMixCommunityOracle(y, site_idx, species_idx,
                                          X_lambda, X_p,
                                          n_sites, n_species, K_max),
        true);
}

// Set the shared per-site abundance offset (= sigma * f). The XPtr is the same
// REGroupOracle the AGHQ engine drives; we downcast to the concrete type to
// reach the offset setter.
// [[Rcpp::export]]
void cpp_nmix_spatial_community_set_offset(SEXP oracle_ptr,
                                           Rcpp::NumericVector z) {
    Rcpp::XPtr<tulpa::REGroupOracle> xp(oracle_ptr);
    tulpaObs::NMixCommunityOracle* orc =
        dynamic_cast<tulpaObs::NMixCommunityOracle*>(xp.get());
    if (orc == nullptr)
        Rcpp::stop("cpp_nmix_spatial_community_set_offset: pointer is not a "
                   "community N-mixture oracle.");
    orc->set_offset(z);
}
