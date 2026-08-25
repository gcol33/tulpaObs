// long_form_group.h
// Grouping of a long-form (species, site) response into per-(species, site)
// row lists, shared by the two community N-mixture entry points that build one:
// the AGHQ oracle (nmix_community_oracle.cpp) and the shared-field EM
// (nmix_community_field.cpp).
//
// The grouping itself is three lines; what it needs is the guard in front of
// it. The buffer is sized from n_species / n_sites, which arrive as plain int
// arguments, and it is indexed with values taken from species_idx / site_idx,
// which arrive separately. Neither R caller derives one from the other -- both
// forward the counts from further up -- so an out-of-range code is a push_back
// on a std::vector object read from arbitrary memory. Relating the two here
// means every caller inherits the check rather than restating it.

#ifndef TULPAOBS_LONG_FORM_GROUP_H
#define TULPAOBS_LONG_FORM_GROUP_H

#include "tobs_shape.h"
#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

// rows[s][i] = the long-form row positions of species s at site i, in input
// order. `n_obs` is the response length every index vector must match.
inline std::vector<std::vector<std::vector<int>>>
group_rows_by_species_site(const Rcpp::IntegerVector& species_idx,
                           const Rcpp::IntegerVector& site_idx,
                           int n_species, int n_sites, R_xlen_t n_obs) {
    namespace sh = tulpaObs::shape;
    sh::check_dim_arg(n_species, "n_species");
    sh::check_dim_arg(n_sites, "n_sites");
    sh::check_len(species_idx, n_obs, "species_idx");
    sh::check_len(site_idx, n_obs, "site_idx");
    sh::check_index1(species_idx, n_species, "species_idx");
    sh::check_index1(site_idx, n_sites, "site_idx");

    std::vector<std::vector<std::vector<int>>> rows(
        n_species, std::vector<std::vector<int>>(n_sites));
    for (R_xlen_t r = 0; r < n_obs; ++r)
        rows[species_idx[r] - 1][site_idx[r] - 1].push_back((int) r);
    return rows;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_LONG_FORM_GROUP_H
