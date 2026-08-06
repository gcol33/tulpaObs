// community_grid_pack.h
// Packing the per-outer-grid-point community Laplace-EM results into the list
// the R side reads, shared by the community areal drivers (community occupancy,
// ms_occu_spatial.cpp; community N-mixture, nmix_community_spatial.cpp).
//
// Every driver produces the same per-point record: the marginal and data
// log-likelihood, the mode (community means followed by the field), the
// iteration record, the community-mean covariance, and per-arm community
// covariances and BLUPs. Two things differ. What the state arm is called
// ("psi" for occupancy, "lambda" for abundance), which is a naming choice on
// the returned list, so the members are reached through pointers-to-member and
// the names are passed in. And whether the family carries a per-point boundary
// diagnostic, which is a null member pointer where it does not.
//
// Include AFTER RcppEigen.h: Rcpp::wrap on an Eigen matrix is looked up when
// this header is parsed.

#ifndef TULPAOBS_COMMUNITY_GRID_PACK_H
#define TULPAOBS_COMMUNITY_GRID_PACK_H

#include <Rcpp.h>
#include <vector>

namespace tulpaObs {

template <class Result, class Mat>
Rcpp::List community_pack_grid(int d, int field_len, int n_grid,
                               const std::vector<Result>& results,
                               const Rcpp::NumericMatrix& theta_grid_out,
                               Mat Result::*sigma_state,
                               Mat Result::*blup_state,
                               const char* sigma_state_name,
                               const char* blup_state_name,
                               const char* p_state_name,
                               int p_state, int p_p, int n_spatial,
                               double Result::*boundary = nullptr) {
    Rcpp::NumericVector log_marginals(n_grid), log_liks(n_grid);
    Rcpp::NumericVector boundaries(boundary ? n_grid : 0);
    Rcpp::IntegerVector n_iters(n_grid);
    Rcpp::LogicalVector convergeds(n_grid);
    Rcpp::NumericMatrix modes(n_grid, d + field_len);   // (mu, field) per grid point
    Rcpp::List vcov_mu(n_grid), Sigma_st(n_grid), Sigma_p(n_grid);
    Rcpp::List blup_st(n_grid), blup_p(n_grid);
    for (int k = 0; k < n_grid; ++k) {
        const Result& rr = results[k];
        log_marginals[k] = rr.log_marginal;
        log_liks[k]      = rr.log_lik;
        n_iters[k]       = rr.n_iter;
        convergeds[k]    = rr.converged;
        if (boundary) boundaries[k] = rr.*boundary;
        for (int j = 0; j < d; ++j)         modes(k, j) = rr.mu(j);
        for (int j = 0; j < field_len; ++j) modes(k, d + j) = rr.field(j);
        vcov_mu[k]  = Rcpp::wrap(rr.vcov_mu);
        Sigma_st[k] = Rcpp::wrap(rr.*sigma_state);
        Sigma_p[k]  = Rcpp::wrap(rr.Sigma_p);
        blup_st[k]  = Rcpp::wrap(rr.*blup_state);
        blup_p[k]   = Rcpp::wrap(rr.blup_p);
    }
    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("theta_grid")     = theta_grid_out,
        Rcpp::Named("log_marginal")   = log_marginals,
        Rcpp::Named("log_lik")        = log_liks,
        Rcpp::Named("modes")          = modes,
        Rcpp::Named("vcov_mu")        = vcov_mu,
        Rcpp::Named(sigma_state_name) = Sigma_st,
        Rcpp::Named("Sigma_p")        = Sigma_p,
        Rcpp::Named(blup_state_name)  = blup_st,
        Rcpp::Named("b_p")            = blup_p,
        Rcpp::Named("n_iter")         = n_iters,
        Rcpp::Named("converged")      = convergeds);
    if (boundary) out["boundary_max"] = boundaries;
    out[p_state_name] = p_state;
    out["p_p"]        = p_p;
    out["n_spatial"]  = n_spatial;
    out["n_grid"]     = n_grid;
    return out;
}

}  // namespace tulpaObs

#endif  // TULPAOBS_COMMUNITY_GRID_PACK_H
