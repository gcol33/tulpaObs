// nmix_community_spatial_oracle.h
// Native compiled REGroupOracle for the SPATIAL community / multispecies
// N-mixture (spAbundance sfMsNMix): the community N-mixture (per-species Royle
// 2004 N-mixture with Gaussian community hyperpriors on the per-species
// abundance / detection coefficients) plus a SHARED areal field offset on the
// abundance arm,
//
//   N_{s,i}        ~ Poisson(lambda_{s,i})
//   y_{s,i,j} | N  ~ Binomial(N_{s,i}, p_{s,i,j})
//   log lambda_{s,i} = X_lambda_i . (mu_lambda + b_lambda_s) + offset_i
//   logit p_{s,i,j}  = X_p_{ij}   . (mu_p     + b_p_s)
//   b_lambda_s ~ N(0, Sigma_lambda),  b_p_s ~ N(0, Sigma_p)
//
// where offset_i = sigma * f_i is the per-site value of the shared ICAR / BYM2
// field, SAME across all species. The grouping factor is the species; the
// per-group RE vector is b_s = (b_lambda_s, b_p_s) (dimension d = p_lambda +
// p_p), entering the abundance / detection linear predictors as coef = mu + b_s.
//
// This is the spatial sibling of NMixCommunityOracle: identical per-species
// marginal / score / observed-info assembly (the per-site marginal kernel is
// the shared nmix_kernel.h), with the field offset added to eta_lambda BEFORE
// the per-site marginal. The shared field f / its hyperparameter (tau for ICAR,
// sigma + rho for BYM2) is conditioned at the value the OUTER nested-Laplace
// driver is integrating; the offset is set on the oracle (set_offset) before
// each tulpa::tulpa_re_aghq() community solve at that node. The community means
// + RE covariance are then integrated by the engine GIVEN the offset, and the
// field is integrated by the outer grid (the nested-approx + debias split: the
// shared field block is an outer-integrated latent Gaussian, the community RE
// block is the inner AGHQ-integrated quantity).
//
// Poisson only (the abundance mixing distribution of the field model). The
// global-NB extension is the non-spatial community path's territory.
//
// Bodies are out-of-line in nmix_community_spatial_oracle.cpp (the Eigen
// per-site assembly is heavy template code; keeping it out of the header stops
// a calling TU from re-instantiating it inline, which overflows MinGW g++ -O2).

#ifndef TULPAOBS_NMIX_COMMUNITY_SPATIAL_ORACLE_H
#define TULPAOBS_NMIX_COMMUNITY_SPATIAL_ORACLE_H

#include "tulpa/aghq_oracle.h"
#include "nmix_kernel.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

namespace tulpaObs {

struct SpatialNMixCommunityOracle : tulpa::REGroupOracle {
    int p_lam = 0, p_p = 0;
    Eigen::MatrixXd Xlam;             // n_sites x p_lambda (shared across species)
    Eigen::VectorXd mu;               // active community means (theta), length d
    Eigen::VectorXd offset;           // per-site abundance offset (= sigma * f), n_sites
    int n_sites = 0;

    // Per species, per site: the cached Poisson marginal (lgamma precompute) and
    // the detection design rows for that site's visits, in input order.
    struct SiteRec {
        int site = 0;                 // 0-based row into Xlam / offset
        NMixSiteCache cache;          // eta-independent lgamma terms
        Eigen::MatrixXd Xp;           // J_i x p_p
    };
    std::vector<std::vector<SiteRec>> sp_sites;   // [n_species][n_sites]

    struct SpeciesEval {
        double logL = 0.0;
        Eigen::VectorXd grad;         // d ell_g / db  (== d ell_g / d mu)
        Eigen::MatrixXd negH;         // -d^2 ell_g / db^2 (marginal observed info)
        Eigen::MatrixXd fisher;       // complete-data Fisher (PSD)
    };

    SpatialNMixCommunityOracle(const Rcpp::IntegerVector& y,
                               const Rcpp::IntegerVector& site_idx,
                               const Rcpp::IntegerVector& species_idx,
                               const Rcpp::NumericMatrix& X_lambda,
                               const Rcpp::NumericMatrix& X_p,
                               int n_sites, int n_species, int K_max);

    static inline double clamp30(double e) {
        return e < -30.0 ? -30.0 : (e > 30.0 ? 30.0 : e);
    }

    // Set the shared per-site abundance offset (length n_sites). Called by the
    // R driver before each community solve at a field node.
    void set_offset(const Rcpp::NumericVector& z);

    SpeciesEval eval_species(int g, const double* b,
                             bool want_negH = true, bool want_fisher = true) const;

    void rebind(const double* theta) override {
        for (int i = 0; i < d; ++i) mu(i) = theta[i];
    }
    void grad_hess(int g, const double* b, double& logL,
                   double* grad, double* negH) const override;
    void node_ll(int g, const double* B, int n_nodes, double* out) const override;
    void theta_score(int g, const double* b, double* dl_dtheta) const override;
    bool newton_hess(int g, const double* b, double* H) const override;
    bool thread_safe() const override { return true; }
};

}  // namespace tulpaObs

#endif  // TULPAOBS_NMIX_COMMUNITY_SPATIAL_ORACLE_H
