// nmix_community_spatial_oracle.cpp
// Out-of-line bodies for SpatialNMixCommunityOracle (declared in the header).
// The per-site marginal math is nmix_kernel.h (the single source); this file is
// the design-sandwiched per-species b-space assembly with the shared field
// offset added to the abundance predictor.

#include "nmix_community_spatial_oracle.h"
#include "nmix_oracle_emit.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <cstddef>
#include <vector>

namespace tulpaObs {

SpatialNMixCommunityOracle::SpatialNMixCommunityOracle(
    const Rcpp::IntegerVector& y,
    const Rcpp::IntegerVector& site_idx,
    const Rcpp::IntegerVector& species_idx,
    const Rcpp::NumericMatrix& X_lambda,
    const Rcpp::NumericMatrix& X_p,
    int n_sites_, int n_species, int K_max) {
    p_lam    = X_lambda.ncol();
    p_p      = X_p.ncol();
    d        = p_lam + p_p;
    n_theta  = d;                    // Poisson: no global dispersion entry
    n_groups = n_species;
    n_sites  = n_sites_;
    mu       = Eigen::VectorXd::Zero(d);
    offset   = Eigen::VectorXd::Zero(n_sites);

    Xlam.resize(n_sites, p_lam);
    for (int i = 0; i < n_sites; ++i)
        for (int c = 0; c < p_lam; ++c) Xlam(i, c) = X_lambda(i, c);

    const int n_obs = y.size();
    std::vector<std::vector<std::vector<int>>> rows(
        n_species, std::vector<std::vector<int>>(n_sites));
    for (int r = 0; r < n_obs; ++r)
        rows[species_idx[r] - 1][site_idx[r] - 1].push_back(r);

    sp_sites.assign(n_species, std::vector<SiteRec>());
    for (int s = 0; s < n_species; ++s) {
        sp_sites[s].reserve(n_sites);
        for (int i = 0; i < n_sites; ++i) {
            SiteRec rec;
            rec.site = i;
            const std::vector<int>& rr = rows[s][i];
            const int J = (int)rr.size();
            std::vector<int> yv(J);
            rec.Xp.resize(J, p_p);
            for (int j = 0; j < J; ++j) {
                yv[j] = y[rr[j]];
                for (int c = 0; c < p_p; ++c) rec.Xp(j, c) = X_p(rr[j], c);
            }
            rec.cache = nmix_precompute_site(yv.data(), J, K_max);
            sp_sites[s].push_back(std::move(rec));
        }
    }
}

void SpatialNMixCommunityOracle::set_offset(const Rcpp::NumericVector& z) {
    if ((int)z.size() != n_sites)
        Rcpp::stop("set_offset: length(z) must equal n_sites.");
    for (int i = 0; i < n_sites; ++i) offset(i) = z[i];
}

SpatialNMixCommunityOracle::SpeciesEval
SpatialNMixCommunityOracle::eval_species(int g, const double* b,
                                         bool want_negH, bool want_fisher) const {
    SpeciesEval e;
    e.grad = Eigen::VectorXd::Zero(d);
    if (want_negH)   e.negH   = Eigen::MatrixXd::Zero(d, d);
    if (want_fisher) e.fisher = Eigen::MatrixXd::Zero(d, d);

    Eigen::VectorXd coef(d);
    for (int i = 0; i < d; ++i) coef(i) = mu(i) + b[i];

    std::vector<double> eta_p;
    for (const SiteRec& rec : sp_sites[g]) {
        const int J = rec.cache.n_visits;
        double eta_lam = offset(rec.site);          // shared field offset
        for (int c = 0; c < p_lam; ++c) eta_lam += Xlam(rec.site, c) * coef(c);
        eta_lam = clamp30(eta_lam);

        eta_p.assign(J, 0.0);
        for (int j = 0; j < J; ++j) {
            double v = 0.0;
            for (int c = 0; c < p_p; ++c) v += rec.Xp(j, c) * coef(p_lam + c);
            eta_p[j] = clamp30(v);
        }

        const NMixSiteResult res =
            compute_nmix_site_cached(rec.cache, eta_p.data(), eta_lam);
        e.logL += res.log_lik;

        for (int c = 0; c < p_lam; ++c)
            e.grad(c) += Xlam(rec.site, c) * res.grad_eta_lambda;
        for (int j = 0; j < J; ++j)
            for (int c = 0; c < p_p; ++c)
                e.grad(p_lam + c) += rec.Xp(j, c) * res.grad_eta_p[j];

        if (!want_negH && !want_fisher) continue;

        const int dd = 1 + J;
        Eigen::MatrixXd Bobs, Bfis;
        if (want_negH)   { Bobs = Eigen::MatrixXd::Zero(dd, dd); Bobs(0, 0) = res.info_eta_lambda; }
        if (want_fisher) { Bfis = Eigen::MatrixXd::Zero(dd, dd); Bfis(0, 0) = res.info_eta_lambda; }
        for (int j = 0; j < J; ++j) {
            if (want_negH)   Bobs(1 + j, 1 + j) = res.info_eta_p[j];
            if (want_fisher) Bfis(1 + j, 1 + j) = res.info_eta_p[j];
        }
        if (want_negH && J > 0) {
            Eigen::VectorXd vv(dd);
            vv(0) = -res.score_wt_lambda;
            for (int j = 0; j < J; ++j) {
                const double pj = (eta_p[j] > 0.0)
                    ? 1.0 / (1.0 + std::exp(-eta_p[j]))
                    : std::exp(eta_p[j]) / (1.0 + std::exp(eta_p[j]));
                vv(1 + j) = pj;
            }
            Bobs.noalias() -= res.var_N * (vv * vv.transpose());
        }

        Eigen::MatrixXd Zi = Eigen::MatrixXd::Zero(dd, d);
        for (int c = 0; c < p_lam; ++c) Zi(0, c) = Xlam(rec.site, c);
        for (int j = 0; j < J; ++j)
            for (int c = 0; c < p_p; ++c) Zi(1 + j, p_lam + c) = rec.Xp(j, c);

        const Eigen::MatrixXd Zt = Zi.transpose();
        if (want_negH) {
            const Eigen::MatrixXd BZ = Bobs * Zi;
            e.negH.noalias() += Zt * BZ;
        }
        if (want_fisher) {
            const Eigen::MatrixXd BZ = Bfis * Zi;
            e.fisher.noalias() += Zt * BZ;
        }
    }
    return e;
}

void SpatialNMixCommunityOracle::grad_hess(int g, const double* b, double& logL,
                                           double* grad, double* negH) const {
    emit_grad_hess(eval_species(g, b, /*want_negH=*/true, /*want_fisher=*/false),
                   d, logL, grad, negH);
}

void SpatialNMixCommunityOracle::node_ll(int g, const double* B, int n_nodes,
                                         double* out) const {
    Eigen::VectorXd coef(d);
    std::vector<double> eta_p;
    for (int k = 0; k < n_nodes; ++k) {
        const double* bk = B + (std::size_t)k * d;
        for (int i = 0; i < d; ++i) coef(i) = mu(i) + bk[i];
        double ll = 0.0;
        for (const SiteRec& rec : sp_sites[g]) {
            const int J = rec.cache.n_visits;
            double eta_lam = offset(rec.site);
            for (int c = 0; c < p_lam; ++c) eta_lam += Xlam(rec.site, c) * coef(c);
            eta_lam = clamp30(eta_lam);
            eta_p.assign(J, 0.0);
            for (int j = 0; j < J; ++j) {
                double v = 0.0;
                for (int c = 0; c < p_p; ++c) v += rec.Xp(j, c) * coef(p_lam + c);
                eta_p[j] = clamp30(v);
            }
            ll += compute_nmix_site_cached(rec.cache, eta_p.data(), eta_lam).log_lik;
        }
        out[k] = ll;
    }
}

void SpatialNMixCommunityOracle::theta_score(int g, const double* b,
                                             double* dl_dtheta) const {
    const SpeciesEval e = eval_species(g, b, /*want_negH=*/false, /*want_fisher=*/false);
    for (int i = 0; i < d; ++i) dl_dtheta[i] = e.grad(i);
}

bool SpatialNMixCommunityOracle::newton_hess(int g, const double* b, double* H) const {
    emit_fisher(eval_species(g, b, /*want_negH=*/false, /*want_fisher=*/true), d, H);
    return true;
}

}  // namespace tulpaObs

// Rcpp factory: build the native spatial community N-mixture oracle and return
// it as an XPtr<tulpa::REGroupOracle>. tulpa::tulpa_re_aghq() consumes the XPtr;
// the per-species marginal / score / observed-info assembly runs in tulpaObs.
// The shared field offset is set separately via cpp_nmix_spatial_community_set_offset.
// [[Rcpp::export]]
SEXP cpp_nmix_spatial_community_oracle(Rcpp::IntegerVector y,
                                       Rcpp::IntegerVector site_idx,
                                       Rcpp::IntegerVector species_idx,
                                       Rcpp::NumericMatrix X_lambda,
                                       Rcpp::NumericMatrix X_p,
                                       int n_sites, int n_species, int K_max) {
    return Rcpp::XPtr<tulpa::REGroupOracle>(
        new tulpaObs::SpatialNMixCommunityOracle(y, site_idx, species_idx,
                                                 X_lambda, X_p,
                                                 n_sites, n_species, K_max),
        true);
}

// Set the shared per-site abundance offset (= sigma * f) on a spatial community
// oracle. The XPtr is the same REGroupOracle the AGHQ engine drives; we
// downcast to the concrete type to reach the offset setter.
// [[Rcpp::export]]
void cpp_nmix_spatial_community_set_offset(SEXP oracle_ptr,
                                           Rcpp::NumericVector z) {
    Rcpp::XPtr<tulpa::REGroupOracle> xp(oracle_ptr);
    tulpaObs::SpatialNMixCommunityOracle* orc =
        dynamic_cast<tulpaObs::SpatialNMixCommunityOracle*>(xp.get());
    if (orc == nullptr)
        Rcpp::stop("cpp_nmix_spatial_community_set_offset: pointer is not a "
                   "SpatialNMixCommunityOracle.");
    orc->set_offset(z);
}
