// nmix_re_oracle.cpp
// Out-of-line bodies for NMixGroupedOracle (declared in nmix_re_oracle.h).
// The per-site marginal math is nmix_kernel.h (the single source) -- this
// file is the Z-sandwich that maps per-site eta-space derivatives back to the
// per-group RE coordinate b through the site-level RE design row Z_i.
// Out-of-line so a translation unit that only calls the oracle does not
// re-instantiate the per-site Eigen assembly inline, matching the layout of
// nmix_community_oracle.cpp.

#include "nmix_re_oracle.h"
#include "nmix_oracle_emit.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <cstddef>
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
    arm      = arm_;
    p_lam    = X_lambda.ncol();
    p_p      = X_p.ncol();
    d        = Z_site.ncol();
    is_nb    = nb;
    n_theta  = p_lam + p_p + (nb ? 1 : 0);
    n_groups = n_groups_;
    n_sites  = n_sites_;

    Xlam.resize(n_sites, p_lam);
    for (int i = 0; i < n_sites; ++i)
        for (int c = 0; c < p_lam; ++c) Xlam(i, c) = X_lambda(i, c);

    Zsite.resize(n_sites, d);
    for (int i = 0; i < n_sites; ++i)
        for (int c = 0; c < d; ++c) Zsite(i, c) = Z_site(i, c);

    // Group long-form rows by site, preserving input visit order. n_sites is
    // the design granularity; sites without observed visits get integer(0)
    // (they carry no marginal information; the kernel returns log_lik = 0).
    const int n_obs = y.size();
    std::vector<std::vector<int>> rows_at_site(n_sites);
    for (int r = 0; r < n_obs; ++r)
        rows_at_site[site_idx[r] - 1].push_back(r);

    Xp_site.assign(n_sites, Eigen::MatrixXd());
    site_cache.assign(n_sites, NMixSiteCache());
    for (int i = 0; i < n_sites; ++i) {
        const std::vector<int>& rr = rows_at_site[i];
        const int J = (int)rr.size();
        std::vector<int> yv(J);
        Xp_site[i].resize(J, p_p);
        for (int j = 0; j < J; ++j) {
            yv[j] = y[rr[j]];
            for (int c = 0; c < p_p; ++c) Xp_site[i](j, c) = X_p(rr[j], c);
        }
        site_cache[i] = nmix_precompute_site(yv.data(), J, K_max);
    }

    // Per-group site lists (0-based site indices).
    sites_by_group.assign(n_groups, std::vector<int>());
    for (int i = 0; i < n_sites; ++i) {
        const int g = site_group[i] - 1;
        if (g >= 0 && g < n_groups) sites_by_group[g].push_back(i);
    }

    eta_lambda_base = Eigen::VectorXd::Zero(n_sites);
    eta_p_site.assign(n_sites, std::vector<double>());
    for (int i = 0; i < n_sites; ++i) eta_p_site[i].assign((std::size_t)Xp_site[i].rows(), 0.0);
}

void NMixGroupedOracle::rebind(const double* theta) {
    // theta layout: [beta_lambda (p_lam) | beta_p (p_p) | log_r if NB].
    // eta_lambda_base = clamp30(Xlam * beta_lambda) -- site-level abundance
    // fixed-effect predictor at every site; the per-group call adds the
    // RE shift on the lambda arm (or evaluates this as-is for the p arm).
    // eta_p_site[i] = clamp30(Xp_i * beta_p) per visit.
    for (int i = 0; i < n_sites; ++i) {
        double e = 0.0;
        for (int c = 0; c < p_lam; ++c) e += Xlam(i, c) * theta[c];
        eta_lambda_base(i) = clamp30(e);
        const int J = (int)Xp_site[i].rows();
        for (int j = 0; j < J; ++j) {
            double v = 0.0;
            for (int c = 0; c < p_p; ++c) v += Xp_site[i](j, c) * theta[p_lam + c];
            eta_p_site[i][j] = clamp30(v);
        }
    }
    if (is_nb) r = std::exp(theta[p_lam + p_p]);
    else       r = std::numeric_limits<double>::infinity();
}

NMixGroupedOracle::GroupEval
NMixGroupedOracle::eval_group(int g, const double* b,
                              bool want_negH, bool want_fisher,
                              bool want_theta_grad) const {
    GroupEval e;
    e.grad = Eigen::VectorXd::Zero(d);
    if (want_negH)        e.negH       = Eigen::MatrixXd::Zero(d, d);
    if (want_fisher)      e.fisher     = Eigen::MatrixXd::Zero(d, d);
    if (want_theta_grad)  e.theta_grad = Eigen::VectorXd::Zero(n_theta);

    std::vector<double> eta_p_shift;
    for (const int i : sites_by_group[g]) {
        // Z_i . b -- the per-site shift the engine applies to the active arm.
        double shift = 0.0;
        for (int c = 0; c < d; ++c) shift += Zsite(i, c) * b[c];

        const int J = (int)Xp_site[i].rows();
        double eta_lam;
        const double* eta_p_ptr;
        if (arm == 0) {
            // lambda arm: shifted predictor enters at the site level; the R
            // closure clamps eta_lambda_base[i] + shift, the kernel sees the
            // clamped value. Detection stays at the FE base.
            eta_lam = clamp30(eta_lambda_base(i) + shift);
            eta_p_ptr = eta_p_site[i].empty() ? nullptr : eta_p_site[i].data();
        } else {
            // p arm: the shift is a per-site scalar applied uniformly to every
            // visit at that site, matching the R closure's apply_p_shift. The
            // R path clamps the shift, then adds (no second clamp on the sum).
            eta_lam = eta_lambda_base(i);
            const double shift_cl = clamp30(shift);
            eta_p_shift.assign(J, 0.0);
            for (int j = 0; j < J; ++j) eta_p_shift[j] = eta_p_site[i][j] + shift_cl;
            eta_p_ptr = eta_p_shift.empty() ? nullptr : eta_p_shift.data();
        }

        const NMixSiteResult res =
            compute_nmix_site_cached(site_cache[i], eta_p_ptr, eta_lam, r);
        e.logL += res.log_lik;

        // ---- b-space score ----
        // arm == lambda: d ell_i / d shift_i = grad_eta_lambda
        // arm == p:      d ell_i / d shift_i = sum_j grad_eta_p[j]
        // Chain to b via shift = Z_i . b -> d ell / db = (d ell / d shift) Z_i.
        double s_d1;
        if (arm == 0) {
            s_d1 = res.grad_eta_lambda;
        } else {
            s_d1 = 0.0;
            for (int j = 0; j < J; ++j) s_d1 += res.grad_eta_p[j];
        }
        for (int c = 0; c < d; ++c) e.grad(c) += s_d1 * Zsite(i, c);

        // ---- b-space curvatures (rank-1 in Z_i: outer product * scalar) ----
        // The eta-space per-site marginal block contracts to a scalar through
        // the ones vector on the p arm. Observed info (negH) carries the
        // Var(N|y) abundance-detection coupling; the Newton Fisher does not.
        if (want_negH || want_fisher) {
            double b_obs = 0.0, b_fis = 0.0;
            if (arm == 0) {
                // (1,1) entry of the per-site B_i; var_N coupling via score_wt^2.
                b_fis = res.info_eta_lambda;
                b_obs = res.info_eta_lambda
                      - res.var_N * res.score_wt_lambda * res.score_wt_lambda;
            } else {
                double sum_info = 0.0;
                for (int j = 0; j < J; ++j) sum_info += res.info_eta_p[j];
                b_fis = sum_info;
                if (J > 0) {
                    double sum_p = 0.0;
                    for (int j = 0; j < J; ++j) {
                        const double ej = eta_p_shift[j];
                        const double pj = (ej > 0.0)
                            ? 1.0 / (1.0 + std::exp(-ej))
                            : std::exp(ej) / (1.0 + std::exp(ej));
                        sum_p += pj;
                    }
                    b_obs = sum_info - res.var_N * sum_p * sum_p;
                } else {
                    b_obs = sum_info;  // = 0
                }
            }
            // negH_i = -d^2 ell_i / db^2 = b_obs * Z_i Z_i^T (rank-1).
            // fisher_i = b_fis * Z_i Z_i^T (PSD Newton curvature).
            for (int c1 = 0; c1 < d; ++c1) {
                const double z1 = Zsite(i, c1);
                if (want_negH) {
                    const double v = b_obs * z1;
                    for (int c2 = 0; c2 < d; ++c2) e.negH(c1, c2) += v * Zsite(i, c2);
                }
                if (want_fisher) {
                    const double v = b_fis * z1;
                    for (int c2 = 0; c2 < d; ++c2) e.fisher(c1, c2) += v * Zsite(i, c2);
                }
            }
        }

        // ---- theta-space data score (length n_theta) ----
        // beta_lambda block: chain grad_eta_lambda through Xlam.row(i).
        // beta_p block: chain per-visit grad_eta_p[j] through Xp_site[i].row(j).
        // log_r block (NB): the per-site dispersion score grad_theta_i.
        if (want_theta_grad) {
            for (int c = 0; c < p_lam; ++c)
                e.theta_grad(c) += Xlam(i, c) * res.grad_eta_lambda;
            for (int j = 0; j < J; ++j)
                for (int c = 0; c < p_p; ++c)
                    e.theta_grad(p_lam + c) += Xp_site[i](j, c) * res.grad_eta_p[j];
            if (is_nb) e.theta_grad(p_lam + p_p) += res.grad_theta;
        }
    }
    return e;
}

void NMixGroupedOracle::grad_hess(int g, const double* b, double& logL,
                                  double* grad, double* negH) const {
    emit_grad_hess(eval_group(g, b, /*want_negH=*/true, /*want_fisher=*/false,
                              /*want_theta_grad=*/false),
                   d, logL, grad, negH);
}

void NMixGroupedOracle::node_ll(int g, const double* B, int n_nodes,
                                double* out) const {
    std::vector<double> eta_p_shift;
    for (int k = 0; k < n_nodes; ++k) {
        const double* bk = B + (std::size_t)k * d;
        double ll = 0.0;
        for (const int i : sites_by_group[g]) {
            double shift = 0.0;
            for (int c = 0; c < d; ++c) shift += Zsite(i, c) * bk[c];
            const int J = (int)Xp_site[i].rows();
            double eta_lam;
            const double* eta_p_ptr;
            if (arm == 0) {
                eta_lam = clamp30(eta_lambda_base(i) + shift);
                eta_p_ptr = eta_p_site[i].empty() ? nullptr : eta_p_site[i].data();
            } else {
                eta_lam = eta_lambda_base(i);
                const double shift_cl = clamp30(shift);
                eta_p_shift.assign(J, 0.0);
                for (int j = 0; j < J; ++j) eta_p_shift[j] = eta_p_site[i][j] + shift_cl;
                eta_p_ptr = eta_p_shift.empty() ? nullptr : eta_p_shift.data();
            }
            ll += compute_nmix_site_cached(site_cache[i], eta_p_ptr, eta_lam, r).log_lik;
        }
        out[k] = ll;
    }
}

void NMixGroupedOracle::theta_score(int g, const double* b,
                                    double* dl_dtheta) const {
    const GroupEval e = eval_group(g, b, /*want_negH=*/false,
                                   /*want_fisher=*/false,
                                   /*want_theta_grad=*/true);
    for (int i = 0; i < n_theta; ++i) dl_dtheta[i] = e.theta_grad(i);
}

bool NMixGroupedOracle::newton_hess(int g, const double* b, double* H) const {
    emit_fisher(eval_group(g, b, /*want_negH=*/false, /*want_fisher=*/true,
                           /*want_theta_grad=*/false), d, H);
    return true;
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
