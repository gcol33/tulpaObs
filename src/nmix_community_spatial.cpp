// nmix_community_spatial.cpp
// Nested Laplace-EM for the SPATIAL community / multispecies N-mixture
// (spAbundance sfMsNMix analogue): a per-species Royle (2004) N-mixture with
// Gaussian community hyperpriors on the per-species coefficients AND one shared
// areal field f on the abundance arm.
//
//   N_{s,i}        ~ Poisson(lambda_{s,i})            (or NB, global size r)
//   y_{s,i,j} | N  ~ Binomial(N_{s,i}, p_{s,i,j})
//   log lambda_{s,i} = X_lambda_i . (mu_lambda + b_lambda_s) + f_{u(i)}
//   logit p_{s,i,j}  = X_p_{ij}   . (mu_p      + b_p_s)
//   b_lambda_s ~ N(0, Sigma_lambda),  b_p_s ~ N(0, Sigma_p)
//   f ~ ICAR(tau) | BYM2(sigma,rho) | CAR(tau,rho) | SPDE (Matern)
//
// The Laplace-EM (field geometry, arrowhead Newton, EM covariance update,
// observed-info final pass, outer hyperparameter-grid driver) is the SAME
// algorithm the community occupancy family uses (ms_occu_spatial.cpp) and lives
// once in community_spatial_em.h; the only thing specific to this family is the
// per-(species, site) marginal cell below (`site_blocks`), which reads the
// Royle N-mixture kernel (nmix_kernel.h) via the pre-grouped, lgamma-cached
// oracle (NMixCommunityOracle) its constructor built.

#include "nmix_kernel.h"
#include "nmix_community_oracle.h"
#include "community_spatial_em.h"
#include "tobs_math.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using tulpaObs::clamp_eta;

namespace {

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;
using tulpaObs::NMixCommunityOracle;
using tulpaObs::NMixSiteResult;
using tulpaObs::compute_nmix_site_cached;
using tulpaObs::CommFieldKind;
using tulpaObs::CommSpdeCtx;
using tulpaObs::CommGridPoint;

// Per-site augmented assembly. na = d + nfield, with the d coefficient coords
// first and the nfield field loadings last (eta coord 0 = lambda, 1..J =
// visits). Fills grad_aug (na) and, when want_block, the symmetric small block
// (na x na): small = Z_aug' B Z_aug with B the per-site eta-space curvature
// (complete-data Fisher diagonal, plus -Var[N|y] rank-1 when want_obs). The
// d x d top-left of small IS the per-species coefficient curvature (the
// nmix_community_em block); the field rows/cols and crosses extend it to the
// shared field. Returns the per-site marginal log-lik. A no-visit (J = 0)
// species-site contributes nothing (interpolated by the field / prior), matching
// the single-species held-out convention.
double site_blocks(const NMixCommunityOracle::SiteRec& rec,
                   const Map<MatrixXd>& Xlam, int p_lam, int p_p,
                   const VectorXd& coef, double field_offset,
                   const double* load, int nfield, double r,
                   bool want_block, bool want_obs,
                   VectorXd& grad_aug, MatrixXd& small,
                   double* boundary_out = nullptr) {
    const int d  = p_lam + p_p;
    const int na = d + nfield;
    const int J  = rec.cache.n_visits;
    grad_aug.setZero(na);
    if (want_block) small.setZero(na, na);
    if (J == 0) return 0.0;

    double eta_lam = field_offset;
    for (int c = 0; c < p_lam; ++c) eta_lam += Xlam(rec.site, c) * coef(c);
    eta_lam = clamp_eta(eta_lam);
    std::vector<double> eta_p(J);
    for (int j = 0; j < J; ++j) {
        double v = 0.0;
        for (int c = 0; c < p_p; ++c) v += rec.Xp(j, c) * coef(p_lam + c);
        eta_p[j] = clamp_eta(v);
    }
    const NMixSiteResult res = compute_nmix_site_cached(rec.cache, eta_p.data(),
                                                        eta_lam, r);
    if (boundary_out) *boundary_out = res.boundary_weight;
    // Gradient: eta coord 0 -> mu_lambda cols (Xlam) and field cols (loadings);
    // eta coord 1+j -> mu_p cols (Xp).
    for (int c = 0; c < p_lam; ++c)
        grad_aug(c) += Xlam(rec.site, c) * res.grad_eta_lambda;
    for (int f = 0; f < nfield; ++f)
        grad_aug(d + f) += load[f] * res.grad_eta_lambda;
    for (int j = 0; j < J; ++j)
        for (int c = 0; c < p_p; ++c)
            grad_aug(p_lam + c) += rec.Xp(j, c) * res.grad_eta_p[j];

    if (!want_block) return res.log_lik;

    // Augmented design Z_aug ((1+J) x na): row 0 = lambda, rows 1+j = visits.
    MatrixXd Z = MatrixXd::Zero(1 + J, na);
    for (int c = 0; c < p_lam; ++c) Z(0, c) = Xlam(rec.site, c);
    for (int f = 0; f < nfield; ++f) Z(0, d + f) = load[f];
    for (int j = 0; j < J; ++j)
        for (int c = 0; c < p_p; ++c) Z(1 + j, p_lam + c) = rec.Xp(j, c);

    // Per-site eta-space block B ((1+J)^2): complete-data Fisher diagonal, plus
    // the -Var[N|y] rank-1 coupling (Louis 1982) when want_obs.
    MatrixXd B = MatrixXd::Zero(1 + J, 1 + J);
    B(0, 0) = res.info_eta_lambda;
    for (int j = 0; j < J; ++j) B(1 + j, 1 + j) = res.info_eta_p[j];
    if (want_obs && J > 0) {
        VectorXd vv(1 + J);
        vv(0) = -res.score_wt_lambda;
        for (int j = 0; j < J; ++j) {
            const double e = eta_p[j];
            vv(1 + j) = (e > 0.0) ? 1.0 / (1.0 + std::exp(-e))
                                  : std::exp(e) / (1.0 + std::exp(e));
        }
        B.noalias() -= res.var_N * (vv * vv.transpose());
    }
    small.noalias() = Z.transpose() * B * Z;
    return res.log_lik;
}

// Adapt site_blocks() to the shared driver's uniform SiteBlockFn signature: a
// (species, site) pair -> the pre-grouped SiteRec is orc.sp_sites[s][site]
// (sp_sites has exactly n_sites entries per species, index == site, built by
// the oracle constructor), so no lookup beyond direct indexing is needed.
auto make_nmix_site_block_fn(const NMixCommunityOracle& orc, const Map<MatrixXd>& Xl,
                             int p_lam, int p_p, double r) {
    return [&orc, &Xl, p_lam, p_p, r](
               int s, int site, const VectorXd& coef, double field_offset,
               const double* load, int nfield, bool want_block, bool want_obs,
               VectorXd& grad_aug, MatrixXd& small, double* boundary_out) -> double {
        return site_blocks(orc.sp_sites[s][site], Xl, p_lam, p_p, coef,
                           field_offset, load, nfield, r, want_block, want_obs,
                           grad_aug, small, boundary_out);
    };
}

// Resolve the oracle XPtr to a concrete NMixCommunityOracle (the pre-grouped,
// lgamma-cached data container its constructor built).
const NMixCommunityOracle& as_community_oracle(SEXP oracle) {
    Rcpp::XPtr<tulpa::REGroupOracle> base(oracle);
    NMixCommunityOracle* orcp = dynamic_cast<NMixCommunityOracle*>(base.get());
    if (orcp == nullptr)
        Rcpp::stop("oracle is not an NMixCommunityOracle.");
    return *orcp;
}

std::vector<int> resolve_map(const Rcpp::IntegerVector& map_R, int n_sites,
                             int n_spatial) {
    if ((int)map_R.size() != n_sites)
        Rcpp::stop("length(map_site_to_unit) must equal n_sites.");
    std::vector<int> map(n_sites);
    for (int s = 0; s < n_sites; ++s) {
        const int u = map_R[s] - 1;
        if (u < 0 || u >= n_spatial)
            Rcpp::stop("map_site_to_unit out of range [1, n_spatial].");
        map[s] = u;
    }
    return map;
}

// Shared entry-point plumbing: unpack the oracle / designs, dimension-check the
// warm start, and hand off to the shared areal driver with an nmix site-block
// factory (closes over r, the per-grid-point NB size / +Inf).
Rcpp::List run_nmix_spatial_grid(
    CommFieldKind kind, int field_len,
    SEXP oracle,
    const Rcpp::IntegerVector& map_site_to_unit_R,
    const Rcpp::NumericMatrix& X_lambda_R,
    int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const std::vector<CommGridPoint>& plan,
    const Rcpp::NumericMatrix& theta_grid_out,
    const Rcpp::NumericVector& mu_init,
    const Rcpp::NumericMatrix& Sigma_lambda_init,
    const Rcpp::NumericMatrix& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose,
    bool progress, int progress_every, double progress_throttle,
    const std::string& progress_file) {

    const NMixCommunityOracle& orc = as_community_oracle(oracle);
    const int n_sites = X_lambda_R.nrow();
    const int p_lam = orc.p_lam, p_p = orc.p_p, d = p_lam + p_p;
    Map<MatrixXd> Xl(REAL(X_lambda_R), n_sites, p_lam);
    std::vector<int> map = resolve_map(map_site_to_unit_R, n_sites, n_spatial);
    if ((int)mu_init.size() != d)
        Rcpp::stop("mu_init length must equal p_lam + p_p.");
    if (Sigma_lambda_init.nrow() != p_lam || Sigma_lambda_init.ncol() != p_lam)
        Rcpp::stop("Sigma_lambda_init must be p_lam x p_lam.");
    if (Sigma_p_init.nrow() != p_p || Sigma_p_init.ncol() != p_p)
        Rcpp::stop("Sigma_p_init must be p_p x p_p.");
    VectorXd mu0 = Map<VectorXd>(REAL(mu_init), d);
    MatrixXd Sl0 = Map<MatrixXd>(REAL(Sigma_lambda_init), p_lam, p_lam);
    MatrixXd Sp0 = Map<MatrixXd>(REAL(Sigma_p_init), p_p, p_p);

    return tulpaObs::run_community_spatial_grid(
        kind, field_len, orc.n_groups, n_sites, p_lam, p_p, map, n_spatial,
        adj_row_ptr, adj_col_idx, n_neighbors, plan, theta_grid_out,
        mu0, Sl0, Sp0, max_iter_em, tol_em, inner_max, inner_tol,
        verbose, "em", "Sigma_lambda", "b_lambda", "p_lambda",
        [&orc, &Xl, p_lam, p_p](double r) {
            return make_nmix_site_block_fn(orc, Xl, p_lam, p_p, r);
        },
        progress, progress_every, progress_throttle, progress_file,
        "nmix-spatial", /*report_boundary=*/true);
}

}  // namespace

// ---------------------------------------------------------------------------
// Entry points (one per field kind). `oracle` is the XPtr from
// cpp_nmix_community_oracle(); the outer grid integrates the field hyperparameter
// (tau / sigma, rho) and, under NB, the size r (outermost axis). Each grid point
// is a cold-restart Laplace-EM (the (lambda, p) identifiability ridge shifts as
// the field absorbs different variance, so warm-starting confounds the inner
// Newton -- matching the single-species spatial path).
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_spatial_icar(
    SEXP oracle,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector r_grid,
    Rcpp::NumericVector mu_init,
    Rcpp::NumericMatrix Sigma_lambda_init,
    Rcpp::NumericMatrix Sigma_p_init,
    int max_iter_em = 100, double tol_em = 1e-4,
    int inner_max = 50, double inner_tol = 1e-6,
    double sigma_beta = 100.0, bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = "") {
    (void) sigma_beta;

    const int n_tau = tau_grid.size(), n_r = r_grid.size();
    const int n_grid = n_tau * n_r;
    Rcpp::NumericMatrix theta_grid_out(n_grid, 2);
    std::vector<tulpaObs::CommGridPoint> plan(n_grid);
    int k = 0;
    for (int ri = 0; ri < n_r; ++ri)
        for (int t = 0; t < n_tau; ++t, ++k) {
            theta_grid_out(k, 0) = tau_grid[t];
            theta_grid_out(k, 1) = r_grid[ri];
            plan[k] = tulpaObs::CommGridPoint{ tau_grid[t], /*rho=*/1.0,
                                               /*log_det_Q_rho=*/0.0,
                                               /*a=*/0.0, /*b=*/0.0, r_grid[ri] };
        }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("tau", "r");
    return run_nmix_spatial_grid(
        tulpaObs::CommFieldKind::ICAR, /*field_len=*/n_spatial, oracle,
        map_site_to_unit_R, X_lambda_R, n_spatial, adj_row_ptr, adj_col_idx,
        n_neighbors, plan, theta_grid_out, mu_init, Sigma_lambda_init, Sigma_p_init,
        max_iter_em, tol_em, inner_max, inner_tol, verbose,
        progress, progress_every, progress_throttle, progress_file);
}

// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_spatial_car_proper(
    SEXP oracle,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector rho_grid,
    Rcpp::NumericVector log_det_Q_rho,   // one per rho grid point
    Rcpp::NumericVector r_grid,
    Rcpp::NumericVector mu_init,
    Rcpp::NumericMatrix Sigma_lambda_init,
    Rcpp::NumericMatrix Sigma_p_init,
    int max_iter_em = 100, double tol_em = 1e-4,
    int inner_max = 50, double inner_tol = 1e-6,
    double sigma_beta = 100.0, bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = "") {
    (void) sigma_beta;

    const int n_tau = tau_grid.size(), n_rho = rho_grid.size(), n_r = r_grid.size();
    if ((int)log_det_Q_rho.size() != n_rho)
        Rcpp::stop("length(log_det_Q_rho) must equal length(rho_grid).");
    const int n_grid = n_tau * n_rho * n_r;
    Rcpp::NumericMatrix theta_grid_out(n_grid, 3);
    std::vector<tulpaObs::CommGridPoint> plan(n_grid);
    int k = 0;
    for (int ri = 0; ri < n_r; ++ri)
        for (int rh = 0; rh < n_rho; ++rh)
            for (int t = 0; t < n_tau; ++t, ++k) {
                theta_grid_out(k, 0) = tau_grid[t];
                theta_grid_out(k, 1) = rho_grid[rh];
                theta_grid_out(k, 2) = r_grid[ri];
                plan[k] = tulpaObs::CommGridPoint{ tau_grid[t], rho_grid[rh],
                                                   log_det_Q_rho[rh],
                                                   /*a=*/0.0, /*b=*/0.0, r_grid[ri] };
            }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("tau", "rho", "r");
    return run_nmix_spatial_grid(
        tulpaObs::CommFieldKind::CAR_PROPER, /*field_len=*/n_spatial, oracle,
        map_site_to_unit_R, X_lambda_R, n_spatial, adj_row_ptr, adj_col_idx,
        n_neighbors, plan, theta_grid_out, mu_init, Sigma_lambda_init, Sigma_p_init,
        max_iter_em, tol_em, inner_max, inner_tol, verbose,
        progress, progress_every, progress_throttle, progress_file);
}

// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_spatial_bym2(
    SEXP oracle,
    Rcpp::IntegerVector map_site_to_unit_R,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::IntegerVector adj_row_ptr,
    Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors,
    int n_spatial,
    Rcpp::NumericVector sigma_grid,
    Rcpp::NumericVector rho_grid,
    double scale_factor,
    Rcpp::NumericVector r_grid,
    Rcpp::NumericVector mu_init,
    Rcpp::NumericMatrix Sigma_lambda_init,
    Rcpp::NumericMatrix Sigma_p_init,
    int max_iter_em = 100, double tol_em = 1e-4,
    int inner_max = 50, double inner_tol = 1e-6,
    double sigma_beta = 100.0, bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = "") {
    (void) sigma_beta;

    const int n_sig = sigma_grid.size(), n_rho = rho_grid.size(), n_r = r_grid.size();
    const int n_grid = n_sig * n_rho * n_r;
    Rcpp::NumericMatrix theta_grid_out(n_grid, 3);
    std::vector<tulpaObs::CommGridPoint> plan(n_grid);
    int k = 0;
    for (int ri = 0; ri < n_r; ++ri)
        for (int rh = 0; rh < n_rho; ++rh)
            for (int sg = 0; sg < n_sig; ++sg, ++k) {
                const double sigma = sigma_grid[sg], rho = rho_grid[rh];
                const double a = sigma * std::sqrt(rho / scale_factor);
                const double b = sigma * std::sqrt(1.0 - rho);
                theta_grid_out(k, 0) = sigma;
                theta_grid_out(k, 1) = rho;
                theta_grid_out(k, 2) = r_grid[ri];
                plan[k] = tulpaObs::CommGridPoint{ /*tau=*/1.0, rho,
                                                   /*log_det_Q_rho=*/0.0,
                                                   a, b, r_grid[ri] };
            }
    Rcpp::colnames(theta_grid_out) = Rcpp::CharacterVector::create("sigma", "rho", "r");
    return run_nmix_spatial_grid(
        tulpaObs::CommFieldKind::BYM2, /*field_len=*/2 * n_spatial, oracle,
        map_site_to_unit_R, X_lambda_R, n_spatial, adj_row_ptr, adj_col_idx,
        n_neighbors, plan, theta_grid_out, mu_init, Sigma_lambda_init, Sigma_p_init,
        max_iter_em, tol_em, inner_max, inner_tol, verbose,
        progress, progress_every, progress_throttle, progress_file);
}

// Continuous Matern (SPDE) shared field on the abundance arm. The field lives
// at n_mesh FEM nodes; the dense projection A (n_sites x n_mesh) maps mesh nodes
// onto sites, shared across species exactly as the areal field is. Per grid
// point the proper Matern precision Q(range, sigma) (carrying tau_spde^2) and
// its log|Q| are built once on the R side and passed in; the outer grid axes are
// (range, sigma) [, NB size r]. n_spatial is the mesh node count here.
// [[Rcpp::export]]
Rcpp::List cpp_nmix_community_spatial_spde(
    SEXP oracle,
    Rcpp::NumericMatrix X_lambda_R,
    Rcpp::NumericMatrix A_R,            // dense [n_sites x n_mesh]
    Rcpp::List Q_list,                  // per-grid-point precision [n_mesh x n_mesh]
    Rcpp::NumericVector log_det_Q,      // per-grid-point log|Q|
    Rcpp::NumericMatrix theta_grid_R,   // [n_grid x n_theta] (range, sigma, r)
    Rcpp::NumericVector r_grid,         // NB size per grid point (or +Inf)
    Rcpp::NumericVector mu_init,
    Rcpp::NumericMatrix Sigma_lambda_init,
    Rcpp::NumericMatrix Sigma_p_init,
    int max_iter_em = 100, double tol_em = 1e-4,
    int inner_max = 50, double inner_tol = 1e-6,
    double sigma_beta = 100.0, bool verbose = false,
    bool progress = false, int progress_every = 0,
    double progress_throttle = 0.0, std::string progress_file = "") {
    (void) sigma_beta;

    const NMixCommunityOracle& orc = as_community_oracle(oracle);
    const int n_sites = X_lambda_R.nrow();
    const int n_mesh  = A_R.ncol();
    const int p_lam = orc.p_lam, p_p = orc.p_p, d = p_lam + p_p;
    if (A_R.nrow() != n_sites)
        Rcpp::stop("nrow(A) must equal nrow(X_lambda).");
    Map<MatrixXd> Xl(REAL(X_lambda_R), n_sites, p_lam);
    const MatrixXd A = Map<MatrixXd>(REAL(A_R), n_sites, n_mesh);
    if ((int)mu_init.size() != d)
        Rcpp::stop("mu_init length must equal p_lam + p_p.");
    if (Sigma_lambda_init.nrow() != p_lam || Sigma_lambda_init.ncol() != p_lam)
        Rcpp::stop("Sigma_lambda_init must be p_lam x p_lam.");
    if (Sigma_p_init.nrow() != p_p || Sigma_p_init.ncol() != p_p)
        Rcpp::stop("Sigma_p_init must be p_p x p_p.");
    VectorXd mu0 = Map<VectorXd>(REAL(mu_init), d);
    MatrixXd Sl0 = Map<MatrixXd>(REAL(Sigma_lambda_init), p_lam, p_lam);
    MatrixXd Sp0 = Map<MatrixXd>(REAL(Sigma_p_init), p_p, p_p);

    return tulpaObs::run_community_spatial_grid_spde(
        orc.n_groups, A, n_mesh, p_lam, p_p, Q_list, log_det_Q, theta_grid_R,
        r_grid, mu0, Sl0, Sp0, max_iter_em, tol_em, inner_max, inner_tol,
        verbose, "em", "Sigma_lambda", "b_lambda", "p_lambda",
        [&orc, &Xl, p_lam, p_p](double r) {
            return make_nmix_site_block_fn(orc, Xl, p_lam, p_p, r);
        },
        progress, progress_every, progress_throttle, progress_file,
        "nmix-spatial", /*report_boundary=*/true);
}
