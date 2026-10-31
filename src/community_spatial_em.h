// community_spatial_em.h
// Shared field-structure layer + Laplace-EM driver for a per-species marginal
// with Gaussian community hyperpriors on the per-species coefficients AND one
// shared ICAR / BYM2 / proper-CAR / SPDE field on the STATE arm (the sfMsNMix /
// svcMsPGOcc pattern). Two response families instantiate this: the community
// Royle N-mixture (nmix_community_spatial.cpp) and the community two-state
// occupancy marginal (ms_occu_spatial.cpp). The algorithm is identical between
// them -- a Laplace-EM over species random effects with a shared areal field,
// solved by a Schur-complement Newton step with step-halving, then an EM update
// of the community covariances, wrapped in an outer hyperparameter grid; only
// the per-species marginal (Royle N-mixture vs occupancy two-state) differs.
//
// The family-specific piece is ONE callable, passed as a CommSiteBlockRef:
//
//   double site_block_fn(int s, int site, const VectorXd& coef,
//                        double field_offset, const double* load, int nfield,
//                        bool want_block, bool want_obs,
//                        VectorXd& grad_aug, MatrixXd& small,
//                        double* boundary_out)
//
// evaluating one (species, site) cell: `coef` is that species' full
// (state, detection) coefficient vector, `field_offset` the shared field's
// current contribution to eta_state at that site, `load`/`nfield` the field
// loading geometry (see FieldGeom below). It fills `grad_aug` (length
// d + nfield: coefficient columns then field columns) and, when `want_block`,
// the symmetric curvature `small` ((d + nfield) x (d + nfield)) via the
// design-sandwiched per-site eta-space block; returns the site's marginal
// log-lik. Each family builds this from its own per-site kernel
// (nmix_kernel.h's compute_nmix_site_cached, ms_occu_kernel.h's
// ms_occu_site_cell) and adapts it to this signature with a small lambda that
// captures its own data container + design; that adapter is the ONLY
// family-specific code left outside this header.
//
// Everything else -- field geometry / offset / log-prior / prior-gradient /
// sum-to-zero centering (areal + SPDE), the arrowhead block-elimination Newton,
// the EM covariance update, the observed-info final pass + community-mean
// covariance, and the outer hyperparameter-grid driver (areal and SPDE) -- is
// family-agnostic and is compiled once, in community_spatial_em.cpp.

#ifndef TULPAOBS_COMMUNITY_SPATIAL_EM_H
#define TULPAOBS_COMMUNITY_SPATIAL_EM_H

#include "function_ref.h"
#include "nmix_progress.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <string>
#include <vector>

namespace tulpaObs {

using Eigen::MatrixXd;
using Eigen::VectorXd;

enum class CommFieldKind { ICAR, CAR_PROPER, BYM2, SPDE };

// SPDE field context: the continuous Matern field lives at n_mesh FEM nodes;
// the dense projection A (n_sites x n_mesh) maps mesh nodes onto sites, and Q
// (n_mesh x n_mesh, carrying tau_spde^2) is the proper Matern precision built
// once per grid point on the R side. The areal kinds carry a single field unit
// per site; SPDE carries the whole mesh row, so the field loading and prior are
// threaded through this struct rather than the (u, a, b) areal triple.
struct CommSpdeCtx {
    const MatrixXd* A = nullptr;     // n_sites x n_mesh
    const MatrixXd* Q = nullptr;     // n_mesh x n_mesh (with tau_spde^2)
    double log_det_Q = 0.0;
    int n_mesh = 0;
};

struct CommSpatialResult {
    VectorXd mu;             // d community means
    VectorXd field;          // n_field_total = nfield * n_spatial (ICAR/CAR/SPDE: f; BYM2: [v; w])
    MatrixXd Sigma_state, Sigma_p;
    MatrixXd blup_state, blup_p;   // S x p_state, S x p_p
    MatrixXd vcov_mu;              // d x d community-mean covariance (b- and field-folded)
    std::vector<MatrixXd> b_cov;   // S of d x d: Cov(b_s | y, mu, field) at the mode
    double log_marginal = R_NegInf;
    double log_lik = R_NegInf;
    bool   converged = false;
    int    n_iter = 0;
    double boundary_max = 0.0;
};

// The per-(species, site) evaluator the driver calls; the header comment above
// describes each argument.
using CommSiteBlockSig = double(int s, int site, const VectorXd& coef,
                                double field_offset, const double* load, int nfield,
                                bool want_block, bool want_obs,
                                VectorXd& grad_aug, MatrixXd& small,
                                double* boundary_out);
using CommSiteBlockRef = FunctionRef<CommSiteBlockSig>;

// One outer-grid-point Laplace-EM fit, over the marginal `site_block_fn`
// evaluates. `field` is laid out [f] (ICAR/CAR/SPDE) or [v; w] (BYM2); `a`, `b`
// are the BYM2 loadings (unused otherwise); `map_site_to_unit` is unused (may be
// empty) for SPDE, which reads its geometry from `spde` instead.
CommSpatialResult community_spatial_em(
    CommFieldKind kind,
    int S, int n_sites, int p_state, int p_p,
    const std::vector<int>& map_site_to_unit,
    int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    double tau, double rho, double log_det_Q_rho, double a, double b,
    const VectorXd& mu_init, const MatrixXd& Sigma_state_init,
    const MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose, const char* verbose_label,
    CommSiteBlockRef site_block_fn,
    const CommSpdeCtx* spde = nullptr);

// Pack the per-grid-point results into the list the R side reads.
Rcpp::List community_spatial_pack(int d, int field_len, int n_grid,
                                  const std::vector<CommSpatialResult>& results,
                                  const Rcpp::NumericMatrix& theta_grid_out,
                                  const char* sigma_state_name,
                                  const char* blup_state_name,
                                  const char* p_state_name,
                                  int p_state, int p_p, int n_spatial,
                                  bool report_boundary);

// One outer-grid point: the field hyperparameters (tau / rho / log|Q(rho)| and
// the BYM2 mixing loadings a, b) plus, for count families, a per-point nuisance
// (NB size r) the inner Laplace-EM conditions on. `r` is ignored by families
// with no such nuisance (occupancy passes a placeholder). Each entry builds a
// plan of these in its own grid-nesting order; the driver below just walks it.
struct CommGridPoint { double tau, rho, log_det_Q_rho, a, b, r; };

// Shared outer-grid driver for an AREAL community spatial nested-Laplace-EM
// (icar / car_proper / bym2). `make_site_block_fn(r)` builds the family-
// specific per-site evaluator closure for one grid point's nuisance `r`
// (ignored by families with none). Progress reporting is opt-in per call
// (occupancy has none today; count families report under "nmix-spatial").
template <class SiteBlockFactory>
Rcpp::List run_community_spatial_grid(
    CommFieldKind kind, int field_len,
    int S, int n_sites, int p_state, int p_p,
    const std::vector<int>& map_site_to_unit, int n_spatial,
    const Rcpp::IntegerVector& adj_row_ptr,
    const Rcpp::IntegerVector& adj_col_idx,
    const Rcpp::IntegerVector& n_neighbors,
    const std::vector<CommGridPoint>& plan,
    const Rcpp::NumericMatrix& theta_grid_out,
    const VectorXd& mu_init, const MatrixXd& Sigma_state_init,
    const MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose, const char* verbose_label,
    const char* sigma_state_name, const char* blup_state_name,
    const char* p_state_name,
    SiteBlockFactory make_site_block_fn,
    bool progress, int progress_every, double progress_throttle,
    const std::string& progress_file, const char* progress_grid_label,
    bool report_boundary) {

    const int n_grid = (int)plan.size();
    std::vector<CommSpatialResult> results(n_grid);
    auto gp = tulpaObs::make_grid_progress(progress_grid_label, n_grid, progress,
                                           progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        const CommGridPoint& g = plan[k];
        auto site_block_fn = make_site_block_fn(g.r);
        results[k] = community_spatial_em(
            kind, S, n_sites, p_state, p_p, map_site_to_unit, n_spatial,
            adj_row_ptr, adj_col_idx, n_neighbors,
            g.tau, g.rho, g.log_det_Q_rho, g.a, g.b,
            mu_init, Sigma_state_init, Sigma_p_init,
            max_iter_em, tol_em, inner_max, inner_tol, verbose, verbose_label,
            CommSiteBlockRef(site_block_fn));
        if (gp) gp->tick();
    }
    if (gp) gp->finish();
    return community_spatial_pack(
        p_state + p_p, field_len, n_grid, results, theta_grid_out,
        sigma_state_name, blup_state_name, p_state_name, p_state, p_p, n_spatial,
        report_boundary);
}

// Shared outer-grid driver for the continuous Matern (SPDE) shared field. The
// field lives at n_mesh FEM nodes; the dense projection A (n_sites x n_mesh)
// maps mesh nodes onto sites, shared across species exactly as the areal field
// is. Per grid point the proper Matern precision Q(range, sigma) (carrying
// tau_spde^2) and its log|Q| are built once on the R side and passed in via
// `Q_list` / `log_det_Q`; the outer grid axes are (range, sigma) [, a per-point
// nuisance `r_grid`].
template <class SiteBlockFactory>
Rcpp::List run_community_spatial_grid_spde(
    int S, const MatrixXd& A, int n_mesh, int p_state, int p_p,
    const Rcpp::List& Q_list, const Rcpp::NumericVector& log_det_Q,
    const Rcpp::NumericMatrix& theta_grid_R, const Rcpp::NumericVector& r_grid,
    const VectorXd& mu_init, const MatrixXd& Sigma_state_init,
    const MatrixXd& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose, const char* verbose_label,
    const char* sigma_state_name, const char* blup_state_name,
    const char* p_state_name,
    SiteBlockFactory make_site_block_fn,
    bool progress, int progress_every, double progress_throttle,
    const std::string& progress_file, const char* progress_grid_label,
    bool report_boundary) {

    const int n_grid = Q_list.size();
    if ((int)log_det_Q.size() != n_grid)
        Rcpp::stop("length(log_det_Q) must equal length(Q_list).");
    if (theta_grid_R.nrow() != n_grid)
        Rcpp::stop("nrow(theta_grid) must equal length(Q_list).");
    if ((int)r_grid.size() != n_grid)
        Rcpp::stop("length(r_grid) must equal length(Q_list).");

    // Empty CSR adjacency (unused on the SPDE path).
    Rcpp::IntegerVector empty_int(0);
    const int n_sites = (int)A.rows();

    std::vector<CommSpatialResult> results(n_grid);
    std::vector<MatrixXd> Qmats(n_grid);   // keep alive for the CommSpdeCtx pointers
    auto gp = tulpaObs::make_grid_progress(progress_grid_label, n_grid, progress,
                                           progress_every, progress_throttle, progress_file);
    for (int k = 0; k < n_grid; ++k) {
        Rcpp::NumericMatrix Qk_R = Q_list[k];
        if (Qk_R.nrow() != n_mesh || Qk_R.ncol() != n_mesh)
            Rcpp::stop("Q_list[[%d]] must be n_mesh x n_mesh.", k + 1);
        Qmats[k] = Eigen::Map<MatrixXd>(REAL(Qk_R), n_mesh, n_mesh);
        CommSpdeCtx sp;
        sp.A = &A; sp.Q = &Qmats[k]; sp.log_det_Q = log_det_Q[k]; sp.n_mesh = n_mesh;
        auto site_block_fn = make_site_block_fn(r_grid[k]);
        results[k] = community_spatial_em(
            CommFieldKind::SPDE, S, n_sites, p_state, p_p,
            /*map=*/std::vector<int>(), /*n_spatial=*/n_mesh,
            empty_int, empty_int, empty_int,
            /*tau=*/1.0, /*rho=*/1.0, /*log_det_Q_rho=*/0.0, /*a=*/0.0, /*b=*/0.0,
            mu_init, Sigma_state_init, Sigma_p_init,
            max_iter_em, tol_em, inner_max, inner_tol, verbose, verbose_label,
            CommSiteBlockRef(site_block_fn), &sp);
        if (gp) gp->tick();
    }
    if (gp) gp->finish();
    return community_spatial_pack(
        p_state + p_p, n_mesh, n_grid, results, theta_grid_R,
        sigma_state_name, blup_state_name, p_state_name, p_state, p_p, n_mesh,
        report_boundary);
}

}  // namespace tulpaObs

#endif  // TULPAOBS_COMMUNITY_SPATIAL_EM_H
