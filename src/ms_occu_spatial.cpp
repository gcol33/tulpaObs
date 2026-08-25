// ms_occu_spatial.cpp
// Areal-spatial community single-season occupancy (ms_occu() + shared field;
// the occupancy analogue of sfMsNMix). A per-species two-state occupancy
// model with Gaussian community hyperpriors on the per-species coefficients
// AND one shared ICAR / BYM2 / proper-CAR / SPDE field on the OCCUPANCY arm:
//
//   logit psi_{s,i} = X_psi_i . (mu_psi + b_psi_s) + f_{u(i)}
//   logit p_{s,i}   = X_p_i   . (mu_p   + b_p_s)
//   b_psi_s ~ N(0, Sigma_psi),  b_p_s ~ N(0, Sigma_p),  one f shared
//
// The latent z integrates out per species-site in closed form (the occupancy
// two-state marginal, ms_occu_kernel.h, site-level detection). The Laplace-EM
// (field geometry, arrowhead Newton, EM covariance update, observed-info final
// pass, outer hyperparameter-grid driver) is the SAME algorithm the community
// N-mixture family uses (nmix_community_spatial.cpp) and lives once in
// community_spatial_em.h; the only thing specific to this family is the
// per-(species, site) marginal cell below (`occ_site_blocks`), which comes from
// ms_occu_site_cell.

#include "tobs_math.h"
#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include "ms_occu_kernel.h"
#include "community_spatial_em.h"

using namespace Rcpp;
using tulpaObs::clamp_eta;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

namespace tulpaObs {
namespace {

// Per-(species) site records: the site index, occupancy / detection design rows
// are shared (site-level) across species, so only the summary (n_valid, n_det,
// any_det) is per species.
struct OccCommData {
    int n_sites = 0, n_species = 0, p_psi = 0, p_p = 0;
    Map<MatrixXd> X_psi;     // n_sites x p_psi
    Map<MatrixXd> X_p;       // n_sites x p_p
    // per species: n_valid / n_det / any_det over sites (length n_sites)
    std::vector<std::vector<int>> n_valid, n_det;
    std::vector<std::vector<char>> any_det;
    std::vector<int> map_site_to_unit;   // length n_sites (0-based field unit)

    OccCommData(const NumericMatrix& Xpsi_R, const NumericMatrix& Xp_R)
        : X_psi(REAL(Xpsi_R), Xpsi_R.nrow(), Xpsi_R.ncol()),
          X_p(REAL(Xp_R), Xp_R.nrow(), Xp_R.ncol()) {}
};

// One site's contribution: writes grad_aug (length na = d + nfield) and, when
// want_block, the curvature `small` (na x na) via the Z' B Z sandwich. The
// occupancy arm (eta_psi) carries the field loadings + X_psi; the detection arm
// (eta_p) carries X_p. Mirrors site_blocks in nmix_community_spatial.cpp but
// over the occupancy cell. Returns the per-site marginal log-lik.
inline double occ_site_blocks(const OccCommData& d_, int s, int site,
                              const VectorXd& coef, double field_offset,
                              const double* load, int nfield,
                              bool want_block, bool want_obs,
                              VectorXd& grad_aug, MatrixXd& small) {
    const int p_psi = d_.p_psi, p_p = d_.p_p, d = p_psi + p_p;
    const int na = d + nfield;
    grad_aug.setZero(na);
    if (want_block) small.setZero(na, na);
    const int nv = d_.n_valid[s][site];
    if (nv == 0) return 0.0;

    double eta_psi = field_offset;
    for (int c = 0; c < p_psi; ++c) eta_psi += d_.X_psi(site, c) * coef(c);
    eta_psi = clamp_eta(eta_psi);
    double eta_p = 0.0;
    for (int c = 0; c < p_p; ++c) eta_p += d_.X_p(site, c) * coef(p_psi + c);
    eta_p = clamp_eta(eta_p);

    const MsOccuSiteCell cell = ms_occu_site_cell(
        eta_psi, eta_p, nv, d_.n_det[s][site], d_.any_det[s][site] != 0, want_obs);

    // gradient: eta_psi -> mu_psi cols (X_psi) + field cols (loadings);
    //           eta_p   -> mu_p cols (X_p).
    for (int c = 0; c < p_psi; ++c) grad_aug(c) += d_.X_psi(site, c) * cell.g_psi;
    for (int f = 0; f < nfield; ++f) grad_aug(d + f) += load[f] * cell.g_psi;
    for (int c = 0; c < p_p; ++c) grad_aug(p_psi + c) += d_.X_p(site, c) * cell.g_p;

    if (!want_block) return cell.log_lik;

    // Augmented design Z (2 x na): row 0 = psi (mu_psi + field), row 1 = p.
    MatrixXd Z = MatrixXd::Zero(2, na);
    for (int c = 0; c < p_psi; ++c) Z(0, c) = d_.X_psi(site, c);
    for (int f = 0; f < nfield; ++f) Z(0, d + f) = load[f];
    for (int c = 0; c < p_p; ++c) Z(1, p_psi + c) = d_.X_p(site, c);

    // Per-site eta-space 2x2 negative-Hessian B.
    MatrixXd B(2, 2);
    B(0, 0) = cell.B_pp_psi;
    B(1, 1) = cell.B_p_p;
    B(0, 1) = cell.B_cross;
    B(1, 0) = cell.B_cross;
    small.noalias() = Z.transpose() * B * Z;
    return cell.log_lik;
}

// Adapt occ_site_blocks() to the shared driver's uniform SiteBlockFn signature.
// `r` is unused (occupancy has no per-grid-point nuisance).
auto make_occ_site_block_fn(const OccCommData& d_) {
    return [&d_](int s, int site, const VectorXd& coef, double field_offset,
                const double* load, int nfield, bool want_block, bool want_obs,
                VectorXd& grad_aug, MatrixXd& small, double* boundary_out) -> double {
        (void) boundary_out;
        return occ_site_blocks(d_, s, site, coef, field_offset, load, nfield,
                               want_block, want_obs, grad_aug, small);
    };
}

// Build the occupancy community data container from R inputs.
OccCommData build_occ_data(const NumericMatrix& X_psi_R, const NumericMatrix& X_p_R,
                           const IntegerMatrix& n_valid_R, const IntegerMatrix& n_det_R,
                           const IntegerVector& map_R, int n_spatial) {
    OccCommData d_(X_psi_R, X_p_R);
    d_.n_sites   = X_psi_R.nrow();
    d_.p_psi     = X_psi_R.ncol();
    d_.p_p       = X_p_R.ncol();
    d_.n_species = n_valid_R.ncol();
    if (n_valid_R.nrow() != d_.n_sites || n_det_R.nrow() != d_.n_sites ||
        n_det_R.ncol() != d_.n_species)
        Rcpp::stop("n_valid / n_det must be n_sites x n_species.");
    d_.n_valid.assign(d_.n_species, std::vector<int>(d_.n_sites, 0));
    d_.n_det.assign(d_.n_species, std::vector<int>(d_.n_sites, 0));
    d_.any_det.assign(d_.n_species, std::vector<char>(d_.n_sites, 0));
    for (int s = 0; s < d_.n_species; ++s)
        for (int i = 0; i < d_.n_sites; ++i) {
            d_.n_valid[s][i] = n_valid_R(i, s);
            d_.n_det[s][i]   = n_det_R(i, s);
            d_.any_det[s][i] = (n_det_R(i, s) > 0) ? 1 : 0;
        }
    if ((int)map_R.size() != d_.n_sites)
        Rcpp::stop("length(map_site_to_unit) must equal n_sites.");
    d_.map_site_to_unit.assign(d_.n_sites, 0);
    for (int i = 0; i < d_.n_sites; ++i) {
        const int u = map_R[i] - 1;
        if (u < 0 || u >= n_spatial)
            Rcpp::stop("map_site_to_unit out of range [1, n_spatial].");
        d_.map_site_to_unit[i] = u;
    }
    return d_;
}

Rcpp::List run_occ_spatial_grid(
    CommFieldKind kind, int field_len,
    const NumericMatrix& X_psi_R, const NumericMatrix& X_p_R,
    const IntegerMatrix& n_valid_R, const IntegerMatrix& n_det_R,
    const IntegerVector& map_R, int n_spatial,
    const IntegerVector& adj_row_ptr, const IntegerVector& adj_col_idx,
    const IntegerVector& n_neighbors,
    const std::vector<CommGridPoint>& plan, const NumericMatrix& theta_grid_out,
    const NumericVector& mu_init, const NumericMatrix& Sigma_psi_init,
    const NumericMatrix& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose) {

    OccCommData d_ = build_occ_data(X_psi_R, X_p_R, n_valid_R, n_det_R, map_R, n_spatial);
    const int p_psi = d_.p_psi, p_p = d_.p_p, d = p_psi + p_p;
    if ((int)mu_init.size() != d)
        Rcpp::stop("mu_init length must equal p_psi + p_p.");
    if (Sigma_psi_init.nrow() != p_psi || Sigma_psi_init.ncol() != p_psi)
        Rcpp::stop("Sigma_psi_init must be p_psi x p_psi.");
    if (Sigma_p_init.nrow() != p_p || Sigma_p_init.ncol() != p_p)
        Rcpp::stop("Sigma_p_init must be p_p x p_p.");
    VectorXd mu0 = Map<VectorXd>(REAL(mu_init), d);
    MatrixXd Spsi0 = Map<MatrixXd>(REAL(Sigma_psi_init), p_psi, p_psi);
    MatrixXd Sp0   = Map<MatrixXd>(REAL(Sigma_p_init), p_p, p_p);

    return tulpaObs::run_community_spatial_grid(
        kind, field_len, d_.n_species, d_.n_sites, p_psi, p_p,
        d_.map_site_to_unit, n_spatial,
        adj_row_ptr, adj_col_idx, n_neighbors, plan, theta_grid_out,
        mu0, Spsi0, Sp0, max_iter_em, tol_em, inner_max, inner_tol,
        verbose, "occ-em", "Sigma_psi", "b_psi", "p_psi",
        [&d_](double /*r*/) { return make_occ_site_block_fn(d_); },
        /*progress=*/false, 0, 0.0, "", "ms-occu-spatial",
        /*report_boundary=*/false);
}

// SPDE entry driver: builds OccCommData (no map_site_to_unit -- SPDE reads its
// geometry from the A projection) and runs the shared SPDE grid driver.
Rcpp::List run_occ_spatial_grid_spde(
    const NumericMatrix& X_psi_R, const NumericMatrix& X_p_R,
    const IntegerMatrix& n_valid_R, const IntegerMatrix& n_det_R,
    const MatrixXd& A, int n_mesh,
    const Rcpp::List& Q_list, const Rcpp::NumericVector& log_det_Q,
    const Rcpp::NumericMatrix& theta_grid_R,
    const NumericVector& mu_init, const NumericMatrix& Sigma_psi_init,
    const NumericMatrix& Sigma_p_init,
    int max_iter_em, double tol_em, int inner_max, double inner_tol,
    bool verbose) {

    OccCommData d_(X_psi_R, X_p_R);
    d_.n_sites   = X_psi_R.nrow();
    d_.p_psi     = X_psi_R.ncol();
    d_.p_p       = X_p_R.ncol();
    d_.n_species = n_valid_R.ncol();
    if (n_valid_R.nrow() != d_.n_sites || n_det_R.nrow() != d_.n_sites ||
        n_det_R.ncol() != d_.n_species)
        Rcpp::stop("n_valid / n_det must be n_sites x n_species.");
    d_.n_valid.assign(d_.n_species, std::vector<int>(d_.n_sites, 0));
    d_.n_det.assign(d_.n_species, std::vector<int>(d_.n_sites, 0));
    d_.any_det.assign(d_.n_species, std::vector<char>(d_.n_sites, 0));
    for (int s = 0; s < d_.n_species; ++s)
        for (int i = 0; i < d_.n_sites; ++i) {
            d_.n_valid[s][i] = n_valid_R(i, s);
            d_.n_det[s][i]   = n_det_R(i, s);
            d_.any_det[s][i] = (n_det_R(i, s) > 0) ? 1 : 0;
        }
    if (A.rows() != d_.n_sites)
        Rcpp::stop("nrow(A) must equal nrow(X_psi).");

    const int p_psi = d_.p_psi, p_p = d_.p_p, d = p_psi + p_p;
    if ((int)mu_init.size() != d)
        Rcpp::stop("mu_init length must equal p_psi + p_p.");
    if (Sigma_psi_init.nrow() != p_psi || Sigma_psi_init.ncol() != p_psi)
        Rcpp::stop("Sigma_psi_init must be p_psi x p_psi.");
    if (Sigma_p_init.nrow() != p_p || Sigma_p_init.ncol() != p_p)
        Rcpp::stop("Sigma_p_init must be p_p x p_p.");
    VectorXd mu0 = Map<VectorXd>(REAL(mu_init), d);
    MatrixXd Spsi0 = Map<MatrixXd>(REAL(Sigma_psi_init), p_psi, p_psi);
    MatrixXd Sp0   = Map<MatrixXd>(REAL(Sigma_p_init), p_p, p_p);
    Rcpp::NumericVector r_grid(Q_list.size(), R_PosInf);   // no dispersion nuisance

    return tulpaObs::run_community_spatial_grid_spde(
        d_.n_species, A, n_mesh, p_psi, p_p, Q_list, log_det_Q, theta_grid_R,
        r_grid, mu0, Spsi0, Sp0, max_iter_em, tol_em, inner_max, inner_tol,
        verbose, "occ-em", "Sigma_psi", "b_psi", "p_psi",
        [&d_](double /*r*/) { return make_occ_site_block_fn(d_); },
        /*progress=*/false, 0, 0.0, "", "ms-occu-spatial",
        /*report_boundary=*/false);
}

}  // namespace
}  // namespace tulpaObs


// ---------------------------------------------------------------------------
// Cross-check entry: per-site occupancy cell score + curvature vs the R oracle
// (FD-validated). Returns log-lik, score (g_psi, g_p), and the 2x2 -Hessian.
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_site_cell(double eta_psi, double eta_p, int n_valid,
                                 int n_det, bool observed) {
    const tulpaObs::MsOccuSiteCell c = tulpaObs::ms_occu_site_cell(
        eta_psi, eta_p, n_valid, n_det, n_det > 0, observed);
    Rcpp::NumericMatrix B(2, 2);
    B(0, 0) = c.B_pp_psi; B(1, 1) = c.B_p_p;
    B(0, 1) = c.B_cross;  B(1, 0) = c.B_cross;
    return Rcpp::List::create(
        Rcpp::Named("log_lik") = c.log_lik,
        Rcpp::Named("grad") = Rcpp::NumericVector::create(c.g_psi, c.g_p),
        Rcpp::Named("neg_hess") = B);
}

// ICAR
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_spatial_icar(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p,
    Rcpp::IntegerMatrix n_valid, Rcpp::IntegerMatrix n_det,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector mu_init, Rcpp::NumericMatrix Sigma_psi_init,
    Rcpp::NumericMatrix Sigma_p_init, int max_iter_em, bool verbose) {
    std::vector<tulpaObs::CommGridPoint> plan;
    Rcpp::NumericMatrix tg(tau_grid.size(), 1);
    for (int i = 0; i < tau_grid.size(); ++i) {
        plan.push_back({tau_grid[i], 1.0, 0.0, 0.0, 0.0, 0.0});
        tg(i, 0) = tau_grid[i];
    }
    colnames(tg) = Rcpp::CharacterVector::create("tau");
    return tulpaObs::run_occ_spatial_grid(
        tulpaObs::CommFieldKind::ICAR, n_spatial, X_psi, X_p, n_valid, n_det,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, mu_init, Sigma_psi_init, Sigma_p_init,
        max_iter_em, 1e-4, 50, 1e-7, verbose);
}

// Proper CAR
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_spatial_car_proper(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p,
    Rcpp::IntegerMatrix n_valid, Rcpp::IntegerMatrix n_det,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector tau_grid,
    Rcpp::NumericVector rho_grid, Rcpp::NumericVector log_det_Q_rho,
    Rcpp::NumericVector mu_init, Rcpp::NumericMatrix Sigma_psi_init,
    Rcpp::NumericMatrix Sigma_p_init, int max_iter_em, bool verbose) {
    std::vector<tulpaObs::CommGridPoint> plan;
    std::vector<double> tg_tau, tg_rho;
    for (int ri = 0; ri < rho_grid.size(); ++ri)
        for (int ti = 0; ti < tau_grid.size(); ++ti) {
            plan.push_back({tau_grid[ti], rho_grid[ri], log_det_Q_rho[ri], 0.0, 0.0, 0.0});
            tg_tau.push_back(tau_grid[ti]); tg_rho.push_back(rho_grid[ri]);
        }
    Rcpp::NumericMatrix tg(plan.size(), 2);
    for (std::size_t i = 0; i < plan.size(); ++i) { tg(i, 0) = tg_tau[i]; tg(i, 1) = tg_rho[i]; }
    colnames(tg) = Rcpp::CharacterVector::create("tau", "rho");
    return tulpaObs::run_occ_spatial_grid(
        tulpaObs::CommFieldKind::CAR_PROPER, n_spatial, X_psi, X_p, n_valid, n_det,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, mu_init, Sigma_psi_init, Sigma_p_init,
        max_iter_em, 1e-4, 50, 1e-7, verbose);
}

// BYM2
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_spatial_bym2(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p,
    Rcpp::IntegerMatrix n_valid, Rcpp::IntegerMatrix n_det,
    Rcpp::IntegerVector map_site_to_unit, int n_spatial,
    Rcpp::IntegerVector adj_row_ptr, Rcpp::IntegerVector adj_col_idx,
    Rcpp::IntegerVector n_neighbors, Rcpp::NumericVector sigma_grid,
    Rcpp::NumericVector rho_grid, double scale_factor,
    Rcpp::NumericVector mu_init, Rcpp::NumericMatrix Sigma_psi_init,
    Rcpp::NumericMatrix Sigma_p_init, int max_iter_em, bool verbose) {
    std::vector<tulpaObs::CommGridPoint> plan;
    std::vector<double> tg_sig, tg_rho;
    for (int ri = 0; ri < rho_grid.size(); ++ri)
        for (int si = 0; si < sigma_grid.size(); ++si) {
            const double sigma = sigma_grid[si], rho = rho_grid[ri];
            const double a = sigma * std::sqrt(rho / scale_factor);
            const double b = sigma * std::sqrt(1.0 - rho);
            plan.push_back({1.0, rho, 0.0, a, b, 0.0});
            tg_sig.push_back(sigma); tg_rho.push_back(rho);
        }
    Rcpp::NumericMatrix tg(plan.size(), 2);
    for (std::size_t i = 0; i < plan.size(); ++i) { tg(i, 0) = tg_sig[i]; tg(i, 1) = tg_rho[i]; }
    colnames(tg) = Rcpp::CharacterVector::create("sigma", "rho");
    return tulpaObs::run_occ_spatial_grid(
        tulpaObs::CommFieldKind::BYM2, 2 * n_spatial, X_psi, X_p, n_valid, n_det,
        map_site_to_unit, n_spatial, adj_row_ptr, adj_col_idx, n_neighbors,
        plan, tg, mu_init, Sigma_psi_init, Sigma_p_init,
        max_iter_em, 1e-4, 50, 1e-7, verbose);
}

// Continuous Matern (SPDE) shared field on the occupancy arm. The field lives
// at n_mesh FEM nodes; the dense projection A (n_sites x n_mesh) maps mesh nodes
// onto sites, shared across species. Per grid point the proper Matern precision
// Q(range, sigma) (carrying tau_spde^2) and its log|Q| are built once on the R
// side and passed in; the outer grid axes are (range, sigma). Mirrors
// cpp_nmix_community_spatial_spde (#239), minus the NB size nuisance occupancy
// has none of.
// [[Rcpp::export]]
Rcpp::List cpp_ms_occu_spatial_spde(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p,
    Rcpp::IntegerMatrix n_valid, Rcpp::IntegerMatrix n_det,
    Rcpp::NumericMatrix A_R,            // dense [n_sites x n_mesh]
    Rcpp::List Q_list,                  // per-grid-point precision [n_mesh x n_mesh]
    Rcpp::NumericVector log_det_Q,      // per-grid-point log|Q|
    Rcpp::NumericMatrix theta_grid_R,   // [n_grid x n_theta] (range, sigma)
    Rcpp::NumericVector mu_init, Rcpp::NumericMatrix Sigma_psi_init,
    Rcpp::NumericMatrix Sigma_p_init, int max_iter_em, bool verbose) {
    const int n_sites = A_R.nrow();
    const int n_mesh  = A_R.ncol();
    const MatrixXd A = Map<MatrixXd>(REAL(A_R), n_sites, n_mesh);
    return tulpaObs::run_occ_spatial_grid_spde(
        X_psi, X_p, n_valid, n_det, A, n_mesh, Q_list, log_det_Q, theta_grid_R,
        mu_init, Sigma_psi_init, Sigma_p_init,
        max_iter_em, 1e-4, 50, 1e-7, verbose);
}
