// cell_coupling_occu_multiscale_cover.cpp
// Per-fit registration of the stateful `OccuMultiscaleCover*Coupling` specs
// against tulpa's CellCouplingSpec registry, plus Rcpp-export'd direct
// single-cell evaluators used by
// tests/testthat/test-occu-multiscale-cover-coupling.R to FD-check every
// closed-form derivative against numerical derivatives of the cell density.
//
// Unlike the 2-level occu_cover specs (stateless, registered once at .onLoad),
// the multiscale spec carries the fit's per-cell plot structure, so the R
// driver re-registers it under a fixed name immediately before each joint fit
// (last-writer-wins; the previous fit's spec is released). See
// R/occu_multiscale_cover_joint.R.

#include "cell_coupling_occu_multiscale_cover.h"

#include <tulpa/cell_coupling.h>
#include <R_ext/Rdynload.h>
#include <Rcpp.h>
#include <memory>
#include <string>
#include <vector>

namespace {

inline tulpa::RegisterCellCouplingFn lookup_registrar() {
    auto fp = (tulpa::RegisterCellCouplingFn) R_GetCCallable(
        "tulpa", "tulpa_register_cell_coupling");
    if (!fp) {
        Rcpp::stop("tulpaObs: R_GetCCallable('tulpa', 'tulpa_register_cell_coupling') "
                   "returned NULL -- tulpa not loaded or ABI mismatch.");
    }
    return fp;
}

// Rebuild the ragged per-cell plot-size structure from the flat R inputs:
// `n_plots_per_cell[c]` plots for cell c, their visit counts read in order
// from `plot_sizes_flat` (cell-major, plot-major).
std::vector<std::vector<int>> build_cell_plot_sizes(
    Rcpp::IntegerVector n_plots_per_cell,
    Rcpp::IntegerVector plot_sizes_flat) {
    const int n_cells = n_plots_per_cell.size();
    std::vector<std::vector<int>> out(n_cells);
    int pos = 0;
    for (int c = 0; c < n_cells; c++) {
        const int mp = n_plots_per_cell[c];
        out[c].resize(mp);
        for (int j = 0; j < mp; j++) {
            if (pos >= plot_sizes_flat.size()) {
                Rcpp::stop("occu_multiscale_cover: plot_sizes_flat shorter than "
                           "sum(n_plots_per_cell).");
            }
            out[c][j] = plot_sizes_flat[pos++];
        }
    }
    return out;
}

// Shared body for the direct single-cell evaluators. Builds 4-arm
// CellEtas / CellResponse / CellDerivs views over one synthetic cell and
// dispatches into `spec.evaluate_cell()`.
template <class Spec>
Rcpp::List eval_one_cell_ms_(double                eta_psi,
                             Rcpp::NumericVector   eta_theta,
                             Rcpp::NumericVector   eta_p,
                             Rcpp::NumericVector   eta_pos,
                             Rcpp::IntegerVector   y_det,
                             Rcpp::NumericVector   y_pos,
                             Rcpp::IntegerVector   plot_sizes,
                             double                phi_pos,
                             const char*           fam_pos_str,
                             tulpa::CurvatureMode  curv) {
    const int M  = eta_theta.size();
    const int Jc = eta_p.size();
    if (eta_pos.size() != Jc || y_det.size() != Jc || y_pos.size() != Jc) {
        Rcpp::stop("cpp_eval_occu_multiscale_cover_*_cell: visit-length mismatch "
                   "(Jc=%d).", Jc);
    }
    if (plot_sizes.size() != M) {
        Rcpp::stop("cpp_eval_occu_multiscale_cover_*_cell: plot_sizes length %d "
                   "!= n plots %d.", (int)plot_sizes.size(), M);
    }
    int sum_sz = 0;
    for (int j = 0; j < M; j++) sum_sz += plot_sizes[j];
    if (sum_sz != Jc) {
        Rcpp::stop("cpp_eval_occu_multiscale_cover_*_cell: sum(plot_sizes)=%d "
                   "!= Jc=%d.", sum_sz, Jc);
    }

    std::vector<double> eta_psi_buf(1, eta_psi);
    std::vector<double> eta_theta_buf(eta_theta.begin(), eta_theta.end());
    std::vector<double> eta_p_buf(eta_p.begin(), eta_p.end());
    std::vector<double> eta_pos_buf(eta_pos.begin(), eta_pos.end());

    std::vector<const double*> arm_eta_ptr(4);
    arm_eta_ptr[0] = eta_psi_buf.data();
    arm_eta_ptr[1] = eta_theta_buf.data();
    arm_eta_ptr[2] = eta_p_buf.data();
    arm_eta_ptr[3] = eta_pos_buf.data();

    std::vector<int> rows_psi(1, 0);
    std::vector<int> rows_theta(M);   for (int j = 0; j < M;  j++) rows_theta[j] = j;
    std::vector<int> rows_p(Jc);      for (int v = 0; v < Jc; v++) rows_p[v] = v;
    std::vector<int> rows_pos(Jc);    for (int v = 0; v < Jc; v++) rows_pos[v] = v;
    std::vector<const int*> arm_rows(4);
    arm_rows[0] = rows_psi.data();
    arm_rows[1] = rows_theta.data();
    arm_rows[2] = rows_p.data();
    arm_rows[3] = rows_pos.data();
    std::vector<int> arm_row_count = {1, M, Jc, Jc};

    std::vector<double> y_det_buf(Jc);
    for (int v = 0; v < Jc; v++) y_det_buf[v] = (double) y_det[v];
    std::vector<double> y_pos_buf(y_pos.begin(), y_pos.end());
    std::vector<const double*> arm_y_ptr(4);
    arm_y_ptr[0] = nullptr;            // psi arm carries no data
    arm_y_ptr[1] = nullptr;            // theta arm carries no data
    arm_y_ptr[2] = y_det_buf.data();
    arm_y_ptr[3] = y_pos_buf.data();

    std::vector<int> n_trials_psi(1, 0);
    std::vector<int> n_trials_theta(M, 0);
    std::vector<int> n_trials_p(Jc, 1);
    std::vector<int> n_trials_pos(Jc, 0);
    std::vector<const int*> arm_n_trials_ptr(4);
    arm_n_trials_ptr[0] = n_trials_psi.data();
    arm_n_trials_ptr[1] = n_trials_theta.data();
    arm_n_trials_ptr[2] = n_trials_p.data();
    arm_n_trials_ptr[3] = n_trials_pos.data();

    std::string fam_psi = "binomial", fam_theta = "binomial", fam_p = "binomial";
    std::string fam_pos = fam_pos_str;
    std::vector<const char*> arm_family_ptr = {
        fam_psi.c_str(), fam_theta.c_str(), fam_p.c_str(), fam_pos.c_str()
    };
    std::vector<double> arm_phi = {1.0, 1.0, 1.0, phi_pos};

    tulpa::CellEtas etas_view;
    etas_view.arm_eta_ptr   = arm_eta_ptr.data();
    etas_view.arm_rows      = arm_rows.data();
    etas_view.arm_row_count = arm_row_count.data();
    etas_view.n_arms_       = 4;

    tulpa::CellResponse y_view;
    y_view.arm_y         = arm_y_ptr.data();
    y_view.arm_n_trials  = arm_n_trials_ptr.data();
    y_view.arm_family    = arm_family_ptr.data();
    y_view.arm_phi       = arm_phi.data();
    y_view.arm_rows      = arm_rows.data();
    y_view.arm_row_count = arm_row_count.data();
    y_view.n_arms_       = 4;

    std::vector<double> grad_psi_buf(1, 0.0), grad_theta_buf(M, 0.0),
                        grad_p_buf(Jc, 0.0), grad_pos_buf(Jc, 0.0);
    std::vector<double> nh_psi_buf(1, 0.0), nh_theta_buf(M, 0.0),
                        nh_p_buf(Jc, 0.0), nh_pos_buf(Jc, 0.0);
    std::vector<double*> arm_grad_ptr = {grad_psi_buf.data(), grad_theta_buf.data(),
                                         grad_p_buf.data(),   grad_pos_buf.data()};
    std::vector<double*> arm_nh_ptr = {nh_psi_buf.data(), nh_theta_buf.data(),
                                       nh_p_buf.data(),   nh_pos_buf.data()};

    // Cross-Hessian buffers for every kk <= ll pair (row counts {1, M, Jc, Jc}).
    auto rc = [&](int k) { return arm_row_count[k]; };
    std::vector<std::vector<std::vector<double>>> cross(4,
        std::vector<std::vector<double>>(4));
    for (int k = 0; k < 4; k++)
        for (int l = k; l < 4; l++)
            cross[k][l].assign((std::size_t)rc(k) * rc(l), 0.0);
    std::vector<std::vector<double*>> cross_inner(4, std::vector<double*>(4, nullptr));
    for (int k = 0; k < 4; k++)
        for (int l = k; l < 4; l++)
            cross_inner[k][l] = cross[k][l].data();
    std::vector<double* const*> cross_outer(4);
    for (int k = 0; k < 4; k++) cross_outer[k] = cross_inner[k].data();

    tulpa::CellDerivs out;
    out.arm_grad          = arm_grad_ptr.data();
    out.arm_neg_hess_diag = arm_nh_ptr.data();
    out.arm_cross_hess    = cross_outer.data();
    out.arm_row_count     = arm_row_count.data();
    out.n_arms_           = 4;
    out.curvature         = curv;

    std::vector<std::vector<int>> cps(1);
    cps[0].assign(plot_sizes.begin(), plot_sizes.end());
    Spec spec(std::move(cps));
    const double cell_ll = spec.evaluate_cell(0, etas_view, y_view, out);

    auto to_mat = [](const std::vector<double>& buf, int nr, int nc) {
        Rcpp::NumericMatrix mm(nr, nc);
        for (int i = 0; i < nr; i++)
            for (int j = 0; j < nc; j++)
                mm(i, j) = buf[(std::size_t)i * nc + j];
        return mm;
    };

    return Rcpp::List::create(
        Rcpp::Named("cell_ll")          = cell_ll,
        Rcpp::Named("grad_psi")         = grad_psi_buf[0],
        Rcpp::Named("grad_theta")       = Rcpp::NumericVector(grad_theta_buf.begin(), grad_theta_buf.end()),
        Rcpp::Named("grad_p")           = Rcpp::NumericVector(grad_p_buf.begin(), grad_p_buf.end()),
        Rcpp::Named("grad_pos")         = Rcpp::NumericVector(grad_pos_buf.begin(), grad_pos_buf.end()),
        Rcpp::Named("neg_hess_psi")     = nh_psi_buf[0],
        Rcpp::Named("neg_hess_theta")   = Rcpp::NumericVector(nh_theta_buf.begin(), nh_theta_buf.end()),
        Rcpp::Named("neg_hess_p")       = Rcpp::NumericVector(nh_p_buf.begin(), nh_p_buf.end()),
        Rcpp::Named("neg_hess_pos")     = Rcpp::NumericVector(nh_pos_buf.begin(), nh_pos_buf.end()),
        Rcpp::Named("cross_psi_theta")  = Rcpp::NumericVector(cross[0][1].begin(), cross[0][1].end()),
        Rcpp::Named("cross_psi_p")      = Rcpp::NumericVector(cross[0][2].begin(), cross[0][2].end()),
        Rcpp::Named("cross_theta_theta")= to_mat(cross[1][1], M, M),
        Rcpp::Named("cross_theta_p")    = to_mat(cross[1][2], M, Jc),
        Rcpp::Named("cross_p_p")        = to_mat(cross[2][2], Jc, Jc)
    );
}

static tulpa::CurvatureMode parse_curvature_(const std::string& s) {
    return (s == "expected" || s == "fisher")
             ? tulpa::CurvatureMode::Expected
             : tulpa::CurvatureMode::Observed;
}

} // namespace


// Register the stateful multiscale spec for the requested positive family,
// carrying this fit's per-cell plot structure. Returns the registry name the
// R driver passes as `cell_coupling`.
// [[Rcpp::export]]
std::string cpp_register_occu_multiscale_cover_coupling(
    std::string         positive,
    Rcpp::IntegerVector n_plots_per_cell,
    Rcpp::IntegerVector plot_sizes_flat) {
    auto cps = build_cell_plot_sizes(n_plots_per_cell, plot_sizes_flat);
    auto fp  = lookup_registrar();
    if (positive == "beta") {
        std::string nm = tulpaObs::BetaPositive::multiscale_spec_name();
        fp(nm.c_str(),
           std::make_shared<tulpaObs::OccuMultiscaleCoverBetaCoupling>(std::move(cps)));
        return nm;
    } else if (positive == "lognormal") {
        std::string nm = tulpaObs::LognormalPositive::multiscale_spec_name();
        fp(nm.c_str(),
           std::make_shared<tulpaObs::OccuMultiscaleCoverLognormalCoupling>(std::move(cps)));
        return nm;
    }
    Rcpp::stop("cpp_register_occu_multiscale_cover_coupling: positive must be "
               "\"lognormal\" or \"beta\".");
}


// Direct single-cell evaluators (FD test harness). `plot_sizes` partitions the
// length-Jc visit vectors into the M plots (sum(plot_sizes) == Jc); arm-1
// (theta) is length M. `sigma_pos` / `phi_pos` is the pos-arm dispersion.
// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_multiscale_cover_lognormal_cell(
    double                eta_psi,
    Rcpp::NumericVector   eta_theta,
    Rcpp::NumericVector   eta_p,
    Rcpp::NumericVector   eta_pos,
    Rcpp::IntegerVector   y_det,
    Rcpp::NumericVector   y_pos,
    Rcpp::IntegerVector   plot_sizes,
    double                sigma_pos,
    std::string           curvature = "observed") {
    return eval_one_cell_ms_<tulpaObs::OccuMultiscaleCoverLognormalCoupling>(
        eta_psi, eta_theta, eta_p, eta_pos, y_det, y_pos, plot_sizes,
        sigma_pos, "lognormal", parse_curvature_(curvature));
}

// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_multiscale_cover_beta_cell(
    double                eta_psi,
    Rcpp::NumericVector   eta_theta,
    Rcpp::NumericVector   eta_p,
    Rcpp::NumericVector   eta_pos,
    Rcpp::IntegerVector   y_det,
    Rcpp::NumericVector   y_pos,
    Rcpp::IntegerVector   plot_sizes,
    double                phi_pos,
    std::string           curvature = "observed") {
    return eval_one_cell_ms_<tulpaObs::OccuMultiscaleCoverBetaCoupling>(
        eta_psi, eta_theta, eta_p, eta_pos, y_det, y_pos, plot_sizes,
        phi_pos, "beta", parse_curvature_(curvature));
}
