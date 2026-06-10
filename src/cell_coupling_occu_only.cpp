// cell_coupling_occu_only.cpp
// Registration of `OccuOnlyCoupling` against tulpa's CellCouplingSpec registry
// via the `tulpa_register_cell_coupling` registered C callable, plus an
// Rcpp-export'd direct evaluator used by tests/testthat/test-occu-only-coupling.R
// to FD-check every closed-form derivative against numerical derivatives of the
// cell log-density (gcol33/tulpaObs#81).
//
// Registration is invoked from R via `.onLoad` (R/tulpaObs-package.R), the same
// pattern the occu_cover specs use.

#include "cell_coupling_occu_only.h"

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

// Map the test-facing curvature string to the engine enum. "expected" /
// "fisher" select the complete-data Fisher curvature; anything else is the
// observed Hessian (the default the inner solver uses at the mode-pass).
static tulpa::CurvatureMode parse_curvature_(const std::string& s) {
    return (s == "expected" || s == "fisher")
             ? tulpa::CurvatureMode::Expected
             : tulpa::CurvatureMode::Observed;
}

} // namespace


// [[Rcpp::export]]
void cpp_register_occu_only_coupling() {
    auto fp = lookup_registrar();
    fp("occu_only", std::make_shared<tulpaObs::OccuOnlyCoupling>());
}


// Direct evaluator returning the spec's closed-form log-density and every
// nonzero derivative buffer for a single synthetic cell. The R-side FD test
// rebuilds the same cell density by varying each eta numerically and compares.
//
// `y_det` is a length-J 0/1 integer vector; `eta_p` is length J. `eta_psi` is
// the occupancy logit. There is no cover arm.
//
// Returns a list with `cell_ll` (scalar) plus `grad_psi` (scalar), `grad_p`
// (length J), `neg_hess_psi` (scalar), `neg_hess_p` (length J), `cross_psi_p`
// (length J; zero in the det case) and `cross_p_p` (J x J row-major; zero in the
// det case and along its diagonal).
// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_only_cell(
    double                     eta_psi,
    Rcpp::NumericVector        eta_p,
    Rcpp::IntegerVector        y_det,
    std::string                curvature = "observed"
) {
    const int Jc = eta_p.size();
    if (y_det.size() != Jc) {
        Rcpp::stop("cpp_eval_occu_only_cell: length mismatch (Jc=%d, y_det=%d).",
                   Jc, (int) y_det.size());
    }
    const tulpa::CurvatureMode curv = parse_curvature_(curvature);

    std::vector<double> eta_psi_buf(1, eta_psi);
    std::vector<double> eta_p_buf(eta_p.begin(), eta_p.end());
    std::vector<const double*> arm_eta_ptr(2);
    arm_eta_ptr[0] = eta_psi_buf.data();
    arm_eta_ptr[1] = eta_p_buf.data();

    std::vector<int> rows_psi(1, 0);
    std::vector<int> rows_p(Jc); for (int v = 0; v < Jc; v++) rows_p[v] = v;
    std::vector<const int*> arm_rows(2);
    arm_rows[0] = rows_psi.data();
    arm_rows[1] = rows_p.data();
    std::vector<int> arm_row_count = {1, Jc};

    std::vector<double> y_det_buf(Jc);
    for (int v = 0; v < Jc; v++) y_det_buf[v] = (double) y_det[v];
    std::vector<const double*> arm_y_ptr(2);
    arm_y_ptr[0] = nullptr;
    arm_y_ptr[1] = y_det_buf.data();

    std::vector<int> n_trials_psi(1, 0);
    std::vector<int> n_trials_p(Jc, 1);
    std::vector<const int*> arm_n_trials_ptr(2);
    arm_n_trials_ptr[0] = n_trials_psi.data();
    arm_n_trials_ptr[1] = n_trials_p.data();

    std::string fam_psi = "binomial", fam_p = "binomial";
    std::vector<const char*> arm_family_ptr = { fam_psi.c_str(), fam_p.c_str() };
    std::vector<double> arm_phi = {1.0, 1.0};

    tulpa::CellEtas etas_view;
    etas_view.arm_eta_ptr   = arm_eta_ptr.data();
    etas_view.arm_rows      = arm_rows.data();
    etas_view.arm_row_count = arm_row_count.data();
    etas_view.n_arms_       = 2;

    tulpa::CellResponse y_view;
    y_view.arm_y         = arm_y_ptr.data();
    y_view.arm_n_trials  = arm_n_trials_ptr.data();
    y_view.arm_family    = arm_family_ptr.data();
    y_view.arm_phi       = arm_phi.data();
    y_view.arm_rows      = arm_rows.data();
    y_view.arm_row_count = arm_row_count.data();
    y_view.n_arms_       = 2;

    std::vector<double> grad_psi_buf(1, 0.0);
    std::vector<double> grad_p_buf(Jc, 0.0);
    std::vector<double> neg_hess_psi_buf(1, 0.0);
    std::vector<double> neg_hess_p_buf(Jc, 0.0);
    std::vector<double*> arm_grad_ptr = {grad_psi_buf.data(), grad_p_buf.data()};
    std::vector<double*> arm_neg_hess_diag_ptr = {neg_hess_psi_buf.data(),
                                                   neg_hess_p_buf.data()};

    std::vector<double> cross_00(1 * 1, 0.0);
    std::vector<double> cross_01(1 * Jc, 0.0);
    std::vector<double> cross_11((std::size_t)Jc * Jc, 0.0);
    std::vector<double*> cross_row0 = {cross_00.data(), cross_01.data()};
    std::vector<double*> cross_row1 = {nullptr,         cross_11.data()};
    std::vector<double* const*> cross_outer = {cross_row0.data(),
                                               cross_row1.data()};

    tulpa::CellDerivs out;
    out.arm_grad           = arm_grad_ptr.data();
    out.arm_neg_hess_diag  = arm_neg_hess_diag_ptr.data();
    out.arm_cross_hess     = cross_outer.data();
    out.arm_row_count      = arm_row_count.data();
    out.n_arms_            = 2;
    out.curvature          = curv;

    tulpaObs::OccuOnlyCoupling spec;
    const double cell_ll = spec.evaluate_cell(0, etas_view, y_view, out);

    Rcpp::NumericVector grad_p_r(grad_p_buf.begin(), grad_p_buf.end());
    Rcpp::NumericVector neg_hess_p_r(neg_hess_p_buf.begin(), neg_hess_p_buf.end());
    Rcpp::NumericVector cross_psi_p_r(cross_01.begin(), cross_01.end());
    Rcpp::NumericMatrix cross_p_p_r(Jc, Jc);
    for (int j = 0; j < Jc; j++) {
        for (int m = 0; m < Jc; m++) {
            cross_p_p_r(j, m) = cross_11[(std::size_t)j * Jc + m];
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("cell_ll")      = cell_ll,
        Rcpp::Named("grad_psi")     = grad_psi_buf[0],
        Rcpp::Named("grad_p")       = grad_p_r,
        Rcpp::Named("neg_hess_psi") = neg_hess_psi_buf[0],
        Rcpp::Named("neg_hess_p")   = neg_hess_p_r,
        Rcpp::Named("cross_psi_p")  = cross_psi_p_r,
        Rcpp::Named("cross_p_p")    = cross_p_p_r
    );
}
