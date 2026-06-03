// cell_coupling_occu_cover.cpp
// Registration of `OccuCoverLognormalCoupling` + `OccuCoverBetaCoupling`
// against tulpa's CellCouplingSpec registry via the
// `tulpa_register_cell_coupling` registered C callable, plus Rcpp-export'd
// direct evaluators used by tests/testthat/test-occu-cover-coupling.R to
// FD-check every closed-form derivative against numerical derivatives of
// the cell log-density.
//
// Registration is invoked from R via `.onLoad` (R/tulpaObs-package.R)
// rather than from a `R_init_tulpaObs` C entry point -- Rcpp's
// compileAttributes already owns R_init, and a single `tulpa_register_*`
// call on package load is the standard pattern across the tulpa
// ecosystem.

#include "cell_coupling_occu_cover.h"

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

// Shared body: builds CellEtas / CellResponse / CellDerivs views over a
// single synthetic cell and dispatches into `spec.evaluate_cell()`. The
// caller picks the spec (lognormal or beta). Used by both direct
// evaluators below.
template <class Spec>
Rcpp::List eval_one_cell_(double                eta_psi,
                          Rcpp::NumericVector   eta_p,
                          Rcpp::NumericVector   eta_pos,
                          Rcpp::IntegerVector   y_det,
                          Rcpp::NumericVector   y_pos,
                          double                phi_pos,
                          const char*           fam_pos_str,
                          tulpa::CurvatureMode  curv) {
    const int Jc = eta_p.size();
    // The aggregated specs carry one pos row per cell (the mean / median cover),
    // so eta_pos / y_pos are length 1 there; the per-visit specs carry Jc pos
    // rows. Size the pos arm from eta_pos rather than assuming it equals Jc.
    const int n_pos = eta_pos.size();
    if (y_det.size() != Jc || y_pos.size() != n_pos ||
        (n_pos != Jc && n_pos != 1)) {
        Rcpp::stop("cpp_eval_occu_cover_*_cell: length mismatch "
                   "(Jc=%d, n_pos=%d).", Jc, n_pos);
    }

    std::vector<double> eta_psi_buf(1, eta_psi);
    std::vector<double> eta_p_buf(eta_p.begin(), eta_p.end());
    std::vector<double> eta_pos_buf(eta_pos.begin(), eta_pos.end());

    std::vector<const double*> arm_eta_ptr(3);
    arm_eta_ptr[0] = eta_psi_buf.data();
    arm_eta_ptr[1] = eta_p_buf.data();
    arm_eta_ptr[2] = eta_pos_buf.data();

    std::vector<int> rows_psi(1, 0);
    std::vector<int> rows_p(Jc); for (int v = 0; v < Jc; v++) rows_p[v] = v;
    std::vector<int> rows_pos(n_pos); for (int v = 0; v < n_pos; v++) rows_pos[v] = v;
    std::vector<const int*> arm_rows(3);
    arm_rows[0] = rows_psi.data();
    arm_rows[1] = rows_p.data();
    arm_rows[2] = rows_pos.data();
    std::vector<int> arm_row_count = {1, Jc, n_pos};

    std::vector<double> y_dummy_psi(1, 0.0);
    std::vector<double> y_det_buf(Jc);
    for (int v = 0; v < Jc; v++) y_det_buf[v] = (double) y_det[v];
    std::vector<double> y_pos_buf(y_pos.begin(), y_pos.end());
    std::vector<const double*> arm_y_ptr(3);
    arm_y_ptr[0] = nullptr;
    arm_y_ptr[1] = y_det_buf.data();
    arm_y_ptr[2] = y_pos_buf.data();

    std::vector<int> n_trials_psi(1, 0);
    std::vector<int> n_trials_p(Jc, 1);
    std::vector<int> n_trials_pos(n_pos, 0);
    std::vector<const int*> arm_n_trials_ptr(3);
    arm_n_trials_ptr[0] = n_trials_psi.data();
    arm_n_trials_ptr[1] = n_trials_p.data();
    arm_n_trials_ptr[2] = n_trials_pos.data();

    std::string fam_psi = "binomial", fam_p = "binomial";
    std::string fam_pos = fam_pos_str;
    std::vector<const char*> arm_family_ptr = {
        fam_psi.c_str(), fam_p.c_str(), fam_pos.c_str()
    };
    std::vector<double> arm_phi = {1.0, 1.0, phi_pos};

    tulpa::CellEtas etas_view;
    etas_view.arm_eta_ptr   = arm_eta_ptr.data();
    etas_view.arm_rows      = arm_rows.data();
    etas_view.arm_row_count = arm_row_count.data();
    etas_view.n_arms_       = 3;

    tulpa::CellResponse y_view;
    y_view.arm_y         = arm_y_ptr.data();
    y_view.arm_n_trials  = arm_n_trials_ptr.data();
    y_view.arm_family    = arm_family_ptr.data();
    y_view.arm_phi       = arm_phi.data();
    y_view.arm_rows      = arm_rows.data();
    y_view.arm_row_count = arm_row_count.data();
    y_view.n_arms_       = 3;

    std::vector<double> grad_psi_buf(1, 0.0);
    std::vector<double> grad_p_buf(Jc, 0.0);
    std::vector<double> grad_pos_buf(n_pos, 0.0);
    std::vector<double> neg_hess_psi_buf(1, 0.0);
    std::vector<double> neg_hess_p_buf(Jc, 0.0);
    std::vector<double> neg_hess_pos_buf(n_pos, 0.0);
    std::vector<double*> arm_grad_ptr = {grad_psi_buf.data(),
                                          grad_p_buf.data(),
                                          grad_pos_buf.data()};
    std::vector<double*> arm_neg_hess_diag_ptr = {neg_hess_psi_buf.data(),
                                                   neg_hess_p_buf.data(),
                                                   neg_hess_pos_buf.data()};

    std::vector<double> cross_00(1 * 1, 0.0);
    std::vector<double> cross_01(1 * Jc, 0.0);
    std::vector<double> cross_02((std::size_t)1 * n_pos, 0.0);
    std::vector<double> cross_11((std::size_t)Jc * Jc, 0.0);
    std::vector<double> cross_12((std::size_t)Jc * n_pos, 0.0);
    std::vector<double> cross_22((std::size_t)n_pos * n_pos, 0.0);

    std::vector<double*> cross_row0 = {cross_00.data(), cross_01.data(), cross_02.data()};
    std::vector<double*> cross_row1 = {nullptr,         cross_11.data(), cross_12.data()};
    std::vector<double*> cross_row2 = {nullptr,         nullptr,         cross_22.data()};
    std::vector<double* const*> cross_outer = {cross_row0.data(),
                                                cross_row1.data(),
                                                cross_row2.data()};

    tulpa::CellDerivs out;
    out.arm_grad           = arm_grad_ptr.data();
    out.arm_neg_hess_diag  = arm_neg_hess_diag_ptr.data();
    out.arm_cross_hess     = cross_outer.data();
    out.arm_row_count      = arm_row_count.data();
    out.n_arms_            = 3;
    out.curvature          = curv;

    Spec spec;
    const double cell_ll = spec.evaluate_cell(0, etas_view, y_view, out);

    Rcpp::NumericVector grad_p_r(grad_p_buf.begin(), grad_p_buf.end());
    Rcpp::NumericVector grad_pos_r(grad_pos_buf.begin(), grad_pos_buf.end());
    Rcpp::NumericVector neg_hess_p_r(neg_hess_p_buf.begin(), neg_hess_p_buf.end());
    Rcpp::NumericVector neg_hess_pos_r(neg_hess_pos_buf.begin(), neg_hess_pos_buf.end());
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
        Rcpp::Named("grad_pos")     = grad_pos_r,
        Rcpp::Named("neg_hess_psi") = neg_hess_psi_buf[0],
        Rcpp::Named("neg_hess_p")   = neg_hess_p_r,
        Rcpp::Named("neg_hess_pos") = neg_hess_pos_r,
        Rcpp::Named("cross_psi_p")  = cross_psi_p_r,
        Rcpp::Named("cross_p_p")    = cross_p_p_r
    );
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
void cpp_register_occu_cover_lognormal_coupling() {
    auto fp = lookup_registrar();
    fp("occu_cover_lognormal",
       std::make_shared<tulpaObs::OccuCoverLognormalCoupling>());
}

// [[Rcpp::export]]
void cpp_register_occu_cover_beta_coupling() {
    auto fp = lookup_registrar();
    fp("occu_cover_beta",
       std::make_shared<tulpaObs::OccuCoverBetaCoupling>());
}

// Cell-aggregated cover variants (tulpaObs#33): the pos arm carries one row per
// detected occupancy unit (the mean / median cover), evaluated once per cell.
// [[Rcpp::export]]
void cpp_register_occu_cover_lognormal_agg_coupling() {
    auto fp = lookup_registrar();
    fp("occu_cover_lognormal_agg",
       std::make_shared<tulpaObs::OccuCoverLognormalAggCoupling>());
}

// [[Rcpp::export]]
void cpp_register_occu_cover_beta_agg_coupling() {
    auto fp = lookup_registrar();
    fp("occu_cover_beta_agg",
       std::make_shared<tulpaObs::OccuCoverBetaAggCoupling>());
}


// Direct evaluator returning the spec's closed-form log-density and every
// nonzero derivative buffer for a single synthetic cell. R-side FD test
// rebuilds the same cell density by varying each eta numerically and
// compares.
//
// `y_det` is a length-J 0/1 integer vector; `y_pos` is the raw lognormal
// data at detected visits (0 ignored at undetected). `eta_p` and
// `eta_pos` are length-J. `sigma_pos` is the lognormal SD on the log
// scale (= the pos arm's phi).
//
// Returns a list with `cell_ll` (scalar) plus `grad_psi` (scalar),
// `grad_p` (length J), `grad_pos` (length J), `neg_hess_psi` (scalar),
// `neg_hess_p` (length J), `neg_hess_pos` (length J), `cross_psi_p`
// (length J; zero in det case), `cross_p_p` (J x J row-major; zero in
// det case and along its diagonal).
// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_cover_lognormal_cell(
    double                     eta_psi,
    Rcpp::NumericVector        eta_p,
    Rcpp::NumericVector        eta_pos,
    Rcpp::IntegerVector        y_det,
    Rcpp::NumericVector        y_pos,
    double                     sigma_pos,
    std::string                curvature = "observed"
) {
    return eval_one_cell_<tulpaObs::OccuCoverLognormalCoupling>(
        eta_psi, eta_p, eta_pos, y_det, y_pos, sigma_pos, "lognormal",
        parse_curvature_(curvature)
    );
}

// Beta-arm twin of `cpp_eval_occu_cover_lognormal_cell`. `y_pos` is in
// (0, 1) at detected visits and 0 at undetected visits; `phi_pos` is the
// beta precision (mean is sigmoid(eta_pos), variance is mu(1-mu)/(1+phi)).
// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_cover_beta_cell(
    double                     eta_psi,
    Rcpp::NumericVector        eta_p,
    Rcpp::NumericVector        eta_pos,
    Rcpp::IntegerVector        y_det,
    Rcpp::NumericVector        y_pos,
    double                     phi_pos,
    std::string                curvature = "observed"
) {
    return eval_one_cell_<tulpaObs::OccuCoverBetaCoupling>(
        eta_psi, eta_p, eta_pos, y_det, y_pos, phi_pos, "beta",
        parse_curvature_(curvature)
    );
}

// Cell-aggregated twins (tulpaObs#33). `eta_pos` and `y_pos` are length 1 (the
// cell's aggregated cover predictor / mean-or-median observation); `eta_p` /
// `y_det` stay length J. The returned `grad_pos` / `neg_hess_pos` are length 1.
// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_cover_lognormal_agg_cell(
    double                     eta_psi,
    Rcpp::NumericVector        eta_p,
    Rcpp::NumericVector        eta_pos,
    Rcpp::IntegerVector        y_det,
    Rcpp::NumericVector        y_pos,
    double                     sigma_pos,
    std::string                curvature = "observed"
) {
    return eval_one_cell_<tulpaObs::OccuCoverLognormalAggCoupling>(
        eta_psi, eta_p, eta_pos, y_det, y_pos, sigma_pos, "lognormal",
        parse_curvature_(curvature)
    );
}

// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_cover_beta_agg_cell(
    double                     eta_psi,
    Rcpp::NumericVector        eta_p,
    Rcpp::NumericVector        eta_pos,
    Rcpp::IntegerVector        y_det,
    Rcpp::NumericVector        y_pos,
    double                     phi_pos,
    std::string                curvature = "observed"
) {
    return eval_one_cell_<tulpaObs::OccuCoverBetaAggCoupling>(
        eta_psi, eta_p, eta_pos, y_det, y_pos, phi_pos, "beta",
        parse_curvature_(curvature)
    );
}
