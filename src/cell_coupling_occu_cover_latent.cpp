// cell_coupling_occu_cover_latent.cpp
// Per-fit registration of the stateful latent-cover occu_cover specs
// (`occu_cover_lognormal_latent` / `occu_cover_beta_latent`) plus Rcpp-export'd
// direct evaluators used by tests/testthat/test-occu-cover-latent.R to FD-check
// the per-unit latent marginal and its eta-derivatives.
//
// Unlike the stateless per-visit / aggregated specs (registered once at
// .onLoad), the latent spec captures the per-unit cover data + the fixed
// within-unit dispersion at construction, so it is (re)registered from R per
// fit with that fit's data. The registry is last-writer-wins and hands the
// joint driver a shared_ptr that outlives the lookup, so sequential fits sharing
// the name are safe.

#include "cell_coupling_occu_cover_latent.h"

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

// Convert an R list of per-unit cover-value vectors to the C++ ragged layout
// the stateful spec captures (indexed by global pos-arm row).
std::vector<std::vector<double>> as_site_data(const Rcpp::List& cover) {
    std::vector<std::vector<double>> out;
    out.reserve(cover.size());
    for (R_xlen_t i = 0; i < cover.size(); ++i) {
        Rcpp::NumericVector v = cover[i];
        out.emplace_back(v.begin(), v.end());
    }
    return out;
}

static tulpa::CurvatureMode parse_curvature_(const std::string& s) {
    return (s == "expected" || s == "fisher")
             ? tulpa::CurvatureMode::Expected
             : tulpa::CurvatureMode::Observed;
}

// Build a stateful latent spec over a single synthetic cell holding ONE
// occupancy unit's cover data, drive evaluate_cell, and return cell_ll + the
// nonzero derivative buffers. Templated on the latent policy.
template <class PosLatent>
Rcpp::List eval_one_latent_cell_(double               eta_psi,
                                 Rcpp::NumericVector  eta_p,
                                 double               eta_pos,
                                 Rcpp::IntegerVector  y_det,
                                 Rcpp::NumericVector  y_pos_vals,
                                 double               disp2,
                                 double               sigma_u,
                                 int                  n_quad,
                                 tulpa::CurvatureMode curv) {
    const int Jc = eta_p.size();
    if (y_det.size() != Jc) {
        Rcpp::stop("eval_one_latent_cell_: y_det length %d != Jc %d.",
                   (int) y_det.size(), Jc);
    }

    // One site, captured by global pos row 0.
    std::vector<std::vector<double>> site_data(
        1, std::vector<double>(y_pos_vals.begin(), y_pos_vals.end()));
    tulpaObs::OccuCoverLatentCoupling<PosLatent> spec(site_data, disp2, n_quad);

    std::vector<double> eta_psi_buf(1, eta_psi);
    std::vector<double> eta_p_buf(eta_p.begin(), eta_p.end());
    std::vector<double> eta_pos_buf(1, eta_pos);
    std::vector<const double*> arm_eta_ptr = {
        eta_psi_buf.data(), eta_p_buf.data(), eta_pos_buf.data() };

    std::vector<int> rows_psi(1, 0);
    std::vector<int> rows_p(Jc); for (int v = 0; v < Jc; v++) rows_p[v] = v;
    std::vector<int> rows_pos(1, 0);
    std::vector<const int*> arm_rows = {
        rows_psi.data(), rows_p.data(), rows_pos.data() };
    std::vector<int> arm_row_count = {1, Jc, 1};

    std::vector<double> y_det_buf(Jc);
    for (int v = 0; v < Jc; v++) y_det_buf[v] = (double) y_det[v];
    std::vector<double> y_pos_dummy(1, 0.0);   // pos arm y unused by the latent spec
    std::vector<const double*> arm_y_ptr = {
        nullptr, y_det_buf.data(), y_pos_dummy.data() };

    std::vector<int> n_trials_psi(1, 0), n_trials_p(Jc, 1), n_trials_pos(1, 0);
    std::vector<const int*> arm_n_trials_ptr = {
        n_trials_psi.data(), n_trials_p.data(), n_trials_pos.data() };

    std::string fam_psi = "binomial", fam_p = "binomial",
                fam_pos = PosLatent::spec_name();
    std::vector<const char*> arm_family_ptr = {
        fam_psi.c_str(), fam_p.c_str(), fam_pos.c_str() };
    std::vector<double> arm_phi = {1.0, 1.0, sigma_u};   // pos phi slot = sigma_u

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

    std::vector<double> grad_psi_buf(1, 0.0), grad_p_buf(Jc, 0.0), grad_pos_buf(1, 0.0);
    std::vector<double> nh_psi_buf(1, 0.0), nh_p_buf(Jc, 0.0), nh_pos_buf(1, 0.0);
    std::vector<double*> arm_grad_ptr = {
        grad_psi_buf.data(), grad_p_buf.data(), grad_pos_buf.data() };
    std::vector<double*> arm_nh_ptr = {
        nh_psi_buf.data(), nh_p_buf.data(), nh_pos_buf.data() };

    // Cross-Hessian buffers (zero in the det branch; populated in nodet).
    std::vector<double> cross_00(1, 0.0), cross_01((std::size_t) Jc, 0.0),
                        cross_02(1, 0.0), cross_11((std::size_t) Jc * Jc, 0.0),
                        cross_12((std::size_t) Jc, 0.0), cross_22(1, 0.0);
    std::vector<double*> cross_row0 = {cross_00.data(), cross_01.data(), cross_02.data()};
    std::vector<double*> cross_row1 = {nullptr, cross_11.data(), cross_12.data()};
    std::vector<double*> cross_row2 = {nullptr, nullptr, cross_22.data()};
    std::vector<double* const*> cross_outer = {
        cross_row0.data(), cross_row1.data(), cross_row2.data() };

    tulpa::CellDerivs out;
    out.arm_grad          = arm_grad_ptr.data();
    out.arm_neg_hess_diag = arm_nh_ptr.data();
    out.arm_cross_hess    = cross_outer.data();
    out.arm_row_count     = arm_row_count.data();
    out.n_arms_           = 3;
    out.curvature         = curv;

    const double cell_ll = spec.evaluate_cell(0, etas_view, y_view, out);

    return Rcpp::List::create(
        Rcpp::Named("cell_ll")      = cell_ll,
        Rcpp::Named("grad_psi")     = grad_psi_buf[0],
        Rcpp::Named("grad_p")       = Rcpp::NumericVector(grad_p_buf.begin(), grad_p_buf.end()),
        Rcpp::Named("grad_pos")     = grad_pos_buf[0],
        Rcpp::Named("neg_hess_psi") = nh_psi_buf[0],
        Rcpp::Named("neg_hess_p")   = Rcpp::NumericVector(nh_p_buf.begin(), nh_p_buf.end()),
        Rcpp::Named("neg_hess_pos") = nh_pos_buf[0]
    );
}

} // namespace


// ---------------------------------------------------------------------------
// Per-fit registration.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
void cpp_register_occu_cover_lognormal_latent_coupling(Rcpp::List cover_values,
                                                       double     sigma_eps,
                                                       int        n_quad = 1) {
    auto fp = lookup_registrar();
    fp("occu_cover_lognormal_latent",
       std::make_shared<tulpaObs::OccuCoverLognormalLatentCoupling>(
           as_site_data(cover_values), sigma_eps, n_quad));
}

// [[Rcpp::export]]
void cpp_register_occu_cover_beta_latent_coupling(Rcpp::List cover_values,
                                                  double     phi_prec,
                                                  int        n_quad = 15) {
    auto fp = lookup_registrar();
    fp("occu_cover_beta_latent",
       std::make_shared<tulpaObs::OccuCoverBetaLatentCoupling>(
           as_site_data(cover_values), phi_prec, n_quad));
}


// ---------------------------------------------------------------------------
// Direct FD evaluators (single synthetic cell, one occupancy unit).
// ---------------------------------------------------------------------------

// `y_pos_vals` is the unit's detected cover values (length m, on the natural
// scale). `sigma_eps` is the fixed within-unit log-scale SD; `sigma_u` the
// cover-latent SD. Returns cell_ll + grad/neg_hess for psi, p, and the single
// pos row.
// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_cover_lognormal_latent_cell(
    double               eta_psi,
    Rcpp::NumericVector  eta_p,
    double               eta_pos,
    Rcpp::IntegerVector  y_det,
    Rcpp::NumericVector  y_pos_vals,
    double               sigma_eps,
    double               sigma_u,
    std::string          curvature = "observed") {
    return eval_one_latent_cell_<tulpaObs::LognormalLatent>(
        eta_psi, eta_p, eta_pos, y_det, y_pos_vals, sigma_eps, sigma_u,
        1, parse_curvature_(curvature));
}

// Beta twin. `y_pos_vals` in (0, 1); `phi_prec` is the fixed beta precision;
// `sigma_u` the cover-latent SD. `n_quad` is the adaptive GH node count.
// [[Rcpp::export]]
Rcpp::List cpp_eval_occu_cover_beta_latent_cell(
    double               eta_psi,
    Rcpp::NumericVector  eta_p,
    double               eta_pos,
    Rcpp::IntegerVector  y_det,
    Rcpp::NumericVector  y_pos_vals,
    double               phi_prec,
    double               sigma_u,
    int                  n_quad   = 15,
    std::string          curvature = "observed") {
    return eval_one_latent_cell_<tulpaObs::BetaLatent>(
        eta_psi, eta_p, eta_pos, y_det, y_pos_vals, phi_prec, sigma_u,
        n_quad, parse_curvature_(curvature));
}
