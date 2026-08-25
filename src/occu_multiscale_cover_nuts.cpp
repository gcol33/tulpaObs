// occu_multiscale_cover_nuts.cpp
// NUTS target for the non-spatial three-level occupancy + cover hurdle
// (occu_multiscale_cover()). The flat parameter vector is
//   theta = (beta_psi, beta_theta, beta_p[site, visit], beta_pos[site, visit],
//            log_disp)
// and the joint log-posterior is the exact three-level marginal (z over cells
// and a over plots both summed in closed form -- the same density the
// non-spatial Laplace path .occu_mscale_cover_nonspatial_ll() optimises) plus weak
// Gaussian coefficient priors. The cover hurdle, the no-detection occupancy /
// availability mixtures, and their eta-derivatives reuse the cover policies and
// nodet_mixture_block from occu_coupling_shared.h (single source of truth with
// the cell-coupling spec) -- no new likelihood math here.
//
// The shared engine (nuts_engine.h) drives tulpa's NUTS. The R reference
// .tobs_occu_mscale_cover_nuts_logpost (R/occu_multiscale_cover_nuts.R) is the
// oracle that cpp_occu_mscale_cover_nuts_joint_logpost is cross-checked against.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "occu_coupling_shared.h"   // sigmoid_ / Pos policies / nodet_mixture_block
#include "tobs_shape.h"
#include "nuts_engine.h"

namespace tulpaObs {

// Per-cell, per-plot, per-visit data + designs for the three-level marginal.
// The eta of each arm is a design-times-coefficient sum; the site-level and
// optional visit-level blocks of the p / pos arms are stored separately so the
// flat coefficient gradient is the design-sandwiched eta-gradient (mirroring the
// R oracle's two-block .occu_ms_eta_visit()).
struct MscaleCoverNutsModel {
    int n_cells = 0, n_plots = 0, max_visits = 0;
    int p_psi = 0, p_theta = 0, p_p_site = 0, p_p_visit = 0,
        p_pos_site = 0, p_pos_visit = 0;
    int o_psi = 0, o_theta = 0, o_p = 0, o_pos = 0, o_disp = 0, total = 0;
    int positive = 0;   // 0 = lognormal, 3 = beta, 4 = gaussian
    double sigma_beta = 5.0;
    double sigma_logdisp = 5.0;

    Rcpp::NumericMatrix X_psi;        // n_cells   x p_psi
    Rcpp::NumericMatrix X_theta;      // n_plots   x p_theta
    Rcpp::NumericMatrix X_p_site;     // n_plots   x p_p_site
    Rcpp::NumericMatrix X_pos_site;   // n_plots   x p_pos_site
    Rcpp::NumericMatrix X_p_visit;    // (n_plots*J) x p_p_visit (visit-major rows), may be empty
    Rcpp::NumericMatrix X_pos_visit;  // (n_plots*J) x p_pos_visit, may be empty
    std::vector<int>    y;            // n_plots*J, plot-major (row p*J + v)
    std::vector<double> y_pos;        // n_plots*J
    std::vector<char>   valid;        // n_plots*J
    std::vector<int>    plot_cell;    // 0-based cell per plot
    std::vector<std::vector<int>> plots_of_cell;   // plot indices per cell
};

inline MscaleCoverNutsModel mscale_cover_nuts_build(const Rcpp::List& spec) {
    MscaleCoverNutsModel m;
    m.X_psi      = Rcpp::as<Rcpp::NumericMatrix>(spec["X_psi"]);
    m.X_theta    = Rcpp::as<Rcpp::NumericMatrix>(spec["X_theta"]);
    m.X_p_site   = Rcpp::as<Rcpp::NumericMatrix>(spec["X_p_site"]);
    m.X_pos_site = Rcpp::as<Rcpp::NumericMatrix>(spec["X_pos_site"]);
    m.n_cells    = Rcpp::as<int>(spec["n_cells"]);
    m.max_visits = Rcpp::as<int>(spec["max_visits"]);
    m.positive   = Rcpp::as<int>(spec["positive"]);

    m.n_plots    = m.X_theta.nrow();
    m.p_psi      = m.X_psi.ncol();
    m.p_theta    = m.X_theta.ncol();
    m.p_p_site   = m.X_p_site.ncol();
    m.p_pos_site = m.X_pos_site.ncol();

    if (spec.containsElementNamed("X_p_visit") &&
        !Rf_isNull(spec["X_p_visit"])) {
        m.X_p_visit = Rcpp::as<Rcpp::NumericMatrix>(spec["X_p_visit"]);
        m.p_p_visit = m.X_p_visit.ncol();
    } else { m.p_p_visit = 0; }
    if (spec.containsElementNamed("X_pos_visit") &&
        !Rf_isNull(spec["X_pos_visit"])) {
        m.X_pos_visit = Rcpp::as<Rcpp::NumericMatrix>(spec["X_pos_visit"]);
        m.p_pos_visit = m.X_pos_visit.ncol();
    } else { m.p_pos_visit = 0; }

    const int p_p   = m.p_p_site   + m.p_p_visit;
    const int p_pos = m.p_pos_site + m.p_pos_visit;
    m.o_psi   = 0;
    m.o_theta = m.p_psi;
    m.o_p     = m.o_theta + m.p_theta;
    m.o_pos   = m.o_p + p_p;
    m.o_disp  = m.o_pos + p_pos;
    m.total   = m.o_disp + 1;

    Rcpp::IntegerVector yv  = spec["y"];        // n_plots*J, plot-major
    Rcpp::NumericVector yp  = spec["y_pos"];
    Rcpp::LogicalVector vd  = spec["valid"];
    Rcpp::IntegerVector pc  = spec["plot_cell"];
    const int N = m.n_plots * m.max_visits;
    if ((int) yv.size() != N) Rcpp::stop("y length != n_plots * max_visits");
    m.y.assign(yv.begin(), yv.end());
    m.y_pos.assign(yp.begin(), yp.end());
    m.valid.assign(N, 0);
    for (int i = 0; i < N; ++i) m.valid[i] = (vd[i] == TRUE) ? 1 : 0;
    m.plot_cell.resize(m.n_plots);
    for (int i = 0; i < m.n_plots; ++i) m.plot_cell[i] = pc[i] - 1;
    m.plots_of_cell.assign(m.n_cells, std::vector<int>());
    for (int i = 0; i < m.n_plots; ++i) {
        const int c = m.plot_cell[i];
        if (c < 0 || c >= m.n_cells) Rcpp::stop("plot_cell out of range");
        m.plots_of_cell[c].push_back(i);
    }
    return m;
}

// Build the [n_plots x J] eta matrices (site predictor broadcast across a
// plot's visits + optional visit-varying part, visit-major rows). Returns flat
// length n_plots*J, plot-major (row p*J + v), matching `y` / `y_pos` / `valid`.
inline void ms_cover_eta_visit(const MscaleCoverNutsModel& m, const double* theta,
                               int o, int p_site, int p_visit,
                               const Rcpp::NumericMatrix& Xs,
                               const Rcpp::NumericMatrix& Xv,
                               std::vector<double>& eta) {
    const int J = m.max_visits, np = m.n_plots;
    eta.assign((std::size_t) np * J, 0.0);
    for (int i = 0; i < np; ++i) {
        double e_site = 0.0;
        for (int k = 0; k < p_site; ++k) e_site += Xs(i, k) * theta[o + k];
        for (int v = 0; v < J; ++v) eta[(std::size_t) i * J + v] = e_site;
    }
    if (p_visit > 0) {
        // X_*_visit rows are ordered (plot-1)*J + visit (visit-major), the same
        // flat plot-major index used here.
        for (int i = 0; i < np; ++i)
            for (int v = 0; v < J; ++v) {
                const int row = i * J + v;
                double e_v = 0.0;
                for (int k = 0; k < p_visit; ++k)
                    e_v += Xv(row, k) * theta[o + p_site + k];
                eta[(std::size_t) row] += e_v;
            }
    }
}

// Joint log-posterior + full coefficient gradient of the three-level marginal.
inline double mscale_cover_nuts_eval(const MscaleCoverNutsModel& m, const double* theta,
                                 double* grad) {
    const int J = m.max_visits, np = m.n_plots;
    const int p_p = m.p_p_site + m.p_p_visit;
    for (int j = 0; j < m.total; ++j) grad[j] = 0.0;
    const double disp = std::exp(theta[m.o_disp]);

    // eta_psi (n_cells), eta_theta (n_plots), eta_p / eta_pos (n_plots*J flat).
    std::vector<double> eta_psi(m.n_cells, 0.0), eta_theta(np, 0.0);
    for (int c = 0; c < m.n_cells; ++c) {
        double e = 0.0;
        for (int k = 0; k < m.p_psi; ++k) e += m.X_psi(c, k) * theta[m.o_psi + k];
        eta_psi[c] = e;
    }
    for (int i = 0; i < np; ++i) {
        double e = 0.0;
        for (int k = 0; k < m.p_theta; ++k) e += m.X_theta(i, k) * theta[m.o_theta + k];
        eta_theta[i] = e;
    }
    std::vector<double> eta_p, eta_pos;
    ms_cover_eta_visit(m, theta, m.o_p,   m.p_p_site,   m.p_p_visit,
                       m.X_p_site,   m.X_p_visit,   eta_p);
    ms_cover_eta_visit(m, theta, m.o_pos, m.p_pos_site, m.p_pos_visit,
                       m.X_pos_site, m.X_pos_visit, eta_pos);

    // Eta-gradient accumulators (one per arm element); sandwiched at the end.
    std::vector<double> ge_psi(m.n_cells, 0.0), ge_theta(np, 0.0);
    std::vector<double> ge_p((std::size_t) np * J, 0.0),
                        ge_pos((std::size_t) np * J, 0.0);
    double ge_logdisp = 0.0;
    double lp = 0.0;

    // Reused per plot: the plot's valid visit rows in the padded grid.
    std::vector<int> det_rows;
    det_rows.reserve(J);

    for (int c = 0; c < m.n_cells; ++c) {
        const std::vector<int>& plots = m.plots_of_cell[c];
        const int M_c = (int) plots.size();
        const double psi = sigmoid_(eta_psi[c]);

        // Per-plot any-detection flag + cell-wide flag.
        bool any_det_cell = false;
        std::vector<char> plot_det(M_c, 0);
        for (int jj = 0; jj < M_c; ++jj) {
            const int pi = plots[jj];
            for (int v = 0; v < J; ++v) {
                const int row = pi * J + v;
                if (m.valid[row] && m.y[row] == 1) { plot_det[jj] = 1; any_det_cell = true; }
            }
        }

        if (any_det_cell) {
            // Branch A: z = 1 forced. log L = log psi + sum_j log Q_j.
            lp += log_safe(psi);
            ge_psi[c] += 1.0 - psi;
            for (int jj = 0; jj < M_c; ++jj) {
                const int pi = plots[jj];
                const double theta_j = sigmoid_(eta_theta[pi]);
                if (plot_det[jj]) {
                    // Detected plot: factorised. log theta + per-visit GLM
                    // terms, from the shared block the cell-coupling spec also
                    // evaluates; this layout gathers the plot's valid visits
                    // out of the padded grid and scatters the scores back.
                    det_rows.clear();
                    for (int v = 0; v < J; ++v) {
                        const int row = pi * J + v;
                        if (m.valid[row]) det_rows.push_back(row);
                    }
                    const int nv = (int) det_rows.size();
                    double g_th = 0.0, nh_th = 0.0;
                    lp += mscale_det_plot_block(
                        PosCodeAccess(m.positive), eta_theta[pi], nv,
                        [&](int t) { return eta_p[det_rows[t]]; },
                        [&](int t) { return eta_pos[det_rows[t]]; },
                        [&](int t) { return (double) m.y[det_rows[t]]; },
                        [&](int t) { return m.y_pos[det_rows[t]]; },
                        disp, /*want_hess=*/false, /*want_logdisp=*/true,
                        g_th, nh_th,
                        [&](int t, double g, double) { ge_p[det_rows[t]] += g; },
                        [&](int t, double g, double) { ge_pos[det_rows[t]] += g; },
                        ge_logdisp);
                    ge_theta[pi] += g_th;
                } else {
                    // Non-detected plot inside an occupied cell: within-plot
                    // availability mixture m_j = theta_j P0_j + (1 - theta_j).
                    // Reuse the shared nodet block over the plot's valid visits.
                    std::vector<int> vrows;
                    for (int v = 0; v < J; ++v) {
                        const int row = pi * J + v;
                        if (m.valid[row]) vrows.push_back(row);
                    }
                    const int nv = (int) vrows.size();
                    std::vector<double> eta_p_buf(nv), g_p_buf(nv), nh_p_buf(nv);
                    for (int t = 0; t < nv; ++t) eta_p_buf[t] = eta_p[vrows[t]];
                    double g_th = 0.0, nh_th = 0.0;
                    lp += nodet_mixture_block(
                        theta_j, eta_p_buf.data(), nv, /*want_hess=*/false,
                        /*expected=*/false, g_th, nh_th,
                        g_p_buf.data(), nh_p_buf.data(), nullptr, nullptr);
                    ge_theta[pi] += g_th;
                    for (int t = 0; t < nv; ++t) ge_p[vrows[t]] += g_p_buf[t];
                }
            }
        } else {
            // Branch B: no detection anywhere in the cell. The mixture and
            // its scores are the same derivation the cell-coupling spec
            // evaluates, so both read it from occu_coupling_shared.h; this
            // layout gathers its padded [n_plots x J] grid into the CSR the
            // block takes, and scatters the scores back through the mask.
            std::vector<int> visit_off(M_c), visit_size(M_c), rows;
            rows.reserve((std::size_t) M_c * J);
            for (int jj = 0; jj < M_c; ++jj) {
                const int pi = plots[jj];
                visit_off[jj] = (int) rows.size();
                for (int v = 0; v < J; ++v) {
                    const int row = pi * J + v;
                    if (m.valid[row]) rows.push_back(row);
                }
                visit_size[jj] = (int) rows.size() - visit_off[jj];
            }

            const MscaleNoDetCell cell = mscale_nodet_cell(
                psi, M_c,
                [&](int jj) { return eta_theta[plots[jj]]; },
                visit_off.data(), visit_size.data(),
                [&](int t)  { return eta_p[rows[t]]; });

            lp += cell.log_lik;
            ge_psi[c] += cell.s_psi;
            for (int jj = 0; jj < M_c; ++jj)
                ge_theta[plots[jj]] += cell.s_theta[jj];
            for (int t = 0; t < cell.n_visits; ++t)
                ge_p[rows[t]] += cell.s_p[t];
        }
    }

    // Sandwich eta-gradients through the designs into coefficient gradients.
    for (int c = 0; c < m.n_cells; ++c)
        for (int k = 0; k < m.p_psi; ++k)
            grad[m.o_psi + k] += ge_psi[c] * m.X_psi(c, k);
    for (int i = 0; i < np; ++i)
        for (int k = 0; k < m.p_theta; ++k)
            grad[m.o_theta + k] += ge_theta[i] * m.X_theta(i, k);
    for (int i = 0; i < np; ++i) {
        double s_p = 0.0, s_pos = 0.0;
        for (int v = 0; v < J; ++v) { s_p += ge_p[i * J + v]; s_pos += ge_pos[i * J + v]; }
        for (int k = 0; k < m.p_p_site; ++k)   grad[m.o_p + k]   += s_p   * m.X_p_site(i, k);
        for (int k = 0; k < m.p_pos_site; ++k) grad[m.o_pos + k] += s_pos * m.X_pos_site(i, k);
    }
    if (m.p_p_visit > 0)
        for (int i = 0; i < np; ++i)
            for (int v = 0; v < J; ++v) {
                const int row = i * J + v;
                for (int k = 0; k < m.p_p_visit; ++k)
                    grad[m.o_p + m.p_p_site + k] += ge_p[row] * m.X_p_visit(row, k);
            }
    if (m.p_pos_visit > 0)
        for (int i = 0; i < np; ++i)
            for (int v = 0; v < J; ++v) {
                const int row = i * J + v;
                for (int k = 0; k < m.p_pos_visit; ++k)
                    grad[m.o_pos + m.p_pos_site + k] += ge_pos[row] * m.X_pos_visit(row, k);
            }
    grad[m.o_disp] += ge_logdisp;

    // Weak Gaussian priors: N(0, sigma_beta^2) on every coefficient (matching
    // the non-spatial Laplace path's default) and a broad N(0, sigma_logdisp^2)
    // on log_disp to keep the dispersion proper -- the two sibling NUTS targets
    // (occu_cover_nuts.cpp, cover_nuts.cpp) both sample log_disp under this
    // prior; this target used to leave it flat.
    const double ib2 = 1.0 / (m.sigma_beta * m.sigma_beta);
    const int n_beta = m.o_disp;   // all coords before log_disp
    for (int k = 0; k < n_beta; ++k) {
        lp      -= 0.5 * ib2 * theta[k] * theta[k];
        grad[k] -= ib2 * theta[k];
    }
    const double ild2 = 1.0 / (m.sigma_logdisp * m.sigma_logdisp);
    lp             -= 0.5 * ild2 * theta[m.o_disp] * theta[m.o_disp];
    grad[m.o_disp] -= ild2 * theta[m.o_disp];
    return lp;
}

inline void mscale_cover_nuts_full_grad(const std::vector<double>& params,
                                    const tulpa::ModelData& data,
                                    const tulpa::ParamLayout& /*layout*/,
                                    std::vector<double>& grad, double* log_post_out) {
    const MscaleCoverNutsModel* m =
        static_cast<const MscaleCoverNutsModel*>(data.model_response_data);
    grad.assign((std::size_t) m->total, 0.0);
    const double lp = mscale_cover_nuts_eval(*m, params.data(), grad.data());
    if (log_post_out) *log_post_out = lp;
}

}  // namespace tulpaObs

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List cpp_occu_mscale_cover_nuts_joint_logpost(Rcpp::List spec,
                                                Rcpp::NumericVector theta,
                                                double sigma_beta,
                                                double sigma_logdisp) {
    tulpaObs::MscaleCoverNutsModel m = tulpaObs::mscale_cover_nuts_build(spec);
    m.sigma_beta = sigma_beta;
    m.sigma_logdisp = sigma_logdisp;
    if ((int) theta.size() != m.total)
        Rcpp::stop("theta length %d != expected %d", (int) theta.size(), m.total);
    Rcpp::NumericVector grad(m.total);
    const double lp = tulpaObs::mscale_cover_nuts_eval(m, theta.begin(), grad.begin());
    return Rcpp::List::create(Rcpp::Named("lp") = lp, Rcpp::Named("grad") = grad);
}

// [[Rcpp::export]]
Rcpp::List cpp_occu_mscale_cover_nuts(Rcpp::List spec, Rcpp::NumericVector theta0,
                                  double sigma_beta, double sigma_logdisp,
                                  Rcpp::Nullable<Rcpp::NumericVector> inv_metric,
                                  int n_iter, int n_warmup, int max_treedepth,
                                  double adapt_delta, int seed, bool verbose) {
    tulpaObs::MscaleCoverNutsModel m = tulpaObs::mscale_cover_nuts_build(spec);
    m.sigma_beta = sigma_beta;
    m.sigma_logdisp = sigma_logdisp;
    return tulpaObs::run_tulpa_nuts(&tulpaObs::mscale_cover_nuts_full_grad, &m, m.total,
                                    theta0, sigma_beta, tulpaObs::shape::optional_numeric(inv_metric.get(), "inv_metric"), n_iter, n_warmup,
                                    max_treedepth, adapt_delta, seed, verbose);
}
