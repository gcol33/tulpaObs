// occu_mscale_cover_ploglik.cpp
// Batched per-cell pointwise log-likelihood for the three-level multiscale
// occupancy + cover family (occu_multiscale_cover, non-spatial). The R reference
// .occu_mscale_cover_nonspatial_ll (R/occu_multiscale_cover.R, per_cell = TRUE) is
// the oracle; this port mirrors it draw for draw and parallelises over draws.
// The marginal is cell (occupancy z) over plot (availability a) over visit
// (detection), plus the per-detected-visit cover density. Designs arrive as
// matrices; the visit blocks are site-major [n_plots x J] (X_*_visit has one row
// per (plot, visit)). Byte-close to the R marginal (~1e-13).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "tobs_shape.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using tulpaObs::stable_plogis;

namespace {
inline double clp(double x) {
  if (x < 1e-12) return 1e-12;
  if (x > 1.0 - 1e-12) return 1.0 - 1e-12;
  return x;
}
}  // namespace

// idx_* are 0-based coefficient offsets into each draw row; p_* the block widths.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_occu_mscale_cover_ploglik(
    Rcpp::NumericMatrix draws,        // [S x total]
    Rcpp::NumericMatrix X_psi,        // [n_cells x p_psi]
    Rcpp::NumericMatrix X_theta,      // [n_plots x p_theta]
    Rcpp::NumericMatrix X_p_site,     // [n_plots x p_p_site]
    Rcpp::NumericMatrix X_p_visit,    // [n_plots*J x p_p_visit] (0 cols if none)
    Rcpp::NumericMatrix X_pos_site,   // [n_plots x p_pos_site]
    Rcpp::NumericMatrix X_pos_visit,  // [n_plots*J x p_pos_visit] (0 cols if none)
    Rcpp::IntegerMatrix y,            // [n_plots x J], 1 = detected, else 0/NA
    Rcpp::NumericMatrix y_pos,        // [n_plots x J] cover at detected visits
    Rcpp::LogicalMatrix valid,        // [n_plots x J]
    Rcpp::IntegerVector plot_cell,    // [n_plots], 1-based cell of each plot
    int positive,                     // 0 = lognormal, 3 = beta, 4 = gaussian

    int idx_psi, int p_psi, int idx_theta, int p_theta,
    int idx_p_site, int p_p_site, int idx_p_visit, int p_p_visit,
    int idx_pos_site, int p_pos_site, int idx_pos_visit, int p_pos_visit,
    int idx_disp, int n_threads
) {
  const int S = draws.nrow();
  const int n_plots = X_theta.nrow();
  const int n_cells = X_psi.nrow();
  const int J = y.ncol();
  const bool has_pv = p_p_visit > 0, has_posv = p_pos_visit > 0;

  // Every design row count and every coefficient offset comes from a different
  // argument than the buffer it indexes: the block widths and offsets are
  // packed in R against `draws`, the visit blocks are site-major [n_plots x J],
  // and plot_cell addresses the per-cell accumulator. None of that is checked
  // downstream -- the sibling multiscale NUTS builder checks its cell index,
  // this kernel checked nothing -- so it is related here, before any .begin().
  namespace sh = tulpaObs::shape;
  const R_xlen_t total = draws.ncol();
  sh::check_block(idx_psi, p_psi, total, "psi");
  sh::check_block(idx_theta, p_theta, total, "theta");
  sh::check_block(idx_p_site, p_p_site, total, "p_site");
  sh::check_block(idx_p_visit, p_p_visit, total, "p_visit");
  sh::check_block(idx_pos_site, p_pos_site, total, "pos_site");
  sh::check_block(idx_pos_visit, p_pos_visit, total, "pos_visit");
  sh::check_block(idx_disp, 1, total, "disp");

  sh::check_dim(X_psi, n_cells, p_psi, "X_psi");
  sh::check_dim(X_theta, n_plots, p_theta, "X_theta");
  sh::check_dim(X_p_site, n_plots, p_p_site, "X_p_site");
  sh::check_dim(X_pos_site, n_plots, p_pos_site, "X_pos_site");
  if (has_pv)
    sh::check_dim(X_p_visit, (R_xlen_t) n_plots * J, p_p_visit, "X_p_visit");
  if (has_posv)
    sh::check_dim(X_pos_visit, (R_xlen_t) n_plots * J, p_pos_visit,
                  "X_pos_visit");
  sh::check_dim(y, n_plots, J, "y");
  sh::check_dim(y_pos, n_plots, J, "y_pos");
  sh::check_dim(valid, n_plots, J, "valid");
  sh::check_len(plot_cell, n_plots, "plot_cell");
  sh::check_index1(plot_cell, n_cells, "plot_cell");

  Rcpp::NumericMatrix ll(S, n_cells);
  const double* pd  = draws.begin();
  const int*    pc  = plot_cell.begin();
  double* pll = ll.begin();

#ifdef _OPENMP
  #pragma omp parallel num_threads(n_threads > 0 ? n_threads : 1)
#endif
  {
    std::vector<double> sum_logpj(n_cells);
    std::vector<unsigned char> det_cell(n_cells);
#ifdef _OPENMP
    #pragma omp for schedule(static)
#endif
    for (int d = 0; d < S; ++d) {
      // draws is column-major [S x p]: coefficient (off + k) for draw d.
      auto coef = [&](int off, int k) { return pd[(std::size_t) (off + k) * S + d]; };

      double log_disp = coef(idx_disp, 0);
      double disp = std::exp(log_disp);

      for (int c = 0; c < n_cells; ++c) { sum_logpj[c] = 0.0; det_cell[c] = 0; }

      // Per-plot three-level term, aggregated into its cell.
      for (int i = 0; i < n_plots; ++i) {
        // theta (availability) predictor for this plot.
        double eta_theta = 0.0;
        for (int k = 0; k < p_theta; ++k)
          eta_theta += X_theta((std::size_t) i, k) * coef(idx_theta, k);
        double theta = clp(stable_plogis(eta_theta));

        // Site-level detection / cover predictors (broadcast across visits).
        double eta_p_site = 0.0, eta_pos_site = 0.0;
        for (int k = 0; k < p_p_site; ++k)
          eta_p_site += X_p_site((std::size_t) i, k) * coef(idx_p_site, k);
        for (int k = 0; k < p_pos_site; ++k)
          eta_pos_site += X_pos_site((std::size_t) i, k) * coef(idx_pos_site, k);

        double sum_hdet = 0.0, sum_1mp = 0.0, sum_cover = 0.0;
        bool det_plot = false;
        for (int j = 0; j < J; ++j) {
          if (!valid((std::size_t) i, j)) continue;
          int vrow = i * J + j;                       // site-major visit row
          double eta_p = eta_p_site;
          if (has_pv) {
            for (int k = 0; k < p_p_visit; ++k)
              eta_p += X_p_visit((std::size_t) vrow, k) * coef(idx_p_visit, k);
          }
          double p = clp(stable_plogis(eta_p));
          double lp = std::log(p), l1mp = std::log(1.0 - p);
          sum_1mp += l1mp;
          int yij = y((std::size_t) i, j);
          if (yij == 1) {
            det_plot = true;
            sum_hdet += lp;
            // Cover term at a detected plot with an observed cover; a missing
            // (NA -> non-finite) cover drops out (missing-at-random cover), the
            // detection mixture above still counts it. Without the guard an NA
            // cover poisons sum_cover with NaN, unlike the 2-level sibling
            // occu_cover_ploglik.cpp
            // which already guards it.
            double cv = y_pos((std::size_t) i, j);
            if (std::isfinite(cv)) {
              double eta_pos = eta_pos_site;
              if (has_posv) {
                double add = 0.0;
                for (int k = 0; k < p_pos_visit; ++k)
                  add += X_pos_visit((std::size_t) vrow, k) * coef(idx_pos_visit, k);
                eta_pos += add;
              }
              double dens;
              if (positive == 3) {
                double mu = clp(stable_plogis(eta_pos));
                double cvc = cv < 1e-9 ? 1e-9 : (cv > 1.0 - 1e-9 ? 1.0 - 1e-9 : cv);
                double a = mu * disp, b = (1.0 - mu) * disp;
                dens = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                       (a - 1.0) * std::log(cvc) + (b - 1.0) * std::log(1.0 - cvc);
              } else if (positive == 4) {
                // identity-Gaussian: mu = eta, residual on the raw response.
                double r = (cv - eta_pos) / disp;
                dens = -0.5 * std::log(2.0 * M_PI) - std::log(disp) - 0.5 * r * r;
              } else {
                double ly = std::log(cv), r = (ly - eta_pos) / disp;
                dens = -0.5 * std::log(2.0 * M_PI) - std::log(disp) - 0.5 * r * r - ly;
              }
              sum_cover += dens;
            }
          } else {
            sum_hdet += l1mp;
          }
        }
        // Plot log-prob given z = 1 (availability a marginalised when no detection).
        double log_theta = std::log(theta), log_1mtheta = std::log(1.0 - theta);
        double log_pj;
        if (det_plot) {
          log_pj = log_theta + sum_hdet + sum_cover;
        } else {
          log_pj = tulpaObs::logsumexp2(log_theta + sum_1mp, log_1mtheta);
        }
        int cell = pc[i] - 1;
        sum_logpj[cell] += log_pj;
        if (det_plot) det_cell[cell] = 1;
      }

      // Cell-level occupancy marginal.
      for (int c = 0; c < n_cells; ++c) {
        double eta_psi = 0.0;
        for (int k = 0; k < p_psi; ++k)
          eta_psi += X_psi((std::size_t) c, k) * coef(idx_psi, k);
        double psi = clp(stable_plogis(eta_psi));
        double log_psi = std::log(psi), log_1mpsi = std::log(1.0 - psi);
        double val;
        if (det_cell[c]) {
          val = log_psi + sum_logpj[c];
        } else {
          val = tulpaObs::logsumexp2(log_psi + sum_logpj[c], log_1mpsi);
        }
        pll[(std::size_t) c * S + d] = val;
      }
    }
  }
  return ll;
}
