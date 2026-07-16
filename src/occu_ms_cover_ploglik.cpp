// occu_ms_cover_ploglik.cpp
// Batched per-cell pointwise log-likelihood for the three-level multiscale
// occupancy + cover family (occu_multiscale_cover, non-spatial). The R reference
// .occu_ms_cover_nonspatial_ll (R/occu_multiscale_cover.R, per_cell = TRUE) is
// the oracle; this port mirrors it draw for draw and parallelises over draws.
// The marginal is cell (occupancy z) over plot (availability a) over visit
// (detection), plus the per-detected-visit cover density. Designs arrive as
// matrices; the visit blocks are site-major [n_plots x J] (X_*_visit has one row
// per (plot, visit)). Byte-close to the R marginal (~1e-13).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace {
inline double clp(double x) {
  if (x < 1e-12) return 1e-12;
  if (x > 1.0 - 1e-12) return 1.0 - 1e-12;
  return x;
}
inline double plogis(double x) {
  if (x >= 0.0) { double z = std::exp(-x); return 1.0 / (1.0 + z); }
  double z = std::exp(x); return z / (1.0 + z);
}
// log(exp(a)+exp(b)) as the R code writes it: m + log(exp(a-m)+exp(b-m)).
inline double lae2(double a, double b) {
  double m = a > b ? a : b;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}
}  // namespace

// idx_* are 0-based coefficient offsets into each draw row; p_* the block widths.
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_occu_ms_cover_ploglik(
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
        double theta = clp(plogis(eta_theta));

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
          double p = clp(plogis(eta_p));
          double lp = std::log(p), l1mp = std::log(1.0 - p);
          sum_1mp += l1mp;
          int yij = y((std::size_t) i, j);
          if (yij == 1) {
            det_plot = true;
            sum_hdet += lp;
            double eta_pos = eta_pos_site;
            if (has_posv) {
              double add = 0.0;
              for (int k = 0; k < p_pos_visit; ++k)
                add += X_pos_visit((std::size_t) vrow, k) * coef(idx_pos_visit, k);
              eta_pos += add;
            }
            double cv = y_pos((std::size_t) i, j);
            double dens;
            if (positive == 3) {
              double mu = clp(plogis(eta_pos));
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
          log_pj = lae2(log_theta + sum_1mp, log_1mtheta);
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
        double psi = clp(plogis(eta_psi));
        double log_psi = std::log(psi), log_1mpsi = std::log(1.0 - psi);
        double val;
        if (det_cell[c]) {
          val = log_psi + sum_logpj[c];
        } else {
          val = lae2(log_psi + sum_logpj[c], log_1mpsi);
        }
        pll[(std::size_t) c * S + d] = val;
      }
    }
  }
  return ll;
}
