// simulate_ms.cpp
// Family generators for the community (multispecies) families. Each draws ONE
// replicate at the per-species, per-site parameters handed in (occupancy /
// abundance / detection probabilities, cover predictors), from R's RNG stream
// via the R:: samplers in a fixed order. simulate_<family>() hands them the
// stated truth, simulate(fit) the fitted per-species values.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "simulate_helpers.h"
using namespace Rcpp;
using tulpaObs::clamp_eta;
using tulpaObs::draw_latent_N;

// Community N-mixture (ms_abun). lambda / p are [n_sites x n_species] fitted
// values; obs_mask flags observed (site, visit, species); size_s the per-species
// NB size (NA = Poisson); omega_s the per-species structural-zero share (NA =
// not zero-inflated). Per species: latent N (n_sites), the structural-zero draws
// when zero-inflated, then the observed visits' binomial detections, site by
// site. Returns an [n_sites x max_v x n_species] array.
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_simulate_ms_nmix(
    Rcpp::NumericMatrix lambda, Rcpp::NumericMatrix p,
    Rcpp::NumericVector size_s, Rcpp::NumericVector omega_s,
    Rcpp::IntegerVector obs_mask, int n_sites, int max_v, int n_species
) {
  Rcpp::RNGScope scope;
  const int* pm = obs_mask.begin();
  const std::size_t sp_stride = (std::size_t) n_sites * max_v;
  Rcpp::IntegerVector ys((std::size_t) n_sites * max_v * n_species);
  std::fill(ys.begin(), ys.end(), NA_INTEGER);
  int* base = ys.begin();
  Rcpp::IntegerMatrix Nlat(n_sites, n_species);
  for (int sp = 0; sp < n_species; ++sp) {
    std::vector<int> N(n_sites);
    for (int i = 0; i < n_sites; ++i) N[i] = draw_latent_N(lambda(i, sp), size_s[sp]);
    if (R_finite(omega_s[sp])) {
      for (int i = 0; i < n_sites; ++i)
        if ((int) R::rbinom(1.0, omega_s[sp]) == 1) N[i] = 0;
    }
    for (int i = 0; i < n_sites; ++i) Nlat(i, sp) = N[i];
    for (int i = 0; i < n_sites; ++i)
      for (int j = 0; j < max_v; ++j) {
        std::size_t off = (std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i;
        if (pm[off] != 0)
          base[off] = (int) R::rbinom((double) N[i], p(i, sp));
      }
  }
  ys.attr("dim") = Rcpp::IntegerVector::create(n_sites, max_v, n_species);
  ys.attr("latent") = Nlat;
  return ys;
}

// Community single-season occupancy (ms_occu). psi / p are [n_sites x n_species]
// fitted values; per species draw z ~ Bernoulli(psi) (all sites first), then the
// observed visits' detections y ~ Bernoulli(z p), site by site. Returns an
// [n_sites x max_v x n_species] array.
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_simulate_ms_occu(
    Rcpp::NumericMatrix psi, Rcpp::NumericMatrix p, Rcpp::IntegerVector valid,
    int n_sites, int max_v, int n_species
) {
  Rcpp::RNGScope scope;
  const int* pvv = valid.begin();
  const std::size_t sp_stride = (std::size_t) n_sites * max_v;
  Rcpp::IntegerVector ys((std::size_t) n_sites * max_v * n_species);
  std::fill(ys.begin(), ys.end(), NA_INTEGER);
  int* base = ys.begin();
  Rcpp::IntegerMatrix zlat(n_sites, n_species);
  for (int sp = 0; sp < n_species; ++sp) {
    std::vector<int> z(n_sites);
    for (int i = 0; i < n_sites; ++i) zlat(i, sp) = z[i] = (int) R::rbinom(1.0, psi(i, sp));
    for (int i = 0; i < n_sites; ++i)
      for (int j = 0; j < max_v; ++j) {
        std::size_t off = (std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i;
        if (pvv[off] != 0) base[off] = (int) R::rbinom(1.0, z[i] * p(i, sp));
      }
  }
  ys.attr("dim") = Rcpp::IntegerVector::create(n_sites, max_v, n_species);
  ys.attr("latent") = zlat;
  return ys;
}

// Community occupancy + cover (ms_occu_cover; the spatial variant folds its field
// into psi / ep_mat in R). psi [n_sites x n_species]; p_mat / ep_mat / valid are
// flat [n_sites x max_v x n_species]. Per species: z ~ Bernoulli(psi); at an
// occupied site the visit detections, then the cover values at its detected
// visits, drawn as .tobs_draw_positive_cover() draws them (beta held inside
// [1e-6, 1 - 1e-6]). `disp` is one value or one per species. Returns list(y,
// y_pos).
// [[Rcpp::export]]
Rcpp::List cpp_simulate_ms_occu_cover(
    Rcpp::NumericMatrix psi, Rcpp::NumericVector p_mat, Rcpp::NumericVector ep_mat,
    Rcpp::IntegerVector valid, Rcpp::NumericVector disp, int positive,
    int n_sites, int max_v, int n_species
) {
  // positive: 0 = lognormal, 3 = beta, 4 = gaussian (see .occu_cover_pos_code /
  // src/occu_coupling_shared.h).
  Rcpp::RNGScope scope;
  const int* pvv = valid.begin();
  const double* pp = p_mat.begin(); const double* pe = ep_mat.begin();
  const std::size_t sp_stride = (std::size_t) n_sites * max_v;
  const bool disp_vec = disp.size() == n_species;
  Rcpp::IntegerVector y((std::size_t) n_sites * max_v * n_species);
  Rcpp::NumericVector yp((std::size_t) n_sites * max_v * n_species);
  std::fill(y.begin(), y.end(), NA_INTEGER);
  std::fill(yp.begin(), yp.end(), NA_REAL);
  int* by = y.begin(); double* byp = yp.begin();
  Rcpp::IntegerMatrix zlat(n_sites, n_species);
  for (int sp = 0; sp < n_species; ++sp) {
    double d = disp_vec ? disp[sp] : disp[0];
    std::vector<int> z(n_sites);
    for (int i = 0; i < n_sites; ++i) zlat(i, sp) = z[i] = (int) R::rbinom(1.0, psi(i, sp));
    for (int i = 0; i < n_sites; ++i) {
      std::vector<int> vis;
      for (int j = 0; j < max_v; ++j)
        if (pvv[(std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i]) vis.push_back(j);
      if (vis.empty()) continue;
      if (z[i] == 0) { for (int j : vis) by[(std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i] = 0; continue; }
      std::vector<int> det;
      for (int j : vis) {
        std::size_t off = (std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i;
        int dd = (int) R::rbinom(1.0, pp[off]);
        by[off] = dd; if (dd == 1) det.push_back(j);
      }
      for (int j : det) {
        std::size_t off = (std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i;
        double eta = pe[off];
        double v;
        if (positive == 3) {
          double mu = tulpa::math::inv_logit(clamp_eta(eta));
          v = R::rbeta(mu * d, (1.0 - mu) * d);
          v = v < 1e-6 ? 1e-6 : (v > 1.0 - 1e-6 ? 1.0 - 1e-6 : v);
        } else if (positive == 4) {
          v = R::rnorm(eta, d);
        } else {
          v = std::exp(R::rnorm(eta, d));
        }
        byp[off] = v;
      }
    }
  }
  Rcpp::IntegerVector dim = Rcpp::IntegerVector::create(n_sites, max_v, n_species);
  y.attr("dim") = dim; yp.attr("dim") = dim;
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("y") = y, Rcpp::Named("y_pos") = yp);
  out.attr("latent") = zlat;
  return out;
}

// Community multi-season occupancy (ms_dyn_occu). psi1 / p [n_sites x n_species];
// gamma / eps [n_sites]; valid flat [n_sites x max_v x n_seasons x n_species].
// Per species: the season-1 state at every site, then interval by interval an
// extinction and a colonisation draw at every site (the site's state picks which
// applies), then the detections site by site and season by season. Returns an
// [n_sites x max_v x n_seasons x n_species] array.
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_simulate_ms_dyn_occu(
    Rcpp::NumericMatrix psi1, Rcpp::NumericMatrix p, Rcpp::NumericVector gamma,
    Rcpp::NumericVector eps, Rcpp::IntegerVector valid,
    int n_sites, int max_v, int n_seasons, int n_species
) {
  Rcpp::RNGScope scope;
  const std::size_t sp_stride = (std::size_t) n_sites * max_v * n_seasons;
  Rcpp::IntegerVector out((std::size_t) sp_stride * n_species);
  std::fill(out.begin(), out.end(), NA_INTEGER);
  const int* pvv = valid.begin();
  int* base = out.begin();
  std::vector<int> z((std::size_t) n_sites * n_seasons), ext(n_sites), col(n_sites);
  Rcpp::IntegerVector zlat((std::size_t) n_sites * n_seasons * n_species);
  for (int sp = 0; sp < n_species; ++sp) {
    for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, psi1(i, sp));
    for (int t = 1; t < n_seasons; ++t) {
      for (int i = 0; i < n_sites; ++i) ext[i] = (int) R::rbinom(1.0, eps[i]);
      for (int i = 0; i < n_sites; ++i) col[i] = (int) R::rbinom(1.0, gamma[i]);
      for (int i = 0; i < n_sites; ++i) {
        const int zp = z[(std::size_t) (t - 1) * n_sites + i];
        z[(std::size_t) t * n_sites + i] = zp * (1 - ext[i]) + (1 - zp) * col[i];
      }
    }
    std::copy(z.begin(), z.end(), zlat.begin() + (std::size_t) sp * n_sites * n_seasons);
    for (int i = 0; i < n_sites; ++i)
      for (int t = 0; t < n_seasons; ++t)
        for (int j = 0; j < max_v; ++j) {
          std::size_t off = (std::size_t) sp * sp_stride +
            (std::size_t) t * n_sites * max_v + (std::size_t) j * n_sites + i;
          if (pvv[off] != 0)
            base[off] = (int) R::rbinom(1.0, z[(std::size_t) t * n_sites + i] * p(i, sp));
        }
  }
  out.attr("dim") = Rcpp::IntegerVector::create(n_sites, max_v, n_seasons, n_species);
  zlat.attr("dim") = Rcpp::IntegerVector::create(n_sites, n_seasons, n_species);
  out.attr("latent") = zlat;
  return out;
}

// Community integrated multi-source occupancy (ms_int_occu). psi [n_sites x
// n_species]; pd a list of D matrices [n_sites x n_species]; valid a list of D
// flat arrays [n_sites x J_d x n_species]. Per species: z ~ Bernoulli(psi); per
// source, the observed visits' detections at occupied sites, site by site.
// Returns a list of D arrays [n_sites x J_d x n_species].
// [[Rcpp::export]]
Rcpp::List cpp_simulate_ms_int_occu(
    Rcpp::NumericMatrix psi, Rcpp::List pd, Rcpp::List valid,
    Rcpp::IntegerVector J_d, int n_sites, int n_species, int D
) {
  Rcpp::RNGScope scope;
  Rcpp::List srcs(D);
  std::vector<Rcpp::IntegerVector> arr(D);
  for (int d = 0; d < D; ++d) {
    arr[d] = Rcpp::IntegerVector((std::size_t) n_sites * J_d[d] * n_species);
    std::fill(arr[d].begin(), arr[d].end(), NA_INTEGER);
  }
  Rcpp::IntegerMatrix zlat(n_sites, n_species);
  for (int sp = 0; sp < n_species; ++sp) {
    std::vector<int> z(n_sites);
    for (int i = 0; i < n_sites; ++i) zlat(i, sp) = z[i] = (int) R::rbinom(1.0, psi(i, sp));
    for (int d = 0; d < D; ++d) {
      Rcpp::NumericMatrix pdd = pd[d];
      Rcpp::IntegerVector vd = valid[d];
      int Jd = J_d[d];
      std::size_t sp_stride = (std::size_t) n_sites * Jd;
      int* base = arr[d].begin();
      const int* pvd = vd.begin();
      for (int i = 0; i < n_sites; ++i)
        for (int j = 0; j < Jd; ++j) {
          std::size_t off = (std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i;
          if (pvd[off] != 0)
            base[off] = (z[i] == 0) ? 0 : (int) R::rbinom(1.0, pdd(i, sp));
        }
    }
  }
  for (int d = 0; d < D; ++d) {
    arr[d].attr("dim") = Rcpp::IntegerVector::create(n_sites, J_d[d], n_species);
    srcs[d] = arr[d];
  }
  srcs.attr("latent") = zlat;
  return srcs;
}
