// simulate_ms.cpp
// C++ generators for the community (multispecies) simulate() methods. Unlike the
// single-species simulators, these draw from the POSTERIOR-MEAN fitted values
// (computed in R, deterministic, no draw selection), so there is no sample.int
// here -- only the RNG data generation, run from R's RNG stream via the R::
// samplers in the SAME order as the former R loops (byte-identical under a seed).

#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

namespace {
const double EB = 30.0;
inline double clampe(double e) { return e > EB ? EB : (e < -EB ? -EB : e); }
inline double plg(double x) {
  if (x >= 0.0) { double z = std::exp(-x); return 1.0 / (1.0 + z); }
  double z = std::exp(x); return z / (1.0 + z);
}
inline int draw_N(double lambda, double size) {   // NB(mu) or Poisson
  if (R_finite(size)) {
    if (lambda <= 0.0) return 0;
    return (int) R::rpois(R::rgamma(size, lambda / size));
  }
  return (int) R::rpois(lambda);
}
}  // namespace

// Community N-mixture (ms_abun). lambda / p are [n_sites x n_species] fitted
// values; obs_mask flags observed (site, visit, species); size_s the per-species
// NB size (NA = Poisson). Per species: latent N (n_sites), then the observed
// visits' binomial detections. Returns a list of nsim [n_sites x max_v x
// n_species] arrays.
// [[Rcpp::export]]
Rcpp::List cpp_simulate_ms_nmix(
    Rcpp::NumericMatrix lambda, Rcpp::NumericMatrix p,
    Rcpp::NumericVector size_s, Rcpp::IntegerVector obs_mask,
    int n_sites, int max_v, int n_species, int nsim
) {
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  const int* pm = obs_mask.begin();
  const std::size_t sp_stride = (std::size_t) n_sites * max_v;
  for (int s = 0; s < nsim; ++s) {
    Rcpp::IntegerVector ys((std::size_t) n_sites * max_v * n_species);
    std::fill(ys.begin(), ys.end(), NA_INTEGER);
    int* base = ys.begin();
    for (int sp = 0; sp < n_species; ++sp) {
      std::vector<int> N(n_sites);
      for (int i = 0; i < n_sites; ++i) N[i] = draw_N(lambda(i, sp), size_s[sp]);
      for (int i = 0; i < n_sites; ++i)
        for (int j = 0; j < max_v; ++j) {
          std::size_t off = (std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i;
          if (pm[off] != 0)
            base[off] = (int) R::rbinom((double) N[i], p(i, sp));
        }
    }
    ys.attr("dim") = Rcpp::IntegerVector::create(n_sites, max_v, n_species);
    out[s] = ys;
  }
  return out;
}

// Community single-season occupancy (ms_occu). psi / p are [n_sites x n_species]
// fitted values; per species draw z ~ Bernoulli(psi) (all sites first), then the
// observed visits' detections y ~ Bernoulli(z p). RNG order matches the R loop.
// [[Rcpp::export]]
Rcpp::List cpp_simulate_ms_occu(
    Rcpp::NumericMatrix psi, Rcpp::NumericMatrix p, Rcpp::IntegerVector valid,
    int n_sites, int max_v, int n_species, int nsim
) {
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  const int* pvv = valid.begin();
  const std::size_t sp_stride = (std::size_t) n_sites * max_v;
  for (int s = 0; s < nsim; ++s) {
    Rcpp::IntegerVector ys((std::size_t) n_sites * max_v * n_species);
    std::fill(ys.begin(), ys.end(), NA_INTEGER);
    int* base = ys.begin();
    for (int sp = 0; sp < n_species; ++sp) {
      std::vector<int> z(n_sites);
      for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, psi(i, sp));
      for (int i = 0; i < n_sites; ++i)
        for (int j = 0; j < max_v; ++j) {
          std::size_t off = (std::size_t) sp * sp_stride + (std::size_t) j * n_sites + i;
          if (pvv[off] != 0) base[off] = (int) R::rbinom(1.0, z[i] * p(i, sp));
        }
    }
    ys.attr("dim") = Rcpp::IntegerVector::create(n_sites, max_v, n_species);
    out[s] = ys;
  }
  return out;
}

// Community occupancy + cover (ms_occu_cover; the spatial variant folds its field
// into psi / ep_mat in R). psi [n_sites x n_species]; p_mat / ep_mat / valid are
// flat [n_sites x max_v x n_species]. Per species: z ~ Bernoulli(psi); at
// occupied sites the visit detections, and the cover value at detected visits
// (beta / lognormal). Returns a list of nsim list(y, y_pos). RNG order matches.
// [[Rcpp::export]]
Rcpp::List cpp_simulate_ms_occu_cover(
    Rcpp::NumericMatrix psi, Rcpp::NumericVector p_mat, Rcpp::NumericVector ep_mat,
    Rcpp::IntegerVector valid, Rcpp::NumericVector disp, int positive,
    int n_sites, int max_v, int n_species, int nsim
) {
  // positive: 0 = lognormal, 3 = beta, 4 = gaussian (see .occu_cover_pos_code /
  // src/occu_coupling_shared.h).
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  const int* pvv = valid.begin();
  const double* pp = p_mat.begin(); const double* pe = ep_mat.begin();
  const std::size_t sp_stride = (std::size_t) n_sites * max_v;
  const bool disp_vec = disp.size() == n_species;
  for (int s = 0; s < nsim; ++s) {
    Rcpp::IntegerVector y((std::size_t) n_sites * max_v * n_species);
    Rcpp::NumericVector yp((std::size_t) n_sites * max_v * n_species);
    std::fill(y.begin(), y.end(), NA_INTEGER);
    std::fill(yp.begin(), yp.end(), NA_REAL);
    int* by = y.begin(); double* byp = yp.begin();
    for (int sp = 0; sp < n_species; ++sp) {
      double d = disp_vec ? disp[sp] : disp[0];
      std::vector<int> z(n_sites);
      for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, psi(i, sp));
      for (int i = 0; i < n_sites; ++i) {
        // observed visits (valid) of this (site, species)
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
          byp[off] = (positive == 3)
                       ? (R::rbeta(plg(clampe(eta)) * d, (1.0 - plg(clampe(eta))) * d))
                       : (positive == 4)
                         ? R::rnorm(eta, d)
                         : std::exp(R::rnorm(eta, d));
        }
      }
    }
    Rcpp::IntegerVector dim = Rcpp::IntegerVector::create(n_sites, max_v, n_species);
    y.attr("dim") = dim; yp.attr("dim") = dim;
    out[s] = Rcpp::List::create(Rcpp::Named("y") = y, Rcpp::Named("y_pos") = yp);
  }
  return out;
}

// Community multi-season occupancy (ms_dyn_occu). psi1 / p [n_sites x n_species];
// gamma / eps [n_sites]; valid flat [n_sites x max_v x n_seasons x n_species].
// Per (species, site): the season-1 state, the survival/colonisation transitions
// (eps then gamma each interval), then the per-season detections. RNG matches.
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_simulate_ms_dyn_occu(
    Rcpp::NumericMatrix psi1, Rcpp::NumericMatrix p, Rcpp::NumericVector gamma,
    Rcpp::NumericVector eps, Rcpp::IntegerVector valid,
    int n_sites, int max_v, int n_seasons, int n_species, int nsim
) {
  Rcpp::RNGScope scope;
  const std::size_t sp_stride = (std::size_t) n_sites * max_v * n_seasons;
  Rcpp::IntegerVector out((std::size_t) sp_stride * n_species * nsim);
  std::fill(out.begin(), out.end(), NA_INTEGER);
  const int* pvv = valid.begin();
  const std::size_t sim_stride = (std::size_t) sp_stride * n_species;
  for (int s = 0; s < nsim; ++s) {
    int* base = out.begin() + (std::size_t) s * sim_stride;
    for (int sp = 0; sp < n_species; ++sp) {
      for (int i = 0; i < n_sites; ++i) {
        std::vector<int> z(n_seasons);
        z[0] = (int) R::rbinom(1.0, psi1(i, sp));
        for (int t = 1; t < n_seasons; ++t) {
          int surv = (int) R::rbinom(1.0, eps[i]);
          int col  = (int) R::rbinom(1.0, gamma[i]);
          z[t] = z[t - 1] * (1 - surv) + (1 - z[t - 1]) * col;
        }
        for (int t = 0; t < n_seasons; ++t)
          for (int j = 0; j < max_v; ++j) {
            std::size_t off = (std::size_t) sp * sp_stride +
              (std::size_t) t * n_sites * max_v + (std::size_t) j * n_sites + i;
            if (pvv[off] != 0) base[off] = (int) R::rbinom(1.0, z[t] * p(i, sp));
          }
      }
    }
  }
  out.attr("dim") = Rcpp::IntegerVector::create(n_sites, max_v, n_seasons, n_species, nsim);
  return out;
}

// Community integrated multi-source occupancy (ms_int_occu). psi [n_sites x
// n_species]; pd a list of D matrices [n_sites x n_species]; valid a list of D
// flat arrays [n_sites x J_d x n_species]. Per species: z ~ Bernoulli(psi); per
// source, the observed visits' detections at occupied sites. Returns a list of
// nsim, each a list of D arrays. RNG order matches the R loop.
// [[Rcpp::export]]
Rcpp::List cpp_simulate_ms_int_occu(
    Rcpp::NumericMatrix psi, Rcpp::List pd, Rcpp::List valid,
    Rcpp::IntegerVector J_d, int n_sites, int n_species, int D, int nsim
) {
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  for (int s = 0; s < nsim; ++s) {
    Rcpp::List srcs(D);
    std::vector<Rcpp::IntegerVector> arr(D);
    for (int d = 0; d < D; ++d) {
      arr[d] = Rcpp::IntegerVector((std::size_t) n_sites * J_d[d] * n_species);
      std::fill(arr[d].begin(), arr[d].end(), NA_INTEGER);
    }
    for (int sp = 0; sp < n_species; ++sp) {
      std::vector<int> z(n_sites);
      for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, psi(i, sp));
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
    out[s] = srcs;
  }
  return out;
}
