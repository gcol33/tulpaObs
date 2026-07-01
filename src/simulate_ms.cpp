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
