// simulate_single.cpp
// C++ generator for the single-season occupancy simulate() method. The posterior
// draw is selected with R_unif_index (the exact primitive sample.int(n, 1) uses,
// exposed by <R_ext/Random.h>), and the latent state z ~ Bernoulli(psi) and
// detection replicate y ~ Bernoulli(z p) draw from R's RNG stream via R::rbinom,
// in the SAME order as the former R loop, so under a fixed seed the result is
// byte-identical. Returns a list of nsim [n_sites x max_visits] integer matrices
// (< 0 marks the missing visits carried from the observed design).

#include <Rcpp.h>
#include <R_ext/Random.h>   // R_unif_index
#include <vector>
#include <cmath>
using namespace Rcpp;

namespace {
inline double plg(double x) {
  if (x >= 0.0) { double z = std::exp(-x); return 1.0 / (1.0 + z); }
  double z = std::exp(x); return z / (1.0 + z);
}
}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_simulate_single(
    Rcpp::NumericMatrix X_occ,   // [n_sites x p_occ]
    Rcpp::NumericMatrix X_det,   // [n_sites x p_det]
    Rcpp::NumericMatrix draws,   // [ndraws x (p_occ + p_det + ...)]
    Rcpp::IntegerMatrix y,       // [n_sites x max_visits], < 0 = missing
    int p_occ, int p_det, int nsim
) {
  const int n_sites = X_occ.nrow(), max_v = y.ncol();
  const int ndr = draws.nrow();
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  const double* pXo = X_occ.begin(); const double* pXd = X_det.begin();
  const double* pd = draws.begin(); const int* py = y.begin();
  std::vector<double> psi(n_sites), p(n_sites);

  for (int s = 0; s < nsim; ++s) {
    int idx = (int) R_unif_index((double) ndr);   // sample.int(ndraws, 1) - 1
    for (int i = 0; i < n_sites; ++i) {
      double eo = 0.0, ed = 0.0;
      for (int k = 0; k < p_occ; ++k) eo += pXo[(std::size_t) k * n_sites + i] * pd[(std::size_t) k * ndr + idx];
      for (int k = 0; k < p_det; ++k) ed += pXd[(std::size_t) k * n_sites + i] * pd[(std::size_t) (p_occ + k) * ndr + idx];
      psi[i] = plg(eo); p[i] = plg(ed);
    }
    std::vector<int> z(n_sites);
    for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, psi[i]);
    Rcpp::IntegerMatrix ysim(n_sites, max_v);
    std::fill(ysim.begin(), ysim.end(), NA_INTEGER);
    // (i outer, j inner) draw order, matching the R double loop.
    for (int i = 0; i < n_sites; ++i)
      for (int j = 0; j < max_v; ++j)
        if (py[(std::size_t) j * n_sites + i] >= 0)
          ysim(i, j) = (int) R::rbinom(1.0, z[i] * p[i]);
    out[s] = ysim;
  }
  return out;
}
