// simulate_single.cpp
// Family generator for single-season occupancy: ONE replicate at the parameter
// vector `theta` (occupancy block, then detection block) plus the realised
// offsets -- per site on occupancy, per (site, visit) on detection -- drawing
// the latent state z ~ Bernoulli(psi) and the detection replicate
// y ~ Bernoulli(z p) from R's RNG stream via R::rbinom in a fixed order.
// simulate_occu(), simulate(fit) and the SBC replicate all call it. Returns an
// [n_sites x max_visits] integer matrix, NA at the visits the observed design
// marks missing (< 0).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "simulate_helpers.h"
using namespace Rcpp;
using tulpaObs::arm_eta;

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_simulate_single(
    Rcpp::NumericMatrix X_occ,   // [n_sites x p_occ]
    Rcpp::NumericMatrix X_det,   // [n_sites x p_det]
    Rcpp::NumericVector theta,   // [p_occ + p_det + ...]
    Rcpp::NumericVector off_occ, // [n_sites], or empty
    Rcpp::NumericMatrix off_det, // [n_sites x max_visits], or 0 x 0
    Rcpp::IntegerMatrix y,       // [n_sites x max_visits], < 0 = missing
    int p_occ, int p_det
) {
  const int n_sites = X_occ.nrow(), max_v = y.ncol();
  const bool has_off_det = off_det.nrow() > 0;
  const Rcpp::NumericVector none(0);
  Rcpp::RNGScope scope;
  const double* pXo = X_occ.begin(); const double* pXd = X_det.begin();
  const double* th = theta.begin(); const int* py = y.begin();
  std::vector<double> psi(n_sites), eta_det(n_sites);
  for (int i = 0; i < n_sites; ++i) {
    psi[i] = tulpa::math::inv_logit(arm_eta(pXo, n_sites, i, th, 0, p_occ, off_occ));
    eta_det[i] = arm_eta(pXd, n_sites, i, th, p_occ, p_det, none);
  }
  std::vector<int> z(n_sites);
  for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, psi[i]);
  Rcpp::IntegerMatrix ysim(n_sites, max_v);
  std::fill(ysim.begin(), ysim.end(), NA_INTEGER);
  // (i outer, j inner) draw order, matching the R double loop.
  for (int i = 0; i < n_sites; ++i)
    for (int j = 0; j < max_v; ++j)
      if (py[(std::size_t) j * n_sites + i] >= 0) {
        double p = tulpa::math::inv_logit(
            has_off_det ? eta_det[i] + off_det(i, j) : eta_det[i]);
        ysim(i, j) = (int) R::rbinom(1.0, z[i] * p);
      }
  ysim.attr("latent") = Rcpp::wrap(z);
  return ysim;
}
