// simulate_distance.cpp
// Family generator for distance sampling: ONE replicate at the parameter vector
// `theta` plus each arm's realised offset, called by simulate_distance(),
// simulate(fit) and, per species, the ms_distance() routes. The per-bin
// detection probabilities pi_b reuse the SAME Gauss-Legendre quadrature the
// distance likelihood integrates against (distance_quad.h: dist_build_quad +
// dist_key_deriv), rather than a separate stats::integrate path -- so the
// simulator draws from exactly the pi the model was fit against (one source of
// truth). N is rpois / rnbinom(mu); the bin counts are a multinomial drawn by
// R's own sequential-binomial algorithm.
//
// pi comes from the engine's fixed high-order rule (order 64, the engine
// default) rather than adaptive QUADPACK: it is the model-consistent pi, and
// the detection functions are smooth, so the rule is effectively exact.

#include <Rcpp.h>
#include "distance_quad.h"
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "simulate_helpers.h"
using namespace Rcpp;
using tulpaObs::arm_eta;
using tulpaObs::draw_latent_N;

namespace {
// R's rmultinom(1, N, prob) via sequential conditional binomials (matches the
// stats C routine's draw sequence).
inline void rmultinom_R(int N, const std::vector<double>& prob, std::vector<int>& out) {
  const int K = (int) prob.size();
  double p_tot = 0.0; for (double p : prob) p_tot += p;
  out.assign(K, 0);
  double rem_p = p_tot; int rem_N = N;
  for (int k = 0; k < K - 1 && rem_N > 0; ++k) {
    if (prob[k] > 0.0) {
      double pp = prob[k] / rem_p; if (pp > 1.0) pp = 1.0;
      out[k] = (pp < 1.0) ? (int) R::rbinom((double) rem_N, pp) : rem_N;
      rem_N -= out[k];
    }
    rem_p -= prob[k];
  }
  out[K - 1] = rem_N;
}
}  // namespace

// [[Rcpp::export]]
Rcpp::IntegerMatrix cpp_simulate_distance(
    Rcpp::NumericMatrix X_lambda, Rcpp::NumericMatrix X_sigma,
    Rcpp::NumericVector theta, Rcpp::NumericVector off_lambda,
    Rcpp::NumericVector off_sigma, Rcpp::NumericVector cutpoints,
    int key, int transect, int quad_order, double b_shape,
    int n_sites, int n_bins, int p_lam, int p_sig,
    bool is_nb, double r_size
) {
  std::vector<double> cut(cutpoints.begin(), cutpoints.end());
  // The rule is (cutpoints, transect, quad_order); the order is the caller's
  // `model$quad_order`, not a constant, or a fit at any other order is
  // simulated from per-bin probabilities it was never fit against.
  if (quad_order < 1) Rcpp::stop("quad_order must be a positive integer.");
  tulpaObs::DistQuad quad = tulpaObs::dist_build_quad(cut, transect, quad_order);
  Rcpp::RNGScope scope;
  const double* pXl = X_lambda.begin(); const double* pXs = X_sigma.begin();
  const double* th = theta.begin();
  std::vector<int> N(n_sites);
  for (int i = 0; i < n_sites; ++i)
    N[i] = draw_latent_N(std::exp(arm_eta(pXl, n_sites, i, th, 0, p_lam, off_lambda)),
                         is_nb ? r_size : R_PosInf);
  Rcpp::IntegerMatrix ys(n_sites, n_bins);
  std::vector<double> prob(n_bins + 1);
  std::vector<int> counts;
  for (int i = 0; i < n_sites; ++i) {
    double sigma = std::exp(arm_eta(pXs, n_sites, i, th, p_lam, p_sig, off_sigma));
    double psum = 0.0;
    for (int b = 0; b < n_bins; ++b) {
      double pi_b = 0.0;
      for (int qi = 0; qi < quad.order; ++qi)
        pi_b += quad.base_w[b][qi] *
                tulpaObs::dist_key_deriv(quad.x[b][qi], key, sigma, b_shape).g;
      prob[b] = pi_b; psum += pi_b;
    }
    prob[n_bins] = (1.0 - psum > 0.0) ? (1.0 - psum) : 0.0;
    if (N[i] > 0) {
      rmultinom_R(N[i], prob, counts);
      for (int b = 0; b < n_bins; ++b) ys(i, b) = counts[b];
    }
    // N[i] == 0 leaves the row at 0 (IntegerMatrix default-initialises to 0)
  }
  ys.attr("latent") = Rcpp::wrap(N);
  return ys;
}
