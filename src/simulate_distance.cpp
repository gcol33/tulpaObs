// simulate_distance.cpp
// C++ generator for the distance-sampling simulate() method. The per-bin
// detection probabilities pi_b reuse the SAME Gauss-Legendre quadrature the
// distance likelihood integrates against (distance_quad.h: dist_build_quad +
// dist_key_deriv), rather than a separate stats::integrate path -- so the
// simulator draws from exactly the pi the model was fit against (one source of
// truth). The draw selection uses R_unif_index; N is rpois / rnbinom(mu); the
// bin counts are a multinomial drawn by R's own sequential-binomial algorithm.
//
// Because pi now comes from the engine's fixed high-order rule (order 64, the
// engine default) instead of adaptive QUADPACK, simulate() output is not
// byte-identical to the former R path; it is the model-consistent pi and the
// detection functions are smooth, so the rule is effectively exact.

#include <Rcpp.h>
#include <R_ext/Random.h>
#include "distance_quad.h"
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "simulate_helpers.h"
using namespace Rcpp;
using tulpaObs::row_draw_dot;
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
Rcpp::List cpp_simulate_distance(
    Rcpp::NumericMatrix X_lambda, Rcpp::NumericMatrix X_sigma,
    Rcpp::NumericMatrix draws, Rcpp::NumericVector cutpoints,
    int key, int transect, double b_shape,
    int n_sites, int n_bins, int p_lam, int p_sig,
    bool is_nb, double r_size, int nsim
) {
  const int ndr = draws.nrow();
  std::vector<double> cut(cutpoints.begin(), cutpoints.end());
  tulpaObs::DistQuad quad = tulpaObs::dist_build_quad(cut, transect, 64);
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  const double* pXl = X_lambda.begin(); const double* pXs = X_sigma.begin();
  const double* pd = draws.begin();
  for (int s = 0; s < nsim; ++s) {
    int idx = (int) R_unif_index((double) ndr);
    std::vector<int> N(n_sites);
    for (int i = 0; i < n_sites; ++i)
      N[i] = draw_latent_N(std::exp(row_draw_dot(pXl, n_sites, i, pd, ndr, idx, 0, p_lam)), is_nb ? r_size : R_PosInf);
    Rcpp::IntegerMatrix ys(n_sites, n_bins);
    std::vector<double> prob(n_bins + 1);
    std::vector<int> counts;
    for (int i = 0; i < n_sites; ++i) {
      double sigma = std::exp(row_draw_dot(pXs, n_sites, i, pd, ndr, idx, p_lam, p_sig));
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
    out[s] = ys;
  }
  return out;
}
