// cover_hurdle_ploglik.cpp
// Parallel pointwise log-likelihood for the cover() hurdle fit (the standalone
// two-part cover model), the WAIC / PSIS-LOO / stacking input. The R reference
// is .tobs_cover_hurdle_ll (R/cover_hurdle_diag.R); this port mirrors it draw
// for draw and is cross-checked byte-close against it (test-cover-ploglik-cpp.R).
//
// Per observation the latent occurrence is a hurdle: absent sites (occur = 0)
// score log(1 - p); present sites (occur = 1) score log p + the positive-arm
// density of the observed cover at that draw's cover predictor. The lognormal,
// gaussian and beta arms evaluate the engine's own per-observation density
// (tulpa/family_density.h), the function the positive-arm fit maximised:
//   lognormal        dnorm(logy, eta, sigma) - logy         (log-cover + Jacobian)
//   lognormal_trunc  the above - log Phi((u - eta)/sigma)   (upper-truncated)
//   ordinal          log( Phi(zu) - Phi(zl) )               (interval class mass)
//   gaussian         dnorm(y, eta, sigma)                   (raw response)
//   beta             dbeta(y, mu phi, (1 - mu) phi)         (mu = plogis(eta),
//                                                            clamped to the fit's
//                                                            [1e-15, 1 - 1e-15])
// Each draw's [N] row is independent, so the draw loop parallelises with no
// shared writes -- the axis WAIC scales on (n.draws x N).
// Family codes: 0 lognormal, 1 lognormal_trunc, 2 ordinal, 3 beta, 4 gaussian.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <tulpa/family_density.h>
#include "cover_hurdle_shape.h"
#include "tobs_math.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using tulpaObs::log_plogis;
using tulpaObs::log_1m_plogis;

// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_cover_hurdle_ploglik(
    Rcpp::NumericMatrix eta_occ,    // [S x N]
    Rcpp::NumericMatrix eta_pos,    // [S x N_pos]
    Rcpp::NumericVector disp,       // [S] (expanded by the caller)
    Rcpp::IntegerVector occur,      // [N], 0/1
    Rcpp::NumericVector y_pos,      // [N_pos]
    Rcpp::IntegerVector pos_col,    // [N], 1-based eta_pos column, 0 if absent
    int positive,
    Rcpp::NumericVector lower,        // [N_pos] ordinal (empty if unused)
    Rcpp::NumericVector upper,        // [N_pos] ordinal
    Rcpp::NumericVector trunc_upper,  // [N_pos] lognormal_trunc ceiling
    int n_threads
) {
  const int S = eta_occ.nrow();
  const int N = eta_occ.ncol();
  tulpaObs::cover_hurdle::check_arms(eta_occ, eta_pos, occur, y_pos, pos_col,
                                     disp, positive, lower, upper, trunc_upper);

  Rcpp::NumericMatrix ll(S, N);
  const double* peo  = eta_occ.begin();
  const double* pep  = eta_pos.begin();
  const int Spos     = eta_pos.nrow();
  const double* pd   = disp.begin();
  const int*    poc  = occur.begin();
  const double* pyp  = y_pos.begin();
  const int*    ppc  = pos_col.begin();
  const double* plo  = lower.begin();
  const double* pup  = upper.begin();
  const double* ptu  = trunc_upper.begin();
  double* pll        = ll.begin();

#ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(n_threads > 0 ? n_threads : 1)
#endif
  for (int d = 0; d < S; ++d) {
    for (int i = 0; i < N; ++i) {
      double eo = peo[(std::size_t) i * S + d];
      if (poc[i] == 0) {
        pll[(std::size_t) i * S + d] = log_1m_plogis(eo);
        continue;
      }
      int j = ppc[i] - 1;                       // eta_pos column for this site
      double e   = pep[(std::size_t) j * Spos + d];
      double sd  = pd[d];
      double y   = pyp[j];
      double dens;
      switch (positive) {
        case 0:  // lognormal (y is log-cover)
          dens = tulpa::log_lik_lognormal_logy(y, e, sd);
          break;
        case 1:  // lognormal_trunc
          dens = tulpa::log_lik_lognormal_logy(y, e, sd)
                 - tulpa::math::portable_pnorm_log((ptu[j] - e) / sd);
          break;
        case 2: { // ordinal interval class mass
          double zl = (plo[j] - e) / sd;
          double zu = (pup[j] - e) / sd;
          double m  = tulpa::math::portable_pnorm(zu) - tulpa::math::portable_pnorm(zl);
          dens = std::log(m > 1e-300 ? m : 1e-300);
          break;
        }
        case 4:  // identity-Gaussian (y is the raw response, no Jacobian)
          dens = tulpa::log_lik_gaussian(y, e, sd);
          break;
        default:  // beta
          dens = tulpa::log_lik_beta_logit(y, e, sd);
          break;
      }
      pll[(std::size_t) i * S + d] = log_plogis(eo) + dens;
    }
  }
  return ll;
}
