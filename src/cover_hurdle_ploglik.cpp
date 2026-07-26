// cover_hurdle_ploglik.cpp
// Parallel pointwise log-likelihood for the cover() hurdle fit (the standalone
// two-part cover model), the WAIC / PSIS-LOO / stacking input. The R reference
// is .tobs_cover_hurdle_ll (R/family_cover_hurdle.R); this port mirrors it draw
// for draw and is cross-checked byte-close against it (test-cover-ploglik-cpp.R).
//
// Per observation the latent occurrence is a hurdle: absent sites (occur = 0)
// score log(1 - p); present sites (occur = 1) score log p + the positive-arm
// density of the observed cover at that draw's cover predictor. The four
// positive families mirror the R kernel exactly:
//   lognormal        dnorm(logy, eta, sigma) - logy         (log-cover + Jacobian)
//   lognormal_trunc  the above - log Phi((u - eta)/sigma)   (upper-truncated)
//   ordinal          log( Phi(zu) - Phi(zl) )               (interval class mass)
//   beta             dbeta(y, mu phi, (1 - mu) phi)         (mu = plogis(eta))
// Each draw's [N] row is independent, so the draw loop parallelises with no
// shared writes -- the axis WAIC scales on (n.draws x N).

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using tulpaObs::stable_plogis;

namespace {

// log(plogis(x)) = log(1/(1+exp(-x))), the stable two-branch form R's
// plogis(x, log.p = TRUE) uses.
inline double log_plogis(double x) {
  if (x >= 0.0) return -std::log1p(std::exp(-x));
  return x - std::log1p(std::exp(x));
}
inline double log_1m_plogis(double x) { return log_plogis(-x); }

inline double dnorm_log(double x, double mean, double sd) {
  double r = (x - mean) / sd;
  return -0.5 * std::log(2.0 * M_PI) - std::log(sd) - 0.5 * r * r;
}

inline double std_pnorm(double z) { return 0.5 * std::erfc(-z * M_SQRT1_2); }

// log Phi(z); deep left tail via the asymptotic expansion so it never logs 0.
inline double log_pnorm(double z) {
  double p = std_pnorm(z);
  if (p > 0.0) return std::log(p);
  return -0.5 * z * z - std::log(-z) - 0.5 * std::log(2.0 * M_PI);
}

// Family codes: 0 lognormal, 1 lognormal_trunc, 2 ordinal, 3 beta.
}  // namespace

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
  if (disp.size() != S) Rcpp::stop("disp must be length S (expand before call).");
  if (occur.size() != N || pos_col.size() != N)
    Rcpp::stop("occur / pos_col must be length N.");

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
          dens = dnorm_log(y, e, sd) - y;
          break;
        case 1:  // lognormal_trunc
          dens = dnorm_log(y, e, sd) - y - log_pnorm((ptu[j] - e) / sd);
          break;
        case 2: { // ordinal interval class mass
          double zl = (plo[j] - e) / sd;
          double zu = (pup[j] - e) / sd;
          double m  = std_pnorm(zu) - std_pnorm(zl);
          dens = std::log(m > 1e-300 ? m : 1e-300);
          break;
        }
        case 4:  // identity-Gaussian (y is the raw response, no Jacobian)
          dens = dnorm_log(y, e, sd);
          break;
        default: { // beta
          double mu = stable_plogis(e);
          double a  = mu * sd;
          double b  = (1.0 - mu) * sd;
          dens = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) +
                 (a - 1.0) * std::log(y) + (b - 1.0) * std::log(1.0 - y);
          break;
        }
      }
      pll[(std::size_t) i * S + d] = log_plogis(eo) + dens;
    }
  }
  return ll;
}
