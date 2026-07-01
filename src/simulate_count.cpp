// simulate_count.cpp
// C++ generators for the count / multistate family simulate() methods (N-mixture,
// removal, false-positive occupancy). Each selects a posterior draw with
// R_unif_index (the sample.int primitive) and draws the latent abundance / state
// and the observations from R's RNG stream via the R:: samplers, in the SAME
// order as the former R loops, so under a fixed seed each is byte-identical.

#include <Rcpp.h>
#include <R_ext/Random.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

namespace {
inline double plg(double x) {
  if (x >= 0.0) { double z = std::exp(-x); return 1.0 / (1.0 + z); }
  double z = std::exp(x); return z / (1.0 + z);
}
inline double rowdot(const double* X, int nrow, int i, const double* dr,
                     int ndr, int idx, int off, int p) {
  double a = 0.0;
  for (int k = 0; k < p; ++k)
    a += X[(std::size_t) k * nrow + i] * dr[(std::size_t) (off + k) * ndr + idx];
  return a;
}
// Latent abundance: NB (mu parameterisation) or Poisson. R's rnbinom(size, mu)
// is rpois(rgamma(size, mu / size)) internally, so replicating that with the
// exposed R:: samplers is byte-identical (same two-draw stream).
inline int draw_N(double lambda, bool is_nb, double r_size) {
  if (is_nb && R_finite(r_size)) {
    if (lambda <= 0.0) return 0;
    return (int) R::rpois(R::rgamma(r_size, lambda / r_size));
  }
  return (int) R::rpois(lambda);
}
}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_simulate_nmix(
    Rcpp::NumericMatrix X_lambda, Rcpp::NumericMatrix X_p, Rcpp::NumericMatrix draws,
    Rcpp::IntegerVector site_idx, Rcpp::IntegerVector visit_idx,
    int n_sites, int max_visits, int p_lam, int p_p,
    bool is_nb, double r_size, int nsim
) {
  const int ndr = draws.nrow(), n_obs = site_idx.size();
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  const double* pXl = X_lambda.begin(); const double* pXp = X_p.begin();
  const double* pd = draws.begin();
  for (int s = 0; s < nsim; ++s) {
    int idx = (int) R_unif_index((double) ndr);
    std::vector<int> N(n_sites);
    for (int i = 0; i < n_sites; ++i)
      N[i] = draw_N(std::exp(rowdot(pXl, n_sites, i, pd, ndr, idx, 0, p_lam)), is_nb, r_size);
    Rcpp::IntegerMatrix ys(n_sites, max_visits);
    std::fill(ys.begin(), ys.end(), NA_INTEGER);
    for (int k = 0; k < n_obs; ++k) {
      double po = plg(rowdot(pXp, n_obs, k, pd, ndr, idx, p_lam, p_p));
      ys(site_idx[k] - 1, visit_idx[k] - 1) = (int) R::rbinom((double) N[site_idx[k] - 1], po);
    }
    out[s] = ys;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_simulate_removal(
    Rcpp::NumericMatrix X_lambda, Rcpp::NumericMatrix X_p, Rcpp::NumericMatrix draws,
    Rcpp::IntegerVector site_idx, Rcpp::IntegerVector visit_idx,
    int n_sites, int n_pass, int p_lam, int p_p,
    bool is_nb, double r_size, int nsim
) {
  const int ndr = draws.nrow(), n_obs = site_idx.size();
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  const double* pXl = X_lambda.begin(); const double* pXp = X_p.begin();
  const double* pd = draws.begin();
  for (int s = 0; s < nsim; ++s) {
    int idx = (int) R_unif_index((double) ndr);
    std::vector<int> N(n_sites);
    for (int i = 0; i < n_sites; ++i)
      N[i] = draw_N(std::exp(rowdot(pXl, n_sites, i, pd, ndr, idx, 0, p_lam)), is_nb, r_size);
    // Per-site-pass detection prob (long form scattered to [n_sites x n_pass]).
    std::vector<double> pmat((std::size_t) n_sites * n_pass, NA_REAL);
    for (int k = 0; k < n_obs; ++k)
      pmat[(std::size_t) (visit_idx[k] - 1) * n_sites + (site_idx[k] - 1)] =
        plg(rowdot(pXp, n_obs, k, pd, ndr, idx, p_lam, p_p));
    Rcpp::IntegerMatrix ys(n_sites, n_pass);
    for (int i = 0; i < n_sites; ++i) {
      int rem = N[i];
      for (int k = 0; k < n_pass; ++k) {
        int yk = (int) R::rbinom((double) rem, pmat[(std::size_t) k * n_sites + i]);
        ys(i, k) = yk; rem -= yk;
      }
    }
    out[s] = ys;
  }
  return out;
}

// Open N-mixture (dyn_abun): per site, initial abundance then a survival +
// recruitment HMM across seasons, observed by binomial detection. Arms are
// per-site (constant across intervals). RNG order is site-major (draw N, then
// per season the transition draws + the J detections), matching the R loop.
// [[Rcpp::export]]
Rcpp::IntegerVector cpp_simulate_dyn_abun(
    Rcpp::NumericMatrix X_lambda, Rcpp::NumericMatrix X_p,
    Rcpp::NumericMatrix X_omega, Rcpp::NumericMatrix X_gamma,
    Rcpp::NumericMatrix draws, int n_sites, int T, int J,
    int p_lam, int p_p, int p_om, int p_gm,
    bool is_nb, double r_disp, int nsim
) {
  const int ndr = draws.nrow();
  const int o_p = p_lam, o_om = p_lam + p_p, o_gm = p_lam + p_p + p_om;
  Rcpp::RNGScope scope;
  // Output: nsim arrays [n_sites x J x T]; return one flat vector when nsim == 1
  // (dim set by R), else the caller reshapes. Here we always return nsim == 1's
  // array (R wraps multi-sim in a list via repeated calls upstream).
  Rcpp::IntegerVector out((std::size_t) n_sites * J * T * (nsim < 1 ? 1 : nsim));
  const double* pXl = X_lambda.begin(); const double* pXp = X_p.begin();
  const double* pXo = X_omega.begin(); const double* pXg = X_gamma.begin();
  const double* pd = draws.begin();
  const std::size_t sim_stride = (std::size_t) n_sites * J * T;
  for (int s = 0; s < nsim; ++s) {
    int idx = (int) R_unif_index((double) ndr);
    int* base = out.begin() + (std::size_t) s * sim_stride;
    for (int i = 0; i < n_sites; ++i) {
      double lambda = std::exp(rowdot(pXl, n_sites, i, pd, ndr, idx, 0, p_lam));
      double pdet = plg(rowdot(pXp, n_sites, i, pd, ndr, idx, o_p, p_p));
      double omega = plg(rowdot(pXo, n_sites, i, pd, ndr, idx, o_om, p_om));
      double gamma = std::exp(rowdot(pXg, n_sites, i, pd, ndr, idx, o_gm, p_gm));
      int N = draw_N(lambda, is_nb, r_disp);
      for (int t = 0; t < T; ++t) {
        if (t > 0) N = (int) R::rbinom((double) N, omega) + (int) R::rpois(gamma);
        for (int j = 0; j < J; ++j)
          base[(std::size_t) i + (std::size_t) j * n_sites + (std::size_t) t * n_sites * J] =
            (int) R::rbinom((double) N, pdet);
      }
    }
  }
  if (nsim >= 1) out.attr("dim") = Rcpp::IntegerVector::create(n_sites, J, T,
                                                              nsim == 1 ? 1 : nsim);
  return out;
}

// [[Rcpp::export]]
Rcpp::List cpp_simulate_fp_occu(
    Rcpp::NumericMatrix X_psi, Rcpp::NumericMatrix X_p11, Rcpp::NumericMatrix X_p10,
    Rcpp::NumericMatrix X_b, Rcpp::NumericMatrix draws,
    int n_sites, int J, int p_psi, int p_p11, int p_p10, int p_b, int nsim
) {
  const int ndr = draws.nrow();
  const int o_p11 = p_psi, o_p10 = p_psi + p_p11, o_b = p_psi + p_p11 + p_p10;
  Rcpp::RNGScope scope;
  Rcpp::List out(nsim);
  const double* pd = draws.begin();
  for (int s = 0; s < nsim; ++s) {
    int idx = (int) R_unif_index((double) ndr);
    std::vector<double> psi(n_sites), p11(n_sites), p10(n_sites), b(n_sites);
    for (int i = 0; i < n_sites; ++i) {
      psi[i] = plg(rowdot(X_psi.begin(), n_sites, i, pd, ndr, idx, 0, p_psi));
      p11[i] = plg(rowdot(X_p11.begin(), n_sites, i, pd, ndr, idx, o_p11, p_p11));
      p10[i] = plg(rowdot(X_p10.begin(), n_sites, i, pd, ndr, idx, o_p10, p_p10));
      b[i]   = plg(rowdot(X_b.begin(),   n_sites, i, pd, ndr, idx, o_b,   p_b));
    }
    std::vector<int> z(n_sites);
    for (int i = 0; i < n_sites; ++i) z[i] = (int) R::rbinom(1.0, psi[i]);
    Rcpp::IntegerMatrix ys(n_sites, J);
    for (int i = 0; i < n_sites; ++i) {
      if (z[i] == 1) {
        for (int j = 0; j < J; ++j) {
          int det = (int) R::rbinom(1.0, p11[i]);   // R draws det (J), then cert (J)
          ys(i, j) = det;                            // placeholder; cert applied below
        }
        for (int j = 0; j < J; ++j) {
          int cert = (int) R::rbinom(1.0, b[i]);
          ys(i, j) = (ys(i, j) == 1) ? (cert == 1 ? 2 : 1) : 0;
        }
      } else {
        for (int j = 0; j < J; ++j) ys(i, j) = (int) R::rbinom(1.0, p10[i]);
      }
    }
    out[s] = ys;
  }
  return out;
}
