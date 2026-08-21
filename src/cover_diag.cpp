// cover_diag.cpp
// C++ kernels for the cover() hurdle posterior diagnostics whose per-draw /
// per-observation loops were pure R:
//   * cpp_cover_pit_cdf -- the deterministic predictive-CDF limits at occupied
//     sites (.tobs_pit_cover), the randomized-PIT building block. All four
//     positive families.
//   * cpp_cover_ppc -- the posterior predictive check (.tobs_ppc_cover). The
//     occurrence and cover replicates are drawn from R's RNG stream via the R::
//     samplers in the SAME order as the former R loop, so under a fixed seed the
//     discrepancy is byte-identical. Serial (the RNG stream is ordered).
// Family codes: 0 lognormal, 1 lognormal_trunc, 2 ordinal, 3 beta.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include "tobs_math.h"
#include "tobs_shape.h"
using namespace Rcpp;
using tulpaObs::stable_plogis;
using tulpaObs::ppc_stat;
namespace shape = tulpaObs::shape;

namespace {
inline double pnorm1(double z) { return 0.5 * std::erfc(-z * M_SQRT1_2); }
}  // namespace

// Deterministic predictive-CDF limits [S x N]. Absent sites: Fl = 0,
// Fu = 1 - p. Occupied sites (ordinal): a class jump [1-p + p F(lo), 1-p + p
// F(up)]; other families: a continuous point 1 - p + p F_pos (Fl = Fu).
// [[Rcpp::export]]
Rcpp::List cpp_cover_pit_cdf(
    Rcpp::NumericMatrix eta_occ,   // [S x N]
    Rcpp::NumericMatrix eta_pos,   // [S x N_pos]
    Rcpp::IntegerVector occur,     // [N] 0/1
    Rcpp::NumericVector y_pos,     // [N_pos] (log-cover for logn/trunc/ordinal)
    Rcpp::IntegerVector pos_col,   // [N] 1-based eta_pos column, 0 if absent
    Rcpp::NumericVector disp,      // [S]
    int positive,
    Rcpp::NumericVector lower, Rcpp::NumericVector upper,
    Rcpp::NumericVector trunc_upper
) {
  const int S = eta_occ.nrow(), N = eta_occ.ncol();
  const int N_pos = eta_pos.ncol();
  shape::check_nrow(eta_pos, S, "eta_pos");
  shape::check_len(occur, N, "occur");
  shape::check_len(pos_col, N, "pos_col");
  shape::check_len(disp, S, "disp");
  shape::check_len(y_pos, N_pos, "y_pos");
  if (positive == 2) {                           // ordinal class bounds
    shape::check_len(lower, N_pos, "lower");
    shape::check_len(upper, N_pos, "upper");
  } else if (positive == 1) {                    // lognormal_trunc ceiling
    shape::check_len(trunc_upper, N_pos, "trunc_upper");
  }
  // Every plot with occurrence 1 carries a positive-arm row, and pos_col holds
  // the 1-based eta_pos column of that row. The 0 marking a plot with no such
  // row would read the column before the first.
  for (int i = 0; i < N; ++i) {
    if (occur[i] != 1) continue;
    const int j = pos_col[i] - 1;
    if (j < 0 || j >= N_pos) {
      Rcpp::stop("pos_col[%d] = %d at a plot with occurrence 1; expected an "
                 "eta_pos column in [1, %d].", i + 1, (int) pos_col[i], N_pos);
    }
  }
  Rcpp::NumericMatrix Fl(S, N), Fu(S, N);
  for (int i = 0; i < N; ++i) {
    for (int s = 0; s < S; ++s) {
      double p = stable_plogis(eta_occ(s, i));
      Fl(s, i) = 0.0; Fu(s, i) = 1.0 - p;
    }
    if (occur[i] != 1) continue;
    int j = pos_col[i] - 1;
    for (int s = 0; s < S; ++s) {
      double p = stable_plogis(eta_occ(s, i)), one_mp = 1.0 - p;
      double e = eta_pos(s, j), sd = disp[s];
      if (positive == 2) {                         // ordinal: class jump
        double Flp = pnorm1((lower[j] - e) / sd), Fup = pnorm1((upper[j] - e) / sd);
        Fl(s, i) = one_mp + p * Flp; Fu(s, i) = one_mp + p * Fup;
      } else {
        double Fpos;
        // codes 0 (lognormal, y on log-scale) and 4 (identity-Gaussian, y raw)
        // both reduce to the Gaussian CDF at the stored y (#112).
        if (positive == 0 || positive == 4) Fpos = pnorm1((y_pos[j] - e) / sd);
        else if (positive == 1)
          Fpos = pnorm1((y_pos[j] - e) / sd) / pnorm1((trunc_upper[j] - e) / sd);
        else { double mu = stable_plogis(e); Fpos = R::pbeta(y_pos[j], mu * sd, (1.0 - mu) * sd, 1, 0); }
        double val = one_mp + p * Fpos;
        Fl(s, i) = val; Fu(s, i) = val;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("cdf_lower") = Fl,
                            Rcpp::Named("cdf_upper") = Fu);
}

// Posterior predictive check. Per draw s: occurrence replicate ~ Bernoulli(p_s)
// (N draws), then the positive replicate over the n_pos occupied plots in the
// family's order -- matching the R loop's RNG order exactly.
// [[Rcpp::export]]
Rcpp::List cpp_cover_ppc(
    Rcpp::NumericMatrix eta_occ,   // [S x N]
    Rcpp::NumericMatrix eta_pos,   // [S x N_pos]
    Rcpp::IntegerVector occur,     // [N]
    Rcpp::NumericVector y_pos_nat, // [N_pos] observed cover, natural scale
    Rcpp::NumericVector disp,      // [S]
    Rcpp::NumericVector trunc_upper, // [N_pos] (family 1 only)
    int positive, bool freeman
) {
  const int S = eta_occ.nrow(), N = eta_occ.ncol(), n_pos = y_pos_nat.size();
  shape::check_dim(eta_pos, S, n_pos, "eta_pos");
  shape::check_len(occur, N, "occur");
  shape::check_len(disp, S, "disp");
  if (positive == 1) {                           // lognormal_trunc ceiling
    shape::check_len(trunc_upper, n_pos, "trunc_upper");
  }
  Rcpp::RNGScope scope;
  Rcpp::NumericVector fit_y(S), fit_rep(S);
  auto stat = [&](double o, double e) { return ppc_stat(o, e, freeman); };
  for (int s = 0; s < S; ++s) {
    double occ_obs = 0.0, occ_rp = 0.0;
    std::vector<int> occ_rep(N);
    for (int i = 0; i < N; ++i) {
      double p = stable_plogis(eta_occ(s, i));
      occ_rep[i] = (int) R::rbinom(1.0, p);
    }
    for (int i = 0; i < N; ++i) {
      double p = stable_plogis(eta_occ(s, i));
      occ_obs += stat((double) occur[i], p);
      occ_rp  += stat((double) occ_rep[i], p);
    }
    double pos_obs = 0.0, pos_rp = 0.0;
    if (n_pos > 0) {
      double sd = disp[s];
      std::vector<double> Epos(n_pos), yrep(n_pos);
      for (int j = 0; j < n_pos; ++j) {
        double e = eta_pos(s, j);
        if (positive == 1) {                        // lognormal_trunc
          double za = (trunc_upper[j] - e) / sd, Phi_za = pnorm1(za);
          Epos[j] = std::exp(e + sd * sd / 2.0) * pnorm1(za - sd) / Phi_za;
          yrep[j] = std::exp(e + sd * R::qnorm(R::unif_rand() * Phi_za, 0, 1, 1, 0));
        } else if (positive == 0 || positive == 2) { // lognormal / ordinal
          Epos[j] = std::exp(e + sd * sd / 2.0);
          yrep[j] = std::exp(R::rnorm(e, sd));
        } else if (positive == 4) {                  // identity-Gaussian (#112)
          Epos[j] = e;                               // mu = eta, response scale
          yrep[j] = R::rnorm(e, sd);
        } else {                                     // beta
          double mu = stable_plogis(e); Epos[j] = mu;
          yrep[j] = R::rbeta(mu * sd, (1.0 - mu) * sd);
        }
      }
      for (int j = 0; j < n_pos; ++j) {
        pos_obs += stat(y_pos_nat[j], Epos[j]);
        pos_rp  += stat(yrep[j], Epos[j]);
      }
    }
    fit_y[s] = occ_obs + pos_obs;
    fit_rep[s] = occ_rp + pos_rp;
  }
  return Rcpp::List::create(Rcpp::Named("fit.y") = fit_y,
                            Rcpp::Named("fit.y.rep") = fit_rep);
}
