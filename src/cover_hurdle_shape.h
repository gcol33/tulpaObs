// cover_hurdle_shape.h
// Argument-shape validation for the two cover-hurdle kernels that read the same
// arm layout: the pointwise log-likelihood behind WAIC / LOO / CPO
// (cover_hurdle_ploglik.cpp) and the PIT / LOO-PIT CDF limits (cover_diag.cpp).
//
// Both take [S x N] occurrence draws, [S x N_pos] positive-arm draws, and a
// per-plot pos_col mapping the second onto the first. The mapping is what makes
// a shared validator worth its own file: pos_col holds 0 at a plot with no
// positive-arm row, so a kernel that follows it at an occupied plot without
// checking reads the column before the first. The PIT kernel checked it and the
// WAIC kernel did not, which is the drift this removes.

#ifndef TULPAOBS_COVER_HURDLE_SHAPE_H
#define TULPAOBS_COVER_HURDLE_SHAPE_H

#include "tobs_shape.h"
#include <Rcpp.h>

namespace tulpaObs {
namespace cover_hurdle {

// `positive` selects which of the optional bound vectors is read: 1 =
// lognormal_trunc (trunc_upper), 2 = ordinal (lower / upper). The others carry
// no bounds and their vectors are not inspected.
inline void check_arms(const Rcpp::NumericMatrix& eta_occ,
                       const Rcpp::NumericMatrix& eta_pos,
                       const Rcpp::IntegerVector& occur,
                       const Rcpp::NumericVector& y_pos,
                       const Rcpp::IntegerVector& pos_col,
                       const Rcpp::NumericVector& disp,
                       int positive,
                       const Rcpp::NumericVector& lower,
                       const Rcpp::NumericVector& upper,
                       const Rcpp::NumericVector& trunc_upper) {
  namespace sh = tulpaObs::shape;
  const R_xlen_t S = eta_occ.nrow(), N = eta_occ.ncol();
  const R_xlen_t N_pos = eta_pos.ncol();

  sh::check_nrow(eta_pos, S, "eta_pos");
  sh::check_len(occur, N, "occur");
  sh::check_len(pos_col, N, "pos_col");
  sh::check_len(disp, S, "disp");
  sh::check_len(y_pos, N_pos, "y_pos");
  if (positive == 2) {                           // ordinal class bounds
    sh::check_len(lower, N_pos, "lower");
    sh::check_len(upper, N_pos, "upper");
  } else if (positive == 1) {                    // lognormal_trunc ceiling
    sh::check_len(trunc_upper, N_pos, "trunc_upper");
  }

  // Every plot with occurrence 1 carries a positive-arm row, and pos_col holds
  // the 1-based eta_pos column of that row. The 0 marking a plot with no such
  // row would read the column before the first.
  for (R_xlen_t i = 0; i < N; ++i) {
    if (occur[i] != 1) continue;
    const R_xlen_t j = (R_xlen_t) pos_col[i] - 1;
    if (j < 0 || j >= N_pos) {
      Rcpp::stop("pos_col[%d] = %d at a plot with occurrence 1; expected an "
                 "eta_pos column in [1, %d].", i + 1, (int) pos_col[i], N_pos);
    }
  }
}

}  // namespace cover_hurdle
}  // namespace tulpaObs

#endif  // TULPAOBS_COVER_HURDLE_SHAPE_H
