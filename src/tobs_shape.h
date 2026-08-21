// tobs_shape.h
// Argument-shape validation shared by the diagnostic and pointwise-likelihood
// kernels.
//
// Those kernels read their arguments through .begin() pointers and through
// Rcpp's operator(), neither of which bounds-checks in a release build, so an
// argument whose shape differs from the one the kernel assumes reads heap bytes
// into a fitted statistic rather than erroring. The shapes to assert are the
// same handful everywhere -- a matrix of stated dimensions, a vector of stated
// length, an index vector that has to land inside another argument -- so they
// live here once and every kernel calls them before it takes a pointer.
//
// Each message names the argument and both shapes, so a mismatch reports what
// was expected and what arrived. Call these BEFORE entering an OpenMP region:
// Rcpp::stop() throws, and an exception must not cross a parallel boundary.

#ifndef TULPAOBS_TOBS_SHAPE_H
#define TULPAOBS_TOBS_SHAPE_H

#include <Rcpp.h>

namespace tulpaObs {
namespace shape {

// A scalar dimension the caller passes alongside the data. Rejected before it
// enters the size arithmetic the length checks below are built from.
inline void check_dim_arg(int n, const char* what) {
  if (n < 0) {
    Rcpp::stop("%s must be non-negative; got %d.", what, n);
  }
}

// Matrix of exactly [nrow x ncol].
template <typename M>
inline void check_dim(const M& m, R_xlen_t nrow, R_xlen_t ncol,
                      const char* what) {
  if (m.nrow() != nrow || m.ncol() != ncol) {
    Rcpp::stop("%s must be [%d x %d]; got [%d x %d].", what, nrow, ncol,
               (R_xlen_t) m.nrow(), (R_xlen_t) m.ncol());
  }
}

// Row count alone, for a design whose column count is free.
template <typename M>
inline void check_nrow(const M& m, R_xlen_t nrow, const char* what) {
  if (m.nrow() != nrow) {
    Rcpp::stop("%s must have %d rows; got %d.", what, nrow,
               (R_xlen_t) m.nrow());
  }
}

// At least `ncol` columns, for a draw matrix a kernel reads a prefix of.
template <typename M>
inline void check_ncol_min(const M& m, R_xlen_t ncol, const char* what) {
  if (m.ncol() < ncol) {
    Rcpp::stop("%s must have at least %d columns; got %d.", what, ncol,
               (R_xlen_t) m.ncol());
  }
}

// Vector of exactly length `n`.
template <typename V>
inline void check_len(const V& v, R_xlen_t n, const char* what) {
  if (v.size() != n) {
    Rcpp::stop("%s must have length %d; got %d.", what, n, (R_xlen_t) v.size());
  }
}

// Every entry of a 0-based index vector inside [0, n). The offending entry is
// reported at its 1-based position, the one the R caller sees.
inline void check_index0(const Rcpp::IntegerVector& v, R_xlen_t n,
                         const char* what) {
  for (R_xlen_t k = 0; k < v.size(); ++k) {
    if (v[k] < 0 || (R_xlen_t) v[k] >= n) {
      Rcpp::stop("%s[%d] = %d is outside [0, %d).", what, k + 1, (int) v[k], n);
    }
  }
}

// Every entry of a 1-based index vector inside [1, n].
inline void check_index1(const Rcpp::IntegerVector& v, R_xlen_t n,
                         const char* what) {
  for (R_xlen_t k = 0; k < v.size(); ++k) {
    if (v[k] < 1 || (R_xlen_t) v[k] > n) {
      Rcpp::stop("%s[%d] = %d is outside [1, %d].", what, k + 1, (int) v[k], n);
    }
  }
}

}  // namespace shape
}  // namespace tulpaObs

#endif  // TULPAOBS_TOBS_SHAPE_H
