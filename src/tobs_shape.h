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
#include <vector>

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
// reported at its 1-based position, the one the R caller sees. NA is rejected
// by name: it reaches C++ as INT_MIN, which reads as an arbitrary index rather
// than as the missing value it is.
inline void check_index0(const Rcpp::IntegerVector& v, R_xlen_t n,
                         const char* what) {
  for (R_xlen_t k = 0; k < v.size(); ++k) {
    if (Rcpp::IntegerVector::is_na(v[k])) {
      Rcpp::stop("%s[%d] is NA.", what, k + 1);
    }
    if (v[k] < 0 || (R_xlen_t) v[k] >= n) {
      Rcpp::stop("%s[%d] = %d is outside [0, %d).", what, k + 1, (int) v[k], n);
    }
  }
}

// Every entry of a 1-based index vector inside [1, n].
inline void check_index1(const Rcpp::IntegerVector& v, R_xlen_t n,
                         const char* what) {
  for (R_xlen_t k = 0; k < v.size(); ++k) {
    if (Rcpp::IntegerVector::is_na(v[k])) {
      Rcpp::stop("%s[%d] is NA.", what, k + 1);
    }
    if (v[k] < 1 || (R_xlen_t) v[k] > n) {
      Rcpp::stop("%s[%d] = %d is outside [1, %d].", what, k + 1, (int) v[k], n);
    }
  }
}

// A coefficient block [off, off + width) inside a draw row of `total` columns.
// The offsets and widths are packed in R against the draw matrix and arrive as
// separate scalar arguments, so a kernel reading one is reading two arguments
// against a third. A zero-width block is an absent arm: nothing reads its
// offset, so it passes whatever that offset says.
inline void check_block(int off, int width, R_xlen_t total, const char* what) {
  if (width < 0) {
    Rcpp::stop("%s block width must be non-negative; got %d.", what, width);
  }
  if (width == 0) return;
  if (off < 0 || (R_xlen_t) off + (R_xlen_t) width > total) {
    Rcpp::stop("%s block [%d, %d) does not fit a draw row of %d columns.",
               what, off, off + width, total);
  }
}

// An optional numeric argument, unwrapped at the .cpp export boundary.
// Rcpp::Nullable<T> must not cross into a header helper -- it crashes MinGW --
// so an export keeps the Nullable formal, hands this the SEXP behind it, and
// passes the returned buffer inward. R_NilValue yields an empty vector, which
// is how a header spells "argument absent"; `n`, when non-negative, is the
// length a present argument is required to have.
inline std::vector<double> optional_numeric(SEXP x, const char* what,
                                            R_xlen_t n = -1) {
  if (Rf_isNull(x)) return std::vector<double>();
  Rcpp::NumericVector v(x);
  if (n >= 0 && v.size() != n) {
    Rcpp::stop("%s must have length %d; got %d.", what, n, (R_xlen_t) v.size());
  }
  return std::vector<double>(v.begin(), v.end());
}

}  // namespace shape
}  // namespace tulpaObs

#endif  // TULPAOBS_TOBS_SHAPE_H
