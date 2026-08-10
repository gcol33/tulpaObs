// occu_cover_ragged.h
// Predictor assembly for the compact (ragged) occu_cover() kernels.
//
// A compact occu_cover model stores one row per VALID visit: `site_of_visit`
// maps each visit to its site, the visit-level designs carry V rows, and the
// site-level designs plus the shared-field contributions carry n_sites rows.
// Three kernels read exactly that layout -- the pointwise log-likelihood
// (occu_cover_ploglik.cpp), the PIT / LOO-PIT CDF limits and the posterior
// predictive check (occu_cover_diag.cpp) -- and each needs the same three
// predictors: the occupancy logit with the shared field folded in, the
// detection logit (site block + visit block), and the cover predictor (site
// block + copied field + visit block). They are assembled here once, so the
// three kernels cannot drift in how those blocks combine.
//
// The detection and cover arms also take an optional per-visit offset [V x S],
// one value per (visit, draw): what an observation-arm random effect adds to
// that arm's linear predictor at that draw. It enters the same place the
// visit-level design block does, so a kernel scores the predictor the fit was
// evaluated at without knowing where the offset came from.
//
// A dense (padded [n_sites x max_visits]) model reaches the same kernels
// through .occu_cover_visit_view() in R, which flattens its grid to this
// one-row-per-valid-visit form in site-major, visit-ascending order. Dense and
// compact builds of the same data therefore hand these kernels byte-identical
// input.

#ifndef TULPAOBS_OCCU_COVER_RAGGED_H
#define TULPAOBS_OCCU_COVER_RAGGED_H

#include <Rcpp.h>
#include <cstddef>
#include "tobs_math.h"

namespace tulpaObs {
namespace occu_cover_ragged {

// Column-major raw views of the arm designs, the per-draw coefficient matrices
// [S x p], and the per-site field contributions [n_sites x S]. The cover-arm
// members stay null for a kernel that scores only the detection summary.
struct Arms {
  const double* X_occ = nullptr;
  const double* X_det_site = nullptr;
  const double* X_det_visit = nullptr;
  const double* X_pos_site = nullptr;
  const double* X_pos_visit = nullptr;
  const double* b_occ = nullptr;
  const double* b_det = nullptr;
  const double* b_pos = nullptr;
  const double* field_occ = nullptr;
  const double* field_pos = nullptr;
  const double* off_det_visit = nullptr;      // [V x S] or null
  const double* off_pos_visit = nullptr;      // [V x S] or null
  const int*    site_of_visit = nullptr;      // [V], 1-based
  int n_sites = 0, V = 0, S = 0;
  int p_occ = 0, p_det_site = 0, p_det_vis = 0;
  int p_pos_site = 0, p_pos_vis = 0;
  double eta_bound = kEtaClampBound;

  bool has_det_visit() const { return p_det_vis > 0; }
  bool has_pos_visit() const { return p_pos_vis > 0; }

  // 0-based site of visit v.
  int site(int v) const { return site_of_visit[v] - 1; }

  // Occupancy logit at site i under draw d, shared field folded in.
  double eta_psi(int i, int d) const {
    return row_draw_dot(X_occ, n_sites, i, b_occ, S, d, 0, p_occ) +
           field_occ[(std::size_t) d * n_sites + i];
  }
  double psi(int i, int d) const {
    return stable_plogis(clamp_eta(eta_psi(i, d), eta_bound));
  }

  // Site-level detection / cover blocks; the visit-level block is added per
  // visit so a site's block is computed once and reused across its visits.
  double eta_p_site(int i, int d) const {
    return row_draw_dot(X_det_site, n_sites, i, b_det, S, d, 0, p_det_site);
  }
  double eta_pos_site(int i, int d) const {
    return row_draw_dot(X_pos_site, n_sites, i, b_pos, S, d, 0, p_pos_site) +
           field_pos[(std::size_t) d * n_sites + i];
  }

  // Visit-level blocks: the visit-level design contribution (zero when that arm
  // carries no visit-level design) plus the per-visit offset an observation-arm
  // random effect contributes at this draw (zero when the arm carries none).
  double eta_p_visit(int v, int d) const {
    double e = has_det_visit()
      ? row_draw_dot(X_det_visit, V, v, b_det, S, d, p_det_site, p_det_vis)
      : 0.0;
    if (off_det_visit) e += off_det_visit[(std::size_t) d * V + v];
    return e;
  }
  double eta_pos_visit(int v, int d) const {
    double e = has_pos_visit()
      ? row_draw_dot(X_pos_visit, V, v, b_pos, S, d, p_pos_site, p_pos_vis)
      : 0.0;
    if (off_pos_visit) e += off_pos_visit[(std::size_t) d * V + v];
    return e;
  }
};

// One [V x S] per-visit offset matrix, or null when the caller passes a
// zero-column matrix (the same "arm carries none" convention as the visit-level
// designs). A present matrix must carry one row per valid visit and one column
// per draw.
inline const double* visit_offset(const Rcpp::NumericMatrix& off, int V, int S,
                                  const char* what) {
  if (off.ncol() == 0) return nullptr;
  if (off.nrow() != V || off.ncol() != S) {
    Rcpp::stop("%s must be [V x S].", what);
  }
  return off.begin();
}

// Occupancy + detection arms: what the CDF-limits kernel needs.
inline Arms make_arms(const Rcpp::NumericMatrix& X_occ,
                      const Rcpp::NumericMatrix& X_det_site,
                      const Rcpp::NumericMatrix& X_det_visit,
                      const Rcpp::IntegerVector& site_of_visit,
                      const Rcpp::NumericMatrix& b_occ,
                      const Rcpp::NumericMatrix& b_det,
                      const Rcpp::NumericMatrix& field_occ,
                      double eta_bound) {
  Arms a;
  a.n_sites = X_occ.nrow();
  a.V       = site_of_visit.size();
  a.S       = b_occ.nrow();
  a.p_occ      = X_occ.ncol();
  a.p_det_site = X_det_site.ncol();
  a.p_det_vis  = X_det_visit.ncol();
  if (X_det_site.nrow() != a.n_sites) {
    Rcpp::stop("X_det_site must have one row per site.");
  }
  if (a.p_det_vis > 0 && X_det_visit.nrow() != a.V) {
    Rcpp::stop("X_det_visit must have one row per valid visit.");
  }
  if (b_det.ncol() != a.p_det_site + a.p_det_vis || b_occ.ncol() != a.p_occ) {
    Rcpp::stop("coefficient draws do not match the arm designs.");
  }
  if (field_occ.nrow() != a.n_sites || field_occ.ncol() != a.S) {
    Rcpp::stop("field_occ must be [n_sites x S].");
  }
  a.X_occ = X_occ.begin();
  a.X_det_site = X_det_site.begin();
  a.X_det_visit = X_det_visit.begin();
  a.b_occ = b_occ.begin();
  a.b_det = b_det.begin();
  a.field_occ = field_occ.begin();
  a.site_of_visit = site_of_visit.begin();
  a.eta_bound = eta_bound;
  return a;
}

// Add the cover arm: what the pointwise log-likelihood and the PPC need on top.
inline void attach_cover(Arms& a, const Rcpp::NumericMatrix& X_pos_site,
                         const Rcpp::NumericMatrix& X_pos_visit,
                         const Rcpp::NumericMatrix& b_pos,
                         const Rcpp::NumericMatrix& field_pos) {
  a.p_pos_site = X_pos_site.ncol();
  a.p_pos_vis  = X_pos_visit.ncol();
  if (X_pos_site.nrow() != a.n_sites) {
    Rcpp::stop("X_pos_site must have one row per site.");
  }
  if (a.p_pos_vis > 0 && X_pos_visit.nrow() != a.V) {
    Rcpp::stop("X_pos_visit must have one row per valid visit.");
  }
  if (b_pos.ncol() != a.p_pos_site + a.p_pos_vis) {
    Rcpp::stop("coefficient draws do not match the cover-arm design.");
  }
  if (field_pos.nrow() != a.n_sites || field_pos.ncol() != a.S) {
    Rcpp::stop("field_pos must be [n_sites x S].");
  }
  a.X_pos_site = X_pos_site.begin();
  a.X_pos_visit = X_pos_visit.begin();
  a.b_pos = b_pos.begin();
  a.field_pos = field_pos.begin();
}

}  // namespace occu_cover_ragged
}  // namespace tulpaObs

#endif  // TULPAOBS_OCCU_COVER_RAGGED_H
