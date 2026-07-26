# =============================================================================
# grid_summary.R - summarising a nested-Laplace outer grid
#
# Every nested-Laplace fit in the package returns, per outer-grid point k, a
# weight w_k (the normalised marginal), the mode m_k of the inner Gaussian, and
# -- where the kernel keeps them -- the within-grid covariance C_k at that mode.
# Two summaries of that grid recur across every areal / joint / community
# fitter, and both live here.
# =============================================================================


# Posterior mean and sd of a scalar hyperparameter carried as one column of the
# outer grid: the first two weighted moments of `v` under the grid weights `w`.
.tobs_weighted_moment <- function(w, v) {
  m <- sum(w * v, na.rm = TRUE)
  s <- sqrt(max(0, sum(w * v^2, na.rm = TRUE) - m^2))
  c(mean = m, sd = s)
}


# Marginal covariance of a coefficient (or coefficient + field) block by the law
# of total covariance over the outer grid:
#
#   V = sum_k w_k [ C_k + (m_k - mbar)(m_k - mbar)' ],   mbar = sum_k w_k m_k
#
# the within-grid Laplace covariance plus the between-grid spread of the modes.
# This is the marginal object -- the one that gives calibrated standard errors
# once the hyperparameters are integrated out, rather than the conditional
# covariance at a single grid point.
#
# `modes` is [n_grid x p], one mode per row; `weights` the normalised grid
# weights; `blocks` the per-grid within covariances, or NULL to keep only the
# between-grid term. `center` supplies mbar when the caller already has it.
# Grid points with a non-positive or non-finite weight are dropped. A NULL / NA
# entry of `blocks` is handled per `on_missing`: "skip" drops that grid point
# entirely, "zero" keeps its between-grid term with a zero within term.
.tobs_grid_vcov <- function(modes, weights, blocks = NULL, center = NULL,
                            on_missing = c("skip", "zero"), symmetrize = FALSE) {
  on_missing <- match.arg(on_missing)
  modes <- as.matrix(modes)
  p     <- ncol(modes)
  mbar  <- if (is.null(center)) as.numeric(crossprod(weights, modes)) else center

  V <- matrix(0, p, p)
  for (k in seq_len(nrow(modes))) {
    wk <- weights[k]
    if (!isTRUE(is.finite(wk) && wk > 0)) next
    Ck      <- if (is.null(blocks)) NULL else blocks[[k]]
    missing <- !is.null(blocks) && (is.null(Ck) || anyNA(Ck))
    if (missing && identical(on_missing, "skip")) next
    within  <- if (is.null(Ck) || missing) matrix(0, p, p) else as.matrix(Ck)
    dk      <- modes[k, ] - mbar
    V <- V + wk * (within + tcrossprod(dk))
  }
  if (symmetrize) V <- (V + t(V)) / 2
  V
}
