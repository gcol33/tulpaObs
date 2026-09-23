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


# The grid as the mixture `.tobs_grid_mixture_draws()` samples, with the same
# cell handling `.tobs_grid_vcov()` applied to the covariance it reports: under
# `on_missing = "skip"` a cell with no usable within covariance carries no
# weight, under "zero" it keeps its weight and is drawn at its mode.
.tobs_grid_mixture <- function(weights, modes, covs,
                               on_missing = c("skip", "zero")) {
  on_missing <- match.arg(on_missing)
  w <- as.numeric(weights)
  w[!is.finite(w) | w < 0] <- 0
  if (identical(on_missing, "skip"))
    w <- w * !vapply(covs, function(C) is.null(C) || anyNA(C), logical(1))
  list(weights = w, modes = as.matrix(modes), covs = covs)
}


# Draws from the outer-grid posterior itself: the mixture
#
#   x ~ sum_k w_k N(m_k, C_k)
#
# over the leading `ncol(modes)` coordinates, a cell picked by its weight and the
# draw taken from that cell's own inner Gaussian. Its first two moments are the
# mean and the `.tobs_grid_vcov()` covariance of the same grid, so every summary
# that read the collapsed Gaussian keeps its mean and SD; what the mixture adds
# is the shape the weighted cells describe (skew, a long tail along a ridge),
# which one Gaussian at the grid mean cannot carry into a quantile.
#
# `covs = NULL` means the fit kept no per-cell covariance at all (an engine
# without stored per-cell precisions), and the draw is the Gaussian at `means`
# with covariance `V`, the only posterior such a fit describes.
#
# A grid wider than `means` is cut to its first `length(means)` coordinates: a
# fit may report a leading run of the coordinates its grid integrated.
#
# `covs[[k]]` is cell k's within covariance; a NULL / non-finite entry draws that
# cell at its mode, the same "zero within term" `.tobs_grid_vcov(on_missing =
# "zero")` gives it. `V` and `means`, when `V` has more coordinates than
# `modes`, supply the trailing ones: coordinates the grid carries no inner
# Gaussian for (the outer hyperparameters, whose reported SD is the engine's
# per-axis read rather than the grid spread). Those are drawn from their
# Gaussian conditional under `V` given the leading draw, which reproduces `V`'s
# trailing block and its cross-covariance with the leading block exactly, since
# the mixture's own covariance is `V`'s leading block.
# `lead` short-circuits the mixture draw for the fixed-effect block: when the
# engine has already sampled those coordinates exactly -- its subspace debias
# returns a corrected `$draws` -- re-drawing them from the per-cell Gaussians
# throws that sample away and reports the approximation the correction was
# asked to replace (gcol33/tulpa#862). The conditional hyperparameter tail is
# built by the same helper either way.
.tobs_grid_mixture_draws <- function(n, weights, modes, covs, means = NULL,
                                     V = NULL, lead = NULL) {
  if (!is.null(lead)) {
    lead <- as.matrix(lead)
    if (nrow(lead) != n) {
      lead <- lead[sample.int(nrow(lead), n, replace = TRUE), , drop = FALSE]
    }
    return(.tobs_grid_draw_tail(lead, means, V, ncol(lead), n))
  }
  if (is.null(covs)) return(.rmvn(n, means, V))
  modes <- as.matrix(modes)
  if (!is.null(means) && ncol(modes) > length(means)) {
    keep  <- seq_along(means)
    modes <- modes[, keep, drop = FALSE]
    covs  <- lapply(covs, function(C)
      if (is.null(C)) NULL else as.matrix(C)[keep, keep, drop = FALSE])
  }
  p <- ncol(modes)
  w <- as.numeric(weights)
  w[!is.finite(w) | w < 0] <- 0
  if (!any(w > 0)) stop("grid mixture draw: no cell carries positive weight.",
                        call. = FALSE)
  counts <- as.integer(stats::rmultinom(1L, size = n, prob = w / sum(w)))
  lead <- matrix(0, n, p)
  pos <- 0L
  for (k in which(counts > 0L)) {
    rows <- pos + seq_len(counts[k])
    Ck <- covs[[k]]
    lead[rows, ] <- if (is.null(Ck) || !all(is.finite(Ck)))
                      matrix(modes[k, ], counts[k], p, byrow = TRUE)
                    else .rmvn(counts[k], modes[k, ], as.matrix(Ck))
    pos <- pos + counts[k]
  }
  lead <- lead[sample.int(n), , drop = FALSE]
  .tobs_grid_draw_tail(lead, means, V, p, n)
}

# The hyperparameter block, drawn from its Gaussian conditional on the
# fixed-effect draws it is handed. Split out so the mixture path and the
# engine-corrected path cannot build it differently.
.tobs_grid_draw_tail <- function(lead, means, V, p, n) {
  if (is.null(V) || ncol(V) <= p) return(lead)

  b <- seq_len(p); h <- (p + 1L):ncol(V)
  mu_b <- as.numeric(means[b]); mu_h <- as.numeric(means[h])
  if (!all(is.finite(V))) {
    return(cbind(lead, matrix(mu_h, n, length(h), byrow = TRUE)))
  }
  eb <- eigen((V[b, b, drop = FALSE] + t(V[b, b, drop = FALSE])) / 2,
              symmetric = TRUE)
  keep <- eb$values > max(eb$values, 0) * 1e-10
  Vbb_inv <- eb$vectors[, keep, drop = FALSE] %*%
             (t(eb$vectors[, keep, drop = FALSE]) / eb$values[keep])
  K  <- V[h, b, drop = FALSE] %*% Vbb_inv
  Vc <- V[h, h, drop = FALSE] - K %*% V[b, h, drop = FALSE]
  ec <- eigen((Vc + t(Vc)) / 2, symmetric = TRUE)
  Lc <- t(ec$vectors %*% diag(sqrt(pmax(ec$values, 0)), nrow = length(h)))
  dev  <- sweep(lead, 2L, mu_b, "-")
  tail <- sweep(dev %*% t(K) + matrix(stats::rnorm(n * length(h)), n) %*% Lc,
                2L, mu_h, "+")
  cbind(lead, tail)
}
