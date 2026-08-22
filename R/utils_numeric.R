# Shared numerical-stability helpers and family-neutral numerical primitives.
# Nothing here knows about a model type, an arm layout, or a likelihood.

# Linear-predictor clamp bound. A logit clamped to [-30, 30] gives detection /
# occupancy probabilities within ~9.4e-14 of {0, 1}, far inside double
# precision, while keeping exp()/plogis() away from overflow. The same bound
# guards log-mean linear predictors on the count arms before exp().
.TOBS_ETA_BOUND <- 30

# Clamp a linear predictor (logit or log-mean) to [-bound, bound] before it
# reaches plogis()/exp(). Elementwise; preserves the dim/names of `e` (pmin/pmax
# copy the attributes of their first argument), so it is safe on vectors and
# matrices alike.
.tobs_clamp_eta <- function(e, bound = .TOBS_ETA_BOUND) {
  pmin(pmax(e, -bound), bound)
}

# log(x) guarded at x <= 0, returning -1e300 rather than -Inf / NaN -- the R twin
# of src/occu_coupling_shared.h::log_safe_. Used by the cover positive-arm density
# so the WAIC / LOO pointwise density and the fit kernel agree at the cover
# boundary (cover exactly 0 or 1) instead of emitting -Inf.
.tobs_log_safe <- function(x) {
  r <- rep(-1e300, length(x))
  p <- !is.na(x) & x > 0
  r[p] <- log(x[p])
  r
}

# Elementwise log(exp(a) + exp(b)), max-shifted for stability. The shifted term
# `pmin(a, b) - m` is <= 0, so the other term contributes exp(0) = 1 exactly and
# the sum is log1p(exp(pmin - m)) -- no catastrophic cancellation. `a` and `b`
# recycle/broadcast as for the arithmetic operators.
.tobs_logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log1p(exp(pmin(a, b) - m))
}

# Draw `n` rows from MVN(mu, sigma) via the Cholesky factor of `sigma`.
#
# Two degenerate cases, both reachable from a non-converged observed-information
# Hessian: a `sigma` carrying non-finite entries yields `mu` repeated (a point
# mass, so downstream draws stay finite), and a finite but non-PD `sigma` falls
# back to independent normals on its diagonal.
.rmvn <- function(n, mu, sigma) {
  p <- length(mu)
  if (any(!is.finite(sigma))) {
    return(matrix(rep(mu, each = n), n, p, byrow = FALSE))
  }
  L <- tryCatch(chol(sigma), error = function(e) NULL)
  z <- matrix(stats::rnorm(n * p), n, p)
  if (is.null(L)) {
    sds <- sqrt(pmax(diag(sigma), 1e-8))
    return(sweep(z * rep(sds, each = n), 2L, mu, "+"))
  }
  sweep(z %*% L, 2L, mu, "+")
}

# Symmetrised central-difference Jacobian of a vector-valued function of `x`.
# Applied to an analytic negative-log-likelihood gradient it IS the observed
# information at the mode, so it stands in wherever a family exposes a gradient
# but no analytic Hessian. The marginal Hessian is symmetric, so the FD estimate
# is symmetrised.
.tobs_fd_jacobian <- function(fn, x, h = 1e-5) {
  p <- length(x)
  J <- matrix(0, p, p)
  for (i in seq_len(p)) {
    xp <- x; xm <- x; xp[i] <- xp[i] + h; xm[i] <- xm[i] - h
    J[, i] <- (fn(xp) - fn(xm)) / (2 * h)
  }
  0.5 * (J + t(J))
}
