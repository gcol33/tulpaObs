# Shared numerical-stability helpers.

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
# boundary (cover exactly 0 or 1) instead of emitting -Inf (gcol33/tulpaObs#133).
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
