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

# Elementwise log(exp(a) + exp(b)), max-shifted for stability. The shifted term
# `pmin(a, b) - m` is <= 0, so the other term contributes exp(0) = 1 exactly and
# the sum is log1p(exp(pmin - m)) -- no catastrophic cancellation. `a` and `b`
# recycle/broadcast as for the arithmetic operators.
.tobs_logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log1p(exp(pmin(a, b) - m))
}
