# Regression tests for the interleaved-missing-visit fix in the C++ NUTS
# occupancy likelihoods (occ_likelihood.h / dyn_occ_likelihood.h /
# integrated_occ_likelihood.h).
#
# The `y` response is stored non-compacted: missing visits keep their column
# position as a -1 sentinel, and `n_visits[i]` holds the COUNT of valid visits.
# The buggy loops ran `j` over `[0, n_visits)` while indexing a `max_visits`-
# strided array, so a missing visit BEFORE a valid one terminated the loop
# early and silently dropped the trailing valid visits. The fix iterates the
# full `max_visits` stride and skips the -1 sentinels via the existing guard.
#
# Invariant under test: two encodings of the SAME information -- the same real
# visits with the unsurveyed slots placed at the back (A) vs the front (B) --
# must produce identical estimates. The Laplace path (R kernel) always honoured
# this; pre-fix the NUTS path did not (on the single-season reproducer it gave
# |A - B| ~ 0.31 on the detection / occupancy probabilities, with the
# leading-NA fit collapsing to the truncated two-visit dataset).

# Both methods must agree to numerical precision: A and B carry byte-identical
# information, so the masked log-likelihood is identical and a single-chain
# fit at a fixed seed is deterministic. Tolerance well below the ~0.31 bug.
.equal_tol <- 1e-3

# tobs validates control options per method, so Laplace cannot carry the
# NUTS sampler knobs. For the occupancy NUTS path `n.iter` is the TOTAL
# iteration count including warmup (see occu_fit.R), so keep n.iter > n.warmup.
.na_test_ctrl <- function(method) {
  if (identical(method, "nuts"))
    list(n.iter = 800L, n.warmup = 400L, n.chains = 1L, seed = 1L,
         verbose = FALSE)
  else list(verbose = FALSE)
}

test_that("single-season occu: leading vs trailing NA give identical fits", {
  skip_if_fast()
  set.seed(20260612)
  N <- 200L; psi_t <- 0.70; p_t <- 0.45
  z  <- rbinom(N, 1L, psi_t)
  mk <- function() ifelse(z == 1L, rbinom(N, 1L, p_t), 0L)
  v1 <- mk(); v2 <- mk(); v3 <- mk(); v4 <- mk()
  NAc <- rep(NA_integer_, N)

  yA <- cbind(v1, v2, v3, v4, NAc, NAc)   # trailing NA: valid in slots 1..4
  yB <- cbind(NAc, NAc, v1, v2, v3, v4)   # leading  NA: valid in slots 3..6
  dimnames(yA) <- dimnames(yB) <- NULL
  dat <- data.frame(dummy = rep(1, N))

  fit_one <- function(y, method) {
    f <- tobs(formula = ~ 1, data = dat, family = occu(),
              detection = ~ 1, y = y, method = method,
              control = .na_test_ctrl(method))
    c(psi = plogis(unname(f$means[["psi_(Intercept)"]])),
      p   = plogis(unname(f$means[["p_(Intercept)"]])))
  }

  for (m in c("laplace", "nuts")) {
    eA <- fit_one(yA, m)
    eB <- fit_one(yB, m)
    expect_equal(eA[["psi"]], eB[["psi"]], tolerance = .equal_tol,
                 info = sprintf("psi A vs B under method = %s", m))
    expect_equal(eA[["p"]],   eB[["p"]],   tolerance = .equal_tol,
                 info = sprintf("p A vs B under method = %s", m))
  }
})


test_that("dynamic occu: leading vs trailing NA give identical fits", {
  skip_if_fast()
  set.seed(20260612)
  N <- 150L; T_seasons <- 2L
  psi1 <- 0.6; gam <- 0.3; eps <- 0.2; p_t <- 0.5

  z <- matrix(NA_integer_, N, T_seasons)
  z[, 1] <- rbinom(N, 1L, psi1)
  for (t in 2:T_seasons)
    z[, t] <- z[, t - 1] * (1 - rbinom(N, 1L, eps)) +
              (1 - z[, t - 1]) * rbinom(N, 1L, gam)

  # Four real visits per (site, season); embed in a 6-visit stride.
  mk <- function(t) {
    sapply(1:4, function(j) ifelse(z[, t] == 1L, rbinom(N, 1L, p_t), 0L))
  }
  to_array <- function(pad_front) {
    a <- array(NA_integer_, dim = c(N, 6L, T_seasons))
    for (t in seq_len(T_seasons)) {
      vis <- mk(t)
      cols <- if (pad_front) 3:6 else 1:4
      a[, cols, t] <- vis
    }
    a
  }
  set.seed(99); yA <- to_array(pad_front = FALSE)   # trailing NA
  set.seed(99); yB <- to_array(pad_front = TRUE)    # leading  NA
  dat <- data.frame(dummy = rep(1, N))

  fit_one <- function(y, method) {
    f <- tobs(formula = ~ 1, data = dat, family = dyn_occu(),
              detection = ~ 1, y = y, colonization = ~ 1, extinction = ~ 1,
              method = method, control = .na_test_ctrl(method))
    plogis(unname(f$means[["p_(Intercept)"]]))
  }

  for (m in c("laplace", "nuts")) {
    pA <- fit_one(yA, m)
    pB <- fit_one(yB, m)
    expect_equal(pA, pB, tolerance = .equal_tol,
                 info = sprintf("p A vs B under method = %s (dynamic)", m))
  }
})
