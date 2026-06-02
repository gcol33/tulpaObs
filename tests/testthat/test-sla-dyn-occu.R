# Tier A regression tests for the SLA / finite-difference d3 path on
# the dyn_occu (HMM) family. See R/sla_dyn_occu.R for the math and
# dev_notes/simplified_laplace_derivation.md §2 for the broader framework.


# ---------------------------------------------------------------------------
# Helper: build a small dyn_occu tobs_model directly so we can probe internals
# without paying for a full EM fit.
# ---------------------------------------------------------------------------
make_small_dyn_model <- function(N = 60, J = 3, n_seasons = 3,
                                  beta_psi1 = c(0.3, 0.5),
                                  beta_p    = c(0.2),
                                  beta_gam  = c(-1.5),
                                  beta_eps  = c(-1.0),
                                  seed = 11) {
  set.seed(seed)
  elev <- rnorm(N)
  data <- data.frame(elev = elev)
  X_psi <- model.matrix(~ elev, data)
  X_p   <- model.matrix(~ 1, data)
  X_col <- model.matrix(~ 1, data)
  X_ext <- model.matrix(~ 1, data)

  psi1 <- plogis(as.numeric(X_psi %*% beta_psi1))
  p    <- plogis(as.numeric(X_p   %*% beta_p))
  gam  <- plogis(as.numeric(X_col %*% beta_gam))
  eps  <- plogis(as.numeric(X_ext %*% beta_eps))

  z <- matrix(0L, N, n_seasons)
  z[, 1] <- rbinom(N, 1, psi1)
  for (t in 2:n_seasons) {
    surv <- z[, t - 1] * (1 - rbinom(N, 1, eps))
    col  <- (1 - z[, t - 1]) * rbinom(N, 1, gam)
    z[, t] <- surv + col
  }

  y <- array(NA_integer_, dim = c(N, J, n_seasons))
  for (i in seq_len(N)) {
    for (t in seq_len(n_seasons)) {
      y[i, , t] <- if (z[i, t] == 1L) rbinom(J, 1, p[i]) else 0L
    }
  }

  list(y = y, data = data,
       beta_true = c(beta_psi1, beta_p, beta_gam, beta_eps))
}


# ---------------------------------------------------------------------------
# Smoke test: end-to-end fit with approx = "simplified_laplace"
# ---------------------------------------------------------------------------
test_that("dyn_occu with simplified_laplace attaches a numeric skew vector", {
  skip_if_fast()
  set.seed(201)
  N <- 60; J <- 3; n_seasons <- 3
  data <- data.frame(x = rnorm(N))
  psi1 <- 0.55; p <- 0.5; gam <- 0.2; eps <- 0.1
  z <- matrix(0L, N, n_seasons)
  z[, 1] <- rbinom(N, 1, psi1)
  for (t in 2:n_seasons) {
    z[, t] <- z[, t - 1] * (1 - rbinom(N, 1, eps)) +
              (1 - z[, t - 1]) * rbinom(N, 1, gam)
  }
  y <- array(NA_integer_, dim = c(N, J, n_seasons))
  for (i in seq_len(N)) {
    for (t in seq_len(n_seasons)) {
      y[i, , t] <- if (z[i, t] == 1L) rbinom(J, 1, p) else 0L
    }
  }

  fit <- tobs(
    formula     = ~ 1,
    data        = data,
    family      = dyn_occu(),
    detection   = ~ 1,
    y           = y,
    col_formula = ~ 1,
    ext_formula = ~ 1,
    method      = "laplace_sla",
    control     = list(verbose = FALSE)
  )

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$sla_status, "simplified_laplace")
  expect_true(is.numeric(fit$skew))
  expect_equal(length(fit$skew), fit$n_params)
  expect_named(fit$skew, names(fit$means))
  expect_true(all(is.finite(fit$skew)))
})


# ---------------------------------------------------------------------------
# Log-likelihood correctness: compare against an independent reference HMM
# forward implementation in linear (non-log) space. Both should agree to
# numerical precision on finite-precision-safe inputs (small N, T, J).
#
# This is a stronger correctness check than "is the EM mode a local max"
# because the EM-Laplace mode is the M-step pseudo-binomial-encoded mode,
# which is a *surrogate* for the marginal HMM MLE — they generally don't
# coincide on finite samples, and we deliberately don't want to silently
# rely on coincidence.
# ---------------------------------------------------------------------------
test_that(".loglik_dyn_occu matches an independent reference HMM forward", {
  set.seed(202)
  N <- 30; J <- 3; n_seasons <- 3
  data <- data.frame(x = rnorm(N))
  psi1 <- 0.55; p <- 0.5; gam <- 0.2; eps <- 0.1
  z <- matrix(0L, N, n_seasons)
  z[, 1] <- rbinom(N, 1, psi1)
  for (t in 2:n_seasons) {
    z[, t] <- z[, t - 1] * (1 - rbinom(N, 1, eps)) +
              (1 - z[, t - 1]) * rbinom(N, 1, gam)
  }
  y <- array(NA_integer_, dim = c(N, J, n_seasons))
  for (i in seq_len(N)) {
    for (t in seq_len(n_seasons)) {
      y[i, , t] <- if (z[i, t] == 1L) rbinom(J, 1, p) else 0L
    }
  }
  model <- tulpaObs:::.tobs_build_model(
    occ_formula = ~ 1, det_formula = ~ 1,
    data        = data, y = y,
    col_formula = ~ 1, ext_formula = ~ 1
  )

  # Independent reference forward in linear space (no log-sum-exp trick).
  # Numerically OK for the small (N, T, J) used here. Uses the same y-array
  # layout the user supplied (3D R array y[i, j, t]) directly, NOT the
  # flattened y_flat — this lets the test catch any indexing mismatch in
  # `.loglik_dyn_occu` against y_flat's column-major flatten.
  ref_loglik <- function(beta, y_array, n_sites, n_seasons, max_visits) {
    psi1 <- plogis(beta[1])
    p    <- plogis(beta[2])
    gam  <- plogis(beta[3])
    eps  <- plogis(beta[4])
    ll <- 0
    for (i in seq_len(n_sites)) {
      a_occ <- psi1; a_un <- 1 - psi1
      site_ll <- 0
      for (t in seq_len(n_seasons)) {
        raw <- y_array[i, , t]
        valid <- !is.na(raw) & raw >= 0
        nv_it <- sum(valid)
        if (nv_it > 0L) {
          y_valid <- raw[valid]
          n_det_t <- sum(y_valid == 1L)
          n_obs_t <- nv_it
          any_det <- n_det_t > 0L
          prob_y_occ <- p^n_det_t * (1 - p)^(n_obs_t - n_det_t)
          if (any_det) {
            # z must be 1
            site_ll <- site_ll + log(a_occ * prob_y_occ)
            a_occ <- 1; a_un <- 0
          } else {
            mix <- a_occ * prob_y_occ + a_un * 1
            site_ll <- site_ll + log(mix)
            a_occ <- (a_occ * prob_y_occ) / mix
            a_un  <- (a_un * 1) / mix
          }
        }
        if (t < n_seasons) {
          new_occ <- a_occ * (1 - eps) + a_un * gam
          new_un  <- a_occ * eps + a_un * (1 - gam)
          a_occ <- new_occ; a_un <- new_un
        }
      }
      ll <- ll + site_ll
    }
    ll
  }

  # Evaluate both at a grid of beta values inside the well-behaved logit
  # range. They must agree to floating-point precision (modulo log-sum-exp
  # vs linear differences in the renormalization step, which are O(1e-12)).
  set.seed(303)
  for (k in seq_len(5)) {
    beta_try <- runif(4, -1.5, 1.5)
    ll_a <- tulpaObs:::.loglik_dyn_occu(beta_try, model)
    ll_b <- ref_loglik(beta_try, y, N, n_seasons, J)
    expect_equal(ll_a, ll_b, tolerance = 1e-8,
                 label = sprintf("loglik vs reference at trial %d (beta=%s)",
                                 k, paste(round(beta_try, 3), collapse = ",")))
  }
})


# ---------------------------------------------------------------------------
# Sigma sanity: gamma is finite on a moderate dyn_occu simulation
# (N = 80, T = 3, J = 4). This is the spec's "Sigma sanity test".
# ---------------------------------------------------------------------------
test_that("SLA gamma is finite on N=80, T=3, J=4 dyn_occu", {
  skip_if_fast()
  set.seed(203)
  N <- 80; J <- 4; n_seasons <- 3
  elev <- rnorm(N)
  data <- data.frame(elev = elev)
  beta_psi1 <- c(0.2, 0.6)   # intercept + elev effect
  beta_p    <- 0.1
  beta_gam  <- -1.2
  beta_eps  <- -0.8

  psi1 <- plogis(beta_psi1[1] + beta_psi1[2] * elev)
  p    <- plogis(beta_p)
  gam  <- plogis(beta_gam)
  eps  <- plogis(beta_eps)
  z <- matrix(0L, N, n_seasons)
  z[, 1] <- rbinom(N, 1, psi1)
  for (t in 2:n_seasons) {
    z[, t] <- z[, t - 1] * (1 - rbinom(N, 1, eps)) +
              (1 - z[, t - 1]) * rbinom(N, 1, gam)
  }
  y <- array(NA_integer_, dim = c(N, J, n_seasons))
  for (i in seq_len(N)) {
    for (t in seq_len(n_seasons)) {
      y[i, , t] <- if (z[i, t] == 1L) rbinom(J, 1, p) else 0L
    }
  }

  fit <- tobs(
    formula     = ~ elev,
    data        = data,
    family      = dyn_occu(),
    detection   = ~ 1,
    y           = y,
    col_formula = ~ 1,
    ext_formula = ~ 1,
    method      = "laplace_sla",
    control     = list(verbose = FALSE)
  )

  expect_identical(fit$sla_status, "simplified_laplace")
  expect_true(is.numeric(fit$skew))
  expect_equal(length(fit$skew), fit$n_params)
  expect_true(all(is.finite(fit$skew)))
  # Each per-coefficient skew should be on a sensible scale; > 5 in magnitude
  # would indicate a Sigma-block scaling bug.
  expect_true(all(abs(fit$skew) < 5))
})


# ---------------------------------------------------------------------------
# Helper: the dyn_occu EM-Laplace log-marginal-likelihood reported by the
# fitter is a *surrogate* (the pseudo-binomial-encoded M-step likelihood, not
# the marginal HMM forward), so we don't expect exact agreement. The
# parent-spec's "cross-check 4" requires order-of-magnitude agreement only;
# we capture that here as a smoke check that the R-side HMM evaluator gives
# a *reasonable* number at the converged mode (negative, finite, not crazy).
# ---------------------------------------------------------------------------
test_that(".loglik_dyn_occu returns a reasonable number at the EM mode", {
  skip_if_fast()
  set.seed(204)
  N <- 50; J <- 3; n_seasons <- 3
  data <- data.frame(x = rnorm(N))
  psi1 <- 0.5; p <- 0.45; gam <- 0.2; eps <- 0.15
  z <- matrix(0L, N, n_seasons)
  z[, 1] <- rbinom(N, 1, psi1)
  for (t in 2:n_seasons) {
    z[, t] <- z[, t - 1] * (1 - rbinom(N, 1, eps)) +
              (1 - z[, t - 1]) * rbinom(N, 1, gam)
  }
  y <- array(NA_integer_, dim = c(N, J, n_seasons))
  for (i in seq_len(N)) {
    for (t in seq_len(n_seasons)) {
      y[i, , t] <- if (z[i, t] == 1L) rbinom(J, 1, p) else 0L
    }
  }

  fit <- tobs(
    formula     = ~ 1,
    data        = data,
    family      = dyn_occu(),
    detection   = ~ 1,
    y           = y,
    col_formula = ~ 1,
    ext_formula = ~ 1,
    method      = "laplace",
    control     = list(verbose = FALSE)
  )
  model <- tulpaObs:::.tobs_build_model(
    occ_formula = ~ 1, det_formula = ~ 1,
    data        = data, y = y,
    col_formula = ~ 1, ext_formula = ~ 1
  )
  ll <- tulpaObs:::.loglik_dyn_occu(as.numeric(fit$means), model)
  expect_true(is.finite(ll))
  # log P(y) <= 0 always; for N=50, T=3, J=3, with ~half occupied and
  # moderate detection, |ll| ~ O(100) (the data per site-season is at
  # worst log(2) ~ 0.7 nats per observation; total cap ~ 150 * log(2) ~ 100).
  expect_lt(ll, 0)
  expect_gt(ll, -1500)
})
