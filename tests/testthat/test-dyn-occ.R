test_that("dynamic occupancy model runs", {
  set.seed(42)
  n_sites <- 30
  n_seasons <- 3
  max_visits <- 3

  y_array <- array(0L, dim = c(n_sites, max_visits, n_seasons))
  z <- matrix(0, n_sites, n_seasons)
  z[, 1] <- rbinom(n_sites, 1, 0.6)
  for (t in 2:n_seasons) {
    z[, t] <- z[, t-1] * (1 - rbinom(n_sites, 1, 0.1)) +
              (1 - z[, t-1]) * rbinom(n_sites, 1, 0.2)
  }
  for (i in seq_len(n_sites))
    for (t in seq_len(n_seasons))
      if (z[i, t] == 1) y_array[i, , t] <- rbinom(max_visits, 1, 0.5)

  fit <- tobs(
    formula     = ~ 1,
    data        = data.frame(x = rnorm(n_sites)),
    family      = dyn_occu(),
    detection   = ~ 1,
    y           = y_array,
    col_formula = ~ 1,
    ext_formula = ~ 1,
    method      = "laplace",
    control     = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_equal(fit$n_params, 4)
  expect_true(fit$intercepts$psi1 > 0 && fit$intercepts$psi1 < 1)
})

test_that("dynamic tobs validates inputs", {
  expect_error(
    tobs(~ 1, data.frame(x = 1:5), family = dyn_occu(),
         detection = ~ 1, y = matrix(0, 3, 2),
         col_formula = ~ 1, ext_formula = ~ 1),
    "3D array"
  )
})


test_that("dyn_occu recovers (psi1, gamma, epsilon, p) within bias tolerance", {
  # Regression test for the y_flat indexing fix: prior to commit fixing
  # R/occu.R::.tobs_build_dynamic, `as.integer(y)` produced column-major
  # flat while src/dyn_occ_likelihood.h and build_dynamic_callbacks read
  # site-major — scrambling visits across (site, season) cells.
  # On N=200 / T=4 / J=4 with p_true=0.5, the buggy version recovered
  # p~0.32 (60% bias). After the fix, all four parameters recover within
  # ~0.05 of truth on a single seed; we assert < 0.1 to cover MC noise.
  set.seed(2026)
  n_sites <- 200; n_seasons <- 4; n_visits <- 4
  psi1 <- 0.5; gam <- 0.3; eps <- 0.2; p_true <- 0.5

  z <- matrix(NA_integer_, n_sites, n_seasons)
  z[, 1] <- rbinom(n_sites, 1, psi1)
  for (t in 2:n_seasons) {
    z[, t] <- ifelse(z[, t-1] == 1,
                     rbinom(n_sites, 1, 1 - eps),
                     rbinom(n_sites, 1, gam))
  }
  y <- array(0L, dim = c(n_sites, n_visits, n_seasons))
  for (i in seq_len(n_sites)) {
    for (t in seq_len(n_seasons)) {
      y[i, , t] <- if (z[i, t]) rbinom(n_visits, 1, p_true) else 0L
    }
  }

  fit <- tobs(
    formula     = ~ 1,
    data        = data.frame(idx = seq_len(n_sites)),
    family      = dyn_occu(),
    detection   = ~ 1,
    y           = y,
    col_formula = ~ 1,
    ext_formula = ~ 1,
    control     = list(verbose = FALSE)
  )

  rec <- c(
    psi1 = plogis(fit$means[["psi1_(Intercept)"]]),
    p    = plogis(fit$means[["p_(Intercept)"]]),
    gam  = plogis(fit$means[["gamma_(Intercept)"]]),
    eps  = plogis(fit$means[["epsilon_(Intercept)"]])
  )
  truth <- c(psi1 = psi1, p = p_true, gam = gam, eps = eps)
  expect_true(all(abs(rec - truth) < 0.1),
              info = sprintf("recovered: %s; truth: %s",
                             paste(round(rec, 3), collapse = ","),
                             paste(round(truth, 3), collapse = ",")))
})

test_that("fitted()$z is the forward-backward smoothed state for dynamic models", {
  # tulpaObs#18: dynamic fitted()$z must be the HMM smoothing posterior
  # P(z_t=1 | y_{1:T}), not the marginal occupancy psi_t. We verify (a) shape
  # [n_sites x n_seasons], (b) any detected (site, season) smooths to 1, and
  # (c) recovery against simulated truth beats the marginal-psi plug-in baseline.
  set.seed(7)
  n_sites <- 250; n_seasons <- 4; n_visits <- 4
  psi1 <- 0.5; gam <- 0.3; eps <- 0.2; p_true <- 0.5

  ztrue <- matrix(NA_integer_, n_sites, n_seasons)
  ztrue[, 1] <- rbinom(n_sites, 1, psi1)
  for (t in 2:n_seasons)
    ztrue[, t] <- ifelse(ztrue[, t-1] == 1, rbinom(n_sites, 1, 1 - eps),
                         rbinom(n_sites, 1, gam))
  y <- array(0L, dim = c(n_sites, n_visits, n_seasons))
  for (i in seq_len(n_sites)) for (t in seq_len(n_seasons))
    y[i, , t] <- if (ztrue[i, t]) rbinom(n_visits, 1, p_true) else 0L

  fit <- tobs(~ 1, data = data.frame(idx = seq_len(n_sites)), family = dyn_occu(),
              detection = ~ 1, y = y, col_formula = ~ 1, ext_formula = ~ 1,
              method = "laplace", control = list(verbose = FALSE))
  z <- fitted(fit)$z
  expect_true(is.matrix(z))
  expect_equal(dim(z), c(n_sites, n_seasons))
  expect_true(all(z >= 0 & z <= 1))

  # Any season with a detection is certainly occupied -> smoothed z == 1.
  det_mask <- apply(y, c(1, 3), function(v) any(v == 1))
  expect_true(all(abs(z[det_mask] - 1) < 1e-6))

  # Smoothed posterior beats the marginal-psi plug-in (the old behavior).
  psi1_hat <- plogis(fit$means[["psi1_(Intercept)"]])
  gam_hat  <- plogis(fit$means[["gamma_(Intercept)"]])
  eps_hat  <- plogis(fit$means[["epsilon_(Intercept)"]])
  psi_marg <- numeric(n_seasons); psi_marg[1] <- psi1_hat
  for (t in 2:n_seasons)
    psi_marg[t] <- psi_marg[t-1] * (1 - eps_hat) + (1 - psi_marg[t-1]) * gam_hat
  zmarg <- matrix(rep(psi_marg, each = n_sites), n_sites, n_seasons)

  brier_sm   <- mean((z     - ztrue)^2)
  brier_marg <- mean((zmarg - ztrue)^2)
  expect_lt(brier_sm, brier_marg)              # strictly better calibration
  expect_gt(mean((z > 0.5) == (ztrue == 1)), 0.9)  # high classification acc
})
