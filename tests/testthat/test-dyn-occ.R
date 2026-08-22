test_that("dynamic occupancy model runs", {
  skip_if_fast()
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
    colonization = ~ 1,
    extinction = ~ 1,
    method      = "laplace",
    control     = list(verbose = FALSE)
  )
  expect_s3_class(fit, "tobs_fit")
  expect_equal(fit$n_params, 4)
  expect_true(fit$intercepts$psi1 > 0 && fit$intercepts$psi1 < 1)
})

test_that("dynamic tobs validates inputs", {
  skip_if_fast()
  expect_error(
    tobs(~ 1, data.frame(x = 1:5), family = dyn_occu(),
         detection = ~ 1, y = matrix(0, 3, 2),
         colonization = ~ 1, extinction = ~ 1),
    "3D array"
  )
})


test_that("dyn_occu recovers (psi1, gamma, epsilon, p) within bias tolerance", {
  skip_if_fast()
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
    colonization = ~ 1,
    extinction = ~ 1,
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
  skip_if_fast()
  # dynamic fitted()$z must be the HMM smoothing posterior P(z_t=1 | y_{1:T}),
  # not the marginal occupancy psi_t. We verify (a) shape [n_sites x n_seasons],
  # (b) any detected (site, season) smooths to 1, and (c) recovery against
  # simulated truth beats the marginal-psi plug-in baseline.
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
              detection = ~ 1, y = y, colonization = ~ 1, extinction = ~ 1,
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


# --- season-varying colonization / extinction ----------

test_that("a season-varying rate covariate routes to the interval path", {
  # A [n_sites x (T-1)] matrix covariate on colonization / extinction triggers
  # the interval-indexed design; a plain site-level covariate does NOT (it stays
  # on the byte-identical constant-rate path).
  sv <- simulate_dyn_occu(N = 40, J = 3, n_seasons = 5,
                          beta_gamma = c(-1, 0.8), seed = 1)
  m_sv <- tulpaObs:::.tobs_build_model(
    occ_formula = ~ 1, det_formula = ~ 1, data = sv$data, y = sv$y,
    col_formula = ~ gamma_cov, ext_formula = ~ 1)
  expect_true(m_sv$col_season_varying)
  expect_false(m_sv$ext_season_varying)
  # long-form design: (n_sites * (T-1)) rows
  expect_equal(nrow(m_sv$X_processes[[3]]), 40 * 4)
  expect_equal(nrow(m_sv$X_processes[[4]]), 40)          # extinction site-level

  # a site-level covariate is NOT season-varying
  d2 <- sv$data; d2$sitecov <- rnorm(40)
  m_c <- tulpaObs:::.tobs_build_model(
    occ_formula = ~ 1, det_formula = ~ 1, data = d2, y = sv$y,
    col_formula = ~ sitecov, ext_formula = ~ 1)
  expect_false(m_c$col_season_varying)
  expect_equal(nrow(m_c$X_processes[[3]]), 40)

  # a matrix covariate with the wrong number of columns errors
  d3 <- sv$data; d3$bad <- matrix(rnorm(40 * 2), 40, 2)
  expect_error(
    tulpaObs:::.tobs_build_model(occ_formula = ~ 1, det_formula = ~ 1,
      data = d3, y = sv$y, col_formula = ~ bad, ext_formula = ~ 1),
    "one column per transition interval")

  # season-varying rates are gated under NUTS (the C++ forward reads one rate
  # per site)
  expect_error(
    tobs(~ 1, data = sv$data, family = dyn_occu(), y = sv$y, detection = ~ 1,
         colonization = ~ gamma_cov, extinction = ~ 1, method = "nuts"),
    "season-varying")
})

test_that("dyn_occu recovers season-varying gamma/epsilon + ~95% coverage", {
  skip_if_fast()
  skip_on_cran()
  # Interval-indexed colonization / extinction driven by a per-(site, interval)
  # covariate (the dyn_abun #80 recipe ported to the colext forward). The
  # weighted-logistic transition M-step plus the season-varying exact HMM-forward
  # marginal refine recover the four transition coefficients with calibrated
  # Wald coverage.
  truth  <- c(-1.0, 0.9, -1.5, -0.7)   # gamma0, gamma1, eps0, eps1
  n_seed <- 20L
  est   <- matrix(NA_real_, n_seed, 4)
  cover <- matrix(FALSE, n_seed, 4)
  for (s in seq_len(n_seed)) {
    sv <- simulate_dyn_occu(N = 400, J = 5, n_seasons = 7,
                            beta_occ = 0.3, beta_det = 0.4,
                            beta_gamma = c(-1.0, 0.9),
                            beta_epsilon = c(-1.5, -0.7), seed = 300 + s)
    fit <- tobs(~ 1, data = sv$data, family = dyn_occu(), y = sv$y,
                detection = ~ 1, colonization = ~ gamma_cov,
                extinction = ~ eps_cov, method = "laplace",
                control = list(progress = FALSE, verbose = FALSE))
    cc <- unlist(coef(fit)); se <- sqrt(diag(vcov(fit)))
    idx <- c(which(names(cc) == "gamma.(Intercept)"),
             which(names(cc) == "gamma.gamma_cov"),
             which(names(cc) == "epsilon.(Intercept)"),
             which(names(cc) == "epsilon.eps_cov"))
    b <- cc[idx]; s_e <- se[idx]
    est[s, ]   <- b
    cover[s, ] <- (truth >= b - 1.96 * s_e) & (truth <= b + 1.96 * s_e)
  }
  expect_equal(colMeans(est), truth, tolerance = 0.08)
  expect_true(all(colMeans(cover) >= 0.85))
})

# --- season-varying detection --------------------------

test_that("a season-varying detection covariate routes to the season path", {
  # A [n_sites x T] matrix covariate on detection triggers the season-indexed
  # (T-column) design; a plain site-level covariate does NOT.
  sv <- simulate_dyn_occu(N = 40, J = 3, n_seasons = 5,
                          beta_det_season = c(0.2, -0.8), seed = 1)
  m_sv <- tulpaObs:::.tobs_build_model(
    occ_formula = ~ 1, det_formula = ~ det_cov, data = sv$data, y = sv$y,
    col_formula = ~ 1, ext_formula = ~ 1)
  expect_true(m_sv$det_season_varying)
  # long-form detection design: (n_sites * T) rows
  expect_equal(nrow(m_sv$X_processes[[2]]), 40 * 5)

  # a site-level detection covariate is NOT season-varying (byte-identical path)
  d2 <- sv$data; d2$sitecov <- rnorm(40)
  m_c <- tulpaObs:::.tobs_build_model(
    occ_formula = ~ 1, det_formula = ~ sitecov, data = d2, y = sv$y,
    col_formula = ~ 1, ext_formula = ~ 1)
  expect_false(isTRUE(m_c$det_season_varying))
  expect_equal(nrow(m_c$X_processes[[2]]), 40)

  # a matrix detection covariate with the wrong number of columns errors
  d3 <- sv$data; d3$bad <- matrix(rnorm(40 * 3), 40, 3)   # T = 5, so 3 is wrong
  expect_error(
    tulpaObs:::.tobs_build_model(occ_formula = ~ 1, det_formula = ~ bad,
      data = d3, y = sv$y, col_formula = ~ 1, ext_formula = ~ 1),
    "one column per primary season")

  # gated under NUTS (the C++ forward reads one detection predictor per site)
  expect_error(
    tobs(~ 1, data = sv$data, family = dyn_occu(), y = sv$y,
         detection = ~ det_cov, colonization = ~ 1, extinction = ~ 1,
         method = "nuts"),
    "season-varying")
})

test_that("dyn_occu recovers season-varying detection + ~95% coverage", {
  skip_if_fast()
  skip_on_cran()
  # Per-(site, season) detection driven by a [n_sites x T] covariate: the E-step
  # emission reads the season's detection probability and the M-step encodes one
  # detection row per (site, season) at the season's covariate, then the exact
  # season-varying HMM-forward marginal refine calibrates the coefficients.
  truth  <- c(0.2, -0.8)               # p intercept, p slope on det_cov
  n_seed <- 20L
  est   <- matrix(NA_real_, n_seed, 2)
  cover <- matrix(FALSE, n_seed, 2)
  for (s in seq_len(n_seed)) {
    sv <- simulate_dyn_occu(N = 350, J = 4, n_seasons = 6, beta_occ = 0.4,
                            beta_det_season = c(0.2, -0.8),
                            gamma = 0.25, epsilon = 0.15, seed = 400 + s)
    fit <- tobs(~ 1, data = sv$data, family = dyn_occu(), y = sv$y,
                detection = ~ det_cov, colonization = ~ 1, extinction = ~ 1,
                method = "laplace", control = list(progress = FALSE, verbose = FALSE))
    cc <- unlist(coef(fit)); se <- sqrt(diag(vcov(fit)))
    idx <- c(which(names(cc) == "p.(Intercept)"),
             which(names(cc) == "p.det_cov"))
    b <- cc[idx]; s_e <- se[idx]
    est[s, ]   <- b
    cover[s, ] <- (truth >= b - 1.96 * s_e) & (truth <= b + 1.96 * s_e)
  }
  expect_equal(colMeans(est), truth, tolerance = 0.08)
  expect_true(all(colMeans(cover) >= 0.85))
})
