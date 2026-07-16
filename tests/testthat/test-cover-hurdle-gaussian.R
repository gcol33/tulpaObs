# Tests for cover(response = "gaussian") -- the identity-Gaussian (delta-normal)
# positive arm (gcol33/tulpaObs#112). It is the lognormal arm on the raw response
# (no log transform, no Jacobian), for a magnitude that lives on a real,
# unbounded scale. Presence is the nonzero sentinel (y != 0), not y > 0.

test_that("cover(response = 'gaussian') constructor is wired through", {
  fam <- cover(response = "gaussian")
  expect_s3_class(fam, "tobs_family")
  expect_equal(fam$name, "cover")
  expect_equal(fam$status, "working")
  expect_equal(fam$default_engine, "laplace")
  expect_equal(fam$params$positive, "gaussian")
  expect_equal(fam$observation, "binomial_plus_gaussian")
})

test_that("simulate_cover(response = 'gaussian') round-trips an unbounded response", {
  sim <- simulate_cover(N = 500, beta_pos = c(2.0, 0.3), sigma_pos = 0.5,
                        response = "gaussian", seed = 42)
  expect_named(sim, c("data", "y", "coords", "truth"))
  expect_equal(sim$truth$response, "gaussian")
  # Absent sites are the exact 0 sentinel; present magnitudes are unbounded.
  expect_true(mean(sim$truth$occur) > 0.1 && mean(sim$truth$occur) < 0.95)
  present <- sim$y[sim$truth$occur == 1L]
  expect_true(any(present != 0))
  # mu on the response scale is eta (identity), not exp(eta + sigma^2/2).
  expect_equal(sim$truth$mu, as.numeric(cbind(1, sim$data$x) %*% sim$truth$beta_pos))
})

test_that("single fit recovers truth (gaussian arm, identity mean)", {
  sim <- simulate_cover(
    N         = 800,
    beta_occ  = c(-0.3, 0.7),
    beta_pos  = c(2.0, 0.4),
    sigma_pos = 0.5,
    response  = "gaussian",
    seed      = 2026
  )
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover(response = "gaussian"),
    y       = sim$y
  )
  expect_s3_class(fit, "cover_fit")
  expect_true(fit$converged)
  expect_equal(fit$positive, "gaussian")
  expect_true(is.na(fit$phi_pos))
  expect_true(is.finite(fit$sigma_pos) && fit$sigma_pos > 0)

  expect_lt(abs(fit$beta_occ[1] - sim$truth$beta_occ[1]), 0.4)
  expect_lt(abs(fit$beta_occ[2] - sim$truth$beta_occ[2]), 0.4)
  expect_lt(abs(fit$beta_pos[1] - sim$truth$beta_pos[1]), 0.2)
  expect_lt(abs(fit$beta_pos[2] - sim$truth$beta_pos[2]), 0.2)
  expect_lt(abs(fit$sigma_pos / sim$truth$sigma_pos - 1), 0.2)

  # predict() on the response scale: conditional mean is the linear predictor.
  newdata <- data.frame(x = sim$data$x)
  p_hat  <- predict(fit, newdata, type = "occupancy")
  mu_hat <- predict(fit, newdata, type = "conditional")
  e_hat  <- predict(fit, newdata, type = "expected")
  expect_equal(e_hat, p_hat * mu_hat, tolerance = 1e-8)
  eta_hat <- as.numeric(cbind(1, sim$data$x) %*% fit$beta_pos)
  expect_equal(as.numeric(mu_hat), eta_hat, tolerance = 1e-8)
})

test_that("gaussian presence is y != 0 (negative magnitudes are present)", {
  # A present site with a negative magnitude must be counted as positive, not
  # collapsed into the absent 0 sentinel.
  x <- rnorm(200)
  set.seed(3)
  occur <- rbinom(200, 1L, plogis(0.2 + 0.5 * x))
  mag   <- rnorm(200, -1.0 + 0.3 * x, 0.5)          # centred below 0
  y     <- ifelse(occur == 1L, mag, 0)
  fit <- tobs(formula = ~ x, data = data.frame(x = x),
              family = cover(response = "gaussian"), y = y)
  expect_equal(fit$n_positive, sum(occur == 1L & mag != 0))
  # A lognormal arm would reject y <= 0; the gaussian arm accepts it.
  expect_true(is.finite(fit$sigma_pos))
})

test_that("cover(gaussian): WAIC works, PPC and NUTS are gated", {
  sim <- simulate_cover(N = 400, beta_pos = c(2.0, 0.4), sigma_pos = 0.5,
                        response = "gaussian", seed = 21)
  fit <- tobs(formula = ~ x, data = sim$data,
              family = cover(response = "gaussian"), y = sim$y)
  w <- tobs_waic(fit)
  expect_true(is.finite(w$waic) && is.finite(w$p_waic))
  expect_error(tobs_ppc(fit), "not defined for cover.*gaussian")
  expect_error(
    tobs(formula = ~ x, data = sim$data, family = cover(response = "gaussian"),
         y = sim$y, method = "nuts"),
    "not yet wired for method = 'nuts'")
})


test_that("repeat fits recover truth in aggregate (gaussian, 20 seeds)", {
  skip_on_cran()
  skip_if_fast()
  truth <- list(beta_occ = c(-0.3, 0.7), beta_pos = c(2.0, 0.4), sigma_pos = 0.5)
  n_seeds <- 20L
  est_occ <- matrix(NA_real_, n_seeds, 2L)
  est_pos <- matrix(NA_real_, n_seeds, 2L)
  se_occ  <- matrix(NA_real_, n_seeds, 2L)
  se_pos  <- matrix(NA_real_, n_seeds, 2L)
  sigma_diffs <- numeric(n_seeds)
  conv <- logical(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- simulate_cover(N = 600, beta_occ = truth$beta_occ,
                          beta_pos = truth$beta_pos, sigma_pos = truth$sigma_pos,
                          response = "gaussian", seed = 500L + r)
    fit <- tryCatch(
      tobs(formula = ~ x, data = sim$data,
           family = cover(response = "gaussian"), y = sim$y),
      error = function(e) NULL)
    if (is.null(fit)) next
    conv[r] <- isTRUE(fit$converged)
    est_occ[r, ] <- fit$beta_occ; se_occ[r, ] <- fit$se_occ
    est_pos[r, ] <- fit$beta_pos; se_pos[r, ] <- fit$se_pos
    sigma_diffs[r] <- abs(fit$sigma_pos - truth$sigma_pos) / truth$sigma_pos
  }
  ok <- conv & is.finite(est_pos[, 1L])
  expect_gte(mean(ok), 0.80)

  # Point recovery (bias below tolerance).
  expect_lt(abs(mean(est_occ[ok, 1L]) - truth$beta_occ[1L]), 0.20)
  expect_lt(abs(mean(est_occ[ok, 2L]) - truth$beta_occ[2L]), 0.20)
  expect_lt(abs(mean(est_pos[ok, 1L]) - truth$beta_pos[1L]), 0.10)
  expect_lt(abs(mean(est_pos[ok, 2L]) - truth$beta_pos[2L]), 0.10)
  expect_lt(mean(sigma_diffs[ok]), 0.10)

  # 95% Wald CI pooled coverage at the 0.85 working-family floor.
  cov_cells <- c(
    abs(est_occ[ok, 1L] - truth$beta_occ[1L]) < 1.96 * se_occ[ok, 1L],
    abs(est_occ[ok, 2L] - truth$beta_occ[2L]) < 1.96 * se_occ[ok, 2L],
    abs(est_pos[ok, 1L] - truth$beta_pos[1L]) < 1.96 * se_pos[ok, 1L],
    abs(est_pos[ok, 2L] - truth$beta_pos[2L]) < 1.96 * se_pos[ok, 2L]
  )
  expect_gte(mean(cov_cells), 0.85)
})
