# Tests for cover(positive = "lognormal") -- Phase 1a.

test_that("cover(positive = 'lognormal') flips to working", {
  fam <- cover(positive = "lognormal")
  expect_s3_class(fam, "tobs_family")
  expect_equal(fam$name, "cover")
  expect_equal(fam$status, "working")
  expect_equal(fam$default_engine, "laplace")
  expect_equal(fam$params$positive, "lognormal")
})

test_that("cover(positive = 'beta') is wired through to the beta Laplace engine", {
  fam <- cover("beta")
  expect_equal(fam$status, "working")
  expect_equal(fam$params$positive, "beta")

  sim <- simulate_cover(N = 300, seed = 1)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("beta"),
    y       = sim$y
  )
  expect_s3_class(fit, "cover_fit")
  expect_equal(fit$positive, "beta")
  expect_true(is.finite(fit$phi_pos) && fit$phi_pos > 0)
  expect_true(is.na(fit$sigma_pos))
})

test_that("simulator round-trips a coefficient prior", {
  sim <- simulate_cover(N = 500, seed = 42)
  expect_named(sim, c("data", "y", "coords", "truth"))
  expect_equal(length(sim$y), 500)
  expect_true(all(sim$y >= 0 & sim$y <= 1))
  # Some sites should be zero, some positive.
  expect_true(mean(sim$truth$occur) > 0.1)
  expect_true(mean(sim$truth$occur) < 0.95)
})

test_that("single fit recovers truth within tolerance and prediction identity holds", {
  sim <- simulate_cover(
    N         = 800,
    beta_occ  = c(-0.5, 0.8),
    beta_pos  = c(-1.0, 0.3),
    sigma_pos = 0.4,
    seed      = 2026
  )
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover(positive = "lognormal"),
    y       = sim$y
  )
  expect_s3_class(fit, "cover_fit")
  expect_true(fit$converged)

  expect_lt(abs(fit$beta_occ[1] - sim$truth$beta_occ[1]), 0.4)
  expect_lt(abs(fit$beta_occ[2] - sim$truth$beta_occ[2]), 0.4)
  expect_lt(abs(fit$beta_pos[1] - sim$truth$beta_pos[1]), 0.2)
  expect_lt(abs(fit$beta_pos[2] - sim$truth$beta_pos[2]), 0.2)
  expect_lt(abs(fit$sigma_pos / sim$truth$sigma_pos - 1), 0.2)

  newdata <- data.frame(x = sim$data$x)
  p_hat  <- predict(fit, newdata, type = "occupancy")
  mu_hat <- predict(fit, newdata, type = "conditional")
  e_hat  <- predict(fit, newdata, type = "expected")
  expect_equal(e_hat, p_hat * mu_hat, tolerance = 1e-8)

  z_occ <- (fit$beta_occ - sim$truth$beta_occ) / fit$se_occ
  z_pos <- (fit$beta_pos - sim$truth$beta_pos) / fit$se_pos
  expect_true(all(abs(z_occ) < 4))
  expect_true(all(abs(z_pos) < 4))
})

test_that("Gaussian arm uses only positive-cover rows", {
  sim <- simulate_cover(N = 300, seed = 7)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover(positive = "lognormal"),
    y       = sim$y
  )
  n_pos_obs <- sum(sim$y > 0)
  expect_equal(fit$n_positive, n_pos_obs)
  expect_equal(length(fit$encoding$pos_data$y), n_pos_obs)
})

test_that("repeat fits recover truth in aggregate (light sanity, 10 reps)", {
  skip_on_cran()
  truth <- list(beta_occ = c(-0.5, 0.8), beta_pos = c(-1.0, 0.3),
                sigma_pos = 0.4)
  hits_occ <- integer(2)
  hits_pos <- integer(2)
  sigma_diffs <- numeric(10)
  for (r in seq_len(10)) {
    sim <- simulate_cover(
      N         = 600,
      beta_occ  = truth$beta_occ,
      beta_pos  = truth$beta_pos,
      sigma_pos = truth$sigma_pos,
      seed      = 100 + r
    )
    fit <- tobs(
      formula = ~ x,
      data    = sim$data,
      family  = cover("lognormal"),
      y       = sim$y
    )
    hits_occ <- hits_occ +
      as.integer(abs(fit$beta_occ - truth$beta_occ) <= 2 * fit$se_occ)
    hits_pos <- hits_pos +
      as.integer(abs(fit$beta_pos - truth$beta_pos) <= 2 * fit$se_pos)
    sigma_diffs[r] <- abs(fit$sigma_pos - truth$sigma_pos) / truth$sigma_pos
  }
  expect_gte(min(hits_occ), 7)
  expect_gte(min(hits_pos), 7)
  expect_lt(mean(sigma_diffs), 0.15)
})

test_that("predict requires newdata with the same columns as the formula", {
  sim <- simulate_cover(N = 200, seed = 11)
  fit <- tobs(
    formula = ~ x,
    data    = sim$data,
    family  = cover("lognormal"),
    y       = sim$y
  )
  expect_error(
    predict(fit, newdata = data.frame(z = rnorm(5))),
    "object 'x' not found",
    fixed = FALSE
  )
})

test_that("a temporal() term requires engine = 'nested_laplace'", {
  sim <- simulate_cover(N = 50, seed = 13)
  sim$data$year <- sample.int(3, nrow(sim$data), replace = TRUE)
  expect_error(
    tobs(
      formula  = ~ x + temporal(year, type = "ar1"),
      data     = sim$data,
      family   = cover("lognormal"),
      y        = sim$y,
      engine   = "laplace"
    ),
    "require engine = 'nested_laplace'"
  )
})

test_that("detection = errors (cover has no detection layer)", {
  sim <- simulate_cover(N = 50, seed = 17)
  expect_error(
    tobs(
      formula   = ~ x,
      data      = sim$data,
      family    = cover("lognormal"),
      y         = sim$y,
      detection = ~ 1
    ),
    "does not use a detection formula"
  )
})
