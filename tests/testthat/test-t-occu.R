# =============================================================================
# test-t-occu.R - multi-season occupancy with an AR1 year random effect
# (t_occu(); spOccupancy tPGOcc; gcol33/tulpaObs#124).
#
# Distinct from dyn_occu() (colext): NO colonization / extinction transition --
# per-(site, season) occupancy is a Bernoulli GLMM with a shared AR1 year effect
# on the logit, fit by an exact Polya-Gamma Gibbs sampler. The primary
# deliverables are the year-effect surface eta_t, the occupancy / detection
# coefficients, and the AR1 innovation SD sigma. The AR1 CORRELATION rho is a
# weakly-identified hyperparameter of a short time series -- given even the TRUE
# eta series it recovers only as the number of seasons grows (rho_hat ~0.08 at
# T=8, ~0.43 at T=20, ~0.59 at T=200 for truth 0.6; dev_notes/_probe_tpg_rho.R),
# so these assert on the surface + coefficients + sigma tightly and on rho only
# directionally, matching how spOccupancy documents the AR1 temporal parameters.
# =============================================================================

test_that("t_occu() gates + S3 surface", {
  sim <- simulate_t_occu(N = 80, T_seasons = 6, J = 3, beta_occ = c(0.2),
                         p = 0.4, rho = 0.6, sigma = 0.7, seed = 1)

  # Requires a detection formula and y.
  expect_error(tobs(~ 1, family = t_occu(), y = sim$y, data = sim$data,
                    method = "pg_gibbs"),
               "detection")
  expect_error(tobs(~ 1, family = t_occu(), detection = ~ 1, data = sim$data,
                    method = "pg_gibbs"),
               "y")
  # pg_gibbs only.
  expect_error(tobs(~ 1, family = t_occu(), detection = ~ 1, y = sim$y,
                    data = sim$data, method = "laplace"),
               "pg_gibbs")
  # >= 2 seasons.
  sim1 <- simulate_t_occu(N = 40, T_seasons = 2, J = 3, seed = 1)
  y1 <- sim1$y[, 1L, , drop = FALSE]
  expect_error(.tobs_build_t_occu(~ 1, ~ 1, sim1$data,
                                  array(y1, dim = c(40L, 1L, 3L))),
               ">= 2 seasons")

  fit <- tobs(~ 1, family = t_occu(), detection = ~ 1, y = sim$y,
              data = sim$data, method = "pg_gibbs",
              control = list(n.iter = 800L, n.warmup = 400L, n.chains = 2L,
                             seed = 1))
  expect_s3_class(fit, "tobs_fit")
  expect_true(all(c("psi_(Intercept)", "p_(Intercept)", "log_sigma_ar1",
                    "rho_ar1") %in% names(fit$means)))
  expect_length(fit$temporal_field, 6L)          # one year effect per season
  expect_true(is.finite(coef(fit)$psi[["(Intercept)"]]))
  expect_true(all(is.finite(fit$rhat)))
})

test_that("t_occu() accepts a list of per-season matrices", {
  sim <- simulate_t_occu(N = 50, T_seasons = 5, J = 3, seed = 2)
  ylist <- lapply(seq_len(5L), function(t) sim$y[, t, ])
  fit <- tobs(~ 1, family = t_occu(), detection = ~ 1, y = ylist,
              data = sim$data, method = "pg_gibbs",
              control = list(n.iter = 600L, n.warmup = 300L, n.chains = 2L,
                             seed = 1))
  expect_s3_class(fit, "tobs_fit")
  expect_length(fit$temporal_field, 5L)
})

test_that("t_occu() recovers the year-effect surface, coefficients + sigma", {
  skip_if_fast()
  skip_on_cran()

  n_seed <- 12L
  eta_cor <- psi0 <- p_hat <- sigma_hat <- rho_hat <- max_rhat <-
    numeric(n_seed)
  psi0_cov <- logical(n_seed)
  b0 <- 0.2

  for (s in seq_len(n_seed)) {
    sim <- simulate_t_occu(N = 220, T_seasons = 12, J = 4, beta_occ = c(b0),
                           p = 0.4, rho = 0.6, sigma = 0.7, seed = 100 + s)
    fit <- tobs(~ 1, family = t_occu(), detection = ~ 1, y = sim$y,
                data = sim$data, method = "pg_gibbs",
                control = list(n.iter = 2000L, n.warmup = 1000L,
                               n.chains = 2L, seed = 7L))
    m <- fit$means
    eta_cor[s]   <- suppressWarnings(cor(fit$temporal_field, sim$truth$eta))
    psi0[s]      <- m[["psi_(Intercept)"]]
    p_hat[s]     <- stats::plogis(m[["p_(Intercept)"]])
    sigma_hat[s] <- exp(m[["log_sigma_ar1"]])
    rho_hat[s]   <- m[["rho_ar1"]]
    max_rhat[s]  <- max(fit$rhat, na.rm = TRUE)
    ci <- confint(fit)["psi_(Intercept)", ]
    psi0_cov[s]  <- ci[1L] <= b0 && b0 <= ci[2L]
  }

  # Year-effect surface = the primary output: recovers essentially exactly.
  expect_gt(median(eta_cor), 0.90)
  expect_gt(min(eta_cor), 0.80)

  # Occupancy + detection intercepts unbiased; nominal interval coverage.
  expect_lt(abs(median(psi0) - b0), 0.12)
  expect_lt(abs(median(p_hat) - 0.4), 0.05)
  expect_gte(mean(psi0_cov), 0.83)

  # AR1 innovation SD recovers (mild small-T inflation tolerated).
  expect_lt(abs(median(sigma_hat) - 0.7), 0.25)

  # AR1 correlation is weakly identified at T=12 -- assert only that it is
  # positive on average (the surface, not rho, is the deliverable).
  expect_gt(median(rho_hat), 0.15)

  # Converged.
  expect_lt(median(max_rhat), 1.1)
})
