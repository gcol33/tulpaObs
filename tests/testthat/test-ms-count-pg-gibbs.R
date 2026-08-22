# =============================================================================
# test-ms-count-pg-gibbs.R - Polya-Gamma Gibbs for the community Bernoulli /
# binomial GLMM (jsdm() and ms_count(response="binomial") method = "pg_gibbs";
# spOccupancy msPGOcc-family / community svcPGBinom).
#
# The community logistic GLMM has no latent state, so the PG sampler is the pure
# per-species conjugate coefficient update + conjugate community mean +
# Inverse-Gamma community variance. Samples the exact community posterior, so the
# community variance recovers where the Laplace-EM attenuates. Only the logistic
# responses (bernoulli / binomial) are routed here.
# =============================================================================

test_that("jsdm()/ms_count() method='pg_gibbs' gates + S3", {
  sim <- simulate_jsdm(N = 60, n_species = 8, seed = 1)
  fit <- tobs(~ x, data = sim$data, family = jsdm(), y = sim$y,
              species = paste0("sp", seq_len(8)), method = "pg_gibbs",
              control = list(n.iter = 800L, n.warmup = 400L, n.chains = 2L,
                             seed = 1, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "pg_gibbs")
  expect_true(all(is.finite(fit$rhat)) && all(fit$rhat < 1.1))
  expect_true(all(c("sd_mu", "coef_mu") %in% names(fit$ms_community)))
  # A Poisson community count rejects pg_gibbs (not a logistic response).
  sp <- simulate_ms_count(N = 60, n_species = 6, response = "poisson", seed = 2)
  expect_error(
    tobs(~ x, data = sp$data, family = ms_count("poisson"), y = sp$y,
         species = paste0("sp", seq_len(6)), method = "pg_gibbs"),
    "logistic")
})

test_that("jsdm() pg_gibbs recovers community means + SD (not attenuated)", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 10L
  mu_t <- c(0.2, 0.7); sd_t <- c(0.7, 0.5)
  mu_e <- sd_e <- matrix(NA_real_, n_seed, 2L); cor_e <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_jsdm(N = 160, n_species = 20, beta_comm_mean = mu_t,
                         beta_comm_sd = sd_t, seed = 300 + s)
    fg <- tryCatch(
      tobs(~ x, data = sim$data, family = jsdm(), y = sim$y,
           species = paste0("sp", seq_len(20)), method = "pg_gibbs",
           control = list(n.iter = 1500L, n.warmup = 750L, n.chains = 2L,
                          seed = s, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fg)) next
    mu_e[s, ] <- fg$means
    sd_e[s, ] <- fg$ms_community$sd_mu
    cor_e[s]  <- min(diag(cor(fg$ms_community$coef_mu, sim$truth$beta_species)))
  }
  expect_lt(abs(mean(mu_e[, 1], na.rm = TRUE) - mu_t[1]), 0.10)
  expect_lt(abs(mean(mu_e[, 2], na.rm = TRUE) - mu_t[2]), 0.10)
  expect_lt(abs(mean(sd_e[, 1], na.rm = TRUE) - sd_t[1]), 0.10)
  expect_lt(abs(mean(sd_e[, 2], na.rm = TRUE) - sd_t[2]), 0.10)
  expect_gt(mean(cor_e, na.rm = TRUE), 0.85)
})

test_that("ms_count('binomial') pg_gibbs recovers the community slope", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 8L
  b1 <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_count(N = 140, n_species = 14, response = "binomial",
                             beta_comm_mean = c(0.3, 0.6), beta_comm_sd = c(0.6, 0.4),
                             trials = 6, seed = 500 + s)
    fb <- tryCatch(
      tobs(~ x, data = sim$data, family = ms_count("binomial"), y = sim$y,
           species = paste0("sp", seq_len(14)), trials = 6, method = "pg_gibbs",
           control = list(n.iter = 1200L, n.warmup = 600L, n.chains = 2L,
                          seed = s, verbose = FALSE)),
      error = function(e) NULL)
    if (!is.null(fb)) b1[s] <- fb$means[[2L]]
  }
  # The community slope recovers (the k-of-n likelihood identifies it cleanly).
  expect_lt(abs(mean(b1, na.rm = TRUE) - 0.6), 0.08)
})
