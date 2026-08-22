# =============================================================================
# test-occu-multi.R - multi-species co-occurrence occupancy (occu_multi();
# Rota et al. 2016; unmarked occuMulti).
#
# Joint occupancy z in {0,1}^S from a log-linear model with first-order (per
# species) and second-order (per pair) natural parameters; the interaction is
# the second-order term (positive = co-occur, negative = avoid). The latent
# state integrates out by enumerating the 2^S states; the exact marginal is
# maximised (optim BFGS) with an observed-information vcov. Non-spatial laplace.
#
# The log-linear natural-parameter slopes trade off (only the MARGINAL occupancy
# and the interaction are cleanly identified), so recovery targets the
# interaction natural parameter, the per-species marginal occupancy, and the
# detection -- not the individual first-order covariate slopes.
# =============================================================================

test_that("occu_multi() constructor + gates", {
  f <- occu_multi()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "occu_multi")
  sim <- simulate_occu_multi(S = 2, N = 40, J = 3, seed = 1)
  expect_error(
    tobs(~ scov1, data = sim$data, family = occu_multi(), detection = ~ 1,
         y = sim$y, species = sim$species, method = "nuts"),
    "laplace")
})

test_that("occu_multi() fits + full S3 surface (S = 2)", {
  sim <- simulate_occu_multi(S = 2, N = 300, J = 5, seed = 3)
  fit <- tobs(~ scov1, data = sim$data, family = occu_multi(), detection = ~ 1,
              y = sim$y, species = sim$species,
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true(isTRUE(fit$convergence$converged))
  expect_true("f_sp1_sp2_(Intercept)" %in% names(fit$means))   # the interaction

  fv <- fitted(fit)
  expect_named(fv, c("psi", "p"))
  expect_equal(dim(fv$psi), c(300L, 2L))
  expect_true(all(fv$psi > 0 & fv$psi < 1))
  expect_equal(predict(fit, type = "state"), fv$psi)
  w <- waic(fit, n.draws = 100L)
  expect_true(is.finite(w$waic))
  s2 <- simulate(fit, nsim = 1)
  expect_length(s2, 2L)
  expect_equal(dim(residuals(fit)$occ), c(300L, 2L))

  # nobs() counts every surveyed (site, visit) cell, over all species.
  expect_identical(nobs(fit),
                   sum(vapply(sim$y, function(m) sum(!is.na(m)), integer(1))))
})

test_that("occu_multi() recovers the interaction sign + marginal occupancy", {
  skip_on_cran()
  skip_if_fast()
  # Intercept-only natural parameters, so f1 / f2 / f12 are cleanly identified.
  n_seed <- 15L
  f12_pos <- f12_neg <- rep(NA_real_, n_seed)
  # Fitted-vs-realized marginal occupancy gap per species (the log-linear model
  # inflates the marginal above plogis(first-order) via the interaction, so the
  # target is the REALIZED colMeans(z), not plogis(natural param)).
  psi_gap1 <- psi_gap2 <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    # Positive interaction (co-occurrence).
    simp <- simulate_occu_multi(
      S = 2, N = 350, J = 5, n_state_covs = 0,
      beta_first = list(c(stats::qlogis(0.45)), c(stats::qlogis(0.5))),
      beta_second = list(c(1.2)),
      beta_p = list(c(stats::qlogis(0.55)), c(stats::qlogis(0.55))),
      seed = 800 + s)
    fp <- tryCatch(tobs(~ 1, data = simp$data, family = occu_multi(),
                        detection = ~ 1, y = simp$y, species = simp$species,
                        control = list(verbose = FALSE, progress = FALSE)),
                   error = function(e) NULL)
    if (!is.null(fp) && isTRUE(fp$convergence$converged)) {
      f12_pos[s] <- fp$means[["f_sp1_sp2_(Intercept)"]]
      fv <- fitted(fp)
      psi_gap1[s] <- mean(fv$psi[, 1]) - mean(simp$truth$z[, 1])
      psi_gap2[s] <- mean(fv$psi[, 2]) - mean(simp$truth$z[, 2])
    }
    # Negative interaction (avoidance).
    simn <- simulate_occu_multi(
      S = 2, N = 350, J = 5, n_state_covs = 0,
      beta_first = list(c(stats::qlogis(0.5)), c(stats::qlogis(0.5))),
      beta_second = list(c(-1.2)),
      beta_p = list(c(stats::qlogis(0.55)), c(stats::qlogis(0.55))),
      seed = 8000 + s)
    fn <- tryCatch(tobs(~ 1, data = simn$data, family = occu_multi(),
                        detection = ~ 1, y = simn$y, species = simn$species,
                        control = list(verbose = FALSE, progress = FALSE)),
                   error = function(e) NULL)
    if (!is.null(fn) && isTRUE(fn$convergence$converged))
      f12_neg[s] <- fn$means[["f_sp1_sp2_(Intercept)"]]
  }
  # The interaction natural parameter recovers (mean near truth) and never
  # flips sign in aggregate.
  expect_lt(abs(mean(f12_pos, na.rm = TRUE) - 1.2), 0.3)
  expect_gt(mean(f12_pos, na.rm = TRUE), 0.5)
  expect_lt(abs(mean(f12_neg, na.rm = TRUE) - (-1.2)), 0.3)
  expect_lt(mean(f12_neg, na.rm = TRUE), -0.5)
  # Fitted marginal per-species occupancy tracks the realized occupancy rate.
  expect_lt(abs(mean(psi_gap1, na.rm = TRUE)), 0.05)
  expect_lt(abs(mean(psi_gap2, na.rm = TRUE)), 0.05)
})
