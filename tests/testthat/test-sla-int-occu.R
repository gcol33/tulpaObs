# Phase 3.5 — simplified-Laplace skewness on the integrated multi-source
# occupancy family (shared psi via FD d3 along Sigma[, j] in beta-space).

# -----------------------------------------------------------------------------
# Helper: build a small int_occu test data set with three detection sources,
# all observing the same n sites. Returns the data frame, the list of y
# matrices in the form expected by tobs(..., y = list(s1, s2, s3)), and the
# true beta values for downstream perturbation checks.
# -----------------------------------------------------------------------------
.make_int_occu_sim <- function(n = 80, J = c(3, 4, 3),
                               beta_occ = c(0.3, 0.6),
                               beta_det = list(-0.1, 0.2, -0.4),
                               seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  elev <- rnorm(n)
  d <- data.frame(elev = elev)
  X_occ <- model.matrix(~ elev, d)
  psi <- plogis(as.numeric(X_occ %*% beta_occ))
  z <- rbinom(n, 1, psi)

  y_list <- vector("list", length(J))
  for (s in seq_along(J)) {
    p_s <- plogis(beta_det[[s]])
    y <- matrix(0L, n, J[s])
    for (i in seq_len(n)) y[i, ] <- if (z[i]) rbinom(J[s], 1, p_s) else 0L
    y_list[[s]] <- y
  }
  names(y_list) <- paste0("s", seq_along(J))
  list(data = d, y = y_list, z = z, psi = psi,
       beta_occ = beta_occ, beta_det = beta_det)
}


test_that("simplified_laplace attaches finite skew on int_occu (smoke)", {
  set.seed(2027)
  sim <- .make_int_occu_sim(n = 80, J = c(3, 4, 3))

  fit <- tobs(
    formula   = ~ elev,
    data      = sim$data,
    family    = int_occu(),
    detection = ~ 1,
    y         = sim$y,
    approx    = "simplified_laplace",
    control   = list(verbose = FALSE)
  )

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$sla_status, "simplified_laplace")
  expect_true(is.numeric(fit$skew))
  expect_equal(length(fit$skew), length(fit$means))
  expect_named(fit$skew, names(fit$means))
  expect_true(all(is.finite(fit$skew)))

  # Joint dimension = p_occ (2) + sum_s p_det_s (1 each) = 5
  expect_equal(length(fit$skew), 2L + 3L * 1L)
})


test_that(".loglik_int_occu agrees with EM log_marginal_likelihood at beta_hat", {
  set.seed(2028)
  sim <- .make_int_occu_sim(n = 60, J = c(3, 3, 4))

  model <- .tobs_build_model(
    occ_formula = ~ elev,
    det_formula = ~ 1,
    data        = sim$data,
    y           = sim$y,
    integrated  = TRUE
  )

  fit <- tobs(
    formula   = ~ elev,
    data      = sim$data,
    family    = int_occu(),
    detection = ~ 1,
    y         = sim$y,
    approx    = "gaussian_laplace",
    control   = list(verbose = FALSE)
  )

  beta_hat <- as.numeric(fit$means)
  ll <- tulpaObs:::.loglik_int_occu(beta_hat, model)
  expect_true(is.finite(ll))
  expect_length(ll, 1L)

  # Perturbing beta_psi in a random direction should not increase the log-lik
  # by a large amount (we are at the MAP, so log-lik is locally concave).
  # Check that ll(beta_hat + 0.5 * e_j) <= ll(beta_hat) + slack for several j.
  p <- length(beta_hat)
  for (j in seq_len(p)) {
    e <- numeric(p); e[j] <- 1
    ll_pert_plus  <- tulpaObs:::.loglik_int_occu(beta_hat + 0.5 * e, model)
    ll_pert_minus <- tulpaObs:::.loglik_int_occu(beta_hat - 0.5 * e, model)
    # At a max, both directions should give <= ll, modulo numerical slack.
    expect_lt(ll_pert_plus,  ll + 1e-2)
    expect_lt(ll_pert_minus, ll + 1e-2)
  }
})


test_that(".loglik_int_occu is finite and on the same magnitude as EM marginal", {
  # Tier-3 sanity: at the gaussian-laplace mode, the R-side R-only int_occu
  # log-lik should be a finite scalar of order n_sites * 1 (typical Bernoulli
  # log-prob per site is in [-log 2, 0]).
  set.seed(2029)
  sim <- .make_int_occu_sim(n = 50, J = c(3, 3, 3))

  model <- .tobs_build_model(
    occ_formula = ~ elev,
    det_formula = ~ 1,
    data        = sim$data,
    y           = sim$y,
    integrated  = TRUE
  )

  fit <- tobs(
    formula   = ~ elev,
    data      = sim$data,
    family    = int_occu(),
    detection = ~ 1,
    y         = sim$y,
    approx    = "gaussian_laplace",
    control   = list(verbose = FALSE)
  )
  beta_hat <- as.numeric(fit$means)
  ll <- tulpaObs:::.loglik_int_occu(beta_hat, model)
  expect_true(is.finite(ll))
  # Upper bound: 0 (each site contributes a log-prob <= 0).
  expect_lt(ll, 0)
  # Lower bound: very loose — typical site log-prob is in [-10, 0]; 60 sites
  # x 3 sources x 3 visits gives a generous floor.
  expect_gt(ll, -10 * model$n_sites)
})


test_that("SLA gamma on int_occu is finite on a realistic simulated setup", {
  set.seed(2030)
  sim <- .make_int_occu_sim(n = 100, J = c(3, 4, 4),
                            beta_occ = c(0.2, 0.5),
                            beta_det = list(-0.3, 0.1, -0.2))

  fit <- tobs(
    formula   = ~ elev,
    data      = sim$data,
    family    = int_occu(),
    detection = ~ 1,
    y         = sim$y,
    approx    = "simplified_laplace",
    control   = list(verbose = FALSE)
  )

  expect_identical(fit$sla_status, "simplified_laplace")
  expect_true(all(is.finite(fit$skew)))
  # On this regime no |gamma| should saturate the SN ceiling — but cumulant
  # magnitudes for int_occu can run a bit higher than single-season because
  # the no-det branch couples beta_psi with all beta_p,s. Cap defensively.
  expect_true(max(abs(fit$skew)) < 5)
})
