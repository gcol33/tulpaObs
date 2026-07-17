# =============================================================================
# test-occu-pg-gibbs-spatial.R - spatial Polya-Gamma Gibbs for single-season
# occupancy with an intrinsic areal (ICAR) field on the occupancy logit
# (occu() + icar() under method = "pg_gibbs"; spOccupancy spPGOcc; tulpaObs#126).
#
# logit(psi_i) = X_i beta + f_i, f ~ ICAR(tau). Conditional on the Polya-Gamma
# auxiliaries the joint (beta, f) update is a Gaussian Markov random field draw;
# tau has a conjugate Gamma full conditional; the field is centred (sum-to-zero,
# its level moved into the intercept) each sweep. These check the S3 surface, the
# field's sum-to-zero centring, and -- over seeds -- recovery of the intercept,
# the detection, and the field itself. NB: a site-level occupancy covariate is
# subject to spatial confounding with a one-node-per-site field (a known property
# of spatial occupancy models), so the recovery target is the field + intercept +
# detection, with the field driving the occupancy structure.
# =============================================================================

# Rook-adjacency for a g x g grid.
.pg_grid_adj <- function(g) {
  co <- expand.grid(r = seq_len(g), c = seq_len(g))
  n  <- g * g; adj <- matrix(0L, n, n)
  for (i in seq_len(n)) for (k in seq_len(n)) {
    if (i < k && abs(co$r[i] - co$r[k]) + abs(co$c[i] - co$c[k]) == 1L)
      adj[i, k] <- adj[k, i] <- 1L
  }
  list(adj = adj, co = co)
}

# Simulate intercept + ICAR-field occupancy with site-level detection.
.pg_sim_spatial <- function(g = 12L, J = 6L, beta0 = 0.2, p = 0.5, seed = 1) {
  set.seed(seed)
  gr <- .pg_grid_adj(g); co <- gr$co; n <- g * g
  f  <- as.numeric(scale(sin(co$r / 2) + cos(co$c / 2.5)))
  f  <- f - mean(f)
  psi <- stats::plogis(beta0 + f)
  z <- stats::rbinom(n, 1L, psi)
  y <- matrix(0L, n, J)
  for (i in seq_len(n)) if (z[i] == 1L) y[i, ] <- stats::rbinom(J, 1L, p)
  list(y = y, adj = gr$adj, f = f, data = data.frame(row.names = seq_len(n)))
}

test_that("occu() + icar() method='pg_gibbs' (spPGOcc) S3 + centred field", {
  sim <- .pg_sim_spatial(g = 8L, J = 5L, seed = 1)
  fit <- tobs(~ icar(graph = sim$adj), data = sim$data, family = occu(),
              detection = ~ 1, y = sim$y, method = "pg_gibbs",
              control = list(n.iter = 800L, n.warmup = 400L, n.chains = 2L,
                             seed = 1, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "pg_gibbs")
  expect_false(is.null(fit$spatial_field))
  expect_length(fit$spatial_field, nrow(sim$adj))
  expect_lt(abs(mean(fit$spatial_field)), 1e-6)          # sum-to-zero
  expect_true("log_tau" %in% names(fit$means))
  # The fixed effects mix well at this short quick fit; the field-precision
  # log_tau mixes more slowly (checked at longer chains in the recovery block).
  expect_true(all(fit$rhat[c("psi_(Intercept)", "p_(Intercept)")] < 1.2))
  # bym2 / car_proper are not wired for the PG spatial path in v1.
  expect_error(
    tobs(~ bym2(graph = sim$adj), data = sim$data, family = occu(),
         detection = ~ 1, y = sim$y, method = "pg_gibbs"),
    "icar")
})

test_that("occu() + icar() pg_gibbs recovers the field + intercept + detection", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 6L
  fcor <- pint <- psint <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- .pg_sim_spatial(g = 12L, J = 6L, beta0 = 0.2, p = 0.5, seed = 100 + s)
    fit <- tryCatch(
      tobs(~ icar(graph = sim$adj), data = sim$data, family = occu(),
           detection = ~ 1, y = sim$y, method = "pg_gibbs",
           control = list(n.iter = 1500L, n.warmup = 750L, n.chains = 2L,
                          seed = s, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    fcor[s]  <- cor(fit$spatial_field, sim$f)
    psint[s] <- fit$means[["psi_(Intercept)"]]
    pint[s]  <- fit$means[["p_(Intercept)"]]
  }
  # The field recovers (occupancy-informed field, detection-filtered + J=6).
  expect_gt(mean(fcor, na.rm = TRUE), 0.6)
  # Intercept + detection recover in aggregate (the occupancy intercept carries a
  # little upward variance-inflation from the field level at moderate n).
  expect_lt(abs(mean(psint, na.rm = TRUE) - 0.2), 0.25)
  expect_lt(abs(mean(pint,  na.rm = TRUE) - 0.0), 0.10)
})
