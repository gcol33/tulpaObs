# Correlated (`|`) / separable-MCAR spatial coefficient fields on the cover
# hurdle (gcol33/tulpaObs#64). A single bar `~ 1 + time | cell` declares the
# intercept field u_cell and the time-slope field s_cell as CORRELATED Besag
# fields with a free 2x2 cross-covariance Sigma (the within-arm covariance among
# the fields), then copies the whole correlated field onto the positive arm with
# one estimated amplitude alpha (the cross-arm transfer). The free Sigma is
# integrated over the outer mode-centred CCD in log-Cholesky coordinates; the
# recovered Sigma SDs and cross-correlation are weighted-quantile summaries of
# the marginalized posterior, never a plug-in of the modal cell.
#
# These tests are parameter-recovery against simulated truth (a KNOWN Sigma:
# both field SDs AND the cross-correlation rho, plus the cross-arm copy alpha),
# the recovery-test bar this statistical fitter must clear -- not a smoke test.

# Rook-adjacency on an nx x ny grid (self-contained so the file runs in isolation).
.mcar_grid_adj <- function(nx, ny) {
  n <- nx * ny
  W <- matrix(0L, n, n)
  id <- function(i, j) (j - 1L) * nx + i
  for (i in seq_len(nx)) for (j in seq_len(ny)) {
    if (i < nx) { a <- id(i, j); b <- id(i + 1L, j); W[a, b] <- W[b, a] <- 1L }
    if (j < ny) { a <- id(i, j); b <- id(i, j + 1L); W[a, b] <- W[b, a] <- 1L }
  }
  W
}

# Draw correlated fields (u, s) with cov Sigma (x) Q^-1 via the LMC factorization
# u = L11 z1, s = L21 z1 + L22 z2 (z1, z2 iid near-intrinsic CAR), L = chol(Sigma).
.mcar_sim_fields <- function(adj, Sigma) {
  n <- nrow(adj)
  Qp <- diag(rowSums(adj)) - 0.99 * adj
  U  <- chol(Qp)
  z1 <- backsolve(U, stats::rnorm(n)); z1 <- z1 - mean(z1)
  z2 <- backsolve(U, stats::rnorm(n)); z2 <- z2 - mean(z2)
  L  <- t(chol(Sigma))
  u <- L[1, 1] * z1
  s <- L[2, 1] * z1 + L[2, 2] * z2
  list(u = u - mean(u), s = s - mean(s))
}

# One correlated-MCAR cover-hurdle dataset (lognormal positive arm). High
# occupancy + lognormal cover so both fields are well identified.
.mcar_sim_cover <- function(adj, Sigma, alpha, seed,
                            n_per = 25L, psi0 = 1.2, pos0 = -0.3, sd_pos = 0.4) {
  set.seed(seed)
  n_s <- nrow(adj)
  fld <- .mcar_sim_fields(adj, Sigma)
  cell <- rep(seq_len(n_s), each = n_per)
  N    <- length(cell)
  time.sc <- stats::rnorm(N)
  eta_psi <- psi0 + fld$u[cell] + time.sc * fld$s[cell]
  z <- stats::rbinom(N, 1, stats::plogis(eta_psi))
  eta_pos <- pos0 + alpha * (fld$u[cell] + time.sc * fld$s[cell])
  cover <- numeric(N)
  pos <- z == 1L
  cover[pos] <- pmin(exp(stats::rnorm(sum(pos), eta_pos[pos], sd_pos)), 0.999)
  list(data = data.frame(cover = cover, time.sc = time.sc, cell = cell),
       y = cover, fld = fld)
}

test_that("correlated `|` bar recovers Sigma (SDs + cross-correlation) and alpha", {
  skip_on_cran()
  skip_if_fast()

  nx <- ny <- 10L
  adj <- .mcar_grid_adj(nx, ny)
  sig_u <- 1.0; sig_s <- 0.7; rho_true <- 0.5; alpha_true <- 1.3
  Sigma <- matrix(c(sig_u^2,             rho_true * sig_u * sig_s,
                    rho_true * sig_u * sig_s, sig_s^2), 2, 2)

  n_seeds <- 4L
  sd_u <- sd_s <- rho_hat <- alpha_hat <- numeric(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- .mcar_sim_cover(adj, Sigma, alpha_true, seed = r)
    fit <- suppressWarnings(tobs(
      formula = ~ time.sc + spatial(~ 1 + time.sc | cell, graph = adj,
                                    to = c("presence", "positive")),
      data = sim$data, y = sim$y, family = cover(positive = "lognormal"),
      method = "nested_laplace",
      control = list(max.iter = 60L, progress = FALSE, verbose = FALSE)))

    expect_s3_class(fit, "cover_fit")
    expect_true(isTRUE(fit$mcar))
    expect_length(fit$sigma_mcar, 2L)
    expect_length(fit$rho_mcar, 1L)
    expect_true(is.finite(fit$alpha_mcar) && fit$alpha_mcar > 0)
    # The integrated outer grid is the mode-centred CCD over the log-Cholesky
    # coordinates of Sigma + the alpha copy axis.
    expect_identical(fit$joint$integration, "ccd")

    sd_u[r]      <- fit$sigma_mcar[[1L]]
    sd_s[r]      <- fit$sigma_mcar[[2L]]
    rho_hat[r]   <- fit$rho_mcar[[1L]]
    alpha_hat[r] <- fit$alpha_mcar
  }

  # Field SDs recover with low bias every seed (the diagonal of Sigma is the
  # best-identified part of the cross-covariance).
  expect_lt(abs(mean(sd_u) - sig_u) / sig_u, 0.20)
  expect_lt(abs(mean(sd_s) - sig_s) / sig_s, 0.25)
  # The free cross-correlation -- the quantity `|` adds over `||` -- recovers a
  # clearly positive value (true 0.5), and the per-seed estimates are positive
  # in the majority of seeds (not collapsed to 0 or pinned to +-1).
  expect_gt(mean(rho_hat), 0.30)
  expect_gte(sum(rho_hat > 0.15), 3L)
  expect_true(all(rho_hat > -0.9 & rho_hat < 0.99))
  # The cross-arm copy amplitude is positive and within a factor of ~2 of truth
  # (alpha = sigma_pos / sigma_presence is weakly identified at small n_pos).
  expect_gt(mean(alpha_hat), 0.5)
  expect_lt(mean(alpha_hat), 2.6)
})

test_that("independent `||` bar does not recover a spurious cross-correlation", {
  skip_on_cran()
  skip_if_fast()

  nx <- ny <- 10L
  adj <- .mcar_grid_adj(nx, ny)
  # Data simulated WITH a correlated truth, but fit with the INDEPENDENT `||`
  # spelling: the independent path has separate per-field precisions and NO
  # cross-covariance parameter at all, so it structurally cannot manufacture a
  # cross-correlation. This is the `||`-vs-`|` distinction (#61 vs #64).
  sig_u <- 1.0; sig_s <- 0.7; rho_true <- 0.5; alpha_true <- 1.3
  Sigma <- matrix(c(sig_u^2,             rho_true * sig_u * sig_s,
                    rho_true * sig_u * sig_s, sig_s^2), 2, 2)
  sim <- .mcar_sim_cover(adj, Sigma, alpha_true, seed = 1L)

  fit <- suppressWarnings(tobs(
    formula = ~ time.sc + spatial(~ 1 + time.sc || cell, graph = adj,
                                  to = c("presence", "positive")),
    data = sim$data, y = sim$y, family = cover(positive = "lognormal"),
    method = "nested_laplace",
    control = list(max.iter = 60L, progress = FALSE, verbose = FALSE)))

  expect_s3_class(fit, "cover_fit")
  # The independent fit is NOT an MCAR fit and carries no cross-correlation.
  expect_false(isTRUE(fit$mcar))
  expect_null(fit$rho_mcar)
  expect_null(fit$sigma_mcar)
  # It does fit the two independent fields (intercept + trend) with their own
  # amplitudes -- the `||` deliverable.
  expect_false(is.null(fit$sigma_trend) || is.na(fit$sigma_trend))
})


# ---- single-arm correlated `|` field on the occurrence arm alone (#109) ------
#
# A correlated intercept + slope field on `to = "presence"` only: a free 2x2
# cross-coefficient Sigma on the occupancy arm with NO cross-arm copy (no alpha).
# The occupancy intercept and time-slope fields covary; the positive arm carries
# none. This is the single-arm counterpart of the both-arm copy path above.

# One dataset with the correlated field on the OCCUPANCY arm only (the positive
# arm's cover has no spatial structure).
.mcar_sim_cover_occ <- function(adj, Sigma, seed,
                                n_per = 25L, psi0 = 1.2, pos0 = -0.3, sd_pos = 0.4) {
  set.seed(seed)
  n_s <- nrow(adj)
  fld <- .mcar_sim_fields(adj, Sigma)
  cell <- rep(seq_len(n_s), each = n_per)
  N    <- length(cell)
  time.sc <- stats::rnorm(N)
  eta_psi <- psi0 + fld$u[cell] + time.sc * fld$s[cell]
  z <- stats::rbinom(N, 1, stats::plogis(eta_psi))
  cover <- numeric(N); pos <- z == 1L
  cover[pos] <- pmin(exp(stats::rnorm(sum(pos), pos0, sd_pos)), 0.999)  # no field
  list(data = data.frame(cover = cover, time.sc = time.sc, cell = cell), y = cover)
}

test_that("single-arm `|` (to = 'presence') recovers Sigma, no cross-arm copy", {
  skip_on_cran()
  skip_if_fast()

  nx <- ny <- 10L
  adj <- .mcar_grid_adj(nx, ny)
  sig_u <- 1.0; sig_s <- 0.7; rho_true <- 0.5
  Sigma <- matrix(c(sig_u^2,                 rho_true * sig_u * sig_s,
                    rho_true * sig_u * sig_s, sig_s^2), 2, 2)

  n_seeds <- 4L
  sd_u <- sd_s <- rho_hat <- numeric(n_seeds)
  for (r in seq_len(n_seeds)) {
    sim <- .mcar_sim_cover_occ(adj, Sigma, seed = r)
    fit <- suppressWarnings(tobs(
      formula = ~ time.sc + spatial(~ 1 + time.sc | cell, graph = adj,
                                    to = "presence"),
      data = sim$data, y = sim$y, family = cover(positive = "lognormal"),
      method = "nested_laplace",
      control = list(max.iter = 60L, progress = FALSE, verbose = FALSE)))

    expect_s3_class(fit, "cover_fit")
    expect_true(isTRUE(fit$mcar))
    expect_length(fit$sigma_mcar, 2L)
    expect_length(fit$rho_mcar, 1L)
    # Single arm: no copy, so no alpha amplitude.
    expect_true(is.na(fit$alpha_mcar))

    sd_u[r]    <- fit$sigma_mcar[[1L]]
    sd_s[r]    <- fit$sigma_mcar[[2L]]
    rho_hat[r] <- fit$rho_mcar[[1L]]
  }

  # Field SDs recover with low bias (the diagonal of Sigma is best identified).
  expect_lt(abs(mean(sd_u) - sig_u) / sig_u, 0.20)
  expect_lt(abs(mean(sd_s) - sig_s) / sig_s, 0.25)
  # The free cross-correlation -- the quantity the single `|` adds over `||` --
  # recovers a clearly positive value (true 0.5), positive in the majority of
  # seeds and never pinned to the boundary.
  expect_gt(mean(rho_hat), 0.30)
  expect_gte(sum(rho_hat > 0.15), 3L)
  expect_true(all(rho_hat > -0.9 & rho_hat < 0.99))
})
