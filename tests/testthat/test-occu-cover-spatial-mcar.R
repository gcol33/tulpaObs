# Correlated (`|`) / separable-MCAR spatial coefficient fields on occu_cover()
# (gcol33/tulpaObs#63). A single bar `~ 1 + x | cell` on the occupancy (psi)
# formula declares the intercept field u_cell and the x-slope field s_cell as
# CORRELATED Besag fields with a free 2x2 cross-covariance Sigma, copied onto the
# cover arm with one amplitude alpha (the cross-arm transfer). The free Sigma is
# integrated over the outer mode-centred CCD in log-Cholesky coordinates; the
# reported Sigma SDs and cross-correlation are weighted-quantile summaries of the
# marginalized posterior, never a plug-in of the modal cell. Unlike the cover
# hurdle the field rides the coupled occupancy mixture, so it exercises the
# coupled per-cell scatter's INDEXED_MULTI path (tulpa engine).

.mcar_grid_adj <- function(nx, ny) {
  n <- nx * ny; W <- matrix(0L, n, n); id <- function(i, j) (j - 1L) * nx + i
  for (i in seq_len(nx)) for (j in seq_len(ny)) {
    if (i < nx) { a <- id(i, j); b <- id(i + 1L, j); W[a, b] <- W[b, a] <- 1L }
    if (j < ny) { a <- id(i, j); b <- id(i, j + 1L); W[a, b] <- W[b, a] <- 1L }
  }
  W
}

# Correlated (u, s) with cov Sigma (x) Q^-1 via the LMC factorization
# u = L11 z1, s = L21 z1 + L22 z2 (z1, z2 iid near-intrinsic CAR), L = chol(Sigma).
.mcar_fields <- function(adj, Sigma) {
  n <- nrow(adj); Qp <- diag(rowSums(adj)) - 0.99 * adj; U <- chol(Qp)
  z1 <- backsolve(U, stats::rnorm(n)); z1 <- z1 - mean(z1)
  z2 <- backsolve(U, stats::rnorm(n)); z2 <- z2 - mean(z2)
  L <- t(chol(Sigma)); u <- L[1, 1] * z1; s <- L[2, 1] * z1 + L[2, 2] * z2
  list(u = u - mean(u), s = s - mean(s))
}

# One correlated-MCAR occu_cover dataset (lognormal positive arm): a per-cell
# correlated (intercept, x-slope) field on occupancy, copied to cover with alpha.
.mcar_occu_data <- function(adj, Sigma, alpha, seed, J = 8L,
                            psi0 = 1.3, p_det = 0.7, pos0 = -0.3, sd_pos = 0.4) {
  set.seed(seed); n <- nrow(adj)
  fld <- .mcar_fields(adj, Sigma); x <- stats::rnorm(n)
  z <- stats::rbinom(n, 1, stats::plogis(psi0 + fld$u + x * fld$s))
  y <- matrix(0L, n, J)
  for (i in seq_len(n)) if (z[i] == 1L) y[i, ] <- stats::rbinom(J, 1, p_det)
  eta_pos <- pos0 + alpha * (fld$u + x * fld$s)
  y_pos <- matrix(0, n, J)
  for (i in seq_len(n)) for (j in seq_len(J))
    if (y[i, j] == 1L) y_pos[i, j] <- exp(stats::rnorm(1, eta_pos[i], sd_pos))
  long <- data.frame(site_id = rep(seq_len(n), each = J),
                     visit = rep(seq_len(J), times = n),
                     y = as.vector(t(y)),
                     det_cov1 = stats::rnorm(n * J),
                     pos_cov1 = stats::rnorm(n * J))
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  list(od = od, cell_dat = data.frame(site_id = seq_len(n), x = x),
       y_pos = y_pos, fld = fld, adj = adj)
}

.mcar_fit <- function(d, max.iter = 80L) {
  suppressWarnings(tobs(
    formula = ~ spatial(~ 1 + x | site_id, graph = d$adj),
    data = d$cell_dat, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1 + copy(spatial()),
    y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
    method = "nested_laplace",
    control = list(max.iter = max.iter, progress = FALSE, verbose = FALSE)))
}


# --- Structural / scope gates (fast, parse-level) ------------------------------

test_that("a correlated `|` occupancy bar routes to the MCAR field", {
  adj <- .mcar_grid_adj(4L, 4L)
  dat <- data.frame(site_id = seq_len(16), x = rnorm(16))
  si  <- tulpaObs:::.occu_cover_spatial_fields(
    ~ spatial(~ 1 + x | site_id, graph = adj), dat)
  expect_true(isTRUE(si$correlated))
  expect_length(si$fields, 2L)            # intercept + x-slope coefficient field
})

test_that("an independent `||` occupancy bar is NOT flagged correlated", {
  adj <- .mcar_grid_adj(4L, 4L)
  dat <- data.frame(site_id = seq_len(16), x = rnorm(16))
  si  <- tulpaObs:::.occu_cover_spatial_fields(
    ~ spatial(~ 1 + x || site_id, graph = adj), dat)
  expect_false(isTRUE(si$correlated))
  expect_length(si$fields, 2L)
})

test_that("a correlated `|` bar enforces its scope gates", {
  adj <- .mcar_grid_adj(4L, 4L)
  dat <- data.frame(site_id = seq_len(16), x = rnorm(16), g = rep(1:4, 4))
  # Single field (intercept only) has no cross-covariance.
  expect_error(
    tulpaObs:::.occu_cover_spatial_fields(
      ~ spatial(~ 1 | site_id, graph = adj), dat),
    "at least one coefficient")
  # Must be the only spatial term.
  expect_error(
    tulpaObs:::.occu_cover_spatial_fields(
      ~ spatial(~ 1 + x | site_id, graph = adj) + icar(graph = adj), dat),
    "only spatial term|MCAR field structure")
  # Does not compose with a per-group RE.
  expect_error(
    tulpaObs:::.occu_cover_spatial_fields(
      ~ spatial(~ 1 + x | site_id, graph = adj) + (1 | g), dat),
    "random effect")
})


# --- Parameter recovery (the recovery-test bar this fitter must clear) ---------

test_that("correlated `|` occupancy bar recovers Sigma (SDs + rho) and alpha", {
  skip_on_cran()
  skip_if_fast()

  adj <- .mcar_grid_adj(10L, 10L)
  sig_u <- 1.0; sig_s <- 0.7; rho_true <- 0.5; alpha_true <- 1.3
  Sigma <- matrix(c(sig_u^2,                  rho_true * sig_u * sig_s,
                    rho_true * sig_u * sig_s, sig_s^2), 2, 2)

  n_seeds <- 4L
  sd_u <- sd_s <- rho_hat <- alpha_hat <- numeric(n_seeds)
  field_cor <- numeric(n_seeds)
  for (r in seq_len(n_seeds)) {
    d   <- .mcar_occu_data(adj, Sigma, alpha_true, seed = r)
    fit <- .mcar_fit(d)

    expect_s3_class(fit, "tobs_fit")
    expect_identical(fit$spatial$type, "mcar")
    expect_length(fit$spatial$sigma_mcar, 2L)
    expect_length(fit$spatial$rho_mcar, 1L)
    expect_true(is.finite(fit$spatial$alpha_mcar) && fit$spatial$alpha_mcar > 0)
    # The intercept field is recovered (the field DOES couple into the occupancy
    # mixture -- the INDEXED_MULTI coupled-scatter path, tulpa engine).
    expect_gt(stats::sd(fit$spatial_field), 1e-6)

    sd_u[r]      <- fit$spatial$sigma_mcar[[1L]]
    sd_s[r]      <- fit$spatial$sigma_mcar[[2L]]
    rho_hat[r]   <- fit$spatial$rho_mcar[[1L]]
    alpha_hat[r] <- fit$spatial$alpha_mcar
    field_cor[r] <- stats::cor(fit$spatial_field, d$fld$u)
  }

  # The intercept field tracks the simulated truth every seed.
  expect_gt(mean(field_cor), 0.6)
  # Field SDs recover with modest bias (the occupancy field is harder to identify
  # than the cover hurdle's, and the Laplace under-disperses small-group SDs).
  expect_lt(abs(mean(sd_u) - sig_u) / sig_u, 0.30)
  expect_lt(abs(mean(sd_s) - sig_s) / sig_s, 0.40)
  # The free cross-correlation -- the quantity `|` adds over `||` -- recovers a
  # clearly positive value (true 0.5) and is not pinned to +-1.
  expect_gt(mean(rho_hat), 0.30)
  expect_true(all(rho_hat > -0.9 & rho_hat < 0.99))
  # The cross-arm copy amplitude is recovered within a factor of ~2 of truth.
  expect_gt(mean(alpha_hat), 0.6)
  expect_lt(mean(alpha_hat), 2.6)
})

test_that("independent `||` occupancy bar carries no cross-correlation", {
  skip_on_cran()
  skip_if_fast()

  adj <- .mcar_grid_adj(10L, 10L)
  sig_u <- 1.0; sig_s <- 0.7; rho_true <- 0.5; alpha_true <- 1.3
  Sigma <- matrix(c(sig_u^2,                  rho_true * sig_u * sig_s,
                    rho_true * sig_u * sig_s, sig_s^2), 2, 2)
  d <- .mcar_occu_data(adj, Sigma, alpha_true, seed = 1L)

  fit <- suppressWarnings(tobs(
    formula = ~ spatial(~ 1 + x || site_id, graph = d$adj),
    data = d$cell_dat, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
    method = "nested_laplace",
    control = list(max.iter = 80L, progress = FALSE, verbose = FALSE)))

  expect_s3_class(fit, "tobs_fit")
  # The independent fit is NOT an MCAR fit and carries no cross-correlation.
  expect_identical(fit$spatial$type, "icar")
  expect_null(fit$spatial$rho_mcar)
  # It does fit the two independent fields (intercept + trend) with their own
  # amplitudes -- the `||` deliverable (gcol33/tulpaObs#61).
  expect_false(is.null(fit$spatial$sigma_trend_mean))
})
