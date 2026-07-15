# Parameter recovery for an SPDE mesh field on the STATE arm of the integrated,
# JSDM, community, and dynamic occupancy paths (gcol33/tulpaObs#21). Each fits a
# model whose state (psi / season-1 psi1) linear predictor carries a smooth
# Matern-like signal `u_true`, through the SPDE-Laplace EM, and asserts
#   (a) state fixed-effect (slope) recovery within tolerance, and
#   (b) field-shape recovery cor(A u_hat, u_true) above a per-arm bound.
#
# Field-shape bounds differ by arm because the SPDE field is recovered through
# each arm's latent state:
#   * jsdm / community pool several independent presence observations per site
#     (one per species, sharing a single site-level field), so the field is
#     sharply identified -- cor comfortably > 0.7.
#   * integrated / dynamic carry ONE latent occupancy state per site (sources
#     share it; the dynamic field enters season-1 psi1 only), so the field is
#     recovered exactly as well as the proven single-season state arm and caps
#     around 0.6-0.7 even at large N -- the same ceiling test-spde-occ.R sees
#     for the single-season state field (its bound is 0.35-0.40). The bounds
#     below (0.5 per seed, 0.6 aggregate) sit above noise and catch structural
#     regressions (field not wired, sign flip, the broadcast dropped) without
#     flagging a correct fit.

`%||%` <- function(a, b) if (is.null(a)) b else a

.spde_field_term <- function() {
  ~ occ_cov + spde(lon, lat, max_edge = c(0.3, 0.6), nu = 1,
                   prior_range = c(0.3, 0.5), prior_sigma = c(0.8, 0.5))
}
.spde_field_cor <- function(fit, u_true) {
  cor(as.numeric(fit$spatial$tulpa_spec$A %*% fit$spatial_field), u_true)
}
.spde_sfield <- function(coords, amp = 0.9) {
  amp * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
}

# --- JSDM: shared site-level field across species, no detection -------------
# Calibration over 6 seeds (dev_notes/probe_spde_calib.R, N = 600, 6 species,
# field amp 0.9): x-slope 0.535-0.629 (truth 0.6), field_cor min 0.765, mean
# 0.806. The field is shared across species.
#
# Since gcol33/tulpaObs#121 jsdm() is the COMMUNITY GLMM (per-species
# coefficients under a Gaussian community covariance), so its spde() field is
# fit by the shared block-coordinate driver (R/community_latent.R) rather than
# the single-block SPDE-Laplace EM: mesh nodes carry the field, the barycentric
# projector A maps them onto sites, the precision is
# Q(kappa) = kappa^4 C0 + 2 kappa^2 G1 + G1 C0^-1 G1 scaled by the driver's tau,
# and kappa (the Matern range) is gridded by the field marginal. That makes it
# method = "nested_laplace", alongside the areal community fields, and the
# community-mean coefficient is `mu_*`.
.sim_spde_jsdm <- function(seed, n_sites = 600, n_species = 6) {
  set.seed(seed)
  coords <- cbind(runif(n_sites), runif(n_sites))
  u_true <- .spde_sfield(coords)
  x_cov  <- rnorm(n_sites)
  beta0  <- c(-0.3, 0.1, 0.4, -0.1, 0.2, -0.2)[seq_len(n_species)]
  y <- matrix(0L, n_sites, n_species)
  for (k in seq_len(n_species)) {
    y[, k] <- rbinom(n_sites, 1, plogis(beta0[k] + 0.6 * x_cov + u_true))
  }
  colnames(y) <- paste0("sp", seq_len(n_species))
  list(data = data.frame(occ_cov = x_cov, lon = coords[, 1], lat = coords[, 2]),
       y = y, u_true = u_true)
}

test_that("tobs() + jsdm() state-arm spde() Laplace recovers beta and the field", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")

  sim <- .sim_spde_jsdm(seed = 1)
  fit <- tobs(.spde_field_term(), data = sim$data, family = jsdm(),
              y = sim$y, species = TRUE, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))

  expect_lt(abs(fit$means["mu_occ_cov"] - 0.6), 0.20)
  expect_false(is.null(fit$spatial_field))
  expect_equal(length(fit$spatial_field), fit$spatial$n_units)
  expect_gt(.spde_field_cor(fit, sim$u_true), 0.7)
  # the Matern range is estimated on the kappa grid
  expect_true(is.finite(fit$spatial_hyper$range))
})

test_that("jsdm() state-arm spde() recovery holds across seeds", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")

  seeds <- c(2L, 3L, 4L)
  cors <- vapply(seeds, function(s) {
    sim <- .sim_spde_jsdm(seed = s)
    fit <- tobs(.spde_field_term(), data = sim$data, family = jsdm(),
                y = sim$y, species = TRUE, method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    .spde_field_cor(fit, sim$u_true)
  }, numeric(1))
  for (c_k in cors) expect_gt(c_k, 0.7)
  expect_gt(mean(cors), 0.72)
})

# --- INTEGRATED: shared psi field, two sources ------------------------------
# Calibration over 6 seeds (N = 1200, 2 sources, field amp 1.0): psi-slope
# 0.649-0.816 (truth 0.7), field_cor min 0.613, mean 0.667 -- the single-season
# state-arm ceiling (the field is on the shared latent occupancy state).
.sim_spde_integrated <- function(seed, n_sites = 1200) {
  set.seed(seed)
  coords <- cbind(runif(n_sites), runif(n_sites))
  u_true <- .spde_sfield(coords, amp = 1.0)
  x_cov  <- rnorm(n_sites)
  z <- rbinom(n_sites, 1, plogis(0.3 + 0.7 * x_cov + u_true))
  mk_src <- function(p0, J) {
    y <- matrix(0L, n_sites, J)
    for (i in seq_len(n_sites)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p0)
    y
  }
  list(data = data.frame(occ_cov = x_cov, lon = coords[, 1], lat = coords[, 2]),
       y = list(src1 = mk_src(0.65, 10), src2 = mk_src(0.55, 8)),
       u_true = u_true)
}

test_that("tobs() + int_occu() state-arm spde() Laplace recovers beta and the field", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")

  sim <- .sim_spde_integrated(seed = 1)
  fit <- tobs(.spde_field_term(), data = sim$data, family = int_occu(),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))

  expect_lt(abs(fit$means["psi_occ_cov"] - 0.7), 0.25)
  expect_false(is.null(fit$spatial_field))
  expect_equal(length(fit$spatial_field), fit$spatial$n_units)
  expect_gt(.spde_field_cor(fit, sim$u_true), 0.5)
})

test_that("int_occu() state-arm spde() recovery holds across seeds", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")

  seeds <- c(2L, 3L, 4L)
  cors <- vapply(seeds, function(s) {
    sim <- .sim_spde_integrated(seed = s)
    fit <- tobs(.spde_field_term(), data = sim$data, family = int_occu(),
                detection = ~ 1, y = sim$y, method = "laplace",
                control = list(verbose = FALSE))
    .spde_field_cor(fit, sim$u_true)
  }, numeric(1))
  for (c_k in cors) expect_gt(c_k, 0.5)
  expect_gt(mean(cors), 0.6)
})

# --- DYNAMIC: psi1 (season-1) field -----------------------------------------
# Calibration over 6 seeds (N = 1200, T = 5, J = 8, field amp 1.0): psi1-slope
# 0.649-0.816 (truth 0.7), field_cor min 0.620, mean 0.667 -- identical to the
# integrated / single-season state-arm ceiling (the field enters psi1 only).
.sim_spde_dynamic <- function(seed, n_sites = 1200, n_seasons = 5, J = 8) {
  set.seed(seed)
  coords <- cbind(runif(n_sites), runif(n_sites))
  u_true <- .spde_sfield(coords, amp = 1.0)
  x_cov  <- rnorm(n_sites)
  gamma <- 0.25; eps <- 0.15; p0 <- 0.6
  z <- matrix(0L, n_sites, n_seasons)
  z[, 1] <- rbinom(n_sites, 1, plogis(0.3 + 0.7 * x_cov + u_true))
  for (t in 2:n_seasons) {
    for (i in seq_len(n_sites)) {
      z[i, t] <- if (z[i, t - 1] == 1L) rbinom(1, 1, 1 - eps)
                 else rbinom(1, 1, gamma)
    }
  }
  y <- array(0L, dim = c(n_sites, J, n_seasons))
  for (t in seq_len(n_seasons)) {
    for (i in seq_len(n_sites)) {
      if (z[i, t] == 1L) y[i, , t] <- rbinom(J, 1, p0)
    }
  }
  list(data = data.frame(occ_cov = x_cov, lon = coords[, 1], lat = coords[, 2]),
       y = y, u_true = u_true)
}

test_that("tobs() + dyn_occu() psi1 spde() Laplace recovers beta and the field", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")

  sim <- .sim_spde_dynamic(seed = 1)
  fit <- tobs(.spde_field_term(), data = sim$data, family = dyn_occu(),
              detection = ~ 1, y = sim$y, colonization = ~ 1, extinction = ~ 1,
              method = "laplace", control = list(verbose = FALSE))

  expect_lt(abs(fit$means["psi1_occ_cov"] - 0.7), 0.25)
  expect_false(is.null(fit$spatial_field))
  expect_equal(length(fit$spatial_field), fit$spatial$n_units)
  expect_gt(.spde_field_cor(fit, sim$u_true), 0.5)
})

test_that("dyn_occu() psi1 spde() recovery holds across seeds", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("tulpaMesh")

  seeds <- c(2L, 3L, 4L)
  cors <- vapply(seeds, function(s) {
    sim <- .sim_spde_dynamic(seed = s)
    fit <- tobs(.spde_field_term(), data = sim$data, family = dyn_occu(),
                detection = ~ 1, y = sim$y, colonization = ~ 1, extinction = ~ 1,
                method = "laplace", control = list(verbose = FALSE))
    .spde_field_cor(fit, sim$u_true)
  }, numeric(1))
  for (c_k in cors) expect_gt(c_k, 0.5)
  expect_gt(mean(cors), 0.6)
})
