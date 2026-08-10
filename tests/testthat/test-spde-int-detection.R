# Detection-arm SPDE mesh field on an INTEGRATED (multi-source) occupancy model
# (gcol33/tulpaObs#216). A term's process is the formula it is written in, so an
# spde() on the `detection` formula loads a continuous Matern field on the
# detection logit; the shared latent occupancy state stays whatever the
# occupancy formula says. The detection arm is one arm observed through S
# sources, so the field is fit once per source block -- each source recovers its
# own realization at its own sites -- and `fit$spatial_field_det` is the
# per-source named list of mesh fields.
#
# Simulator: psi non-spatial (occ_cov slope 0.7); both sources' detection logits
# carry the SAME smooth signal `u_true` plus a det_cov slope (0.4) on top of
# per-source intercepts. Two sources at every site, 8 and 6 visits.

`%||%` <- function(a, b) if (is.null(a)) b else a

.spde_int_det_term <- function() {
  ~ det_cov + spde(lon, lat, max_edge = c(0.25, 0.5), nu = 1,
                   prior_range = c(0.3, 0.5), prior_sigma = c(0.8, 0.5))
}

.sim_spde_int_det <- function(seed, n_sites = 600, J1 = 8, J2 = 6,
                              beta_occ = c(0.2, 0.7), beta_det = 0.4) {
  set.seed(seed)
  coords  <- cbind(runif(n_sites), runif(n_sites))
  u_true  <- 0.9 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
  x_cov   <- rnorm(n_sites)
  det_cov <- rnorm(n_sites)
  z <- rbinom(n_sites, 1, plogis(beta_occ[1] + beta_occ[2] * x_cov))
  mk <- function(J, p0) {
    p <- plogis(p0 + beta_det * det_cov + u_true)
    y <- matrix(0L, n_sites, J)
    for (i in seq_len(n_sites)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
    y
  }
  list(data = data.frame(occ_cov = x_cov, det_cov = det_cov,
                         lon = coords[, 1], lat = coords[, 2]),
       y = list(src1 = mk(J1, -0.2), src2 = mk(J2, -0.5)),
       u_true = u_true)
}

.fit_spde_int_det <- function(sim, detection = .spde_int_det_term(),
                              method = "laplace") {
  tobs(~ occ_cov, data = sim$data, family = int_occu(),
       detection = detection, y = sim$y, method = method,
       control = list(verbose = FALSE, progress = FALSE))
}

.int_det_field_cor <- function(fit, u_true) {
  A <- fit$spatial$tulpa_spec$A
  vapply(fit$spatial_field_det,
         function(u) cor(as.numeric(A %*% u), u_true), numeric(1))
}

# --- the term reaches the fitter branch -------------------------------------
# Structural, no fit: the builder tags the detection formula's term with the
# detection process, so `.tobs_structures_from_model()` hands the fitter a spec
# with shared = c(FALSE, TRUE) and `build_integrated_callbacks()` resolves a
# per-source projection. This is the wiring the fitter branch was written for
# and that no formula used to reach.
test_that("an integrated detection formula carries its structured term to the arm", {
  skip_if_no_tulpamesh()

  sim   <- .sim_spde_int_det(seed = 1, n_sites = 40, J1 = 3, J2 = 2)
  model <- .tobs_build_model(occ_formula = ~ occ_cov,
                             det_formula = .spde_int_det_term(),
                             data = sim$data, y = sim$y, integrated = TRUE)

  structs <- .tobs_structures_from_model(model)
  expect_false(is.null(structs$spatial))
  expect_identical(structs$spatial$type, "spde")
  expect_identical(structs$spatial$shared, c(FALSE, TRUE))
  expect_null(.spatial_for_arm(structs$spatial, 1L))
  expect_false(is.null(.spatial_for_arm(structs$spatial, 2L)))

  # The per-source detection blocks carry the field, broadcast onto each
  # source's own sites; the occupancy block does not.
  cb <- build_integrated_callbacks(model, structs$spatial)
  specs <- cb$m_step_encode(rep(0.5, model$n_sites))
  expect_null(specs$occ$spatial)
  for (s in seq_len(model$n_sources)) {
    blk <- specs[[paste0("det", s)]]
    expect_false(is.null(blk$spatial))
    expect_equal(nrow(blk$spatial$A), nrow(sim$y[[s]]))
  }
})

# --- recovery ---------------------------------------------------------------
# Calibration at N = 600 (seeds 1-3, field amp 0.9): per-source field cor
# 0.872-0.926, occupancy slope 0.613-0.692 (truth 0.7), detection slope
# 0.335-0.455 (truth 0.4). The bounds below sit outside that spread and catch
# structural regressions (the field not wired, a sign flip, the per-source
# broadcast dropped).
test_that("integrated detection-arm spde() recovers the per-source field", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  sim <- .sim_spde_int_det(seed = 1)
  fit <- .fit_spde_int_det(sim)

  # The field lives on the detection arm, one realization per source; the
  # occupancy arm carries none.
  expect_true(is.list(fit$spatial_field_det))
  expect_identical(names(fit$spatial_field_det), c("src1", "src2"))
  expect_null(fit$spatial_field)
  for (u in fit$spatial_field_det) expect_equal(length(u), fit$spatial$n_units)

  cors <- .int_det_field_cor(fit, sim$u_true)
  for (k in seq_along(cors)) expect_gt(cors[k], 0.75)

  expect_lt(abs(fit$means["psi_occ_cov"] - 0.7), 0.25)
  expect_lt(abs(fit$means["src1_det_cov"] - 0.4), 0.20)
  expect_lt(abs(fit$means["src2_det_cov"] - 0.4), 0.20)
  expect_true(is.finite(fit$sds["src1_det_cov"]))
})

test_that("integrated detection-arm spde() recovery holds across seeds", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  seeds <- c(2L, 3L)
  cors <- unlist(lapply(seeds, function(s) {
    sim <- .sim_spde_int_det(seed = s)
    .int_det_field_cor(.fit_spde_int_det(sim), sim$u_true)
  }))
  for (c_k in cors) expect_gt(c_k, 0.75)
  expect_gt(mean(cors), 0.80)
})

# The field is doing work, not riding along: dropping the term from the same
# data moves the fit. Dropping it leaves the per-site detection signal in the
# residual, and both sources' detection slopes come back attenuated below truth
# (seed 1: 0.24 / 0.37 against 0.40); the field lifts both. It lifts src2 past
# truth (0.46), so the assertion is on the direction and size of the shift, not
# on the fit with the field being nearer truth on every coefficient.
test_that("the detection field changes the fit it is dropped from", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  sim      <- .sim_spde_int_det(seed = 1)
  fit      <- .fit_spde_int_det(sim)
  fit_flat <- .fit_spde_int_det(sim, detection = ~ det_cov)

  expect_null(fit_flat$spatial_field_det)
  slopes <- c("src1_det_cov", "src2_det_cov")
  for (nm in slopes) {
    expect_lt(fit_flat$means[[nm]], 0.4)
    expect_gt(fit$means[[nm]] - fit_flat$means[[nm]], 0.05)
  }
  # The per-source detection intercepts move too -- the field is not absorbed
  # into a single arm-wide level.
  ints <- c("src1_(Intercept)", "src2_(Intercept)")
  expect_gt(min(abs(fit$means[ints] - fit_flat$means[ints])), 0.02)
})

# --- gates ------------------------------------------------------------------
# Everything the detection arm does NOT fit errors with a pointer. It is fit by
# the single-Laplace EM, which attaches a continuous mesh field per M-step
# block; the areal kinds and the temporal / re / svc / latent classes are built
# into the multi-block latent prior the nested-Laplace path attaches to the
# STATE block, so on a detection formula they would be fit against the wrong
# arm.
.int_gate_data <- function(n = 60) {
  set.seed(7)
  d <- data.frame(occ_cov = rnorm(n), det_cov = rnorm(n),
                  lon = runif(n), lat = runif(n),
                  grp = factor(rep(seq_len(6), length.out = n)),
                  yr  = rep(seq_len(5), length.out = n))
  z <- rbinom(n, 1, 0.6)
  mk <- function(J, p) {
    y <- matrix(0L, n, J)
    for (i in seq_len(n)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p)
    y
  }
  list(data = d, y = list(src1 = mk(5, 0.6), src2 = mk(4, 0.5)))
}

.int_gate_fit <- function(g, detection, occ = ~ occ_cov, method = "laplace") {
  tobs(occ, data = g$data, family = int_occu(), detection = detection,
       y = g$y, method = method,
       control = list(verbose = FALSE, progress = FALSE))
}

test_that("non-spde structured terms on an integrated detection formula error", {
  g <- .int_gate_data()
  expect_error(.int_gate_fit(g, ~ det_cov + (1 | grp)),
               "detection formula is not supported")
  expect_error(.int_gate_fit(g, ~ det_cov + temporal(yr, type = "ar1")),
               "detection formula is not supported")
  expect_error(.int_gate_fit(g, ~ det_cov + latent(2)),
               "detection formula is not supported")
})

test_that("an areal field on an integrated detection formula errors with a pointer", {
  g <- .int_gate_data()
  adj <- matrix(0L, nrow(g$data), nrow(g$data))
  for (i in seq_len(nrow(g$data) - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  expect_error(.int_gate_fit(g, ~ det_cov + icar(graph = adj)),
               "areal kinds are grid-integrated on the state arm")
})

test_that("a field on both integrated arms at once errors", {
  skip_if_no_tulpamesh()
  g <- .int_gate_data()
  expect_error(
    .int_gate_fit(g, .spde_int_det_term(),
                  occ = ~ occ_cov + spde(lon, lat, max_edge = c(0.3, 0.6),
                                         nu = 1, prior_range = c(0.3, 0.5),
                                         prior_sigma = c(0.8, 0.5))),
    "carries one spatial field")
})

# The detection arm is one arm read by S sources, so it is one detection
# formula. Per-source formulas that differ carry fixed effects only.
test_that("a structured term on differing per-source formulas errors", {
  skip_if_no_tulpamesh()
  g <- .int_gate_data()
  expect_error(
    .int_gate_fit(g, list(.spde_int_det_term(), ~ det_cov)),
    "same detection formula at every source")
  # Per-source formulas written the same way ARE the shared arm: the same
  # detection model at every source, spelled out once per source.
  m <- .tobs_build_model(occ_formula = ~ occ_cov,
                         det_formula = list(.spde_int_det_term(),
                                            .spde_int_det_term()),
                         data = g$data, y = g$y, integrated = TRUE)
  expect_identical(.tobs_structures_from_model(m)$spatial$shared,
                   c(FALSE, TRUE))
})

# The nested-Laplace route builds its latent blocks upstream and attaches the
# whole multi-block prior to the STATE block, with no arm channel -- a
# detection-arm term there would be fit against occupancy. It stops instead.
test_that("nested_laplace rejects a detection-arm term rather than moving it", {
  skip_if_no_tulpamesh()
  g <- .int_gate_data()
  expect_error(.int_gate_fit(g, .spde_int_det_term(), method = "nested_laplace"),
               "attaches its latent blocks to the state arm only")
  # Same reject on the single-season path, which shares the builder.
  expect_error(
    tobs(~ occ_cov, data = g$data, family = occu(),
         detection = .spde_int_det_term(), y = g$y[[1]],
         method = "nested_laplace",
         control = list(verbose = FALSE, progress = FALSE)),
    "attaches its latent blocks to the state arm only")
})

# The per-source list is an integrated-only shape; a single-season detection
# field stays the plain mesh vector its consumers read.
test_that("a single-season detection field keeps its vector shape", {
  skip_on_cran()
  skip_if_fast()
  skip_if_no_tulpamesh()

  sim <- .sim_spde_int_det(seed = 1, n_sites = 300, J1 = 8, J2 = 6)
  fit <- tobs(~ occ_cov, data = sim$data, family = occu(),
              detection = .spde_int_det_term(), y = sim$y[[1]],
              method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_true(is.numeric(fit$spatial_field_det))
  expect_equal(length(fit$spatial_field_det), fit$spatial$n_units)
})
