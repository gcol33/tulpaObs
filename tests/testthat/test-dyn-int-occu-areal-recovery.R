# =============================================================================
# test-dyn-int-occu-areal-recovery.R - shared areal ICAR field on the first-
# season occupancy (psi1) arm of the multi-season integrated occupancy model
# (spOccupancy stIntPGOcc; gcol33/tulpaObs#122).
#
# psi1 sets only the initial mixing weight of each site's HMM, so the per-site
# field gradient is the Fisher-identity score w1 - psi1 (the smoothed season-1
# occupancy), which drives the shared areal-BFGS nested-Laplace fit. These check
# the routing / gates and, over seeds, recovery of the interior field (by cor,
# never sigma -- the null-field trap) plus the transition rates.
# =============================================================================

.dio_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}
.dio_field <- function(side, sd = 0.9) {
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f  <- sd * scale(sin(co$r / side * pi) + cos(co$c / side * pi))[, 1]
  f - mean(f)
}

test_that("dyn_int_occu() nested_laplace field is registered and gated", {
  expect_true("nested_laplace" %in% tulpaObs:::.tobs_family_methods$dyn_int_occu)
  side <- 5L; A <- .dio_grid_graph(side)
  sim <- simulate_dyn_int_occu(N = side * side, T_seasons = 4, S = 2, J = 3,
                               field = .dio_field(side), seed = 1)
  # a field needs nested_laplace; plain laplace with a field errors with a pointer
  expect_error(
    tobs(~ 1 + icar(graph = A), data = sim$data, family = dyn_int_occu(),
         detection = ~ 1, colonization = ~ 1, extinction = ~ 1, y = sim$y,
         sources = sim$sources, method = "laplace"),
    "nested_laplace")
  # bym2 / car_proper are gated to icar in v1
  expect_error(
    tobs(~ 1 + bym2(graph = A), data = sim$data, family = dyn_int_occu(),
         detection = ~ 1, colonization = ~ 1, extinction = ~ 1, y = sim$y,
         sources = sim$sources, method = "nested_laplace"),
    "icar")
})

test_that("dyn_int_occu() + icar recovers the shared field + transitions", {
  skip_if_fast()
  skip_on_cran()
  # stIntPGOcc: the field is the new object; assert field recovery by cor on an
  # INTERIOR field, plus the transition rates in aggregate.
  side <- 9L; N <- side * side; A <- .dio_grid_graph(side)
  ftrue <- .dio_field(side, sd = 0.9)
  n_seed <- 12L
  fcor <- rep(NA_real_, n_seed)
  gm <- ep <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_dyn_int_occu(N = N, T_seasons = 4, S = 2, J = 3,
                                 psi1 = 0.5, gamma = 0.3, eps = 0.2,
                                 p = c(0.45, 0.6), field = ftrue, seed = 200 + s)
    fit <- tryCatch(
      tobs(~ 1 + icar(graph = A), data = sim$data, family = dyn_int_occu(),
           detection = ~ 1, colonization = ~ 1, extinction = ~ 1, y = sim$y,
           sources = sim$sources, method = "nested_laplace",
           control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    fcor[s] <- stats::cor(fit$spatial_field, ftrue)
    gm[s]   <- stats::plogis(fit$means[["gamma_(Intercept)"]])
    ep[s]   <- stats::plogis(fit$means[["eps_(Intercept)"]])
  }
  # interior field recovery (cor), pooled across all sources/seasons.
  expect_gt(stats::median(fcor[!is.na(fcor)]), 0.7)
  # transition rates recover in aggregate.
  expect_lt(abs(mean(gm, na.rm = TRUE) - 0.30), 0.06)
  expect_lt(abs(mean(ep, na.rm = TRUE) - 0.20), 0.06)
})

# --- varying-coefficient bar on the psi1 arm (svcTIntPGOcc, #122) ---

# A second interior field, shifted from .dio_field so the intercept and trend
# surfaces are not collinear (which would make them unidentifiable).
.dio_field2 <- function(side, sd = 0.9) {
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f  <- sd * scale(cos(co$r / side * pi * 1.3) + sin(co$c / side * pi * 0.8))[, 1]
  f - mean(f)
}

test_that("dyn_int_occu() + SVC bar recovers intercept + trend fields", {
  skip_if_fast()
  skip_on_cran()
  # svcTIntPGOcc: an intercept field PLUS a covariate-weighted (varying-coefficient)
  # field on the first-season occupancy arm, both shared. The weighted areal block
  # flows through the same multi-block areal-BFGS driver as stIntPGOcc; the psi1
  # score w1 - psi1 scatters to each block (weighted by the covariate for the trend
  # block). Assert BOTH interior fields recover by cor (never sigma; null-field trap).
  side <- 9L; N <- side * side; A <- .dio_grid_graph(side)
  f0 <- .dio_field(side, sd = 0.9); f1 <- .dio_field2(side, sd = 0.9)
  n_seed <- 12L
  c0 <- c1 <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_dyn_int_occu(N = N, T_seasons = 4, S = 2, J = 3,
                                 psi1 = 0.5, gamma = 0.3, eps = 0.2,
                                 p = c(0.45, 0.6), field = f0, trend = f1,
                                 seed = 600 + s)
    fit <- tryCatch(
      tobs(~ spatial(~ 1 + w || cell, graph = A), data = sim$data,
           family = dyn_int_occu(), detection = ~ 1, colonization = ~ 1,
           extinction = ~ 1, y = sim$y, sources = sim$sources,
           method = "nested_laplace", control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    c0[s] <- stats::cor(fit$spatial_field, f0)
    c1[s] <- if (!is.null(fit$trend_field)) stats::cor(fit$trend_field, f1) else NA_real_
  }
  # Both the intercept and the varying-coefficient surface recover; the SVC field
  # is more weakly identified (weighted loading, jointly separated from the
  # intercept), so a looser floor than the intercept field.
  expect_gt(stats::median(c0[!is.na(c0)]), 0.75)
  expect_gt(stats::median(c1[!is.na(c1)]), 0.55)
})
