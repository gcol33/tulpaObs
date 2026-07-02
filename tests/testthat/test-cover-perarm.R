# =============================================================================
# test-cover-perarm.R - cover() per-arm formulas (arm = formula).
#
# `cover(presence = ~ ..., positive = ~ ...)` gives each hurdle arm its own fixed
# effects (two independent designs), while the single shared `formula` stays the
# back-compat spelling. Structured terms in a per-arm formula are rejected in this
# first cut (fields go on the shared formula with `to =`).
# =============================================================================

.cp_sim <- function(seed = 1L, N = 3000L) {
  set.seed(seed)
  x1 <- stats::rnorm(N); x2 <- stats::rnorm(N)
  pres <- stats::rbinom(N, 1L, stats::plogis(-0.2 + 1.0 * x1))
  mu   <- stats::plogis(-0.5 + 0.8 * x2)
  y <- numeric(N); pos <- pres == 1L
  y[pos] <- stats::rbeta(sum(pos), mu[pos] * 20, (1 - mu[pos]) * 20)
  y[y >= 1] <- 1 - 1e-6
  list(data = data.frame(x1 = x1, x2 = x2), y = y)
}

test_that("per-arm formulas give each arm its own fixed effects and recover", {
  skip_on_cran()
  s <- .cp_sim()
  fit <- tobs(presence = ~ x1, positive = ~ x2,
              family = cover(positive = "beta"), data = s$data, y = s$y,
              method = "laplace")

  # Each arm carries ONLY its own covariate.
  expect_true("x1" %in% names(fit$beta_occ) && !("x2" %in% names(fit$beta_occ)))
  expect_true("x2" %in% names(fit$beta_pos) && !("x1" %in% names(fit$beta_pos)))
  # ... and recovers its truth.
  expect_lt(abs(fit$beta_occ[["x1"]] - 1.0), 0.15)
  expect_lt(abs(fit$beta_pos[["x2"]] - 0.8), 0.15)
})

test_that("the shared single formula stays back-compat (both arms share the FE)", {
  skip_on_cran()
  s <- .cp_sim()
  fit <- tobs(~ x1 + x2, family = cover(positive = "beta"),
              data = s$data, y = s$y, method = "laplace")
  expect_true(all(c("x1", "x2") %in% names(fit$beta_occ)))
  expect_true(all(c("x1", "x2") %in% names(fit$beta_pos)))
})

test_that("a spatial field placed in an arm formula == the to= spelling", {
  skip_on_cran()
  side <- 6L; nc <- side * side
  adj <- matrix(0L, nc, nc)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    if (r > 1L)   adj[idx(r, c), idx(r - 1L, c)] <- 1L
    if (r < side) adj[idx(r, c), idx(r + 1L, c)] <- 1L
    if (c > 1L)   adj[idx(r, c), idx(r, c - 1L)] <- 1L
    if (c < side) adj[idx(r, c), idx(r, c + 1L)] <- 1L
  }
  set.seed(3); reps <- 10L; cell <- rep(seq_len(nc), each = reps); N <- length(cell)
  x <- stats::rnorm(N); tt <- stats::rnorm(N)
  pres <- stats::rbinom(N, 1L, stats::plogis(0.2 + 0.5 * x))
  mu <- stats::plogis(0.3 * tt)
  y <- numeric(N); pos <- pres == 1L
  y[pos] <- stats::rbeta(sum(pos), mu[pos] * 15, (1 - mu[pos]) * 15); y[y >= 1] <- 1 - 1e-6
  dat <- data.frame(x = x, t = tt, cell = cell)
  ctrl <- list(progress = FALSE, integration = "ccd")

  fit_pa <- suppressWarnings(tobs(
    presence = ~ x + t, positive = ~ x + t + spatial(~ 1 || cell, graph = adj),
    family = cover(positive = "beta"), data = dat, y = y,
    method = "nested_laplace", control = ctrl))
  fit_to <- suppressWarnings(tobs(
    ~ x + t + spatial(~ 1 || cell, graph = adj, to = "positive"),
    family = cover(positive = "beta"), data = dat, y = y,
    method = "nested_laplace", control = ctrl))

  expect_equal(unname(fit_pa$beta_occ), unname(fit_to$beta_occ), tolerance = 1e-8)
  expect_equal(unname(fit_pa$beta_pos), unname(fit_to$beta_pos), tolerance = 1e-8)
})

test_that("a temporal() / re() term in a per-arm formula is rejected with a pointer", {
  s <- .cp_sim(N = 200L)
  expect_error(
    tobs(presence = ~ x1 + (1 | x2), positive = ~ x2,
         family = cover(positive = "beta"), data = s$data, y = s$y,
         method = "laplace"),
    "not routed per-arm")
})

test_that("only one per-arm formula errors (need both, or the shared one)", {
  s <- .cp_sim(N = 200L)
  expect_error(
    tobs(presence = ~ x1, family = cover(positive = "beta"),
         data = s$data, y = s$y, method = "laplace"),
    "BOTH")
})
