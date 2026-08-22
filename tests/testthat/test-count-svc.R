# Single-species spatially-varying-coefficient count -- count() + a weighted areal
# bar (the spAbundance svcAbund analogue).
#
#   log mu_i = X_i beta + f0_{cell(i)} + w_i * f1_{cell(i)}
#
# An intercept field plus one varying-coefficient field per covariate on the
# abundance formula. No new engine: .tobs_to_multi_block_prior already emits one
# weighted latent block per resolved field for model_type "count" (the weighted
# ones carry a per-site svc_weight), and .tobs_nested_attach_field_summary
# already loops the areal blocks into spatial_field + trend_fields. The fit
# records the per-site field contribution sum_k W[i,k] f_k[i] as
# model$count_field_offset so fitted() adds the FULL contribution rather than the
# intercept field alone.
#
# Recovery-grade (per the "statistical code needs recovery tests" rule): both
# latent surfaces correlate with the simulated truth across seeds. Thresholds are
# set from a measured run (dev_notes/probe_svc_abund.R): intercept field
# 0.93-0.96, trend field 0.85-0.91 at side = 12 (144 sites).

# --- fixtures --------------------------------------------------------------

.svca_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

# Two distinct smooth sum-to-zero surfaces: f0 (intercept) and f1 (the
# coefficient on w, so it enters eta weighted by w_i).
.svca_sim <- function(side = 12L, beta = c(0.6, 0.4), f0_sd = 0.6, f1_sd = 0.5,
                      seed = 1L) {
  set.seed(seed)
  A  <- .svca_grid_graph(side); N <- nrow(A)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f0 <- f0_sd * scale(sin(co$r / side * pi) + cos(co$c / side * pi))[, 1]
  f0 <- f0 - mean(f0)
  f1 <- f1_sd * scale(cos(co$r / side * 2 * pi) - sin(co$c / side * pi))[, 1]
  f1 <- f1 - mean(f1)
  d  <- data.frame(x = stats::rnorm(N), w = stats::rnorm(N), cell = seq_len(N))
  X  <- stats::model.matrix(~ x, d)
  d$y <- stats::rpois(N, exp(as.numeric(X %*% beta) + f0 + d$w * f1))
  list(data = d, graph = A, f0 = f0, f1 = f1, beta = beta, N = N)
}

.svca_fit <- function(d) {
  A <- d$graph
  tobs(y ~ x + spatial(~ 1 + w || cell, graph = A), data = d$data,
       family = count(), method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}


# --- (1) dispatch gates ----------------------------------------------------

test_that("count() SVC gates the unsupported field kinds", {
  d <- .svca_sim(side = 6L, seed = 11L)
  A <- d$graph

  # bym2 mixes a structured + unstructured component with hyperparameter-
  # dependent scales, so its per-cell field is not reconstructed on this path
  expect_error(
    tobs(y ~ x + spatial(~ 1 + w || cell, graph = A, type = "bym2"),
         data = d$data, family = count(), method = "nested_laplace"),
    "icar|car_proper|bym2")
  # a varying-coefficient field still needs nested_laplace
  expect_error(
    tobs(y ~ x + spatial(~ 1 + w || cell, graph = A), data = d$data,
         family = count(), method = "laplace"),
    "nested_laplace")
  # areal count is Poisson-only (the field absorbs all overdispersion)
  expect_error(
    tobs(y ~ x + spatial(~ 1 + w || cell, graph = A), data = d$data,
         family = count("negbin"), method = "nested_laplace"),
    "identifiable|Poisson|areal")
})


# --- (2) recovery ----------------------------------------------------------

test_that("svcAbund recovers the intercept and the varying-coefficient field", {
  skip_on_cran()
  d   <- .svca_sim(seed = 1L)
  fit <- .svca_fit(d)

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")

  # both surfaces present and recovered
  expect_length(fit$spatial_field, d$N)
  expect_false(is.null(fit$trend_field))
  expect_length(fit$trend_fields, 1L)
  expect_gt(stats::cor(fit$spatial_field, d$f0), 0.85)
  expect_gt(stats::cor(fit$trend_field,   d$f1), 0.75)

  # one SD hyperparameter per field
  expect_true(all(c("sigma", "sigma_trend") %in% names(fit$means)))

  # fitted() adds the FULL weighted field contribution, not the intercept field
  # alone -- so the recorded offset must differ from spatial_field.
  off <- fit$model$count_field_offset
  expect_length(off, d$N)
  expect_false(isTRUE(all.equal(as.numeric(off),
                                as.numeric(fit$spatial_field))))
  expect_gt(stats::cor(fitted(fit)$mu, d$data$y), 0.6)
})

test_that("svcAbund recovers both fields across seeds", {
  skip_if_fast()
  skip_on_cran()
  n_seed <- 5L
  c0 <- numeric(n_seed); c1 <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d   <- .svca_sim(seed = 100L + s)
    fit <- .svca_fit(d)
    c0[s] <- stats::cor(fit$spatial_field, d$f0)
    c1[s] <- stats::cor(fit$trend_field,   d$f1)
  }
  expect_gt(stats::median(c0), 0.90)
  expect_gt(stats::median(c1), 0.80)
})
