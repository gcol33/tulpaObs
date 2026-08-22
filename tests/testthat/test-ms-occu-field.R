# Community-spatial / SVC occupancy -- ms_occu() + a shared areal field or a
# varying-coefficient bar on the occupancy arm (the spOccupancy sfMsPGOcc /
# svcMsPGOcc analogue). A plain intercept field keeps the in-tree C++
# community-spatial path; a bar spatial(~ 1 + w || cell, graph) routes to the
# block-coordinate occupancy-field fitter (R/ms_occu_field.R): the community
# occupancy Laplace-EM with the field as a psi offset, alternated with a
# two-state-marginal field update. The field is informed by every species at
# each site, so it recovers alongside the community means.

.msof_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

# logit psi_{s,i} = X_i (mu + b_s) + f0_i + w_i f1_i; logit p = alpha_s.
.msof_sim <- function(side = 10L, S = 16L, J = 4L, svc = TRUE, seed = 1L) {
  set.seed(seed)
  A <- .msof_grid_graph(side); Ns <- nrow(A)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f0 <- 0.8 * scale(sin(co$r/side*pi) + cos(co$c/side*pi))[, 1]; f0 <- f0 - mean(f0)
  f1 <- 0.8 * scale(cos(co$r/side*pi*1.3) + sin(co$c/side*pi*0.7))[, 1]; f1 <- f1 - mean(f1)
  d <- data.frame(x = stats::rnorm(Ns), w = stats::rnorm(Ns), cell = seq_len(Ns))
  X <- stats::model.matrix(~ x, d)
  bocc <- vapply(1:2, function(j) stats::rnorm(S, c(0, 0.6)[j], c(0.4, 0.3)[j]),
                 numeric(S))
  alp  <- stats::rnorm(S, 0.4, 0.4)
  yarr <- array(NA_integer_, c(Ns, J, S))
  for (s in seq_len(S)) {
    eta <- as.numeric(X %*% bocc[s, ]) + f0 + if (svc) d$w * f1 else 0
    psi <- stats::plogis(eta); p <- stats::plogis(alp[s]); z <- stats::rbinom(Ns, 1, psi)
    for (i in seq_len(Ns)) yarr[i, , s] <- stats::rbinom(J, 1, z[i] * p)
  }
  dimnames(yarr) <- list(NULL, NULL, paste0("sp", seq_len(S)))
  list(y = yarr, data = d, graph = A, f0 = f0, f1 = f1, Ns = Ns)
}

test_that("svcMsPGOcc (occupancy SVC bar) recovers the intercept + trend fields", {
  skip_on_cran()
  d <- .msof_sim(side = 10L, S = 16L, seed = 3L)
  fit <- tobs(~ x + spatial(~ 1 + w || cell, graph = d$graph), data = d$data,
              family = ms_occu(), detection = ~ 1, y = d$y,
              species = paste0("sp", seq_len(16L)),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_false(is.null(fit$spatial_field))
  expect_false(is.null(fit$trend_field))
  expect_gt(stats::cor(fit$spatial_field, d$f0), 0.8)   # intercept field
  expect_gt(stats::cor(fit$trend_field,   d$f1), 0.8)   # varying-coefficient field
})

# Smoke coverage of the routing decision itself: a plain intercept field must not
# be diverted onto the block-coordinate R path, since the dedicated C++ fitter is
# the faster and already recovery-tested one. Asserted on a small grid, with the
# field-recovery threshold left to the gated block.
test_that("a plain intercept ms_occu field routes to the C++ spatial fitter", {
  # The cost here is the nested-Laplace outer grid, not the design, so the fixture
  # is trimmed on species and visits rather than on sites.
  d <- .msof_sim(side = 4L, S = 3L, J = 2L, svc = FALSE, seed = 4L)
  fit <- tobs(~ x + icar(graph = d$graph), data = d$data, family = ms_occu(),
              detection = ~ 1, y = d$y, species = paste0("sp", seq_len(3L)),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_length(fit$spatial_field, length(d$f0))
  # The C++ path reports no block-coordinate settings; the R driver always does.
  expect_null(fit$latent_control)
})

test_that("plain intercept ms_occu field keeps the C++ community-spatial path", {
  skip_if_fast()
  skip_on_cran()
  d <- .msof_sim(side = 10L, S = 16L, svc = FALSE, seed = 4L)
  fit <- tobs(~ x + icar(graph = d$graph), data = d$data, family = ms_occu(),
              detection = ~ 1, y = d$y, species = paste0("sp", seq_len(16L)),
              method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_gt(stats::cor(fit$spatial_field, d$f0), 0.8)
})

test_that("svcMsPGOcc recovers both occupancy fields over seeds", {
  skip_if_fast()
  skip_on_cran()
  n_seed <- 8L
  c0 <- numeric(n_seed); c1 <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d <- .msof_sim(side = 10L, S = 16L, seed = 100 + s)
    fit <- tobs(~ x + spatial(~ 1 + w || cell, graph = d$graph), data = d$data,
                family = ms_occu(), detection = ~ 1, y = d$y,
                species = paste0("sp", seq_len(16L)), method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    c0[s] <- stats::cor(fit$spatial_field, d$f0)
    c1[s] <- stats::cor(fit$trend_field,   d$f1)
  }
  expect_gt(stats::median(c0), 0.8)
  expect_gt(stats::median(c1), 0.8)
})
