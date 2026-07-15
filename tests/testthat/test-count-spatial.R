# Areal count() -- a plain areal field (icar / bym2 / car_proper) on the
# abundance formula, routed to nested-Laplace (the spAbundance spAbund analogue).
# The field is a latent GMRF prior on the count GLMM block, integrated over its
# hyperparameters by the shared nested-Laplace EM machinery (no new C++).
# Poisson only: with one field node per site the negbin size / gaussian residual
# variance and the latent field are not jointly identified under the fixed-phi
# outer loop (gcol33/tulpaObs#117), so those are gated.
#
# Recovery-grade (per the "statistical code needs recovery tests" rule): the
# fixed effects recover to near-nominal 95% coverage across seeds and the latent
# field correlates with the simulated truth. The structural tests check the
# method registry and the dispatch gates.

# --- fixtures --------------------------------------------------------------

.count_grid_graph <- function(side) {
  N <- side * side
  A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

# Poisson counts with a smooth sum-to-zero areal field: log mu_i = X_i beta + f_i.
.count_sim_areal <- function(side = 10L, beta = c(0.5, 0.5), field_sd = 0.7,
                             seed = 1L) {
  set.seed(seed)
  A <- .count_grid_graph(side); N <- nrow(A)
  coord <- expand.grid(r = seq_len(side), c = seq_len(side))
  f <- field_sd * scale(sin(coord$r / side * pi) + cos(coord$c / side * pi))[, 1]
  f <- f - mean(f)
  data <- data.frame(x = stats::rnorm(N))
  X <- stats::model.matrix(~ x, data)
  y <- stats::rpois(N, exp(as.numeric(X %*% beta) + f))
  list(y = y, data = data, graph = A, field = f, beta = beta, N = N)
}


# --- (1) registry ----------------------------------------------------------

test_that("count() reports nested_laplace as a supported backend", {
  expect_true("nested_laplace" %in% tulpaObs:::.tobs_family_methods$count)
})


# --- (2) dispatch gates ----------------------------------------------------

test_that("areal count() gates the unsupported forms", {
  d <- .count_sim_areal(side = 6L, seed = 11L)

  # non-Poisson + areal: dispersion / field not jointly identified
  expect_error(
    tobs(~ x + icar(graph = d$graph), data = d$data, y = d$y,
         family = count("negbin"), method = "nested_laplace"),
    "not yet supported|not.*jointly identified|areal")
  expect_error(
    tobs(~ x + icar(graph = d$graph), data = d$data, y = d$y,
         family = count("gaussian"), method = "nested_laplace"),
    "not yet supported|not.*jointly identified|areal")

  # A varying-coefficient / weighted field IS wired (the svcAbund analogue,
  # gcol33/tulpaObs#120); its recovery lives in test-count-svc.R. A group_var
  # naming the field node is the identity map when the graph has one node per
  # site, so it fits rather than erroring.
  expect_s3_class(
    tobs(~ x + icar(graph = d$graph, weight = x, group_var = "cell"),
         data = cbind(d$data, cell = seq_len(d$N)), y = d$y,
         family = count(), method = "nested_laplace",
         control = list(verbose = FALSE, progress = FALSE)),
    "tobs_fit")

  # An AGGREGATING group_var (sites > cells) is not yet reconstructed here.
  expect_error(
    tobs(~ x + icar(graph = .count_grid_graph(3L), group_var = "cell"),
         data = cbind(d$data, cell = rep(seq_len(9L), length.out = d$N)),
         y = d$y, family = count(), method = "nested_laplace"),
    "group_var|sites > cells|one field node per site")

  # bym2 (mixed structured/unstructured field) is not reconstructed on the
  # generic nested path -> pointer to icar()/car_proper()
  expect_error(
    tobs(~ x + bym2(graph = d$graph), data = d$data, y = d$y,
         family = count(), method = "nested_laplace"),
    "icar.*car_proper|bym2")

  # engine mismatch either way
  expect_error(
    tobs(~ x, data = d$data, y = d$y, family = count(),
         method = "nested_laplace"),
    "needs a plain areal field")
  expect_error(
    tobs(~ x + icar(graph = d$graph), data = d$data, y = d$y,
         family = count(), method = "laplace"),
    "needs method = .nested_laplace.")
})


# --- (3) single fit + S3 ---------------------------------------------------

test_that("a Poisson areal count fit recovers the field and wires S3", {
  skip_on_cran()
  d   <- .count_sim_areal(side = 10L, seed = 7L)
  fit <- tobs(~ x + icar(graph = d$graph), data = d$data, y = d$y,
              family = count(), method = "nested_laplace",
              control = list(progress = FALSE, verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")
  expect_equal(length(fit$spatial_field), d$N)
  # the latent field tracks the simulated truth
  expect_gt(stats::cor(fit$spatial_field, d$field), 0.75)
  # the field SD hyperparameter is finite and positive
  expect_true(is.finite(fit$means[["sigma"]]) && fit$means[["sigma"]] > 0)

  # S3 surface runs
  expect_length(unlist(coef(fit)), 2L)
  expect_true(all(is.finite(diag(vcov(fit)))))
  expect_true(all(is.finite(confint(fit))))
  # fitted() is field-aware (in-sample), so it correlates with the counts
  ft <- fitted(fit)$mu
  expect_length(ft, d$N)
  expect_gt(stats::cor(ft, d$y), 0.5)
  # predict(newdata) is fixed-effect only (no field node at a new site)
  pr <- predict(fit, newdata = d$data)
  expect_length(pr, d$N)

  # WAIC is field-aware: the shared field is part of the log-mean, so the areal
  # fit scores better than the fixed-effect-only model on the same data (a wrong
  # FE-only WAIC would ignore the field and not improve).
  w_sp <- tobs_waic(fit)
  expect_true(is.finite(w_sp$waic))
  fit_fe <- tobs(~ x, data = d$data, y = d$y, family = count(),
                 control = list(progress = FALSE, verbose = FALSE))
  expect_lt(w_sp$waic, tobs_waic(fit_fe)$waic)
})


# --- (4) fixed-effect recovery + coverage (Poisson + ICAR) -----------------

test_that("Poisson areal count recovers coefficients with ~95% coverage", {
  skip_if_fast()
  skip_on_cran()
  beta   <- c(0.5, 0.5)
  n_seed <- 20L
  cover  <- matrix(FALSE, n_seed, length(beta))
  est    <- matrix(NA_real_, n_seed, length(beta))
  fcor   <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d   <- .count_sim_areal(side = 10L, beta = beta, seed = 400 + s)
    fit <- tobs(~ x + icar(graph = d$graph), data = d$data, y = d$y,
                family = count(), method = "nested_laplace",
                control = list(progress = FALSE, verbose = FALSE))
    b  <- unname(unlist(coef(fit)))
    se <- sqrt(diag(vcov(fit)))
    est[s, ]   <- b
    cover[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
    fcor[s]    <- stats::cor(fit$spatial_field, d$field)
  }
  expect_equal(colMeans(est), beta, tolerance = 0.08)
  expect_true(all(colMeans(cover) >= 0.85))
  expect_gt(stats::median(fcor), 0.8)
})


# --- (5) car_proper field kind ---------------------------------------------

test_that("areal count recovers the field under car_proper", {
  skip_if_fast()
  skip_on_cran()
  d <- .count_sim_areal(side = 10L, seed = 21L)

  fit_cp <- tobs(~ x + car_proper(graph = d$graph), data = d$data, y = d$y,
                 family = count(), method = "nested_laplace",
                 control = list(progress = FALSE, verbose = FALSE))
  expect_identical(fit_cp$method, "nested_laplace")
  expect_equal(length(fit_cp$spatial_field), d$N)
  expect_gt(stats::cor(fit_cp$spatial_field, d$field), 0.7)
})
