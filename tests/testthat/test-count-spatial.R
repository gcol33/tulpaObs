# Areal count() -- a plain areal field (icar / bym2 / car_proper) on the
# abundance formula, routed to nested-Laplace (the spAbundance spAbund analogue).
# The field is a latent GMRF prior on the count GLMM block, integrated over its
# hyperparameters by the shared nested-Laplace EM machinery (no new C++). Poisson
# only: with one field node per site the negbin size / gaussian residual variance
# and the latent field are not jointly identified under the fixed-phi outer loop,
# so those are gated.
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

  # A varying-coefficient / weighted field IS wired (the svcAbund analogue);
  # its recovery lives in test-count-svc.R. A group_var naming the field node
  # is the identity map when the graph has one node per site, so it fits
  # rather than erroring.
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
  w_sp <- waic(fit)
  expect_true(is.finite(w_sp$waic))
  fit_fe <- tobs(~ x, data = d$data, y = d$y, family = count(),
                 control = list(progress = FALSE, verbose = FALSE))
  expect_lt(w_sp$waic, waic(fit_fe)$waic)
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

test_that("areal count recovers the field + slope under bym2", {
  skip_if_fast()
  skip_on_cran()
  # bym2 reconstructs the rho-mixed unit field z = sqrt(rho/scale) * phi +
  # sqrt(1 - rho) * theta on the generic nested-Laplace field summary (a distinct
  # path from icar/car_proper, which carry a single per-site block).
  fcor <- slopes <- numeric(3)
  for (s in seq_len(3)) {
    d <- .count_sim_areal(side = 10L, beta = c(0.5, 0.5), seed = 400 + s)
    fit <- tobs(~ x + bym2(graph = d$graph), data = d$data, y = d$y,
                family = count(), method = "nested_laplace",
                control = list(progress = FALSE, verbose = FALSE))
    if (s == 1L) {
      expect_identical(fit$method, "nested_laplace")
      expect_equal(length(fit$spatial_field), d$N)
      expect_true("sigma" %in% names(fit$means))
    }
    fcor[s]   <- stats::cor(fit$spatial_field, d$field)
    slopes[s] <- unname(unlist(coef(fit)))[2L]
  }
  # The bym2 field tracks the simulated (structured) truth, and the abundance
  # slope recovers on average.
  expect_gt(mean(fcor), 0.75)
  expect_lt(abs(mean(slopes) - 0.5), 0.15)
})


# --- (6) binomial areal field (spOccupancy svcPGBinom) ---------------------
# Unlike negbin / gaussian, a binomial areal count IS identified against a
# per-node field: the variance is pinned by the trial count n, so there is no
# free dispersion for the field to absorb.

# Binomial successes with a smooth sum-to-zero areal field:
# logit p_i = X_i beta + f_i, y_i ~ Binom(n_i, p_i).
.count_sim_areal_binom <- function(side = 12L, beta = c(0.2, 0.6),
                                    field_sd = 0.8, trials = 10L, seed = 1L) {
  set.seed(seed)
  A <- .count_grid_graph(side); N <- nrow(A)
  coord <- expand.grid(r = seq_len(side), c = seq_len(side))
  f <- field_sd * scale(sin(coord$r / side * pi) + cos(coord$c / side * pi))[, 1]
  f <- f - mean(f)
  data <- data.frame(x = stats::rnorm(N))
  X <- stats::model.matrix(~ x, data)
  p <- stats::plogis(as.numeric(X %*% beta) + f)
  y <- stats::rbinom(N, size = trials, prob = p)
  list(y = y, data = data, graph = A, field = f, beta = beta, N = N,
       trials = trials)
}

test_that("count('binomial') is identified against an areal field (not gated)", {
  d <- .count_sim_areal_binom(side = 8L, seed = 11L)
  # binomial does NOT hit the non-Poisson dispersion gate
  fit <- tobs(~ x + icar(graph = d$graph), data = d$data, y = d$y,
              family = count("binomial"), trials = d$trials,
              method = "nested_laplace",
              control = list(progress = FALSE, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")
  expect_equal(length(fit$spatial_field), d$N)
})

test_that("binomial areal count recovers the field + coefficients (trials>1)", {
  skip_if_fast()
  skip_on_cran()
  beta   <- c(0.2, 0.6)
  n_seed <- 20L
  cover  <- matrix(FALSE, n_seed, length(beta))
  est    <- matrix(NA_real_, n_seed, length(beta))
  fcor   <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    d   <- .count_sim_areal_binom(side = 12L, beta = beta, trials = 10L,
                                  seed = 600 + s)
    fit <- tobs(~ x + icar(graph = d$graph), data = d$data, y = d$y,
                family = count("binomial"), trials = d$trials,
                method = "nested_laplace",
                control = list(progress = FALSE, verbose = FALSE))
    b  <- unname(unlist(coef(fit)))
    se <- sqrt(diag(vcov(fit)))
    est[s, ]   <- b
    cover[s, ] <- (beta >= b - 1.96 * se) & (beta <= b + 1.96 * se)
    fcor[s]    <- stats::cor(fit$spatial_field, d$field)
  }
  expect_equal(colMeans(est), beta, tolerance = 0.08)
  expect_true(all(colMeans(cover) >= 0.85))
  # interior field: assert on field cor, never on sigma (the null-field trap)
  expect_gt(stats::median(fcor), 0.8)
})

test_that("Bernoulli areal count recovers the field (trials = 1, svcPGBinom)", {
  skip_if_fast()
  skip_on_cran()
  # One Bernoulli per node identifies the field but the fixed effects are weakly
  # identified against it (as in spOccupancy's svcPGBinom); assert field
  # recovery, which is the point of the family.
  fcor <- numeric(12L)
  for (s in seq_len(12L)) {
    d   <- .count_sim_areal_binom(side = 14L, trials = 1L, field_sd = 0.6,
                                  seed = 700 + s)
    fit <- tobs(~ x + icar(graph = d$graph), data = d$data, y = d$y,
                family = count("binomial"), trials = 1L,
                method = "nested_laplace",
                control = list(progress = FALSE, verbose = FALSE))
    fcor[s] <- stats::cor(fit$spatial_field, d$field)
  }
  expect_gt(stats::median(fcor), 0.6)
})

test_that("areal count recovers a continuous SPDE field + slope", {
  skip_if_fast()
  skip_on_cran()
  skip_if_no_tulpamesh()
  # A continuous-mesh spde() field on the count formula: the latent lives on the
  # mesh nodes (fit$spatial_field, length n_mesh) and the barycentric projector
  # fit$spatial$tulpa_spec$A maps it onto sites. The field summary reconstructs the
  # mesh field (grid-averaged block modes) and .count_spatial_field_offset projects
  # it to sites (A %*% mesh_field) for fitted() / predict().
  set.seed(42)
  n <- 300L
  coords <- cbind(runif(n), runif(n))
  u  <- 0.8 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])
  x  <- rnorm(n)
  y  <- rpois(n, exp(0.5 + 0.5 * x + u))
  dat <- data.frame(x = x, lon = coords[, 1], lat = coords[, 2])

  fit <- tobs(~ x + spde(lon, lat, max_edge = c(0.3, 0.6), nu = 1,
                         prior_range = c(0.3, 0.5), prior_sigma = c(0.7, 0.5)),
              data = dat, y = y, family = count(), method = "nested_laplace",
              control = list(progress = FALSE, verbose = FALSE))

  expect_identical(fit$method, "nested_laplace")
  # The mesh field lives on n_mesh nodes; the projector maps it to sites.
  expect_false(is.null(fit$spatial_field))
  expect_false(is.null(fit$spatial$tulpa_spec$A))
  expect_equal(length(fit$spatial_field), fit$spatial$n_units)
  field_at_sites <- as.numeric(fit$spatial$tulpa_spec$A %*% fit$spatial_field)
  expect_gt(cor(field_at_sites, u), 0.7)

  # Abundance slope recovers.
  expect_lt(abs(unname(unlist(coef(fit)))[2] - 0.5), 0.15)

  # fitted() projects the mesh field to sites (length n_sites, tracks y).
  ft <- fitted(fit)
  expect_equal(length(ft$mu), n)
  expect_gt(cor(ft$mu, y), 0.5)
})


test_that("areal count recovers a slope under a continuous NNGP gp() field", {
  skip_if_fast()
  skip_on_cran()
  # A continuous NNGP gp() field on the count formula routes to tulpa's single-
  # block `nngp` nested-Laplace kernel ( follow-up): the GP marginal variance and
  # range are integrated on the kernel's own outer grid and the field is
  # Schur-folded out, so the fit reports grid-integrated fixed effects plus the
  # GP hyperparameter posterior (fit$gp_hyper). The per-cell field itself is not
  # reconstructed on this path -- spde() is the route for a continuous field map.
  set.seed(7)
  n <- 150L
  coords <- cbind(runif(n), runif(n))
  D <- as.matrix(dist(coords))
  u <- as.numeric(t(chol(exp(-D / 0.25) + diag(1e-6, n))) %*% rnorm(n)); u <- u - mean(u)
  x <- rnorm(n)
  y <- rpois(n, exp(0.6 + 0.6 * x + u))
  dat <- data.frame(x = x, lon = coords[, 1], lat = coords[, 2])

  fit <- tobs(~ x + gp(lon, lat, prior_range = c(0.1, 0.05)),
              data = dat, y = y, family = count("poisson"),
              method = "nested_laplace",
              control = list(progress = FALSE, verbose = FALSE))

  expect_identical(fit$method, "nested_laplace")
  # Slope recovers under the field-integrated marginal.
  expect_lt(abs(unname(unlist(coef(fit)))[2] - 0.6), 0.15)
  # The GP hyperparameter posterior (marginal variance + range) is surfaced.
  expect_true(all(c("sigma2", "phi_gp") %in% fit$gp_hyper$theta_names))
  expect_true(all(is.finite(fit$gp_hyper$mean)))
  # The field is integrated out (no per-cell map on this path).
  expect_null(fit$spatial_field)
  # S3 surface: fitted / predict / WAIC.
  expect_equal(length(fitted(fit)$mu), n)
  expect_length(predict(fit), n)
  expect_true(is.finite(waic(fit)$waic))
})

test_that("count() spatial-field gates: multiscale_gp and negbin + gp", {
  skip_on_cran()
  set.seed(9); n <- 50L
  dat <- data.frame(x = rnorm(n), lon = runif(n), lat = runif(n))
  y <- rpois(n, exp(0.5 + 0.4 * dat$x))
  expect_error(
    tobs(~ x + multiscale_gp(lon, lat), data = dat, y = y,
         family = count("poisson"), method = "nested_laplace",
         control = list(verbose = FALSE)),
    "multiscale_gp")
  expect_error(
    tobs(~ x + gp(lon, lat, prior_range = c(0.1, 0.05)), data = dat, y = y,
         family = count("negbin"), method = "nested_laplace",
         control = list(verbose = FALSE)),
    "identifiable|confounded")
})
