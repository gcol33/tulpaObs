# =============================================================================
# test-jsdm-spatial.R - areal-spatial joint species distribution model
# (jsdm() + shared field; tulpaObs#76).
#
# The JSDM observes presence/absence directly (no detection process), with shared
# fixed-effect coefficients, a scalar per-species random intercept, and one shared
# ICAR / BYM2 / proper-CAR field on the latent occupancy, fit by the in-tree
# JSDM-spatial nested Laplace-EM (R/jsdm_spatial.R, src/jsdm_spatial.cpp).
#
# Coverage: (1) the per-(species, site) Bernoulli cell score + negative Hessian
# byte-exact vs the closed form (FD); (2) fixed-effect + shared-field recovery
# (ICAR) over several seeds; (3) the proper-CAR + BYM2 field kinds run + recover
# the field; (4) S3 (coef / vcov / confint / ranef / spatial field); (5) gates
# (engine routing, field kind, one-unit-per-site).
# =============================================================================


# --- shared fixtures -------------------------------------------------------

.jsds_grid_graph <- function(side) {
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

# JSDM presence/absence with a smooth shared field on the latent occupancy and a
# per-species random intercept. logit psi_{s,i} = X_i . beta + b_s + f_{u(i)}.
.jsds_sim <- function(side = 8L, n_species = 14L, beta = c(-0.2, 0.9),
                      sd_re = 0.5, field_sd = 0.9, seed = 1L) {
  set.seed(seed)
  A <- .jsds_grid_graph(side); N <- nrow(A)
  coord <- expand.grid(r = seq_len(side), c = seq_len(side))
  f <- field_sd * scale(sin(coord$r / side * pi) + cos(coord$c / side * pi))[, 1]
  f <- f - mean(f)
  data <- data.frame(x = stats::rnorm(N))
  X <- stats::model.matrix(~ x, data)
  b_s <- stats::rnorm(n_species, 0, sd_re)
  y <- matrix(NA_integer_, N, n_species)
  for (s in seq_len(n_species)) {
    psi <- stats::plogis(as.numeric(X %*% beta) + b_s[s] + f)
    y[, s] <- stats::rbinom(N, 1, psi)
  }
  colnames(y) <- paste0("sp", seq_len(n_species))
  list(y = y, data = data, graph = A, field = f, beta = beta, b_s = b_s,
       species = paste0("sp", seq_len(n_species)))
}


# --- (1) per-cell Bernoulli score + negative Hessian vs FD -----------------

test_that("jsdm spatial cell score + curvature match finite differences", {
  skip_on_cran()
  ll <- function(eta, y) {
    psi <- pmin(pmax(stats::plogis(eta), 1e-12), 1 - 1e-12)
    if (y == 1L) log(psi) else log1p(-psi)
  }
  max_g <- 0; max_h <- 0
  for (eta in c(-2.5, -0.7, 0, 0.6, 1.8)) for (y in c(0L, 1L)) {
    cell <- cpp_jsdm_site_cell(eta, y)
    h <- 1e-6
    g_fd <- (ll(eta + h, y) - ll(eta - h, y)) / (2 * h)
    hh <- 1e-4
    H_fd <- -(ll(eta + hh, y) - 2 * ll(eta, y) + ll(eta - hh, y)) / hh^2
    max_g <- max(max_g, abs(cell$grad - g_fd))
    max_h <- max(max_h, abs(cell$neg_hess - H_fd))
  }
  expect_lt(max_g, 1e-6)
  expect_lt(max_h, 1e-4)
})


# --- (2) fixed-effect + shared-field recovery (ICAR), several seeds --------

test_that("jsdm + icar() recovers fixed effects + the shared field", {
  skip_on_cran()
  skip_if_fast()
  cors <- numeric(0); b0 <- numeric(0); b1 <- numeric(0)
  for (sd in 1:4) {
    sim <- .jsds_sim(side = 8L, n_species = 14L, seed = sd)
    fit <- tobs(~ x + icar(graph = sim$graph), data = sim$data, family = jsdm(),
                y = sim$y, species = sim$species, method = "nested_laplace",
                control = list(verbose = FALSE))
    expect_equal(fit$method, "nested_laplace")
    expect_false(is.null(fit$spatial_field))
    cf <- coef(fit)$psi
    cors <- c(cors, cor(fit$spatial_field, sim$field))
    b0 <- c(b0, cf[["(Intercept)"]]); b1 <- c(b1, cf[["x"]])
  }
  # Shared latent field recovered (up to the sum-to-zero constraint).
  expect_gt(mean(cors), 0.80)
  expect_true(all(cors > 0.70))
  # Fixed-effect community means recovered (small finite-sample / attenuation
  # bias on the intercept; the x slope is the discriminating quantity).
  expect_lt(abs(mean(b0) - (-0.2)), 0.25)
  expect_lt(abs(mean(b1) - 0.9), 0.20)
})


# --- (3) proper-CAR + BYM2 field kinds run + recover the field -------------

test_that("jsdm + car_proper()/bym2() recover the shared field", {
  skip_on_cran()
  skip_if_fast()
  sim <- .jsds_sim(side = 8L, n_species = 14L, seed = 2L)
  for (term in c("car_proper", "bym2")) {
    f <- stats::as.formula(sprintf("~ x + %s(graph = sim$graph)", term))
    fit <- tobs(f, data = sim$data, family = jsdm(), y = sim$y,
                species = sim$species, method = "nested_laplace",
                control = list(verbose = FALSE))
    expect_equal(fit$method, "nested_laplace")
    expect_gt(cor(fit$spatial_field, sim$field), 0.70)
    expect_lt(abs(coef(fit)$psi[["x"]] - 0.9), 0.25)
  }
})


# --- (4) S3 + per-species random intercept ---------------------------------

test_that("jsdm spatial S3 works (coef / vcov / confint / ranef / field)", {
  skip_on_cran()
  skip_if_fast()
  sim <- .jsds_sim(side = 7L, n_species = 12L, sd_re = 0.6, seed = 5L)
  fit <- tobs(~ x + icar(graph = sim$graph), data = sim$data, family = jsdm(),
              y = sim$y, species = sim$species, method = "nested_laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))

  cf <- coef(fit)
  expect_named(cf, "psi")
  expect_named(cf$psi, c("(Intercept)", "x"))

  V <- vcov(fit)
  expect_equal(dim(V), c(2L, 2L))
  expect_true(all(diag(V) > 0))

  ci <- confint(fit)
  expect_equal(nrow(ci), 2L)

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 12L)

  # Per-species random intercepts track the simulated deviations.
  expect_gt(cor(fit$jsdm_re$blup, sim$b_s), 0.6)
  expect_length(fit$spatial_field, nrow(sim$graph))
})


# --- (5) gates -------------------------------------------------------------

test_that("jsdm spatial gates: engine routing, field kind, one-unit-per-site", {
  skip_on_cran()
  sim <- .jsds_sim(side = 5L, n_species = 4L, seed = 1L)

  # A field on the formula needs method = "nested_laplace".
  expect_error(
    tobs(~ x + icar(graph = sim$graph), data = sim$data, family = jsdm(),
         y = sim$y, species = sim$species, method = "laplace"),
    "nested_laplace")

  # NUTS + a field points to nested_laplace.
  expect_error(
    tobs(~ x + icar(graph = sim$graph), data = sim$data, family = jsdm(),
         y = sim$y, species = sim$species, method = "nuts"),
    "nested_laplace")

  # nested_laplace without a field errors with a pointer to the non-spatial path.
  expect_error(
    tobs(~ x, data = sim$data, family = jsdm(), y = sim$y,
         species = sim$species, method = "nested_laplace"),
    "areal field")

  # car() is the improper non-intrinsic CAR; the areal JSDM path takes
  # icar / bym2 / car_proper only.
  expect_error(
    tobs(~ x + car(graph = sim$graph), data = sim$data, family = jsdm(),
         y = sim$y, species = sim$species, method = "nested_laplace"),
    "icar")
})
