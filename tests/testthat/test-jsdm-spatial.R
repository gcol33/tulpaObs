# =============================================================================
# test-jsdm-spatial.R - areal-spatial joint species distribution model
# (jsdm() + a shared field; gcol33/tulpaObs#76, #121).
#
# jsdm() observes presence/absence directly (no detection process). Since #121 it
# is the COMMUNITY GLMM -- per-species coefficients under a Gaussian community
# covariance (the spOccupancy lfJSDM / sfJSDM model class), i.e. ms_count() with a
# logit link -- so a shared ICAR / proper-CAR / BYM2 field is fit by the block
# coordinate driver (R/community_latent.R): the community Laplace-EM with the
# field as a per-site offset, alternated with the areal Newton.
#
# Coverage here: the three areal field kinds recover the shared field alongside
# the community means, plus S3 and the dispatch gates. The latent-factor
# compositions (lfJSDM / sfJSDM), the community-model structure, and the NUTS
# target live in test-jsdm-community.R.
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

# JSDM presence/absence with a smooth shared field on the latent occupancy and
# per-species coefficient deviations under a community covariance:
# logit psi_{s,i} = X_i . (mu + b_s) + f_{u(i)}.
.jsds_sim <- function(side = 10L, n_species = 14L, beta = c(-0.2, 0.9),
                      sd_re = c(0.5, 0.3), field_sd = 0.9, seed = 1L) {
  set.seed(seed)
  A <- .jsds_grid_graph(side); N <- nrow(A)
  coord <- expand.grid(r = seq_len(side), c = seq_len(side))
  f <- field_sd * scale(sin(coord$r / side * pi) + cos(coord$c / side * pi))[, 1]
  f <- f - mean(f)
  data <- data.frame(x = stats::rnorm(N))
  X <- stats::model.matrix(~ x, data)
  bs <- vapply(1:2, function(j) stats::rnorm(n_species, beta[j], sd_re[j]),
               numeric(n_species))
  y <- matrix(NA_integer_, N, n_species)
  for (s in seq_len(n_species)) {
    y[, s] <- stats::rbinom(N, 1, stats::plogis(as.numeric(X %*% bs[s, ]) + f))
  }
  colnames(y) <- paste0("sp", seq_len(n_species))
  list(y = y, data = data, graph = A, field = f, beta = beta, bs = bs,
       species = paste0("sp", seq_len(n_species)))
}


# --- (1) shared-field recovery (ICAR), several seeds -----------------------

test_that("jsdm + icar() recovers community means + the shared field", {
  skip_on_cran()
  skip_if_fast()
  cors <- numeric(0); b0 <- numeric(0); b1 <- numeric(0)
  for (sd in 1:4) {
    sim <- .jsds_sim(seed = sd)
    fit <- tobs(~ x + icar(graph = sim$graph), data = sim$data, family = jsdm(),
                y = sim$y, species = sim$species, method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    expect_equal(fit$method, "nested_laplace")
    expect_false(is.null(fit$spatial_field))
    cors <- c(cors, cor(fit$spatial_field, sim$field))
    b0 <- c(b0, unname(fit$means[1L])); b1 <- c(b1, unname(fit$means[2L]))
  }
  # Shared latent field recovered (up to the sum-to-zero constraint).
  expect_gt(mean(cors), 0.80)
  expect_true(all(cors > 0.70))
  # Community-mean coefficients recovered (a small finite-sample bias on the
  # intercept, which trades against the field level; the slope discriminates).
  expect_lt(abs(mean(b0) - (-0.2)), 0.30)
  expect_lt(abs(mean(b1) - 0.9), 0.25)
})


# --- (2) proper-CAR + BYM2 field kinds run + recover the field -------------

test_that("jsdm + car_proper()/bym2() recover the shared field", {
  skip_on_cran()
  skip_if_fast()
  sim <- .jsds_sim(seed = 2L)
  for (term in c("car_proper", "bym2")) {
    f <- stats::as.formula(sprintf("~ x + %s(graph = sim$graph)", term))
    fit <- tobs(f, data = sim$data, family = jsdm(), y = sim$y,
                species = sim$species, method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    expect_equal(fit$method, "nested_laplace")
    expect_gt(cor(fit$spatial_field, sim$field), 0.70)
    expect_lt(abs(unname(fit$means[2L]) - 0.9), 0.30)
  }
})


# --- (3) S3 + per-species coefficient deviations ---------------------------

test_that("jsdm spatial S3 works (coef / vcov / confint / ranef / field)", {
  skip_on_cran()
  skip_if_fast()
  sim <- .jsds_sim(side = 9L, n_species = 12L, seed = 5L)
  fit <- tobs(~ x + icar(graph = sim$graph), data = sim$data, family = jsdm(),
              y = sim$y, species = sim$species, method = "nested_laplace",
              control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))

  V <- vcov(fit)
  expect_equal(dim(V), c(2L, 2L))
  expect_true(all(diag(V) > 0))

  ci <- confint(fit)
  expect_equal(nrow(ci), 2L)

  # ranef(): per-species deviations, one row per (species, term).
  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 12L * 2L)

  # The per-species coefficients track the simulated per-species truth.
  expect_gt(cor(fit$ms_community$coef_mu[, "x"], sim$bs[, 2L]), 0.5)
  expect_length(fit$spatial_field, nrow(sim$graph))
})


# --- (4) gates -------------------------------------------------------------

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

  # car() is the improper non-intrinsic CAR; the areal community path takes
  # icar / car_proper / bym2 only.
  expect_error(
    tobs(~ x + car(graph = sim$graph), data = sim$data, family = jsdm(),
         y = sim$y, species = sim$species, method = "nested_laplace"),
    "icar|car_proper|bym2")
})
