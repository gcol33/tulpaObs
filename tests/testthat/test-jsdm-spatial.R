# =============================================================================
# test-jsdm-spatial.R - areal-spatial joint species distribution model
# (jsdm() + a shared field).
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


# JSDM presence/absence with a smooth shared field on the latent occupancy and
# per-species coefficient deviations under a community covariance:
# logit psi_{s,i} = X_i . (mu + b_s) + f_{u(i)}.
.jsds_sim <- function(side = 10L, n_species = 14L, beta = c(-0.2, 0.9),
                      sd_re = c(0.5, 0.3), field_sd = 0.9, seed = 1L) {
  set.seed(seed)
  A <- rook_adj(side); N <- nrow(A)
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
  cors <- numeric(0); dev0 <- numeric(0); dev1 <- numeric(0)
  for (sd in 1:4) {
    sim <- .jsds_sim(seed = sd)
    fit <- tobs(~ x + icar(graph = sim$graph), data = sim$data, family = jsdm(),
                y = sim$y, species = sim$species, method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    expect_equal(fit$method, "nested_laplace")
    expect_false(is.null(fit$spatial_field))
    cors <- c(cors, cor(fit$spatial_field, sim$field))
    real <- colMeans(sim$bs)
    dev0 <- c(dev0, unname(fit$means[1L]) - real[1L])
    dev1 <- c(dev1, unname(fit$means[2L]) - real[2L])
  }
  # Shared latent field recovered (up to the sum-to-zero constraint).
  expect_gt(mean(cors), 0.80)
  expect_true(all(cors > 0.70))
  # Community means against the seed's REALIZED mean (colMeans(sim$bs)), not the
  # nominal c(-0.2, 0.9): this loop already averaged over seeds before comparing
  # to the nominal, so calls it statistically valid as it stood; retargeting is
  # a power improvement here, not a bug fix. Budget is 5x the SE of a 4-seed
  # mean, from a fresh 16-seed measurement of this exact fixture (dev0/dev1 =
  # fit$means[1:2] - colMeans(sim$bs), seeds 1-16): intercept sd 0.054 -> SE_4
  # 0.027 -> budget 0.135; slope sd 0.052 -> SE_4 0.026 -> budget 0.130.
  #
  # The slope carries a real, one-sided finite-sample bias, not just draw
  # noise: over the 16-seed measurement the mean deviation is -0.076 (se 0.013,
  # ~5.8 se from zero, 15 of 16 seeds negative) -- the shared field absorbs part
  # of the slope's signal under the sum-to-zero constraint. On the seeds this
  # loop actually runs (1-4) the mean deviation is -0.031, comfortably inside
  # the budget; the budget is sized with margin above the larger 16-seed figure
  # rather than tuned to only this loop's own smaller sample, so it does not
  # silently paper over the known bias.
  expect_lt(abs(mean(dev0)), 0.135)
  expect_lt(abs(mean(dev1)), 0.130)
})


# --- (2) proper-CAR + BYM2 field kinds run + recover the field -------------

test_that("jsdm + car_proper()/bym2() recover the shared field", {
  skip_on_cran()
  skip_if_fast()
  sim <- .jsds_sim(seed = 2L)
  real_slope <- colMeans(sim$bs)[2L]
  for (term in c("car_proper", "bym2")) {
    f <- stats::as.formula(sprintf("~ x + %s(graph = sim$graph)", term))
    fit <- tobs(f, data = sim$data, family = jsdm(), y = sim$y,
                species = sim$species, method = "nested_laplace",
                control = list(verbose = FALSE, progress = FALSE))
    expect_equal(fit$method, "nested_laplace")
    expect_gt(cor(fit$spatial_field, sim$field), 0.70)
    # Single fit vs the seed's REALIZED slope mean (#155), not the nominal 0.9:
    # this was a single-seed comparison to the nominal population constant, the
    # actual defect pattern. Budget = 3 sd of the icar-route deviation over 16
    # seeds of the sibling fixture above (sd 0.052, max 0.130), reused here as a
    # proxy since the field TYPE (icar/car_proper/bym2) does not change what
    # identifies the community slope against the shared field -- all three
    # trade the same sum-to-zero constraint. Seed 2 itself sits at -0.071
    # (within budget); the multi-seed test above documents this fixture also
    # carries a real ~-0.076 population-level slope bias, so the budget is set
    # with margin above that, not just this one draw.
    expect_lt(abs(unname(fit$means[2L]) - real_slope), 0.16)
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
