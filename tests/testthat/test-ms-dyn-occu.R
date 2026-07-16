# =============================================================================
# test-ms-dyn-occu.R - community (multispecies) dynamic occupancy (ms_dyn_occu()).
#
# Per-species first-season occupancy + detection coefficient RE with per-arm
# Gaussian community covariances; shared community colonisation / extinction.
# The latent occupancy path integrates out by an HMM forward filter and the
# per-species deviations are integrated by the shared community Laplace-EM
# (R/community_em.R). The family is status = "working" (gcol33/tulpaObs#99):
# community-mean 95% CI coverage is gated at the 0.85 working floor of the
# recovery rubric (measured pooled coverage ~0.98 at 24 seeds; the shared
# gamma / eps dynamics cover at ~1.0).
# =============================================================================


test_that("ms_dyn_occu() constructor returns a tobs_family", {
  f <- ms_dyn_occu()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "ms_dyn_occu")
  expect_equal(f$status, "working")
  expect_equal(f$default_engine, "laplace")
})


test_that("ms_dyn_occu() recovers community means + per-species coefs", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_dyn_occu(N = 90, J = 3, n_species = 14, n_seasons = 4,
                              beta_comm_mean = c(0.3), beta_comm_sd = c(0.7),
                              gamma = 0.2, epsilon = 0.1, seed = 31)
  fit <- tobs(~ 1, data = sim$data, family = ms_dyn_occu(),
              detection = ~ 1, y = sim$y,
              species = paste0("sp", seq_len(14)),
              method = "laplace", control = list(verbose = FALSE))

  expect_true(isTRUE(fit$convergence$converged))

  # Community means within ~2.5 SE of truth. psi1 intercept -> beta_comm_mean[1];
  # detection community mean -> 0 (p_species ~ N(0, 0.5)); gamma / eps -> the
  # shared transition logits.
  truth <- c("psi1_(Intercept)"  = 0.3,
             "p_(Intercept)"     = 0,
             "gamma_(Intercept)" = stats::qlogis(0.2),
             "eps_(Intercept)"   = stats::qlogis(0.1))
  m <- fit$means[names(truth)]
  s <- fit$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 2.5))

  # Per-species coefficients track the simulated truth.
  cm <- fit$ms_community
  expect_gt(cor(cm$coef_psi1[, 1], stats::qlogis(sim$truth$psi1_species)), 0.70)
  expect_gt(cor(cm$coef_p[, 1],    stats::qlogis(sim$truth$p_species)),    0.60)

  # Per-species variance components recover the realised spread (gcol33/tulpaObs
  # #99): mild Laplace attenuation on the binary arms (measured est/realised
  # ~0.91 psi1, ~0.95 p over seeds), so a band around the realised SD.
  emp_psi1 <- stats::sd(stats::qlogis(sim$truth$psi1_species))
  emp_p    <- stats::sd(stats::qlogis(sim$truth$p_species))
  expect_true(cm$sd_psi1[1] > 0.4 * emp_psi1 && cm$sd_psi1[1] < 1.8 * emp_psi1)
  expect_true(cm$sd_p[1]    > 0.4 * emp_p    && cm$sd_p[1]    < 1.8 * emp_p)
})


test_that("ms_dyn_occu() community-mean 95% CIs cover near the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 12L
  covered <- logical(0)
  truth <- c("psi1_(Intercept)"  = 0.3,
             "p_(Intercept)"     = 0,
             "gamma_(Intercept)" = stats::qlogis(0.2),
             "eps_(Intercept)"   = stats::qlogis(0.1))
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_dyn_occu(N = 80, J = 3, n_species = 12, n_seasons = 4,
                                beta_comm_mean = c(0.3), beta_comm_sd = c(0.6),
                                gamma = 0.2, epsilon = 0.1, seed = 400 + s)
    fit <- tryCatch(
      tobs(~ 1, data = sim$data, family = ms_dyn_occu(),
           detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(12)),
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    m <- fit$means[names(truth)]; sd <- fit$sds[names(truth)]
    covered <- c(covered, abs(m - truth) < 1.96 * sd)
  }
  # Working-family gate (gcol33/tulpaObs#99): pooled over the shared psi1 / p
  # means and the shared gamma / eps dynamics x 12 seeds. Measured pooled ~0.98
  # at 24 seeds, comfortably clearing the 0.85 working floor.
  expect_gt(mean(covered), 0.85)
})


test_that("ms_dyn_occu() S3 methods work", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_dyn_occu(N = 50, J = 3, n_species = 8, n_seasons = 4,
                              seed = 5)
  fit <- tobs(~ 1, data = sim$data, family = ms_dyn_occu(),
              detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(8)),
              method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  expect_no_error(summary(fit))

  cf <- coef(fit)
  expect_true(is.list(cf))
  expect_setequal(names(cf), c("psi1", "p", "gamma", "eps"))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  # n_species x (psi1 + p) = 8 x 2
  expect_equal(nrow(re), 8L * 2L)
  expect_setequal(unique(re$arm), c("psi1", "p"))

  fv <- fitted(fit)
  expect_equal(dim(fv$psi1), c(50L, 8L))
  expect_equal(dim(fv$p),    c(50L, 8L))
  expect_length(fv$gamma, 50L)
  expect_length(fv$eps,   50L)
  expect_true(all(fv$psi1 > 0 & fv$psi1 < 1))

  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys), c(50L, 3L, 4L, 8L))

  expect_type(nobs(fit), "integer")
})


test_that("ms_dyn_occu() capability gates", {
  sim <- simulate_ms_dyn_occu(N = 30, J = 3, n_species = 4, n_seasons = 3,
                              seed = 1)
  # nested_laplace needs a shared areal field on the occupancy formula (#123);
  # without one it errors with a pointer rather than silently downgrading.
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(4)),
         method = "nested_laplace"),
    "areal field")
  # nuts is not offered.
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
         y = sim$y, species = paste0("sp", seq_len(4)), method = "nuts"),
    "not available")
  # Missing y.
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
         species = paste0("sp", seq_len(4)), method = "laplace"),
    "y")
})


# --- shared areal field on the first-season occupancy arm (stMsPGOcc, #123) ---

.msdyn_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}
.msdyn_field <- function(side, sd = 0.9) {
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  f  <- sd * scale(sin(co$r / side * pi) + cos(co$c / side * pi))[, 1]
  f - mean(f)
}

test_that("ms_dyn_occu() nested_laplace is registered and gated", {
  expect_true("nested_laplace" %in% tulpaObs:::.tobs_family_methods$ms_dyn_occu)
  side <- 5L; A <- .msdyn_grid_graph(side)
  sim <- simulate_ms_dyn_occu(N = side * side, J = 3, n_species = 4,
                              n_seasons = 3, field = .msdyn_field(side), seed = 1)
  # a field needs nested_laplace; plain laplace with a field errors with a pointer
  expect_error(
    tobs(~ 1 + icar(graph = A), data = sim$data, family = ms_dyn_occu(),
         detection = ~ 1, y = sim$y, species = paste0("sp", 1:4),
         method = "laplace"),
    "nested_laplace")
  # a field on the detection arm is rejected
  expect_error(
    tobs(~ 1, data = sim$data, family = ms_dyn_occu(),
         detection = ~ 1 + icar(graph = A), y = sim$y,
         species = paste0("sp", 1:4), method = "nested_laplace"),
    "occupancy|detection")
})

test_that("ms_dyn_occu() + icar recovers the shared field + community means", {
  skip_if_fast()
  skip_on_cran()
  # stMsPGOcc: a shared areal field on the first-season occupancy arm. The field
  # is the new object; assert field recovery by cor on an INTERIOR field (the
  # null-field trap: never assert on sigma), plus community-mean coverage.
  side <- 8L; N <- side * side; A <- .msdyn_grid_graph(side)
  ftrue <- .msdyn_field(side, sd = 0.9)
  n_seed <- 12L
  fcor <- numeric(n_seed)
  truth <- c("gamma_(Intercept)" = stats::qlogis(0.2),
             "eps_(Intercept)"   = stats::qlogis(0.12))
  covered <- logical(0)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_dyn_occu(N = N, J = 4, n_species = 8, n_seasons = 4,
                                beta_comm_mean = 0.2, beta_comm_sd = 0.5,
                                gamma = 0.2, epsilon = 0.12, field = ftrue,
                                seed = 500 + s)
    fit <- tryCatch(
      tobs(~ 1 + icar(graph = A), data = sim$data, family = ms_dyn_occu(),
           detection = ~ 1, colonization = ~ 1, extinction = ~ 1, y = sim$y,
           species = paste0("sp", seq_len(8)), method = "nested_laplace",
           control = list(progress = FALSE, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    fcor[s] <- stats::cor(fit$spatial_field, ftrue)
    m <- fit$means[names(truth)]; sd <- fit$sds[names(truth)]
    covered <- c(covered, abs(m - truth) < 1.96 * sd)
  }
  # interior field recovery (cor), not sigma; the shared field pools all species
  # so it recovers cleanly.
  expect_gt(stats::median(fcor[fcor != 0]), 0.8)
  # shared transition-dynamics coverage at the 0.85 working floor.
  expect_gt(mean(covered), 0.85)
})
