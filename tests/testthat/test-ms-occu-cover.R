# =============================================================================
# test-ms-occu-cover.R - community (multispecies) joint occupancy-detection +
# cover (ms_occu_cover()).
#
# Per-species coefficient RE with Gaussian community covariances across the
# psi / p / cover arms; the latent presence z integrates out in closed form
# (the occu_cover marginal) and the per-species deviations are integrated by a
# Laplace-EM with Louis/Schur community-mean SEs (read at the natural scale).
#
# These cover:
#   - family constructor + capability gates
#   - community-mean point recovery + per-species coefficient recovery
#   - community-mean 95% CI coverage over seeds
#   - the S3 method surface
#   - a beta-arm recovery smoke
#
# The family is status = "experimental", so the coverage gate is the softer
# >= 0.80 floor (per the recovery rubric in tulpaObs CLAUDE.md). The community
# MEANS are at the natural scale and unbiased; the binary-detection community
# VARIANCE component carries the documented Laplace small-cluster attenuation
# (the EM n_quad = 1 regime), so its recovery is checked loosely.
# =============================================================================


# Build the visit-level tobs_data() detection design shared across species.
.msoc_visits <- function(N, J, visit_data) {
  long <- data.frame(
    site_id = rep(seq_len(N), each = J),
    visit   = rep(seq_len(J), times = N),
    yy      = 0L,
    det_cov1 = visit_data$det_cov1,
    pos_cov1 = visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "yy", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  od$det.covs
}


test_that("ms_occu_cover() constructor returns a tobs_family", {
  f <- ms_occu_cover()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "ms_occu_cover")
  expect_equal(f$status, "experimental")
  expect_equal(f$params$positive, "beta")
  expect_equal(ms_occu_cover("lognormal")$params$positive, "lognormal")
})


test_that("ms_occu_cover() enforces its capability gates", {
  sim <- simulate_ms_occu_cover(n_species = 4, N = 30, J = 3,
                                positive = "lognormal", seed = 1)
  vis <- .msoc_visits(30, 3, sim$visit_data)

  # Missing species.
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = sim$y, y_pos = sim$y_pos, visits = vis, method = "laplace"),
    "species"
  )
  # Missing y_pos.
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = sim$y, visits = vis, species = sim$species, method = "laplace"),
    "y_pos"
  )
  # nested_laplace is not offered (community coupled-field engine is upstream-
  # pending).
  expect_error(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
         method = "nested_laplace"),
    "not available"
  )
  # A structured term on the detection arm is rejected: Stage 1 supports a
  # shared icar() field on the occupancy arm only (gcol33/tulpa#67).
  adj <- matrix(0L, 30, 30)
  for (i in seq_len(29)) adj[i, i + 1L] <- adj[i + 1L, i] <- 1L
  expect_error(
    tobs(~ occ_cov1, data = sim$data,
         family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1 + icar(graph = adj),
         positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos, visits = vis,
         species = sim$species, method = "laplace"),
    "occupancy arm only"
  )
})


test_that("ms_occu_cover() recovers community means + per-species coefs (lognormal)", {
  skip_on_cran()
  skip_if_fast()
  set.seed(21)
  sim <- simulate_ms_occu_cover(
    n_species = 14, N = 90, J = 4,
    mu_occ = c(stats::qlogis(0.45), 0.7), mu_p = c(0.2, -0.4),
    mu_pos = c(log(0.12), 0.5), sd_occ = 0.5, sd_p = 0.4, sd_pos = 0.4,
    positive = "lognormal", sigma_pos = 0.4, seed = 21)
  vis <- .msoc_visits(90, 4, sim$visit_data)

  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE))

  expect_true(isTRUE(fit$convergence$converged))

  # Community means: within ~2.5 SE of truth on every coordinate (beta arms).
  beta_means <- fit$means[seq_len(6)]
  truth <- c(sim$truth$mu_occ, sim$truth$mu_p, sim$truth$mu_pos)
  z <- abs(beta_means - truth) / fit$sds[seq_len(6)]
  expect_true(all(z < 2.5))

  # Shared dispersion recovered.
  expect_lt(abs(exp(fit$means[["log_sigma_pos"]]) - 0.4), 0.1)

  # Per-species coefficients: high correlation with simulated truth on the
  # occupancy and cover arms; detection moderate (binary-data RE attenuation).
  # Occupancy BLUPs are detection-filtered and shrunk toward the community
  # mean, so they recover less sharply than the count-informed coefficients of
  # ms_abun (correct Bayesian shrinkage, not bias).
  cm <- fit$ms_community
  expect_gt(min(diag(cor(cm$coef_occ, sim$truth$beta_occ))), 0.78)
  expect_gt(min(diag(cor(cm$coef_pos, sim$truth$beta_pos))), 0.80)
  expect_gt(min(diag(cor(cm$coef_p,   sim$truth$beta_p))),   0.50)

  # Community SDs positive and finite; the occupancy / cover arms track the
  # realized spread, the detection arm is attenuated but non-degenerate.
  expect_true(all(cm$sd_occ > 0 & cm$sd_pos > 0 & cm$sd_p > 0))
  emp_occ <- apply(sim$truth$beta_occ, 2, sd)
  emp_pos <- apply(sim$truth$beta_pos, 2, sd)
  expect_true(all(cm$sd_occ > 0.35 * emp_occ & cm$sd_occ < 1.8 * emp_occ))
  expect_true(all(cm$sd_pos > 0.35 * emp_pos & cm$sd_pos < 1.8 * emp_pos))
})


test_that("ms_occu_cover() community-mean 95% CIs cover near the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 15L
  covered <- logical(0)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_occu_cover(
      n_species = 12, N = 70, J = 4,
      mu_occ = c(stats::qlogis(0.45), 0.6), mu_p = c(0.2, -0.4),
      mu_pos = c(log(0.12), 0.4), sd_occ = 0.5, sd_p = 0.4, sd_pos = 0.4,
      positive = "lognormal", sigma_pos = 0.4, seed = 300 + s)
    vis <- .msoc_visits(70, 4, sim$visit_data)
    fit <- tryCatch(
      tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
           detection = ~ det_cov1, positive = ~ pos_cov1,
           y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    bm    <- fit$means[seq_len(6)]
    bsd   <- fit$sds[seq_len(6)]
    truth <- c(sim$truth$mu_occ, sim$truth$mu_p, sim$truth$mu_pos)
    covered <- c(covered, abs(bm - truth) < 1.96 * bsd)
  }
  # Nominal 95%; experimental floor at 0.80 with Monte-Carlo slack on 15 x 6.
  expect_gt(mean(covered), 0.80)
})


test_that("ms_occu_cover() S3 methods work", {
  skip_on_cran()
  skip_if_fast()
  set.seed(5)
  sim <- simulate_ms_occu_cover(n_species = 8, N = 45, J = 3,
                                positive = "lognormal", seed = 5)
  vis <- .msoc_visits(45, 3, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1,
              y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  expect_no_error(summary(fit))

  # coef(): per-process list (community means by arm).
  cf <- coef(fit)
  expect_true(is.list(cf))
  expect_setequal(names(cf), c("psi", "p", "pos"))
  expect_setequal(names(cf$psi), c("(Intercept)", "occ_cov1"))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  # n_species x (p_occ + p_p + p_pos) = 8 x (2 + 2 + 2)
  expect_equal(nrow(re), 8L * 6L)
  expect_true(all(c("species", "arm", "term", "estimate") %in% names(re)))
  expect_setequal(unique(re$arm), c("psi", "p", "pos"))

  fv <- fitted(fit)
  expect_equal(dim(fv$psi),   c(45L, 8L))
  expect_equal(dim(fv$p),     c(45L, 8L))
  expect_equal(dim(fv$cover), c(45L, 8L))
  expect_true(all(fv$psi > 0 & fv$psi < 1))
  expect_true(all(fv$cover > 0))

  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys$y),     c(45L, 3L, 8L))
  expect_equal(dim(ys$y_pos), c(45L, 3L, 8L))

  expect_type(nobs(fit), "integer")
})


test_that("ms_occu_cover() recovers community means (beta arm, smoke)", {
  skip_on_cran()
  skip_if_fast()
  n_seeds <- 8L
  occ_x <- p_x <- pos_x <- phi_e <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    sim <- simulate_ms_occu_cover(
      n_species = 12, N = 70, J = 4,
      mu_occ = c(stats::qlogis(0.5), 0.6), mu_p = c(0.2, -0.4),
      mu_pos = c(stats::qlogis(0.3), 0.3), sd_occ = 0.5, sd_p = 0.4,
      sd_pos = 0.3, positive = "beta", phi = 25, seed = 500 + s)
    vis <- .msoc_visits(70, 4, sim$visit_data)
    fit <- tryCatch(
      tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("beta"),
           detection = ~ det_cov1, positive = ~ pos_cov1,
           y = sim$y, y_pos = sim$y_pos, visits = vis, species = sim$species,
           method = "laplace", control = list(verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    occ_x[s] <- fit$means[["psi_occ_cov1"]]
    p_x[s]   <- fit$means[["p_det_cov1"]]
    pos_x[s] <- fit$means[["pos_pos_cov1"]]
    phi_e[s] <- exp(fit$means[["log_phi"]])
  }
  expect_lt(abs(mean(occ_x, na.rm = TRUE) - 0.6),  0.25)
  expect_lt(abs(mean(p_x,   na.rm = TRUE) - (-0.4)), 0.25)
  expect_lt(abs(mean(pos_x, na.rm = TRUE) - 0.3),  0.25)
  expect_gt(mean(phi_e, na.rm = TRUE), 8)   # well above the boundary
})
