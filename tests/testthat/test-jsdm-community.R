# jsdm() as the community GLMM on an observed presence/absence response -- the
# spOccupancy lfJSDM / sfJSDM model class:
#
#   y_{s,i}         ~ Bernoulli(psi_{s,i})            (observed, no detection)
#   logit psi_{s,i} = X_i (mu + b_s) [+ f_{u(i)}] [+ sum_q lambda_{s,q} zeta_{q,i}]
#   b_s ~ N(0, Sigma)                                 (community covariance)
#
# This is exactly the ms_count() community model with a logit link, so jsdm()
# shares its binder, community Laplace-EM, latent-structure driver
# (R/community_latent.R), NUTS target and S3 methods. latent() factors give
# lfJSDM; a shared areal field alongside them gives sfJSDM.
#
# Thresholds from a measured run (dev_notes/probe_jsdm_community.R): lfJSDM
# residual correlation 0.87-0.93 at N=300/S=16; sfJSDM field 0.98, residual
# 0.82-0.92 at 256 sites / S=20. Presence/absence observed directly carries more
# information per (site, species) than a detection history, so these sit above
# the ms_occu factor thresholds.

.jsdmc_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

.jsdmc_sim <- function(N = 300L, S = 16L, Q = 2L, load_sd = 0.8, field = NULL,
                       seed = 1L) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(N))
  X <- stats::model.matrix(~ x, d)
  bs <- vapply(1:2, function(j) stats::rnorm(S, c(0.2, 0.8)[j], c(0.4, 0.3)[j]),
               numeric(S))
  eta <- X %*% t(bs)
  cor_res <- NULL
  if (Q > 0L) {
    lam <- matrix(stats::rnorm(S * Q, 0, load_sd), S, Q)
    if (!is.null(field)) lam <- scale(lam, scale = FALSE)
    zeta <- matrix(stats::rnorm(N * Q), N, Q)
    eta  <- eta + zeta %*% t(lam)
    cor_res <- stats::cov2cor(tcrossprod(lam) + diag(1e-8, S))
  }
  if (!is.null(field)) eta <- eta + matrix(field, N, S)
  y <- matrix(stats::rbinom(N * S, 1, stats::plogis(eta)), N, S,
              dimnames = list(NULL, paste0("sp", seq_len(S))))
  list(y = y, data = d, S = S, beta = c(0.2, 0.8),
       beta_real = colMeans(bs), cor_res = cor_res)
}


# The community mean is a POPULATION constant, c(0.2, 0.8), but each seed draws
# its own 16 species around it, so the mean this fixture actually realizes sits
# an SD of beta_sd / sqrt(S) -- 0.10 on the intercept, 0.075 on the slope --
# away from that constant. Scoring a fit against the population constant
# therefore spends most of its tolerance on draw noise, and what is left is too
# coarse to see estimator bias: that is how survived, a slope inflated 1.44x
# sitting inside `tolerance = 0.35`.
#
# So the assertions below score against `beta_real`, the mean of that seed's own
# draw, which removes the draw noise and leaves a pure estimator budget. They
# also compare ABSOLUTELY. testthat's numeric `tolerance` is relative only while
# the target is larger than the tolerance itself and switches to absolute below
# it (all.equal.numeric), so one `tolerance = 0.35` meant +-0.28 on the 0.8 slope
# and +-0.35 on the 0.2 intercept -- a budget that reads uniform and is not.
#
# Budgets from a 16-seed measurement (dev_notes/probe_153_intercept.R, seeds
# 11-26). Deviation from the realized mean: lfJSDM intercept sd 0.035 /
# max 0.093, slope sd 0.052 / max 0.151; the factor-free community GLMM
# intercept max 0.074, slope max 0.059.
# `expect_community_mean()` is shared across the community families that have
# been retargeted -- tests/testthat/helper-community-mean.R.


# --- (1) registry + gates ---------------------------------------------------

test_that("jsdm() reports the community backends and gates the rest", {
  m <- tulpaObs:::.tobs_family_methods$jsdm
  expect_true(all(c("laplace", "nested_laplace", "nuts") %in% m))
  # the single-block correction routes belonged to the former shared-FE model
  expect_false(any(c("laplace_sla", "laplace_gibbs", "laplace_mi") %in% m))

  d <- .jsdmc_sim(N = 60L, S = 6L, Q = 0L, seed = 3L)
  expect_error(
    tobs(~ x, data = d$data, family = jsdm(), y = d$y,
         species = colnames(d$y), detection = ~ 1, method = "laplace"),
    "no detection process")
  # a factor-only model is the block-coordinate laplace
  expect_error(
    tobs(~ x + latent(2), data = d$data, family = jsdm(), y = d$y,
         species = colnames(d$y), method = "nested_laplace"),
    "block-coordinate|laplace")
  # nested_laplace needs a field
  expect_error(
    tobs(~ x, data = d$data, family = jsdm(), y = d$y,
         species = colnames(d$y), method = "nested_laplace"),
    "needs a shared areal field")
  # y must be binary
  bad <- d$y; bad[1, 1] <- 3
  expect_error(
    tobs(~ x, data = d$data, family = jsdm(), y = bad,
         species = colnames(d$y), method = "laplace"),
    "0, 1, or NA")
})


# --- (2) the community model ------------------------------------------------

test_that("jsdm() recovers community means with per-species coefficients", {
  skip_on_cran()
  d <- .jsdmc_sim(Q = 0L, seed = 1L)
  f <- tobs(~ x, data = d$data, family = jsdm(), y = d$y,
            species = colnames(d$y), method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(f, "tobs_fit")
  # Seed 1 realizes a community intercept of 0.340 against the population's
  # 0.200, so scoring this against the population constant spends more than half
  # the old 0.25 on the draw. Against the realized mean the factor-free route is
  # unbiased (-0.006 over 16 seeds) and agrees with lme4's glmer on the same
  # fixture to 0.001, which is what licenses the tighter budget here.
  expect_community_mean(f, d$beta_real, c(0.12, 0.12))
  # per-species coefficients under a community covariance (not a scalar intercept)
  expect_equal(dim(f$ms_community$coef_mu), c(16L, 2L))
  expect_equal(dim(f$ms_community$Sigma_mu), c(2L, 2L))
  expect_true(all(f$ms_community$sd_mu > 0))
  # fitted() is on the probability scale via the logit link
  fv <- fitted(f)$mu
  expect_equal(dim(fv), c(300L, 16L))
  expect_true(all(fv > 0 & fv < 1))
  expect_true(is.finite(waic(f)$waic))
  expect_equal(dim(ranef(f)), c(32L, 4L))
})


# --- (3) lfJSDM: latent factors ---------------------------------------------

# Smoke coverage of the lfJSDM path: the factor block and the S3 surface on a
# small fixture, no threshold tied to the size.
test_that("lfJSDM wires the factor block and S3", {
  d <- .jsdmc_sim(N = 60L, S = 5L, Q = 1L, seed = 11L)
  f <- tobs(~ x + latent(1), data = d$data, family = jsdm(), y = d$y,
            species = colnames(d$y), method = "laplace",
            control = list(max.outer = 2L, factor.starts = 1L,
                           verbose = FALSE, progress = FALSE))
  expect_s3_class(f, "tobs_fit")
  expect_identical(f$ms_factor$n_factors, 1L)
  expect_equal(dim(f$ms_factor$residual_cov), c(5L, 5L))
  expect_equal(dim(f$ms_factor$loadings), c(5L, 1L))
  expect_false(is.null(f$model$count_factor_offset))
  expect_true(is.finite(waic(f)$waic))
})

test_that("lfJSDM recovers residual species co-occurrence", {
  skip_if_fast()
  skip_on_cran()
  d <- .jsdmc_sim(Q = 2L, seed = 11L)
  f <- tobs(~ x + latent(2), data = d$data, family = jsdm(), y = d$y,
            species = colnames(d$y), method = "laplace",
            control = list(verbose = FALSE, progress = FALSE))
  expect_identical(f$ms_factor$n_factors, 2L)
  expect_equal(dim(f$ms_factor$residual_cov), c(16L, 16L))
  off <- upper.tri(d$cor_res)
  expect_gt(stats::cor(f$ms_factor$residual_cor[off], d$cor_res[off]), 0.8)
  # The community coefficients used to be inflated here -- slope 1.44x truth,
  # one-sided over 16 seeds -- because the factor magnitude was a joint-mode
  # estimate with nothing bounding it. With the magnitude set by the joint site
  # marginal the slope is unbiased: 0.796 vs 0.800 over seeds 11-26, 7 seeds
  # above truth and 9 below.
  #
  # Scored per coefficient against this seed's realized mean, so neither a
  # one-sided shift nor a draw offset can hide in the other. On seed 11 the fit
  # deviates 0.039 (intercept) and 0.034 (slope); the pre-fix slope deviated
  # 0.166 on the same seed, so the slope budget separates them.
  expect_community_mean(f, d$beta_real, c(0.10, 0.12))
  expect_true(is.finite(waic(f)$waic))
})

test_that("lfJSDM recovers the residual correlation over seeds", {
  skip_if_fast()
  skip_on_cran()
  rc <- numeric(4L)
  dev <- matrix(NA_real_, 4L, 2L)
  for (s in seq_len(4L)) {
    d <- .jsdmc_sim(Q = 2L, seed = 10L + s)
    f <- tobs(~ x + latent(2), data = d$data, family = jsdm(), y = d$y,
              species = colnames(d$y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
    off <- upper.tri(d$cor_res)
    rc[s] <- stats::cor(f$ms_factor$residual_cor[off], d$cor_res[off])
    dev[s, ] <- unname(f$means[1:2]) - d$beta_real
  }
  expect_gt(stats::median(rc), 0.85)
  # A one-sided coefficient shift is a property of the MEAN deviation, not of
  # any single fit, so it is asserted here rather than on the single-seed fit
  # above: averaging 4 seeds cuts the per-seed spread in half and makes the
  # budget an estimator budget. Measured over these 4 seeds the mean deviation
  # is +0.036 (intercept) and +0.017 (slope); before was fixed the slope's was
  # +0.325.
  expect_lt(abs(mean(dev[, 1])), 0.09)
  expect_lt(abs(mean(dev[, 2])), 0.09)
})


# --- (4) sfJSDM: shared field + factors -------------------------------------

test_that("sfJSDM recovers BOTH the shared field and the factors", {
  skip_if_fast()
  skip_on_cran()
  side <- 16L
  A  <- .jsdmc_grid_graph(side)
  co <- expand.grid(r = seq_len(side), c = seq_len(side))
  fl <- scale(sin(co$r / side * pi) + cos(co$c / side * pi))[, 1]
  fl <- fl - mean(fl)
  d  <- .jsdmc_sim(N = nrow(A), S = 20L, Q = 2L, load_sd = 0.7, field = fl,
                   seed = 21L)
  f  <- tobs(~ x + icar(graph = A) + latent(2), data = d$data, family = jsdm(),
             y = d$y, species = colnames(d$y), method = "nested_laplace",
             control = list(verbose = FALSE, progress = FALSE))
  expect_identical(f$method, "nested_laplace")
  expect_false(is.null(f$spatial_field))
  expect_false(is.null(f$ms_factor))
  expect_gt(stats::cor(f$spatial_field, fl), 0.9)
  off <- upper.tri(d$cor_res)
  expect_gt(stats::cor(f$ms_factor$residual_cor[off], d$cor_res[off]), 0.75)
  expect_true(is.finite(waic(f)$waic))
})


# --- (5) NUTS: the target must match the R oracle byte-for-byte -------------

test_that("jsdm NUTS log-posterior + gradient match the R oracle", {
  skip_on_cran()
  d <- .jsdmc_sim(N = 60L, S = 6L, Q = 0L, seed = 2L)
  model <- tulpaObs:::.tobs_build_ms_count(
    formula = ~ x, data = d$data, y = d$y, species = colnames(d$y),
    response = "bernoulli")
  expect_identical(model$link, "logit")
  P   <- model$process_info[[1L]]$p; S <- model$n_species
  lay <- tulpaObs:::.tobs_ms_count_nuts_layout(P, S, "bernoulli")
  # bernoulli carries no dispersion, so it shares the Poisson layout
  expect_false(lay$is_nb); expect_false(lay$is_gauss)
  pri  <- tulpaObs:::.tobs_ms_count_nuts_priors()
  Y    <- matrix(as.numeric(model$y), model$n_sites, S)
  spec <- list(X = model$X, y = Y, family = "bernoulli")
  set.seed(102)
  for (rep in 1:4) {
    theta <- stats::rnorm(lay$total, 0, 0.4)
    theta[lay$chol_beta] <- theta[lay$chol_beta] * 0.3
    r_ora <- tulpaObs:::.tobs_ms_count_nuts_logpost(
      theta, model$X, Y, lay, pri, sigma.beta = 10, sigma.logr = 1.5)
    cpp <- cpp_ms_count_nuts_joint_logpost(spec, theta, pri, 10, 1.5)
    expect_equal(cpp$lp, r_ora$lp, tolerance = 1e-8)
    expect_lt(max(abs(cpp$grad - r_ora$grad)), 1e-7)
  }
})

test_that("jsdm NUTS recovers community means and agrees with Laplace", {
  skip_if_fast()
  skip_on_cran()
  d <- .jsdmc_sim(N = 200L, S = 10L, Q = 0L, seed = 5L)
  lap <- tobs(~ x, data = d$data, family = jsdm(), y = d$y,
              species = colnames(d$y), method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  nut <- tobs(~ x, data = d$data, family = jsdm(), y = d$y,
              species = colnames(d$y), method = "nuts",
              control = list(n.iter = 600L, n.warmup = 600L, seed = 1,
                             verbose = FALSE, progress = FALSE))
  expect_identical(nut$method, "nuts")
  expect_equal(unname(nut$means[1:2]), d$beta, tolerance = 0.3)
  # NUTS agrees with the Laplace-EM mode
  expect_equal(unname(nut$means[1:2]), unname(lap$means[1:2]), tolerance = 0.2)
  expect_true(all(is.finite(nut$sds[1:2])))
})
