# Parameter-recovery tests for the prior-regularised occupancy Laplace path.
#
# Background (2026-05-15). The unpenalised EM-Laplace MAP sits on the psi-p
# identifiability ridge at small J: small revisit counts mean "site is
# unoccupied" and "site is occupied but detection failed" are nearly
# indistinguishable from data alone. NUTS escapes via priors; the unpenalised
# Laplace doesn't. The default prior (see `occu_priors()`, attached per M-step
# block as a tulpa `beta_prior`) adds a
# weakly-informative Normal prior on the detection intercept (sd = 1.5),
# detection slopes (sd = 2.5), psi intercept (sd = 2), and psi slopes
# (sd = 5).
#
# These tests deliberately use small N and J so the ridge is sharp; the
# point of the prior is to break it without imparting a bias of its own at
# moderate-information settings.

test_that("penalised default recovers beta_p across multiple seeds", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 15L
  truth <- list(
    beta_psi_int   = 0.5,
    beta_psi_slope = 0.5,
    beta_p_int     = 0.0,
    beta_p_slope   = 0.8
  )
  N <- 200L; J <- 4L

  hat_pen <- vector("list", n_seeds)
  for (s in seq_len(n_seeds)) {
    set.seed(2000L + s)
    x_occ <- rnorm(N)
    x_det <- rnorm(N)
    psi   <- plogis(truth$beta_psi_int + truth$beta_psi_slope * x_occ)
    p_tru <- plogis(truth$beta_p_int   + truth$beta_p_slope   * x_det)
    z     <- rbinom(N, 1L, psi)
    y     <- matrix(0L, N, J)
    for (i in seq_len(N)) {
      if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p_tru[i])
    }
    dat <- data.frame(x_occ = x_occ, x_det = x_det)

    fit <- tryCatch(
      tobs(formula = ~ x_occ, data = dat, family = occu(),
           detection = ~ x_det, y = y, method = "laplace",
           control = list(verbose = FALSE)),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    hat_pen[[s]] <- list(
      beta_p_int   = unname(fit$means["p_(Intercept)"]),
      beta_p_slope = unname(fit$means["p_x_det"]),
      se_p_int     = unname(fit$sds[["p_(Intercept)"]]),
      se_p_slope   = unname(fit$sds[["p_x_det"]])
    )
  }
  hat_pen <- Filter(Negate(is.null), hat_pen)
  expect_true(length(hat_pen) >= floor(0.8 * n_seeds))

  est_int   <- vapply(hat_pen, `[[`, numeric(1), "beta_p_int")
  est_slope <- vapply(hat_pen, `[[`, numeric(1), "beta_p_slope")
  # `|estimate - truth| < 0.2` for the detection slope across seeds.
  expect_lt(median(abs(est_slope - truth$beta_p_slope)), 0.2)
  # Loose tolerance on the detection intercept: 0.4 reflects residual
  # finite-sample slack even with the prior on board.
  expect_lt(median(abs(est_int - truth$beta_p_int)), 0.4)
})


test_that("disabling the prior makes the detection slope bias worse", {
  skip_on_cran()
  skip_if_fast()

  # Probe the ridge using the detection-slope bias. At N = 200, J = 4 the
  # unpenalised MAP drifts along the psi-p ridge and *under*-shrinks the
  # detection slope toward zero (smoking-gun diagnostic 2026-05-15: mean
  # `p_x_det` was 0.537 vs truth 0.8 across 10 seeds, vs 0.815 penalised).
  n_seeds <- 15L
  truth_slope <- 0.8
  N <- 200L; J <- 4L

  est_pen   <- numeric(0)
  est_unpen <- numeric(0)

  for (s in seq_len(n_seeds)) {
    set.seed(3000L + s)
    x_occ <- rnorm(N)
    x_det <- rnorm(N)
    psi   <- plogis(0.5 + 0.5 * x_occ)
    p_tru <- plogis(0.0 + truth_slope * x_det)
    z     <- rbinom(N, 1L, psi)
    y     <- matrix(0L, N, J)
    for (i in seq_len(N)) {
      if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p_tru[i])
    }
    dat <- data.frame(x_occ = x_occ, x_det = x_det)

    fit_pen <- tryCatch(
      tobs(formula = ~ x_occ, data = dat, family = occu(),
           detection = ~ x_det, y = y, method = "laplace",
           control = list(verbose = FALSE)),
      error = function(e) NULL
    )
    fit_unp <- tryCatch(
      tobs(formula = ~ x_occ, data = dat, family = occu(),
           detection = ~ x_det, y = y, method = "laplace",
           priors = FALSE,
           control = list(verbose = FALSE)),
      error = function(e) NULL
    )
    if (!is.null(fit_pen)) {
      est_pen <- c(est_pen, unname(fit_pen$means[["p_x_det"]]))
    }
    if (!is.null(fit_unp)) {
      est_unpen <- c(est_unpen, unname(fit_unp$means[["p_x_det"]]))
    }
  }

  expect_true(length(est_pen)   >= floor(0.8 * n_seeds))
  expect_true(length(est_unpen) >= floor(0.8 * n_seeds))

  bias_pen   <- mean(est_pen)   - truth_slope
  bias_unpen <- mean(est_unpen) - truth_slope
  # Penalised mean bias must be materially smaller in magnitude than the
  # unpenalised bias — this is the sanity check that the prior is doing
  # something at the ridge.
  expect_lt(abs(bias_pen), abs(bias_unpen))
  expect_lt(abs(bias_pen), 0.1)
})


test_that("95% Wald CI covers the detection slope at the nominal rate", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 20L
  truth_slope <- 0.8
  N <- 200L; J <- 4L

  covered <- logical(0)
  for (s in seq_len(n_seeds)) {
    set.seed(4000L + s)
    x_occ <- rnorm(N)
    x_det <- rnorm(N)
    psi   <- plogis(0.5 + 0.5 * x_occ)
    p_tru <- plogis(0.0 + truth_slope * x_det)
    z     <- rbinom(N, 1L, psi)
    y     <- matrix(0L, N, J)
    for (i in seq_len(N)) {
      if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p_tru[i])
    }
    dat <- data.frame(x_occ = x_occ, x_det = x_det)

    fit <- tryCatch(
      tobs(formula = ~ x_occ, data = dat, family = occu(),
           detection = ~ x_det, y = y, method = "laplace",
           control = list(verbose = FALSE)),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    est <- unname(fit$means[["p_x_det"]])
    se  <- unname(fit$sds[["p_x_det"]])
    if (!is.finite(est) || !is.finite(se) || se <= 0) next
    lo <- est - 1.96 * se
    hi <- est + 1.96 * se
    covered <- c(covered, lo <= truth_slope && truth_slope <= hi)
  }
  expect_true(length(covered) >= floor(0.8 * n_seeds))

  # Nominal 95% CI; allow 85% floor as per the recovery-tests-not-smoke-tests
  # rubric in CLAUDE.md (small-N + EM-conditional Wald CI is conservative
  # but not strictly nominal).
  expect_gte(mean(covered), 0.85)
})


test_that("occu_priors() validates its arguments", {
  expect_s3_class(occu_priors(), "occu_priors")

  expect_error(
    occu_priors(p_intercept = list(mean = 0, sd = -1)),
    "must be positive"
  )
  expect_error(
    occu_priors(p_intercept = list(mean = NA, sd = 1)),
    "must be finite"
  )
  expect_error(
    occu_priors(p_intercept = list(sd = 1)),
    "must be a list with `mean` and `sd`"
  )

  # Inf sd is allowed (no penalty on that bucket).
  pr <- occu_priors(p_slope = list(mean = 0, sd = Inf))
  expect_true(is.infinite(pr$p_slope$sd))
})


test_that("priors = FALSE disables the penalty entirely", {
  skip_on_cran()
  skip_if_fast()

  set.seed(123)
  N <- 200L; J <- 4L
  x_occ <- rnorm(N); x_det <- rnorm(N)
  psi   <- plogis(0.5 + 0.5 * x_occ)
  p_tru <- plogis(0.0 + 0.8 * x_det)
  z     <- rbinom(N, 1L, psi)
  y     <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    if (z[i] == 1L) y[i, ] <- rbinom(J, 1L, p_tru[i])
  }
  dat <- data.frame(x_occ = x_occ, x_det = x_det)

  fit_off <- tobs(formula = ~ x_occ, data = dat, family = occu(),
                  detection = ~ x_det, y = y, method = "laplace",
                  priors = FALSE,
                  control = list(verbose = FALSE))
  expect_null(fit_off$priors)
  # priors=FALSE routes through tulpa::tulpa_em_laplace -> tulpa_laplace
  expect_s3_class(fit_off, "tobs_fit")
})
