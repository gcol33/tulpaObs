# =============================================================================
# test-sbc-registry.R -- the SBC family registry (gcol33/tulpaObs#207)
#
# `test-sbc.R` measures the coupled `occu_cover` route. This file covers the
# OTHER half of the registry: the families whose site marginals multiply, and
# the registry contract itself.
#
# Three tiers:
#
#   STRUCTURE (always on, no fitting). Every entry supplies the callbacks the
#   driver reads, and every callback slot holds a function. A family cannot be
#   added half-registered and discovered at the first run.
#
#   CONTRACT (skip_on_cran). Per family: the callbacks compose, the pooled data
#   set carries both blocks under disjoint site labels, and refitting the
#   observed data ALONE reproduces the observed fit -- which is what says the
#   rebuilt formulas, family and method are the ones it came from. A fit
#   carrying a structured term is refused rather than scored on its
#   coefficients alone.
#
#   ACCEPTANCE (skip_if_fast + skip_on_cran). The calibration measurement: the
#   reported posterior uniform, a deliberately mis-scaled control not.
# =============================================================================

.sbc_reg_ctl <- list(verbose = FALSE, progress = FALSE)

# One fixture per registered non-coupled family, small enough that a refit costs
# a fraction of a second. Each returns the fitted model SBC then rechecks.
.SBC_REG_FIXTURES <- list(
  occu = function(N = 60L) {
    sim <- simulate_occu(N = N, J = 4L, seed = 11L)
    suppressWarnings(tobs(~ occ_cov1, data = sim$data, family = occu(),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  count = function(N = 120L) {
    sim <- simulate_count(N = N, beta = c(1, 0.5), seed = 12L)
    suppressWarnings(tobs(~ x, data = sim$data, family = count("poisson"),
                          y = sim$y, method = "laplace",
                          control = .sbc_reg_ctl))
  },
  abun = function(N = 60L) {
    sim <- simulate_abun(N = N, J = 4L, seed = 13L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data, family = abun(),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  royle_nichols = function(N = 100L) {
    sim <- simulate_royle_nichols(N = N, J = 5L, seed = 14L)
    suppressWarnings(tobs(~ x, data = sim$data, family = royle_nichols(),
                          detection = ~ x, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  occu_ttd = function(N = 100L) {
    sim <- simulate_occu_ttd(N = N, J = 4L, seed = 15L)
    suppressWarnings(tobs(~ psi_cov1, data = sim$data,
                          family = occu_ttd(surveyLength = sim$Tmax),
                          detection = ~ rate_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  fp_occu = function(N = 120L) {
    sim <- simulate_fp_occu(N = N, J = 5L, seed = 16L)
    suppressWarnings(tobs(~ occ_cov1, data = sim$data, family = fp_occu(),
                          detection = ~ occ_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  removal = function(N = 80L) {
    sim <- simulate_removal(N = N, K = 4L, seed = 17L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data, family = removal(),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  distance = function(N = 80L) {
    sim <- simulate_distance(N = N, seed = 18L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data,
                          family = distance(cutpoints = sim$cutpoints,
                                            key = "halfnorm",
                                            transect = "line"),
                          detection = ~ sigma_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  double_observer = function(N = 150L) {
    sim <- simulate_double_observer(N = N, beta_lambda = c(log(8), 0.4),
                                    beta_p1 = c(stats::qlogis(0.5), 0.2),
                                    beta_p2 = c(stats::qlogis(0.45), -0.1),
                                    seed = 19L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data,
                          family = double_observer(),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  }
)


# ---------------------------------------------------------------------------
# Structure
# ---------------------------------------------------------------------------

test_that("every registry entry supplies the callbacks the driver reads", {
  reg <- tulpaObs:::.TOBS_SBC_REGISTRY
  req <- tulpaObs:::.TOBS_SBC_REQUIRED
  expect_gt(length(reg), 0L)

  for (fam in names(reg)) {
    e <- reg[[fam]]
    info <- paste("family", fam)
    expect_true(all(req %in% names(e)), info = info)
    for (slot in req) expect_true(is.function(e[[slot]]), info = info)
    # Optional slots, each a function when present. A rank arm needs one of the
    # two log-likelihood readers; the rest fall back to the shared driver's.
    for (slot in c("spec", "data", "pool", "loglik", "loglik_many")) {
      if (!is.null(e[[slot]])) expect_true(is.function(e[[slot]]), info = info)
    }
    expect_true(!is.null(e$loglik) || !is.null(e$loglik_many), info = info)
    # The key IS the family name `tobs()` dispatches on, so an entry filed under
    # a name no fit reports would be unreachable.
    expect_true(exists(fam, envir = asNamespace("tulpaObs"),
                       mode = "function"), info = info)
  }

  expect_identical(tulpaObs:::.tobs_sbc_registered(), sort(names(reg)))
  # The fixtures below cover every registered family bar the coupled one, which
  # test-sbc.R owns. A new entry lands here without a fixture and this fails.
  expect_setequal(names(reg),
                  c("occu_cover", names(.SBC_REG_FIXTURES)))
})


test_that("the roster in the error message names what is registered", {
  expect_error(sbc(list()), "No sbc method is registered")
  fake <- structure(list(model = list()), class = "tobs_fit",
                    tobs_family = list(name = "not_a_family"))
  expect_error(sbc(fake), "not registered for family")
  expect_error(sbc(fake), "occu_cover")
})


# ---------------------------------------------------------------------------
# Contract
# ---------------------------------------------------------------------------

test_that("each registered family composes its callbacks end to end", {
  skip_on_cran()

  for (fam in names(.SBC_REG_FIXTURES)) {
    fit <- .SBC_REG_FIXTURES[[fam]]()
    expect_identical(attr(fit, "tobs_family")$name, fam)

    m <- sbc(fit, model.only = TRUE)
    expect_true(all(c("data_obs", "fit", "draw_theta", "simulate", "pool",
                      "arms", "group_ids") %in% names(m)), info = fam)
    # Every fixed effect the fit reports is scored; nothing is silently fixed.
    expect_setequal(attr(m, "quantities"), names(fit$means))
    expect_length(attr(m, "fixed"), 0L)

    # Refitting the observed data ALONE reproduces the observed fit: the
    # rebuilt formula, family, method and control are the ones it came from.
    again <- m$fit(m$data_obs)
    expect_equal(again$means, fit$means, tolerance = 1e-8, info = fam)

    # `draw_theta` is reproducible from its seed alone and moves with it, which
    # is what the driver's RNG split relies on.
    t1 <- m$draw_theta(fit, 11L)
    expect_identical(t1, m$draw_theta(fit, 11L), info = fam)
    expect_false(isTRUE(all.equal(unname(t1),
                                  unname(m$draw_theta(fit, 12L)))), info = fam)

    # Pooling keeps BOTH data sets under disjoint site labels: fitting the
    # replicate alone would be ordinary SBC under a hand-made prior, not the
    # posterior experiment, and shared labels would deny the premise
    # `tulpa::sbc()` verifies.
    rep <- m$simulate(t1, 4400L)
    pl  <- m$pool(m$data_obs, rep)
    n_o <- length(m$group_ids(m$data_obs)); n_r <- length(m$group_ids(rep))
    expect_length(unique(m$group_ids(pl)), n_o + n_r)
    expect_identical(nrow(pl$y), n_o + n_r)
    expect_identical(nrow(pl$cells), n_o + n_r)
    expect_equal(pl$y[seq_len(n_o), , drop = FALSE], m$data_obs$y, info = fam)

    pooled <- m$fit(pl)
    expect_s3_class(pooled, "tobs_fit")
    expect_true(all(is.finite(pooled$means)), info = fam)

    # The rank arm reads the family's own marginal, so a reference set and the
    # truth score through one kernel and come back finite.
    a <- m$arms(pooled, list(theta = t1))
    expect_true("log_lik" %in% names(a$posterior), info = fam)
    expect_setequal(setdiff(names(a$posterior), "log_lik"),
                    attr(m, "quantities"))
  }
})


test_that("a structured term is refused, not scored on the coefficients", {
  skip_on_cran()
  fit <- .SBC_REG_FIXTURES$occu(N = 40L)

  # The field is a latent quantity shared across sites that theta does not
  # hold; scoring the coefficients alone would report a calibration the
  # experiment did not measure.
  spatial <- fit
  spatial$spatial <- list(type = "icar")
  expect_error(sbc(spatial, model.only = TRUE), "structured term")
  expect_error(sbc(spatial, model.only = TRUE), "fresh cells")

  re <- fit
  re$re_effects <- list(g = 1)
  expect_error(sbc(re, model.only = TRUE), "structured term")
})


test_that("a visit-level observation design is refused, not rebuilt", {
  skip_on_cran()

  # The replicate comes from the family's own simulate() kernel, which builds
  # detection from the SITE-level design; a visit-level column would be scored
  # by the refit and missing from the data it sees, and the ranks would not say
  # so. Refused at build time instead.
  N <- 40L; J <- 3L
  sim <- simulate_occu(N = N, J = J, seed = 31L)
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit = rep(seq_len(J), times = N),
                     y = as.vector(t(sim$y)),
                     det_v = stats::rnorm(N * J))
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = "det_v")
  fit <- suppressWarnings(tobs(~ occ_cov1, data = sim$data, family = occu(),
                               detection = ~ det_v, y = od$y,
                               visits = od$det.covs, method = "laplace",
                               control = .sbc_reg_ctl))
  expect_true(ncol(fit$model$X_det_visit) > 0L)
  expect_error(sbc(fit, model.only = TRUE), "visit-level")
})


test_that("the replicate generator draws at the theta it is handed", {
  skip_on_cran()

  # count(): the mean response has a closed form in theta alone, so the
  # generator can be checked against it rather than against itself. The band is
  # absolute and sized on the estimator's own Monte-Carlo error: at 200
  # replicates of 120 Poisson sites the standard error of the pooled mean is
  # under 0.02, so 0.08 is four of them.
  fit <- .SBC_REG_FIXTURES$count()
  m <- sbc(fit, model.only = TRUE)
  th <- m$draw_theta(fit, 3L)
  mu <- exp(as.vector(fit$model$X_occ %*% th[seq_len(ncol(fit$model$X_occ))]))

  reps <- vapply(seq_len(200L),
                 function(i) mean(m$simulate(th, 77000L + i)$y), numeric(1))
  expect_lt(abs(mean(reps) - mean(mu)), 0.08)
})


# ---------------------------------------------------------------------------
# Acceptance
# ---------------------------------------------------------------------------

test_that("occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  fit <- .SBC_REG_FIXTURES$occu()
  # `bad.factor` is 1.5 rather than the default 1.25 because 100 simulations on
  # four coefficients do not resolve a 20% mis-scale: measured over two seeds
  # the 1.25 control landed at 1.8e-3 and 2.1e-2, so an assertion on it would
  # be reporting the seed. At 1.5 it lands at 7.7e-7 and 2.2e-6.
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.5, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- names(fit$means)
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  # The gate is the smallest uniformity p-value over the four coefficients, at
  # a threshold well below any per-quantity band level: each band holds at 0.95
  # SEPARATELY, so requiring four at once fails on a calibrated algorithm about
  # a fifth of the time. This is a multiplicity-aware regression guard.
  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)

  # A deliberately mis-scaled posterior has to fail the same read, or the band
  # is not measuring anything. Same simulations, same fits, narrower report.
  expect_lt(min(pu("narrow")), 1e-3)
})


test_that("double_observer posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220. bad.factor = 1.5, the same tuning occu() needed: 100
  # simulations on six coefficients do not resolve the default 20% mis-scale
  # cleanly. Measured (seed = 0): posterior min p_unif 0.057, narrow 1.3e-8.
  fit <- .SBC_REG_FIXTURES$double_observer()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.5, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- names(fit$means)
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})
