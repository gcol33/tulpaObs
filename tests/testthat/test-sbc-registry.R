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
  },
  dyn_occu = function(N = 80L) {
    sim <- simulate_dyn_occu(N = N, J = 4L, n_seasons = 5L,
                             beta_occ = c(0.2, 0.6), beta_det = c(0.4),
                             gamma = 0.25, epsilon = 0.15, seed = 21L)
    suppressWarnings(tobs(~ x, data = sim$data, family = dyn_occu(),
                          detection = ~ 1, colonization = ~ 1,
                          extinction = ~ 1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  int_occu = function(N = 60L) {
    set.seed(22L)
    x_cov <- rnorm(N); det_cov <- rnorm(N)
    z <- rbinom(N, 1, plogis(0.2 + 0.7 * x_cov))
    mk <- function(J, p0) {
      p <- plogis(p0 + 0.4 * det_cov)
      y <- matrix(0L, N, J)
      for (i in seq_len(N)) if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p[i])
      y
    }
    dat <- data.frame(occ_cov = x_cov, det_cov = det_cov)
    yy <- list(src1 = mk(4L, -0.2), src2 = mk(3L, -0.5))
    suppressWarnings(tobs(~ occ_cov, data = dat, family = int_occu(),
                          detection = ~ det_cov, y = yy,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  dyn_abun = function(N = 80L) {
    sim <- simulate_dyn_abun(N = N, T = 4L, J = 3L,
                             beta_lambda = c(log(5), 0.3), p = 0.5,
                             omega = 0.6, gamma = 1.0, seed = 23L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data, family = dyn_abun(),
                          detection = ~ 1, omega = ~ 1, gamma = ~ 1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  gdistremoval = function(N = 150L) {
    cutp <- c(0, 25, 50, 75, 100)
    sim <- simulate_gdistremoval(N = N, cutpoints = cutp, n_periods = 4L,
                                 beta_lambda = c(log(30), 0.3),
                                 beta_sigma = c(log(18), 0.1),
                                 beta_r = c(stats::qlogis(0.4), -0.2), seed = 11L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data, y = sim$y,
                          y_rem = sim$y_rem,
                          family = gdistremoval(cutpoints = cutp),
                          detection = ~ det_cov1, removal = ~ rem_cov1,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  occu_categorical = function(N = 300L) {
    sim <- simulate_occu_categorical(N = N, seed = 11L)
    suppressWarnings(tobs(~ x, data = sim$data, family = occu_categorical(),
                          y = sim$y, method = "laplace", control = .sbc_reg_ctl))
  },
  distsamp_open = function(N = 100L) {
    cutp <- c(0, 10, 20, 30, 40)
    sim <- simulate_distsamp_open(N = N, cutpoints = cutp, n_seasons = 4L,
                                  beta_lambda = c(log(15), 0.3),
                                  beta_sigma = c(log(15), 0.1),
                                  omega = 0.7, gamma = 2.5, seed = 7L)
    suppressWarnings(tobs(~ abund_cov1, data = sim$data,
                          family = distsamp_open(cutpoints = cutp),
                          detection = ~ det_cov1, y = sim$y,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  occu_multi = function(N = 150L) {
    sim <- simulate_occu_multi(S = 2L, N = N, J = 4L, seed = 5L)
    suppressWarnings(tobs(~ scov1, data = sim$data, family = occu_multi(),
                          detection = ~ 1, y = sim$y, species = sim$species,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  dyn_int_occu = function(N = 60L) {
    sim <- simulate_dyn_int_occu(N = N, T_seasons = 4L, S = 2L, J = 3L,
                                 psi1 = 0.5, gamma = 0.3, eps = 0.2,
                                 p = c(0.4, 0.6), seed = 25L)
    suppressWarnings(tobs(~ 1, data = sim$data, family = dyn_int_occu(),
                          detection = ~ 1, colonization = ~ 1,
                          extinction = ~ 1, y = sim$y, sources = sim$sources,
                          method = "laplace", control = .sbc_reg_ctl))
  },
  t_occu = function(N = 100L) {
    sim <- simulate_t_occu(N = N, T_seasons = 6L, J = 3L, beta_occ = c(0.2, 0.6),
                           p = 0.4, rho = 0.6, sigma = 0.7, seed = 31L)
    suppressWarnings(tobs(~ x, data = sim$data, family = t_occu(),
                          detection = ~ 1, y = sim$y,
                          method = "pg_gibbs", control = .sbc_reg_ctl))
  },
  # ms_occu: not registered -- see gcol33/tulpaObs#226 (R/sbc.R section 6j).
  cover = function(N = 200L) {
    sim <- simulate_cover(N = N, beta_occ = c(-0.5, 0.8), beta_pos = c(-1.0, 0.3),
                          sigma_pos = 0.4, response = "lognormal", seed = 51L)
    suppressWarnings(tobs(~ x, data = sim$data, family = cover("lognormal"),
                          y = sim$y, method = "laplace", control = .sbc_reg_ctl))
  },
  # ms_int_occu: not registered -- see gcol33/tulpaObs#226 (R/sbc.R section 6l).
  # Confirmed (not just suspected) to share ms_occu's failure mode.
  occu_multiscale_cover = function(n_cells = 40L) {
    sim <- simulate_occu_multiscale_cover(
      n_cells = n_cells, plots_per_cell = 4L, visits_per_plot = 2L,
      beta_psi = c(0.4, 0.6), beta_theta = c(0.2, 0.5),
      beta_p = c(0.0, 0.5), beta_pos = c(log(0.10), -0.4),
      positive = "lognormal", phi = 0.35, seed = 61L)
    suppressWarnings(tobs(~ x_cell + icar(graph = sim$adj, group_var = "cell"),
                          data = sim$data,
                          family = occu_multiscale_cover(response = "lognormal"),
                          detection = ~ x_pdet, availability = ~ x_plot,
                          positive = ~ x_cov, y = sim$y, y_pos = sim$y_pos,
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

# occu_categorical carries no fit$means -- two independent Laplace-Gaussian
# blocks (presence, class), not tulpa's usual joint MVN summary -- so this
# flattens it with the SAME arm-prefix scheme .tobs_sbc_draws_occu_categorical()
# uses, letting the generic checks below read one name per reported
# coefficient regardless of family.
.sbc_reg_means <- function(fit) {
  if (!is.null(fit$means)) return(fit$means)
  if (inherits(fit, "occu_categorical_fit")) {
    occ <- fit$beta_occ
    names(occ) <- paste0("occ_", names(occ))
    cls <- as.numeric(fit$beta_class)
    names(cls) <- as.vector(outer(rownames(fit$beta_class),
                                  colnames(fit$beta_class),
                                  function(r, cl) paste0("class_", cl, "_", r)))
    return(c(occ, cls))
  }
  if (inherits(fit, "cover_fit")) {
    occ <- fit$beta_occ; names(occ) <- paste0("occ_", names(occ))
    pos <- fit$beta_pos; names(pos) <- paste0("pos_", names(pos))
    return(c(occ, pos))
  }
  stop("no means accessor for this fit class (", paste(class(fit), collapse = "/"),
       ")", call. = FALSE)
}

test_that("each registered family composes its callbacks end to end", {
  skip_on_cran()

  for (fam in names(.SBC_REG_FIXTURES)) {
    fit <- .SBC_REG_FIXTURES[[fam]]()
    expect_identical(attr(fit, "tobs_family")$name, fam)

    m <- sbc(fit, model.only = TRUE)
    expect_true(all(c("data_obs", "fit", "draw_theta", "simulate", "pool",
                      "arms", "group_ids") %in% names(m)), info = fam)
    # Every fixed effect the fit reports is scored; nothing is silently fixed.
    expect_setequal(attr(m, "quantities"), names(.sbc_reg_means(fit)))
    expect_length(attr(m, "fixed"), 0L)

    # Refitting the observed data ALONE reproduces the observed fit: the
    # rebuilt formula, family, method and control are the ones it came from.
    again <- m$fit(m$data_obs)
    expect_equal(.sbc_reg_means(again), .sbc_reg_means(fit), tolerance = 1e-8,
                info = fam)

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
    expect_identical(nrow(pl$cells), n_o + n_r)
    # `y` is a plain 2D matrix for most families, a 3D [site x visit x season]
    # array for the multi-season group (pooled on the site axis alone), a list
    # of one matrix (or, for the product shape, one 3D array) per source for
    # the multi-source group (each pooled on the site axis independently), or
    # a plain length-N vector for a family with one observation per unit and
    # no visit/season axis (occu_categorical).
    slice_site <- function(z) if (length(dim(z)) == 3L)
      z[seq_len(n_o), , , drop = FALSE] else z[seq_len(n_o), , drop = FALSE]
    if (is.list(pl$y) && !is.data.frame(pl$y)) {
      expect_identical(vapply(pl$y, nrow, integer(1)),
                       stats::setNames(rep(n_o + n_r, length(pl$y)), names(pl$y)),
                       info = fam)
      obs_slice <- lapply(pl$y, slice_site)
    } else if (is.null(dim(pl$y))) {
      expect_identical(length(pl$y), n_o + n_r)
      obs_slice <- pl$y[seq_len(n_o)]
    } else {
      expect_identical(nrow(pl$y), n_o + n_r)
      obs_slice <- slice_site(pl$y)
    }
    # Values only -- a family's own y may carry dimnames (site/bin/season
    # labels) the freshly-built pooled array never does; that carries no
    # information this check is about.
    strip_dn <- function(z) if (is.list(z)) lapply(z, unname) else unname(z)
    expect_equal(strip_dn(obs_slice), strip_dn(m$data_obs$y), info = fam)

    pooled <- m$fit(pl)
    expect_s3_class(pooled, "tobs_fit")
    expect_true(all(is.finite(.sbc_reg_means(pooled))), info = fam)

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


test_that("dyn_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (multi-season group). bad.factor = 1.75 -- the extra
  # tick past double_observer's 1.5 that epsilon's weak identification at 5
  # seasons needs to clear 1e-3 cleanly. Measured (seed = 0): posterior min
  # p_unif 0.090 (epsilon), narrow max 1.5e-4 (epsilon).
  fit <- .SBC_REG_FIXTURES$dyn_occu()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

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


test_that("int_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#225 (via #220). Root-caused to two compounding bugs, not
  # an SBC adapter issue: (1) int_occu()'s per-source detection design
  # matrix lost its column names when padded to full-site width
  # (R/occu.R), so the autoscale unscale step couldn't find the intercept
  # column and left the detection intercept in standardized-covariate units
  # while correctly unscaling the slope; (2) model_type == "integrated"
  # never got the exact-marginal Newton debiasing single-season occu() and
  # dyn_occu() already have (R/int_occu_marginal.R). Verified: int_occu()
  # with one source now matches occu() on the same data to full optim
  # precision (all four coefficients and both SDs bit-identical). bad.factor
  # = 1.75, the same tuning dyn_occu() needed. Measured (seed = 0): posterior
  # min p_unif 0.171, narrow max 1.0e-5.
  fit <- .SBC_REG_FIXTURES$int_occu()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- names(fit$means)
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok_int_occu <- pu("posterior")
  expect_length(ok_int_occu, length(qs))
  expect_gt(min(ok_int_occu), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})


test_that("dyn_abun posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (multi-season group). dyn_abun shares dyn_occu's 3D
  # response and site-axis pooling but has its own working simulate(), so the
  # replicate is the shared simple-family route, not a bespoke forward
  # simulator. Measured (seed = 0): posterior min p_unif 0.087, narrow max
  # 1.6e-7.
  fit <- .SBC_REG_FIXTURES$dyn_abun()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- names(fit$means)
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok_dyn_abun <- pu("posterior")
  expect_length(ok_dyn_abun, length(qs))
  expect_gt(min(ok_dyn_abun), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})


test_that("gdistremoval posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (multi-response group). bad.factor = 2.0 -- the
  # single time-step, two-response-matrix shape needed a bigger tick than
  # the other registered families to clear 1e-3 on both arms cleanly.
  # Measured (seed = 0): posterior min p_unif 0.449, narrow max 2.0e-9.
  fit <- .SBC_REG_FIXTURES$gdistremoval()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 2.0, seed = 0L)

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


test_that("distsamp_open posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (multi-season group). Constant-dynamics, Poisson only
  # for v1 -- shares dyn_abun's 3D response/pooling and fit$means/fit$draws
  # shape. bad.factor tuned below.
  fit <- .SBC_REG_FIXTURES$distsamp_open()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

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


test_that("occu_multi posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (multi-response group). Same list-of-matrices response
  # shape int_occu() pools, but simulate() is custom (joint multi-species
  # state, not independent per-source arms). Measured (seed = 0): posterior
  # min p_unif 0.019, narrow max 4.3e-6.
  fit <- .SBC_REG_FIXTURES$occu_multi()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

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


test_that("dyn_int_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (the multi-season x multi-source product shape,
  # section 6h). A named list of S per-source 3D arrays, pooled with
  # `.tobs_sbc_pool_named_3d`; simulate() wraps the family's own
  # `.tobs_simulate_dyn_int_occu()` handler. Measured (seed = 0): posterior
  # min p_unif 0.033, narrow max 2.5e-5.
  fit <- .SBC_REG_FIXTURES$dyn_int_occu()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

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


test_that("t_occu posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (section 6i). A pg_gibbs family whose fit$draws is
  # already the real pooled posterior sample; loglik_many is a Laplace
  # approximation to the AR1 year effect's marginal (FD-validated gradient
  # and Hessian, brute-force-cross-checked against a dense grid at T = 2).
  # `qs` only scores the reported psi/p/AR1-hyperparameter coefficients
  # (matching every other family's test); log_lik's own presence/finiteness
  # is covered by the generic cross-family CONTRACT test.
  # Measured (seed = 0): posterior min p_unif 0.081, narrow max 1.2e-4.
  fit <- .SBC_REG_FIXTURES$t_occu()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

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


test_that("cover posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (multiarm-S3 group). Same two-independent-block
  # shape as occu_categorical (presence, positive); positive = "lognormal"
  # only for v1, dispersion held fixed (no SE anywhere in the package for
  # it). Checked at both N=200 and N=600 during development -- consistent,
  # no anomaly like the ms_occu near-miss (#226). Measured (seed = 0, N=200):
  # posterior min p_unif 0.124, narrow max 1.8e-8.
  fit <- .SBC_REG_FIXTURES$cover()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$cover$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})


test_that("occu_categorical posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (multiarm-S3 group). occu_categorical has no
  # fit$means/fit$draws -- two independent Laplace-Gaussian blocks (presence,
  # class) instead of a joint MVN draw matrix -- so the scored quantity names
  # come off the registry's own draws() rather than fit$means. Measured
  # (seed = 0): posterior min p_unif 0.057, narrow max 2.8e-5.
  fit <- .SBC_REG_FIXTURES$occu_categorical()
  res <- sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L)

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs <- colnames(tulpaObs:::.TOBS_SBC_REGISTRY$occu_categorical$draws(fit, 2L))
  pu <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok <- pu("posterior")
  expect_length(ok, length(qs))
  expect_gt(min(ok), 1e-3)
  expect_lt(min(pu("narrow")), 1e-3)
})


test_that("occu_multiscale_cover posterior SBC: correct fit uniform, mis-scaled is not", {
  skip_on_cran()
  skip_if_fast()

  # gcol33/tulpaObs#220 (multiarm-S3 group, section 6m). The standard
  # single-block fit shape (unlike cover()); the exchangeable unit is the
  # CELL, not the plot, so pooling/site track cell indices. Checked at three
  # configurations (n_cells=40 seeds 0/1, n_cells=120 seed 0) before
  # registering -- all consistent, no anomaly like the ms_occu/ms_int_occu
  # near-misses (#226). Measured (n_cells=40, seed=0): posterior min p_unif
  # 0.093 (psi_(Intercept)), narrow max 5.5e-6.
  fit <- .SBC_REG_FIXTURES$occu_multiscale_cover()
  res <- suppressMessages(sbc(fit, n.sim = 100L, n.draws = 1000L, n.ref = 200L,
                  controls = "narrow", bad.factor = 1.75, seed = 0L))

  expect_s3_class(res, "sbc")
  expect_identical(res$premises$pooling, "verified")
  expect_identical(res$premises$fresh_groups, "verified (disjoint group labels)")

  rp <- res$report
  qs2 <- names(fit$means)
  pu2 <- function(arm) {
    r <- rp[rp$arm == arm & rp$quantity %in% qs2, ]
    stats::setNames(r$p_unif, r$quantity)
  }

  ok2 <- pu2("posterior")
  expect_length(ok2, length(qs2))
  expect_gt(min(ok2), 1e-3)
  expect_lt(min(pu2("narrow")), 1e-3)
})
