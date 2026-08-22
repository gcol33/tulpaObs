# =============================================================================
# test-occu-cover-nuts-re.R - observation-arm random intercepts under
# occu_cover(method = "nuts").
#
# A `(1 | g)` on the detection or positive-cover formula becomes one shared
# non-centered ReBlock (src/nuts_re_block.h) in the occu_cover NUTS target, with
# its own SAMPLED log sigma_re. Tests:
#   - C++ FullGradFn == R oracle over the new coordinates, and the R oracle's
#     analytic gradient == a central finite difference of its own lp
#   - byte-identity: a block that is off (z = 0, log_sigma_re = 0) leaves lp and
#     the coefficient gradient bit for bit unchanged, and so does an empty block
#     list -- the no-RE sampler is the same target it was
#   - the prior on the SAMPLED coordinate is exactly N(0, sigma.logre^2): at
#     z = 0 the data term does not depend on log sigma_re, so the target reduces
#     to the prior alone (a missing or stray change-of-variables term would show
#     up here and nowhere in a gradient check)
#   - the sampled group SD covers a known truth over seeds, on both arms
#   - dispatch: what still errors, and with which limitation named
# =============================================================================


# One non-spatial occu_cover simulation carrying a per-visit `habitat` grouping
# with a random intercept on the requested arm. `arm = "p"` puts the truth on
# detection, `arm = "pos"` on the cover magnitude.
.ocnre_sim <- function(seed, arm = "p", N = 120L, J = 5L, n_g = 8L,
                       sigma_re = 0.8, beta_occ = NULL, beta_p = NULL,
                       beta_pos = NULL, sigma_pos = 0.4) {
  simulate_occu_cover(
    N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive = "lognormal", beta_occ = beta_occ, beta_p = beta_p,
    beta_pos = beta_pos, sigma_pos = sigma_pos,
    re_det_groups = if (identical(arm, "p")) n_g else NULL, sigma_re_p = sigma_re,
    re_pos_groups = if (identical(arm, "pos")) n_g else NULL,
    sigma_re_pos = sigma_re, seed = seed)
}

.ocnre_fit <- function(sim, detection = ~ det_cov1 + (1 | habitat),
                       positive = ~ pos_cov1, method = "nuts", control = list(),
                       occurrence = ~ occ_cov1) {
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  tobs(formula = occurrence, data = sim$data, family = occu_cover("lognormal"),
       detection = detection, positive = positive, y = sim$y, y_pos = y_pos,
       visits = sim$visit_data, method = method,
       control = utils::modifyList(list(verbose = FALSE, progress = FALSE),
                                   control))
}

# The bound model, reached through the front door with a two-iteration run (the
# RE designs are resolved by the dispatcher, so there is no shortcut past it).
.ocnre_model <- function(sim, detection, positive) {
  .ocnre_fit(sim, detection, positive,
             control = list(n.iter = 2L, n.warmup = 2L, n.chains = 1L,
                            seed = 1L))$model
}

# Truth BLUPs in the fit's own level order. A fit sorts its levels
# lexicographically ("hab10" lands before "hab2"), so a positional comparison
# against the simulator's numeric order silently scrambles the pairing.
.ocnre_truth_blup <- function(sim, re, arm = "p") {
  lev <- if (identical(arm, "p")) sim$truth$re_det_levels
         else sim$truth$re_pos_levels
  b   <- as.numeric(if (identical(arm, "p")) sim$truth$b_p_re
                    else sim$truth$b_pos_re)
  out <- b[match(re$levels, lev)]
  stopifnot(!anyNA(out))
  out
}

# The C++ spec + R block list for a model, exactly as the fitter assembles them.
.ocnre_spec <- function(model) {
  blocks <- tulpaObs:::.occu_cover_nuts_re_blocks(model)
  spec   <- tulpaObs:::.tobs_occu_cover_nuts_spec(model)
  if (!is.null(blocks)) {
    spec$re_blocks <- lapply(blocks, function(b)
      list(re_arm = b$arm, re_group = b$group, n_re_groups = b$n_groups,
           sigma_re_lsd = b$sigma_lsd))
  }
  pin <- model$process_info
  list(spec = spec, blocks = blocks,
       n_par = pin[[1L]]$p + pin[[2L]]$p + pin[[3L]]$p + 1L,
       n_re = if (is.null(blocks)) 0L
              else sum(vapply(blocks, function(b) b$n_groups + 1L, numeric(1))))
}


test_that("occu_cover NUTS RE: C++ target == R oracle, and both == a finite difference", {
  skip_on_cran()
  sim <- .ocnre_sim(11L, arm = "p", N = 50L, J = 4L, n_g = 4L)
  # A second per-visit grouping, so the detection and cover arms can each carry
  # their own factor (two blocks, laid out back to back).
  set.seed(111L)
  sim$visit_data$region <- factor(paste0("r", sample.int(3L, 50L * 4L, TRUE)))

  cases <- list(
    det  = list(d = ~ det_cov1 + (1 | habitat), p = ~ pos_cov1),
    pos  = list(d = ~ det_cov1,                 p = ~ pos_cov1 + (1 | region)),
    both = list(d = ~ det_cov1 + (1 | habitat), p = ~ pos_cov1 + (1 | region)))

  for (nm in names(cases)) {
    cs    <- cases[[nm]]
    model <- .ocnre_model(sim, cs$d, cs$p)
    sp    <- .ocnre_spec(model)
    expect_false(is.null(sp$blocks))

    set.seed(7)
    theta <- c(stats::rnorm(sp$n_par, 0, 0.3), stats::rnorm(sp$n_re, 0, 0.3))
    theta[sp$n_par] <- log(0.4)                       # a sane log-dispersion

    rr <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, model, 5, 5,
                                                   re = sp$blocks)
    cc <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(sp$spec, theta, 5, 5)
    expect_equal(rr$lp, cc$lp, tolerance = 1e-10, info = nm)
    expect_equal(rr$grad, as.numeric(cc$grad), tolerance = 1e-10, info = nm)
    expect_length(rr$grad, sp$n_par + sp$n_re)

    # Every coordinate, the RE block's z and log_sigma_re included.
    fd <- numeric(length(theta)); h <- 1e-5
    for (k in seq_along(theta)) {
      tp <- theta; tp[k] <- tp[k] + h
      tm <- theta; tm[k] <- tm[k] - h
      fd[k] <- (tulpaObs:::.tobs_occu_cover_nuts_logpost(tp, model, 5, 5, re = sp$blocks)$lp -
                tulpaObs:::.tobs_occu_cover_nuts_logpost(tm, model, 5, 5, re = sp$blocks)$lp) / (2 * h)
    }
    expect_equal(rr$grad, fd, tolerance = 1e-4, info = nm)
  }
})


test_that("occu_cover NUTS RE: an off block leaves the no-RE target bit for bit", {
  skip_on_cran()
  sim   <- .ocnre_sim(12L, arm = "p", N = 50L, J = 4L, n_g = 4L)
  model <- .ocnre_model(sim, ~ det_cov1 + (1 | habitat), ~ pos_cov1)
  sp    <- .ocnre_spec(model)
  spec0 <- tulpaObs:::.tobs_occu_cover_nuts_spec(model)   # no re_blocks entry

  set.seed(3)
  th0 <- c(stats::rnorm(sp$n_par - 1L, 0, 0.3), log(0.4))
  G   <- sp$blocks[[1L]]$n_groups

  cc0 <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec0, th0, 5, 5)
  # z = 0 makes every offset 0; log_sigma_re = 0 makes its prior term exactly 0.
  cc1 <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(
    sp$spec, c(th0, rep(0, G), 0), 5, 5)
  expect_identical(cc1$lp, cc0$lp)
  expect_identical(as.numeric(cc1$grad)[seq_len(sp$n_par)], as.numeric(cc0$grad))
  # The log-SD score is d eta / d log_sigma = the offset itself, so it too is
  # exactly 0 here. (The whitened z score is NOT: sigma_re = 1 at log_sigma = 0,
  # so each z_g still carries its group's summed arm score -- the block is off,
  # not inert.)
  expect_identical(as.numeric(cc1$grad)[sp$n_par + G + 1L], 0)

  # An empty block list is the same target as no entry at all.
  spec_empty <- spec0; spec_empty$re_blocks <- list()
  cce <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec_empty, th0, 5, 5)
  expect_identical(cce$lp, cc0$lp)
  expect_identical(as.numeric(cce$grad), as.numeric(cc0$grad))

  rr0 <- tulpaObs:::.tobs_occu_cover_nuts_logpost(th0, model, 5, 5)
  rre <- tulpaObs:::.tobs_occu_cover_nuts_logpost(th0, model, 5, 5, re = list())
  expect_identical(rr0, rre)
})


test_that("occu_cover NUTS RE: the group SD's prior is N(0, sigma.logre^2) on the sampled coordinate", {
  skip_on_cran()
  sim   <- .ocnre_sim(13L, arm = "p", N = 50L, J = 4L, n_g = 4L)
  model <- .ocnre_model(sim, ~ det_cov1 + (1 | habitat), ~ pos_cov1)
  sp    <- .ocnre_spec(model)
  G     <- sp$blocks[[1L]]$n_groups
  lsd   <- sp$blocks[[1L]]$sigma_lsd
  expect_equal(lsd, 1.5)

  set.seed(5)
  th0  <- c(stats::rnorm(sp$n_par - 1L, 0, 0.3), log(0.4))
  at   <- function(ls) tulpaObs:::cpp_occu_cover_nuts_joint_logpost(
    sp$spec, c(th0, rep(0, G), ls), 5, 5)
  base <- at(0)

  # At z = 0 the data term does not depend on log sigma_re, so the whole
  # dependence is the prior. Both the density and its score must be exactly the
  # Gaussian -- a missing (or double-counted) change-of-variables term biases the
  # SD posterior while leaving every gradient check green.
  for (ls in c(-2, -1, -0.5, 0.5, 1, 2)) {
    cur <- at(ls)
    expect_equal(cur$lp - base$lp, -0.5 * ls^2 / lsd^2, tolerance = 1e-10)
    expect_equal(as.numeric(cur$grad)[sp$n_par + G + 1L], -ls / lsd^2,
                 tolerance = 1e-10)
  }
})


test_that("occu_cover NUTS RE: the fit reports the block and ranef() reads it", {
  skip_on_cran()
  sim <- .ocnre_sim(14L, arm = "p", N = 80L, J = 5L, n_g = 5L)
  fit <- .ocnre_fit(sim, control = list(n.iter = 300L, n.warmup = 300L,
                                        n.chains = 2L, seed = 2L))
  expect_identical(fit$method, "nuts")
  expect_named(fit$re, "p")
  expect_identical(fit$re$p$arm, "p")
  expect_identical(fit$re$p$var, "habitat")
  expect_identical(fit$re$p$levels, paste0("hab", 1:5))
  expect_length(fit$re$p$blup, 5L)
  expect_gt(fit$re$p$sigma, 0)
  expect_true(is.finite(fit$re$p$sigma_median))
  expect_length(fit$re$p$sigma_draws, fit$n_samples)
  # The SD is SAMPLED, so it has a posterior spread and its own convergence row.
  expect_gt(fit$re$p$sigma_sd, 0)
  expect_true(is.finite(fit$nuts$re_sigma_rhat[["p"]]))
  expect_equal(fit$nuts$sigma_logre, 1.5)

  # The reported draw / vcov surface stays the coefficient block, so the
  # coefficient inference surface is unchanged by the RE.
  expect_equal(ncol(fit$draws), fit$n_params)
  expect_identical(colnames(fit$draws), names(fit$means))
  expect_equal(ncol(fit$re_draws), 6L)              # 5 whitened z + log sigma
  expect_true(is.finite(as.numeric(logLik(fit))))
  sm <- summary(fit)
  expect_equal(nrow(sm), fit$n_params)

  rf <- ranef(fit)
  expect_s3_class(rf, "data.frame")
  expect_equal(nrow(rf), 5L)
  expect_true(all(rf$arm == "p"))
  expect_identical(rf$group, paste0("hab", 1:5))

  # A no-RE fit on the same data carries no RE surface at all.
  fit0 <- .ocnre_fit(sim, detection = ~ det_cov1,
                     control = list(n.iter = 100L, n.warmup = 100L,
                                    n.chains = 1L, seed = 2L))
  expect_null(fit0$re)
  expect_null(fit0$re_draws)
})


test_that("occu_cover NUTS RE: crossed groupings on both arms each get their own SD", {
  skip_on_cran()
  skip_if_fast()
  sim <- .ocnre_sim(15L, arm = "p", N = 80L, J = 5L, n_g = 5L)
  set.seed(151L)
  sim$visit_data$region <- factor(paste0("r", sample.int(4L, 80L * 5L, TRUE)))
  fit <- .ocnre_fit(sim, detection = ~ det_cov1 + (1 | habitat) + (1 | region),
                    positive = ~ pos_cov1 + (1 | region),
                    control = list(n.iter = 400L, n.warmup = 400L,
                                   n.chains = 1L, seed = 4L))
  # Two terms share the detection arm, so both are keyed by their grouping var;
  # the lone cover-arm term keeps the bare arm key.
  expect_named(fit$re, c("p:habitat", "p:region", "pos"))
  expect_identical(vapply(fit$re, `[[`, character(1), "arm"),
                   c("p:habitat" = "p", "p:region" = "p", "pos" = "pos"))
  expect_true(all(vapply(fit$re, `[[`, numeric(1), "sigma") > 0))
  expect_equal(ncol(fit$re_draws), 5L + 1L + 4L + 1L + 4L + 1L)
  rf <- ranef(fit)
  expect_equal(nrow(rf), 5L + 4L + 4L)
})


test_that("occu_cover NUTS samples a detection-arm group SD that covers the truth", {
  skip_on_cran()
  skip_if_fast()

  seeds    <- 1:8
  sigma_re <- 0.8
  res <- t(vapply(seeds, function(s) {
    sim <- .ocnre_sim(7000L + s, arm = "p", N = 200L, J = 6L, n_g = 10L,
                      sigma_re = sigma_re, beta_occ = c(stats::qlogis(0.55), 0.8),
                      beta_p = c(0.2, 0.6), beta_pos = c(log(0.12), -0.4))
    fit <- .ocnre_fit(sim, control = list(n.iter = 1000L, n.warmup = 1000L,
                                          n.chains = 2L, seed = 1L))
    q <- stats::quantile(fit$re$p$sigma_draws, c(0.025, 0.5, 0.975))
    c(med = unname(q[2L]),
      cover = as.numeric(q[1L] <= sigma_re && sigma_re <= q[3L]),
      blup_cor = stats::cor(fit$re$p$blup,
                            .ocnre_truth_blup(sim, fit$re$p, "p")),
      div = fit$nuts$divergent_total,
      rhat = unname(fit$nuts$re_sigma_rhat[["p"]]))
  }, numeric(5)))

  # The variance component is right-skewed at ten groups, so the seed-averaged
  # POSTERIOR MEDIAN is the summary compared against the truth (a one-sided shift
  # is a property of the mean over seeds, not of one fit). Measured 2026-08-10:
  # median-of-medians 0.784 against 0.8, coverage 8/8, BLUP cor 0.94, no
  # divergences.
  expect_lt(abs(mean(res[, "med"]) - sigma_re), 0.30)
  expect_gte(mean(res[, "cover"]), 0.75)         # 95% CrI, >= 6/8 seeds
  expect_gt(mean(res[, "blup_cor"]), 0.7)
  expect_lte(max(res[, "div"]), 20)
  expect_lt(max(res[, "rhat"]), 1.05)
})


test_that("occu_cover NUTS samples a cover-arm group SD that covers the truth", {
  skip_on_cran()
  skip_if_fast()

  seeds    <- 1:8
  sigma_re <- 0.5
  res <- t(vapply(seeds, function(s) {
    sim <- .ocnre_sim(7500L + s, arm = "pos", N = 200L, J = 6L, n_g = 10L,
                      sigma_re = sigma_re, beta_occ = c(stats::qlogis(0.55), 0.8),
                      beta_p = c(0.2, 0.6), beta_pos = c(log(0.12), -0.4))
    fit <- .ocnre_fit(sim, detection = ~ det_cov1,
                      positive = ~ pos_cov1 + (1 | habitat),
                      control = list(n.iter = 1000L, n.warmup = 1000L,
                                     n.chains = 2L, seed = 1L))
    q <- stats::quantile(fit$re$pos$sigma_draws, c(0.025, 0.5, 0.975))
    c(med = unname(q[2L]),
      cover = as.numeric(q[1L] <= sigma_re && sigma_re <= q[3L]),
      blup_cor = stats::cor(fit$re$pos$blup,
                            .ocnre_truth_blup(sim, fit$re$pos, "pos")),
      div = fit$nuts$divergent_total)
  }, numeric(4)))

  # The cover arm observes a continuous response at every detected visit, so its
  # variance component is far better identified than the detection arm's: the
  # seed-averaged posterior median lands on the truth to ~0.01 and the BLUPs to
  # cor ~0.99 (measured 2026-08-10).
  expect_lt(abs(mean(res[, "med"]) - sigma_re), 0.12)
  expect_gt(mean(res[, "blup_cor"]), 0.9)
  expect_lte(max(res[, "div"]), 20)
  # Coverage measured 6/8. Both misses are boundary near-misses (CrI upper 0.485
  # and lower 0.503 against a truth of 0.5), which is what ten groups buys: the
  # simulator CENTRES its ten drawn deviations, so the realised spread the fit
  # sees is itself a draw around the population SD the interval is scored
  # against. The gate sits below the measurement rather than on it, so a real
  # regression still shows while the estimand's own noise does not.
  expect_gte(mean(res[, "cover"]), 0.6)
})


test_that("occu_cover NUTS RE: the remaining limits error, each naming itself", {
  skip_on_cran()
  sim <- .ocnre_sim(16L, arm = "p", N = 40L, J = 4L, n_g = 4L)
  ctl <- list(n.iter = 2L, n.warmup = 2L, n.chains = 1L, seed = 1L)

  # The plain Laplace route is a coefficient-marginal fit with no latent block.
  expect_error(.ocnre_fit(sim, method = "laplace",
                          control = list(max.iter = 5L)),
               "nested_laplace")

  # A random SLOPE needs the grid-integrated weighted / free-Sigma blocks.
  expect_error(
    .ocnre_fit(sim, detection = ~ det_cov1 + (0 + det_cov1 | habitat),
               control = ctl),
    "random SLOPE")
  expect_error(
    .ocnre_fit(sim, detection = ~ det_cov1 + (1 + det_cov1 | habitat),
               control = ctl),
    "random SLOPE")

  # Composed with the coupled areal field, both blocks go through the joint
  # nested-Laplace engine.
  adj <- matrix(0L, 40L, 40L)
  for (i in seq_len(39L)) adj[i, i + 1L] <- adj[i + 1L, i] <- 1L
  expect_error(
    .ocnre_fit(sim, occurrence = ~ 1 + car_proper(graph = adj), control = ctl),
    "NON-SPATIAL")
})
