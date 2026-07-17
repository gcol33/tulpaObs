# =============================================================================
# test-dyn-int-occu.R - multi-season integrated occupancy (dyn_int_occu();
# spOccupancy tIntPGOcc; gcol33/tulpaObs#122).
#
# The product of dynamic occupancy (a multi-season HMM: psi1, colonization gamma,
# extinction eps) and integrated occupancy (a per-season emission pooling several
# detection sources). The latent occupancy sequence integrates out by the
# two-state HMM forward; the exact marginal is maximised (optim BFGS) with an
# observed-information vcov. v1: full site/season overlap, constant transitions,
# site-level per-source detection, non-spatial laplace.
# =============================================================================

test_that("dyn_int_occu() constructor + gates", {
  f <- dyn_int_occu()
  expect_s3_class(f, "tobs_family")
  expect_equal(f$name, "dyn_int_occu")
  sim <- simulate_dyn_int_occu(N = 60, T_seasons = 3, S = 2, seed = 1)
  # colonization / extinction are required.
  expect_error(
    tobs(~ 1, data = sim$data, family = dyn_int_occu(), detection = ~ 1,
         y = sim$y, sources = sim$sources),
    "colonization")
  # NUTS is not a supported engine.
  expect_error(
    tobs(~ 1, data = sim$data, family = dyn_int_occu(), detection = ~ 1,
         colonization = ~ 1, extinction = ~ 1, y = sim$y, sources = sim$sources,
         method = "nuts"),
    "laplace")
  # A single source is not integrated.
  expect_error(
    tobs(~ 1, data = sim$data, family = dyn_int_occu(), detection = ~ 1,
         colonization = ~ 1, extinction = ~ 1, y = sim$y[1], sources = "src1"),
    ">= 2")
})

test_that("dyn_int_occu() fits + full S3 surface", {
  sim <- simulate_dyn_int_occu(N = 250, T_seasons = 4, S = 2, J = 3,
                               psi1 = 0.5, gamma = 0.3, eps = 0.2,
                               p = c(0.4, 0.6), seed = 3)
  fit <- tobs(~ 1, data = sim$data, family = dyn_int_occu(), detection = ~ 1,
              colonization = ~ 1, extinction = ~ 1, y = sim$y,
              sources = sim$sources, control = list(verbose = FALSE, progress = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_true(isTRUE(fit$convergence$converged))
  expect_true(all(c("psi1_(Intercept)", "gamma_(Intercept)", "eps_(Intercept)",
                    "p_src1_(Intercept)", "p_src2_(Intercept)") %in% names(fit$means)))
  fv <- fitted(fit)
  expect_named(fv, c("psi1", "gamma", "eps", "p"))
  expect_true(all(fv$psi1 > 0 & fv$psi1 < 1))
  expect_equal(ncol(fv$p), 2L)
  expect_length(predict(fit, type = "state"), 250L)
  expect_true(is.finite(predict(fit, type = "colonization")[1]))
  w <- tobs_waic(fit, n.draws = 100L)
  expect_true(is.finite(w$waic))
  s2 <- simulate(fit, nsim = 1)
  expect_length(s2, 2L)
  expect_equal(dim(s2[[1]]), dim(sim$y[[1]]))
  expect_length(residuals(fit)$occ, 250L)
})

test_that("dyn_int_occu() recovers psi1 / gamma / eps + per-source detection", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 20L
  ps <- gm <- ep <- p1 <- p2 <- rep(NA_real_, n_seed)
  # 95% Wald coverage of each intercept on the logit scale, pooled over params.
  truth_logit <- c(psi1 = stats::qlogis(0.5), gamma = stats::qlogis(0.3),
                   eps = stats::qlogis(0.2), p_src1 = stats::qlogis(0.35),
                   p_src2 = stats::qlogis(0.6))
  cov_hits <- 0L; cov_tot <- 0L
  for (s in seq_len(n_seed)) {
    sim <- simulate_dyn_int_occu(N = 350, T_seasons = 5, S = 2, J = 3,
                                 psi1 = 0.5, gamma = 0.3, eps = 0.2,
                                 p = c(0.35, 0.6), seed = 400 + s)
    fit <- tryCatch(
      tobs(~ 1, data = sim$data, family = dyn_int_occu(), detection = ~ 1,
           colonization = ~ 1, extinction = ~ 1, y = sim$y, sources = sim$sources,
           control = list(verbose = FALSE, progress = FALSE)),
      error = function(e) NULL)
    if (is.null(fit) || !isTRUE(fit$convergence$converged)) next
    m <- fit$means
    ps[s] <- stats::plogis(m[["psi1_(Intercept)"]])
    gm[s] <- stats::plogis(m[["gamma_(Intercept)"]])
    ep[s] <- stats::plogis(m[["eps_(Intercept)"]])
    p1[s] <- stats::plogis(m[["p_src1_(Intercept)"]])
    p2[s] <- stats::plogis(m[["p_src2_(Intercept)"]])
    for (nm in c("psi1_(Intercept)", "gamma_(Intercept)", "eps_(Intercept)",
                 "p_src1_(Intercept)", "p_src2_(Intercept)")) {
      key <- sub("_\\(Intercept\\)", "", nm)
      lo <- m[[nm]] - 1.96 * fit$sds[[nm]]; hi <- m[[nm]] + 1.96 * fit$sds[[nm]]
      cov_tot <- cov_tot + 1L
      if (truth_logit[[key]] >= lo && truth_logit[[key]] <= hi)
        cov_hits <- cov_hits + 1L
    }
  }
  # Unbiased recovery of every arm.
  expect_lt(abs(mean(ps, na.rm = TRUE) - 0.50), 0.05)
  expect_lt(abs(mean(gm, na.rm = TRUE) - 0.30), 0.05)
  expect_lt(abs(mean(ep, na.rm = TRUE) - 0.20), 0.05)
  expect_lt(abs(mean(p1, na.rm = TRUE) - 0.35), 0.05)
  expect_lt(abs(mean(p2, na.rm = TRUE) - 0.60), 0.05)
  # 95% Wald intervals cover at the 0.85 pooled working floor (>= 20 seeds x 5).
  expect_gt(cov_hits / cov_tot, 0.85)
})
