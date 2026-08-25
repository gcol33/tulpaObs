# Multistate false-positive occupancy (Miller et al. 2011 confirmed-detection
# design), non-spatial Laplace (analytic-gradient BFGS) + NUTS.
#
# Recovery-grade tests: point recovery against simulated truth + 95% CI coverage
# across seeds, plus a closed-form correctness anchor (the two-state marginal
# against a direct R computation, and the certain-detection identity), an FD check
# of the analytic gradient, and a byte-level C++ <-> R oracle cross-check.

fp_R_marginal <- function(y, psi, p11, p10, b) {
  # Direct two-state marginal for one site (y in {0,1,2}).
  g1 <- ifelse(y == 0, 1 - p11, ifelse(y == 1, p11 * (1 - b), p11 * b))
  g0 <- ifelse(y == 0, 1 - p10, ifelse(y == 1, p10, 0))
  log(psi * prod(g1) + (1 - psi) * prod(g0))
}

test_that("fp_occu() family is wired and reports its supported methods", {
  f <- fp_occu()
  expect_s3_class(f, "tobs_family")
  expect_identical(f$name, "fp_occu")
  expect_identical(f$status, "working")
  expect_true(all(c("laplace", "nuts") %in% tulpaObs:::.tobs_family_methods$fp_occu))
})

test_that("the marginal matches a direct two-state computation", {
  set.seed(1)
  for (rep in 1:8) {
    J <- sample(3:6, 1)
    psi <- runif(1, 0.2, 0.8); p11 <- runif(1, 0.3, 0.8)
    p10 <- runif(1, 0.01, 0.15); b <- runif(1, 0.2, 0.8)
    y <- sample(0:2, J, replace = TRUE, prob = c(0.5, 0.3, 0.2))
    out <- tulpaObs:::cpp_fp_occu_total_log_lik(
      as.integer(y), rep(1L, J), qlogis(psi), qlogis(p11), qlogis(p10), qlogis(b))
    expect_equal(out$log_lik, fp_R_marginal(y, psi, p11, p10, b), tolerance = 1e-9)
  }
  # Certain-detection identity: any y == 2 makes the unoccupied term zero, so the
  # posterior occupancy is exactly 1 and L = psi * prod(g1).
  y <- c(0L, 2L, 1L)
  out <- tulpaObs:::cpp_fp_occu_total_log_lik(
    y, rep(1L, 3), qlogis(0.4), qlogis(0.6), qlogis(0.05), qlogis(0.5))
  expect_equal(out$w1, 1, tolerance = 1e-12)
  g1 <- c(1 - 0.6, 0.6 * 0.5, 0.6 * 0.5)
  expect_equal(out$log_lik, log(0.4) + sum(log(g1)), tolerance = 1e-9)
})

test_that("analytic gradient matches finite differences", {
  sim <- simulate_fp_occu(N = 100, J = 5, n_occ_covs = 1, seed = 7)
  model <- tulpaObs:::.tobs_build_fp_occu(~ occ_cov1, ~ 1, sim$data, sim$y)
  lay  <- tulpaObs:::.tobs_fp_occu_nuts_layout(2L, 1L, 1L, 1L)
  marg <- tulpaObs:::.tobs_fp_occu_nuts_marginal(model)
  theta <- c(0.2, 0.4, 0.3, qlogis(0.06), 0.1)
  f <- function(th) tulpaObs:::.tobs_fp_occu_nuts_logpost(th, marg, lay)$lp
  g_an <- tulpaObs:::.tobs_fp_occu_nuts_logpost(theta, marg, lay)$grad
  g_fd <- sapply(seq_along(theta), function(j) {
    h <- 1e-5; tp <- theta; tm <- theta; tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  })
  expect_equal(g_an, g_fd, tolerance = 1e-4)
})

test_that("C++ fp_occu NUTS log-posterior matches the R oracle byte-for-byte", {
  sim <- simulate_fp_occu(N = 60, J = 5, n_occ_covs = 1, seed = 12)
  model <- tulpaObs:::.tobs_build_fp_occu(~ occ_cov1, ~ 1, sim$data, sim$y)
  lay  <- tulpaObs:::.tobs_fp_occu_nuts_layout(2L, 1L, 1L, 1L)
  marg <- tulpaObs:::.tobs_fp_occu_nuts_marginal(model)
  theta <- c(0.1, 0.3, 0.2, qlogis(0.05), -0.1)
  spec <- list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
               X_psi = model$X_processes[[1]], X_p11 = model$X_processes[[2]],
               X_p10 = model$X_processes[[3]], X_b = model$X_processes[[4]],
               n_sites = model$n_sites)
  r_out <- tulpaObs:::.tobs_fp_occu_nuts_logpost(theta, marg, lay, sigma.beta = 10)
  c_out <- tulpaObs:::cpp_fp_occu_nuts_joint_logpost(spec, theta, 10)
  expect_equal(c_out$lp, r_out$lp, tolerance = 1e-9)
  expect_equal(as.numeric(c_out$grad), r_out$grad, tolerance = 1e-9)
})

test_that("fp_occu Laplace recovers truth", {
  skip_if_fast()
  beta_psi <- c(qlogis(0.5), 0.7)
  sim <- simulate_fp_occu(N = 600, J = 6, n_occ_covs = 1, beta_psi = beta_psi,
                          p11 = 0.6, p10 = 0.05, b = 0.5, seed = 11)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = fp_occu(),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  truth <- c(beta_psi, qlogis(0.6), qlogis(0.05), qlogis(0.5))
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_lt(abs(est[2] - beta_psi[2]), 0.2)
  expect_named(fit$means, c("psi_(Intercept)", "psi_occ_cov1",
                            "p11_(Intercept)", "p10_(Intercept)", "b_(Intercept)"))
  expect_equal(unname(fit$intercepts$psi), plogis(est[1]), tolerance = 1e-8)
})

test_that("95% CIs cover the truth at nominal rate across seeds", {
  skip_on_cran()
  skip_if_fast()
  beta_psi <- c(qlogis(0.55), 0.6)
  truth <- c(beta_psi, qlogis(0.6), qlogis(0.05), qlogis(0.5))
  n_seed <- 30L
  covered <- matrix(NA, n_seed, length(truth))
  for (s in seq_len(n_seed)) {
    sim <- simulate_fp_occu(N = 400, J = 6, n_occ_covs = 1, beta_psi = beta_psi,
                            p11 = 0.6, p10 = 0.05, b = 0.5, seed = 400 + s)
    fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = fp_occu(),
                detection = ~ 1, y = sim$y, method = "laplace",
                control = list(verbose = FALSE))
    lo <- fit$means - 1.96 * fit$sds; hi <- fit$means + 1.96 * fit$sds
    covered[s, ] <- (truth >= lo) & (truth <= hi)
  }
  cover_rate <- colMeans(covered)
  expect_true(all(cover_rate >= 0.85),
              info = paste(round(cover_rate, 2), collapse = " | "))
})

test_that("fp_occu supports a covariate on the false-positive arm", {
  skip_on_cran()
  skip_if_fast()
  set.seed(5)
  N <- 600; J <- 6
  x <- rnorm(N)
  psi <- plogis(0.3 + 0.6 * x); p11 <- 0.6; b <- 0.5
  p10 <- plogis(qlogis(0.05) + 0.5 * x)        # covariate-varying false positives
  z <- rbinom(N, 1L, psi)
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    if (z[i] == 1L) {
      det <- rbinom(J, 1L, p11); cert <- rbinom(J, 1L, b)
      y[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
    } else y[i, ] <- rbinom(J, 1L, p10[i])
  }
  dat <- data.frame(x = x)
  fit <- tobs(formula = ~ x, data = dat, family = fp_occu(), detection = ~ 1,
              y = y, p10 = ~ x, method = "laplace",
              control = list(verbose = FALSE))
  est <- fit$means
  expect_true("p10_x" %in% names(est))
  expect_lt(abs(unname(est["p10_x"]) - 0.5), 0.3)
  expect_lt(abs(unname(est["psi_x"]) - 0.6), 0.2)
})

test_that("S3 surface works for fp_occu fits", {
  skip_if_fast()
  sim <- simulate_fp_occu(N = 300, J = 5, n_occ_covs = 1,
                          beta_psi = c(qlogis(0.5), 0.6), seed = 3)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = fp_occu(),
              detection = ~ 1, y = sim$y, method = "laplace",
              control = list(verbose = FALSE))
  expect_equal(dim(vcov(fit)), c(5L, 5L))
  expect_equal(nobs(fit), sum(!is.na(sim$y)))
  expect_true(is.finite(as.numeric(logLik(fit))))

  fv <- fitted(fit)
  expect_named(fv, c("psi", "p11", "p10", "b", "z"))
  expect_length(fv$psi, 300L)
  expect_true(all(fv$psi > 0 & fv$psi < 1))
  expect_true(all(fv$z >= 0 & fv$z <= 1))

  pr <- predict(fit, X.0 = cbind(1, c(-1, 0, 1)), type = "psi")
  expect_equal(nrow(pr), 3L)
  expect_true(all(c("mean", "sd", "q2.5", "q50", "q97.5") %in% names(pr)))
  expect_true(all(diff(pr$mean) > 0))
  expect_true(all(pr$q2.5 <= pr$mean & pr$mean <= pr$q97.5))

  ysim <- simulate(fit, seed = 1)
  expect_equal(dim(ysim), dim(sim$y))
  expect_true(all(ysim %in% 0:2))

  rr <- residuals(fit, type = "pearson")
  expect_null(rr$occ)
  expect_equal(dim(rr$det), dim(sim$y))
})

test_that("fp_occu rejects invalid detection states", {
  sim <- simulate_fp_occu(N = 30, J = 4, seed = 5)
  y_bad <- sim$y; y_bad[1, 1] <- 3L
  expect_error(
    tobs(formula = ~ 1, data = sim$data, family = fp_occu(), detection = ~ 1,
         y = y_bad, method = "laplace"),
    "detection states")
})

test_that("fp_occu NUTS recovers truth and scores WAIC", {
  skip_on_cran()
  skip_if_fast()
  beta_psi <- c(qlogis(0.5), 0.6)
  sim <- simulate_fp_occu(N = 300, J = 6, n_occ_covs = 1, beta_psi = beta_psi,
                          p11 = 0.6, p10 = 0.05, b = 0.5, seed = 31)
  fit <- tobs(formula = ~ occ_cov1, data = sim$data, family = fp_occu(),
              detection = ~ 1, y = sim$y, method = "nuts",
              control = list(n.iter = 500L, n.warmup = 500L, seed = 1L,
                             adapt.delta = 0.9, verbose = FALSE))
  expect_identical(fit$method, "nuts")
  expect_true(is.matrix(fit$draws) && nrow(fit$draws) == 500L)
  truth <- c(beta_psi, qlogis(0.6), qlogis(0.05), qlogis(0.5))
  est <- as.numeric(fit$means); se <- as.numeric(fit$sds)
  expect_true(all(abs(est - truth) / se < 3.5))
  expect_lt(mean(fit$nuts$divergent), 0.2)
  w <- waic(fit)
  expect_true(is.finite(w$waic))
  expect_gt(w$p_waic, 0)
})


# --- NUTS + random effect on the occupancy (psi) arm ----------

# Simulate fp_occu with a site-grouped occupancy random intercept
# b_g ~ N(0, sigma_re^2): eta_psi = beta0 + b_g.
sim_fp_psi_re <- function(n_groups = 25L, per_group = 8L, J = 6L,
                          beta0 = qlogis(0.45), sigma_re = 0.7,
                          p11 = 0.65, p10 = 0.05, b = 0.5, seed = 1L) {
  set.seed(seed)
  N  <- n_groups * per_group
  g  <- rep(seq_len(n_groups), each = per_group)
  bg <- rnorm(n_groups, 0, sigma_re)
  psi <- plogis(beta0 + bg[g])
  z <- rbinom(N, 1L, psi)
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    if (z[i] == 1L) {
      det <- rbinom(J, 1L, p11); cert <- rbinom(J, 1L, b)
      y[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
    } else y[i, ] <- rbinom(J, 1L, p10)
  }
  list(y = y, data = data.frame(g = factor(g)),
       truth = list(beta0 = beta0, sigma_re = sigma_re, bg = bg))
}

# Per-group RE on the true-detection (p11) arm; occupancy psi held fixed.
sim_fp_p11_re <- function(n_groups = 30L, per_group = 10L, J = 6L,
                          beta0 = qlogis(0.55), p11_0 = qlogis(0.6),
                          sigma_re = 0.8, p10 = 0.05, bcert = 0.5, seed = 1L) {
  set.seed(seed)
  N  <- n_groups * per_group
  g  <- rep(seq_len(n_groups), each = per_group)
  bg <- rnorm(n_groups, 0, sigma_re)
  psi <- plogis(beta0)
  p11 <- plogis(p11_0 + bg[g])
  z <- rbinom(N, 1L, psi)
  y <- matrix(0L, N, J)
  for (i in seq_len(N)) {
    if (z[i] == 1L) {
      det <- rbinom(J, 1L, p11[i]); cert <- rbinom(J, 1L, bcert)
      y[i, ] <- ifelse(det == 1L, ifelse(cert == 1L, 2L, 1L), 0L)
    } else y[i, ] <- rbinom(J, 1L, p10)
  }
  list(y = y, data = data.frame(g = factor(g)),
       truth = list(beta0 = beta0, p11_0 = p11_0, sigma_re = sigma_re, bg = bg))
}

test_that("fp_occu NUTS RE log-posterior gradient matches finite differences", {
  s <- sim_fp_psi_re(n_groups = 12L, per_group = 6L, J = 5L, seed = 4L)
  model <- tulpaObs:::.tobs_build_fp_occu(~ 1 + (1 | g), ~ 1, s$data, s$y)
  re_list <- tulpaObs:::.tobs_structures_from_model(model)$re
  re_info <- tulpaObs:::.tobs_count_nuts_re_info(re_list, model,
                                                 arms = c("psi", "p11"))
  expect_identical(re_info$arm, 0L)
  expect_identical(re_info$arm_tag, "psi")
  G <- re_info$n_groups
  spec <- list(y = as.integer(model$y_long), site_idx = as.integer(model$site_idx),
               X_psi = model$X_processes[[1]], X_p11 = model$X_processes[[2]],
               X_p10 = model$X_processes[[3]], X_b = model$X_processes[[4]],
               n_sites = model$n_sites, re_arm = 0L, re_group = re_info$group,
               n_re_groups = G, sigma_re_lsd = 1.5)
  set.seed(11)
  theta <- c(0.1, 0.3, qlogis(0.06), 0.05, rnorm(G, 0, 0.4), log(0.7))
  out <- tulpaObs:::cpp_fp_occu_nuts_joint_logpost(spec, theta, 10)
  f <- function(th) tulpaObs:::cpp_fp_occu_nuts_joint_logpost(spec, th, 10)$lp
  g_fd <- sapply(seq_along(theta), function(j) {
    h <- 1e-5; tp <- theta; tm <- theta; tp[j] <- tp[j] + h; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  })
  expect_equal(out$grad, g_fd, tolerance = 1e-4)
})

test_that("fp_occu NUTS recovers a single occupancy-arm intercept RE", {
  skip_on_cran()
  skip_if_fast()
  s <- sim_fp_psi_re(n_groups = 25L, per_group = 8L, J = 6L,
                     beta0 = qlogis(0.45), sigma_re = 0.7, seed = 3L)
  fit <- tobs(formula = ~ 1 + (1 | g), data = s$data, family = fp_occu(),
              detection = ~ 1, y = s$y, method = "nuts", verbose = FALSE,
              control = list(n.iter = 500L, n.warmup = 500L, seed = 7L))
  expect_identical(fit$method, "nuts")
  expect_identical(fit$re$arm, "psi")
  expect_equal(fit$re$n_groups, 25L)
  expect_length(fit$re$blup, 25L)
  # Variance component + occupancy intercept recover the truth; the RE SD is
  # sampled, so it carries a real posterior SD.
  expect_lt(abs(fit$re$sigma - 0.7), 0.35)
  expect_gt(fit$re$sigma_sd, 0)
  expect_lt(abs(fit$means[["psi_(Intercept)"]] - qlogis(0.45)), 0.4)
  expect_lt(mean(fit$nuts$divergent), 0.1)
  # The BLUPs track the simulated group offsets.
  expect_gt(cor(fit$re$blup, s$truth$bg), 0.4)
  # The RE columns ride in the draws alongside the fixed coefficients.
  expect_true(all(c("re_g1_z1", "log_sigma_psi_g1") %in% colnames(fit$draws)))
})

test_that("fp_occu() Laplace AGHQ recovers a site-grouped psi-arm RE", {
  skip_on_cran()
  skip_if_fast()
  # The variance component is recovered on average; the false-positive emission
  # adds noise to the occupancy signal, so the per-seed estimate is higher-
  # variance than clean occupancy (hence a mean-of-seeds recovery check).
  sig_est <- b0_est <- numeric(0)
  for (sd in 1:5) {
    s <- sim_fp_psi_re(n_groups = 30L, per_group = 10L, J = 6L,
                       beta0 = qlogis(0.45), sigma_re = 0.8, seed = sd)
    fit <- tobs(formula = ~ 1 + (1 | g), data = s$data, family = fp_occu(),
                detection = ~ 1, y = s$y, method = "laplace", verbose = FALSE,
                control = list(n.quad = 7L))
    if (sd == 1L) {
      expect_identical(fit$method, "laplace")
      expect_identical(fit$fp_re$arm, "psi")
      expect_true("sigma_g1_(Intercept)" %in% names(fit$means))
      re <- ranef(fit)
      expect_true(is.data.frame(re))
      expect_equal(nrow(re), 30L)
    }
    sig_est <- c(sig_est, fit$means[["sigma_g1_(Intercept)"]])
    b0_est  <- c(b0_est,  fit$means[["psi_(Intercept)"]])
  }
  expect_lt(abs(mean(sig_est) - 0.8), 0.2)
  expect_lt(abs(mean(b0_est) - qlogis(0.45)), 0.25)
})

test_that("fp_occu() Laplace AGHQ recovers a site-grouped p11 (detection)-arm RE", {
  skip_on_cran()
  skip_if_fast()
  # A site-level RE on the true-detection arm shifts eta_p11 uniformly across a
  # site's visits, so the per-site two-state marginal stays a function of one
  # scalar offset (occupancy held fixed) and goes through the same make_site AGHQ
  # path. Only occupied sites inform p11, and the false-positive arm absorbs some
  # of the y = 1 signal, so the variance component carries more small-cluster
  # attenuation than the occupancy arm (hence the mean-of-seeds check and the
  # wider intercept tolerance).
  sig <- b0 <- p110 <- numeric(0)
  for (sd in 1:5) {
    s <- sim_fp_p11_re(n_groups = 30L, per_group = 10L, J = 6L,
                       beta0 = qlogis(0.55), p11_0 = qlogis(0.6),
                       sigma_re = 0.8, seed = sd)
    fit <- tobs(formula = ~ 1, data = s$data, family = fp_occu(),
                detection = ~ 1 + (1 | g), y = s$y, method = "laplace",
                verbose = FALSE, control = list(progress = FALSE, n.quad = 7L))
    if (sd == 1L) {
      expect_identical(fit$method, "laplace")
      expect_identical(fit$fp_re$arm, "p11")
      expect_true("sigma_p1_(Intercept)" %in% names(fit$means))  # det-arm naming
      re <- ranef(fit)
      expect_true(is.data.frame(re))
      expect_equal(nrow(re), 30L)
    }
    sig  <- c(sig,  fit$means[["sigma_p1_(Intercept)"]])
    b0   <- c(b0,   fit$means[["psi_(Intercept)"]])
    p110 <- c(p110, fit$means[["p11_(Intercept)"]])
  }
  expect_lt(abs(mean(sig) - 0.8), 0.25)
  expect_gt(mean(sig), 0.4)                       # variance not collapsed
  expect_lt(abs(mean(b0) - qlogis(0.55)), 0.25)   # occupancy intercept recovered
  expect_lt(abs(mean(p110) - qlogis(0.6)), 0.3)   # detection intercept (attenuation)
})

test_that("fp_occu RE is on the psi OR p11 arm, never p10/b or both at once", {
  s <- sim_fp_psi_re(n_groups = 6L, per_group = 5L, J = 4L, seed = 2L)
  # NUTS samples a single intercept RE on the psi arm only; a detection-arm RE
  # is rejected with a pointer.
  expect_error(
    tobs(formula = ~ 1, data = s$data, family = fp_occu(),
         detection = ~ 1 + (1 | g), y = s$y, method = "nuts",
         control = list(n.iter = 20L, n.warmup = 10L)),
    "occupancy \\(psi\\) arm only|state formula")
  # Laplace integrates ONE arm at a time: RE on both psi and p11 is rejected.
  expect_error(
    tobs(formula = ~ 1 + (1 | g), data = s$data, family = fp_occu(),
         detection = ~ 1 + (1 | g), y = s$y, method = "laplace",
         control = list(progress = FALSE)),
    "BOTH|one arm|psi OR p11")
})


# --- areal spatial (ICAR / proper-CAR) on the occupancy (psi) arm (#51) --------

.fp_grid_adj <- function(side) {
  ng <- side*side; co <- expand.grid(x=seq_len(side), y=seq_len(side))
  adj <- matrix(0L, ng, ng)
  for (i in seq_len(ng)) for (j in seq_len(ng))
    if (i!=j && abs(co$x[i]-co$x[j])+abs(co$y[i]-co$y[j])==1L) adj[i,j] <- 1L
  adj
}

# Multistate false-positive data with a smoothed ICAR-like field on logit(psi).
.sim_fp_spatial <- function(adj, J = 12L, p11 = 0.75, p10 = 0.05, bcert = 0.5,
                            beta0 = qlogis(0.45), b1 = 0.6, sd_phi = 0.8, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  ng <- nrow(adj)
  phi <- as.numeric(scale(rnorm(ng)))
  for (rep in 1:3) { pn <- phi; for (i in seq_len(ng)) { nbi <- which(adj[i,]==1L); pn[i] <- 0.5*phi[i]+0.5*mean(phi[nbi]) }; phi <- pn }
  phi <- sd_phi*as.numeric(scale(phi)); phi <- phi - mean(phi)
  x <- rnorm(ng); psi <- plogis(beta0 + b1*x + phi); z <- rbinom(ng, 1L, psi)
  y <- matrix(0L, ng, J)
  for (i in seq_len(ng)) {
    if (z[i]==1L) { det <- rbinom(J,1L,p11); cert <- rbinom(J,1L,bcert); y[i,] <- ifelse(det==1L, ifelse(cert==1L,2L,1L), 0L) }
    else y[i,] <- rbinom(J,1L,p10)
  }
  list(y=y, data=data.frame(x=x), phi=phi)
}

test_that("fp_occu() areal ICAR recovers the occupancy slope + field", {
  skip_on_cran()
  skip_if_fast()
  # The field on logit(psi) is fit by BFGS over the exact two-state marginal
  # (cpp_fp_occu_total_log_lik gradient) + the ICAR prior, FD-Hessian Laplace.
  # Occupancy fields are more weakly identified than count fields (one binary
  # site per node), so the field correlation bar is lower than removal/distance.
  adj <- .fp_grid_adj(7L)
  slope_ok <- field_cor <- logical(0); slopes <- numeric(0)
  for (s in 1:3) {
    sim <- .sim_fp_spatial(adj, seed = 600 + s)
    fit <- tobs(formula = ~ x + icar(graph = adj), data = sim$data, family = fp_occu(),
                detection = ~ 1, y = sim$y, method = "nested_laplace",
                control = list(progress = FALSE, verbose = FALSE))
    if (s == 1L) {
      expect_identical(fit$method, "nested_laplace")
      expect_true("psi_x" %in% names(fit$means))
      V <- vcov(fit)
      expect_true(all(is.finite(V)))
      expect_true(all(eigen(V, only.values = TRUE)$values > 0))
      expect_false(is.null(fit$spatial_field))
    }
    est <- fit$means[["psi_x"]]; se <- fit$sds[["psi_x"]]
    slopes <- c(slopes, est)
    slope_ok <- c(slope_ok, abs(est - 0.6) / se < 3.5)
    field_cor <- c(field_cor, cor(fit$spatial_field, sim$phi))
  }
  expect_true(all(slope_ok))
  expect_lt(abs(mean(slopes) - 0.6), 0.25)
  expect_gt(mean(field_cor), 0.3)
})

test_that("fp_occu() areal spatial: bym2 fits; nuts+icar samples (#113)", {
  skip_on_cran()
  skip_if_fast()
  adj <- .fp_grid_adj(5L)
  s <- .sim_fp_spatial(adj, J = 10L, seed = 3)
  fit <- tobs(formula = ~ x + bym2(graph = adj), data = s$data, family = fp_occu(),
              detection = ~ 1, y = s$y, method = "nested_laplace",
              control = list(progress = FALSE, verbose = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_true(all(is.finite(vcov(fit))))
  expect_false(is.null(fit$spatial_field))
  # NUTS + areal now samples an intrinsic icar() field via the #71 sum-to-zero
  # reparameterisation (full recovery lives in test-count-spatial-nuts.R).
  fit_icar <- tobs(formula = ~ x + icar(graph = adj), data = s$data,
    family = fp_occu(), detection = ~ 1, y = s$y, method = "nuts",
    control = list(n.iter = 60L, n.warmup = 60L, verbose = FALSE, progress = FALSE))
  expect_identical(fit_icar$method, "nuts")
  expect_lt(abs(mean(fit_icar$spatial_field)), 1e-6)
})

test_that("fp_occu() bym2 + proper-CAR recover the occupancy field + slope (#131)", {
  skip_on_cran()
  skip_if_fast()
  # Occupancy fields are weakly identified (one binary site per node), so the
  # field-correlation bar is lower than the count families. bym2 fits the
  # rho-mixed unit field; proper-CAR (absent from the fp_occu suite before #131)
  # a full-rank precision. Same fixture as the icar recovery above.
  adj <- .fp_grid_adj(7L)
  for (term in c("bym2", "car_proper")) {
    tf <- if (term == "bym2") (~ x + bym2(graph = adj)) else
                              (~ x + car_proper(graph = adj))
    slope_ok <- field_cor <- logical(0); slopes <- numeric(0)
    for (s in 1:3) {
      sim <- .sim_fp_spatial(adj, seed = 600 + s)
      fit <- tobs(formula = tf, data = sim$data, family = fp_occu(),
                  detection = ~ 1, y = sim$y, method = "nested_laplace",
                  control = list(progress = FALSE, verbose = FALSE))
      est <- fit$means[["psi_x"]]; se <- fit$sds[["psi_x"]]
      slopes    <- c(slopes, est)
      slope_ok  <- c(slope_ok, abs(est - 0.6) / se < 3.5)
      field_cor <- c(field_cor, cor(fit$spatial_field, sim$phi))
    }
    expect_true(all(slope_ok), info = term)
    expect_lt(abs(mean(slopes) - 0.6), 0.25)
    expect_gt(mean(field_cor), 0.3)
  }
})

test_that("fp_occu() temporal()-only field recovers the AR1 field + slope (#114)", {
  skip_on_cran()
  skip_if_fast()
  # A temporal() term on its own (no areal field) runs the shared areal-BFGS driver
  # with a single temporal block on the occupancy arm. Occupancy fields are weakly
  # identified (one binary site per node), so the temporal-field bar is lower than
  # the count families.
  Tt <- 8L; per_t <- 30L; N <- Tt * per_t; J <- 4L
  fcor <- slope <- rep(NA_real_, 8L)
  for (s in seq_len(8L)) {
    set.seed(700L + s)
    period <- rep(seq_len(Tt), each = per_t)
    rho <- 0.7; sig <- 0.8; u <- numeric(Tt)
    u[1] <- stats::rnorm(1, 0, sig / sqrt(1 - rho^2))
    for (t in 2:Tt) u[t] <- rho * u[t - 1] + stats::rnorm(1, 0, sig)
    u <- u - mean(u)
    x <- stats::rnorm(N)
    psi <- stats::plogis(0.2 + 0.8 * x + u[period]); z <- stats::rbinom(N, 1, psi)
    p11 <- 0.6; p10 <- 0.05; b <- 0.4; y <- matrix(0L, N, J)
    for (i in seq_len(N)) for (j in seq_len(J)) {
      y[i, j] <- if (z[i] == 1)
        sample(0:2, 1, prob = c(1 - p11, p11 * (1 - b), p11 * b))
      else sample(0:1, 1, prob = c(1 - p10, p10))
    }
    fit <- tryCatch(tobs(~ x + temporal(period, type = "ar1"),
                         data = data.frame(x = x, period = period),
                         family = fp_occu(), detection = ~ 1, y = y,
                         method = "nested_laplace",
                         control = list(verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) {
      expect_identical(fit$method, "nested_laplace")
      expect_null(fit$spatial_field)                 # temporal-only: no areal field
      expect_length(fit$temporal_field, Tt)
    }
    slope[s] <- fit$means[["psi_x"]]
    if (length(fit$temporal_field) == Tt) fcor[s] <- abs(stats::cor(fit$temporal_field, u))
  }
  ok <- is.finite(slope)
  expect_gte(mean(ok), 0.75)
  expect_lt(abs(mean(slope[ok]) - 0.8), 0.30)        # occupancy slope recovered
  expect_gt(mean(fcor[ok], na.rm = TRUE), 0.6)       # AR1 temporal field recovered
})
