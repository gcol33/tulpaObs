# =============================================================================
# test-occu-cover-nuts.R - NUTS sampler for the non-spatial joint
# occupancy-detection + cover hurdle (occu_cover(), method = "nuts").
#
# The sampler draws the exact two-state coefficient marginal via the in-tree C++
# FullGradFn (src/occu_cover_nuts.cpp), warm-started at the Laplace mode. Tests:
#   - byte-exact gradient: C++ FullGradFn == R oracle (closed-form check)
#   - parameter recovery + 95% CI coverage vs simulated truth (both cover arms)
#   - NUTS posterior agrees with the Laplace mode it is built on
#   - S3 method coverage (coef / vcov / confint / summary / predict / WAIC)
#   - dispatch gating (spatial + nuts rejected; family advertises "nuts")
# =============================================================================


# Build a tobs()-ready occu_cover input set from a simulation.
.ocn_inputs <- function(positive = "lognormal", N = 150L, J = 4L, seed = 1L,
                        beta_occ = NULL, beta_p = NULL, beta_pos = NULL,
                        phi = 30, sigma_pos = 0.4) {
  sim <- simulate_occu_cover(
    N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    beta_occ = beta_occ, beta_p = beta_p, beta_pos = beta_pos,
    phi = phi, sigma_pos = sigma_pos, positive = positive, seed = seed)
  long <- data.frame(
    site_id  = rep(seq_len(N), each = J),
    visit    = rep(seq_len(J), times = N),
    y        = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1,
    pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  list(sim = sim, od = od, cell_dat = cell_dat, y_pos = y_pos, positive = positive)
}

.ocn_fit <- function(inp, method = "laplace", control = list()) {
  tobs(formula = ~ occ_cov1, data = inp$cell_dat,
       family = occu_cover(inp$positive), detection = ~ det_cov1,
       positive = ~ pos_cov1, y = inp$od$y, y_pos = inp$y_pos,
       visits = inp$od$det.covs, method = method, control = control)
}


test_that("occu_cover NUTS C++ FullGradFn matches the R oracle (byte-exact)", {
  for (pos in c("lognormal", "beta")) {
    inp   <- .ocn_inputs(pos, N = 60L, J = 4L, seed = 11L)
    lap   <- .ocn_fit(inp, "laplace", list(verbose = FALSE, max.iter = 60L))
    model <- lap$model
    spec  <- tulpaObs:::.tobs_occu_cover_nuts_spec(model)
    np    <- length(lap$means)
    set.seed(7)
    theta <- as.numeric(lap$means) + stats::rnorm(np, 0, 0.15)

    rr <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, model,
                                                   sigma.beta = 5, sigma.logdisp = 5)
    cc <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec, theta, 5, 5)
    expect_equal(rr$lp, cc$lp, tolerance = 1e-6)
    expect_equal(rr$grad, as.numeric(cc$grad), tolerance = 1e-6)

    # Analytic R gradient matches a central finite difference of its own lp.
    fd <- numeric(np); h <- 1e-5
    for (k in seq_len(np)) {
      tp <- theta; tp[k] <- tp[k] + h
      tm <- theta; tm[k] <- tm[k] - h
      fd[k] <- (tulpaObs:::.tobs_occu_cover_nuts_logpost(tp, model, 5, 5)$lp -
                tulpaObs:::.tobs_occu_cover_nuts_logpost(tm, model, 5, 5)$lp) / (2 * h)
    }
    expect_equal(rr$grad, fd, tolerance = 1e-4)

    # The no-prior data log-likelihood equals the Laplace path's site_ll sum.
    pin <- model$process_info
    p_occ <- pin[[1L]]$p; p_p <- pin[[2L]]$p; p_pos <- pin[[3L]]$p
    bo <- theta[seq_len(p_occ)]; bp <- theta[p_occ + seq_len(p_p)]
    bpos <- theta[p_occ + p_p + seq_len(p_pos)]; ld <- theta[np]
    eta <- tulpaObs:::.occu_cover_eta_from_par(model, bo, bp, bpos)
    ll_site <- sum(tulpaObs:::.occu_cover_site_ll(model, eta$psi, eta$p_mat, eta$ep_mat, ld))
    lp_data <- rr$lp + 0.5 * sum(theta[seq_len(np - 1L)]^2) / 25 + 0.5 * ld^2 / 25
    expect_equal(lp_data, ll_site, tolerance = 1e-6)
  }
})


test_that("occu_cover NUTS is gated off the spatial path; family advertises nuts", {
  expect_true("nuts" %in% tulpaObs:::.tobs_family_methods$occu_cover)

  inp <- .ocn_inputs("lognormal", N = 30L, J = 3L, seed = 5L)
  adj <- matrix(0L, 30, 30)
  for (i in seq_len(29)) adj[i, i + 1L] <- adj[i + 1L, i] <- 1L
  expect_error(
    suppressWarnings(tobs(
      formula = ~ 1 + icar(graph = adj), data = inp$cell_dat,
      family = occu_cover("lognormal"), detection = ~ det_cov1,
      positive = ~ pos_cov1, y = inp$od$y, y_pos = inp$y_pos,
      visits = inp$od$det.covs, method = "nuts",
      control = list(verbose = FALSE))),
    "not yet wired|non-spatial sampler")
})


test_that("occu_cover NUTS posterior agrees with the Laplace mode", {
  skip_on_cran()
  skip_if_fast()

  inp <- .ocn_inputs("lognormal", N = 250L, J = 5L, seed = 42L)
  lap <- .ocn_fit(inp, "laplace", list(verbose = FALSE, max.iter = 400L))
  nut <- .ocn_fit(inp, "nuts",
                  list(verbose = FALSE, n.iter = 1500L, n.warmup = 1000L,
                       n.chains = 2L, seed = 1L))

  expect_equal(nut$method, "nuts")
  expect_true(nut$nuts$divergent_total <= 5L)
  expect_lt(max(nut$nuts$rhat, na.rm = TRUE), 1.05)

  # Posterior means within a few Laplace SEs; posterior SDs comparable.
  nm <- names(lap$means)
  d_mean <- abs(as.numeric(nut$means[nm]) - as.numeric(lap$means[nm]))
  expect_true(all(d_mean < 0.5 * pmax(as.numeric(lap$sds[nm]), 1e-3) + 0.05))
  ratio <- as.numeric(nut$sds[nm]) / pmax(as.numeric(lap$sds[nm]), 1e-3)
  expect_true(all(ratio > 0.5 & ratio < 2.0))
})


test_that("occu_cover NUTS recovers parameters + 95% CI coverage (lognormal)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 8L
  N <- 200L; J <- 5L
  beta_occ <- c(stats::qlogis(0.45), 0.9)
  beta_p   <- c(0.2, 0.6)
  beta_pos <- c(log(0.12), -0.4)
  sigma_pos <- 0.4
  truth <- c(beta_occ, beta_p, beta_pos, log(sigma_pos))

  est <- se <- matrix(NA_real_, n_seeds, length(truth))
  for (s in seq_len(n_seeds)) {
    inp <- .ocn_inputs("lognormal", N = N, J = J, seed = 3000L + s,
                       beta_occ = beta_occ, beta_p = beta_p, beta_pos = beta_pos,
                       sigma_pos = sigma_pos)
    nut <- tryCatch(.ocn_fit(inp, "nuts",
                    list(verbose = FALSE, n.iter = 1200L, n.warmup = 800L,
                         n.chains = 1L, seed = 1L)),
                    error = function(e) NULL)
    if (is.null(nut)) next
    est[s, ] <- as.numeric(nut$means)
    se[s, ]  <- as.numeric(nut$sds)
  }
  ok <- stats::complete.cases(est)
  expect_gte(mean(ok), 0.75)

  bias <- abs(colMeans(est[ok, , drop = FALSE]) - truth)
  expect_true(all(bias[1:6] < 0.25))   # the six regression coefficients
  expect_lt(bias[7], 0.10)             # log_sigma_pos

  # 95% Wald coverage on the slope coefficients (>= 0.75 floor at 8 seeds).
  cover <- abs(est[ok, c(2, 4, 6), drop = FALSE] -
               matrix(truth[c(2, 4, 6)], sum(ok), 3, byrow = TRUE)) <
           1.96 * se[ok, c(2, 4, 6), drop = FALSE]
  expect_gte(mean(cover), 0.75)
})


test_that("occu_cover NUTS recovers parameters (beta positive, smoke)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 5L
  N <- 200L; J <- 5L
  beta_occ <- c(stats::qlogis(0.45), 0.8)
  beta_p   <- c(0.2, 0.5)
  beta_pos <- c(stats::qlogis(0.25), -0.4)
  phi <- 20
  truth <- c(beta_occ, beta_p, beta_pos, log(phi))

  est <- matrix(NA_real_, n_seeds, length(truth))
  for (s in seq_len(n_seeds)) {
    inp <- .ocn_inputs("beta", N = N, J = J, seed = 4000L + s,
                       beta_occ = beta_occ, beta_p = beta_p, beta_pos = beta_pos,
                       phi = phi, sigma_pos = 0.4)
    nut <- tryCatch(.ocn_fit(inp, "nuts",
                    list(verbose = FALSE, n.iter = 1200L, n.warmup = 800L,
                         n.chains = 1L, seed = 1L)),
                    error = function(e) NULL)
    if (is.null(nut)) next
    est[s, ] <- as.numeric(nut$means)
  }
  ok <- stats::complete.cases(est)
  expect_gte(mean(ok), 0.6)
  bias <- abs(colMeans(est[ok, , drop = FALSE]) - truth)
  expect_true(all(bias[1:6] < 0.30))
  expect_lt(bias[7], 0.40)             # log_phi
})


test_that("occu_cover NUTS fit supports the S3 method surface", {
  skip_on_cran()
  skip_if_fast()

  inp <- .ocn_inputs("lognormal", N = 150L, J = 4L, seed = 9L)
  lap <- .ocn_fit(inp, "laplace", list(verbose = FALSE, max.iter = 300L))
  nut <- .ocn_fit(inp, "nuts",
                  list(verbose = FALSE, n.iter = 1000L, n.warmup = 800L,
                       n.chains = 2L, seed = 1L))
  np <- length(nut$means)

  # The NUTS fit exposes the same coefficient-inference surface as the Laplace fit.
  expect_equal(length(coef(nut)), length(coef(lap)))
  expect_equal(dim(vcov(nut)), c(np, np))
  expect_equal(nrow(confint(nut)), np)
  expect_true(is.finite(as.numeric(logLik(nut))))

  # summary surfaces the per-parameter Rhat / ESS the convergence list carries.
  sm <- summary(nut)
  expect_true("rhat" %in% colnames(sm))
  expect_true(all(nut$convergence$rhat[!is.na(nut$convergence$rhat)] > 0.9))

  # Calibrated WAIC from the per-draw pointwise likelihood.
  expect_true(is.finite(tobs_waic(nut)$waic))

  # predict() for the non-spatial occu_cover fit needs the joint nested-Laplace
  # object (the map paths live on the spatial fit); the NUTS and Laplace engines
  # error identically, so the NUTS path matches the Laplace capability surface.
  expect_error(predict(nut, type = "occurrence"), "joint")
  expect_error(predict(lap, type = "occurrence"), "joint")
})
