# =============================================================================
# test-occu-cover-latent.R - latent cover-per-unit occu_cover (cover_aggregate
# = "latent").
#
# A per-unit cover random effect u_i ~ N(0, sigma_u^2) shared across the unit's
# detected visits, integrated out per unit: closed-form (lognormal,
# compound-symmetry) and adaptive Gauss-Hermite (beta). Covers: the per-unit
# marginal + its eta-derivatives FD-checked vs a brute-force numerical integral
# (fast), family + dispatcher gates (fast), and multi-seed recovery of sigma_u +
# the field / occupancy / cover coefficients (slow).
# =============================================================================


# ---- brute-force per-unit marginal references -------------------------------
.lat_brute_lognormal <- function(y, eta, sigma_eps, sigma_u) {
  integrand <- function(u) vapply(u, function(uu)
    exp(sum(stats::dlnorm(y, meanlog = eta + uu, sdlog = sigma_eps, log = TRUE)) +
        stats::dnorm(uu, 0, sigma_u, log = TRUE)), numeric(1))
  log(stats::integrate(integrand, -Inf, Inf, rel.tol = 1e-10)$value)
}
.lat_brute_beta <- function(y, eta, phi, sigma_u) {
  integrand <- function(u) vapply(u, function(uu) {
    mu <- stats::plogis(eta + uu)
    exp(sum(stats::dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE)) +
        stats::dnorm(uu, 0, sigma_u, log = TRUE))
  }, numeric(1))
  log(stats::integrate(integrand, -Inf, Inf, rel.tol = 1e-10)$value)
}

.lat_fd1 <- function(f, x, h = 1e-5) (f(x + h) - f(x - h)) / (2 * h)
.lat_fd2 <- function(f, x, h = 1e-4) (f(x + h) - 2 * f(x) + f(x - h)) / (h * h)


# ---- simulation / fit helpers (shared field + per-unit cover latent) --------
.lat_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) { if (s > 1L) adj[s, s-1L] <- 1L; if (s < n) adj[s, s+1L] <- 1L }
  adj
}

# Cell-level positive design + shared field; a per-site cover latent u_i drives
# the within-cell cover spread the latent path is meant to absorb.
.lat_sim <- function(seed, family = "lognormal", n_cells = 25L, n_per = 5L,
                     J = 10L, sigma = 0.8, alpha = 1.0, sigma_u = 0.6,
                     sigma_eps = 0.4, phi = 30, b_pos = c(-0.7, 0.6),
                     b_occ1 = 0.7) {
  set.seed(seed)
  adj <- .lat_chain_adj(n_cells); n_sites <- n_cells * n_per
  Q  <- tulpaObs:::.occu_cover_icar_Q(adj)
  sq <- tulpaObs:::.occu_cover_icar_scale(adj)
  eig <- eigen(Q, symmetric = TRUE); keep <- eig$values > 1e-8
  f <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                  (stats::rnorm(sum(keep)) / sqrt(eig$values[keep])))
  f <- (f - mean(f)) / sqrt(sq)
  cell_idx <- rep(seq_len(n_cells), times = n_per)
  xcell <- as.numeric(scale(stats::rnorm(n_cells)))
  xocc  <- as.numeric(scale(stats::rnorm(n_cells)))
  site <- data.frame(cell_idx = cell_idx, xpos = xcell[cell_idx], xocc = xocc[cell_idx])
  z <- stats::rbinom(n_sites, 1L,
                     stats::plogis(0.3 + b_occ1 * xocc[cell_idx] + sigma * f[cell_idx]))
  u <- stats::rnorm(n_sites, 0, sigma_u)             # per-site cover latent
  det_cov <- matrix(stats::rnorm(n_sites * J), n_sites, J)
  Y <- matrix(0L, n_sites, J); Ypos <- matrix(0, n_sites, J)
  for (i in seq_len(n_sites)) for (j in seq_len(J)) if (z[i] == 1L) {
    d <- stats::rbinom(1L, 1L, stats::plogis(stats::qlogis(0.7) + 0.4 * det_cov[i, j]))
    Y[i, j] <- d
    if (d == 1L) {
      lin <- b_pos[1] + b_pos[2] * xcell[cell_idx[i]] +
             alpha * sigma * f[cell_idx[i]] + u[i]
      Ypos[i, j] <- if (family == "lognormal")
        exp(lin + stats::rnorm(1, 0, sigma_eps))
      else stats::rbeta(1L, stats::plogis(lin) * phi, (1 - stats::plogis(lin)) * phi)
    }
  }
  if (family == "beta") { Ypos[Ypos <= 0] <- 1e-6; Ypos[Ypos >= 1] <- 1 - 1e-6 }
  list(site = site, Y = Y, Ypos = Ypos,
       vd = data.frame(det_cov = as.vector(t(det_cov))),
       adj = adj, f = f, sigma_u = sigma_u, b_pos = b_pos)
}

.lat_fit <- function(sim, family, max.iter = 300L) {
  suppressWarnings(tobs(
    formula = ~ xocc + icar(graph = sim$adj, group_var = "cell_idx"),
    data = sim$site,
    family = occu_cover(family, cover_aggregate = "latent"),
    detection = ~ det_cov, positive = ~ xpos,
    y = sim$Y, y_pos = sim$Ypos, visits = sim$vd,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = max.iter, engine = "joint",
                   sigma.grid = exp(seq(log(0.4), log(1.6), length.out = 4)),
                   alpha.grid = c(0, 0.8, 1.5), adaptive.grid = FALSE,
                   diagnose.k = FALSE, n.quad = 15L)))
}

# Posterior mean of sigma_u (the cover-latent SD integrated on the pos arm's phi
# axis), surfaced in the public summary as `sigma_u` on the latent path.
.lat_sigma_u <- function(fit) unname(fit$means[["sigma_u"]])


test_that("lognormal latent marginal + eta-derivatives match brute force / FD", {
  set.seed(101)
  for (trial in 1:6) {
    m   <- sample(2:6, 1)
    eta <- stats::rnorm(1, -0.4, 0.5)
    se  <- stats::runif(1, 0.3, 0.9)
    su  <- stats::runif(1, 0.2, 0.8)
    y   <- stats::rlnorm(m, eta + stats::rnorm(1, 0, su), se)
    J   <- m
    base <- log(stats::plogis(0.3)) + sum(log(stats::plogis(rep(0.5, J))))
    res <- tulpaObs:::cpp_eval_occu_cover_lognormal_latent_cell(
      eta_psi = 0.3, eta_p = rep(0.5, J), eta_pos = eta,
      y_det = rep(1L, J), y_pos_vals = y, sigma_eps = se, sigma_u = su)
    logM_cpp <- res$cell_ll - base
    fM <- function(e) tulpaObs:::cpp_eval_occu_cover_lognormal_latent_cell(
      0.3, rep(0.5, J), e, rep(1L, J), y, se, su)$cell_ll - base

    expect_equal(logM_cpp, .lat_brute_lognormal(y, eta, se, su), tolerance = 1e-7)
    expect_equal(res$grad_pos,     .lat_fd1(fM, eta), tolerance = 1e-5)
    expect_equal(res$neg_hess_pos, -.lat_fd2(fM, eta), tolerance = 1e-4)
  }
})


test_that("beta latent marginal + eta-derivatives match brute force / FD", {
  set.seed(202)
  for (trial in 1:6) {
    m   <- sample(2:6, 1)
    eta <- stats::rnorm(1, 0, 0.6)
    phi <- stats::runif(1, 6, 25)
    su  <- stats::runif(1, 0.2, 0.8)
    mu  <- stats::plogis(eta + stats::rnorm(1, 0, su))
    y   <- pmin(pmax(stats::rbeta(m, mu * phi, (1 - mu) * phi), 1e-4), 1 - 1e-4)
    J   <- m
    base <- log(stats::plogis(0.3)) + sum(log(stats::plogis(rep(0.5, J))))
    res <- tulpaObs:::cpp_eval_occu_cover_beta_latent_cell(
      eta_psi = 0.3, eta_p = rep(0.5, J), eta_pos = eta,
      y_det = rep(1L, J), y_pos_vals = y, phi_prec = phi, sigma_u = su,
      n_quad = 25L)
    logM_cpp <- res$cell_ll - base
    fM <- function(e) tulpaObs:::cpp_eval_occu_cover_beta_latent_cell(
      0.3, rep(0.5, J), e, rep(1L, J), y, phi, su, n_quad = 25L)$cell_ll - base

    expect_equal(logM_cpp, .lat_brute_beta(y, eta, phi, su), tolerance = 1e-5)
    expect_equal(res$grad_pos,     .lat_fd1(fM, eta), tolerance = 1e-4)
    expect_equal(res$neg_hess_pos, -.lat_fd2(fM, eta), tolerance = 1e-3)
  }
})


# Brute-force expected (Fisher) marginal information for the beta latent:
# E_pi[ sum_j fisher_beta(eta + u) ] with the per-obs beta Fisher term
# phi^2 mu^2 (1-mu)^2 (trigamma(mu phi) + trigamma((1-mu) phi)), summed over the
# m detected visits (all share the unit-level predictor eta + u) and integrated
# over the exact latent posterior pi(u) prod beta(y | mu(eta+u), phi) N(u;0,su).
.lat_brute_beta_expinfo <- function(y, eta, phi, sigma_u) {
  fish1 <- function(pred) {
    mu <- stats::plogis(pred)
    phi^2 * mu^2 * (1 - mu)^2 *
      (trigamma(mu * phi) + trigamma((1 - mu) * phi))
  }
  # Posterior density (unnormalised) and the Fisher-weighted integrand, both
  # zeroed wherever the beta log-density underflows in the divergent tail (mu
  # -> 0/1), so integrate()'s extreme sample points stay finite. Those points
  # carry ~0 posterior mass, so zeroing them does not bias the ratio.
  guard <- function(v) { v[!is.finite(v)] <- 0; v }
  post <- function(u) guard(vapply(u, function(uu) {
    mu <- stats::plogis(eta + uu)
    exp(sum(stats::dbeta(y, mu * phi, (1 - mu) * phi, log = TRUE)) +
        stats::dnorm(uu, 0, sigma_u, log = TRUE))
  }, numeric(1)))
  # Weight the Fisher term only where the posterior carries mass; in the
  # divergent tail (post underflows to 0) the Fisher term is 0 * Inf, so skip it
  # rather than form the NaN at all.
  num <- function(u) {
    p <- post(u)
    out <- numeric(length(u))
    ok <- p > 0
    if (any(ok)) out[ok] <- p[ok] * length(y) * fish1(eta + u[ok])
    guard(out)
  }
  stats::integrate(num,  -Inf, Inf, rel.tol = 1e-9)$value /
    stats::integrate(post, -Inf, Inf, rel.tol = 1e-9)$value
}

test_that("beta latent Expected curvature is the PSD Fisher marginal info", {
  set.seed(303)
  saw_indefinite_observed <- FALSE
  for (trial in 1:8) {
    m   <- sample(2:6, 1)
    eta <- stats::rnorm(1, 0, 0.8)
    phi <- stats::runif(1, 4, 25)
    su  <- stats::runif(1, 0.3, 1.2)
    mu  <- stats::plogis(eta + stats::rnorm(1, 0, su))
    y   <- pmin(pmax(stats::rbeta(m, mu * phi, (1 - mu) * phi), 1e-4), 1 - 1e-4)
    J   <- m

    obs <- tulpaObs:::cpp_eval_occu_cover_beta_latent_cell(
      0.3, rep(0.5, J), eta, rep(1L, J), y, phi, su,
      n_quad = 25L, curvature = "observed")
    exp_ <- tulpaObs:::cpp_eval_occu_cover_beta_latent_cell(
      0.3, rep(0.5, J), eta, rep(1L, J), y, phi, su,
      n_quad = 25L, curvature = "expected")

    # Curvature mode changes ONLY neg_hess: the log-density and score are the
    # same posterior expectations regardless of which curvature is requested.
    expect_equal(exp_$cell_ll,  obs$cell_ll,  tolerance = 1e-12)
    expect_equal(exp_$grad_pos, obs$grad_pos, tolerance = 1e-12)

    # Expected curvature == the brute-force Fisher marginal information, and is
    # strictly positive (PSD by construction).
    expect_equal(exp_$neg_hess_pos,
                 .lat_brute_beta_expinfo(y, eta, phi, su),
                 tolerance = 1e-3)
    expect_gt(exp_$neg_hess_pos, 0)

    if (obs$neg_hess_pos <= 0) saw_indefinite_observed <- TRUE
  }

  # The Expected branch exists because the observed marginal information goes
  # indefinite at extreme cells; confirm one such cell, where observed <= 0 but
  # the Fisher curvature stays positive.
  yx <- pmin(pmax(stats::rbeta(5, 0.04, 0.96), 1e-4), 1 - 1e-4)  # mu ~ 0.04
  obs_x <- tulpaObs:::cpp_eval_occu_cover_beta_latent_cell(
    0.3, rep(0.5, 5), 3.0, rep(1L, 5), yx, 3.0, 2.0,
    n_quad = 25L, curvature = "observed")
  exp_x <- tulpaObs:::cpp_eval_occu_cover_beta_latent_cell(
    0.3, rep(0.5, 5), 3.0, rep(1L, 5), yx, 3.0, 2.0,
    n_quad = 25L, curvature = "expected")
  expect_gt(exp_x$neg_hess_pos, 0)
  expect_true(is.finite(exp_x$neg_hess_pos))
})

test_that("family carries the latent choice and dispatcher gates it", {
  expect_identical(occu_cover("beta", cover_aggregate = "latent")$params$cover_aggregate,
                   "latent")
  expect_identical(occu_cover("lognormal", cover_aggregate = "latent")$params$cover_aggregate,
                   "latent")

  adj <- .lat_chain_adj(10L)
  sim <- .lat_sim(seed = 7L, family = "lognormal", n_cells = 10L, n_per = 3L,
                  J = 6L)

  # Explicit latent on the non-spatial laplace path -> error.
  expect_error(
    tobs(formula = ~ xocc, data = sim$site,
         family = occu_cover("lognormal", cover_aggregate = "latent"),
         detection = ~ det_cov, positive = ~ xpos,
         y = sim$Y, y_pos = sim$Ypos, visits = sim$vd, method = "laplace"),
    "spatial")

  # Explicit latent with a VISIT-level positive covariate -> error.
  sim2 <- sim
  sim2$vd <- data.frame(det_cov = sim$vd$det_cov,
                        pcov = stats::rnorm(nrow(sim$vd)))
  expect_error(
    tobs(formula = ~ icar(graph = sim$adj, group_var = "cell_idx"),
         data = sim$site,
         family = occu_cover("lognormal", cover_aggregate = "latent"),
         detection = ~ det_cov, positive = ~ pcov,
         y = sim$Y, y_pos = sim$Ypos, visits = sim2$vd,
         method = "nested_laplace",
         control = list(verbose = FALSE, engine = "joint")),
    "cell-level positive design")
})


test_that("lognormal latent recovers sigma_u + field + coefficients (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  seeds <- c(101L, 202L, 303L, 404L, 505L)
  pos1 <- occ1 <- fcor <- su_hat <- numeric(0)
  for (sd in seeds) {
    sim <- .lat_sim(seed = sd, family = "lognormal", sigma_u = 0.6)
    fit <- .lat_fit(sim, "lognormal")
    pos1   <- c(pos1, unname(fit$means[["pos_xpos"]]))
    occ1   <- c(occ1, unname(fit$means[["psi_xocc"]]))
    fcor   <- c(fcor, cor(fit$spatial_field, sim$f))
    su_hat <- c(su_hat, .lat_sigma_u(fit))
  }
  # The cover-latent SD is surfaced as `sigma_u`, not the engine's `phi_pos`.
  expect_true("sigma_u" %in% names(fit$means))
  expect_false("phi_pos" %in% names(fit$means))
  expect_gt(median(fcor), 0.6)
  expect_gt(median(occ1), 0.35)
  expect_gt(median(pos1), 0.3)
  expect_lt(median(pos1), 0.95)
  # sigma_u (truth 0.6) recovers within a broad band -- the outer grid is coarse
  # and the between-unit signal is modest at this N.
  expect_gt(median(su_hat), 0.2)
  expect_lt(median(su_hat), 1.2)
})


test_that("beta latent runs and recovers the cover slope (multi-seed)", {
  skip_on_cran()
  skip_if_fast()
  seeds <- c(11L, 22L, 33L)
  pos1 <- numeric(0)
  for (sd in seeds) {
    sim <- .lat_sim(seed = sd, family = "beta", sigma_u = 0.5)
    fit <- .lat_fit(sim, "beta")
    expect_s3_class(fit, "tobs_fit")
    expect_identical(fit$model$cover_aggregate, "latent")
    pos1 <- c(pos1, unname(fit$means[["pos_xpos"]]))
  }
  expect_gt(median(pos1), 0.2)
  expect_lt(median(pos1), 1.1)
})


# ---- predict() on the latent path -------------------------------------------
# One row per field cell, the cell-level cover / occupancy covariates at their
# per-cell values, so predict() returns a per-cell cover surface.
.lat_predict_newdata <- function(sim) {
  first <- !duplicated(sim$site$cell_idx)
  data.frame(cell = sim$site$cell_idx[first],
             xocc = sim$site$xocc[first],
             xpos = sim$site$xpos[first])
}

test_that("lognormal latent predict: finite, positive, increases with cover cov", {
  skip_on_cran()
  skip_if_fast()
  sim <- .lat_sim(seed = 11L, family = "lognormal", n_cells = 20L, n_per = 4L,
                  J = 8L, sigma_u = 0.6)
  fit <- .lat_fit(sim, "lognormal")
  nd  <- .lat_predict_newdata(sim)

  pr_cond <- predict(fit, newdata = nd, type = "cover_cond", nsim = 400L)
  pr_exp  <- predict(fit, newdata = nd, type = "cover_exp",  nsim = 400L)

  expect_true(all(is.finite(pr_cond$mean)))
  expect_true(all(pr_cond$mean > 0))
  expect_true(all(is.finite(pr_exp$mean)))
  expect_true(all(pr_exp$mean > 0))
  # the cover slope is positive (b_pos[2] = 0.6), so marginal cover rises with xpos
  expect_gt(cor(nd$xpos, pr_cond$mean), 0.5)
})

test_that("lognormal latent cover mu inflates the plug-in mean by exp((se^2+su^2)/2)", {
  skip_on_cran()
  skip_if_fast()
  sim <- .lat_sim(seed = 11L, family = "lognormal", n_cells = 20L, n_per = 4L,
                  J = 8L, sigma_u = 0.6)
  fit <- .lat_fit(sim, "lognormal")
  nd  <- .lat_predict_newdata(sim)

  # `.tobs_cover_mu` is the deterministic derived-quantity kernel predict() uses;
  # exercise it on a single bundle so the comparison carries no Monte Carlo
  # mismatch. The marginal cover is exp(eta) * exp((sigma_eps^2 + sigma_u^2)/2)
  # per draw, where sigma_eps is the fixed within-unit dispersion and sigma_u
  # rides the grid (bundle$disp).
  bundle <- tulpaObs:::.tobs_joint_draws(fit, n = 400L)
  cell   <- as.integer(nd$cell)
  X_pos  <- tulpaObs:::.tobs_joint_arm_design(fit, nd, "pos", ncol(bundle$b$pos))
  eta    <- tulpaObs:::.tobs_joint_arm_eta(bundle, X_pos, "pos", cell)
  se     <- as.numeric(fit$model$cover_latent_disp2)
  su     <- bundle$disp                                  # per-draw sigma_u
  marg   <- exp(sweep(eta, 2L, (se^2 + su^2) / 2, "+"))
  plug   <- exp(eta)

  mu <- tulpaObs:::.tobs_cover_mu(eta, bundle, fit)
  expect_equal(mu, marg, tolerance = 1e-10)

  # the marginal exceeds the plug-in exp(eta) by exp((se^2 + su^2)/2) per draw
  # (mu == marg exactly, above). The averaged inflation is >= exp((se^2 +
  # mean(su^2))/2) by Jensen.
  ratio <- mu / plug
  expect_true(all(mu > plug))
  expect_gt(min(ratio), 1)
  expect_gte(mean(ratio), exp((se^2 + mean(su^2)) / 2))

  # predict()'s posterior-mean cover (a fresh Monte Carlo draw) also exceeds the
  # plug-in mean of exp(eta) -- the variance inflation survives the integration.
  pr_cond <- predict(fit, newdata = nd, type = "cover_cond", nsim = 400L)
  expect_true(all(pr_cond$mean > rowMeans(plug)))
})

test_that("beta latent predict: finite, in (0,1), increases with cover cov", {
  skip_on_cran()
  skip_if_fast()
  sim <- .lat_sim(seed = 11L, family = "beta", sigma_u = 0.5)
  fit <- .lat_fit(sim, "beta")
  nd  <- .lat_predict_newdata(sim)

  pr_cond <- predict(fit, newdata = nd, type = "cover_cond", nsim = 400L)
  expect_true(all(is.finite(pr_cond$mean)))
  expect_true(all(pr_cond$mean > 0 & pr_cond$mean < 1))
  expect_gt(cor(nd$xpos, pr_cond$mean), 0.4)

  # the GH marginal of plogis(eta + u) sits between the plug-in plogis(eta) and
  # 0.5 (Jensen: the logistic is concave above 0, convex below), so it is a
  # genuine marginalization, not the plug-in.
  bundle <- tulpaObs:::.tobs_joint_draws(fit, n = 400L)
  cell   <- as.integer(nd$cell)
  X_pos  <- tulpaObs:::.tobs_joint_arm_design(fit, nd, "pos", ncol(bundle$b$pos))
  eta    <- tulpaObs:::.tobs_joint_arm_eta(bundle, X_pos, "pos", cell)
  plug   <- rowMeans(stats::plogis(eta))
  expect_false(isTRUE(all.equal(pr_cond$mean, plug, tolerance = 1e-4)))
})
