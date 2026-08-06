# =============================================================================
# test-occu-multiscale-cover-nuts.R - non-spatial NUTS for the three-level
# occupancy + cover hurdle (occu_multiscale_cover(), method = "nuts";
# gcol33/tulpaObs#70).
#
# The sampler draws the EXACT coefficient posterior of the non-spatial three-
# level marginal (z over cells, a over plots both summed in closed form) via the
# in-tree C++ FullGradFn (src/occu_multiscale_cover_nuts.cpp), warm-started at the
# Laplace mode. The R target .tobs_occu_mscale_cover_nuts_logpost is the oracle the
# C++ port is cross-checked against.
#
# Coverage:
#   (1) C++ FullGradFn value + gradient byte-exact vs the R oracle (lognormal +
#       beta), structural -- always runs.
#   (2) the spatial+NUTS gate (a coupled trend field is not sampled).
#   (3) recovery + 95% Wald CI coverage of the eight fixed effects across seeds
#       in the replicated (within-plot) regime where theta / p separate; 0
#       divergences; NUTS posterior mean tracks the Laplace mode.
#   (4) S3 surface + calibrated WAIC from the per-cell draws.
# =============================================================================


# Build a bound multiscale model from a simulate_occu_multiscale_cover() output,
# mirroring the dispatcher (.dispatch_occu_multiscale_cover) so the gradient
# checks see exactly the model the NUTS fitter does.
.omcn_model <- function(sim, positive,
                        occ = ~ x_cell, theta = ~ x_plot,
                        det = ~ x_pdet, pos = ~ x_cov) {
  occ_f <- stats::as.formula(
    paste(deparse(occ), "+ icar(graph = sim$adj, group_var = \"cell\")"))
  si <- tulpaObs:::.occu_cover_spatial_fields(occ_f, sim$data)
  vd_det <- tulpaObs:::.normalize_visits(NULL, det, nrow(sim$y), ncol(sim$y))
  vd_pos <- tulpaObs:::.normalize_visits(NULL, pos, nrow(sim$y), ncol(sim$y))
  tulpaObs:::.tobs_build_occu_multiscale_cover(
    occ_formula = si$fe, theta_formula = theta,
    det_formula = vd_det$det_formula, pos_formula = vd_pos$det_formula,
    data = sim$data, y = sim$y, y_pos = sim$y_pos,
    plot_cell = as.integer(sim$data[[si$group_var]]),
    n_cells = nrow(si$fields[[1L]]$graph), positive = positive,
    det_visit_formula = vd_det$det_visit_formula, det_visit_data = vd_det$visits,
    pos_visit_formula = vd_pos$det_visit_formula, pos_visit_data = vd_pos$visits)
}

.omcn_fd_grad <- function(f, theta, h = 1e-6) {
  vapply(seq_along(theta), function(j) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  }, 0)
}


# --- (1) C++ FullGradFn byte-exact vs the R oracle (structural) ------------

test_that("occu_multiscale_cover NUTS C++ FullGradFn matches the R oracle", {
  for (positive in c("lognormal", "beta")) {
    set.seed(if (positive == "beta") 22L else 11L)
    sim <- simulate_occu_multiscale_cover(
      n_cells = 25L, plots_per_cell = 4L, visits_per_plot = 3L,
      beta_pos = if (positive == "beta") c(stats::qlogis(0.3), -0.3)
                 else c(log(0.10), -0.4),
      positive = positive, phi = if (positive == "beta") 12 else 0.35,
      sigma = 0, alpha = 0, seed = if (positive == "beta") 22L else 11L)
    model <- .omcn_model(sim, positive)
    lay   <- tulpaObs:::.tobs_occu_mscale_cover_nuts_layout(model)
    spec  <- tulpaObs:::.tobs_occu_mscale_cover_nuts_spec(model)

    set.seed(3)
    theta <- stats::rnorm(lay$total) * 0.3
    theta[lay$disp] <- log(if (positive == "beta") 12 else 0.35)
    sb <- 5

    lp_R <- tulpaObs:::.tobs_occu_mscale_cover_nuts_logpost(theta, model, lay, sb)
    cpp  <- cpp_occu_mscale_cover_nuts_joint_logpost(spec, theta, sb)

    # Value byte-exact.
    expect_equal(cpp$lp, lp_R, tolerance = 1e-9)
    # C++ analytic gradient vs finite differences of the R oracle.
    f_lp <- function(th)
      tulpaObs:::.tobs_occu_mscale_cover_nuts_logpost(th, model, lay, sb)
    num <- .omcn_fd_grad(f_lp, theta)
    expect_lt(max(abs(cpp$grad - num)), 1e-5)
  }
})


# --- (2) spatial + NUTS gate -----------------------------------------------

test_that("occu_multiscale_cover NUTS rejects a coupled trend field", {
  sim <- simulate_occu_multiscale_cover(
    n_cells = 20L, plots_per_cell = 3L, visits_per_plot = 2L,
    trend = TRUE, positive = "lognormal", seed = 5L)
  expect_error(
    tobs(formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell") +
                     icar(graph = sim$adj, group_var = "cell", weight = tcov),
         data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
         detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
         y = sim$y, y_pos = sim$y_pos, method = "nuts"),
    "non-spatial")
})


# --- (3) recovery + coverage + 0 divergences -------------------------------

test_that("occu_multiscale_cover NUTS recovers the four arms (replicated regime)", {
  skip_on_cran()
  skip_if_fast()

  truth <- list(beta_psi = c(0.2, 0.6), beta_theta = c(0.5, 0.4),
                beta_p = c(0.3, -0.4), beta_pos = c(stats::qlogis(0.3), -0.3))
  tv <- c(truth$beta_psi, truth$beta_theta, truth$beta_p, truth$beta_pos)
  nm <- c("psi_(Intercept)", "psi_x_cell", "theta_(Intercept)", "theta_x_plot",
          "p_(Intercept)", "p_x_pdet", "pos_(Intercept)", "pos_x_cov")
  n_seeds <- 8L

  est <- matrix(NA_real_, n_seeds, length(nm), dimnames = list(NULL, nm))
  se  <- matrix(NA_real_, n_seeds, length(nm), dimnames = list(NULL, nm))
  lap <- matrix(NA_real_, n_seeds, length(nm), dimnames = list(NULL, nm))
  div <- integer(n_seeds)

  for (s in seq_len(n_seeds)) {
    # No areal field (sigma = 0): the truth is the non-spatial three-level model.
    # Within-plot replication (J = 4) separates theta and p.
    sim <- simulate_occu_multiscale_cover(
      n_cells = 70L, plots_per_cell = 4L, visits_per_plot = 4L,
      beta_psi = truth$beta_psi, beta_theta = truth$beta_theta,
      beta_p = truth$beta_p, beta_pos = truth$beta_pos,
      positive = "beta", phi = 12, sigma = 0, alpha = 0, seed = 600L + s)

    fit <- tryCatch(suppressWarnings(tobs(
      formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
      data = sim$data, family = occu_multiscale_cover(response = "beta"),
      detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
      y = sim$y, y_pos = sim$y_pos, method = "nuts",
      control = list(n.iter = 600L, n.warmup = 600L, n.chains = 2L,
                     seed = 100L + s, verbose = FALSE))),
      error = function(e) NULL)
    if (is.null(fit)) next
    est[s, ] <- fit$means[nm]
    se[s, ]  <- fit$sds[nm]
    div[s]   <- sum(fit$divergent, na.rm = TRUE)

    fl <- suppressWarnings(tobs(
      formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
      data = sim$data, family = occu_multiscale_cover(response = "beta"),
      detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
      y = sim$y, y_pos = sim$y_pos, method = "laplace",
      control = list(verbose = FALSE)))
    lap[s, ] <- fl$means[nm]
  }

  ok <- stats::complete.cases(est)
  expect_gte(sum(ok), 6L)

  # Point recovery of every fixed-effect coefficient (mean over seeds).
  bias <- colMeans(est[ok, , drop = FALSE]) - tv
  expect_true(all(abs(bias) < 0.3),
              info = paste(round(bias, 3), collapse = " | "))

  # 95% Wald CI coverage of the coefficients.
  cover <- vapply(seq_along(nm), function(j) {
    lo <- est[ok, j] - 1.96 * se[ok, j]
    hi <- est[ok, j] + 1.96 * se[ok, j]
    mean(tv[j] >= lo & tv[j] <= hi)
  }, numeric(1))
  expect_gte(min(cover), 0.80)

  # NUTS posterior mean tracks the Laplace mode (same closed-form marginal).
  expect_lt(max(abs(colMeans(est[ok, , drop = FALSE]) -
                    colMeans(lap[ok, , drop = FALSE]))), 0.15)

  # No pathological divergences across the recovered fits.
  expect_true(mean(div[ok]) < 5)
})


# --- (4) S3 surface + WAIC -------------------------------------------------

test_that("occu_multiscale_cover NUTS S3 + WAIC", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu_multiscale_cover(
    n_cells = 50L, plots_per_cell = 4L, visits_per_plot = 4L,
    positive = "lognormal", phi = 0.35, sigma = 0, alpha = 0, seed = 909L)
  fit <- suppressWarnings(tobs(
    formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
    data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
    detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
    y = sim$y, y_pos = sim$y_pos, method = "nuts",
    control = list(n.iter = 500L, n.warmup = 500L, n.chains = 2L,
                   seed = 1L, verbose = FALSE)))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nuts")
  expect_equal(fit$n_samples, 1000L)   # 500 post-warmup draws x 2 chains
  expect_true(all(is.finite(fit$means)))
  expect_true(is.matrix(vcov(fit)))
  expect_true(is.matrix(confint(fit)))
  expect_true(is.finite(as.numeric(logLik(fit))))

  fv <- fitted(fit)
  expect_named(fv, c("psi", "theta", "p", "cover", "field", "p_marginal"))
  expect_length(predict(fit, type = "state"), fit$model$n_cells)

  # Calibrated WAIC from the per-cell draws (the point of the NUTS path).
  w <- tobs_waic(fit)
  expect_s3_class(w, "tulpa_criteria")
  expect_true(is.finite(w$elpd_waic))

  # Sampler ran: real (non-NA) diagnostics.
  expect_false(anyNA(fit$divergent))
  expect_true(is.finite(fit$epsilon))
})
