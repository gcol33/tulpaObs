# =============================================================================
# test-occu-multiscale-cover-recovery.R
# - three-level occupancy + cover hurdle, joint nested-Laplace
# cell-coupling path.
#
# Structural gates (always run):
#   - rejects a state formula with no cell-declaring areal term / no group_var
# Non-spatial Laplace path (#53): iid cells, no field; recovers the four arms.
# Recovery gate (skip_on_cran): across seeds with moderate occupancy /
# availability / detection rates and enough cells the four arms separate:
#   - point recovery of the 8 fixed-effect coefficients within 0.25 (mean)
#   - 95% Wald CI coverage of the coefficients (working-family gate, #97): pooled
#     >= 0.85, per-coordinate >= 0.65 (measured pooled ~0.95)
#   - field SHAPE recovery: mean cor(z_hat, f_true) > 0.70
#   - sigma within 0.40, alpha within 0.55 (small-field nested-Laplace
#     attenuation, matching the occu_cover spatial gate's loose hyper tols)
# =============================================================================


test_that("occu_multiscale_cover() requires a cell-declaring areal term", {
  sim <- simulate_occu_multiscale_cover(n_cells = 12L, plots_per_cell = 3L,
                                        visits_per_plot = 2L, seed = 1L)
  fam <- occu_multiscale_cover(response = "lognormal")

  # No areal term -> cells are undeclared, on either method.
  expect_error(
    tobs(formula = ~ x_cell, data = sim$data, family = fam,
         detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
         y = sim$y, y_pos = sim$y_pos, method = "nested_laplace"),
    "areal"
  )
  expect_error(
    tobs(formula = ~ x_cell, data = sim$data, family = fam,
         detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
         y = sim$y, y_pos = sim$y_pos, method = "laplace"),
    "areal"
  )

  # Areal field but no group_var naming the cell column.
  expect_error(
    tobs(formula = ~ x_cell + icar(graph = sim$adj), data = sim$data,
         family = fam, detection = ~ x_pdet, availability = ~ x_plot,
         positive = ~ x_cov, y = sim$y, y_pos = sim$y_pos,
         method = "nested_laplace"),
    "group_var"
  )
})


test_that("occu_multiscale_cover() recovers the four arms + field (nested-Laplace)", {
  skip_on_cran()
  skip_if_fast()

  truth <- list(beta_psi = c(0.0, 0.5), beta_theta = c(0.4, 0.4),
                beta_p = c(0.3, 0.5), beta_pos = c(log(0.12), -0.3),
                sigma = 0.6, alpha = 1.0)
  tv <- c(truth$beta_psi, truth$beta_theta, truth$beta_p, truth$beta_pos)
  nm <- c("psi_(Intercept)", "psi_x_cell", "theta_(Intercept)", "theta_x_plot",
          "p_(Intercept)", "p_x_pdet", "pos_(Intercept)", "pos_x_cov")
  n_seeds <- 15L

  est <- matrix(NA_real_, n_seeds, length(nm), dimnames = list(NULL, nm))
  se  <- matrix(NA_real_, n_seeds, length(nm), dimnames = list(NULL, nm))
  hyp <- matrix(NA_real_, n_seeds, 2L, dimnames = list(NULL, c("sigma", "alpha")))
  field_cor <- rep(NA_real_, n_seeds)

  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_multiscale_cover(
      n_cells = 80L, plots_per_cell = 5L, visits_per_plot = 3L,
      beta_psi = truth$beta_psi, beta_theta = truth$beta_theta,
      beta_p = truth$beta_p, beta_pos = truth$beta_pos,
      positive = "lognormal", sigma = truth$sigma, alpha = truth$alpha,
      seed = 5000L + s)

    fit <- tryCatch(suppressWarnings(tobs(
      formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
      data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
      detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
      y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
      control = list(sigma.grid = c(0.3, 0.6, 1.0),
                     alpha.grid = c(0, 0.5, 1, 2),
                     diagnose.k = FALSE, max.iter = 500L))),
      error = function(e) NULL)
    if (is.null(fit)) next

    est[s, ] <- fit$means[nm]
    se[s, ]  <- fit$sds[nm]
    hyp[s, ] <- fit$means[c("sigma", "alpha")]
    f_true <- sim$truth$f - mean(sim$truth$f)
    field_cor[s] <- stats::cor(fit$spatial_field, f_true)
  }

  ok <- stats::complete.cases(est)
  expect_gte(sum(ok), 12L)   # pipeline runs end-to-end across most seeds

  # Point recovery of every fixed-effect coefficient (mean over seeds).
  bias <- colMeans(est[ok, , drop = FALSE]) - tv
  for (j in seq_along(nm)) {
    expect_lt(abs(bias[j]), 0.25, label = paste0("bias ", nm[j]))
  }

  # 95% Wald CI coverage of the coefficients (working-family gate): pooled over
  # the eight coefficients x seeds at the 0.85 floor, 0.65 per-coordinate floor
  # for Monte-Carlo slack on the shared-field psi / p intercepts. Measured
  # pooled coverage ~0.95 at 18 seeds.
  cover <- vapply(seq_along(nm), function(j) {
    lo <- est[ok, j] - 1.96 * se[ok, j]
    hi <- est[ok, j] + 1.96 * se[ok, j]
    mean(tv[j] >= lo & tv[j] <= hi)
  }, numeric(1))
  expect_gte(mean(cover), 0.85)
  expect_gte(min(cover),  0.65)

  # Field shape + hyperparameters (loose: small-field nested-Laplace attenuation).
  expect_gt(mean(field_cor[ok]), 0.70)
  expect_lt(abs(mean(hyp[ok, "sigma"]) - truth$sigma), 0.40)
  expect_lt(abs(mean(hyp[ok, "alpha"]) - truth$alpha), 0.55)
})


test_that("occu_multiscale_cover() beta positive arm fits end-to-end", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu_multiscale_cover(
    n_cells = 50L, plots_per_cell = 4L, visits_per_plot = 3L,
    beta_pos = c(stats::qlogis(0.3), -0.3), positive = "beta", phi = 12,
    sigma = 0.6, alpha = 1.0, seed = 909L)
  fit <- suppressWarnings(tobs(
    formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
    data = sim$data, family = occu_multiscale_cover(response = "beta"),
    detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
    y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
    control = list(sigma.grid = c(0.3, 0.6, 1.0), alpha.grid = c(0, 0.5, 1, 2),
                   diagnose.k = FALSE, max.iter = 500L)))
  expect_s3_class(fit, "tobs_fit")
  expect_true(all(is.finite(fit$means)))
  expect_true(all(c("psi_(Intercept)", "theta_(Intercept)",
                    "p_(Intercept)", "pos_(Intercept)") %in% names(fit$means)))
})

test_that("occu_multiscale_cover(\"gaussian\") fits + recovers (nested_laplace, #127)", {
  skip_on_cran()
  skip_if_fast()
  # Identity-Gaussian cover arm on the shared-field joint engine: mu = eta.
  sim <- simulate_occu_multiscale_cover(
    n_cells = 60L, plots_per_cell = 4L, visits_per_plot = 2L,
    beta_psi = c(0.4, 0.6), beta_theta = c(0.4, 0.5), beta_p = c(0.3, 0.5),
    beta_pos = c(2.0, -0.4), positive = "gaussian", phi = 0.35,
    sigma = 0.6, alpha = 1.0, seed = 3L)
  fit <- suppressWarnings(tobs(
    formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
    data = sim$data, family = occu_multiscale_cover(response = "gaussian"),
    detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
    y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
    control = list(sigma.grid = c(0.3, 0.6, 1.0), alpha.grid = c(0, 0.5, 1, 2),
                   diagnose.k = FALSE, max.iter = 500L)))
  expect_s3_class(fit, "tobs_fit")
  expect_true(all(is.finite(fit$means)))
  expect_lt(abs(fit$means[["pos_(Intercept)"]] - 2.0), 0.3)
  expect_lt(abs(fit$means[["pos_x_cov"]] - (-0.4)),    0.25)
})

test_that("occu_multiscale_cover(\"gaussian\") non-spatial Laplace recovers (#127)", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 8L
  int_e <- slope_e <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_occu_multiscale_cover(
      n_cells = 80L, plots_per_cell = 4L, visits_per_plot = 4L,
      beta_psi = c(0.2, 0.6), beta_theta = c(0.5, 0.4), beta_p = c(0.3, -0.4),
      beta_pos = c(2.0, -0.4), positive = "gaussian", phi = 0.35,
      sigma = 0, alpha = 0, seed = 400L + s)
    fit <- tryCatch(suppressWarnings(tobs(
      formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
      data = sim$data, family = occu_multiscale_cover(response = "gaussian"),
      detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
      y = sim$y, y_pos = sim$y_pos, method = "laplace",
      control = list(verbose = FALSE))), error = function(e) NULL)
    if (is.null(fit)) next
    int_e[s]   <- fit$means[["pos_(Intercept)"]]
    slope_e[s] <- fit$means[["pos_x_cov"]]
  }
  expect_lt(abs(mean(int_e,   na.rm = TRUE) - 2.0),    0.15)
  expect_lt(abs(mean(slope_e, na.rm = TRUE) - (-0.4)), 0.15)
})

test_that("occu_multiscale_cover() fitted() / predict() (#53)", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_occu_multiscale_cover(
    n_cells = 40L, plots_per_cell = 4L, visits_per_plot = 3L,
    sigma = 0.6, alpha = 1.0, positive = "lognormal", seed = 77L)
  fit <- suppressWarnings(tobs(
    formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
    data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
    detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
    y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
    control = list(sigma.grid = c(0.3, 0.6, 1.0), alpha.grid = c(0, 0.5, 1, 2),
                   diagnose.k = FALSE, max.iter = 500L)))

  n_cells <- fit$model$n_cells; n_plots <- fit$model$n_plots
  fv <- fitted(fit)
  expect_named(fv, c("psi", "theta", "p", "cover", "field", "p_marginal"))
  expect_length(fv$psi, n_cells)
  expect_length(fv$theta, n_plots)
  expect_length(fv$cover, n_plots)
  expect_true(all(fv$psi > 0 & fv$psi < 1))
  expect_true(all(fv$theta > 0 & fv$theta < 1))
  expect_true(all(fv$p > 0 & fv$p < 1))
  expect_true(all(fv$cover > 0))                       # lognormal mean > 0
  expect_true(all(fv$p_marginal >= 0 & fv$p_marginal <= 1))

  # predict() routes types to the matching arm, draw-based (mean/sd/interval);
  # in-sample only. Its mean is E[response(eta)] over the coefficient draws,
  # NOT fitted()'s plug-in response(E[eta]) -- Jensen's inequality separates
  # the two under a nonlinear link (most visibly on the lognormal cover mean),
  # so they track each other only up to that gap, not to Monte Carlo noise
  # alone; assert agreement on the (near-linear) logit arms and a loose
  # correlation everywhere.
  pr_psi   <- predict(fit, type = "state")
  pr_theta <- predict(fit, type = "availability")
  pr_p     <- predict(fit, type = "detection")
  pr_cover <- predict(fit, type = "cover")
  expect_true(all(c("cell", "mean", "sd", "lwr", "upr") %in% names(pr_psi)))
  expect_equal(nrow(pr_psi), n_cells)
  expect_equal(nrow(pr_theta), n_plots)
  expect_equal(nrow(pr_p), n_plots)
  expect_equal(nrow(pr_cover), n_plots)
  expect_true(cor(pr_psi$mean, fv$psi) > 0.9)
  expect_true(cor(pr_theta$mean, fv$theta) > 0.9)
  expect_true(cor(pr_p$mean, fv$p) > 0.9)
  expect_true(cor(pr_cover$mean, fv$cover) > 0.9)
  expect_true(all(pr_psi$lwr <= pr_psi$mean & pr_psi$mean <= pr_psi$upr))

  expect_error(predict(fit, newdata = sim$data), "not supported")
  expect_error(predict(fit, type = "bogus"), "not supported")
})

test_that("occu_multiscale_cover() non-spatial Laplace recovers truth (#53)", {
  skip_on_cran()
  skip_if_fast()
  # No areal field (sigma = 0): the data-generating truth is the non-spatial
  # three-level model, so method = "laplace" should recover the four arms.
  truth <- list(beta_psi = c(0.2, 0.6), beta_theta = c(0.5, 0.4),
                beta_p = c(0.3, -0.4), beta_pos = c(stats::qlogis(0.3), -0.3))
  set.seed(303)
  n_seed <- 6L
  est <- matrix(NA_real_, n_seed, 8L)
  for (s in seq_len(n_seed)) {
    sim <- simulate_occu_multiscale_cover(
      n_cells = 80L, plots_per_cell = 4L, visits_per_plot = 4L,
      beta_psi = truth$beta_psi, beta_theta = truth$beta_theta,
      beta_p = truth$beta_p, beta_pos = truth$beta_pos,
      positive = "beta", phi = 12, sigma = 0, alpha = 0, seed = 300L + s)
    fit <- suppressWarnings(tobs(
      formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
      data = sim$data, family = occu_multiscale_cover(response = "beta"),
      detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
      y = sim$y, y_pos = sim$y_pos, method = "laplace",
      control = list(verbose = FALSE)))
    expect_s3_class(fit, "tobs_fit")
    est[s, ] <- as.numeric(fit$means[c(
      "psi_(Intercept)", "psi_x_cell", "theta_(Intercept)", "theta_x_plot",
      "p_(Intercept)", "p_x_pdet", "pos_(Intercept)", "pos_x_cov")])
  }
  tv <- c(truth$beta_psi, truth$beta_theta, truth$beta_p, truth$beta_pos)
  bias <- colMeans(est) - tv
  expect_true(all(abs(bias) < 0.3),
              info = paste(round(bias, 2), collapse = " | "))
})

test_that("occu_multiscale_cover() coupled trend field recovers its shape (#53)", {
  skip_on_cran()
  skip_if_fast()
  # A second ICAR field, weighted per cell by `tcov`, on top of the intercept
  # field: icar(... weight = tcov) in the psi formula. The trend field shape
  # should recover even though its amplitude (sigma_trend) attenuates under the
  # small-field nested-Laplace, matching occu_cover's coupled-trend behaviour.
  n_seeds <- 4L
  field_cor <- rep(NA_real_, n_seeds)
  for (s in seq_len(n_seeds)) {
    sim <- simulate_occu_multiscale_cover(
      n_cells = 60L, plots_per_cell = 5L, visits_per_plot = 3L,
      sigma = 0.6, alpha = 1.0, trend = TRUE, sigma_trend = 0.8,
      alpha_trend = 1.0, positive = "lognormal", seed = 1200L + s)
    fit <- tryCatch(suppressWarnings(tobs(
      formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell") +
                  icar(graph = sim$adj, group_var = "cell", weight = tcov),
      data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
      detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
      y = sim$y, y_pos = sim$y_pos, method = "nested_laplace",
      control = list(sigma.grid = c(0.3, 0.6, 1.0), alpha.grid = c(0, 0.5, 1, 2),
                     alpha.grid.trend = c(0, 0.5, 1, 2), diagnose.k = FALSE,
                     max.iter = 400L))),
      error = function(e) NULL)
    if (is.null(fit)) next
    if (s == 1L) {
      # The trend field is surfaced separately from the intercept field, with
      # its own sigma_trend / alpha_trend hyperparameters.
      expect_length(fit$trend_field, fit$model$n_cells)
      expect_length(fit$spatial_field, fit$model$n_cells)
      expect_true(all(c("sigma_trend", "alpha_trend") %in% names(fit$means)))
      expect_true(is.finite(fit$means[["sigma_trend"]]))
    }
    ft <- sim$truth$f_trend - mean(sim$truth$f_trend)
    field_cor[s] <- abs(stats::cor(fit$trend_field, ft))
  }
  ok <- is.finite(field_cor)
  expect_gte(sum(ok), 3L)
  expect_gt(mean(field_cor[ok]), 0.6)
})


test_that("occu_multiscale_cover() surfaces + reduces without within-plot replication (#97)", {
  # Single releves (one visit per plot) carry no within-plot replication, so the
  # availability theta and detection p are not separately identified -- the fit
  # identifies cell occupancy psi and the product theta * p (it reduces to
  # occu_cover). The dispatcher surfaces the reduction, and the product (not the
  # individual levels) is what recovers.

  # 1) Surfaced: a single-releve fit emits the identifiability note; a fit with
  #    within-plot replicate visits does not.
  sim1 <- simulate_occu_multiscale_cover(
    n_cells = 30L, plots_per_cell = 4L, visits_per_plot = 1L,
    positive = "beta", phi = 12, sigma = 0, alpha = 0, seed = 801L)
  msg1 <- testthat::capture_messages(suppressWarnings(tobs(
    ~ x_cell + icar(graph = sim1$adj, group_var = "cell"), sim1$data,
    family = occu_multiscale_cover("beta"), detection = ~ x_pdet,
    availability = ~ x_plot, positive = ~ x_cov, y = sim1$y, y_pos = sim1$y_pos,
    method = "laplace", control = list(verbose = FALSE))))
  expect_true(any(grepl("within-plot replication", msg1)))

  sim2 <- simulate_occu_multiscale_cover(
    n_cells = 30L, plots_per_cell = 4L, visits_per_plot = 3L,
    positive = "beta", phi = 12, sigma = 0, alpha = 0, seed = 802L)
  msg2 <- testthat::capture_messages(suppressWarnings(tobs(
    ~ x_cell + icar(graph = sim2$adj, group_var = "cell"), sim2$data,
    family = occu_multiscale_cover("beta"), detection = ~ x_pdet,
    availability = ~ x_plot, positive = ~ x_cov, y = sim2$y, y_pos = sim2$y_pos,
    method = "laplace", control = list(verbose = FALSE))))
  expect_false(any(grepl("within-plot replication", msg2)))

  # 2) Reduction: across seeds the cell occupancy psi and the product theta * p
  #    recover, even though theta and p separately do not.
  skip_on_cran()
  skip_if_fast()
  b_psi <- c(0.2, 0.6); b_theta <- c(0.5, 0.0); b_p <- c(0.3, 0.0)
  prod_true <- stats::plogis(b_theta[1]) * stats::plogis(b_p[1])
  n_seed <- 6L
  psi_int <- psi_x <- prod_hat <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_occu_multiscale_cover(
      n_cells = 120L, plots_per_cell = 4L, visits_per_plot = 1L,
      beta_psi = b_psi, beta_theta = b_theta, beta_p = b_p,
      beta_pos = c(stats::qlogis(0.3), 0.0), positive = "beta", phi = 12,
      sigma = 0, alpha = 0, seed = 810L + s)
    fit <- tryCatch(suppressMessages(suppressWarnings(tobs(
      ~ x_cell + icar(graph = sim$adj, group_var = "cell"), sim$data,
      family = occu_multiscale_cover("beta"), detection = ~ x_pdet,
      availability = ~ x_plot, positive = ~ x_cov, y = sim$y, y_pos = sim$y_pos,
      method = "laplace", control = list(verbose = FALSE)))), error = function(e) NULL)
    if (is.null(fit)) next
    psi_int[s]  <- fit$means[["psi_(Intercept)"]]
    psi_x[s]    <- fit$means[["psi_x_cell"]]
    th <- fit$means[["theta_(Intercept)"]]; pp <- fit$means[["p_(Intercept)"]]
    prod_hat[s] <- stats::plogis(th) * stats::plogis(pp)
  }
  ok <- is.finite(prod_hat)
  expect_gte(sum(ok), 4L)
  # psi (cell occupancy) and the identified product theta * p recover.
  expect_lt(abs(mean(psi_int[ok]) - b_psi[1]), 0.15)
  expect_lt(abs(mean(psi_x[ok])   - b_psi[2]), 0.20)
  expect_lt(abs(mean(prod_hat[ok]) - prod_true), 0.05)
})
