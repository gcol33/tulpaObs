# =============================================================================
# test-residuals-contract.R - one return shape for residuals(), across families.
#
# `residuals.tobs_fit()` documents `list(occ = , det = )`: `occ` the state-level
# series (per site, or per site and season), `det` the observation-level one
# (per visit / pass / distance bin), either NULL for a family that has no such
# level. Three shapes used to be in use -- that list, a bare matrix, and
# `list(mu = , det = )` -- and eight model types reached a fallback that
# differenced against a `z` their `fitted()` does not carry, returning an
# all-NA or an empty vector as a success.
# =============================================================================

ctl_res <- list(verbose = FALSE, progress = FALSE, max.iter = 15L)

# Every element the contract allows, and nothing else.
expect_residual_contract <- function(fit) {
  for (ty in c("deviance", "pearson", "response")) {
    r <- residuals(fit, type = ty)
    expect_true(is.list(r))
    expect_true(all(names(r) %in% c("occ", "det")))
    expect_true(!is.null(r$occ) || !is.null(r$det))
  }
  invisible(TRUE)
}


test_that("an unregistered model type is an error naming it, not an empty answer", {
  skip_on_cran()
  sim <- simulate_occu(N = 30, J = 3, seed = 2)
  fit <- tobs(~ occ_cov1, data = sim$data, detection = ~ 1, y = sim$y,
              family = occu(), method = "laplace", control = ctl_res)
  expect_residual_contract(fit)

  bad <- fit
  bad$model$model_type <- "not_a_family"
  expect_error(residuals(bad), "not_a_family")
  expect_error(residuals(bad), "no handler registered")
})


test_that("the single-season residual is unchanged by the handler extraction", {
  skip_on_cran()
  sim <- simulate_occu(N = 60, J = 3, seed = 2)
  fit <- tobs(~ occ_cov1, data = sim$data, detection = ~ 1, y = sim$y,
              family = occu(), method = "laplace", control = ctl_res)
  fv <- fitted(fit)
  y  <- fit$model$y
  z_obs <- apply(y, 1, function(row) as.integer(any(row[row >= 0] == 1)))
  eps <- 1e-10

  for (ty in c("response", "pearson", "deviance")) {
    r <- residuals(fit, type = ty)
    z <- fv$z
    occ_ref <- switch(ty,
      response = z_obs - z,
      pearson  = (z_obs - z) / sqrt(z * (1 - z) + eps),
      deviance = sign(z_obs - z) * sqrt(2 * abs(ifelse(z_obs == 1, -log(z + eps),
                                                       -log(1 - z + eps)))))
    expect_equal(r$occ, occ_ref)

    ex <- matrix(z * fv$p, nrow(y), ncol(y))
    det_ref <- matrix(NA_real_, nrow(y), ncol(y))
    seen <- y >= 0
    det_ref[seen] <- switch(ty,
      response = (y - ex)[seen],
      pearson  = ((y - ex) / sqrt(ex * (1 - ex) + eps))[seen],
      deviance = (sign(y - ex) * sqrt(2 * abs(ifelse(y == 1, -log(ex + eps),
                                                     -log(1 - ex + eps)))))[seen])
    expect_equal(r$det, det_ref)
  }
})


test_that("a count family reports its visit-level series in det, with occ NULL", {
  skip_on_cran()
  sim <- simulate_abun(N = 40, J = 3, seed = 1)
  fit <- tobs(~ abund_cov1, data = sim$data, detection = ~ 1, y = sim$y,
              family = abun(), method = "laplace", control = ctl_res)
  expect_residual_contract(fit)
  r <- residuals(fit)
  expect_null(r$occ)
  expect_equal(dim(r$det), dim(sim$y))

  # count(): one series, one value per response row, in the unit-level slot
  sc <- simulate_count(N = 60, seed = 1)
  fc <- tobs(~ x, data = sc$data, family = count("poisson"), y = sc$y,
             control = ctl_res)
  expect_residual_contract(fc)
  expect_null(residuals(fc)$det)
  expect_length(residuals(fc)$occ, 60L)
})


test_that("int_occu() residuals read the multi-source detection, not NA", {
  skip_on_cran()
  sim <- simulate_int_occu(N = 60, J = c(3, 4), seed = 3)
  fit <- tobs(~ x, data = sim$data, detection = ~ 1, y = sim$y,
              family = int_occu(), method = "laplace", control = ctl_res)
  expect_residual_contract(fit)
  r <- residuals(fit)
  # the fallback this replaces differenced against rep(NA, n_sites)
  expect_length(r$occ, fit$model$n_sites)
  expect_true(all(is.finite(r$occ)))

  # the indicator is "any source detected the site at any visit", the event
  # fitted()'s own z branch conditions on
  z_obs <- numeric(fit$model$n_sites)
  for (s in seq_along(fit$model$y_sources)) {
    ys  <- fit$model$y_sources[[s]]
    map <- fit$model$site_maps[[s]] + 1L
    for (rw in seq_len(nrow(ys))) {
      row <- ys[rw, ]
      if (any(row[row >= 0] == 1)) z_obs[map[[rw]]] <- 1
    }
  }
  expect_equal(residuals(fit, type = "response")$occ, z_obs - fitted(fit)$z)
})


test_that("the community count families score their own marginal", {
  skip_on_cran()
  sim <- simulate_ms_abun(n_species = 3, N = 40, J = 3, seed = 1)
  fit <- suppressWarnings(
    tobs(~ abund_cov1, data = sim$data, family = ms_abun(mixture = "poisson"),
         detection = ~ 1, y = sim$y, species = paste0("sp", 1:3),
         method = "laplace", control = ctl_res))
  expect_residual_contract(fit)
  r <- residuals(fit, type = "response")
  expect_null(r$occ)
  expect_equal(dim(r$det), dim(fit$model$y))
  # response residual is y - lambda * p, per species
  fv <- tulpaObs:::.tobs_fitted_ms_nmix(fit)
  J  <- dim(fit$model$y)[2L]
  for (s in 1:3) {
    mu <- pmax(matrix(fv$lambda[, s] * fv$p[, s], fit$model$n_sites, J), 1e-10)
    expect_equal(unname(r$det[, , s]), unname(fit$model$y[, , s] - mu))
  }

  cutp <- c(0, 25, 50, 75, 100)
  sd2 <- simulate_ms_distance(n_species = 3, N = 40, cutpoints = cutp,
                              transect = "line", key = "halfnorm", seed = 0L)
  fd <- suppressWarnings(
    tobs(~ abund_cov1, data = sd2$data,
         family = ms_distance(key = "halfnorm", transect = "line",
                              cutpoints = cutp),
         detection = ~ 1, y = sd2$y, species = paste0("sp", 1:3),
         method = "laplace", control = ctl_res))
  expect_residual_contract(fd)
  rd <- residuals(fd, type = "response")
  expect_null(rd$occ)
  expect_equal(dim(rd$det), dim(fd$model$y))
  fvd <- tulpaObs:::.tobs_fitted_ms_distance(fd)
  for (s in 1:3) {
    pi_mat <- t(vapply(fvd$sigma[, s], function(sg)
      tulpaObs:::.distance_pi(sg, fd$model$cutpoints, fd$model$key,
                              fd$model$transect, NULL),
      numeric(fd$model$n_bins)))
    mu <- pmax(fvd$lambda[, s] * pi_mat, 1e-10)
    expect_equal(unname(rd$det[, , s]), unname(fd$model$y[, , s] - mu))
  }
})


test_that("ms_occu_cover scores its state arm like its community siblings", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_occu_cover(n_species = 3, N = 30, J = 3,
                                positive = "lognormal", seed = 1)
  fit <- suppressWarnings(
    tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
         detection = ~ det_cov1, positive = ~ pos_cov1,
         y = sim$y, y_pos = sim$y_pos, visits = sim$visit_data,
         species = sim$species, method = "laplace", control = ctl_res))
  expect_residual_contract(fit)
  r <- residuals(fit, type = "response")
  expect_equal(dim(r$occ), c(30L, 3L))
  z_obs <- tulpaObs:::.tobs_community_ever_detected(fit$model)
  expect_equal(unname(r$occ), unname(z_obs - fitted(fit)$psi))
})


test_that("occu_cover() scores its state arm like its community sibling", {
  skip_on_cran()
  sim <- simulate_occu_cover(N = 30, J = 3, n_occ_covs = 1L, n_det_covs = 1L,
                             n_pos_covs = 1L, positive = "lognormal", seed = 1)
  fit <- tobs(occurrence = ~ occ_cov1, data = sim$data,
             family = occu_cover("lognormal"), detection = ~ det_cov1,
             positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
             visits = sim$visit_data, method = "laplace", control = ctl_res)
  expect_residual_contract(fit)
  r <- residuals(fit, type = "response")
  expect_null(r$det)
  expect_length(r$occ, fit$model$n_sites)
  expect_true(all(is.finite(r$occ)))

  any_det <- tulpaObs:::.occu_cover_visit_view(fit$model)$any_det
  expect_length(any_det, fit$model$n_sites)
  expect_true(all(any_det %in% c(0L, 1L)))
})


test_that("occu_multiscale_cover() scores the cell, not the plot", {
  skip_on_cran()
  sim <- simulate_occu_multiscale_cover(n_cells = 12L, plots_per_cell = 3L,
                                        visits_per_plot = 2L, seed = 1L)
  fit <- tobs(formula = ~ x_cell + icar(graph = sim$adj, group_var = "cell"),
             data = sim$data, family = occu_multiscale_cover(response = "lognormal"),
             detection = ~ x_pdet, availability = ~ x_plot, positive = ~ x_cov,
             y = sim$y, y_pos = sim$y_pos, method = "laplace", control = ctl_res)
  expect_residual_contract(fit)
  r <- residuals(fit, type = "response")
  expect_null(r$det)
  expect_length(r$occ, fit$model$n_cells)
  expect_true(all(is.finite(r$occ)))
})


test_that("check_model() drops Moran's I rather than throwing on it", {
  skip_on_cran()
  sim <- simulate_abun(N = 40, J = 3, seed = 1)
  fit <- tobs(~ abund_cov1, data = sim$data, detection = ~ 1, y = sim$y,
              family = abun(), method = "laplace", control = ctl_res)
  xy <- cbind(sim$data$abund_cov1, sim$data$abund_cov2)

  # an N-mixture has no state-level residual, so there is nothing to test for
  # spatial autocorrelation; the report says nothing rather than erroring after
  # printing, and `moran` comes back NULL rather than a swallowed failure
  out <- NULL
  expect_silent(invisible(capture.output(
    out <- check_model(fit, coords = xy, plot = FALSE))))
  expect_null(out$moran)

  # the occupancy family it is keyed on does report one
  so <- simulate_occu(N = 40, J = 3, seed = 2)
  fo <- tobs(~ occ_cov1, data = so$data, detection = ~ 1, y = so$y,
             family = occu(), method = "laplace", control = ctl_res)
  xyo <- cbind(so$data$occ_cov1, so$data$occ_cov2)
  invisible(capture.output(oo <- check_model(fo, coords = xyo, plot = FALSE)))
  expect_false(is.null(oo$moran))
})
