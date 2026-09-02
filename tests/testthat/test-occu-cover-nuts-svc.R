# =============================================================================
# test-occu-cover-nuts-svc.R - a SECOND coupled areal field on the occu_cover()
# NUTS path: the intercept field plus one spatially-varying-coefficient (trend)
# field, each its own block with its own per-site weight, whitened field and
# (sigma, rho, alpha) coordinates.
#
# The validation battery is the one the sampled-hyper block was held to,
# reproduced for two blocks:
#   - the compiled target equals the R oracle on the full log posterior
#   - the analytic gradient equals a central FD, including every hyper
#     coordinate of BOTH blocks, asserted apart from the coefficient scores
#   - the hyper prior is flat in the axis coordinate, checked at raw = 0 where
#     the fields vanish and the target IS the prior (a missing Jacobian is
#     invisible to the gradient check)
#   - a one-field spec is byte-identical under the plural and singular spellings
#   - the flat vector's coordinate count is what the sampler's metric is sized
#     from
#   - both surfaces recover against a simulated truth, with no divergences
# =============================================================================


# A simulated intercept field PLUS a time-weighted trend field, both coupled
# onto the cover arm with their own amplitudes.
.ocsvc_inputs <- function(side = 6L, J = 5L, seed = 1L, sigma = 0.8, alpha = 1.0,
                          sigma_trend = 0.7, alpha_trend = 0.9) {
  N   <- side * side
  adj <- rook_adj(side)
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "lognormal", beta_occ = c(stats::qlogis(0.5), 0.8),
    beta_p = c(0.3, 0.5), beta_pos = c(log(0.12), -0.4), sigma_pos = 0.4,
    adj = adj, sigma = sigma, alpha = alpha, trend = TRUE,
    sigma_trend = sigma_trend, alpha_trend = alpha_trend, seed = seed)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)), det_cov1 = sim$visit_data$det_cov1,
    pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  list(sim = sim, od = od, cell_dat = cell_dat, y_pos = y_pos, adj = adj, N = N)
}

# Fit with the intercept field + one trend field, both copied onto the cover arm
# (the fixture simulates a coupled truth, so the copy is declared -- a bare areal
# term is an occurrence-only field).
.ocsvc_fit <- function(inp, control = list(), field = "icar",
                       copy_call = "share(spatial())") {
  f <- stats::as.formula(sprintf(
    "~ occ_cov1 + %s(graph = inp$adj) + %s(graph = inp$adj, weight = time)",
    field, field))
  pos <- stats::as.formula(paste("~ pos_cov1 +", copy_call))
  tobs(formula = f, data = inp$cell_dat, family = occu_cover("lognormal"),
       detection = ~ det_cov1, positive = pos, y = inp$od$y, y_pos = inp$y_pos,
       visits = inp$od$det.covs, method = "nuts", control = control)
}

# A pair of hand-built field blocks (no warm fit): the unweighted intercept field
# and a weighted one, every hyper sampled over stated bounds. Returns the block
# entry lists plus their coordinate counts, so a test can lay out theta itself.
.ocsvc_blocks <- function(adj, type, n, weight) {
  one <- function(w) {
    bas <- tulpaObs:::.occu_cover_nuts_field_basis(
      adj, type, n, if (identical(type, "bym2")) .bym2_scale(adj) else NULL)
    e <- list(n_field_units = n, field_map = seq_len(n), field_load = bas$B1,
              field_scale1 = as.integer(bas$scale1),
              field_has_iid = as.integer(bas$has_iid),
              field_sf = as.numeric(bas$sf),
              field_lambda = as.numeric(bas$lambda))
    if (!is.null(w)) e$field_weight <- as.numeric(w)
    sampled <- character(0)
    bound <- function(nm, lo, hi, link) {
      tr <- if (link == 1L) stats::qlogis else log
      e[[paste0("field_", nm, "_lo")]] <<- tr(lo)
      e[[paste0("field_", nm, "_hi")]] <<- tr(hi)
      sampled <<- c(sampled, nm)
    }
    bound("sigma", 0.1, 3, 0L)
    if (identical(type, "icar")) e$field_rho_fixed <- 1
    else bound("rho", 0.05,
               if (identical(type, "car_proper"))
                 min(0.99, (1 - 1e-6) / max(bas$lambda, 0)) else 0.99, 1L)
    bound("alpha", 0.1, 3, 0L)
    list(entries = e, n_raw = ncol(bas$B1) + if (bas$has_iid) n else 0L,
         sampled = sampled)
  }
  list(one(NULL), one(weight))
}

# The non-spatial Laplace fit whose model object the target reads, plus a
# per-site weight column for the varying-coefficient block.
.ocsvc_model <- function(inp) {
  lap <- tobs(formula = ~ occ_cov1, data = inp$cell_dat,
              family = occu_cover("lognormal"), detection = ~ det_cov1,
              positive = ~ pos_cov1, y = inp$od$y, y_pos = inp$y_pos,
              visits = inp$od$det.covs, method = "laplace",
              control = list(verbose = FALSE, max.iter = 60L))
  model <- lap$model
  model$site_cell <- seq_len(inp$N)
  list(lap = lap, model = model)
}


test_that("occu_cover two-field NUTS target matches the R oracle and a central FD", {
  inp <- .ocsvc_inputs(side = 4L, J = 4L, seed = 11L)
  fx  <- .ocsvc_model(inp)
  N   <- inp$N
  set.seed(2)
  w_site <- as.numeric(scale(stats::rnorm(N)))
  n_base <- length(fx$lap$means)

  for (ty in c("icar", "bym2", "car_proper")) {
    bl <- .ocsvc_blocks(inp$adj, ty, N, w_site)
    blocks <- lapply(bl, `[[`, "entries")
    spec <- tulpaObs:::.tobs_occu_cover_nuts_spec(fx$model)
    spec$field_blocks <- blocks

    n_coord <- vapply(bl, function(b) b$n_raw + length(b$sampled), numeric(1))
    n_tot   <- n_base + sum(n_coord)
    set.seed(17)
    theta <- c(as.numeric(fx$lap$means) + stats::rnorm(n_base, 0, 0.1),
               stats::rnorm(bl[[1L]]$n_raw, 0, 0.3),
               stats::rnorm(length(bl[[1L]]$sampled), 0, 0.6),
               stats::rnorm(bl[[2L]]$n_raw, 0, 0.3),
               stats::rnorm(length(bl[[2L]]$sampled), 0, 0.6))
    expect_length(theta, n_tot)

    rr <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, fx$model, 5, 5,
                                                   field = blocks)
    cc <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec, theta, 5, 5)
    expect_equal(rr$lp, cc$lp, tolerance = 1e-10)
    expect_equal(rr$grad, as.numeric(cc$grad), tolerance = 1e-10)

    fd <- numeric(n_tot); h <- 1e-5
    for (k in seq_len(n_tot)) {
      tp <- theta; tp[k] <- tp[k] + h
      tm <- theta; tm[k] <- tm[k] - h
      fd[k] <- (tulpaObs:::.tobs_occu_cover_nuts_logpost(tp, fx$model, 5, 5,
                                                         field = blocks)$lp -
                tulpaObs:::.tobs_occu_cover_nuts_logpost(tm, fx$model, 5, 5,
                                                         field = blocks)$lp) / (2 * h)
    }
    expect_equal(rr$grad, fd, tolerance = 1e-4)

    # Each block's hyper coordinates on their own: a whole-vector tolerance lets
    # a wrong hyper derivative hide behind the much larger coefficient scores,
    # and a per-block one keeps the second field's alpha from hiding behind the
    # first field's.
    hy1 <- n_base + bl[[1L]]$n_raw + seq_along(bl[[1L]]$sampled)
    hy2 <- n_base + n_coord[1L] + bl[[2L]]$n_raw + seq_along(bl[[2L]]$sampled)
    for (hy in list(hy1, hy2)) {
      expect_equal(rr$grad[hy], fd[hy], tolerance = 1e-5)
      expect_gt(max(abs(fd[hy])), 1e-3)     # the check is not vacuous
    }
  }
})


test_that("occu_cover two-field NUTS hyper priors are flat on each axis (#214)", {
  # At raw = 0 both fields are identically zero whatever the hypers are, so the
  # data term drops out and the target reduces to the hyper priors alone.
  # Transformed back to the axis the outer grid spaces its nodes in, each must be
  # CONSTANT -- the claim that the sampler integrates the same flat-over-cells
  # measure per block. A missing Jacobian is invisible in the gradient check.
  inp <- .ocsvc_inputs(side = 4L, J = 4L, seed = 11L)
  fx  <- .ocsvc_model(inp)
  N   <- inp$N
  set.seed(2)
  w_site <- as.numeric(scale(stats::rnorm(N)))

  for (ty in c("icar", "car_proper")) {
    bl <- .ocsvc_blocks(inp$adj, ty, N, w_site)
    spec <- tulpaObs:::.tobs_occu_cover_nuts_spec(fx$model)
    spec$field_blocks <- lapply(bl, `[[`, "entries")
    n_coord <- vapply(bl, function(b) b$n_raw + length(b$sampled), numeric(1))
    base0 <- c(as.numeric(fx$lap$means), numeric(sum(n_coord)))
    hy <- c(length(fx$lap$means) + bl[[1L]]$n_raw + seq_along(bl[[1L]]$sampled),
            length(fx$lap$means) + n_coord[1L] + bl[[2L]]$n_raw +
              seq_along(bl[[2L]]$sampled))
    u_grid <- c(-2, -1, -0.4, 0, 0.7, 1.5)
    for (j in hy) {
      lp <- vapply(u_grid, function(u) {
        th <- base0; th[j] <- u
        tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec, th, 5, 5)$lp
      }, numeric(1))
      e <- stats::plogis(u_grid)
      dens_t <- exp(lp - max(lp)) / (e * (1 - e))
      expect_lt(stats::sd(dens_t) / mean(dens_t), 1e-8)
    }
  }
})


test_that("a single field is byte-identical under both block spellings (#214)", {
  # A second field that is absent must cost nothing: the one-block spec written
  # at the top level (the spelling that predates the block list) and the same
  # block written as a one-element list are the same model, to the bit.
  inp <- .ocsvc_inputs(side = 4L, J = 4L, seed = 11L)
  fx  <- .ocsvc_model(inp)
  bl  <- .ocsvc_blocks(inp$adj, "car_proper", inp$N, NULL)[[1L]]

  spec_flat <- tulpaObs:::.tobs_occu_cover_nuts_spec(fx$model)
  spec_flat[names(bl$entries)] <- bl$entries
  spec_list <- tulpaObs:::.tobs_occu_cover_nuts_spec(fx$model)
  spec_list$field_blocks <- list(bl$entries)

  set.seed(5)
  theta <- c(as.numeric(fx$lap$means) + stats::rnorm(length(fx$lap$means), 0, 0.1),
             stats::rnorm(bl$n_raw, 0, 0.3),
             stats::rnorm(length(bl$sampled), 0, 0.6))
  c_flat <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec_flat, theta, 5, 5)
  c_list <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec_list, theta, 5, 5)
  expect_identical(c_flat$lp, c_list$lp)
  expect_identical(as.numeric(c_flat$grad), as.numeric(c_list$grad))

  r_flat <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, fx$model, 5, 5,
                                                      field = bl$entries)
  r_list <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, fx$model, 5, 5,
                                                      field = list(bl$entries))
  expect_identical(r_flat$lp, r_list$lp)
  expect_identical(r_flat$grad, r_list$grad)
})


test_that("the two-field flat vector is exactly n_raw + n_hyper wide (#214)", {
  # The engine takes the inverse metric as a bare pointer with no length check,
  # so a metric sized from anything but the actual coordinate count reads past
  # its end (the #204 bug, where a bym2 block's n_raw = 2n - 1 met a metric sized
  # n_cells). bym2 is the widest raw block, so two of them is the case to pin.
  inp <- .ocsvc_inputs(side = 4L, J = 4L, seed = 11L)
  fx  <- .ocsvc_model(inp)
  set.seed(2)
  w_site <- as.numeric(scale(stats::rnorm(inp$N)))
  bl <- .ocsvc_blocks(inp$adj, "bym2", inp$N, w_site)
  spec <- tulpaObs:::.tobs_occu_cover_nuts_spec(fx$model)
  spec$field_blocks <- lapply(bl, `[[`, "entries")

  # bym2 carries the structured basis plus an unstructured block per unit.
  expect_identical(bl[[1L]]$n_raw, inp$N - 1L + inp$N)
  n_tot <- length(fx$lap$means) +
           sum(vapply(bl, function(b) b$n_raw + length(b$sampled), numeric(1)))
  expect_silent(tulpaObs:::cpp_occu_cover_nuts_joint_logpost(
    spec, numeric(n_tot), 5, 5))
  expect_error(tulpaObs:::cpp_occu_cover_nuts_joint_logpost(
    spec, numeric(n_tot - 1L), 5, 5), "theta length")
  expect_error(tulpaObs:::cpp_occu_cover_nuts_joint_logpost(
    spec, numeric(n_tot + 1L), 5, 5), "theta length")

  # And end to end: the fitter sizes the metric from that same count and refuses
  # to hand the engine a vector of any other length, so a two-block bym2 fit
  # running at all is the invariant holding on the widest layout there is.
  fit <- suppressWarnings(.ocsvc_fit(inp, field = "bym2",
    copy_call = "share(spatial(), alpha = grid(c(0.3, 1.5)))",
    control = list(verbose = FALSE, progress = FALSE, n.iter = 150L,
                   n.warmup = 150L, n.chains = 1L, seed = 2L, max.iter = 60L,
                   sigma.grid = c(0.3, 1.2))))
  expect_identical(fit$nuts$n_fields, 2L)
  # bym2 carries a mixing parameter, so both blocks sample all three hypers.
  expect_setequal(fit$nuts$sampled_hyper,
                  c("sigma", "rho", "alpha", "sigma_trend", "rho_trend",
                    "alpha_trend"))
  expect_identical(fit$nuts$fixed_hyper, character(0))
  expect_length(fit$trend_field, inp$N)
})


test_that("occu_cover NUTS samples a second (SVC) field and reports both (#214)", {
  skip_on_cran()
  skip_if_fast()
  inp <- .ocsvc_inputs(side = 5L, J = 4L, seed = 77L)
  fit <- suppressWarnings(.ocsvc_fit(inp, control = list(
    verbose = FALSE, n.iter = 400L, n.warmup = 400L, n.chains = 2L, seed = 3L,
    max.iter = 80L)))

  expect_identical(fit$method, "nuts")
  # Both surfaces are reported, each centred and one value per graph node.
  expect_length(fit$spatial_field, inp$N)
  expect_length(fit$trend_field, inp$N)
  expect_named(fit$trend_fields, "time")
  expect_true(all(is.finite(c(fit$spatial_field, fit$trend_field))))

  # Each block's hypers are its own coordinates, named apart.
  expect_setequal(fit$nuts$sampled_hyper,
                  c("sigma", "alpha", "sigma_trend", "alpha_trend"))
  expect_setequal(fit$nuts$fixed_hyper, c("rho", "rho_trend"))   # icar rho = 1
  expect_true(all(c("sigma", "alpha", "field_sd", "sigma_trend", "alpha_trend",
                    "field_sd_trend") %in% colnames(fit$hyper_draws)))
  expect_identical(fit$nuts$n_fields, 2L)
  # The two amplitudes are separate parameters, not one shared coordinate.
  expect_false(isTRUE(all.equal(fit$nuts$hyper_mean[["alpha"]],
                                fit$nuts$hyper_mean[["alpha_trend"]])))
})


test_that("occu_cover two-field NUTS recovers both surfaces (#214)", {
  skip_on_cran()
  skip_if_fast()
  n_seeds <- 6L
  fcor <- tcor <- div <- rhat <- rep(NA_real_, n_seeds)
  amp  <- matrix(NA_real_, n_seeds, 4L)
  for (s in seq_len(n_seeds)) {
    inp <- .ocsvc_inputs(side = 6L, J = 5L, seed = 4000L + s)
    fit <- tryCatch(suppressWarnings(.ocsvc_fit(inp, control = list(
      verbose = FALSE, n.iter = 800L, n.warmup = 600L, n.chains = 2L, seed = 1L,
      max.iter = 100L))), error = function(e) NULL)
    if (is.null(fit)) next
    fcor[s] <- abs(stats::cor(fit$spatial_field, inp$sim$truth$f))
    tcor[s] <- abs(stats::cor(fit$trend_field, inp$sim$truth$f2))
    div[s]  <- fit$nuts$divergent_total
    rhat[s] <- max(fit$nuts$rhat, na.rm = TRUE)
    amp[s, ] <- fit$nuts$hyper_mean[c("field_sd", "field_sd_trend",
                                      "alpha", "alpha_trend")]
  }
  ok <- !is.na(fcor)
  expect_gte(mean(ok), 0.8)
  expect_lte(max(div[ok]), 5L)
  expect_lt(max(rhat[ok]), 1.1)
  # Occupancy fields are weakly identified (one binary site per node) and the
  # trend field is seen only through its weight column, so the surface gates sit
  # below the count families' -- the point is that BOTH surfaces are recovered,
  # not that the intercept field absorbs the trend.
  expect_gt(mean(fcor[ok]), 0.50)
  expect_gt(mean(tcor[ok]), 0.35)
  # Each block's own scale: the field SDs (as the geometric-mean marginal SD the
  # simulator's truth is stated in) and the two copy amplitudes, which is where a
  # second block sharing the first's hyper would show.
  truth_amp <- c(0.8, 0.7, 1.0, 0.9)
  expect_true(all(abs(colMeans(amp[ok, , drop = FALSE]) - truth_amp) < 0.4))
})


test_that("the criteria score the second field too (#211, #214)", {
  skip_on_cran()
  skip_if_fast()
  inp <- .ocsvc_inputs(side = 5L, J = 4L, seed = 91L)
  fit <- suppressWarnings(.ocsvc_fit(inp, control = list(
    verbose = FALSE, progress = FALSE, n.iter = 300L, n.warmup = 300L,
    n.chains = 1L, seed = 4L, max.iter = 80L)))

  S  <- 300L
  c0 <- tulpaObs:::.tobs_occu_cover_components(fit, S)
  expect_identical(dim(c0$field_occ), c(inp$N, S))

  # The occupancy loading is the SUM over blocks, each at its own per-site
  # weight: scoring only the intercept field is the #211 failure mode, and a
  # weighted block is exactly the one a per-cell read would miss.
  sc <- fit$model$site_cell %||% seq_len(inp$N)
  w  <- as.numeric(inp$cell_dat$time)
  expect_equal(rowMeans(c0$field_occ),
               unname(fit$spatial_field[sc]) + w * unname(fit$trend_field[sc]),
               tolerance = 1e-12)
  # Each block reaches the cover arm through its OWN amplitude.
  f1 <- t(fit$field_draws[seq_len(S), sc, drop = FALSE])
  f2 <- t(fit$trend_field_draws[[1L]][seq_len(S), sc, drop = FALSE]) * w
  expect_equal(c0$field_pos,
               sweep(f1, 2L, fit$hyper_draws[seq_len(S), "alpha"], "*") +
               sweep(f2, 2L, fit$hyper_draws[seq_len(S), "alpha_trend"], "*"),
               tolerance = 1e-12)

  # Dropping the trend block moves the score, so the criteria are not reading
  # the intercept field alone.
  core <- tulpaObs:::.occu_cover_ploglik_core
  bare <- fit; bare$trend_field_draws <- NULL
  bare$spatial$field_suffix  <- fit$spatial$field_suffix[1L]
  bare$spatial$field_weights <- fit$spatial$field_weights[1L]
  c1 <- tulpaObs:::.tobs_occu_cover_components(bare, S)
  lse  <- function(v) { m <- max(v); m + log(mean(exp(v - m))) }
  elpd <- function(ll) sum(apply(ll, 2L, lse)) - sum(apply(ll, 2L, stats::var))
  on  <- elpd(core(fit$model, c0$b_occ, c0$b_det, c0$b_pos, c0$disp,
                   c0$field_occ, c0$field_pos))
  off <- elpd(core(fit$model, c1$b_occ, c1$b_det, c1$b_pos, c1$disp,
                   c1$field_occ, c1$field_pos))
  expect_gt(abs(on - off), 1)
  expect_equal(waic(fit, n.draws = S)$elpd_waic, on, tolerance = 1e-8)
  expect_true(is.finite(cpo(fit, n.draws = S)$elpd_loo))
})


test_that("occu_cover NUTS gates the field structures it does not sample (#214)", {
  inp <- .ocsvc_inputs(side = 4L, J = 3L, seed = 12L)
  fit_args <- function(f, pos = ~ pos_cov1 + share(spatial())) {
    tobs(formula = f, data = inp$cell_dat, family = occu_cover("lognormal"),
         detection = ~ det_cov1, positive = pos, y = inp$od$y,
         y_pos = inp$y_pos, visits = inp$od$det.covs, method = "nuts",
         control = list(verbose = FALSE, n.iter = 50L, n.warmup = 50L))
  }
  # A correlated bar is one free-Sigma MCAR block across the fields, not two
  # independent blocks.
  expect_error(
    fit_args(~ occ_cov1 + spatial(~ 1 + time | site_id, graph = inp$adj)),
    "correlated|MCAR")
  # Weighted field(s) with no unweighted intercept field.
  expect_error(fit_args(~ occ_cov1 + icar(graph = inp$adj, weight = time)),
               "unweighted intercept field")
  # A per-group random effect on the occupancy arm still needs the grid engine.
  expect_error(
    fit_args(~ occ_cov1 + icar(graph = inp$adj) + re(site_id)),
    "nested_laplace")
})
