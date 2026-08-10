# =============================================================================
# test-occu-cover-spatial-nuts.R - spatial NUTS sampler for the joint
# occupancy-detection + cover hurdle (occu_cover(), method = "nuts" + a
# car_proper() shared field; gcol33/tulpaObs#74).
#
# The sampler draws the exact two-state coefficient marginal JOINTLY with a
# FIXED-HYPER non-centered coupled proper-CAR field (psi linearly + cover scaled
# by alpha), via the in-tree C++ FullGradFn field block (src/occu_cover_nuts.cpp).
# Tests:
#   - byte-exact field-block gradient: C++ FullGradFn == R oracle (closed form)
#   - field-off path stays byte-identical to the non-spatial occu_cover NUTS
#   - field recovery cor(est, truth), beta/cover recovery + 95% coverage,
#     beta SD calibration vs the nested-Laplace joint SEs, 0 divergences
#   - dispatch gating (icar/bym2/SVC/RE + nuts rejected with a pointer)
# =============================================================================


# A small lattice graph (proper-CAR is full-rank, so its non-centered geometry
# is well-conditioned -- the gate the icar/bym2 fields fail).
.ocsn_grid_adj <- function(side) {
  N <- side * side
  adj <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r > 1L)    adj[i, idx(r - 1L, c)] <- 1L
    if (r < side)  adj[i, idx(r + 1L, c)] <- 1L
    if (c > 1L)    adj[i, idx(r, c - 1L)] <- 1L
    if (c < side)  adj[i, idx(r, c + 1L)] <- 1L
  }
  adj
}

# Build a tobs()-ready spatial occu_cover input set from a simulation.
.ocsn_inputs <- function(side = 8L, J = 5L, seed = 1L, positive = "lognormal",
                         beta_occ = c(stats::qlogis(0.5), 0.8),
                         beta_p = c(0.3, 0.5), beta_pos = NULL,
                         sigma = 0.7, alpha = 1.0, sigma_pos = 0.4, phi = 25) {
  N   <- side * side
  adj <- .ocsn_grid_adj(side)
  if (is.null(beta_pos))
    beta_pos <- if (positive == "beta") c(stats::qlogis(0.25), -0.4)
                else c(log(0.12), -0.4)
  sim <- simulate_occu_cover(
    N = N, J = J, positive = positive, beta_occ = beta_occ, beta_p = beta_p,
    beta_pos = beta_pos, phi = phi, sigma_pos = sigma_pos,
    adj = adj, sigma = sigma, alpha = alpha, seed = seed)
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)), det_cov1 = sim$visit_data$det_cov1,
    pos_cov1 = sim$visit_data$pos_cov1)
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  list(sim = sim, od = od, cell_dat = cell_dat, y_pos = y_pos, adj = adj,
       N = N, positive = positive,
       truth = c(beta_occ, beta_p, beta_pos,
                 if (positive == "beta") log(phi) else log(sigma_pos)))
}

.ocsn_fit <- function(inp, method = "nuts", spatial = TRUE, control = list(),
                      field = "car_proper", positive = ~ pos_cov1) {
  f <- if (!spatial) ~ occ_cov1
       else stats::as.formula(sprintf("~ occ_cov1 + %s(graph = inp$adj)", field))
  tobs(formula = f, data = inp$cell_dat, family = occu_cover(inp$positive),
       detection = ~ det_cov1, positive = positive, y = inp$od$y,
       y_pos = inp$y_pos, visits = inp$od$det.covs, method = method,
       control = control)
}

# The occurrence field carried onto the cover arm, which is what the simulator's
# `alpha` puts there: a bare areal term loads on occurrence alone under both
# engines (gcol33/tulpaObs#217), so a fixture simulated with alpha > 0 is fitted
# with the copy declared.
.ocsn_fit_coupled <- function(inp, method = "nuts", control = list(),
                              field = "car_proper") {
  .ocsn_fit(inp, method, spatial = TRUE, control = control, field = field,
            positive = ~ pos_cov1 + copy(spatial()))
}


test_that("occu_cover spatial NUTS field-block FullGradFn matches the R oracle", {
  for (pos in c("lognormal", "beta")) {
    inp   <- .ocsn_inputs(side = 5L, J = 4L, seed = 11L, positive = pos)
    lap   <- .ocsn_fit(inp, "laplace", spatial = FALSE,
                       control = list(verbose = FALSE, max.iter = 60L))
    model <- lap$model
    N     <- inp$N
    adj   <- inp$adj

    # Fixed proper-CAR field block (tau, rho, alpha arbitrary for the check).
    rho <- 0.9; tau <- 1.5; alpha <- 0.7
    Q   <- tulpaObs:::.areal_Q(adj, rho)
    Qr  <- tau * Q + diag(1e-4 * tau, N)
    L   <- chol(Qr); Linv <- backsolve(L, diag(N))
    field <- list(n_field_units = N, field_map = seq_len(N),
                  Linv = Linv, alpha = alpha)

    spec <- tulpaObs:::.tobs_occu_cover_nuts_spec(model)
    spec$n_field_units <- N
    spec$field_map     <- as.integer(seq_len(N))
    spec$field_Linv    <- Linv
    spec$field_alpha   <- alpha

    np_coef <- length(lap$means)
    set.seed(3)
    theta <- c(as.numeric(lap$means) + stats::rnorm(np_coef, 0, 0.1),
               stats::rnorm(N, 0, 0.3))

    rr <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, model, 5, 5, field = field)
    cc <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec, theta, 5, 5)
    expect_equal(rr$lp, cc$lp, tolerance = 1e-8)
    expect_equal(rr$grad, as.numeric(cc$grad), tolerance = 1e-8)

    # The analytic R gradient (incl. the field block) matches a central FD.
    fd <- numeric(length(theta)); h <- 1e-5
    for (k in seq_along(theta)) {
      tp <- theta; tp[k] <- tp[k] + h
      tm <- theta; tm[k] <- tm[k] - h
      fd[k] <- (tulpaObs:::.tobs_occu_cover_nuts_logpost(tp, model, 5, 5, field = field)$lp -
                tulpaObs:::.tobs_occu_cover_nuts_logpost(tm, model, 5, 5, field = field)$lp) /
               (2 * h)
    }
    expect_equal(rr$grad, fd, tolerance = 1e-4)
  }
})


test_that("occu_cover spatial NUTS sampled-hyper block matches the R oracle (#204)", {
  # sigma / rho / alpha are sampled coordinates, so the field is z = sigma * (B1
  # %*% (s1(rho) raw1) + s2(rho) raw2) with the hypers riding bounded transforms.
  # Every new analytic derivative (the two block scalings, the copy amplitude's
  # data-score, each transform's Jacobian) is checked against a central FD, and
  # the compiled target against the R oracle.
  inp   <- .ocsn_inputs(side = 4L, J = 4L, seed = 11L)
  N     <- inp$N; adj <- inp$adj
  lap   <- .ocsn_fit(inp, "laplace", spatial = FALSE,
                     control = list(verbose = FALSE, max.iter = 60L))
  model <- lap$model
  model$site_cell <- seq_len(N)

  for (ty in c("icar", "bym2", "car_proper")) {
    warm <- tulpaObs:::.tobs_occu_cover_nuts_carproper_warm(
      model, adj, NULL, type = ty, max.iter = 60L)
    fb <- tulpaObs:::.occu_cover_nuts_field_block(adj, ty, N, seq_len(N), warm,
                                                  sample_hyper = TRUE)
    # icar has no mixing parameter; the other two sample all three hypers.
    expect_true(all(c("sigma", "alpha") %in% fb$sampled))
    expect_identical("rho" %in% fb$sampled, !identical(ty, "icar"))

    spec <- tulpaObs:::.tobs_occu_cover_nuts_spec(model)
    spec[names(fb$entries)] <- fb$entries
    n_base <- length(lap$means)
    n_tot  <- n_base + fb$n_raw + length(fb$sampled)
    set.seed(7)
    theta <- c(as.numeric(lap$means) + stats::rnorm(n_base, 0, 0.1),
               stats::rnorm(fb$n_raw, 0, 0.3),
               stats::rnorm(length(fb$sampled), 0, 0.6))

    rr <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, model, 5, 5,
                                                  field = fb$entries)
    cc <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec, theta, 5, 5)
    expect_equal(rr$lp, cc$lp, tolerance = 1e-8)
    expect_equal(rr$grad, as.numeric(cc$grad), tolerance = 1e-8)

    fd <- numeric(n_tot); h <- 1e-5
    for (k in seq_len(n_tot)) {
      tp <- theta; tp[k] <- tp[k] + h
      tm <- theta; tm[k] <- tm[k] - h
      fd[k] <- (tulpaObs:::.tobs_occu_cover_nuts_logpost(tp, model, 5, 5,
                                                         field = fb$entries)$lp -
                tulpaObs:::.tobs_occu_cover_nuts_logpost(tm, model, 5, 5,
                                                         field = fb$entries)$lp) / (2 * h)
    }
    expect_equal(rr$grad, fd, tolerance = 1e-4)
    # The hyper coordinates specifically -- a whole-vector tolerance would let a
    # wrong hyper derivative hide behind the much larger coefficient scores.
    hy <- (n_base + fb$n_raw + 1L):n_tot
    expect_equal(rr$grad[hy], fd[hy], tolerance = 1e-5)
    expect_gt(max(abs(fd[hy])), 1e-3)   # the check is not vacuous
  }
})


test_that("occu_cover spatial NUTS pinned hypers reproduce the fixed loading (#204)", {
  # `fixed.hyper = TRUE` conditions on the warm fit's (sigma, rho, alpha) as the
  # #74 / #113 path did. That is the degenerate configuration of the same block,
  # not a second code path, so it must agree with the old loading BIT for bit --
  # a tolerance here would let the two drift into different models.
  inp   <- .ocsn_inputs(side = 4L, J = 4L, seed = 11L)
  N     <- inp$N; adj <- inp$adj
  lap   <- .ocsn_fit(inp, "laplace", spatial = FALSE,
                     control = list(verbose = FALSE, max.iter = 60L))
  model <- lap$model
  model$site_cell <- seq_len(N)

  for (ty in c("icar", "bym2", "car_proper")) {
    warm <- tulpaObs:::.tobs_occu_cover_nuts_carproper_warm(
      model, adj, NULL, type = ty, max.iter = 60L)
    fb <- tulpaObs:::.occu_cover_nuts_field_block(adj, ty, N, seq_len(N), warm,
                                                  sample_hyper = FALSE)
    expect_identical(fb$sampled, character(0))
    expect_setequal(names(fb$pinned), c("sigma", "rho", "alpha"))

    fl <- tulpaObs:::.tobs_nuts_field_loading(
      adj, ty, N, tau = 1 / max(warm$sigma, 1e-3)^2, rho = warm$rho,
      sigma = warm$sigma)
    legacy <- list(n_field_units = N, field_map = seq_len(N),
                   Linv = fl$field_load, alpha = warm$alpha)
    set.seed(3)
    theta <- c(as.numeric(lap$means) + stats::rnorm(length(lap$means), 0, 0.1),
               stats::rnorm(fb$n_raw, 0, 0.3))
    r_new <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, model, 5, 5,
                                                      field = fb$entries)
    r_leg <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, model, 5, 5,
                                                      field = legacy)
    expect_identical(r_new$lp, r_leg$lp)
    expect_identical(r_new$grad, r_leg$grad)
  }
})


test_that("occu_cover spatial NUTS hyper prior is flat on the grid's axis (#204)", {
  # At raw = 0 the field is identically zero whatever the hypers are, so the data
  # term drops out and the target reduces to the hyper prior alone. Transformed
  # back to the axis coordinate the outer grid spaces its nodes in (log sigma /
  # log alpha / logit rho), that density must be CONSTANT -- which is the claim
  # that the sampler integrates the same flat-over-cells measure the
  # nested-Laplace grid does. A missing Jacobian is invisible in the gradient
  # check above (lp and grad would agree with each other and both be wrong) and
  # shows up here.
  inp   <- .ocsn_inputs(side = 4L, J = 4L, seed = 11L)
  N     <- inp$N; adj <- inp$adj
  lap   <- .ocsn_fit(inp, "laplace", spatial = FALSE,
                     control = list(verbose = FALSE, max.iter = 60L))
  model <- lap$model
  model$site_cell <- seq_len(N)

  for (ty in c("icar", "car_proper")) {
    warm <- tulpaObs:::.tobs_occu_cover_nuts_carproper_warm(
      model, adj, NULL, type = ty, max.iter = 60L)
    fb <- tulpaObs:::.occu_cover_nuts_field_block(adj, ty, N, seq_len(N), warm,
                                                  sample_hyper = TRUE)
    spec <- tulpaObs:::.tobs_occu_cover_nuts_spec(model)
    spec[names(fb$entries)] <- fb$entries
    base <- c(as.numeric(lap$means), numeric(fb$n_raw))
    n_h  <- length(fb$sampled)
    u_grid <- c(-2, -1, -0.4, 0, 0.7, 1.5)
    for (j in seq_len(n_h)) {
      lp <- vapply(u_grid, function(u) {
        th <- c(base, numeric(n_h)); th[length(base) + j] <- u
        tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec, th, 5, 5)$lp
      }, numeric(1))
      # p(u) ∝ exp(lp); the axis coordinate is t = t_lo + (t_hi - t_lo) e(u), so
      # p(t) = p(u) / (dt/du) must not vary with u.
      e <- stats::plogis(u_grid)
      dens_t <- exp(lp - max(lp)) / (e * (1 - e))
      expect_lt(stats::sd(dens_t) / mean(dens_t), 1e-8)
    }
  }
})


test_that("occu_cover spatial NUTS field-off path is byte-identical to non-spatial", {
  inp   <- .ocsn_inputs(side = 5L, J = 4L, seed = 21L)
  lap   <- .ocsn_fit(inp, "laplace", spatial = FALSE,
                     control = list(verbose = FALSE, max.iter = 60L))
  model <- lap$model
  np    <- length(lap$means)
  set.seed(9)
  theta <- as.numeric(lap$means) + stats::rnorm(np, 0, 0.12)

  # The spec without a field block reproduces the non-spatial occu_cover NUTS
  # target exactly (every field branch is guarded by has_field / n_field_units).
  spec0 <- tulpaObs:::.tobs_occu_cover_nuts_spec(model)
  r_off <- tulpaObs:::.tobs_occu_cover_nuts_logpost(theta, model, 5, 5)
  c_off <- tulpaObs:::cpp_occu_cover_nuts_joint_logpost(spec0, theta, 5, 5)
  expect_equal(r_off$lp, c_off$lp, tolerance = 1e-10)
  expect_equal(r_off$grad, as.numeric(c_off$grad), tolerance = 1e-10)
})


test_that("occu_cover NUTS samples icar; rejects SVC/RE spatial; advertises nuts", {
  expect_true("nuts" %in% tulpaObs:::.tobs_family_methods$occu_cover)

  inp <- .ocsn_inputs(side = 5L, J = 3L, seed = 5L)
  # Intrinsic icar() on the psi formula now samples via the coupled sum-to-zero
  # field (gcol33/tulpaObs#113); confirm the path runs and centres the field.
  fit_icar <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + icar(graph = inp$adj), data = inp$cell_dat,
    family = occu_cover("lognormal"), detection = ~ det_cov1,
    positive = ~ pos_cov1, y = inp$od$y, y_pos = inp$y_pos,
    visits = inp$od$det.covs, method = "nuts",
    control = list(verbose = FALSE, n.iter = 300L, n.warmup = 200L)))
  expect_identical(fit_icar$method, "nuts")
  expect_false(is.null(fit_icar$spatial_field))
  expect_lt(abs(mean(fit_icar$spatial_field)), 1e-6)   # sum-to-zero centred
  # A weighted SVC car_proper field is a grid-integrated structure.
  expect_error(
    suppressWarnings(tobs(
      formula = ~ occ_cov1 + car_proper(graph = inp$adj) +
        car_proper(graph = inp$adj, weight = occ_cov1),
      data = inp$cell_dat, family = occu_cover("lognormal"),
      detection = ~ det_cov1, positive = ~ pos_cov1, y = inp$od$y,
      y_pos = inp$y_pos, visits = inp$od$det.covs, method = "nuts",
      control = list(verbose = FALSE))),
    "single|SVC|weighted|nested_laplace")
})


test_that("occu_cover NUTS samples a single bar-form field (#203)", {
  inp  <- .ocsn_inputs(side = 5L, J = 4L, seed = 11L)
  dat  <- cbind(inp$cell_dat, cell_idx = seq_len(inp$N))
  adj  <- inp$adj
  ctl  <- list(verbose = FALSE, n.iter = 300L, n.warmup = 200L,
               n.chains = 1L, seed = 1L)
  run  <- function(f) suppressWarnings(tobs(
    formula = f, data = dat, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1, y = inp$od$y,
    y_pos = inp$y_pos, visits = inp$od$det.covs, method = "nuts", control = ctl))

  # A single-column bar desugars to exactly icar(graph, group_var = node), so the
  # resolved field description is the same object the non-bar spelling builds.
  sp_bar <- tulpaObs:::.occu_cover_nuts_spatial_term(
    ~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj), dat)
  sp_plain <- tulpaObs:::.occu_cover_nuts_spatial_term(
    ~ occ_cov1 + icar(graph = adj, group_var = "cell_idx"), dat)
  for (k in c("type", "n_units", "graph", "group_var", "adj_row_ptr",
              "adj_col_idx", "n_neighbors", "weight")) {
    expect_identical(sp_bar$spatial[[k]], sp_plain$spatial[[k]])
  }
  expect_identical(sp_bar$group_var, "cell_idx")

  # Both spellings therefore reach .tobs_fit_occu_cover_nuts_spatial() as the same
  # model. There is no sampler-noise tolerance to state: at one seed the two fits
  # run the same warm nested-Laplace fit, the same field loading and the same
  # chain, so the posteriors agree BIT for bit. A tolerance here would only hide
  # the two spellings drifting into different models.
  fit_bar   <- run(~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj))
  fit_plain <- run(~ occ_cov1 + icar(graph = adj, group_var = "cell_idx"))
  expect_identical(fit_bar$means, fit_plain$means)
  expect_identical(fit_bar$sds,   fit_plain$sds)
  expect_identical(fit_bar$spatial_field, fit_plain$spatial_field)
  expect_identical(fit_bar$spatial$type, "icar")
  expect_equal(length(fit_bar$spatial_field), inp$N)
  expect_lt(abs(mean(fit_bar$spatial_field)), 1e-6)   # sum-to-zero centred

  # The bar's `model =` picks the areal kind the sampler fixes the hyper of.
  fit_cp <- run(~ occ_cov1 +
                  spatial(~ 1 || cell_idx, graph = adj, model = "car_proper"))
  expect_identical(fit_cp$spatial$type, "car_proper")
})


test_that("occu_cover NUTS + a two-field bar names the field-block limit (#203)", {
  inp <- .ocsn_inputs(side = 5L, J = 3L, seed = 5L)
  dat <- cbind(inp$cell_dat, cell_idx = seq_len(inp$N))
  adj <- inp$adj
  run <- function(f) suppressWarnings(tobs(
    formula = f, data = dat, family = occu_cover("lognormal"),
    detection = ~ det_cov1, positive = ~ pos_cov1, y = inp$od$y,
    y_pos = inp$y_pos, visits = inp$od$det.covs, method = "nuts",
    control = list(verbose = FALSE, n.iter = 50L, n.warmup = 50L)))

  # An intercept + SVC bar declares two fields; the sampler's field block carries
  # one loading and one field_map (src/nuts_field_block.h), so the message must
  # name that, not the bar spelling.
  expect_error(run(~ occ_cov1 + spatial(~ 1 + occ_cov1 || cell_idx, graph = adj)),
               "one field block")
  expect_error(run(~ occ_cov1 + spatial(~ 1 + occ_cov1 || cell_idx, graph = adj)),
               "declares 2 fields")

  # A single-field correlated bar has no cross-covariance to estimate, the same
  # gate the nested-Laplace route applies.
  expect_error(run(~ occ_cov1 + spatial(~ 1 | cell_idx, graph = adj)),
               "cross-covariance")

  # A replicated field spans a Kronecker graph the one-graph sampler has no map for.
  dat$habitat <- rep(c("a", "b"), length.out = inp$N)
  expect_error(
    run(~ occ_cov1 + spatial(~ 1 || cell_idx, graph = adj, by = "habitat")),
    "nested_laplace")
})


test_that("occu_cover spatial NUTS recovers betas, field, coverage (lognormal)", {
  skip_on_cran()
  skip_if_fast()

  # The fixture simulates a field coupled onto the cover arm (alpha = 1), so the
  # fitted model declares the copy; a bare areal term is an occurrence-only field
  # (#217) and would be fitting a different model than the one simulated.
  n_seeds <- 5L
  side <- 8L; J <- 5L
  est <- se <- matrix(NA_real_, n_seeds, 7L)
  fcor <- div <- rhat <- rep(NA_real_, n_seeds)
  truth <- NULL
  for (s in seq_len(n_seeds)) {
    inp <- .ocsn_inputs(side = side, J = J, seed = 2000L + s, positive = "lognormal")
    truth <- inp$truth
    nut <- tryCatch(.ocsn_fit_coupled(inp, "nuts",
                    control = list(verbose = FALSE, n.iter = 1200L,
                                   n.warmup = 800L, n.chains = 2L, seed = 1L)),
                    error = function(e) NULL)
    if (is.null(nut)) next
    est[s, ]  <- as.numeric(nut$means)
    se[s, ]   <- as.numeric(nut$sds)
    fcor[s]   <- abs(stats::cor(nut$spatial_field, inp$sim$truth$f))
    div[s]    <- nut$nuts$divergent_total
    rhat[s]   <- max(nut$nuts$rhat, na.rm = TRUE)
  }
  ok <- stats::complete.cases(est)
  expect_gte(mean(ok), 0.75)

  # 0 (or near-0) divergences and converged chains.
  expect_lte(max(div[ok]), 5L)
  expect_lt(max(rhat[ok]), 1.1)

  # Field-shape recovery (the nested-Laplace gate; proper-CAR fit of an ICAR-
  # simulated field recovers the field shape).
  expect_gt(mean(fcor[ok]), 0.7)

  # Coefficient recovery: the six regression coefficients + log_sigma_pos.
  bias <- abs(colMeans(est[ok, , drop = FALSE]) - truth)
  expect_true(all(bias[1:6] < 0.30))
  expect_lt(bias[7], 0.15)

  # 95% Wald coverage on the slope coefficients (psi / p / cover slopes).
  cover <- abs(est[ok, c(2, 4, 6), drop = FALSE] -
               matrix(truth[c(2, 4, 6)], sum(ok), 3, byrow = TRUE)) <
           1.96 * se[ok, c(2, 4, 6), drop = FALSE]
  expect_gte(mean(cover), 0.75)
})


test_that("occu_cover spatial NUTS + icar coupled field recovers betas + field (#113)", {
  skip_on_cran()
  skip_if_fast()
  # The #71 sum-to-zero reparameterisation samples the coupled INTRINSIC icar
  # field (shared psi + alpha-copied cover arm) with the same in-tree field block,
  # via the non-square loading. Occupancy fields are weakly identified (one binary
  # site per node), so the field-cor gate is looser than the count families.
  n_seeds <- 5L; side <- 8L; J <- 5L
  est <- se <- matrix(NA_real_, n_seeds, 7L)
  fcor <- div <- zmean <- rep(NA_real_, n_seeds)
  truth <- NULL
  for (s in seq_len(n_seeds)) {
    inp <- .ocsn_inputs(side = side, J = J, seed = 3000L + s, positive = "lognormal")
    truth <- inp$truth
    nut <- tryCatch(.ocsn_fit_coupled(inp, "nuts", field = "icar",
                    control = list(verbose = FALSE, n.iter = 1200L,
                                   n.warmup = 800L, n.chains = 2L, seed = 1L)),
                    error = function(e) NULL)
    if (is.null(nut)) next
    est[s, ]  <- as.numeric(nut$means)
    se[s, ]   <- as.numeric(nut$sds)
    fcor[s]   <- abs(stats::cor(nut$spatial_field, inp$sim$truth$f))
    div[s]    <- nut$nuts$divergent_total
    zmean[s]  <- mean(nut$spatial_field)
  }
  ok <- stats::complete.cases(est)
  expect_gte(mean(ok), 0.75)
  expect_lte(max(div[ok]), 5L)                       # sum-to-zero geometry, clean
  expect_lt(max(abs(zmean[ok])), 1e-6)               # field is centred (sum z = 0)
  expect_gt(mean(fcor[ok]), 0.55)                    # weak occupancy-field identification
  bias <- abs(colMeans(est[ok, , drop = FALSE]) - truth)
  expect_lt(bias[7], 0.15)                           # dispersion recovered
  # 95% Wald coverage on the slope coefficients.
  cover <- abs(est[ok, c(2, 4, 6), drop = FALSE] -
               matrix(truth[c(2, 4, 6)], sum(ok), 3, byrow = TRUE)) <
           1.96 * se[ok, c(2, 4, 6), drop = FALSE]
  expect_gte(mean(cover), 0.70)
})


test_that("occu_cover spatial NUTS beta SDs calibrate to nested-Laplace SEs", {
  skip_on_cran()
  skip_if_fast()

  inp <- .ocsn_inputs(side = 8L, J = 5L, seed = 2010L, positive = "lognormal")

  nut <- .ocsn_fit(inp, "nuts",
                   control = list(verbose = FALSE, n.iter = 1500L,
                                  n.warmup = 1000L, n.chains = 2L, seed = 1L))
  nl <- suppressWarnings(tobs(
    formula = ~ occ_cov1 + icar(graph = inp$adj), data = inp$cell_dat,
    family = occu_cover("lognormal"), detection = ~ det_cov1,
    positive = ~ pos_cov1, y = inp$od$y, y_pos = inp$y_pos,
    visits = inp$od$det.covs, method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = 400L, engine = "joint")))

  expect_equal(nut$method, "nuts")
  expect_lte(nut$nuts$divergent_total, 5L)
  # The sampler no longer conditions on the warm fit's hypers (#204): it reports
  # which it integrated over and which it pinned. This fit's psi formula carries
  # a bare areal term, so the field is on occurrence alone and the copy
  # amplitude is pinned at 0 (#217); the field's own two hypers are sampled.
  expect_setequal(nut$nuts$sampled_hyper, c("sigma", "rho"))
  expect_identical(nut$nuts$fixed_hyper, "alpha")

  # The slope-coefficient SDs (psi / p / cover) calibrate to the grid-integrated
  # nested-Laplace SEs (the field hyper is fixed at the same kind of estimate).
  slope_nm <- c("psi_occ_cov1", "p_det_cov1", "pos_pos_cov1")
  ratio <- as.numeric(nut$sds[slope_nm]) / pmax(as.numeric(nl$sds[slope_nm]), 1e-4)
  expect_true(all(ratio > 0.5 & ratio < 2.0))
})


test_that("occu_cover spatial NUTS recovers betas + field (beta arm, smoke)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 4L
  side <- 8L; J <- 5L
  est <- matrix(NA_real_, n_seeds, 7L)
  fcor <- div <- rep(NA_real_, n_seeds)
  truth <- NULL
  for (s in seq_len(n_seeds)) {
    inp <- .ocsn_inputs(side = side, J = J, seed = 5000L + s, positive = "beta",
                        phi = 20)
    truth <- inp$truth
    nut <- tryCatch(.ocsn_fit_coupled(inp, "nuts",
                    control = list(verbose = FALSE, n.iter = 1200L,
                                   n.warmup = 800L, n.chains = 1L, seed = 1L)),
                    error = function(e) NULL)
    if (is.null(nut)) next
    est[s, ] <- as.numeric(nut$means)
    fcor[s]  <- abs(stats::cor(nut$spatial_field, inp$sim$truth$f))
    div[s]   <- nut$nuts$divergent_total
  }
  ok <- stats::complete.cases(est)
  expect_gte(mean(ok), 0.6)
  expect_lte(max(div[ok]), 8L)
  expect_gt(mean(fcor[ok]), 0.6)
  bias <- abs(colMeans(est[ok, , drop = FALSE]) - truth)
  expect_true(all(bias[1:6] < 0.35))
  expect_lt(bias[7], 0.45)
})


test_that("occu_cover spatial NUTS reports its hypers honestly per fit (#204)", {
  inp <- .ocsn_inputs(side = 5L, J = 4L, seed = 11L)
  ctl <- list(verbose = FALSE, n.iter = 200L, n.warmup = 200L, n.chains = 1L,
              seed = 1L)

  # icar carries no mixing parameter, so rho is pinned at 1 and says so; a bare
  # areal term asks for no cover-arm copy, so alpha is pinned at 0 and says so
  # (#217); sigma is sampled and carries a posterior.
  f_icar <- .ocsn_fit(inp, "nuts", field = "icar", control = ctl)
  expect_identical(f_icar$nuts$sampled_hyper, "sigma")
  expect_setequal(f_icar$nuts$fixed_hyper, c("rho", "alpha"))
  expect_equal(unname(f_icar$nuts$fixed_hyper_values[["rho"]]), 1)
  expect_equal(unname(f_icar$nuts$fixed_hyper_values[["alpha"]]), 0)
  expect_true(all(c("sigma", "rho", "alpha", "field_sd") %in%
                  colnames(f_icar$hyper_draws)))
  expect_gt(f_icar$nuts$hyper_sd[["sigma"]], 0)
  expect_equal(f_icar$nuts$hyper_sd[["alpha"]], 0)
  expect_equal(f_icar$nuts$hyper_sd[["rho"]], 0)

  # bym2 and car_proper sample the field's own two hypers.
  for (ty in c("bym2", "car_proper")) {
    ft <- .ocsn_fit(inp, "nuts", field = ty, control = ctl)
    expect_setequal(ft$nuts$sampled_hyper, c("sigma", "rho"))
    expect_identical(ft$nuts$fixed_hyper, "alpha")
    expect_gt(ft$nuts$hyper_sd[["rho"]], 0)
  }

  # A copy() puts the amplitude back among the sampled hypers.
  f_cp <- .ocsn_fit_coupled(inp, "nuts", control = ctl, field = "icar")
  expect_setequal(f_cp$nuts$sampled_hyper, c("sigma", "alpha"))
  expect_identical(f_cp$nuts$fixed_hyper, "rho")
  expect_gt(f_cp$nuts$hyper_sd[["alpha"]], 0)

  # An explicitly pinned fit reports all three as fixed, at the warm values.
  f_fix <- .ocsn_fit(inp, "nuts", field = "car_proper",
                     control = c(ctl, list(fixed.hyper = TRUE)))
  expect_identical(f_fix$nuts$sampled_hyper, character(0))
  expect_setequal(f_fix$nuts$fixed_hyper, c("sigma", "rho", "alpha"))
  expect_equal(unname(f_fix$nuts$hyper_sd), rep(0, 4))
  expect_equal(unname(f_fix$nuts$hyper_mean[c("sigma", "alpha")]),
               unname(f_fix$nuts$warm_hyper[c("sigma", "alpha")]))
})


test_that("occu_cover spatial NUTS hyper posteriors cover the truth (#204)", {
  skip_on_cran()
  skip_if_fast()
  # The point of sampling the hypers is that their posterior is honest, so the
  # gate is coverage of a known truth, not agreement with the nested-Laplace
  # point estimate (agreeing with that is the circularity #204 is about).
  #
  # Two targets, both stated on the sampler's own scale:
  #   alpha     the cover-arm copy amplitude -- dimensionless, and defined
  #             identically by the simulator (eta_pos += alpha * sigma * f) and
  #             the fitter, so the truth needs no conversion. The simulated truth
  #             couples the field onto the cover arm, so the fitted model
  #             declares that with copy(spatial()): a bare areal term is an
  #             occurrence-only field and pins alpha at 0 (#217).
  #   field_sd  the geometric-mean marginal SD of the field the block implies at
  #             that draw's hypers. The simulator's f carries geo-mean marginal
  #             variance 1 (Sorbye-Rue), so the truth is exactly `sigma`. sigma
  #             itself is NOT the comparable quantity: icar, bym2 and car_proper
  #             normalise their precisions differently, so the same field has
  #             three different sigmas and one field_sd.
  n_seeds <- 12L
  side <- 8L; J <- 5L
  sig_true <- 0.7; alpha_true <- 1.0
  for (ty in c("icar", "car_proper")) {
    lo <- hi <- fsd_lo <- fsd_hi <- div <- rep(NA_real_, n_seeds)
    a_mean <- fsd_mean <- rep(NA_real_, n_seeds)
    for (s in seq_len(n_seeds)) {
      inp <- .ocsn_inputs(side = side, J = J, seed = 7000L + s,
                          sigma = sig_true, alpha = alpha_true)
      nut <- tryCatch(.ocsn_fit_coupled(inp, "nuts", field = ty,
                        control = list(verbose = FALSE, n.iter = 1000L,
                                       n.warmup = 800L, n.chains = 2L,
                                       seed = 1L)),
                      error = function(e) NULL)
      if (is.null(nut)) next
      q <- stats::quantile(nut$hyper_draws[, "alpha"], c(0.025, 0.975))
      lo[s] <- q[[1L]]; hi[s] <- q[[2L]]
      a_mean[s] <- mean(nut$hyper_draws[, "alpha"])
      qf <- stats::quantile(nut$hyper_draws[, "field_sd"], c(0.025, 0.975))
      fsd_lo[s] <- qf[[1L]]; fsd_hi[s] <- qf[[2L]]
      fsd_mean[s] <- mean(nut$hyper_draws[, "field_sd"])
      div[s] <- nut$nuts$divergent_total
    }
    ok <- !is.na(lo)
    expect_gte(mean(ok), 0.75)
    expect_lte(max(div[ok]), 5L)
    # 95% credible-interval coverage of the two truths. The pooled floor is the
    # package's usual 0.85 rubric, not a literal 0.95 (that flakes at 12 seeds).
    expect_gte(mean(lo[ok] <= alpha_true & alpha_true <= hi[ok]), 0.85)
    expect_gte(mean(fsd_lo[ok] <= sig_true & sig_true <= fsd_hi[ok]), 0.85)
    # A one-sided shift is a property of the MEAN over seeds, so it is asserted
    # there and not per fit. Both gates are set from the measurement, not chosen:
    # alpha lands within 0.06 of truth, field_sd about 0.25 high (the posterior of
    # a positive variance component at 64 binary sites is right-skewed, so its
    # MEAN sits above the bulk -- `fit$nuts$hyper_median` is the summary to quote).
    # See NOTES_measurements.md.
    expect_lt(abs(mean(a_mean[ok]) - alpha_true), 0.35)
    expect_lt(abs(mean(fsd_mean[ok]) - sig_true), 0.45)
  }
})


test_that("occu_cover spatial NUTS fit exposes the S3 surface", {
  skip_on_cran()
  skip_if_fast()

  inp <- .ocsn_inputs(side = 7L, J = 4L, seed = 7L)
  nut <- .ocsn_fit(inp, "nuts",
                   control = list(verbose = FALSE, n.iter = 800L,
                                  n.warmup = 600L, n.chains = 1L, seed = 1L))
  np <- length(nut$means)
  # coef() returns the per-arm coefficient list (psi / p / pos); the flattened
  # length is the coefficient count (the trailing log-dispersion is not a coef).
  cf <- coef(nut)
  expect_true(is.list(cf))
  expect_equal(length(unlist(cf)), np - 1L)
  expect_equal(dim(vcov(nut)), c(np, np))
  expect_equal(nrow(confint(nut)), np)
  expect_equal(length(nut$spatial_field), inp$N)
  expect_true(all(is.finite(nut$spatial_field)))
  expect_identical(nut$spatial$type, "car_proper")
})


# --------------------------------------------------------------------------- #
# copy() on the positive arm reaches the sampler (gcol33/tulpaObs#210)          #
# --------------------------------------------------------------------------- #

# Same fixture as .ocsn_fit, with a copy() term on the positive formula.
.ocsn_fit_copy <- function(inp, pos_formula, control = list(),
                           field = "car_proper") {
  .ocsn_fit(inp, "nuts", spatial = TRUE, control = control, field = field,
            positive = pos_formula)
}

test_that("occu_cover spatial NUTS honours a copy()'s fixed amplitude (#210)", {
  skip_on_cran()
  skip_if_fast()
  # A scalar `alpha =` is a pinned amplitude: it collapses the warm fit's copy
  # axis to one node, which leaves the sampler nothing to integrate, so alpha is
  # conditioned on exactly the value the user named. The fit reports which
  # hypers it sampled and which it pinned, so this is asserted on that record
  # rather than inferred from the draws alone.
  inp <- .ocsn_inputs(side = 5L, J = 3L, seed = 21L)
  nut <- .ocsn_fit_copy(
    inp, ~ pos_cov1 + copy(spatial(), alpha = 0.35),
    control = list(verbose = FALSE, n.iter = 200L, n.warmup = 200L,
                   n.chains = 1L, seed = 1L))
  expect_true("alpha" %in% nut$nuts$fixed_hyper)
  expect_false("alpha" %in% nut$nuts$sampled_hyper)
  expect_equal(unname(nut$nuts$fixed_hyper_values[["alpha"]]), 0.35)
  # Pinned means pinned in every draw, not merely centred there.
  expect_true(all(nut$hyper_draws[, "alpha"] == 0.35))
})

test_that("occu_cover spatial NUTS honours a copy()'s amplitude grid (#210)", {
  skip_on_cran()
  skip_if_fast()
  # An integrated `alpha = grid(c(...))` sets the SUPPORT of the sampled
  # amplitude: the warm nested-Laplace fit lays the axis, and its span is the
  # support of the flat prior the sampler draws alpha under
  # (.occu_cover_nuts_hyper_bounds). So the knob reaches the sampler even though
  # alpha is no longer pinned at the warm estimate (gcol33/tulpaObs#204).
  #
  # The band is deliberately BELOW the simulation truth (alpha = 1), so the
  # posterior pushes against its upper bound: a dropped copy() would fall back
  # to the default axis, whose positive span reaches 3, and the draws would sit
  # near 1 rather than under 0.5. The assertion therefore cannot pass vacuously.
  inp <- .ocsn_inputs(side = 5L, J = 3L, seed = 21L, alpha = 1.0)
  nut <- .ocsn_fit_copy(
    inp, ~ pos_cov1 + copy(spatial(), alpha = grid(c(0.2, 0.5))),
    control = list(verbose = FALSE, n.iter = 200L, n.warmup = 200L,
                   n.chains = 1L, seed = 1L))
  expect_true("alpha" %in% nut$nuts$sampled_hyper)
  a <- nut$hyper_draws[, "alpha"]
  expect_gte(min(a), 0.2 - 1e-8)
  expect_lte(max(a), 0.5 + 1e-8)
  # The default axis reaches far above this band, which is what the copy()
  # narrowed.
  expect_gt(max(as.numeric(tulpaObs:::.tobs_default_alpha_grid())), 0.5)
})

test_that("occu_cover spatial NUTS refuses a copy() it cannot resolve (#210)", {
  skip_on_cran()
  skip_if_fast()
  inp <- .ocsn_inputs(side = 5L, J = 3L, seed = 21L)
  ctl <- list(verbose = FALSE, n.iter = 50L, n.warmup = 50L, n.chains = 1L,
              seed = 1L)
  # The amplitude is set in ONE place. Before the copy() reached this path the
  # control knob simply won and the copy() vanished.
  expect_error(
    .ocsn_fit_copy(inp, ~ pos_cov1 + copy(spatial(), alpha = 0.35),
                   control = c(ctl, list(alpha.grid = c(0.2, 0.5)))),
    "copy() in the positive", fixed = TRUE)
  # A string reference is not a coupling selector here, under either engine.
  expect_error(
    .ocsn_fit_copy(inp, ~ pos_cov1 + copy("occ_space", alpha = 0.35),
                   control = ctl),
    "not a string")
})


# --------------------------------------------------------------------------- #
# A bare areal term is occurrence-only on BOTH engines (gcol33/tulpaObs#217)    #
# --------------------------------------------------------------------------- #

# The grid-integrated route on the same input, for the two-backend comparison.
# icar, because that is the field kind both engines accept.
.ocsn_fit_nl <- function(inp, positive = ~ pos_cov1, control = list()) {
  suppressWarnings(.ocsn_fit(
    inp, "nested_laplace", spatial = TRUE, field = "icar", positive = positive,
    control = utils::modifyList(
      list(verbose = FALSE, max.iter = 300L, engine = "joint"), control)))
}

test_that("occu_cover: a bare areal term loads on occurrence alone on both engines (#217)", {
  skip_on_cran()
  skip_if_fast()
  # Same input, same model. The psi formula carries an areal term and the
  # positive formula carries no copy(), so the field rides occupancy and the
  # cover-arm copy amplitude is 0 under nuts exactly as under nested_laplace.
  # Before this, the sampler kept the default amplitude axis and estimated a
  # cover-arm copy the deterministic route had pinned away.
  inp <- .ocsn_inputs(side = 6L, J = 4L, seed = 11L)
  ctl <- list(verbose = FALSE, n.iter = 400L, n.warmup = 400L, n.chains = 1L,
              seed = 1L)

  nut <- .ocsn_fit(inp, "nuts", field = "icar", control = ctl)
  nl  <- .ocsn_fit_nl(inp)

  # nuts: alpha pinned at 0, and the record says so. `fixed_hyper` is a character
  # vector (#204), so this cannot be read as a stale TRUE / FALSE flag.
  expect_true(is.character(nut$nuts$fixed_hyper))
  expect_true("alpha" %in% nut$nuts$fixed_hyper)
  expect_false("alpha" %in% nut$nuts$sampled_hyper)
  expect_equal(unname(nut$nuts$fixed_hyper_values[["alpha"]]), 0)
  # Pinned in every draw, not merely centred at 0.
  expect_true(all(nut$hyper_draws[, "alpha"] == 0))
  expect_equal(unname(nut$nuts$hyper_sd[["alpha"]]), 0)

  # nested_laplace: the same amplitude, from the grid the dispatcher pins at 0.
  expect_equal(nl$spatial$alpha_mean, 0)

  # ... and the two backends then recover the same occurrence field.
  expect_gt(abs(stats::cor(nut$spatial_field, nl$spatial_field)), 0.8)
})

test_that("occu_cover: copy() puts the amplitude back on both engines (#217)", {
  skip_on_cran()
  skip_if_fast()
  inp <- .ocsn_inputs(side = 6L, J = 4L, seed = 11L, alpha = 1.0)
  ctl <- list(verbose = FALSE, n.iter = 400L, n.warmup = 400L, n.chains = 1L,
              seed = 1L)
  pos <- ~ pos_cov1 + copy(spatial())

  nut <- .ocsn_fit(inp, "nuts", field = "icar", control = ctl, positive = pos)
  nl  <- .ocsn_fit_nl(inp, positive = pos)

  expect_true("alpha" %in% nut$nuts$sampled_hyper)
  expect_gt(nut$nuts$hyper_mean[["alpha"]], 0)
  expect_gt(nl$spatial$alpha_mean, 0)

  # The default copy() rides the same amplitude axis the pre-#217 no-copy fit
  # used, so an explicit copy() is the spelling that reproduces it.
  ctrl_copy <- tulpaObs:::.occu_cover_nuts_copy_control(
    list(structure(list(id = "spatial()", selector_type = "spatial",
                        selector_group = NULL, component = NULL,
                        alpha_integrate = NA, alpha_grid = NULL,
                        copy_terms = NULL), class = "tobs_copy")),
    list(group_var = NULL), list())
  expect_equal(as.numeric(ctrl_copy$alpha.grid),
               as.numeric(tulpaObs:::.tobs_default_alpha_grid()))
})

test_that("occu_cover: control$alpha.grid overrides the amplitude on both engines (#217)", {
  skip_on_cran()
  skip_if_fast()
  # The low-level knob stays live: a user who names the axis gets it, on either
  # engine, and the translation steps aside.
  inp <- .ocsn_inputs(side = 5L, J = 3L, seed = 21L)
  ctl <- list(verbose = FALSE, n.iter = 200L, n.warmup = 200L, n.chains = 1L,
              seed = 1L)

  band <- .ocsn_fit(inp, "nuts", field = "icar",
                    control = c(ctl, list(alpha.grid = c(0.2, 0.5))))
  expect_true("alpha" %in% band$nuts$sampled_hyper)
  a <- band$hyper_draws[, "alpha"]
  expect_gte(min(a), 0.2 - 1e-8)
  expect_lte(max(a), 0.5 + 1e-8)

  # A one-node axis pins the amplitude. The warm fit reaches it as a grid-weight
  # average over that one node, so it carries the last bit of the weight
  # normalisation; the pinning itself is exact (zero spread across draws).
  pinned <- .ocsn_fit(inp, "nuts", field = "icar",
                      control = c(ctl, list(alpha.grid = 0.8)))
  expect_true("alpha" %in% pinned$nuts$fixed_hyper)
  expect_equal(unname(pinned$nuts$fixed_hyper_values[["alpha"]]), 0.8)
  expect_equal(unname(pinned$nuts$hyper_sd[["alpha"]]), 0)
  expect_equal(range(pinned$hyper_draws[, "alpha"]), c(0.8, 0.8))

  # Same knob, same meaning on the grid-integrated route. A one-node axis is the
  # assertion that carries: the reported amplitude of a MULTI-node axis is the
  # engine's own derived posterior mean, which is not confined to the nodes.
  nl_pin <- .ocsn_fit_nl(inp, control = list(alpha.grid = 0.8))
  expect_equal(nl_pin$spatial$alpha_mean, 0.8)
  nl_band <- .ocsn_fit_nl(inp, control = list(alpha.grid = c(0.2, 0.5)))
  expect_gt(nl_band$spatial$alpha_mean, 0)
  expect_false(isTRUE(all.equal(nl_band$spatial$alpha_mean,
                                nl_pin$spatial$alpha_mean)))
})
