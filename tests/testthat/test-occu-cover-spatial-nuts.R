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
                      field = "car_proper") {
  f <- if (!spatial) ~ occ_cov1
       else stats::as.formula(sprintf("~ occ_cov1 + %s(graph = inp$adj)", field))
  tobs(formula = f, data = inp$cell_dat, family = occu_cover(inp$positive),
       detection = ~ det_cov1, positive = ~ pos_cov1, y = inp$od$y,
       y_pos = inp$y_pos, visits = inp$od$det.covs, method = method,
       control = control)
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


test_that("occu_cover spatial NUTS recovers betas, field, coverage (lognormal)", {
  skip_on_cran()
  skip_if_fast()

  n_seeds <- 5L
  side <- 8L; J <- 5L
  est <- se <- matrix(NA_real_, n_seeds, 7L)
  fcor <- div <- rhat <- rep(NA_real_, n_seeds)
  truth <- NULL
  for (s in seq_len(n_seeds)) {
    inp <- .ocsn_inputs(side = side, J = J, seed = 2000L + s, positive = "lognormal")
    truth <- inp$truth
    nut <- tryCatch(.ocsn_fit(inp, "nuts",
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
    nut <- tryCatch(.ocsn_fit(inp, "nuts", field = "icar",
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
  expect_true(isTRUE(nut$nuts$fixed_hyper))

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
    nut <- tryCatch(.ocsn_fit(inp, "nuts",
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
