# temporal() on the occu_cover() occurrence formula, joint nested-Laplace route:
# one psi-arm temporal block (RW1 / RW2 / AR1 / iid over the term's time levels),
# reported, drawn and predicted with the rest of the joint posterior.

.oct_data <- function(n_cells = 10L, n_int = 5L, J = 6L, psi_t = NULL,
                      seed = 2L) {
  set.seed(seed)
  if (is.null(psi_t)) psi_t <- 0.6 * sin(seq_len(n_int))
  adj <- matrix(0L, n_cells, n_cells)
  for (i in seq_len(n_cells - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  rows <- vector("list", n_cells * n_int)
  k <- 0L
  for (cell in seq_len(n_cells)) for (ti in seq_len(n_int)) {
    psi <- stats::plogis(-0.3 + psi_t[ti] + 0.15 * (cell - n_cells / 2))
    z  <- stats::rbinom(1L, 1L, psi)
    x1 <- stats::rnorm(J)
    p  <- stats::plogis(0.2 + 0.4 * x1)
    occur <- if (z == 1L) stats::rbinom(J, 1L, p) else rep(0L, J)
    cover <- ifelse(occur == 1L,
                    pmin(pmax(stats::rbeta(J, 2, 4), 1e-4), 1 - 1e-4), 0)
    k <- k + 1L
    rows[[k]] <- data.frame(
      cell_idx = cell, yr.interval = ti, site_key = paste(cell, ti, sep = "_"),
      visit = seq_len(J), time.sc = ti - (n_int + 1) / 2, x1 = x1,
      grp = (cell %% 3L) + 1L, occur = as.integer(occur), cover.flat = cover)
  }
  list(dd = do.call(rbind, rows), adj = adj, psi_t = psi_t)
}

.oct_ctrl <- list(engine = "joint", n.threads = 1L, n.threads.outer = 1L,
                  max.iter = 200L, sigma.grid = c(0.5, 1.5),
                  phi.grid.pos = c(2, 10), integration = "grid",
                  verbose = FALSE, progress = FALSE)

.oct_fit <- function(sim, occ, det = ~ x1, pos = ~ time.sc,
                     method = "nested_laplace", control = .oct_ctrl,
                     family = occu_cover(response = "beta",
                                         cover.aggregate = "none")) {
  dd  <- sim$dd
  adj <- sim$adj
  environment(occ) <- environment()
  tobs(occurrence = occ, data = dd, family = family, detection = det,
       positive = pos, response = "occur", y.pos = "cover.flat",
       site = "site_key", visit = "visit",
       det.covs = c("x1", "time.sc", "yr.interval"),
       occ.covs = c("cell_idx", "time.sc", "yr.interval", "grp"),
       method = method, control = control)
}

test_that("an rw1 temporal() term on psi joins the joint fit as its own block", {
  skip_on_cran()
  sim <- .oct_data()
  fit <- .oct_fit(sim, ~ temporal(yr.interval, type = "rw1") +
                    spatial(~ 1 || cell_idx, graph = adj))
  expect_true("sigma_temporal" %in% names(fit$means))
  expect_true("b2.tau" %in% colnames(fit$joint_fit$theta_grid))
  te <- fit$re$temporal
  expect_identical(te$kind, "temporal")
  expect_identical(te$arm, "psi")
  expect_identical(te$var, "yr.interval")
  expect_equal(te$levels, as.character(1:5))
  expect_length(te$blup, 5L)
  expect_true(all(te$blup_sd > 0))
  # The level structure of the simulated sin() trajectory: up early, down late.
  expect_gt(cor(te$blup, sim$psi_t), 0.5)
  rf <- ranef(fit)
  expect_true(all(c("1", "5") %in% as.character(rf$level[rf$group == "yr.interval"])))
})

test_that("rw2, ar1 and iid temporal() terms carry their own axes", {
  skip_on_cran()
  sim <- .oct_data()
  axes <- list(rw2 = "b2.tau", ar1 = c("b2.tau", "b2.rho"), iid = "b2.sigma")
  for (ty in names(axes)) {
    occ <- stats::as.formula(sprintf(paste0(
      "~ temporal(yr.interval, type = '%s') + ",
      "spatial(~ 1 || cell_idx, graph = adj)"), ty))
    fit <- .oct_fit(sim, occ)
    expect_true(all(axes[[ty]] %in% colnames(fit$joint_fit$theta_grid)),
                info = ty)
    expect_true("sigma_temporal" %in% names(fit$means), info = ty)
    if (ty == "ar1") expect_true("rho_temporal" %in% names(fit$means))
    # The iid type rides tulpa's non-centred iid block, whose effect is the
    # cell's SD node times the latent; the RW / AR1 latent is the effect.
    expect_identical(fit$re$temporal$scale_axes,
                     if (ty == "iid") "b2.sigma" else NA_character_, info = ty)
  }
})

test_that("temporal() composes with a psi random intercept", {
  skip_on_cran()
  sim <- .oct_data()
  fit <- .oct_fit(sim, ~ temporal(yr.interval, type = "rw1") + (1 | grp) +
                    spatial(~ 1 || cell_idx, graph = adj))
  expect_setequal(names(fit$re), c("psi", "temporal"))
  expect_true(all(c("sigma_re", "sigma_temporal") %in% names(fit$means)))
  expect_true(is.finite(suppressWarnings(waic(fit))$estimates["waic", "Estimate"]))
})

test_that("a non-linear psi trajectory is recovered by rw1 better than by a line", {
  skip_on_cran()
  # Up then down: a linear time slope cannot follow it.
  psi_t <- c(-1.2, 0.2, 1.0, 0.2, -1.2)
  sim <- .oct_data(n_cells = 16L, n_int = 5L, J = 8L, psi_t = psi_t, seed = 11L)
  f_rw <- .oct_fit(sim, ~ temporal(yr.interval, type = "rw1") +
                     spatial(~ 1 || cell_idx, graph = adj))
  f_ln <- .oct_fit(sim, ~ time.sc + spatial(~ 1 || cell_idx, graph = adj))
  nd <- data.frame(cell = 1:16, cell_idx = 1:16, time.sc = 0, yr.interval = 3)
  set.seed(1)
  traj <- function(f) predict(f, newdata = nd, type = "trajectory",
                              times = (1:5) - 3, time.col = "time.sc",
                              nsim = 400L, draws = FALSE, aggregate = TRUE)
  rw <- traj(f_rw); ln <- traj(f_ln)
  truth <- vapply(seq_along(psi_t), function(t)
    mean(stats::plogis(-0.3 + psi_t[t] + 0.15 * (1:16 - 8))), numeric(1))
  err_rw <- sqrt(mean((rw$psi_mean - truth)^2))
  err_ln <- sqrt(mean((ln$psi_mean - truth)^2))
  expect_lt(err_rw, err_ln)
  # The middle interval is the peak of the fitted trajectory as of the truth.
  expect_identical(which.max(rw$psi_mean), 3L)
  # The per-interval effect is on the logit scale of psi itself (the RW latent
  # enters the predictor unscaled), so its centred levels track the simulated
  # deviations in size, not only in shape.
  b  <- f_rw$re$temporal$blup
  tc <- psi_t - mean(psi_t)
  slope <- sum(b * tc) / sum(tc^2)
  expect_gt(slope, 0.5)
  expect_lt(slope, 1.5)
  w_rw <- suppressWarnings(waic(f_rw))$estimates["waic", "Estimate"]
  w_ln <- suppressWarnings(waic(f_ln))$estimates["waic", "Estimate"]
  expect_lt(w_rw, w_ln)
})

test_that("predict() moves psi with the fitted temporal effect", {
  skip_on_cran()
  sim <- .oct_data()
  fit <- .oct_fit(sim, ~ temporal(yr.interval, type = "rw1") +
                    spatial(~ 1 || cell_idx, graph = adj))
  nd <- data.frame(cell = 1:10, cell_idx = 1:10, time.sc = 0, yr.interval = 3)
  set.seed(3)
  by_sc <- predict(fit, newdata = nd, type = "trajectory", times = (1:5) - 3,
                   time.col = "time.sc", nsim = 2000L, draws = FALSE,
                   aggregate = TRUE)
  set.seed(3)
  by_int <- predict(fit, newdata = nd, type = "trajectory", times = 1:5,
                    time.col = "yr.interval", nsim = 2000L, draws = FALSE,
                    aggregate = TRUE)
  # The two time columns pair one-to-one, so both spellings read the same level.
  expect_equal(by_sc$psi_mean, by_int$psi_mean, tolerance = 0.03)
  # Across intervals psi follows the fitted per-level effect, on the logit scale.
  lg <- stats::qlogis(by_int$psi_median)
  expect_gt(cor(diff(lg), diff(fit$re$temporal$blup)), 0.9)
  expect_error(predict(fit, newdata = nd, type = "trajectory",
                       times = c(-2, 2.5), time.col = "time.sc", nsim = 50L),
               "not a fitted time")
  expect_error(predict(fit, newdata = data.frame(cell = 1:3, time.sc = 0),
                       type = "occurrence", nsim = 50L),
               "needs that column")
  expect_error(predict(fit, newdata = transform(nd, yr.interval = 9),
                       type = "occurrence", nsim = 50L),
               "not a fitted level")
})

test_that("temporal() where the route cannot carry it is refused by name", {
  skip_on_cran()
  sim <- .oct_data()
  sp <- "spatial(~ 1 || cell_idx, graph = adj)"
  f  <- function(x) stats::as.formula(paste("~", x))
  tt <- "temporal(yr.interval, type = 'rw1')"
  expect_error(.oct_fit(sim, f(paste(tt, "+", sp)), det = f(paste("x1 +", tt))),
               "temporal")
  expect_error(.oct_fit(sim, f(paste(tt, "+", sp)),
                        pos = f(paste("time.sc +", tt))), "temporal")
  expect_error(.oct_fit(sim, f(tt), method = "laplace", control = list()),
               "temporal")
  expect_error(.oct_fit(sim, f(tt)), "requires a spatial term")
  expect_error(.oct_fit(sim, f(paste(tt, "+", sp)), method = "nuts",
                        control = list(n.iter = 10L)), "temporal")
  expect_error(.oct_fit(sim, f(paste(
    "temporal(yr.interval, type = 'rw1', group = grp) +", sp))),
    "temporal\\(group = \\)")
  expect_error(.oct_fit(sim, f(paste(
    "temporal(yr.interval, type = 'rw1', tau.shape = 2) +", sp))),
    "tau.shape")
  expect_error(.oct_fit(sim, f(paste(tt, "+ temporal(time.sc, type = 'rw1') +",
                                     sp))), "2 temporal\\(\\) terms")
  expect_error(.oct_fit(sim, f(paste(
    tt, "+ spatial(~ 1 + time.sc | cell_idx, graph = adj)"))),
    "does not compose with a temporal")
})
