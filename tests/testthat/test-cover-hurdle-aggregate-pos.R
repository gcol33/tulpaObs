# =============================================================================
# test-cover-hurdle-aggregate-pos.R - control$aggregate.pos: exact grouped-beta
# sufficient-statistic reduction of the positive (cover) arm (tulpaObs#49).
#
# Beta is not a count family, so there is no single-row collapse. Plots sharing
# the positive design row AND every per-observation latent component (cell, trend
# weight) are exchangeable Beta(mu*phi, (1-mu)*phi) draws; the beta log-density
# is linear in log(y) and log(1-y), so one row carrying (n, sum log y,
# sum log(1-y)) reproduces the log-likelihood, gradient and Fisher Hessian
# exactly. tulpa's built-in beta spec reads those sufficient statistics
# (slog_y / slog_1my on the arm). Default ON for the beta arm; the recovery
# sign-off behind the default lives in test-cover-hurdle-aggregate-recovery.R.
# =============================================================================


.aop_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) { if (s > 1L) adj[s, s - 1L] <- 1L; if (s < n) adj[s, s + 1L] <- 1L }
  adj
}

.aop_icar_f <- function(adj) {
  Q  <- tulpaObs:::.occu_cover_icar_Q(adj)
  sc <- tulpaObs:::.occu_cover_icar_scale(adj)
  eig <- eigen(Q, symmetric = TRUE); keep <- eig$values > 1e-8
  fk <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                     (stats::rnorm(sum(keep)) / sqrt(eig$values[keep])))
  (fk - mean(fk)) / sqrt(sc)
}

# Cell-level positive design (x, time at cell level), so the detected plots
# within a cell share the full positive linear predictor and the grouped collapse
# actually fires.
.aop_sim <- function(seed, n_s = 16L, n_per = 12L, trend = TRUE) {
  set.seed(seed); adj <- .aop_chain_adj(n_s)
  f1 <- .aop_icar_f(adj); f2 <- .aop_icar_f(adj)
  cell <- rep(seq_len(n_s), each = n_per); N <- length(cell)
  x    <- as.numeric(scale(stats::rnorm(n_s)))[cell]
  time <- as.numeric(scale(stats::rnorm(n_s)))[cell]
  eta_o <- -0.3 + 0.7 * x + 0.8 * f1[cell] + if (trend) 0.6 * f2[cell] * time else 0
  occur <- stats::rbinom(N, 1L, stats::plogis(eta_o))
  eta_p <- 0.2 - 0.5 * x + 0.8 * f1[cell] + if (trend) 0.9 * 0.6 * f2[cell] * time else 0
  mu <- stats::plogis(eta_p); y <- numeric(N); pos <- occur == 1L
  y[pos] <- pmin(pmax(stats::rbeta(sum(pos), mu[pos] * 18, (1 - mu[pos]) * 18), 1e-6), 1 - 1e-6)
  list(data = data.frame(x = x, time = time, region = factor(cell)), y = y, adj = adj)
}

.aop_fit <- function(s, trend, agg.occ = FALSE, agg.pos = FALSE) {
  ctrl <- list(verbose = FALSE, sigma.grid = c(0.5, 0.8, 1.2), rho.grid = 0.5,
               phi.grid = c(8, 18, 40), adaptive.grid = FALSE, max.iter = 300L,
               aggregate.occ = agg.occ, aggregate.pos = agg.pos)
  if (trend) {
    ctrl$trend <- list(weight = "time")
    ctrl$alpha.grid <- c(0, 0.5, 1.0); ctrl$alpha.grid.trend <- c(0, 0.5, 1.0)
  } else {
    ctrl$sigma.pos.grid <- c(0.4, 0.8)
  }
  suppressWarnings(tobs(
    formula = ~ x + bym2(graph = s$adj, group_var = "region"),
    data = s$data, family = cover("beta"), y = s$y,
    method = "nested_laplace", control = ctrl))
}

.aop_expect_identical_fits <- function(fa, fb, tol = 1e-7) {
  expect_true(fa$converged && fb$converged)
  expect_equal(fa$beta_occ,     fb$beta_occ,     tolerance = tol)
  expect_equal(fa$beta_pos,     fb$beta_pos,     tolerance = tol)
  expect_equal(fa$se_occ,       fb$se_occ,       tolerance = tol)
  expect_equal(fa$se_pos,       fb$se_pos,       tolerance = tol)
  expect_equal(fa$phi_pos,      fb$phi_pos,      tolerance = tol)
  expect_equal(fa$log_marginal, fb$log_marginal, tolerance = tol)
}


test_that(".cover_aggregate_pos collapses exchangeable beta rows to exact sufficient stats", {
  # Three groups by (X row, idx, weight); the beta sufficient statistics are
  # the per-group counts and sums of log(y) / log(1-y).
  y   <- c(0.2, 0.5, 0.8, 0.3, 0.6, 0.1)
  X   <- rbind(c(1, 2), c(1, 2), c(1, 2), c(3, 4), c(3, 4), c(1, 2))
  idx <- c(5L, 5L, 5L, 7L, 7L, 6L)

  ag <- tulpaObs:::.cover_aggregate_pos(y, X, list(idx = idx))

  # Groups (group-id order): (X=1,2;idx=5) rows 1-3, (X=1,2;idx=6) row 6,
  # (X=3,4;idx=7) rows 4-5.
  expect_equal(ag$n, c(3L, 1L, 2L))
  expect_equal(ag$X, rbind(c(1, 2), c(1, 2), c(3, 4)))
  expect_equal(ag$keys$idx, c(5L, 6L, 7L))
  expect_equal(ag$slog_y,   c(sum(log(y[1:3])),       log(y[6]),       sum(log(y[4:5]))))
  expect_equal(ag$slog_1my, c(sum(log1p(-y[1:3])), log1p(-y[6]), sum(log1p(-y[4:5]))))

  # Sufficiency invariants: counts and the two log-sums are conserved.
  expect_equal(sum(ag$n), length(y))
  expect_equal(sum(ag$slog_y),   sum(log(y)))
  expect_equal(sum(ag$slog_1my), sum(log1p(-y)))

  # Same (X, idx) but a different weight must NOT merge.
  ag2 <- tulpaObs:::.cover_aggregate_pos(
    y = c(0.2, 0.7), X = rbind(c(1, 2), c(1, 2)),
    keys = list(idx = c(5L, 5L), w = c(0.5, 0.9)))
  expect_equal(ag2$n, c(1L, 1L))
})


test_that("aggregate.pos reduces and preserves the single-block beta cover() fit", {
  skip_on_cran()
  skip_if_fast()
  s <- .aop_sim(101L, trend = FALSE)
  ff <- .aop_fit(s, trend = FALSE, agg.pos = FALSE)
  fp <- .aop_fit(s, trend = FALSE, agg.pos = TRUE)

  # Reduction fires: the grouped positive arm has strictly fewer rows.
  enc <- tulpaObs:::encode_cover_hurdle(
    ~ x + bym2(graph = s$adj, group_var = "region"), s$data, s$y, positive = "beta")
  expect_gt(ff$n_positive, 0L)

  .aop_expect_identical_fits(fp, ff)

  # And aggregating BOTH arms is still byte-identical to the full fit.
  fb <- .aop_fit(s, trend = FALSE, agg.occ = TRUE, agg.pos = TRUE)
  .aop_expect_identical_fits(fb, ff)
})


test_that("aggregate.pos defaults ON for the beta arm", {
  skip_on_cran()
  skip_if_fast()
  s <- .aop_sim(101L, trend = FALSE)
  # Default control (aggregate.pos unset) must aggregate the beta positive arm:
  # byte-identical both to the explicit full per-plot fit and to the explicit
  # aggregate.pos = TRUE fit -- the property that licenses the flipped default.
  ctrl <- list(verbose = FALSE, sigma.grid = c(0.5, 0.8, 1.2), rho.grid = 0.5,
               phi.grid = c(8, 18, 40), adaptive.grid = FALSE, max.iter = 300L,
               sigma.pos.grid = c(0.4, 0.8))
  fd <- suppressWarnings(tobs(
    formula = ~ x + bym2(graph = s$adj, group_var = "region"),
    data = s$data, family = cover("beta"), y = s$y,
    method = "nested_laplace", control = ctrl))   # default -> pos arm aggregated
  ff <- .aop_fit(s, trend = FALSE, agg.pos = FALSE)
  fp <- .aop_fit(s, trend = FALSE, agg.pos = TRUE)
  .aop_expect_identical_fits(fd, ff)
  .aop_expect_identical_fits(fd, fp)
})


test_that("aggregate.pos reduces and preserves the coupled-trend beta cover() fit", {
  skip_on_cran()
  skip_if_fast()
  for (seed in c(101L, 202L)) {
    s <- .aop_sim(seed, trend = TRUE)
    ff <- .aop_fit(s, trend = TRUE, agg.pos = FALSE)
    fp <- .aop_fit(s, trend = TRUE, agg.pos = TRUE)
    .aop_expect_identical_fits(fp, ff)
    fb <- .aop_fit(s, trend = TRUE, agg.occ = TRUE, agg.pos = TRUE)
    .aop_expect_identical_fits(fb, ff)
  }
})


# Multi-block path (spatial + AR1 temporal + IID RE). The positive linear
# predictor carries the spatial cell, the AR1 year and the IID observer, so the
# grouped collapse fires only when plots share ALL of them; replicate plots per
# (cell, year, obs) combo with a cell-level positive design make that happen.
.aop_sim_multi <- function(seed, n_s = 6L, n_year = 2L, n_obs = 2L, reps = 4L) {
  set.seed(seed); adj <- .aop_chain_adj(n_s); f1 <- .aop_icar_f(adj)
  combos <- expand.grid(cell = seq_len(n_s), year = seq_len(n_year), obs = seq_len(n_obs))
  combos <- combos[rep(seq_len(nrow(combos)), each = reps), ]
  cell <- combos$cell; N <- length(cell)
  x   <- as.numeric(scale(stats::rnorm(n_s)))[cell]              # cell-level
  ar  <- (as.numeric(scale(stats::rnorm(n_year))) * 0.3)[combos$year]
  iid <- stats::rnorm(n_obs, 0, 0.25)[combos$obs]
  eta_o <- -0.2 + 0.6 * x + 0.8 * f1[cell] + ar + iid
  occur <- stats::rbinom(N, 1L, stats::plogis(eta_o))
  eta_p <- 0.2 - 0.4 * x + 0.8 * f1[cell] + ar + iid
  mu <- stats::plogis(eta_p); y <- numeric(N); pos <- occur == 1L
  y[pos] <- pmin(pmax(stats::rbeta(sum(pos), mu[pos] * 18, (1 - mu[pos]) * 18), 1e-6), 1 - 1e-6)
  list(data = data.frame(x = x, region = factor(cell),
                         year = combos$year, obs = combos$obs),
       y = y, adj = adj)
}

.aop_fit_multi <- function(s, agg.pos) suppressWarnings(tobs(
  formula = ~ x + bym2(graph = s$adj, group_var = "region") +
              temporal(year, type = "ar1") + re(obs, type = "iid"),
  data = s$data, family = cover("beta"), y = s$y, method = "nested_laplace",
  control = list(verbose = FALSE, aggregate.occ = FALSE, aggregate.pos = agg.pos,
                 sigma.grid = c(0.5, 1.0), rho.grid = 0.5, sigma.pos.grid = c(0.5, 1.0),
                 tau.temporal.grid = 4, rho.temporal.grid = 0.5, sigma.re.grid = 0.25,
                 phi.grid = c(8, 18, 40), adaptive.grid = FALSE, max.iter = 300L)))

test_that("aggregate.pos reduces and preserves the multi-block beta cover() fit", {
  skip_on_cran()
  skip_if_fast()
  s <- .aop_sim_multi(101L)

  # The collapse genuinely fires: the (cell, year, obs) key on a cell-level
  # positive design has fewer levels than the positive rows.
  n_key <- nrow(unique(s$data[s$y > 0, c("x", "region", "year", "obs")]))
  expect_lt(n_key, sum(s$y > 0))

  ff <- .aop_fit_multi(s, agg.pos = FALSE)
  fp <- .aop_fit_multi(s, agg.pos = TRUE)
  .aop_expect_identical_fits(fp, ff)
})


test_that("aggregate.pos errors for the lognormal positive arm", {
  s <- .aop_sim(101L, trend = FALSE)
  expect_error(
    suppressWarnings(tobs(
      formula = ~ x + bym2(graph = s$adj, group_var = "region"),
      data = s$data, family = cover("lognormal"), y = s$y,
      method = "nested_laplace",
      control = list(verbose = FALSE, aggregate.pos = TRUE,
                     sigma.grid = c(0.5, 1.0), rho.grid = 0.5,
                     sigma.pos.grid = c(0.4, 0.8), adaptive.grid = FALSE))),
    "aggregate.pos"
  )
})
