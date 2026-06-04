# =============================================================================
# test-cover-hurdle-aggregate-occ.R - control$aggregate.occ on every cover()
# dispatch path (single-block, coupled-trend, multi-block).
#
# `aggregate.occ = TRUE` collapses the occurrence (binomial) arm to its exact
# sufficient statistic: observations agreeing on the occurrence design row AND
# every per-observation component of the linear predictor (spatial cell, any
# temporal / RE block index, any trend weight) are exchangeable Bernoulli trials,
# so they reduce to one Binomial row (n = count, y = successes). Because the
# engine's binomial kernel carries no combinatorial constant, this leaves the
# log-likelihood, gradient and Hessian pointwise unchanged -- it is a row-count
# reduction, not an approximation. The decisive test is therefore EQUIVALENCE on
# each path: the aggregated fit must reproduce the unaggregated fit numerically.
# Parameter recovery is inherited from the unaggregated paths (the other cover
# nested-joint / trend / multi-block tests), which the aggregated path equals
# bit-for-bit here.
# =============================================================================

.aoc_chain_adj <- function(N) {
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  adj
}

.aoc_grid2x2_adj <- function() {
  n_s <- 4L; g <- 2L; adj <- matrix(0L, n_s, n_s)
  for (i in 1:g) for (j in 1:g) {
    s <- (i - 1L) * g + j
    if (i > 1L) adj[s, (i - 2L) * g + j]       <- 1L
    if (i < g)  adj[s, i * g + j]              <- 1L
    if (j > 1L) adj[s, (i - 1L) * g + (j - 1L)] <- 1L
    if (j < g)  adj[s, (i - 1L) * g + (j + 1L)] <- 1L
  }
  adj
}

.aoc_icar_f <- function(adj) {
  Q   <- tulpaObs:::.occu_cover_icar_Q(adj)
  eig <- eigen(Q, symmetric = TRUE); keep <- eig$values > 1e-8
  z <- stats::rnorm(sum(keep))
  f <- as.numeric(eig$vectors[, keep, drop = FALSE] %*% (z / sqrt(eig$values[keep])))
  (f - mean(f)) / stats::sd(f)
}

# Equivalence of every comparable summary between two cover_fit objects.
.aoc_expect_equal_fits <- function(f1, f0, extra = character(), tol = 1e-8) {
  expect_true(f0$converged && f1$converged)
  expect_equal(f1$beta_occ, f0$beta_occ, tolerance = tol)
  expect_equal(f1$beta_pos, f0$beta_pos, tolerance = tol)
  expect_equal(f1$se_occ,   f0$se_occ,   tolerance = tol)
  expect_equal(f1$se_pos,   f0$se_pos,   tolerance = tol)
  # The global marginal likelihood is the most sensitive single summary: the
  # binomial kernel carries no lchoose constant, so the aggregated and per-plot
  # log-marginals are equal, not merely shifted by a constant.
  expect_equal(f1$log_marginal, f0$log_marginal, tolerance = tol)
  for (e in extra) expect_equal(f1[[e]], f0[[e]], tolerance = tol)
}


# --- (1) the aggregator itself ----------------------------------------------

test_that(".cover_aggregate_occ collapses exchangeable rows to exact sufficient stats", {
  y   <- c(1, 0, 1, 1, 0, 1)
  X   <- rbind(c(1, 2), c(1, 2), c(1, 2), c(3, 4), c(3, 4), c(1, 2))
  spi <- c(5L, 5L, 5L, 7L, 7L, 6L)
  w   <- c(0.5, 0.5, 0.5, 0.9, 0.9, 0.5)

  ag <- tulpaObs:::.cover_aggregate_occ(y, X, list(spi = spi, w = w))

  # Three groups, ordered by group id: (X=1,2;spi=5;w=.5) rows 1-3,
  # (X=1,2;spi=6;w=.5) row 6, (X=3,4;spi=7;w=.9) rows 4-5.
  expect_equal(ag$y, c(2, 1, 1))
  expect_equal(ag$n, c(3L, 1L, 2L))
  expect_equal(ag$X, rbind(c(1, 2), c(1, 2), c(3, 4)))
  expect_equal(ag$keys$spi, c(5L, 6L, 7L))
  expect_equal(ag$keys$w,   c(0.5, 0.5, 0.9))

  # Sufficiency invariants: successes and trials are conserved.
  expect_equal(sum(ag$y), sum(y))
  expect_equal(sum(ag$n), length(y))

  # Every key is part of the grouping key: same (X, cell) but a different
  # weight must NOT merge (different linear predictor).
  ag2 <- tulpaObs:::.cover_aggregate_occ(
    y = c(1, 0), X = rbind(c(1, 2), c(1, 2)),
    keys = list(spi = c(5L, 5L), w = c(0.5, 0.9)))
  expect_equal(ag2$n, c(1L, 1L))
  expect_equal(nrow(ag2$X), 2L)

  # A second index axis (e.g. a temporal code) refines the key: same (X, cell)
  # but different time does not merge.
  ag3 <- tulpaObs:::.cover_aggregate_occ(
    y = c(1, 1), X = rbind(c(1, 2), c(1, 2)),
    keys = list(spi = c(5L, 5L), t = c(1L, 2L)))
  expect_equal(ag3$n, c(1L, 1L))
  expect_equal(ag3$keys$t, c(1L, 2L))

  # Row order does not change the aggregated result (groups keyed, not positional).
  ord <- c(4L, 1L, 6L, 2L, 5L, 3L)
  agp <- tulpaObs:::.cover_aggregate_occ(
    y[ord], X[ord, , drop = FALSE], list(spi = spi[ord], w = w[ord]))
  key <- function(a) order(a$keys$spi)
  expect_equal(ag$y[key(ag)],        agp$y[key(agp)])
  expect_equal(ag$n[key(ag)],        agp$n[key(agp)])
  expect_equal(ag$keys$spi[key(ag)], agp$keys$spi[key(agp)])
})


# --- (2) single-block path --------------------------------------------------

.aoc_sim_single <- function(n_s = 16L, n_per = 6L, seed = 11L,
                            sigma = 0.7, alpha = 1.0, sd_pos = 0.4) {
  set.seed(seed); adj <- .aoc_chain_adj(n_s); f1 <- .aoc_icar_f(adj)
  cell <- rep(seq_len(n_s), each = n_per); N <- length(cell)
  x <- as.numeric(scale(stats::rnorm(n_s)))[cell]      # cell-level occ covariate
  eta_o <- -0.3 + 0.7 * x + sigma * f1[cell]
  occur <- stats::rbinom(N, 1L, stats::plogis(eta_o))
  eta_p <- 0.4 - 0.5 * x + alpha * sigma * f1[cell]
  y <- ifelse(occur == 1L, pmin(exp(stats::rnorm(N, eta_p, sd_pos)), 1 - 1e-6), 0)
  list(data = data.frame(x = x, region = factor(cell)), y = y, adj = adj, n_s = n_s, N = N)
}

.aoc_fit_single <- function(s, agg) suppressWarnings(tobs(
  formula = ~ x + bym2(graph = s$adj, group_var = "region"),
  data = s$data, family = cover("lognormal"), y = s$y, method = "nested_laplace",
  control = list(verbose = FALSE, aggregate.occ = agg, sigma.grid = c(0.5, 1.0),
                 rho.grid = 0.5, sigma.pos.grid = c(0.4, 0.8), phi.grid = c(0.3, 0.5),
                 adaptive.grid = FALSE)))

test_that("aggregate.occ reduces and preserves the single-block cover() fit", {
  skip_if_fast()
  s <- .aoc_sim_single()

  # Reduction actually fires: cell-level design collapses occ rows to one
  # Binomial row per cell.
  enc <- tulpaObs:::encode_cover_hurdle(
    ~ x + bym2(graph = s$adj, group_var = "region"), s$data, s$y, positive = "lognormal")
  spi <- tulpa::prior_from_spec(
    enc$spatial_spec, s$data[enc$obs_keep, , drop = FALSE])$spatial_idx
  og  <- tulpaObs:::.cover_aggregate_occ(enc$occ_data$y, enc$occ_data$X, list(spi = spi))
  expect_equal(length(og$y), s$n_s)
  expect_lt(length(og$y), nrow(enc$occ_data$X))
  expect_equal(sum(og$y), sum(enc$occ_data$y))
  expect_equal(sum(og$n), nrow(enc$occ_data$X))

  .aoc_expect_equal_fits(.aoc_fit_single(s, TRUE), .aoc_fit_single(s, FALSE),
                         extra = "sigma_pos")
})


# --- (3) coupled-trend path -------------------------------------------------

.aoc_sim_trend <- function(n_s = 16L, n_per = 6L, seed = 7L,
                           sigma = 0.8, alpha = 1.0,
                           sigma_tr = 0.6, alpha_tr = 0.9, sd_pos = 0.4) {
  set.seed(seed); adj <- .aoc_chain_adj(n_s)
  f1 <- .aoc_icar_f(adj); f2 <- .aoc_icar_f(adj)
  cell  <- rep(seq_len(n_s), each = n_per); N <- length(cell)
  x     <- as.numeric(scale(stats::rnorm(n_s)))[cell]
  time  <- as.numeric(scale(stats::rnorm(n_s)))[cell]
  eta_o <- -0.3 + 0.7 * x + sigma * f1[cell] + sigma_tr * f2[cell] * time
  occur <- stats::rbinom(N, 1L, stats::plogis(eta_o))
  eta_p <- 0.4 - 0.5 * x + alpha * sigma * f1[cell] + alpha_tr * sigma_tr * f2[cell] * time
  y <- ifelse(occur == 1L, pmin(exp(stats::rnorm(N, eta_p, sd_pos)), 1 - 1e-6), 0)
  list(data = data.frame(x = x, time = time, region = factor(cell)),
       y = y, adj = adj, n_s = n_s, N = N)
}

.aoc_fit_trend <- function(s, agg) suppressWarnings(tobs(
  formula = ~ x + bym2(graph = s$adj, group_var = "region"),
  data = s$data, family = cover("lognormal"), y = s$y, method = "nested_laplace",
  control = list(verbose = FALSE, trend = list(weight = "time"), aggregate.occ = agg,
                 sigma.grid = c(0.5, 1.0), rho.grid = 0.5,
                 alpha.grid = c(0, 1.0), alpha.grid.trend = c(0, 1.0),
                 phi.grid = c(0.3, 0.5), adaptive.grid = FALSE)))

test_that("aggregate.occ reduces and preserves the coupled-trend cover() fit", {
  skip_if_fast()
  for (seed in c(7L, 19L)) {
    s <- .aoc_sim_trend(seed = seed)
    if (seed == 7L) {
      enc <- tulpaObs:::encode_cover_hurdle(
        ~ x + bym2(graph = s$adj, group_var = "region"), s$data, s$y, positive = "lognormal")
      data_obs <- s$data[enc$obs_keep, , drop = FALSE]
      spi <- tulpa::prior_from_spec(enc$spatial_spec, data_obs)$spatial_idx
      og  <- tulpaObs:::.cover_aggregate_occ(
        enc$occ_data$y, enc$occ_data$X, list(spi = spi, w = data_obs$time))
      expect_equal(length(og$y), s$n_s)
      expect_lt(length(og$y), nrow(enc$occ_data$X))
    }
    .aoc_expect_equal_fits(.aoc_fit_trend(s, TRUE), .aoc_fit_trend(s, FALSE),
                           extra = c("sigma_pos", "sigma_trend", "alpha_trend"))
  }
})


# --- (4) multi-block path (spatial + AR1 temporal + IID RE) ------------------

.aoc_sim_multi <- function(seed = 7001L, n_s = 4L, n_years = 3L, n_obs = 3L, N = 140L,
                           sigma = 0.6, rho = 0.7, alpha = 1.1, sigma_year = 0.3,
                           sigma_obs = 0.25, rho_ar = 0.6, phi_b = 30) {
  set.seed(seed); adj <- .aoc_grid2x2_adj()
  s_idx <- sample.int(n_s, N, TRUE)
  t_idx <- sample.int(n_years, N, TRUE)
  o_idx <- sample.int(n_obs, N, TRUE)
  phi_f <- stats::rnorm(n_s); phi_f <- phi_f - mean(phi_f)
  th_f  <- stats::rnorm(n_s); th_f  <- th_f  - mean(th_f)
  w_u <- sqrt(rho) * phi_f + sqrt(1 - rho) * th_f
  ar <- numeric(n_years); ar[1] <- stats::rnorm(1, 0, sigma_year)
  for (t in 2:n_years) ar[t] <- rho_ar * ar[t - 1] +
      stats::rnorm(1, 0, sigma_year * sqrt(1 - rho_ar^2))
  iota <- stats::rnorm(n_obs, 0, sigma_obs)
  x <- as.numeric(scale(stats::rnorm(n_s)))[s_idx]     # cell-level occ covariate
  eta_o <- 0.2 + 0.7 * x + sigma * w_u[s_idx] + ar[t_idx] + iota[o_idx]
  occur <- stats::rbinom(N, 1L, stats::plogis(eta_o))
  eta_p <- -0.5 - 0.3 * x + alpha * sigma * w_u[s_idx] + ar[t_idx] + iota[o_idx]
  mu <- stats::plogis(eta_p); cp <- stats::rbeta(N, mu * phi_b, (1 - mu) * phi_b)
  cp <- pmin(pmax(cp, 1e-6), 1 - 1e-6)
  y <- ifelse(occur == 1L, cp, 0)
  list(data = data.frame(x = x, region = factor(s_idx, levels = seq_len(n_s)),
                         year = t_idx, obs = o_idx), y = y, adj = adj, N = N)
}

.aoc_fit_multi <- function(s, agg) suppressWarnings(tobs(
  formula = ~ x + bym2(graph = s$adj, group_var = "region") +
              temporal(year, type = "ar1") + re(obs, type = "iid"),
  data = s$data, family = cover("beta"), y = s$y, method = "nested_laplace",
  control = list(verbose = FALSE, aggregate.occ = agg, sigma.grid = c(0.4, 0.8),
                 rho.grid = 0.6, sigma.pos.grid = c(0.5, 1.0), tau.temporal.grid = 4,
                 rho.temporal.grid = 0.5, sigma.re.grid = 0.25, phi.grid = c(10, 40),
                 adaptive.grid = FALSE)))

test_that("aggregate.occ preserves the multi-block cover() fit on the full key", {
  skip_if_fast()
  s <- .aoc_sim_multi()

  # The occurrence linear predictor depends on cell, year AND observer, so the
  # aggregation key is (X-row, region, year, obs). On this design that key has
  # far fewer levels than N -- the occ arm genuinely collapses -- while no two
  # rows differing in time / observer may merge. Equivalence below is the proof
  # that the finer key is respected (a wrong merge would change the likelihood).
  n_key <- nrow(unique(s$data[, c("region", "year", "obs")]))
  expect_lt(n_key, s$N)

  .aoc_expect_equal_fits(.aoc_fit_multi(s, TRUE), .aoc_fit_multi(s, FALSE),
                         extra = "phi_pos")
})
