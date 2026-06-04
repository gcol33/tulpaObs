# =============================================================================
# test-cover-hurdle-aggregate-occ.R - control$aggregate.occ on the coupled-trend
# cover() path.
#
# `aggregate.occ = TRUE` collapses the occurrence (binomial) arm to its exact
# sufficient statistic: plots sharing an occurrence design row, spatial cell, and
# per-observation trend weight are exchangeable Bernoulli trials, so they reduce
# to one Binomial row (n = count, y = successes). Because the engine's binomial
# kernel carries no combinatorial constant, this leaves the log-likelihood,
# gradient and Hessian pointwise unchanged -- it is a row-count reduction, not an
# approximation. The decisive test is therefore EQUIVALENCE: the aggregated fit
# must reproduce the unaggregated fit numerically. Parameter recovery is then
# inherited from the unaggregated trend path (test-occu-cover-trend.R), which the
# aggregated path equals bit-for-bit here.
# =============================================================================

.aoc_chain_adj <- function(N) {
  adj <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) adj[s, s - 1L] <- 1L
    if (s < N)  adj[s, s + 1L] <- 1L
  }
  adj
}

# Cell-level cover-hurdle with a coupled spatially-varying trend. The occurrence
# covariate `x` and the trend weight `time` are constant within a cell, so the
# n_per plots of a cell are exchangeable on the occurrence arm and collapse to
# one Binomial row; the positive (lognormal) arm stays per-plot.
.aoc_sim <- function(n_s = 16L, n_per = 6L, seed = 7L,
                     sigma = 0.8, alpha = 1.0,
                     sigma_tr = 0.6, alpha_tr = 0.9, sd_pos = 0.4) {
  set.seed(seed)
  adj <- .aoc_chain_adj(n_s)
  Q   <- tulpaObs:::.occu_cover_icar_Q(adj)
  eig <- eigen(Q, symmetric = TRUE); keep <- eig$values > 1e-8
  draw_f <- function() {
    z <- stats::rnorm(sum(keep))
    f <- as.numeric(eig$vectors[, keep, drop = FALSE] %*% (z / sqrt(eig$values[keep])))
    (f - mean(f)) / stats::sd(f)
  }
  f1 <- draw_f(); f2 <- draw_f()
  cell  <- rep(seq_len(n_s), each = n_per); N <- length(cell)
  xcell <- as.numeric(scale(stats::rnorm(n_s)))
  tcell <- as.numeric(scale(stats::rnorm(n_s)))
  x <- xcell[cell]; time <- tcell[cell]
  eta_occ <- -0.3 + 0.7 * x + sigma * f1[cell] + sigma_tr * f2[cell] * time
  occur   <- stats::rbinom(N, 1L, stats::plogis(eta_occ))
  eta_pos <- 0.4 - 0.5 * x + alpha * sigma * f1[cell] + alpha_tr * sigma_tr * f2[cell] * time
  y <- ifelse(occur == 1L, pmin(exp(stats::rnorm(N, eta_pos, sd_pos)), 1 - 1e-6), 0)
  list(data = data.frame(x = x, time = time, region = factor(cell)),
       y = y, adj = adj, n_s = n_s, N = N)
}

.aoc_fit <- function(s, agg, max.iter = 200L) {
  suppressWarnings(tobs(
    formula = ~ x + bym2(graph = s$adj, group_var = "region"),
    data = s$data, family = cover("lognormal"), y = s$y,
    method = "nested_laplace",
    control = list(verbose = FALSE, max.iter = max.iter,
                   trend = list(weight = "time"), aggregate.occ = agg,
                   sigma.grid = c(0.5, 1.0), rho.grid = 0.5,
                   alpha.grid = c(0, 1.0), alpha.grid.trend = c(0, 1.0),
                   phi.grid = c(0.3, 0.5), adaptive.grid = FALSE)))
}


test_that(".cover_aggregate_occ collapses exchangeable rows to exact sufficient stats", {
  y   <- c(1, 0, 1, 1, 0, 1)
  X   <- rbind(c(1, 2), c(1, 2), c(1, 2), c(3, 4), c(3, 4), c(1, 2))
  spi <- c(5L, 5L, 5L, 7L, 7L, 6L)
  w   <- c(0.5, 0.5, 0.5, 0.9, 0.9, 0.5)

  ag <- tulpaObs:::.cover_aggregate_occ(y, X, spi, w)

  # Three groups, ordered by group id: (X=1,2;spi=5;w=.5) rows 1-3,
  # (X=1,2;spi=6;w=.5) row 6, (X=3,4;spi=7;w=.9) rows 4-5.
  expect_equal(ag$y,   c(2, 1, 1))
  expect_equal(ag$n,   c(3L, 1L, 2L))
  expect_equal(ag$spi, c(5L, 6L, 7L))
  expect_equal(ag$w,   c(0.5, 0.5, 0.9))
  expect_equal(ag$X,   rbind(c(1, 2), c(1, 2), c(3, 4)))

  # Sufficiency invariants: successes and trials are conserved.
  expect_equal(sum(ag$y), sum(y))
  expect_equal(sum(ag$n), length(y))

  # The weight is part of the grouping key: same (X, cell) but a different
  # weight must NOT merge (different linear predictor).
  ag2 <- tulpaObs:::.cover_aggregate_occ(
    y = c(1, 0), X = rbind(c(1, 2), c(1, 2)), spi = c(5L, 5L), w = c(0.5, 0.9))
  expect_equal(ag2$n, c(1L, 1L))
  expect_equal(nrow(ag2$X), 2L)

  # Row order does not change the aggregated result (groups keyed, not positional).
  ord <- c(4L, 1L, 6L, 2L, 5L, 3L)
  agp <- tulpaObs:::.cover_aggregate_occ(y[ord], X[ord, , drop = FALSE], spi[ord], w[ord])
  key  <- function(a) order(a$spi)
  expect_equal(ag$y[key(ag)],   agp$y[key(agp)])
  expect_equal(ag$n[key(ag)],   agp$n[key(agp)])
  expect_equal(ag$spi[key(ag)], agp$spi[key(agp)])

  # The trend-free (w = NULL) branch groups on (X, cell) alone.
  ag0 <- tulpaObs:::.cover_aggregate_occ(y, X, spi, w = NULL)
  expect_null(ag0$w)
  expect_equal(sum(ag0$n), length(y))
})


test_that("aggregate.occ reduces the occurrence row count", {
  s <- .aoc_sim()
  enc <- tulpaObs:::encode_cover_hurdle(
    ~ x + bym2(graph = s$adj, group_var = "region"),
    s$data, s$y, positive = "lognormal")
  data_obs <- s$data[enc$obs_keep, , drop = FALSE]
  spi <- tulpa::prior_from_spec(enc$spatial_spec, data_obs)$spatial_idx
  og  <- tulpaObs:::.cover_aggregate_occ(
    enc$occ_data$y, enc$occ_data$X, spi, data_obs$time)

  # Cell-level covariates: occurrence rows collapse to one per cell.
  expect_lt(length(og$y), nrow(enc$occ_data$X))
  expect_equal(length(og$y), s$n_s)
  expect_equal(sum(og$y), sum(enc$occ_data$y))
  expect_equal(sum(og$n), nrow(enc$occ_data$X))
})


test_that("aggregate.occ leaves the coupled-trend cover() fit numerically unchanged", {
  skip_if_fast()

  for (seed in c(7L, 19L)) {
    s  <- .aoc_sim(seed = seed)
    f0 <- .aoc_fit(s, agg = FALSE)
    f1 <- .aoc_fit(s, agg = TRUE)

    expect_true(f0$converged && f1$converged)

    tol <- 1e-8
    # Fixed effects and their SEs.
    expect_equal(f1$beta_occ, f0$beta_occ, tolerance = tol)
    expect_equal(f1$beta_pos, f0$beta_pos, tolerance = tol)
    expect_equal(f1$se_occ,   f0$se_occ,   tolerance = tol)
    expect_equal(f1$se_pos,   f0$se_pos,   tolerance = tol)

    # Hyperparameter posterior (both shared fields: sigma/alpha/rho per block)
    # and the positive-arm dispersion.
    expect_equal(f1$hyperpar$spatial, f0$hyperpar$spatial, tolerance = tol)
    expect_equal(f1$sigma_pos,   f0$sigma_pos,   tolerance = tol)
    expect_equal(f1$sigma_trend, f0$sigma_trend, tolerance = tol)
    expect_equal(f1$alpha_trend, f0$alpha_trend, tolerance = tol)

    # The global marginal likelihood is the most sensitive single summary: the
    # binomial kernel carries no lchoose constant, so the aggregated and
    # per-plot log-marginals are equal, not merely shifted by a constant.
    expect_equal(f1$log_marginal, f0$log_marginal, tolerance = tol)
  }
})
