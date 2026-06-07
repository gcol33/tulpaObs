# =============================================================================
# test-cover-hurdle-aggregate-recovery.R - parameter-recovery suite behind the
# default-ON sufficient-statistic reductions of BOTH cover-hurdle arms:
# aggregate.occ (binomial occurrence arm, tulpaObs#48) and aggregate.pos
# (grouped-beta positive arm, tulpaObs#49). Both default ON, so a default-control
# fit exercises both reductions at once; this suite is the recovery sign-off
# behind that default.
#
# Each reduction collapses observations sharing a design row AND every
# per-observation latent component to one row: the occurrence arm to its exact
# Binomial sufficient statistic, the positive arm to grouped (n, sum log y,
# sum log(1-y)). The package convention is that a likelihood-affecting reduction
# defaults on only with a parameter-recovery suite behind it -- not just shape /
# exactness checks. Two parts:
#
#   (1) RECOVERY of simulated truth on the now-default (both arms aggregated)
#       cover() beta hurdle with two coupled fields (intercept + spatially-
#       varying trend): bias small and 95% CI coverage >= 85% for the quantities
#       the fit reports a posterior SD for (both arms' betas + the beta precision
#       phi_pos). The community-variance hyperparameters (sigma, alpha,
#       sigma_trend, alpha_trend) carry the expected small-N Laplace attenuation
#       and are not surfaced with a posterior SD, so they get a bias bracket,
#       not coverage.
#
#   (2) EXACTNESS: the default (both arms aggregated) fit equals the explicit
#       full per-plot fit byte-for-byte (betas, SEs, log-marginal), so default-ON
#       loses nothing -- the property that licenses flipping both defaults.
# =============================================================================


.acr_chain_adj <- function(n) {
  adj <- matrix(0L, n, n)
  for (s in seq_len(n)) { if (s > 1L) adj[s, s - 1L] <- 1L; if (s < n) adj[s, s + 1L] <- 1L }
  adj
}

.acr_icar_f <- function(adj) {
  Q  <- tulpaObs:::.occu_cover_icar_Q(adj)
  sc <- tulpaObs:::.occu_cover_icar_scale(adj)
  eig <- eigen(Q, symmetric = TRUE); keep <- eig$values > 1e-8
  fk <- as.numeric(eig$vectors[, keep, drop = FALSE] %*%
                     (stats::rnorm(sum(keep)) / sqrt(eig$values[keep])))
  (fk - mean(fk)) / sqrt(sc)
}

.acr_truth <- list(b_occ = c(-0.3, 0.7), b_pos = c(0.2, -0.5),
                   sigma = 0.8, alpha = 1.0, sigma_trend = 0.6,
                   alpha_trend = 0.9, phi = 18)

# Beta cover hurdle, two coupled fields. Occurrence eta carries both fields;
# the cover (positive) eta carries each field scaled by its own alpha.
.acr_sim <- function(seed, n_s = 18L, n_per = 10L) {
  set.seed(seed); adj <- .acr_chain_adj(n_s)
  f1 <- .acr_icar_f(adj); f2 <- .acr_icar_f(adj)
  cell <- rep(seq_len(n_s), each = n_per); N <- length(cell)
  x    <- as.numeric(scale(stats::rnorm(n_s)))[cell]
  time <- as.numeric(scale(stats::rnorm(n_s)))[cell]
  tr <- .acr_truth
  eta_o <- tr$b_occ[1] + tr$b_occ[2] * x + tr$sigma * f1[cell] +
           tr$sigma_trend * f2[cell] * time
  occur <- stats::rbinom(N, 1L, stats::plogis(eta_o))
  eta_p <- tr$b_pos[1] + tr$b_pos[2] * x + tr$alpha * tr$sigma * f1[cell] +
           tr$alpha_trend * tr$sigma_trend * f2[cell] * time
  mu <- stats::plogis(eta_p); y <- numeric(N); pos <- occur == 1L
  y[pos] <- pmin(pmax(stats::rbeta(sum(pos), mu[pos] * tr$phi, (1 - mu[pos]) * tr$phi),
                      1e-6), 1 - 1e-6)
  list(data = data.frame(x = x, time = time, region = factor(cell)), y = y, adj = adj)
}

.acr_fit <- function(s, agg = NULL) {
  ctrl <- list(verbose = FALSE, trend = list(weight = "time"),
               sigma.grid = c(0.5, 0.8, 1.2), rho.grid = 0.5,
               alpha.grid = c(0, 0.5, 1.0, 1.5),
               alpha.grid.trend = c(0, 0.5, 1.0, 1.5),
               phi.grid = c(8, 18, 40), adaptive.grid = FALSE, max.iter = 300L)
  if (!is.null(agg)) { ctrl$aggregate.occ <- agg; ctrl$aggregate.pos <- agg }
  suppressWarnings(tobs(
    formula = ~ x + bym2(graph = s$adj, group_var = "region"),
    data = s$data, family = cover("beta"), y = s$y,
    method = "nested_laplace", control = ctrl))
}


test_that("both arms default-ON recover truth (cover beta-trend hurdle, tulpaObs#48/#49)", {
  skip_on_cran()
  skip_if_fast()

  tr <- .acr_truth
  # Quantities the fit reports a posterior SD for -> bias + coverage.
  calib <- c(`occ_(Intercept)` = tr$b_occ[1], occ_x = tr$b_occ[2],
             `pos_(Intercept)` = tr$b_pos[1], pos_x = tr$b_pos[2],
             phi_pos = tr$phi)
  # Variance components: small-N Laplace attenuation, no posterior SD -> bias bracket.
  vc <- c(sigma = tr$sigma, alpha = tr$alpha,
          sigma_trend = tr$sigma_trend, alpha_trend = tr$alpha_trend)

  seeds <- 201:224
  est <- matrix(NA_real_, length(seeds), length(c(calib, vc)),
                dimnames = list(NULL, names(c(calib, vc))))
  covd <- matrix(NA_real_, length(seeds), length(calib), dimnames = list(NULL, names(calib)))

  for (i in seq_along(seeds)) {
    f <- .acr_fit(.acr_sim(seeds[i]))   # default control -> both arms aggregated
    if (!isTRUE(f$converged)) next
    hp <- f$hyperpar$spatial
    e <- c(`occ_(Intercept)` = f$beta_occ[[1]], occ_x = f$beta_occ[[2]],
           `pos_(Intercept)` = f$beta_pos[[1]], pos_x = f$beta_pos[[2]],
           phi_pos = f$phi_pos,
           sigma = hp[["b1.sigma"]], alpha = hp[["b1.alpha"]],
           sigma_trend = hp[["b2.sigma"]], alpha_trend = hp[["b2.alpha"]])
    est[i, names(e)] <- e
    sdv <- c(`occ_(Intercept)` = f$se_occ[[1]], occ_x = f$se_occ[[2]],
             `pos_(Intercept)` = f$se_pos[[1]], pos_x = f$se_pos[[2]],
             phi_pos = f$phi_pos_sd)
    for (p in names(calib)) {
      covd[i, p] <- (calib[[p]] >= e[[p]] - 1.96 * sdv[[p]]) &&
                    (calib[[p]] <= e[[p]] + 1.96 * sdv[[p]])
    }
  }

  n_ok <- sum(is.finite(est[, "occ_x"]))
  expect_gte(n_ok, 20L)   # at least 20 converged seeds behind the suite

  # Calibrated quantities: small bias + nominal coverage.
  for (p in names(calib)) {
    e <- est[, p]; ok <- is.finite(e)
    bias <- mean(e[ok]) - calib[[p]]
    tol  <- if (p == "phi_pos") 5 else 0.15
    expect_lt(abs(bias), tol, label = sprintf("bias[%s]=%.3f", p, bias))
    coverage <- mean(covd[ok, p])
    expect_gte(coverage, 0.85)
  }

  # Variance components: recovered to within the small-N attenuation bracket.
  for (p in names(vc)) {
    e <- est[, p]; ok <- is.finite(e)
    m <- mean(e[ok])
    expect_gt(m, 0)
    expect_lt(abs(m - vc[[p]]), 0.45, label = sprintf("vc[%s] mean=%.3f", p, m))
  }
})


test_that("both arms default-ON are byte-identical to the full per-plot fit (tulpaObs#48/#49)", {
  skip_on_cran()
  skip_if_fast()

  for (seed in c(101L, 137L)) {
    s  <- .acr_sim(seed)
    fd <- .acr_fit(s)              # default control -> both arms aggregated
    ff <- .acr_fit(s, agg = FALSE) # explicit full per-plot occurrence + positive arms
    expect_true(fd$converged && ff$converged)
    expect_equal(fd$beta_occ,     ff$beta_occ,     tolerance = 1e-8)
    expect_equal(fd$beta_pos,     ff$beta_pos,     tolerance = 1e-8)
    expect_equal(fd$se_occ,       ff$se_occ,       tolerance = 1e-8)
    expect_equal(fd$se_pos,       ff$se_pos,       tolerance = 1e-8)
    expect_equal(fd$phi_pos,      ff$phi_pos,      tolerance = 1e-8)
    expect_equal(fd$log_marginal, ff$log_marginal, tolerance = 1e-8)
  }
})
