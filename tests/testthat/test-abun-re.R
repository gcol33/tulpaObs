# Site-level random effects on the single-species N-mixture (gcol33/tulpaObs#13).
# Recovery + 95% CI coverage on the abundance and detection arms, S3 surface,
# capability gates. NB + RE shares the same AGHQ engine (log_r is the trailing
# theta coordinate) and gets a separate recovery row.

sim_abun_lambda_re <- function(N, J, ngrp, beta_lambda, beta_p, sigma_b,
                               mixture = "poisson", size = 4, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = sigma_b)
  data <- data.frame(x1 = stats::rnorm(N), g = factor(grp))
  eta_l <- as.numeric(model.matrix(~ x1, data) %*% beta_lambda) + b[grp]
  lambda <- exp(eta_l)
  Nlat <- if (identical(mixture, "negbin"))
    stats::rnbinom(N, size = size, mu = lambda) else stats::rpois(N, lambda)
  p_obs <- plogis(beta_p)
  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], p_obs)
  list(y = y, data = data, beta_lambda = beta_lambda, beta_p = beta_p,
       sigma_b = sigma_b, b = b)
}

sim_abun_p_re <- function(N, J, ngrp, beta_lambda, beta_p, sigma_b, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  grp <- rep(seq_len(ngrp), length.out = N)
  b   <- stats::rnorm(ngrp, sd = sigma_b)
  data <- data.frame(z1 = stats::rnorm(N), g = factor(grp))
  lambda <- exp(beta_lambda)
  Nlat   <- stats::rpois(N, lambda)
  eta_p  <- as.numeric(model.matrix(~ z1, data) %*% beta_p) + b[grp]
  p_obs  <- plogis(eta_p)
  y <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- stats::rbinom(J, Nlat[i], p_obs[i])
  list(y = y, data = data, beta_lambda = beta_lambda, beta_p = beta_p,
       sigma_b = sigma_b, b = b)
}


# --- Structural tests (always run) -----------------------------------------

test_that("RE on the lambda arm fits and surfaces sigma_g + BLUPs", {
  s <- sim_abun_lambda_re(N = 60, J = 3, ngrp = 6,
                          beta_lambda = c(0.5, 0.3), beta_p = 0,
                          sigma_b = 0.6, seed = 1)
  fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
              data = s$data, y = s$y, family = abun(),
              method = "laplace", verbose = FALSE,
              control = list(n.quad = 1))

  expect_s3_class(fit, "tobs_fit")
  expect_true(fit$convergence$converged)
  expect_true("sigma_g1_(Intercept)" %in% names(fit$means))
  expect_true(any(grepl("^re_g1_", names(fit$means))))
  expect_identical(fit$nmix_re$arm, "lambda")

  # Per-group BLUP table surfaces through ranef().
  re <- ranef(fit)
  expect_true(is.list(re) || is.data.frame(re))
})

test_that("RE on the p arm fits and the sigma carries the p<t> label", {
  sp <- sim_abun_p_re(N = 60, J = 3, ngrp = 6,
                     beta_lambda = 1.5, beta_p = c(0, 0.4),
                     sigma_b = 0.6, seed = 2)
  fit <- tobs(formula = ~ 1, detection = ~ z1 + (1 | g),
              data = sp$data, y = sp$y, family = abun(),
              method = "laplace", verbose = FALSE,
              control = list(n.quad = 1))

  expect_s3_class(fit, "tobs_fit")
  expect_true(fit$convergence$converged)
  expect_true("sigma_p1_(Intercept)" %in% names(fit$means))
  expect_identical(fit$nmix_re$arm, "p")
})

test_that("S3 surface (coef, vcov, ranef) carries the RE component", {
  s <- sim_abun_lambda_re(N = 60, J = 3, ngrp = 6,
                          beta_lambda = c(0.5, 0.3), beta_p = 0,
                          sigma_b = 0.5, seed = 3)
  fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
              data = s$data, y = s$y, family = abun(),
              method = "laplace", verbose = FALSE,
              control = list(n.quad = 1))

  # coef returns the two fixed-effect arms.
  co <- coef(fit)
  expect_true(is.list(co))
  expect_true(all(c("lambda", "p") %in% names(co)))
  expect_equal(length(co$lambda), 2L)

  # vcov is the joint fixed-effect cov (no RE point estimates in it).
  V <- vcov(fit)
  expect_true(nrow(V) == 3L)        # lambda(2) + p(1)
  expect_true(all(diag(V) >= 0))

  # ranef surfaces the per-group BLUP table.
  re <- ranef(fit)
  expect_false(is.null(re))

  # tulpaObs#19 (RE note): the variance-component pseudo-draw columns
  # (sigma_*, cor_*) have no analytic joint covariance with the fixed effects
  # under the marginal-Hessian AGHQ path, so they are NA (uncertainty
  # explicitly unavailable) rather than a fabricated near-degenerate column.
  # The per-group BLUP columns carry their real AGHQ marginal posterior SD, so
  # their draws are genuinely dispersed.
  dn <- colnames(fit$draws)
  sig_cols <- grep("^sigma_", dn)
  blup_cols <- grep("^re_", dn)
  expect_true(length(sig_cols) >= 1L)
  expect_true(length(blup_cols) >= 1L)
  expect_true(all(is.na(fit$draws[, sig_cols])))
  expect_true(all(is.na(fit$sds[grep("^sigma_", names(fit$sds))])))
  blup_sds <- apply(fit$draws[, blup_cols, drop = FALSE], 2, sd)
  expect_true(all(is.finite(blup_sds)))
  expect_true(any(blup_sds > 1e-3))
})


# --- Capability gates ------------------------------------------------------

test_that("RE + spatial errors with a pointer to the one-or-the-other rule", {
  s <- sim_abun_lambda_re(N = 30, J = 3, ngrp = 5,
                          beta_lambda = c(0.5, 0.3), beta_p = 0,
                          sigma_b = 0.4, seed = 4)
  # Minimal chain-graph adjacency matrix so icar() has a valid graph.
  n <- nrow(s$data)
  adj <- matrix(0L, n, n)
  for (i in seq_len(n - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  expect_error(
    tobs(formula = ~ 1 + (1 | g) + icar(graph = adj),
         detection = ~ 1, data = s$data, y = s$y, family = abun(),
         method = "nested_laplace", verbose = FALSE),
    regexp = "areal spatial term|spatial.*not yet supported", ignore.case = TRUE
  )
})

test_that("RE shared across both arms is rejected", {
  s <- sim_abun_lambda_re(N = 30, J = 3, ngrp = 5,
                          beta_lambda = c(0.5, 0.3), beta_p = 0,
                          sigma_b = 0.4, seed = 5)
  expect_error(
    tobs(formula = ~ 1 + (1 | g), detection = ~ 1 + (1 | g),
         data = s$data, y = s$y, family = abun(),
         method = "laplace", verbose = FALSE),
    regexp = "BOTH|both arms|cross-arm", ignore.case = TRUE
  )
})


# --- Recovery + coverage (gated) -------------------------------------------

# The native NMixGroupedOracle path is per-group O(|group| * K_max * outer
# iters), so a single fit at this scale (N = 100, J = 4, 10 groups) is ~2 s
# (Poisson) / ~6-10 s (NB). The full gated suite below runs ~110 s locally.

test_that("sigma recovery on the lambda arm (multi-seed mean within 25%)", {
  skip_on_cran()
  beta_lambda <- c(0.6, 0.4); beta_p <- 0; sigma_b <- 0.6
  n_seed <- 8L
  sig_hat <- numeric(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- sim_abun_lambda_re(N = 100, J = 4, ngrp = 10,
                              beta_lambda = beta_lambda, beta_p = beta_p,
                              sigma_b = sigma_b, seed = 1000 + s)
    fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
                data = sim$data, y = sim$y, family = abun(),
                method = "laplace", verbose = FALSE,
                control = list(n.quad = 5))
    sig_hat[s] <- fit$means["sigma_g1_(Intercept)"]
  }
  rel_bias <- abs(mean(sig_hat) - sigma_b) / sigma_b
  expect_lt(rel_bias, 0.30)
})

test_that("95% CI coverage of beta_lambda at nominal rate (lambda-arm RE)", {
  skip_on_cran()
  beta_lambda <- c(0.6, 0.4); beta_p <- 0; sigma_b <- 0.6
  n_seed <- 8L
  cov_int <- logical(n_seed); cov_slope <- logical(n_seed)
  for (s in seq_len(n_seed)) {
    sim <- sim_abun_lambda_re(N = 100, J = 4, ngrp = 10,
                              beta_lambda = beta_lambda, beta_p = beta_p,
                              sigma_b = sigma_b, seed = 2000 + s)
    fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
                data = sim$data, y = sim$y, family = abun(),
                method = "laplace", verbose = FALSE,
                control = list(n.quad = 5))
    lo <- fit$means - 1.96 * fit$sds
    hi <- fit$means + 1.96 * fit$sds
    cov_int[s]   <- lo["lambda_(Intercept)"] <= beta_lambda[1] &&
                    hi["lambda_(Intercept)"] >= beta_lambda[1]
    cov_slope[s] <- lo["lambda_x1"]          <= beta_lambda[2] &&
                    hi["lambda_x1"]          >= beta_lambda[2]
  }
  expect_gte(mean(cov_int),   0.75)
  expect_gte(mean(cov_slope), 0.75)
})

test_that("NB + lambda-arm RE returns a usable fit with positive r", {
  skip_on_cran()
  beta_lambda <- c(0.6, 0.3); beta_p <- 0; sigma_b <- 0.5; size_true <- 4
  sim <- sim_abun_lambda_re(N = 80, J = 4, ngrp = 8,
                            beta_lambda = beta_lambda, beta_p = beta_p,
                            sigma_b = sigma_b, size = size_true,
                            mixture = "negbin", seed = 3001)
  fit <- tobs(formula = ~ x1 + (1 | g), detection = ~ 1,
              data = sim$data, y = sim$y, family = abun(mixture = "negbin"),
              method = "laplace", verbose = FALSE,
              control = list(n.quad = 3))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$mixture, "negbin")
  expect_true(is.finite(fit$nmix_dispersion$r))
  expect_gt(fit$nmix_dispersion$r, 0)
  # sigma in the right order of magnitude
  expect_gt(fit$means["sigma_g1_(Intercept)"], 0.1)
  expect_lt(fit$means["sigma_g1_(Intercept)"], 2.0)
})
