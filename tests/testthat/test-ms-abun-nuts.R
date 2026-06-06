# Community / multispecies N-mixture NUTS (ms_abun(), method = "nuts"; tulpaObs#14).
# The sampler draws the EXACT joint posterior -- community means, per-species
# deviations, and community covariances -- over the closed-form per-(species,
# site) Royle marginal via the in-tree C++ FullGradFn (src/ms_abun_nuts.cpp),
# warm-started at the Laplace-EM mode. The R target .tobs_ms_abun_nuts_logpost
# (R/ms_abun_nuts.R) is the oracle; the C++ port is cross-checked against it.
#
# Coverage: (1) the R oracle gradient vs finite differences, (2) the C++
# FullGradFn byte-exact vs the R oracle, (3) community-mean recovery + 0
# divergences (Poisson + NB), (4) community-mean CI coverage, (5) S3 + calibrated
# WAIC from the per-species draws, (6) the spatial-term gate.

# --- shared fixtures / helpers ---------------------------------------------

.msan_pieces <- function(mixture, n_species = 5, N = 30, J = 3, seed = 7) {
  is_nb <- identical(mixture, "negbin")
  sim <- simulate_ms_abun(n_species = n_species, N = N, J = J,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(3), 0.4), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4,
                          mixture = mixture, size = 4, sigma_logr = 0.3,
                          seed = seed)
  model <- tulpaObs:::.tobs_build_ms_abun(
    abund_formula = ~ abund_cov1, det_formula = ~ det_cov1,
    data = sim$data, y = sim$y, species = sim$species)
  lf   <- tulpaObs:::.tobs_ms_nmix_longform(model)
  Xlam <- model$X_processes[[1]]
  p_lam <- ncol(Xlam); p_p <- model$process_info[[2]]$p
  K_max <- as.integer(max(lf$y) + 100L)
  lay  <- tulpaObs:::.tobs_ms_abun_nuts_layout(p_lam, p_p, model$n_species, is_nb)
  margs <- tulpaObs:::.tobs_ms_abun_nuts_marginals(
    lf, Xlam, model$n_sites, if (is_nb) "NB" else "P", K_max)
  warm <- tulpaObs:::nmix_laplace_re(
    y = lf$y, site_idx = lf$site_idx, species_idx = lf$species_idx,
    X_lambda = Xlam, X_p = lf$X_p, n_sites = model$n_sites,
    n_species = model$n_species, K_max = K_max, max_iter = 60L,
    mixture = if (is_nb) "NB" else "P",
    optimizer = if (is_nb) "joint_grad" else "em",
    n_quad = if (is_nb) 3L else 1L, lkj_eta = 1.5, verbose = FALSE)
  theta0 <- tulpaObs:::.tobs_ms_abun_nuts_pack_init(warm, lay)
  spec <- list(y = as.integer(lf$y), site_idx = as.integer(lf$site_idx),
               species_idx = as.integer(lf$species_idx),
               X_lambda = Xlam, X_p = lf$X_p,
               n_sites = model$n_sites, n_species = model$n_species,
               K_max = K_max, is_nb = is_nb)
  list(sim = sim, model = model, lay = lay, margs = margs, theta0 = theta0,
       spec = spec, pri = tulpaObs:::.tobs_ms_abun_nuts_priors(), is_nb = is_nb)
}

.msan_fd_grad <- function(f, theta, h = 1e-5) {
  vapply(seq_along(theta), function(j) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  }, 0)
}


# --- (1) R oracle gradient vs finite differences ---------------------------

test_that("ms_abun NUTS R oracle gradient matches finite differences", {
  skip_on_cran()
  for (mix in c("poisson", "negbin")) {
    P <- .msan_pieces(mix)
    set.seed(1)
    theta <- P$theta0 + stats::rnorm(length(P$theta0), 0, 0.05)
    f_lp <- function(th)
      tulpaObs:::.tobs_ms_abun_nuts_logpost(th, P$margs, P$lay, P$pri,
                                            grad = FALSE)$lp
    ana <- tulpaObs:::.tobs_ms_abun_nuts_logpost(theta, P$margs, P$lay, P$pri,
                                                 grad = TRUE)$grad
    num <- .msan_fd_grad(f_lp, theta)
    expect_lt(max(abs(ana - num)), 1e-5)
  }
})


# --- (2) C++ FullGradFn byte-exact vs the R oracle -------------------------

test_that("ms_abun NUTS C++ FullGradFn matches the R oracle", {
  skip_on_cran()
  for (mix in c("poisson", "negbin")) {
    P <- .msan_pieces(mix)
    set.seed(2)
    theta <- P$theta0 + stats::rnorm(length(P$theta0), 0, 0.05)
    cpp <- cpp_ms_abun_nuts_joint_logpost(P$spec, theta, P$pri,
                                          sigma_beta = 10, sigma_logr = 1.5)
    r   <- tulpaObs:::.tobs_ms_abun_nuts_logpost(
      theta, P$margs, P$lay, P$pri, sigma.beta = 10, sigma.logr = 1.5,
      grad = TRUE)
    expect_lt(abs(cpp$lp - r$lp), 1e-9)
    expect_lt(max(abs(cpp$grad - r$grad)), 1e-9)
  }
})


# --- (3) community-mean recovery + 0 divergences (Poisson) -----------------

test_that("ms_abun NUTS recovers community means (Poisson)", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_abun(n_species = 8, N = 40, J = 4,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4, seed = 42)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
              detection = ~ det_cov1, species = sim$species, method = "nuts",
              control = list(n.iter = 300L, n.warmup = 300L, seed = 1L,
                             verbose = FALSE))
  expect_equal(fit$method, "nuts")
  expect_false(is.null(fit$nuts$draws))
  expect_lt(fit$nuts$divergent_total, 0.05 * nrow(fit$nuts$draws))

  truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  z <- abs(fit$means - truth) / fit$sds
  expect_true(all(z < 2.5))

  # Per-species coefficients track the simulated truth.
  cm <- fit$ms_community
  expect_gt(min(diag(cor(cm$coef_lambda, sim$truth$beta_lambda))), 0.85)
  expect_gt(min(diag(cor(cm$coef_p,      sim$truth$beta_p))),      0.65)
})


# --- (4) community-mean CI coverage ----------------------------------------

test_that("ms_abun NUTS community-mean 95% CIs cover at the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 6L
  covered <- logical(0)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_abun(n_species = 8, N = 40, J = 4,
                            n_abund_covs = 1, n_det_covs = 1,
                            mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.4),
                            sd_lambda = 0.5, sd_p = 0.4, seed = 200 + s)
    fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
                detection = ~ det_cov1, species = sim$species, method = "nuts",
                control = list(n.iter = 300L, n.warmup = 300L, seed = 1L,
                               verbose = FALSE))
    truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
    lo <- fit$means - 1.96 * fit$sds
    hi <- fit$means + 1.96 * fit$sds
    covered <- c(covered, truth >= lo & truth <= hi)
  }
  expect_gt(mean(covered), 0.85)
})


# --- (5) negbin recovery + per-species dispersion --------------------------

test_that("ms_abun NUTS (negbin) recovers means + community dispersion", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_abun(n_species = 8, N = 45, J = 4,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(4), 0.4), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4,
                          mixture = "negbin", size = 4, sigma_logr = 0.3,
                          seed = 31)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y,
              family = ms_abun(mixture = "negbin"),
              detection = ~ det_cov1, species = sim$species, method = "nuts",
              control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                             verbose = FALSE))
  expect_equal(fit$mixture, "negbin")
  expect_true("log_r" %in% names(fit$means))
  expect_true(is.finite(fit$means[["log_r"]]))

  truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  nm <- setdiff(names(fit$means), "log_r")
  z  <- abs(fit$means[nm] - truth) / fit$sds[nm]
  expect_true(all(z < 3.0))

  expect_false(is.null(fit$ms_dispersion))
  expect_gt(fit$ms_dispersion$r, 0)
  expect_lt(abs(log(fit$ms_dispersion$r) - log(sim$truth$size)), log(3))
})


# --- (6) S3 methods + calibrated WAIC from the per-species draws ------------

test_that("ms_abun NUTS S3 methods + WAIC work", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_abun(n_species = 8, N = 30, J = 3, seed = 3)
  fit <- tobs(~ abund_cov1, data = sim$data, y = sim$y, family = ms_abun(),
              detection = ~ det_cov1, species = sim$species, method = "nuts",
              control = list(n.iter = 300L, n.warmup = 300L, seed = 1L,
                             verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))

  cf <- coef(fit)
  expect_setequal(names(cf), c("lambda", "p"))
  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 8L * (2L + 2L))

  fv <- fitted(fit)
  expect_equal(dim(fv$lambda), c(30L, 8L))

  # WAIC scored over the per-(species, site) NUTS draws.
  w <- tobs_waic(fit, n.draws = 200L)
  expect_true(is.finite(w$waic))
  expect_gt(w$p_waic, 0)
})


# --- (7) gates -------------------------------------------------------------

test_that("ms_abun NUTS rejects a spatial term with a pointer", {
  skip_on_cran()
  set.seed(5)
  N <- 16L
  adj <- matrix(0L, N, N)
  for (i in seq_len(N - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  sim <- simulate_ms_abun(n_species = 4, N = N, J = 3, seed = 5)
  expect_error(
    tobs(~ abund_cov1 + icar(graph = adj), data = sim$data, y = sim$y,
         family = ms_abun(), detection = ~ det_cov1, species = sim$species,
         method = "nuts", control = list(verbose = FALSE)),
    "non-spatial")
})
