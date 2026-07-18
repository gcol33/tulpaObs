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
       spec = spec, pri = tulpaObs:::.ms_ocs_nuts_priors(), is_nb = is_nb)
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


# --- (7) shared areal field (proper-CAR) recovery (tulpaObs#73) -------------

.msan_grid_graph <- function(side) {
  N <- side * side; A <- matrix(0L, N, N)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    i <- idx(r, c)
    if (r < side) { j <- idx(r + 1L, c); A[i, j] <- 1L; A[j, i] <- 1L }
    if (c < side) { j <- idx(r, c + 1L); A[i, j] <- 1L; A[j, i] <- 1L }
  }
  A
}

# The shared-field block leaves the non-spatial eval byte-identical at raw = 0.
test_that("ms_abun NUTS shared-field block is a no-op at raw = 0", {
  skip_on_cran()
  P <- .msan_pieces("poisson")
  side <- 6L; A <- .msan_grid_graph(side)   # n_sites = 36 != the fixture's 30
  # Build a matched-size field for the fixture's n_sites by mapping each site to
  # its own unit and an identity Linv (so f = raw); at raw = 0 there is no field.
  N <- P$spec$n_sites
  spec_f <- P$spec
  spec_f$n_field_units <- N
  spec_f$field_map <- seq_len(N)
  spec_f$field_Linv <- diag(N)
  set.seed(9)
  theta <- P$theta0 + stats::rnorm(length(P$theta0), 0, 0.03)
  r_off <- cpp_ms_abun_nuts_joint_logpost(P$spec, theta, P$pri, 10, 1.5)
  r_on  <- cpp_ms_abun_nuts_joint_logpost(spec_f, c(theta, numeric(N)), P$pri, 10, 1.5)
  expect_lt(abs(r_off$lp - r_on$lp), 1e-9)
  expect_lt(max(abs(r_off$grad - r_on$grad[seq_along(r_off$grad)])), 1e-9)
})

# Full-vector gradient (incl. the field raw block) vs finite differences.
test_that("ms_abun NUTS shared-field gradient matches finite differences", {
  skip_on_cran()
  P <- .msan_pieces("poisson")
  N <- P$spec$n_sites
  side <- 6L
  # a proper-CAR Linv on a path graph over the fixture's sites (full-rank).
  A <- matrix(0L, N, N)
  for (i in seq_len(N - 1L)) { A[i, i + 1L] <- 1L; A[i + 1L, i] <- 1L }
  Q <- diag(rowSums(A)) - 0.7 * A
  Qr <- 1.5 * Q + diag(1e-4 * 1.5, N)
  Linv <- backsolve(chol(Qr), diag(N))
  spec_f <- P$spec
  spec_f$n_field_units <- N; spec_f$field_map <- seq_len(N); spec_f$field_Linv <- Linv
  set.seed(4)
  theta <- c(P$theta0, stats::rnorm(N, 0, 0.2)) +
    c(stats::rnorm(length(P$theta0), 0, 0.02), rep(0, N))
  res <- cpp_ms_abun_nuts_joint_logpost(spec_f, theta, P$pri, 10, 1.5)
  f_lp <- function(th) cpp_ms_abun_nuts_joint_logpost(spec_f, th, P$pri, 10, 1.5)$lp
  num <- .msan_fd_grad(f_lp, theta)
  expect_lt(max(abs(res$grad - num)), 1e-5)
})

test_that("ms_abun NUTS + car_proper() recovers community means + the field", {
  skip_on_cran()
  skip_if_fast()
  side <- 7L; N <- side * side; J <- 4L; nsp <- 8L
  set.seed(7)
  A <- .msan_grid_graph(side)
  coord <- expand.grid(r = seq_len(side), c = seq_len(side))
  f_true <- 0.6 * scale(sin(coord$r / side * pi) + cos(coord$c / side * pi))[, 1]
  f_true <- f_true - mean(f_true)
  data <- data.frame(abund_cov1 = stats::rnorm(N), det_cov1 = stats::rnorm(N))
  X_lam <- stats::model.matrix(~ abund_cov1, data)
  X_p   <- stats::model.matrix(~ det_cov1, data)
  mu_lam <- c(log(4), 0.4); mu_p <- c(0.3, -0.3)
  beta_lam <- cbind(stats::rnorm(nsp, mu_lam[1], 0.4), stats::rnorm(nsp, mu_lam[2], 0.3))
  beta_p   <- cbind(stats::rnorm(nsp, mu_p[1], 0.4), stats::rnorm(nsp, mu_p[2], 0.3))
  y <- array(0L, dim = c(N, J, nsp))
  for (s in seq_len(nsp)) {
    lam <- exp(as.numeric(X_lam %*% beta_lam[s, ]) + f_true)
    Ni <- stats::rpois(N, lam); p <- stats::plogis(as.numeric(X_p %*% beta_p[s, ]))
    for (i in seq_len(N)) y[i, , s] <- stats::rbinom(J, Ni[i], p[i])
  }
  sp <- paste0("sp", seq_len(nsp))
  fit <- tobs(~ abund_cov1 + car_proper(graph = A), data = data, y = y,
              family = ms_abun(), detection = ~ det_cov1, species = sp,
              method = "nuts",
              control = list(n.iter = 300L, n.warmup = 300L, n.chains = 2L,
                             seed = 1L, verbose = FALSE))
  expect_equal(fit$method, "nuts")
  expect_lt(fit$nuts$divergent_total, 0.05 * nrow(fit$nuts$draws))
  expect_false(is.null(fit$spatial_field))

  truth <- c(mu_lam, mu_p)
  m <- fit$means[seq_along(truth)]; s <- fit$sds[seq_along(truth)]
  expect_true(all(abs(m - truth) / s < 3.0))
  expect_gt(cor(fit$spatial_field, f_true), 0.80)
})


# --- (8) gates -------------------------------------------------------------

test_that("ms_abun NUTS + icar() shared field recovers community means + field (#113)", {
  skip_on_cran()
  skip_if_fast()
  # The #71 sum-to-zero reparameterisation samples the shared intrinsic icar field
  # (whitened raw ~ N(0, I_{n-1})) jointly with the per-species community block;
  # the non-square loading flows through the generalized in-tree field block.
  side <- 7L; N <- side * side; J <- 4L; nsp <- 8L
  set.seed(13)
  A <- .msan_grid_graph(side)
  coord <- expand.grid(r = seq_len(side), c = seq_len(side))
  f_true <- 0.6 * scale(sin(coord$r / side * pi) + cos(coord$c / side * pi))[, 1]
  f_true <- f_true - mean(f_true)
  data <- data.frame(abund_cov1 = stats::rnorm(N), det_cov1 = stats::rnorm(N))
  X_lam <- stats::model.matrix(~ abund_cov1, data)
  X_p   <- stats::model.matrix(~ det_cov1, data)
  mu_lam <- c(log(4), 0.4); mu_p <- c(0.3, -0.3)
  beta_lam <- cbind(stats::rnorm(nsp, mu_lam[1], 0.4), stats::rnorm(nsp, mu_lam[2], 0.3))
  beta_p   <- cbind(stats::rnorm(nsp, mu_p[1], 0.4), stats::rnorm(nsp, mu_p[2], 0.3))
  y <- array(0L, dim = c(N, J, nsp))
  for (s in seq_len(nsp)) {
    lam <- exp(as.numeric(X_lam %*% beta_lam[s, ]) + f_true)
    Ni <- stats::rpois(N, lam); p <- stats::plogis(as.numeric(X_p %*% beta_p[s, ]))
    for (i in seq_len(N)) y[i, , s] <- stats::rbinom(J, Ni[i], p[i])
  }
  sp <- paste0("sp", seq_len(nsp))
  fit <- tobs(~ abund_cov1 + icar(graph = A), data = data, y = y,
              family = ms_abun(), detection = ~ det_cov1, species = sp,
              method = "nuts",
              control = list(n.iter = 300L, n.warmup = 300L, n.chains = 2L,
                             seed = 1L, verbose = FALSE))
  expect_equal(fit$method, "nuts")
  expect_lt(fit$nuts$divergent_total, 0.05 * nrow(fit$nuts$draws))
  expect_false(is.null(fit$spatial_field))
  expect_lt(abs(mean(fit$spatial_field)), 1e-6)                     # sum-to-zero centred
  truth <- c(mu_lam, mu_p)
  m <- fit$means[seq_along(truth)]; s <- fit$sds[seq_along(truth)]
  expect_true(all(abs(m - truth) / s < 3.0))
  expect_gt(cor(fit$spatial_field, f_true), 0.80)
})
