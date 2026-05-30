# Spatial community / multispecies N-mixture (ms_abun() + a shared areal field on
# the abundance arm; the spAbundance sfMsNMix analogue, gcol33/tulpaObs#12). The
# fit is the in-tree nested Laplace-EM: per-species coefficient RE with Gaussian
# community covariances PLUS one shared ICAR / BYM2 / proper-CAR field on
# log lambda, integrated over the field-hyperparameter (and NB size) grid.
#
# These fits loop S species x an EM x an outer grid, so each is ~1-3 min: every
# block here is skip_on_cran(), and the multi-seed / extra-engine blocks add
# skip_if_fast().

# Dense rook (4-neighbour) adjacency for a g x g grid.
rook_adj <- function(g) {
  n <- g * g
  A <- matrix(0L, n, n)
  idx <- function(r, c) (r - 1L) * g + c
  for (r in seq_len(g)) for (c in seq_len(g)) {
    i <- idx(r, c)
    if (r > 1) A[i, idx(r - 1L, c)] <- 1L
    if (r < g) A[i, idx(r + 1L, c)] <- 1L
    if (c > 1) A[i, idx(r, c - 1L)] <- 1L
    if (c < g) A[i, idx(r, c + 1L)] <- 1L
  }
  A
}


test_that("spatial ms_abun (ICAR) recovers community means and the shared field", {
  skip_on_cran()
  adj <- rook_adj(7L)
  # 20 species: a community covariance is a variance component, so it needs
  # enough groups to be identified -- the non-spatial recovery fixture uses 14,
  # and at fewer the Laplace-EM attenuates / collapses the weakest coordinate
  # (the detection intercept here), as for any small-group variance-component fit.
  sim <- simulate_ms_abun(n_species = 20, J = 5,
                          n_abund_covs = 1, n_det_covs = 1,
                          mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.3),
                          sd_lambda = 0.5, sd_p = 0.4,
                          graph = adj, sigma.field = 0.6, seed = 7)
  fit <- tobs(~ abund_cov1 + icar(graph = adj),
              detection = ~ det_cov1, family = ms_abun(),
              data = sim$data, y = sim$y, species = sim$species,
              method = "nested_laplace", control = list(verbose = FALSE))

  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "nested_laplace")
  expect_true(isTRUE(fit$convergence$converged))

  # Community means within ~3 SE of truth on every coordinate.
  truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
  z <- abs(fit$means - truth) / fit$sds
  expect_true(all(z < 3))

  # Field SHAPE recovery: the posterior-mean field tracks the simulated one.
  expect_length(fit$spatial_field, nrow(adj))
  expect_gt(cor(fit$spatial_field, sim$truth$field), 0.8)

  # Community covariances track the realized per-species spread (not -> 0).
  cm <- fit$ms_community
  emp_sd <- c(apply(sim$truth$beta_lambda, 2, sd), apply(sim$truth$beta_p, 2, sd))
  est_sd <- c(cm$sd_lambda, cm$sd_p)
  expect_true(all(est_sd > 0.4 * emp_sd & est_sd < 1.8 * emp_sd))
})


test_that("spatial ms_abun S3 surface carries the field", {
  skip_on_cran()
  adj <- rook_adj(5L)
  sim <- simulate_ms_abun(n_species = 6, J = 4, graph = adj,
                          sigma.field = 0.5, seed = 3)
  fit <- tobs(~ abund_cov1 + icar(graph = adj),
              detection = ~ det_cov1, family = ms_abun(),
              data = sim$data, y = sim$y, species = sim$species,
              method = "nested_laplace", control = list(verbose = FALSE))

  expect_no_error(print(fit))
  cf <- coef(fit)
  expect_true(is.list(cf))
  expect_setequal(names(cf), c("lambda", "p"))
  expect_setequal(names(cf$lambda), c("(Intercept)", "abund_cov1"))
  expect_setequal(names(cf$p), c("(Intercept)", "det_cov1"))
  expect_equal(nrow(vcov(fit)), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 6L * (2L + 2L))

  # fitted lambda includes the shared field offset (one unit per site).
  fv <- fitted(fit)
  expect_equal(dim(fv$lambda), c(nrow(adj), 6L))
  expect_true(all(fv$lambda > 0))

  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys), c(nrow(adj), 4L, 6L))

  # Spatial hyperparameter posterior is reported.
  expect_true(!is.null(fit$ms_hyper$tau))
})


test_that("spatial ms_abun (BYM2) fits and recovers the field shape", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(6L)
  sim <- simulate_ms_abun(n_species = 8, J = 5, graph = adj,
                          sigma.field = 0.6, seed = 21)
  fit <- tobs(~ abund_cov1 + bym2(graph = adj),
              detection = ~ det_cov1, family = ms_abun(),
              data = sim$data, y = sim$y, species = sim$species,
              method = "nested_laplace", control = list(verbose = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_length(fit$spatial_field, nrow(adj))
  expect_gt(cor(fit$spatial_field, sim$truth$field), 0.7)
  expect_true(!is.null(fit$ms_hyper$sigma) && !is.null(fit$ms_hyper$rho))
})


test_that("spatial ms_abun (proper CAR) fits and recovers the field shape", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(6L)
  sim <- simulate_ms_abun(n_species = 8, J = 5, graph = adj,
                          sigma.field = 0.6, seed = 22)
  fit <- tobs(~ abund_cov1 + car_proper(graph = adj),
              detection = ~ det_cov1, family = ms_abun(),
              data = sim$data, y = sim$y, species = sim$species,
              method = "nested_laplace", control = list(verbose = FALSE))
  expect_identical(fit$method, "nested_laplace")
  expect_length(fit$spatial_field, nrow(adj))
  expect_gt(cor(fit$spatial_field, sim$truth$field), 0.7)
  expect_true(!is.null(fit$ms_hyper$rho))
})


test_that("spatial ms_abun (negbin) integrates the size r over the grid", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(6L)
  sim <- simulate_ms_abun(n_species = 8, J = 5, graph = adj, sigma.field = 0.5,
                          mu_lambda = c(log(4), 0.4), mu_p = c(0.3, -0.3),
                          mixture = "negbin", size = 4, seed = 41)
  fit <- tobs(~ abund_cov1 + icar(graph = adj),
              detection = ~ det_cov1, family = ms_abun(mixture = "negbin"),
              data = sim$data, y = sim$y, species = sim$species,
              method = "nested_laplace", control = list(verbose = FALSE))
  expect_equal(fit$mixture, "negbin")
  # r is grid-integrated (a hyperparameter, not a model coordinate): no log_r
  # column, summarized in ms_dispersion / ms_hyper instead.
  expect_false("log_r" %in% names(fit$means))
  expect_false(is.null(fit$ms_dispersion))
  expect_true(is.finite(fit$ms_dispersion$r) && fit$ms_dispersion$r > 0)
  # Within a factor of ~3 of truth (grid-integrated r is coarse but ballpark).
  expect_lt(abs(log(fit$ms_dispersion$r) - log(sim$truth$size)), log(3))
  expect_gt(cor(fit$spatial_field, sim$truth$field), 0.7)
})


test_that("spatial ms_abun community-mean 95% CIs cover near the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  adj <- rook_adj(5L)
  n_seed <- 12L
  covered <- logical(0)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_abun(n_species = 8, J = 5,
                            n_abund_covs = 1, n_det_covs = 1,
                            mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.3),
                            sd_lambda = 0.5, sd_p = 0.4,
                            graph = adj, sigma.field = 0.5, seed = 200 + s)
    fit <- tobs(~ abund_cov1 + icar(graph = adj),
                detection = ~ det_cov1, family = ms_abun(),
                data = sim$data, y = sim$y, species = sim$species,
                method = "nested_laplace", control = list(verbose = FALSE))
    truth <- c(sim$truth$mu_lambda, sim$truth$mu_p)
    lo <- fit$means - 1.96 * fit$sds
    hi <- fit$means + 1.96 * fit$sds
    covered <- c(covered, truth >= lo & truth <= hi)
  }
  # Nominal 95%; Monte-Carlo slack on 12 x 4 = 48 intervals.
  expect_gt(mean(covered), 0.8)
})

test_that("spatial ms_abun interops with spAbundance::sfMsNMix (smoke)", {
  skip_on_cran()
  skip_if_fast()
  skip_if_not_installed("spAbundance")
  # SMOKE only: a properly-converged sfMsNMix chain runs for hours, so this fits
  # a tiny chain purely to confirm the interop (data simulated here feeds
  # sfMsNMix, both return a finite community abundance intercept). It does NOT
  # assert numerical agreement -- that needs the long offline benchmark in
  # dev_notes/probe_ms_abun_spatial_vs_spabundance.R.
  adj <- rook_adj(5L); N <- nrow(adj)
  sim <- simulate_ms_abun(n_species = 6, J = 4, graph = adj, sigma.field = 0.4,
                          mu_lambda = c(log(4), 0.5), mu_p = c(0.3, -0.3), seed = 9)
  J <- dim(sim$y)[2]

  fit <- tobs(~ abund_cov1 + icar(graph = adj), detection = ~ det_cov1,
              family = ms_abun(), data = sim$data, y = sim$y,
              species = sim$species, method = "nested_laplace",
              control = list(verbose = FALSE))
  expect_true(is.finite(fit$means[["lambda_(Intercept)"]]))

  coords   <- cbind((seq_len(N) - 1L) %/% 5L, (seq_len(N) - 1L) %% 5L) + 0.0
  y_sp     <- aperm(sim$y, c(3, 1, 2))                   # species x site x rep
  det_covs <- list(det_cov1 = matrix(rep(sim$data$det_cov1, J), N, J))
  out <- tryCatch(
    spAbundance::sfMsNMix(
      abund.formula = ~ abund_cov1, det.formula = ~ det_cov1,
      data = list(y = y_sp,
                  abund.covs = data.frame(abund_cov1 = sim$data$abund_cov1),
                  det.covs = det_covs, coords = coords),
      n.factors = 1, cov.model = "exponential", NNGP = TRUE, n.neighbors = 5,
      n.batch = 10, batch.length = 25, n.burn = 100, n.thin = 1, n.chains = 1,
      verbose = FALSE),
    error = function(e) skip(paste("sfMsNMix unavailable:", conditionMessage(e))))

  # Plumbing: sfMsNMix returns community abundance-coefficient samples of the
  # right shape (intercept + 1 slope), both finite.
  expect_s3_class(out, "sfMsNMix")
  expect_equal(ncol(out$beta.comm.samples), 2L)
  expect_true(is.finite(mean(out$beta.comm.samples[, 1])))
})
