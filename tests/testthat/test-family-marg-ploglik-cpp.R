# Batched pointwise log-likelihood for the count / multistate families whose
# per-site marginal was already C++ (compute_*_site) but whose per-draw loop was
# R: N-mixture, removal, distance, false-positive occupancy, open N-mixture
# (dyn_abun). The kernels loop draws over the SAME per-site kernel, so they are
# BYTE-IDENTICAL to the former R loop (the oracle here, via the family marginal
# builders) and thread-count invariant -- a divergence would be a bug, both
# evaluate the identical marginal.

.eq_nanaware <- function(a, b) all((a == b) | (is.nan(a) & is.nan(b)) |
                                   (is.na(a) & is.na(b)))

# Reports the magnitude of any divergence. A bare TRUE/FALSE cannot tell a
# last-ULP difference (the two call sites compile to different floating-point
# contractions) from a structural one, and the two call for opposite responses.
expect_bitwise <- function(R, C) {
  if (.eq_nanaware(R, C)) return(testthat::succeed())
  d <- abs(R - C)
  d <- d[is.finite(d)]
  testthat::fail(sprintf(
    "not bitwise equal: %d of %d cells differ, max abs %.3e, max rel %.3e",
    sum(d > 0), length(R), max(d), max(d / pmax(abs(R[is.finite(R - C)]), 1e-300))))
}

test_that("N-mixture batched ploglik == R loop (Poisson + NB)", {
  set.seed(71)
  for (nb in c(FALSE, TRUE)) {
    n_sites <- 35L; S <- 25L
    nv <- pmax(1L, rpois(n_sites, 3)); site_idx <- rep(seq_len(n_sites), nv)
    n_obs <- length(site_idx)
    model <- list(model_type = "nmix", y_long = as.integer(rpois(n_obs, 2)),
                  site_idx = site_idx, n_sites = n_sites,
                  X_processes = list(cbind(1, rnorm(n_sites)), cbind(1, rnorm(n_obs))))
    nc <- 4L + (if (nb) 1L else 0L); draws <- matrix(rnorm(S * nc, 0, .5), S, nc)
    if (nb) colnames(draws) <- c(rep("", 4), "log_r")
    marg <- nmix_site_marginal(model$y_long, model$site_idx, model$X_processes[[1]],
                               model$X_processes[[2]], mixture = if (nb) "NB" else "P",
                               K_max = as.integer(max(model$y_long) + 100L))
    R <- matrix(0, S, n_sites)
    for (s in seq_len(S)) { r <- if (nb) exp(draws[s, 5L]) else Inf
      R[s, ] <- marg$eval_beta(draws[s, 1:2], draws[s, 3:4], r = r)$log_lik_site }
    C1 <- .tobs_ploglik_nmix(model, draws, 1L)
    expect_bitwise(R, C1)
    expect_identical(.tobs_ploglik_nmix(model, draws, 4L), C1)
  }
})

test_that("removal batched ploglik == R loop", {
  set.seed(72)
  n_sites <- 35L; S <- 25L; np <- 3L
  site_idx <- rep(seq_len(n_sites), each = np); n_obs <- length(site_idx)
  model <- list(model_type = "removal", y_long = as.integer(rpois(n_obs, 1.5)),
                site_idx = site_idx, n_sites = n_sites,
                X_processes = list(cbind(1, rnorm(n_sites)), cbind(1, rnorm(n_obs))))
  draws <- matrix(rnorm(S * 4, 0, .5), S, 4)
  marg <- .tobs_removal_nuts_marginal(model, mixture = "P")
  R <- matrix(0, S, n_sites)
  for (s in seq_len(S)) R[s, ] <- marg$eval_beta(draws[s, 1:2], draws[s, 3:4], r = Inf)$log_lik_site
  C1 <- .tobs_ploglik_removal(model, draws, 1L)
  expect_bitwise(R, C1)
  expect_identical(.tobs_ploglik_removal(model, draws, 4L), C1)
})

test_that("false-positive occupancy batched ploglik == R loop", {
  set.seed(73)
  n_sites <- 40L; S <- 25L; J <- 4L
  site_idx <- rep(seq_len(n_sites), each = J); n_obs <- length(site_idx)
  Xs <- replicate(4, cbind(1, rnorm(n_sites)), simplify = FALSE)
  model <- list(model_type = "fp_occu", y_long = as.integer(rbinom(n_obs, 1L, 0.3)),
                site_idx = site_idx, n_sites = n_sites, X_processes = Xs,
                process_info = lapply(Xs, function(X) list(p = ncol(X))))
  draws <- matrix(rnorm(S * 8, 0, .5), S, 8)
  lay <- .tobs_fp_occu_nuts_layout(2L, 2L, 2L, 2L); marg <- .tobs_fp_occu_nuts_marginal(model)
  R <- matrix(0, S, n_sites)
  for (s in seq_len(S)) R[s, ] <- marg$eval_beta(draws[s, lay$psi], draws[s, lay$p11],
                                                 draws[s, lay$p10], draws[s, lay$b])$log_lik_site
  C1 <- .tobs_ploglik_fp_occu(model, draws, 1L)
  expect_bitwise(R, C1)
  expect_identical(.tobs_ploglik_fp_occu(model, draws, 4L), C1)
})

test_that("distance batched ploglik == R loop", {
  set.seed(7)
  n_sites <- 25L; n_bins <- 3L; S <- 15L
  model <- list(model_type = "distance", y = matrix(rpois(n_sites * n_bins, 4), n_sites, n_bins),
                X_processes = list(cbind(1, rnorm(n_sites) * 0.3), cbind(1, rnorm(n_sites) * 0.3)),
                process_info = list(list(p = 2L), list(p = 2L)), n_sites = n_sites,
                key = "halfnorm", transect = "line", quad_order = 20L, cutpoints = c(0, 5, 10, 15))
  draws <- matrix(rnorm(S * 4, 0, 0.25), S, 4); draws[, 3] <- draws[, 3] + 2.3
  marg <- .tobs_distance_nuts_marginal(model, mixture = "P")
  R <- matrix(0, S, n_sites)
  for (s in seq_len(S)) R[s, ] <- marg$eval_beta(draws[s, 1:2], draws[s, 3:4], eta_b = 0, r = Inf)$log_lik_site
  C1 <- .tobs_ploglik_distance(model, draws, 1L)
  expect_bitwise(R, C1)
  expect_identical(.tobs_ploglik_distance(model, draws, 4L), C1)
})

test_that("open N-mixture (dyn_abun) batched ploglik == R loop", {
  set.seed(75)
  n_sites <- 20L; T <- 3L; J <- 2L; S <- 15L; K <- 20L
  Xs <- replicate(4, cbind(1, rnorm(n_sites)), simplify = FALSE)
  model <- list(model_type = "dyn_abun", y_flat = as.integer(rpois(n_sites * T * J, 1)),
                n_sites = n_sites, n_seasons = T, max_visits = J, K_max = K,
                mixture = "poisson", X_processes = Xs,
                process_info = lapply(Xs, function(X) list(p = ncol(X))))
  draws <- matrix(rnorm(S * 8, 0, .3), S, 8)
  lay <- .tobs_dyn_abun_nuts_layout(2L, 2L, 2L, 2L); marg <- .tobs_dyn_abun_nuts_marginal(model)
  R <- matrix(0, S, n_sites)
  for (s in seq_len(S)) R[s, ] <- marg$eval_beta(draws[s, lay$lambda], draws[s, lay$p],
                                                 draws[s, lay$omega], draws[s, lay$gamma])$log_lik_site
  C1 <- .tobs_ploglik_dyn_abun(model, draws, 1L)
  expect_bitwise(R, C1)
  expect_identical(.tobs_ploglik_dyn_abun(model, draws, 4L), C1)
})
