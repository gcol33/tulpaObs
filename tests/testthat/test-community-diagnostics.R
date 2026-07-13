# test-community-diagnostics.R - WAIC / DIC / CPO for the community-occupancy
# families (ms_occu / ms_dyn_occu / ms_int_occu). The per-(species, site)
# marginal is scored over the community-mean pseudo-draws with per-species BLUP
# deviations plugged in (R/community_ploglik.R).

test_that("ms_occu two-state marginal matches the C++ single-season kernel", {
  skip_on_cran()
  sim <- simulate_ms_occu(N = 60, J = 3, n_species = 6, seed = 11)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:6),
              method = "laplace", control = list(verbose = FALSE))

  ll_comm <- .tobs_ploglik_community(fit, fit$draws)
  m  <- fit$model; cm <- fit$ms_community
  p_occ <- m$process_info[[1]]$p; p_p <- m$process_info[[2]]$p
  # Rebuild species 1 as a single-season shim and score with the C++ kernel.
  draws_s <- cbind(
    sweep(fit$draws[, seq_len(p_occ), drop = FALSE], 2, cm$blup_psi[1, ], "+"),
    sweep(fit$draws[, p_occ + seq_len(p_p), drop = FALSE], 2, cm$blup_p[1, ], "+"))
  colnames(draws_s) <- colnames(fit$draws)
  ys <- m$y[, , 1]; ys[!m$valid[, , 1]] <- -1L
  shim <- list(model_type = "single", X_processes = list(m$X_occ, m$X_det),
               process_info = m$process_info, y = ys, n_sites = m$n_sites)
  ll_cpp <- .tobs_ploglik_replicated(shim, draws_s, 1L)

  expect_equal(ll_comm[, seq_len(m$n_sites)], ll_cpp, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(ncol(ll_comm), m$n_sites * m$n_species)
})

test_that("tobs_waic / tobs_dic / tobs_cpo work on ms_occu", {
  skip_on_cran()
  sim <- simulate_ms_occu(N = 60, J = 3, n_species = 6, seed = 21)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:6),
              method = "laplace", control = list(verbose = FALSE))
  w <- tobs_waic(fit); d <- tobs_dic(fit); cpo <- tobs_cpo(fit)
  expect_true(is.finite(w$waic))
  expect_true(is.finite(d$dic))
  expect_true(is.finite(cpo$lpml))
  # WAIC and DIC estimate the same expected deviance; agree to a few percent.
  expect_lt(abs(w$waic - d$dic) / w$waic, 0.05)
  # LPML (sum of per-obs CPO) matches -0.5 * WAIC scale roughly (both log-scores).
  expect_lt(cpo$lpml, 0)
})

test_that("tobs_waic / tobs_dic / tobs_cpo work on ms_int_occu", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_int_occu(N = 80, J = c(3, 4), n_species = 6, seed = 22)
  fit <- tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:6),
              method = "laplace", control = list(verbose = FALSE))
  w <- tobs_waic(fit); d <- tobs_dic(fit); cpo <- tobs_cpo(fit)
  expect_true(is.finite(w$waic))
  expect_true(is.finite(d$dic))
  expect_true(is.finite(cpo$lpml))
})

test_that("tobs_waic / tobs_dic / tobs_cpo work on ms_dyn_occu", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_dyn_occu(N = 50, J = 3, n_species = 6, n_seasons = 4,
                              gamma = 0.2, epsilon = 0.1, seed = 23)
  fit <- tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", 1:6),
              method = "laplace", control = list(verbose = FALSE))
  w <- tobs_waic(fit); d <- tobs_dic(fit); cpo <- tobs_cpo(fit)
  expect_true(is.finite(w$waic))
  expect_true(is.finite(d$dic))
  expect_true(is.finite(cpo$lpml))
})
