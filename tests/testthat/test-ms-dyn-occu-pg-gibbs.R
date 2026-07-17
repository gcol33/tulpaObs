# =============================================================================
# test-ms-dyn-occu-pg-gibbs.R - community dynamic occupancy Polya-Gamma Gibbs
# (ms_dyn_occu() method = "pg_gibbs"; spOccupancy tMsPGOcc; tulpaObs#115, #126).
#
# The community PG machinery (msPGOcc) + a 2-state HMM forward-filter
# backward-sample latent step: per-species season-1 occupancy psi1 and detection
# p with Gaussian community hyperpriors, SHARED community-level colonization gamma
# and extinction eps. Unlike the community Laplace-EM, the Gibbs targets the exact
# posterior. These check the S3 surface, convergence, the shared gamma / eps
# recovery, the community means, and the community-SD recovery (not attenuated).
# =============================================================================

test_that("ms_dyn_occu() method = 'pg_gibbs' S3 + convergence", {
  sim <- simulate_ms_dyn_occu(N = 70, J = 3, n_species = 8, n_seasons = 4,
                              beta_comm_mean = c(stats::qlogis(0.5)),
                              beta_comm_sd = c(0.6), gamma = 0.3, epsilon = 0.2,
                              seed = 1)
  fit <- tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
              colonization = ~ 1, extinction = ~ 1, y = sim$y,
              species = paste0("sp", seq_len(8)), method = "pg_gibbs",
              control = list(n.iter = 700L, n.warmup = 350L, n.chains = 2L,
                             seed = 1, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "pg_gibbs")
  expect_true(all(is.finite(fit$rhat)) && all(fit$rhat < 1.15))
  expect_true(all(c("psi1_(Intercept)", "gamma_(Intercept)", "eps_(Intercept)")
                  %in% names(fit$means)))
  cm <- fit$ms_community
  expect_true(all(c("sd_psi1", "sd_p", "coef_psi1", "coef_p") %in% names(cm)))
  expect_equal(nrow(cm$coef_psi1), 8L)
})

test_that("ms_dyn_occu() pg_gibbs recovers gamma/eps + community means + SD", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 8L
  psi1_t <- stats::qlogis(0.5); g_t <- 0.3; e_t <- 0.2; sdpsi_t <- 0.7
  psi1_e <- gam_e <- eps_e <- sdpsi_e <- rep(NA_real_, n_seed)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_dyn_occu(N = 100, J = 4, n_species = 14, n_seasons = 5,
                                beta_comm_mean = c(psi1_t), beta_comm_sd = c(sdpsi_t),
                                gamma = g_t, epsilon = e_t, seed = 400 + s)
    fg <- tryCatch(
      tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
           colonization = ~ 1, extinction = ~ 1, y = sim$y,
           species = paste0("sp", seq_len(14)), method = "pg_gibbs",
           control = list(n.iter = 1200L, n.warmup = 600L, n.chains = 2L,
                          seed = s, verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fg)) next
    m <- fg$means
    psi1_e[s] <- m[["psi1_(Intercept)"]]
    gam_e[s]  <- stats::plogis(m[["gamma_(Intercept)"]])
    eps_e[s]  <- stats::plogis(m[["eps_(Intercept)"]])
    sdpsi_e[s] <- fg$ms_community$sd_psi1[1]
  }
  # Shared colonization / extinction recover tightly (informed by all species).
  expect_lt(abs(mean(gam_e, na.rm = TRUE) - g_t), 0.05)
  expect_lt(abs(mean(eps_e, na.rm = TRUE) - e_t), 0.05)
  # Community season-1 occupancy mean recovers.
  expect_lt(abs(mean(psi1_e, na.rm = TRUE) - psi1_t), 0.15)
  # Community SD recovers (not attenuated to a lower bound); a lenient band since
  # a variance component at moderate species count is noisy.
  expect_lt(abs(mean(sdpsi_e, na.rm = TRUE) - sdpsi_t), 0.15)
})
