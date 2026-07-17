# =============================================================================
# test-ms-occu-pg-gibbs.R - community occupancy Polya-Gamma Gibbs
# (ms_occu() method = "pg_gibbs"; spOccupancy msPGOcc; gcol33/tulpaObs#115, #126).
#
# The hierarchical PG-Gibbs: per-species PG-augmented conjugate coefficient
# updates + conjugate community means + Inverse-Gamma community variances. Unlike
# the community Laplace-EM (R/community_em.R), whose variance components carry
# small-cluster attenuation (a documented lower bound), the Gibbs targets the
# exact posterior, so the community VARIANCE recovers rather than shrinking. These
# check the S3 surface, convergence, community-mean recovery, per-species
# coefficient recovery, and -- the point of the engine -- community-SD recovery
# that is not attenuated the way the Laplace-EM is.
# =============================================================================

.mspg_sim <- function(seed, ns = 18) {
  simulate_ms_occu(N = 110, J = 5, n_species = ns,
                   beta_comm_mean = c(stats::qlogis(0.45), 0.7),
                   beta_comm_sd = c(0.6, 0.4),
                   alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5), seed = seed)
}

test_that("ms_occu() method = 'pg_gibbs' S3 + convergence", {
  sim <- .mspg_sim(1, ns = 12)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", seq_len(12)), method = "pg_gibbs",
              control = list(n.iter = 1200L, n.warmup = 600L, n.chains = 2L,
                             seed = 1, verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_identical(fit$method, "pg_gibbs")
  expect_true(all(is.finite(fit$rhat)) && all(fit$rhat < 1.1))
  cm <- fit$ms_community
  expect_true(all(c("sd_psi", "sd_p", "coef_psi", "coef_p") %in% names(cm)))
  expect_equal(dim(cm$coef_psi), c(12L, 2L))
  expect_true(all(cm$sd_psi > 0) && all(cm$sd_p > 0))
})

test_that("ms_occu() pg_gibbs recovers community means, per-species coefs, and SD", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 10L
  mu_t <- c(stats::qlogis(0.45), 0.7, 0.2)          # psi int / psi slope / p int
  sdpsi_t <- c(0.6, 0.4); sdp_t <- 0.5
  mu_e <- matrix(NA_real_, n_seed, 3L)
  sdpsi_e <- matrix(NA_real_, n_seed, 2L); sdp_e <- rep(NA_real_, n_seed)
  cor_psi <- rep(NA_real_, n_seed)
  sd_pg <- sd_lap <- rep(NA_real_, n_seed)          # sd_psi slope: pg vs laplace
  for (s in seq_len(n_seed)) {
    sim <- .mspg_sim(300 + s, ns = 18)
    sp <- paste0("sp", seq_len(18))
    fg <- tryCatch(tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
                        y = sim$y, species = sp, method = "pg_gibbs",
                        control = list(n.iter = 1800L, n.warmup = 900L,
                                       n.chains = 2L, seed = s, verbose = FALSE)),
                   error = function(e) NULL)
    fl <- tryCatch(tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
                        y = sim$y, species = sp, method = "laplace",
                        control = list(verbose = FALSE)), error = function(e) NULL)
    if (is.null(fg)) next
    mu_e[s, ] <- fg$means
    sdpsi_e[s, ] <- fg$ms_community$sd_psi; sdp_e[s] <- fg$ms_community$sd_p
    cor_psi[s] <- min(diag(cor(fg$ms_community$coef_psi, sim$truth$beta_species)))
    sd_pg[s] <- fg$ms_community$sd_psi[2]
    if (!is.null(fl)) sd_lap[s] <- fl$ms_community$sd_psi[2]
  }
  # Community means recover in aggregate.
  expect_lt(abs(mean(mu_e[, 1], na.rm = TRUE) - mu_t[1]), 0.12)
  expect_lt(abs(mean(mu_e[, 2], na.rm = TRUE) - mu_t[2]), 0.12)
  expect_lt(abs(mean(mu_e[, 3], na.rm = TRUE) - mu_t[3]), 0.12)
  # Per-species occupancy coefficients recover (the occupancy BLUPs are
  # detection-filtered and shrunk, so a looser floor than a count-informed arm --
  # the same property the community Laplace test floors at 0.78).
  expect_gt(mean(cor_psi, na.rm = TRUE), 0.75)
  # Community SDs recover (NOT attenuated to a lower bound).
  expect_lt(abs(mean(sdpsi_e[, 1], na.rm = TRUE) - sdpsi_t[1]), 0.12)
  expect_lt(abs(mean(sdpsi_e[, 2], na.rm = TRUE) - sdpsi_t[2]), 0.12)
  expect_lt(abs(mean(sdp_e, na.rm = TRUE) - sdp_t), 0.12)
  # The PG community SD is systematically larger than the attenuated Laplace-EM
  # SD (the whole point: the Gibbs removes the small-cluster attenuation).
  expect_gt(mean(sd_pg - sd_lap, na.rm = TRUE), 0)
})
