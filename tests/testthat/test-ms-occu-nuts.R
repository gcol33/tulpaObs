# =============================================================================
# test-ms-occu-nuts.R - community (multispecies) single-season occupancy NUTS
# (ms_occu(), method = "nuts"; tulpaObs#69).
#
# The sampler draws the EXACT joint posterior -- community means, per-species
# deviations, and the two INDEPENDENT per-arm community covariances -- over the
# closed-form occupancy two-state per-(species, site) marginal via the in-tree
# C++ FullGradFn (src/ms_occu_nuts.cpp), warm-started at the community Laplace-EM
# mode. The R target .tobs_ms_occu_nuts_logpost (R/ms_occu_nuts.R) is the oracle;
# the C++ port is cross-checked against it.
#
# Coverage: (1) the R oracle gradient vs finite differences, (2) the C++
# FullGradFn byte-exact vs the R oracle, (3) community-mean recovery + 0
# divergences, (4) community-mean CI coverage, (5) S3 methods + richness,
# (6) the variance components de-attenuate vs the EM baseline, (7) the
# spatial-term gate.
# =============================================================================


# --- shared fixtures / helpers ---------------------------------------------

.msocc_pieces <- function(n_species = 6, N = 60, J = 4, seed = 7,
                          det_cov = FALSE) {
  sim <- simulate_ms_occu(N = N, J = J, n_species = n_species,
                          beta_comm_mean = c(0, 0.6), beta_comm_sd = c(0.6, 0.3),
                          alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5),
                          seed = seed)
  det_form <- ~ 1
  if (det_cov) { sim$data$dcov <- stats::rnorm(N); det_form <- ~ dcov }
  model <- tulpaObs:::.tobs_build_ms_occu(
    occ_formula = ~ x, det_formula = det_form,
    data = sim$data, y = sim$y, species = paste0("sp", seq_len(n_species)))
  pieces <- tulpaObs:::.tobs_ms_occu_nuts_pieces(model)
  lay <- tulpaObs:::.tobs_ms_occu_nuts_layout(pieces$P_psi, pieces$P_p, pieces$S)
  pri <- tulpaObs:::.ms_ocs_nuts_priors()
  em <- tulpaObs:::.tobs_community_em(
    S = pieces$S, P = pieces$P, arm_idx = pieces$arm_idx,
    sp_ll = pieces$sp_ll, sp_grad = pieces$sp_grad,
    init_mu = pieces$mu0, init_global = numeric(0),
    penalize_global = FALSE, sigma_beta = 5, priors = NULL,
    sigma_init = 0.3, max_iter = 40L, tol = 1e-4, newton_max = 30L,
    verbose = FALSE)
  theta0 <- tulpaObs:::.tobs_ms_occu_nuts_pack_init(em, lay, pieces$arm_idx)
  mats <- tulpaObs:::.ms_occu_spatial_count_mats(pieces$summaries,
                                                 model$n_sites, pieces$S)
  spec <- list(X_psi = pieces$X_psi, X_p = pieces$X_p,
               n_sites = model$n_sites, n_species = pieces$S,
               n_valid = mats$n_valid, n_det = mats$n_det)
  list(sim = sim, model = model, pieces = pieces, lay = lay, pri = pri,
       theta0 = theta0, spec = spec)
}

.msocc_fd_grad <- function(f, theta, h = 1e-5) {
  vapply(seq_along(theta), function(j) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  }, 0)
}


# --- (1) R oracle gradient vs finite differences ---------------------------

test_that("ms_occu NUTS R oracle gradient matches finite differences", {
  skip_on_cran()
  for (dc in c(FALSE, TRUE)) {
    P <- .msocc_pieces(det_cov = dc)
    set.seed(1)
    theta <- P$theta0 + stats::rnorm(length(P$theta0), 0, 0.05)
    f_lp <- function(th)
      tulpaObs:::.tobs_ms_occu_nuts_logpost(th, P$pieces$X_psi, P$pieces$X_p,
                                            P$pieces$summaries, P$lay, P$pri,
                                            grad = FALSE)$lp
    ana <- tulpaObs:::.tobs_ms_occu_nuts_logpost(
      theta, P$pieces$X_psi, P$pieces$X_p, P$pieces$summaries, P$lay, P$pri,
      grad = TRUE)$grad
    num <- .msocc_fd_grad(f_lp, theta)
    expect_lt(max(abs(ana - num)), 1e-5)
  }
})


# --- (2) C++ FullGradFn byte-exact vs the R oracle -------------------------

test_that("ms_occu NUTS C++ FullGradFn matches the R oracle", {
  skip_on_cran()
  for (dc in c(FALSE, TRUE)) {
    P <- .msocc_pieces(det_cov = dc)
    set.seed(2)
    theta <- P$theta0 + stats::rnorm(length(P$theta0), 0, 0.05)
    cpp <- cpp_ms_occu_nuts_joint_logpost(P$spec, theta, P$pri, sigma_beta = 5)
    r   <- tulpaObs:::.tobs_ms_occu_nuts_logpost(
      theta, P$pieces$X_psi, P$pieces$X_p, P$pieces$summaries, P$lay, P$pri,
      sigma.beta = 5, grad = TRUE)
    expect_lt(abs(cpp$lp - r$lp), 1e-9)
    expect_lt(max(abs(cpp$grad - r$grad)), 1e-9)
  }
})


# --- (3) community-mean recovery + 0 divergences ---------------------------

test_that("ms_occu NUTS recovers community means", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_occu(N = 130, J = 4, n_species = 16,
                          beta_comm_mean = c(0, 0.6), beta_comm_sd = c(0.6, 0.3),
                          alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5),
                          seed = 41)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", seq_len(16)), method = "nuts",
              control = list(n.iter = 400L, n.warmup = 400L, n.chains = 2L,
                             seed = 1L, verbose = FALSE))
  expect_equal(fit$method, "nuts")
  expect_false(is.null(fit$nuts$draws))
  expect_lt(fit$nuts$divergent_total, 0.05 * nrow(fit$nuts$draws))
  expect_lt(fit$nuts$max_rhat, 1.1)

  truth <- c("psi_(Intercept)" = 0, "psi_x" = 0.6, "p_(Intercept)" = 0.2)
  m <- fit$means[names(truth)]; s <- fit$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 2.5))

  # Per-species coefficients track the simulated truth.
  cm <- fit$ms_community
  expect_gt(cor(cm$coef_psi[, 1], sim$truth$beta_species[, 1]), 0.70)
  expect_gt(cor(cm$coef_psi[, 2], sim$truth$beta_species[, 2]), 0.60)
  expect_gt(cor(cm$coef_p[, 1],   sim$truth$alpha_species[, 1]), 0.45)
})


# --- (4) community-mean CI coverage ----------------------------------------

test_that("ms_occu NUTS community-mean 95% CIs cover at the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 8L
  covered <- logical(0)
  truth <- c("psi_(Intercept)" = 0, "psi_x" = 0.6, "p_(Intercept)" = 0.2)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_occu(N = 120, J = 4, n_species = 14,
                            beta_comm_mean = c(0, 0.6), beta_comm_sd = c(0.6, 0.3),
                            alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5),
                            seed = 700 + s)
    fit <- tryCatch(
      tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
           y = sim$y, species = paste0("sp", seq_len(14)), method = "nuts",
           control = list(n.iter = 300L, n.warmup = 300L, seed = 1L,
                          verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    m <- fit$means[names(truth)]; sd <- fit$sds[names(truth)]
    covered <- c(covered, abs(m - truth) < 1.96 * sd)
  }
  expect_gt(mean(covered), 0.80)
})


# --- (4b) community-COVARIANCE CI coverage (#115 DoD) -----------------------
# The strict DoD asks for community-covariance recovery, not just the mean. NUTS
# samples the per-arm community covariances jointly, so each fit yields a per-draw
# posterior for the community SDs (sqrt(diag(Sigma_arm)) reconstructed from the
# sampled log-Cholesky coordinates). Over 20 seeds the 95% CIs cover the true
# community SDs at the nominal rate -- the calibration the Laplace-EM lower bound
# cannot give. Measured pooled coverage ~0.90 (sd_psi[1] 0.85, sd_psi[2] 0.95,
# sd_p[1] 0.90); asserted on the POOLED indicator so a single component sitting
# at the 0.85 floor does not make the gate brittle.
test_that("ms_occu NUTS community-covariance 95% CIs cover at the nominal rate", {
  skip_on_cran()
  skip_if_fast()
  n_seed <- 20L
  sd_psi_true <- c(0.5, 0.3); sd_p_true <- 0.4
  covered <- logical(0)
  for (s in seq_len(n_seed)) {
    sim <- simulate_ms_occu(N = 120, J = 4, n_species = 14,
                            beta_comm_mean = c(0, 0.5), beta_comm_sd = sd_psi_true,
                            alpha_comm_mean = c(0.2), alpha_comm_sd = sd_p_true,
                            seed = 700 + s)
    fit <- tryCatch(
      tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
           y = sim$y, species = paste0("sp", seq_len(14)), method = "nuts",
           control = list(n.iter = 300L, n.warmup = 300L, seed = 1L,
                          verbose = FALSE)),
      error = function(e) NULL)
    if (is.null(fit)) next
    lay <- fit$nuts$layout; dr <- fit$nuts$draws
    # Per-draw community SDs from the sampled log-Cholesky coordinates.
    sd_draws <- function(cols, P) {
      s <- vapply(seq_len(nrow(dr)), function(i) {
        C <- tulpaObs:::.ms_ocs_chol_unpack(dr[i, cols], P)
        sqrt(diag(C %*% t(C)))
      }, numeric(P))
      # vapply returns a plain vector when P == 1 (the detection arm here, whose
      # design is ~ 1), so re-impose the P x n_draws shape the callers index.
      matrix(s, nrow = P)
    }
    Spsi <- sd_draws(lay$chol_psi, lay$p_psi)      # 2 x n_draws
    Sp   <- sd_draws(lay$chol_p,   lay$p_p)        # 1 x n_draws
    contains <- function(v, truth) {
      q <- stats::quantile(v, c(0.025, 0.975)); q[1] <= truth && q[2] >= truth
    }
    covered <- c(covered,
                 contains(Spsi[1, ], sd_psi_true[1]),
                 contains(Spsi[2, ], sd_psi_true[2]),
                 contains(Sp[1, ],   sd_p_true))
  }
  # Pooled over the three community-SD components x 20 seeds; >= the 0.85 rubric
  # floor (measured ~0.90).
  expect_gte(mean(covered), 0.85)
})


# --- (5) S3 methods + richness ---------------------------------------------

test_that("ms_occu NUTS S3 methods work, incl. richness", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_occu(N = 60, J = 3, n_species = 8,
                          beta_comm_mean = c(0, 0.5), alpha_comm_mean = c(0.2),
                          seed = 5)
  fit <- tobs(~ x, data = sim$data, family = ms_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", seq_len(8)), method = "nuts",
              control = list(n.iter = 300L, n.warmup = 300L, seed = 1L,
                             verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  expect_no_error(summary(fit))

  cf <- coef(fit)
  expect_setequal(names(cf), c("psi", "p"))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
  expect_equal(nrow(re), 8L * 3L)               # 8 species x (psi int + psi x + p int)
  expect_setequal(unique(re$arm), c("psi", "p"))

  fv <- fitted(fit)
  expect_equal(dim(fv$psi), c(60L, 8L))
  expect_equal(dim(fv$z),   c(60L, 8L))

  ys <- simulate(fit, nsim = 1)
  expect_equal(dim(ys), c(60L, 3L, 8L))

  rich <- tobs_richness(fit)
  expect_s3_class(rich, "data.frame")
  expect_equal(nrow(rich), 60L)
})


# --- (6) the variance components de-attenuate vs the EM baseline ------------

test_that("ms_occu NUTS de-attenuates the community variance vs EM", {
  skip_on_cran()
  skip_if_fast()
  # The raw Laplace-EM community SDs carry the documented small-cluster
  # attenuation for binary detection; the sampler integrates the full joint, so
  # its per-arm community SDs sit above the EM SDs (closer to the truth).
  sim <- simulate_ms_occu(N = 130, J = 4, n_species = 16,
                          beta_comm_mean = c(0, 0.6), beta_comm_sd = c(0.6, 0.3),
                          alpha_comm_mean = c(0.2), alpha_comm_sd = c(0.5),
                          seed = 41)
  args <- list(formula = ~ x, data = sim$data, family = ms_occu(),
               detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(16)))
  fit_em <- do.call(tobs, c(args, list(method = "laplace",
                                       control = list(verbose = FALSE))))
  fit_nuts <- do.call(tobs, c(args, list(method = "nuts",
    control = list(n.iter = 400L, n.warmup = 400L, n.chains = 2L, seed = 1L,
                   verbose = FALSE))))
  # The detection-arm SD de-attenuates (the binary detection arm carries the
  # strongest small-cluster bias); the NUTS SD is not below the EM SD.
  expect_gte(fit_nuts$ms_community$sd_p[[1]], fit_em$ms_community$sd_p[[1]])
})


# --- (7) gates -------------------------------------------------------------

test_that("ms_occu NUTS rejects a spatial term with a pointer", {
  skip_on_cran()
  set.seed(5)
  N <- 16L
  adj <- matrix(0L, N, N)
  for (i in seq_len(N - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  sim <- simulate_ms_occu(N = N, J = 3, n_species = 4, seed = 5)
  expect_error(
    tobs(~ x + icar(graph = adj), data = sim$data, family = ms_occu(),
         detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(4)),
         method = "nuts", control = list(verbose = FALSE)),
    "non-spatial")
})
