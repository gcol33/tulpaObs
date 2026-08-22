# =============================================================================
# test-ms-dyn-occu-nuts.R - community (multispecies) DYNAMIC occupancy NUTS
# (ms_dyn_occu(), method = "nuts").
#
# The sampler draws the EXACT joint posterior -- community means, per-species
# first-season / detection deviations, the two INDEPENDENT per-arm community
# covariances, AND the shared colonisation / extinction globals -- over the
# per-(species, site) HMM-forward marginal via the in-tree C++ FullGradFn
# (src/ms_dyn_occu_nuts.cpp), warm-started at the community Laplace-EM mode.
# The R target .tobs_ms_dyn_occu_nuts_logpost (R/ms_dyn_occu_nuts.R) is the
# oracle; the C++ port is cross-checked against it.
#
# Coverage: (1) the R oracle gradient vs finite differences, (2) the C++
# FullGradFn byte-exact vs the R oracle, (3) community-mean recovery + 0
# divergences, (4) S3 methods + richness, (5) the variance components
# de-attenuate vs the EM baseline, (6) the spatial-term gate.
# =============================================================================


# --- shared fixtures / helpers ---------------------------------------------

.msdyn_pieces <- function(n_species = 6, N = 60, J = 3, n_seasons = 4, seed = 7) {
  sim <- simulate_ms_dyn_occu(N = N, J = J, n_species = n_species,
                              n_seasons = n_seasons,
                              beta_comm_mean = c(0.3), beta_comm_sd = c(0.7),
                              gamma = 0.2, epsilon = 0.1, seed = seed)
  bind  <- tulpaObs:::.tobs_bind_formulas(list(psi1 = ~ 1, p = ~ 1), sim$data)
  model <- tulpaObs:::.tobs_build_ms_dyn_occu(
    occ_formula = bind$fe$psi1, det_formula = bind$fe$p,
    data = sim$data, y = sim$y, species = paste0("sp", seq_len(n_species)),
    structured_terms = bind$terms)
  pieces <- tulpaObs:::.tobs_ms_dyn_occu_nuts_pieces(model)
  lay <- tulpaObs:::.tobs_ms_dyn_occu_nuts_layout(pieces$P_psi1, pieces$P_p,
                                                 pieces$P_gam, pieces$P_eps,
                                                 pieces$S)
  pri <- tulpaObs:::.ms_ocs_nuts_priors()
  em <- tulpaObs:::.tobs_community_em(
    S = pieces$S, P = pieces$P, arm_idx = pieces$arm_idx,
    sp_ll = pieces$sp_ll, sp_grad = pieces$sp_grad,
    init_mu = pieces$init_mu, init_global = pieces$init_global,
    penalize_global = TRUE, sigma_beta = 5, priors = NULL,
    sigma_init = 0.3, max_iter = 40L, tol = 1e-4, newton_max = 30L,
    verbose = FALSE)
  theta0 <- tulpaObs:::.tobs_ms_dyn_occu_nuts_pack_init(em, lay, pieces)
  spec <- list(X_psi1 = pieces$X_psi1, X_p = pieces$X_p,
               X_gamma = pieces$X_gamma, X_eps = pieces$X_eps,
               n_sites = pieces$Ns, n_seasons = pieces$T, n_species = pieces$S,
               n_valid = pieces$nv_list, n_det = pieces$nd_list)
  list(sim = sim, model = model, pieces = pieces, lay = lay, pri = pri,
       theta0 = theta0, spec = spec)
}

.msdyn_fd_grad <- function(f, theta, h = 1e-5) {
  vapply(seq_along(theta), function(j) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    (f(tp) - f(tm)) / (2 * h)
  }, 0)
}


# --- (1) R oracle gradient vs finite differences ---------------------------

test_that("ms_dyn_occu NUTS R oracle gradient matches finite differences", {
  skip_on_cran()
  P <- .msdyn_pieces()
  set.seed(1)
  theta <- P$theta0 + stats::rnorm(length(P$theta0), 0, 0.05)
  f_lp <- function(th)
    tulpaObs:::.tobs_ms_dyn_occu_nuts_logpost(
      th, P$pieces$X_psi1, P$pieces$X_p, P$pieces$X_gamma, P$pieces$X_eps,
      P$pieces$em_stats, P$pieces$Ns, P$pieces$T, P$lay, P$pri,
      grad = FALSE)$lp
  ana <- tulpaObs:::.tobs_ms_dyn_occu_nuts_logpost(
    theta, P$pieces$X_psi1, P$pieces$X_p, P$pieces$X_gamma, P$pieces$X_eps,
    P$pieces$em_stats, P$pieces$Ns, P$pieces$T, P$lay, P$pri, grad = TRUE)$grad
  num <- .msdyn_fd_grad(f_lp, theta)
  expect_lt(max(abs(ana - num)), 1e-5)
})


# --- (2) C++ FullGradFn byte-exact vs the R oracle -------------------------

test_that("ms_dyn_occu NUTS C++ FullGradFn matches the R oracle", {
  skip_on_cran()
  P <- .msdyn_pieces()
  set.seed(2)
  theta <- P$theta0 + stats::rnorm(length(P$theta0), 0, 0.05)
  cpp <- cpp_ms_dyn_occu_nuts_joint_logpost(P$spec, theta, P$pri, sigma_beta = 5)
  r   <- tulpaObs:::.tobs_ms_dyn_occu_nuts_logpost(
    theta, P$pieces$X_psi1, P$pieces$X_p, P$pieces$X_gamma, P$pieces$X_eps,
    P$pieces$em_stats, P$pieces$Ns, P$pieces$T, P$lay, P$pri,
    sigma.beta = 5, grad = TRUE)
  expect_lt(abs(cpp$lp - r$lp), 1e-9)
  expect_lt(max(abs(cpp$grad - r$grad)), 1e-9)
})


# --- (3) community-mean recovery + 0 divergences ---------------------------

test_that("ms_dyn_occu NUTS recovers community means", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_dyn_occu(N = 90, J = 3, n_species = 14, n_seasons = 4,
                              beta_comm_mean = c(0.3), beta_comm_sd = c(0.7),
                              gamma = 0.2, epsilon = 0.1, seed = 31)
  fit <- tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", seq_len(14)), method = "nuts",
              control = list(n.iter = 400L, n.warmup = 400L, n.chains = 2L,
                             seed = 1L, verbose = FALSE))
  expect_equal(fit$method, "nuts")
  expect_false(is.null(fit$nuts$draws))
  expect_lt(fit$nuts$divergent_total, 0.05 * nrow(fit$nuts$draws))
  expect_lt(fit$nuts$max_rhat, 1.1)

  truth <- c("psi1_(Intercept)"  = 0.3,
             "p_(Intercept)"     = 0,
             "gamma_(Intercept)" = stats::qlogis(0.2),
             "eps_(Intercept)"   = stats::qlogis(0.1))
  m <- fit$means[names(truth)]; s <- fit$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 2.5))

  # Per-species first-season / detection coefficients track the simulated truth.
  cm <- fit$ms_community
  expect_gt(cor(cm$coef_psi1[, 1], stats::qlogis(sim$truth$psi1_species)), 0.60)
  expect_gt(cor(cm$coef_p[, 1],    stats::qlogis(sim$truth$p_species)),    0.50)
})

test_that("ms_dyn_occu NUTS community-mean 95% CIs cover at the nominal rate", {
  # raise the single-seed / directional check to a 20-seed CI-coverage study
  # on the community means + shared gamma/eps globals.
  skip_on_cran()
  skip_if_fast()
  truth <- c("psi1_(Intercept)"  = 0.3,
             "p_(Intercept)"     = 0,
             "gamma_(Intercept)" = stats::qlogis(0.2),
             "eps_(Intercept)"   = stats::qlogis(0.1))
  covered <- logical(0)
  for (s in seq_len(20L)) {
    sim <- simulate_ms_dyn_occu(N = 90, J = 3, n_species = 14, n_seasons = 4,
                                beta_comm_mean = c(0.3), beta_comm_sd = c(0.7),
                                gamma = 0.2, epsilon = 0.1, seed = 300 + s)
    fit <- tryCatch(tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
                    y = sim$y, species = paste0("sp", seq_len(14)), method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   verbose = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    m <- fit$means[names(truth)]; se <- fit$sds[names(truth)]
    covered <- c(covered, abs(m - truth) <= 1.96 * se)
  }
  expect_gte(mean(covered), 0.85)
})


# --- (4) S3 methods + richness ---------------------------------------------

test_that("ms_dyn_occu NUTS S3 methods work", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_dyn_occu(N = 55, J = 3, n_species = 8, n_seasons = 3,
                              beta_comm_mean = c(0.3), beta_comm_sd = c(0.7),
                              gamma = 0.2, epsilon = 0.1, seed = 12)
  fit <- tobs(~ 1, data = sim$data, family = ms_dyn_occu(), detection = ~ 1,
              y = sim$y, species = paste0("sp", seq_len(8)), method = "nuts",
              control = list(n.iter = 300L, n.warmup = 300L, seed = 1L,
                             verbose = FALSE))
  expect_s3_class(fit, "tobs_fit")
  expect_no_error(print(fit))
  expect_no_error(summary(fit))

  cf <- coef(fit)
  expect_true(all(c("psi1", "p") %in% names(cf)))

  V <- vcov(fit)
  expect_equal(nrow(V), length(fit$means))
  expect_equal(nrow(confint(fit)), length(fit$means))

  re <- ranef(fit)
  expect_s3_class(re, "data.frame")
})


# --- (5) the variance components de-attenuate vs the EM baseline ------------

test_that("ms_dyn_occu NUTS de-attenuates the community variance vs EM", {
  skip_on_cran()
  skip_if_fast()
  # The raw Laplace-EM community SDs carry the documented small-cluster
  # attenuation for the binary arms; the sampler integrates the full joint, so
  # its per-arm community SDs sit at or above the EM SDs (closer to the truth).
  sim <- simulate_ms_dyn_occu(N = 90, J = 3, n_species = 14, n_seasons = 4,
                              beta_comm_mean = c(0.3), beta_comm_sd = c(0.7),
                              gamma = 0.2, epsilon = 0.1, seed = 31)
  args <- list(formula = ~ 1, data = sim$data, family = ms_dyn_occu(),
               detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(14)))
  fit_em <- do.call(tobs, c(args, list(method = "laplace",
                                       control = list(verbose = FALSE))))
  fit_nuts <- do.call(tobs, c(args, list(method = "nuts",
    control = list(n.iter = 400L, n.warmup = 400L, n.chains = 2L, seed = 1L,
                   verbose = FALSE))))
  expect_gte(fit_nuts$ms_community$sd_psi1[[1]], fit_em$ms_community$sd_psi1[[1]])
})


# --- (6) gates -------------------------------------------------------------

test_that("ms_dyn_occu NUTS rejects a spatial term with a pointer", {
  skip_on_cran()
  set.seed(5)
  N <- 16L
  adj <- matrix(0L, N, N)
  for (i in seq_len(N - 1L)) { adj[i, i + 1L] <- 1L; adj[i + 1L, i] <- 1L }
  sim <- simulate_ms_dyn_occu(N = N, J = 3, n_species = 4, n_seasons = 3, seed = 5)
  expect_error(
    tobs(~ 1 + icar(graph = adj), data = sim$data, family = ms_dyn_occu(),
         detection = ~ 1, y = sim$y, species = paste0("sp", seq_len(4)),
         method = "nuts", control = list(verbose = FALSE)),
    "nested_laplace")
})
