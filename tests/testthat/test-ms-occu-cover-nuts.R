# Community joint occupancy + cover NUTS (ms_occu_cover(), method = "nuts";
# #115 part B7). The R target .tobs_ms_occu_cover_nuts_logpost
# (R/ms_occu_cover_nuts.R) is the oracle the C++ FullGradFn
# (src/ms_occu_cover_nuts.cpp) is cross-checked against, the joint-cover analogue
# of the ms_occu / ms_int_occu community NUTS targets: three non-centered
# per-species arms (occ + p + pos) each with a log-Cholesky community covariance,
# plus ONE shared community log-dispersion scalar. It reuses the per-species
# occu_cover marginal + gradient (.occu_cover_sp_ll / _grad) from the laplace path
# verbatim, and the shared C++ cover-density kernels (occu_coupling_shared.h).
# Coverage: (1) the R oracle gradient vs finite differences (incl the shared
# log-dispersion coord), (2) the C++ FullGradFn byte-exact vs the oracle, (3) the
# b_from_z reconstruction, (4) end-to-end community-mean recovery + variance
# de-attenuation vs the Laplace-EM path.

.msoc_nuts_visits <- function(N, J, visit_data) {
  long <- data.frame(site_id = rep(seq_len(N), each = J),
                     visit   = rep(seq_len(J), times = N), yy = 0L,
                     det_cov1 = visit_data$det_cov1, pos_cov1 = visit_data$pos_cov1)
  od <- tobs_data(long, y = "yy", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  od$det.covs
}

test_that("ms_occu_cover NUTS R oracle gradient matches finite differences", {
  skip_on_cran()
  skip_if_fast()
  set.seed(1)
  sim <- simulate_ms_occu_cover(n_species = 4, N = 40, J = 3,
                                positive = "lognormal", seed = 1)
  vis <- .msoc_nuts_visits(40, 3, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1, y = sim$y,
              y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE, max.iter = 5L))
  model <- fit$model
  S <- model$n_species
  pil <- model$process_info
  P_occ <- pil[[1]]$p; P_p <- pil[[2]]$p; P_pos <- pil[[3]]$p
  views <- lapply(seq_len(S), function(s)
    tulpaObs:::.ms_occu_cover_species_view(model, s))
  lay <- tulpaObs:::.tobs_ms_occu_cover_nuts_layout(P_occ, P_p, P_pos, S)
  pri <- tulpaObs:::.ms_ocs_nuts_priors()
  set.seed(7); theta <- rnorm(lay$total, 0, 0.4); theta[lay$log_disp] <- log(0.6)

  o  <- tulpaObs:::.tobs_ms_occu_cover_nuts_logpost(theta, views, lay,
          priors = pri, sigma.beta = 5)
  an <- o$grad; h <- 1e-6; fd <- numeric(lay$total)
  for (j in seq_len(lay$total)) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    lp <- tulpaObs:::.tobs_ms_occu_cover_nuts_logpost(tp, views, lay,
            priors = pri, sigma.beta = 5, grad = FALSE)$lp
    lm <- tulpaObs:::.tobs_ms_occu_cover_nuts_logpost(tm, views, lay,
            priors = pri, sigma.beta = 5, grad = FALSE)$lp
    fd[j] <- (lp - lm) / (2 * h)
  }
  # Relative match: the log-dispersion coordinate can carry a large magnitude far
  # from the mode, so compare on the scale of the gradient.
  expect_lt(max(abs(an - fd)) / max(1, max(abs(fd))), 1e-4)
  expect_gt(cor(an, fd), 0.9999)
})

test_that("ms_occu_cover NUTS C++ FullGradFn matches the R oracle byte-for-byte", {
  skip_on_cran()
  skip_if_fast()
  set.seed(1)
  sim <- simulate_ms_occu_cover(n_species = 4, N = 40, J = 3,
                                positive = "lognormal", seed = 1)
  vis <- .msoc_nuts_visits(40, 3, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1, y = sim$y,
              y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE, max.iter = 5L))
  model <- fit$model
  S <- model$n_species; pil <- model$process_info
  P_occ <- pil[[1]]$p; P_p <- pil[[2]]$p; P_pos <- pil[[3]]$p
  views <- lapply(seq_len(S), function(s)
    tulpaObs:::.ms_occu_cover_species_view(model, s))
  lay <- tulpaObs:::.tobs_ms_occu_cover_nuts_layout(P_occ, P_p, P_pos, S)
  pri <- tulpaObs:::.ms_ocs_nuts_priors()
  set.seed(7); theta <- rnorm(lay$total, 0, 0.4); theta[lay$log_disp] <- log(0.6)

  o_r  <- tulpaObs:::.tobs_ms_occu_cover_nuts_logpost(theta, views, lay,
            priors = pri, sigma.beta = 5)
  spec <- tulpaObs:::.tobs_ms_occu_cover_nuts_spec(model)
  o_c  <- tulpaObs:::cpp_ms_occu_cover_nuts_joint_logpost(spec, theta, pri, 5)

  expect_lt(abs(o_r$lp - o_c$lp), 1e-8)
  expect_lt(max(abs(o_r$grad - o_c$grad)), 1e-8)
})

test_that("ms_occu_cover b_from_z round-trips a whitened deviation matrix", {
  set.seed(3)
  S <- 3L; P_occ <- 2L; P_p <- 1L; P_pos <- 2L
  lay <- tulpaObs:::.tobs_ms_occu_cover_nuts_layout(P_occ, P_p, P_pos, S)
  theta <- numeric(lay$total)
  mk_chol <- function(P) { A <- matrix(rnorm(P * P), P, P); t(chol(crossprod(A) + diag(P))) }
  C_occ <- mk_chol(P_occ); C_p <- mk_chol(P_p); C_pos <- mk_chol(P_pos)
  theta[lay$chol_occ] <- tulpaObs:::.ms_ocs_chol_pack(C_occ)
  theta[lay$chol_p]   <- tulpaObs:::.ms_ocs_chol_pack(C_p)
  theta[lay$chol_pos] <- tulpaObs:::.ms_ocs_chol_pack(C_pos)
  Z <- matrix(rnorm(S * lay$P), S, lay$P)
  for (s in seq_len(S)) theta[tulpaObs:::.ms_ocs_b_idx(lay, s)] <- Z[s, ]

  B <- tulpaObs:::.tobs_ms_occu_cover_nuts_b_from_z(theta, lay)
  exp_B <- matrix(0, S, lay$P)
  for (s in seq_len(S)) {
    exp_B[s, lay$occ] <- as.numeric(C_occ %*% Z[s, lay$occ])
    exp_B[s, lay$p]   <- as.numeric(C_p   %*% Z[s, lay$p])
    exp_B[s, lay$pos] <- as.numeric(C_pos %*% Z[s, lay$pos])
  }
  expect_equal(B, exp_B, tolerance = 1e-12)
})

test_that("ms_occu_cover NUTS recovers community means + de-attenuates the variance", {
  skip_on_cran()
  skip_if_fast()
  set.seed(21)
  sim <- simulate_ms_occu_cover(
    n_species = 12, N = 90, J = 4,
    mu_occ = c(stats::qlogis(0.45), 0.7), mu_p = c(0.2, -0.4),
    mu_pos = c(log(0.12), 0.5), sd_occ = 0.5, sd_p = 0.4, sd_pos = 0.4,
    positive = "lognormal", sigma_pos = 0.4, seed = 21)
  vis <- .msoc_nuts_visits(90, 4, sim$visit_data)
  args0 <- list(occurrence = ~ occ_cov1, data = sim$data,
                family = ms_occu_cover("lognormal"), detection = ~ det_cov1,
                positive = ~ pos_cov1, y = sim$y, y_pos = sim$y_pos,
                visits = vis, species = sim$species)
  lap <- do.call(tobs, c(args0, method = "laplace", list(control = list(verbose = FALSE))))
  nut <- do.call(tobs, c(args0, method = "nuts",
    list(control = list(n.iter = 500L, n.warmup = 500L, seed = 1L, verbose = FALSE))))

  expect_identical(nut$method, "nuts")
  expect_equal(nut$nuts$divergent_total, 0L)

  # Community means recover (within ~3 SE of truth on every beta coordinate).
  beta_means <- nut$means[seq_len(6)]
  truth <- c(sim$truth$mu_occ, sim$truth$mu_p, sim$truth$mu_pos)
  expect_true(all(abs(beta_means - truth) / nut$sds[seq_len(6)] < 3))
  # Shared dispersion recovered.
  expect_lt(abs(exp(nut$means[["log_sigma_pos"]]) - 0.4), 0.12)

  # Per-species coefficients track the simulated truth.
  cm <- nut$ms_community
  expect_gt(min(diag(cor(cm$coef_occ, sim$truth$beta_occ))), 0.65)
  expect_gt(min(diag(cor(cm$coef_pos, sim$truth$beta_pos))), 0.75)

  # NUTS removes the Laplace small-cluster attenuation of the community SD: the
  # sampled occupancy-arm SD recovers toward truth (0.5) and clears the Laplace
  # estimate, which collapses badly at this per-species n.
  expect_gt(mean(nut$ms_community$sd_occ), mean(lap$ms_community$sd_occ))
  expect_gt(mean(nut$ms_community$sd_occ), 0.25)
})

test_that("ms_occu_cover dispersion-RE oracle gradient matches finite differences", {
  # The per-species dispersion-RE variant (#115 B7 follow-up): a fourth 1-D
  # community arm log_disp_s = mu_ld + sigma_ld * z_ld_s. Validates the analytic
  # gradient (incl the mu_ld mean and the sigma_ld = chol_ld coordinate) on a real
  # small community model. The compiled sampler is the follow-up (C++ port).
  skip_on_cran()
  skip_if_fast()
  set.seed(1)
  sim <- simulate_ms_occu_cover(n_species = 4, N = 40, J = 3,
                                positive = "lognormal", seed = 1)
  vis <- .msoc_nuts_visits(40, 3, sim$visit_data)
  fit <- tobs(~ occ_cov1, data = sim$data, family = ms_occu_cover("lognormal"),
              detection = ~ det_cov1, positive = ~ pos_cov1, y = sim$y,
              y_pos = sim$y_pos, visits = vis, species = sim$species,
              method = "laplace", control = list(verbose = FALSE, max.iter = 5L))
  model <- fit$model; S <- model$n_species; pil <- model$process_info
  P_occ <- pil[[1]]$p; P_p <- pil[[2]]$p; P_pos <- pil[[3]]$p
  views <- lapply(seq_len(S), function(s)
    tulpaObs:::.ms_occu_cover_species_view(model, s))
  lay <- tulpaObs:::.tobs_ms_occu_cover_re_disp_layout(P_occ, P_p, P_pos, S)
  pri <- tulpaObs:::.ms_ocs_nuts_priors()
  set.seed(7); theta <- rnorm(lay$total, 0, 0.4); theta[lay$mu_ld] <- log(0.6)

  o  <- tulpaObs:::.tobs_ms_occu_cover_re_disp_logpost(theta, views, lay,
          priors = pri, sigma.beta = 5)
  an <- o$grad; h <- 1e-6; fd <- numeric(lay$total)
  for (j in seq_len(lay$total)) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    lp <- tulpaObs:::.tobs_ms_occu_cover_re_disp_logpost(tp, views, lay,
            priors = pri, sigma.beta = 5, grad = FALSE)$lp
    lm <- tulpaObs:::.tobs_ms_occu_cover_re_disp_logpost(tm, views, lay,
            priors = pri, sigma.beta = 5, grad = FALSE)$lp
    fd[j] <- (lp - lm) / (2 * h)
  }
  expect_lt(max(abs(an - fd)) / max(1, max(abs(fd))), 1e-4)
  expect_gt(cor(an, fd), 0.9999)
})
