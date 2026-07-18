# Community joint occupancy + cover NUTS (ms_occu_cover(), method = "nuts";
# #115 part B7). The R target .tobs_ms_occu_cover_nuts_logpost
# (R/ms_occu_cover_nuts.R) is the oracle a C++ FullGradFn port
# (src/ms_occu_cover_nuts.cpp) will be cross-checked against, the joint-cover
# analogue of the ms_occu / ms_int_occu community NUTS targets: three
# non-centered per-species arms (occ + p + pos) each with a log-Cholesky
# community covariance, plus ONE shared community log-dispersion scalar. It reuses
# the per-species occu_cover marginal + gradient (.occu_cover_sp_ll / _grad) from
# the laplace path verbatim. This block validates the oracle's analytic gradient
# (including the shared log-dispersion coordinate) against finite differences on a
# real small community model, plus the b_from_z reconstruction. The compiled
# sampler + recovery are the follow-up (the C++ port).

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
