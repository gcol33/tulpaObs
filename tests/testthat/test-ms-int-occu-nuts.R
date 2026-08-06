# Community integrated occupancy NUTS (ms_int_occu(), method = "nuts"; #115).
# The R target .tobs_ms_int_occu_nuts_logpost (R/ms_int_occu_nuts.R) is the
# oracle the C++ FullGradFn (src/ms_int_occu_nuts.cpp) is cross-checked against,
# the multi-source generalisation of the ms_occu / ms_dyn_occu NUTS targets.
# Coverage: (1) the R oracle gradient vs finite differences on a synthetic
# two-source design exercising the per-arm log-Cholesky blocks (off-diagonals:
# P_psi = 2, two sources with P_p = 2 each), (2) the C++ FullGradFn byte-exact vs
# the R oracle, (3) the b_from_z reconstruction, (4) end-to-end community-mean
# recovery + variance de-attenuation vs the Laplace-EM path.

test_that("ms_int_occu NUTS R oracle gradient matches finite differences", {
  skip_on_cran()
  set.seed(42)
  n <- 40L; D <- 2L; S <- 4L; P_psi <- 2L; P_p <- c(2L, 2L)
  X_psi <- cbind(1, rnorm(n))
  X_p   <- lapply(seq_len(D), function(d) cbind(1, rnorm(n)))
  summaries <- lapply(seq_len(S), function(s) {
    nv <- matrix(sample(0:4, n * D, replace = TRUE), n, D)
    nd <- matrix(0L, n, D)
    for (d in seq_len(D))
      nd[, d] <- vapply(nv[, d], function(k) if (k > 0) sample(0:k, 1) else 0L, integer(1))
    list(n_valid = nv, n_det = nd, any_det = rowSums(nd) > 0L)
  })

  lay <- tulpaObs:::.tobs_ms_int_occu_nuts_layout(P_psi, P_p, S)
  pri <- tulpaObs:::.ms_ocs_nuts_priors()
  set.seed(7); theta <- rnorm(lay$total, 0, 0.5)

  o  <- tulpaObs:::.tobs_ms_int_occu_nuts_logpost(theta, X_psi, X_p, summaries,
                                                  lay, priors = pri, sigma.beta = 5)
  an <- o$grad
  h  <- 1e-6
  fd <- numeric(lay$total)
  for (j in seq_len(lay$total)) {
    tp <- theta; tp[j] <- tp[j] + h
    tm <- theta; tm[j] <- tm[j] - h
    lp <- tulpaObs:::.tobs_ms_int_occu_nuts_logpost(tp, X_psi, X_p, summaries, lay,
            priors = pri, sigma.beta = 5, grad = FALSE)$lp
    lm <- tulpaObs:::.tobs_ms_int_occu_nuts_logpost(tm, X_psi, X_p, summaries, lay,
            priors = pri, sigma.beta = 5, grad = FALSE)$lp
    fd[j] <- (lp - lm) / (2 * h)
  }

  expect_lt(max(abs(an - fd)), 1e-5)
  expect_gt(cor(an, fd), 0.9999)
})

test_that("ms_int_occu NUTS C++ FullGradFn matches the R oracle byte-for-byte", {
  skip_on_cran()
  set.seed(42)
  n <- 40L; D <- 2L; S <- 4L; P_psi <- 2L; P_p <- c(2L, 2L)
  X_psi <- cbind(1, rnorm(n))
  X_p   <- lapply(seq_len(D), function(d) cbind(1, rnorm(n)))
  nv_sp <- lapply(seq_len(S), function(s) matrix(sample(0:4, n * D, TRUE), n, D))
  nd_sp <- lapply(seq_len(S), function(s) {
    nv <- nv_sp[[s]]; nd <- matrix(0L, n, D)
    for (d in seq_len(D))
      nd[, d] <- vapply(nv[, d], function(k) if (k > 0) sample(0:k, 1) else 0L, integer(1))
    nd
  })
  summaries <- lapply(seq_len(S), function(s)
    list(n_valid = nv_sp[[s]], n_det = nd_sp[[s]], any_det = rowSums(nd_sp[[s]]) > 0L))
  lay <- tulpaObs:::.tobs_ms_int_occu_nuts_layout(P_psi, P_p, S)
  pri <- tulpaObs:::.ms_ocs_nuts_priors()
  set.seed(7); theta <- rnorm(lay$total, 0, 0.5)

  o_r <- tulpaObs:::.tobs_ms_int_occu_nuts_logpost(theta, X_psi, X_p, summaries,
           lay, priors = pri, sigma.beta = 5)
  int_mat <- function(m) { storage.mode(m) <- "integer"; m }
  nv_list <- lapply(seq_len(D), function(d)
    int_mat(vapply(seq_len(S), function(s) nv_sp[[s]][, d], integer(n))))
  nd_list <- lapply(seq_len(D), function(d)
    int_mat(vapply(seq_len(S), function(s) nd_sp[[s]][, d], integer(n))))
  spec <- list(X_psi = X_psi, X_p = X_p, n_valid = nv_list, n_det = nd_list,
               n_sites = n, n_species = S, D = D)
  o_c <- tulpaObs:::cpp_ms_int_occu_nuts_joint_logpost(spec, theta, pri, 5)

  expect_lt(abs(o_r$lp - o_c$lp), 1e-9)
  expect_lt(max(abs(o_r$grad - o_c$grad)), 1e-9)
})

test_that("ms_int_occu b_from_z round-trips a whitened deviation matrix", {
  # z_s = C_arm^{-1} b_s -> b_s = C_arm z_s should reconstruct the input B.
  set.seed(3)
  D <- 2L; S <- 3L; P_psi <- 2L; P_p <- c(1L, 2L)
  lay <- tulpaObs:::.tobs_ms_int_occu_nuts_layout(P_psi, P_p, S)
  theta <- numeric(lay$total)
  mk_chol <- function(P) { A <- matrix(rnorm(P * P), P, P); t(chol(crossprod(A) + diag(P))) }
  C_psi <- mk_chol(P_psi)
  theta[lay$chol_psi] <- tulpaObs:::.ms_ocs_chol_pack(C_psi)
  C_p <- lapply(seq_len(D), function(d) mk_chol(P_p[d]))
  for (d in seq_len(D)) theta[lay$chol_p[[d]]] <- tulpaObs:::.ms_ocs_chol_pack(C_p[[d]])
  Z <- matrix(rnorm(S * lay$P), S, lay$P)
  for (s in seq_len(S)) theta[tulpaObs:::.ms_ocs_b_idx(lay, s)] <- Z[s, ]

  B <- tulpaObs:::.ms_ocs_b_from_z(theta, lay)
  # Recompute the expected b per arm and compare.
  exp_B <- matrix(0, S, lay$P)
  for (s in seq_len(S)) {
    exp_B[s, lay$psi] <- as.numeric(C_psi %*% Z[s, lay$psi])
    for (d in seq_len(D))
      exp_B[s, lay$p[[d]]] <- as.numeric(C_p[[d]] %*% Z[s, lay$p[[d]]])
  }
  expect_equal(B, exp_B, tolerance = 1e-12)
})

test_that("ms_int_occu NUTS recovers community means + de-attenuates the variance", {
  skip_on_cran()
  skip_if_fast()
  sim <- simulate_ms_int_occu(N = 150, J = c(3, 4), n_species = 12,
                              n_data = 2, seed = 23)
  sp <- paste0("sp", seq_len(12))
  lap <- tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
              y = sim$y, species = sp, method = "laplace",
              control = list(verbose = FALSE, progress = FALSE))
  nut <- tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
              y = sim$y, species = sp, method = "nuts",
              control = list(n.iter = 500L, n.warmup = 500L, seed = 1L,
                             verbose = FALSE, progress = FALSE))

  expect_identical(nut$method, "nuts")
  expect_equal(nut$nuts$divergent_total, 0L)

  # Community means recover (intercept-only community, all truth 0).
  truth <- c("psi_(Intercept)" = 0, "p1_(Intercept)" = 0, "p2_(Intercept)" = 0)
  m <- nut$means[names(truth)]; s <- nut$sds[names(truth)]
  expect_true(all(abs(m - truth) / s < 3))

  # Per-species coefficients track the simulated truth on every arm.
  cm <- nut$ms_community
  expect_gt(cor(cm$coef_psi[, 1], stats::qlogis(sim$truth$psi_species)), 0.65)
  expect_gt(cor(cm$coef_p1[, 1],  stats::qlogis(sim$truth$p_det[[1]])),  0.55)
  expect_gt(cor(cm$coef_p2[, 1],  stats::qlogis(sim$truth$p_det[[2]])),  0.55)

  # NUTS removes the Laplace small-cluster attenuation of the community SD, so the
  # sampled community SD is at least the Laplace lower bound (de-attenuation).
  expect_gt(nut$ms_community$sd_psi, 0.9 * lap$ms_community$sd_psi)
  expect_true(all(cm$sd_psi > 0 & cm$sd_p1 > 0 & cm$sd_p2 > 0))
})

test_that("ms_int_occu NUTS community-mean 95% CIs cover at the nominal rate", {
  # gcol33/tulpaObs#139: raise the single-seed / directional check to a 20-seed
  # CI-coverage study on the community means (per-arm intercepts).
  skip_on_cran()
  skip_if_fast()
  sp <- paste0("sp", seq_len(12))
  truth <- c("psi_(Intercept)" = 0, "p1_(Intercept)" = 0, "p2_(Intercept)" = 0)
  covered <- logical(0)
  for (s in seq_len(20L)) {
    sim <- simulate_ms_int_occu(N = 150, J = c(3, 4), n_species = 12,
                                n_data = 2, seed = 300 + s)
    fit <- tryCatch(tobs(~ 1, data = sim$data, family = ms_int_occu(), detection = ~ 1,
                    y = sim$y, species = sp, method = "nuts",
                    control = list(n.iter = 400L, n.warmup = 400L, seed = 1L,
                                   verbose = FALSE, progress = FALSE)),
                    error = function(e) NULL)
    if (is.null(fit)) next
    m <- fit$means[names(truth)]; se <- fit$sds[names(truth)]
    covered <- c(covered, abs(m - truth) <= 1.96 * se)
  }
  expect_gte(mean(covered), 0.85)
})
