# Community integrated occupancy NUTS (ms_int_occu(), method = "nuts"; #115).
# The R target .tobs_ms_int_occu_nuts_logpost (R/ms_int_occu_nuts.R) is the
# oracle a C++ FullGradFn port will be cross-checked against, the multi-source
# generalisation of the ms_occu / ms_dyn_occu NUTS targets. This block validates
# the oracle's analytic gradient against finite differences on a synthetic
# multi-source design that exercises the per-arm log-Cholesky blocks (including
# off-diagonals: P_psi = 2, two sources with P_p = 2 each). The compiled sampler
# + recovery are a follow-up (the C++ port).

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

  B <- tulpaObs:::.tobs_ms_int_occu_nuts_b_from_z(theta, lay)
  # Recompute the expected b per arm and compare.
  exp_B <- matrix(0, S, lay$P)
  for (s in seq_len(S)) {
    exp_B[s, lay$psi] <- as.numeric(C_psi %*% Z[s, lay$psi])
    for (d in seq_len(D))
      exp_B[s, lay$p[[d]]] <- as.numeric(C_p[[d]] %*% Z[s, lay$p[[d]]])
  }
  expect_equal(B, exp_B, tolerance = 1e-12)
})
