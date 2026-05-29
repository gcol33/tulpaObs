# =============================================================================
# probe_bench_v3_vs_jc.R - measure how the v3 / jc wall-clock gap scales
# with N. Issue #32 claims "10-100x on real EVA scale (~1000 cells)";
# probe_recovery measured 161x at N=49. This script extends to N=100, 196.
# Two seeds per N to expose variance.
# =============================================================================

suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})

make_grid_adj <- function(nrow, ncol) {
  N <- nrow * ncol
  adj <- matrix(0L, N, N)
  idx <- function(i, j) (i - 1L) * ncol + j
  for (i in seq_len(nrow)) {
    for (j in seq_len(ncol)) {
      k <- idx(i, j)
      if (i > 1L)    adj[k, idx(i - 1L, j)] <- 1L
      if (i < nrow)  adj[k, idx(i + 1L, j)] <- 1L
      if (j > 1L)    adj[k, idx(i, j - 1L)] <- 1L
      if (j < ncol)  adj[k, idx(i, j + 1L)] <- 1L
    }
  }
  adj
}

fit_one <- function(N_side, J, engine, seed) {
  N <- N_side * N_side
  adj <- make_grid_adj(N_side, N_side)
  sim <- simulate_occu_cover(
    N = N, J = J, n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive = "lognormal", sigma_pos = 0.35,
    adj = adj, sigma = 0.8, alpha = 1.0, seed = seed
  )
  visit_df <- data.frame(
    det_cov1 = rep(sim$visit_data$det_cov1, length.out = N * J),
    pos_cov1 = rep(sim$visit_data$pos_cov1, length.out = N * J)
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    tobs(formula   = ~ occ_cov1 + icar(graph = adj),
         detection = ~ det_cov1,
         data      = sim$data,
         family    = occu_cover(positive = "lognormal"),
         method    = "nested_laplace",
         y         = sim$y,
         y_pos     = sim$y_pos,
         visits    = visit_df,
         positive  = ~ pos_cov1,
         control   = list(engine = engine, verbose = FALSE)),
    error = function(e) list(error = conditionMessage(e))
  )
  dt <- proc.time()[["elapsed"]] - t0
  ok <- is.null(fit$error)
  list(dt = dt, ok = ok, msg = if (!ok) fit$error else NA_character_)
}

main <- function() {
  cells <- c(49L, 100L, 196L)  # 7x7, 10x10, 14x14
  J <- 4L
  seeds <- c(11L, 12L)
  cat("== bench v3_nested vs joint_coupled (J=4) ==\n\n")
  cat(sprintf("  %5s   %6s   %12s   %12s   %10s\n",
              "N", "seed", "v3 (s)", "jc (s)", "speedup"))
  for (N in cells) {
    Ns <- as.integer(sqrt(N))
    for (s in seeds) {
      r_v3 <- fit_one(Ns, J, "v3_nested",      s)
      r_jc <- fit_one(Ns, J, "joint_coupled",  s)
      if (!r_v3$ok) { cat(sprintf("  %5d   %6d   v3 ERROR: %s\n", N, s, r_v3$msg)); next }
      if (!r_jc$ok) { cat(sprintf("  %5d   %6d   jc ERROR: %s\n", N, s, r_jc$msg)); next }
      cat(sprintf("  %5d   %6d   %12.2f   %12.2f   %9.0fx\n",
                  N, s, r_v3$dt, r_jc$dt, r_v3$dt / max(r_jc$dt, 1e-3)))
    }
  }
}

invisible(main())
