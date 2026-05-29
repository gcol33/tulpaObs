# =============================================================================
# probe_recovery_v3_vs_jc.R - compare v3_nested vs joint_coupled by their
# recovery of truth across a 10-seed sweep, rather than against each other.
#
# Each engine runs in its native default mode (v3 jointly estimates log_sigma
# and log_sigma_pos via outer BFGS; jc pins sigma_pos at the empirical SD and
# integrates sigma + alpha on the outer grid). Per-engine, per-seed:
#   * |intercept estimate - truth|
#   * |alpha estimate - truth|
#   * field correlation with truth
#   * wall-clock
# Reports mean and max bias across seeds for each engine. If jc bias is
# comparable or better, the default flip is justified.
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

ALPHA_TRUE     <- 1.0
SIGMA_TRUE     <- 0.8
SIGMA_POS_TRUE <- 0.35

sim_one <- function(seed, N_side = 7L, J = 4L) {
  N <- N_side * N_side
  adj <- make_grid_adj(N_side, N_side)
  sim <- simulate_occu_cover(
    N         = N,
    J         = J,
    n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive  = "lognormal",
    sigma_pos = SIGMA_POS_TRUE,
    adj       = adj,
    sigma     = SIGMA_TRUE,
    alpha     = ALPHA_TRUE,
    seed      = seed
  )
  list(sim = sim, adj = adj)
}

fit_one <- function(sim_pack, engine) {
  sim <- sim_pack$sim
  adj <- sim_pack$adj
  N <- nrow(sim$y); J <- ncol(sim$y)
  visit_df <- data.frame(
    det_cov1 = rep(sim$visit_data$det_cov1, length.out = N * J),
    pos_cov1 = rep(sim$visit_data$pos_cov1, length.out = N * J)
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    tobs(
      formula   = ~ occ_cov1 + icar(graph = adj),
      detection = ~ det_cov1,
      data      = sim$data,
      family    = occu_cover(positive = "lognormal"),
      method    = "nested_laplace",
      y         = sim$y,
      y_pos     = sim$y_pos,
      visits    = visit_df,
      positive  = ~ pos_cov1,
      control   = list(engine = engine, verbose = FALSE)
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  dt <- proc.time()[["elapsed"]] - t0
  list(fit = fit, dt = dt)
}

# Extract a common (alpha, field, intercepts) tuple from either engine.
extract <- function(fit) {
  if (!is.null(fit$error)) return(NULL)
  m <- fit$means
  alpha_hat <- if ("alpha" %in% names(m)) unname(m[["alpha"]]) else NA_real_
  field_hat <- fit$spatial_field
  list(
    psi_int   = unname(m[["psi_(Intercept)"]]),
    p_int     = unname(m[["p_(Intercept)"]]),
    pos_int   = unname(m[["pos_(Intercept)"]]),
    alpha     = alpha_hat,
    field     = field_hat
  )
}

main <- function(n_seeds = 10L, N_side = 7L, J = 4L) {
  cat("== probe recovery v3_nested vs joint_coupled ==\n")
  cat(sprintf("Fixture: %d x %d grid (N=%d), J=%d\n",
              N_side, N_side, N_side * N_side, J))
  cat(sprintf("Truth: alpha=%g, sigma=%g, sigma_pos=%g\n",
              ALPHA_TRUE, SIGMA_TRUE, SIGMA_POS_TRUE))
  cat(sprintf("Seeds: 1..%d\n\n", n_seeds))

  rows <- list()
  for (s in seq_len(n_seeds)) {
    cat(sprintf("seed %d ... ", s))
    sp <- sim_one(s, N_side = N_side, J = J)
    psi_truth <- sp$sim$truth$beta_occ[1L]
    p_truth   <- sp$sim$truth$beta_p[1L]
    pos_truth <- sp$sim$truth$beta_pos[1L]
    f_truth   <- sp$sim$truth$f

    r_v3 <- fit_one(sp, "v3_nested")
    r_jc <- fit_one(sp, "joint_coupled")
    ex_v3 <- extract(r_v3$fit)
    ex_jc <- extract(r_jc$fit)
    if (is.null(ex_v3) || is.null(ex_jc)) {
      cat(" SKIP (error)\n"); next
    }
    field_cor_v3 <- if (!is.null(ex_v3$field) && length(ex_v3$field) == length(f_truth))
                      stats::cor(ex_v3$field, f_truth) else NA_real_
    field_cor_jc <- if (!is.null(ex_jc$field) && length(ex_jc$field) == length(f_truth))
                      stats::cor(ex_jc$field, f_truth) else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(
      seed       = s,
      err_psi_v3 = abs(ex_v3$psi_int - psi_truth),
      err_psi_jc = abs(ex_jc$psi_int - psi_truth),
      err_p_v3   = abs(ex_v3$p_int - p_truth),
      err_p_jc   = abs(ex_jc$p_int - p_truth),
      err_pos_v3 = abs(ex_v3$pos_int - pos_truth),
      err_pos_jc = abs(ex_jc$pos_int - pos_truth),
      err_a_v3   = abs(ex_v3$alpha - ALPHA_TRUE),
      err_a_jc   = abs(ex_jc$alpha - ALPHA_TRUE),
      fcor_v3    = field_cor_v3,
      fcor_jc    = field_cor_jc,
      dt_v3      = r_v3$dt,
      dt_jc      = r_jc$dt
    )
    cat(sprintf("v3=%.1fs jc=%.1fs (%.0fx)\n",
                r_v3$dt, r_jc$dt, r_v3$dt / max(r_jc$dt, 1e-3)))
  }

  if (length(rows) == 0L) { cat("(no rows)\n"); return(invisible(NULL)) }
  df <- do.call(rbind, rows)
  cat("\n== Recovery error (|est - truth|) across seeds ==\n")
  metrics <- list(
    list(name = "psi_int",   v3 = df$err_psi_v3, jc = df$err_psi_jc),
    list(name = "p_int",     v3 = df$err_p_v3,   jc = df$err_p_jc),
    list(name = "pos_int",   v3 = df$err_pos_v3, jc = df$err_pos_jc),
    list(name = "alpha",     v3 = df$err_a_v3,   jc = df$err_a_jc)
  )
  cat(sprintf("  %-10s   %12s   %12s   %12s   %12s\n",
              "param", "v3 mean|err|", "jc mean|err|", "v3 max|err|", "jc max|err|"))
  for (m in metrics) {
    cat(sprintf("  %-10s   %12.4f   %12.4f   %12.4f   %12.4f\n",
                m$name, mean(m$v3), mean(m$jc), max(m$v3), max(m$jc)))
  }
  cat(sprintf("\n  %-10s   v3 = %.3f +/- %.3f   jc = %.3f +/- %.3f\n",
              "field cor",
              mean(df$fcor_v3, na.rm = TRUE), stats::sd(df$fcor_v3, na.rm = TRUE),
              mean(df$fcor_jc, na.rm = TRUE), stats::sd(df$fcor_jc, na.rm = TRUE)))
  cat(sprintf("\n  %-10s   v3 = %.1fs +/- %.1f   jc = %.2fs +/- %.2f   speedup = %.0fx\n",
              "wall-clock",
              mean(df$dt_v3), stats::sd(df$dt_v3),
              mean(df$dt_jc), stats::sd(df$dt_jc),
              mean(df$dt_v3) / mean(df$dt_jc)))
  invisible(df)
}

invisible(main(n_seeds = 10L, N_side = 7L, J = 4L))
