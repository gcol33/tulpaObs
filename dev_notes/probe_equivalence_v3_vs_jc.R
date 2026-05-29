# =============================================================================
# probe_equivalence_v3_vs_jc.R - compare v3_nested (pure R) vs joint_coupled
# (tulpa C++) on the spatial occu_cover() recovery fixture.
#
# Goal: decide whether joint_coupled can become the default for
# `method = "nested_laplace"` and the v3 pure-R fitter can be deleted
# (gcol33/tulpa#32 closing follow-up).
#
# Per seed, fits the same simulation under both engines and reports:
#   * intercept beta estimates (psi, p, pos)
#   * hyperparam estimates (sigma, alpha)
#   * intercept SDs
#   * wall-clock per engine
#   * |delta| / sd for each pair (z-score: <2 = compatible)
#
# Runs across seeds 1..n_seeds for variance, prints a per-seed table and a
# summary mean +/- sd of the absolute z-scores.
# =============================================================================

suppressMessages({
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpa", quiet = TRUE)
  devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE)
})

# ---- fixture builder --------------------------------------------------------

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

sim_one <- function(seed, N_side = 10L, J = 6L,
                    alpha_true = 1.0, sigma_true = 0.8,
                    sigma_pos_true = 0.35) {
  N <- N_side * N_side
  adj <- make_grid_adj(N_side, N_side)
  sim <- simulate_occu_cover(
    N         = N,
    J         = J,
    n_occ_covs = 1L, n_det_covs = 1L, n_pos_covs = 1L,
    positive  = "lognormal",
    sigma_pos = sigma_pos_true,
    adj       = adj,
    sigma     = sigma_true,
    alpha     = alpha_true,
    seed      = seed
  )
  list(sim = sim, adj = adj)
}

# ---- fit harness ------------------------------------------------------------

fit_one <- function(sim_pack, engine, control_extra = list()) {
  sim <- sim_pack$sim
  adj <- sim_pack$adj
  N <- nrow(sim$y); J <- ncol(sim$y)
  visit_df <- data.frame(
    det_cov1 = rep(sim$visit_data$det_cov1, length.out = N * J),
    pos_cov1 = rep(sim$visit_data$pos_cov1, length.out = N * J)
  )
  # jc holds sigma_pos pinned pre-fit unless given a phi.grid.pos axis. v3
  # estimates log_sigma_pos jointly via outer BFGS. To compare them on the
  # same parameter set, integrate sigma_pos on the jc side too.
  if (identical(engine, "joint_coupled") && is.null(control_extra$phi.grid.pos)) {
    control_extra$phi.grid.pos <- seq(0.1, 1.0, length.out = 5)
  }
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
      control   = c(list(engine = engine, verbose = FALSE), control_extra)
    ),
    error = function(e) list(error = conditionMessage(e))
  )
  dt <- proc.time()[["elapsed"]] - t0
  list(fit = fit, dt = dt)
}

# ---- extract comparable summary --------------------------------------------

summarise_fit <- function(fit) {
  if (!is.null(fit$error)) {
    return(list(ok = FALSE, msg = fit$error))
  }
  m <- fit$means
  s <- fit$sds
  pick <- function(nm) c(mean = unname(m[nm]), sd = unname(s[nm]))
  # v3 reports log_sigma / log_sigma_pos; jc reports sigma / phi_pos.
  # Translate to a common sigma / sigma_pos representation via delta method
  # (sd(exp(x)) ~= exp(x_mean) * sd(x)).
  if ("log_sigma" %in% names(m)) {
    lm <- unname(m[["log_sigma"]]); ls <- unname(s[["log_sigma"]])
    sigma_pair <- c(mean = exp(lm), sd = exp(lm) * ls)
  } else if ("sigma" %in% names(m)) {
    sigma_pair <- pick("sigma")
  } else {
    sigma_pair <- c(mean = NA_real_, sd = NA_real_)
  }
  if ("log_sigma_pos" %in% names(m)) {
    lm <- unname(m[["log_sigma_pos"]]); ls <- unname(s[["log_sigma_pos"]])
    sigma_pos_pair <- c(mean = exp(lm), sd = exp(lm) * ls)
  } else if ("phi_pos" %in% names(m)) {
    sigma_pos_pair <- pick("phi_pos")
  } else {
    sigma_pos_pair <- c(mean = NA_real_, sd = NA_real_)
  }
  list(
    ok        = TRUE,
    psi_int   = pick("psi_(Intercept)"),
    p_int     = pick("p_(Intercept)"),
    pos_int   = pick("pos_(Intercept)"),
    sigma     = sigma_pair,
    alpha     = pick("alpha"),
    sigma_pos = sigma_pos_pair
  )
}

# ---- main sweep -------------------------------------------------------------

main <- function(n_seeds = 3L, N_side = 8L, J = 4L) {
  cat("== probe equivalence v3_nested vs joint_coupled ==\n")
  cat(sprintf("Fixture: %d x %d grid (N = %d cells), J = %d visits\n",
              N_side, N_side, N_side * N_side, J))
  cat(sprintf("Truth: sigma = 0.8, alpha = 1.0, lognormal SD = 0.35\n"))
  cat(sprintf("Seeds: 1..%d\n\n", n_seeds))

  rows <- list()
  for (s in seq_len(n_seeds)) {
    cat(sprintf("-- seed %d --\n", s))
    sp <- sim_one(s, N_side = N_side, J = J)
    r_v3 <- fit_one(sp, engine = "v3_nested")
    r_jc <- fit_one(sp, engine = "joint_coupled")

    sm_v3 <- summarise_fit(r_v3$fit)
    sm_jc <- summarise_fit(r_jc$fit)

    if (!sm_v3$ok) { cat("  v3 ERROR:", sm_v3$msg, "\n"); next }
    if (!sm_jc$ok) { cat("  jc ERROR:", sm_jc$msg, "\n"); next }

    params <- c("psi_int", "p_int", "pos_int", "sigma", "alpha", "sigma_pos")
    cat(sprintf("  %-10s  %12s  %12s  %12s  %12s   z\n",
                "param", "v3 mean", "v3 sd", "jc mean", "jc sd"))
    for (p in params) {
      v3 <- sm_v3[[p]]; jc <- sm_jc[[p]]
      pooled_sd <- sqrt((v3["sd"]^2 + jc["sd"]^2) / 2)
      z <- abs(v3["mean"] - jc["mean"]) / max(pooled_sd, 1e-8)
      cat(sprintf("  %-10s  %12.4f  %12.4f  %12.4f  %12.4f  %5.2f\n",
                  p, v3["mean"], v3["sd"], jc["mean"], jc["sd"], z))
      rows[[length(rows) + 1L]] <- list(
        seed = s, param = p, v3_mean = v3["mean"], v3_sd = v3["sd"],
        jc_mean = jc["mean"], jc_sd = jc["sd"], z = z
      )
    }
    cat(sprintf("  wall-clock: v3 = %.1fs   jc = %.1fs   speedup = %.1fx\n\n",
                r_v3$dt, r_jc$dt, r_v3$dt / max(r_jc$dt, 1e-3)))
  }

  if (length(rows) == 0L) { cat("(no rows)\n"); return(invisible(NULL)) }
  df <- do.call(rbind, lapply(rows, as.data.frame))
  cat("\n== Summary across seeds ==\n")
  by_param <- split(df, df$param)
  for (p in names(by_param)) {
    z <- by_param[[p]]$z
    cat(sprintf("  %-10s   mean|z| = %.2f   max|z| = %.2f   ok(<2) = %d/%d\n",
                p, mean(z), max(z), sum(z < 2), length(z)))
  }
  invisible(df)
}

invisible(main(n_seeds = 3L, N_side = 7L, J = 4L))
