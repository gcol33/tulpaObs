# Scaling profile for gcol33/tulpaObs#101.
#
# occu_cover() joint-coupled engine cost vs number of ICAR nodes. Confirms
# whether per-outer-grid-point cost scales super-linearly and isolates the
# Pareto-k diagnostic cost + the R-side post-processing cost from the engine
# inner-Newton grid integration.
#
# Fixed SMALL deterministic outer grid (2 sigma x 2 alpha = 4 cells, no
# adaptive refine, no var-of-means consistency) so per-grid-point cost =
# total / 4 and the cost is comparable across field sizes.

suppressMessages({
  library(tulpa)
  library(tulpaObs)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# 2D rook-adjacency grid on an side x side lattice -> side^2 cells.
grid_adj <- function(side) {
  n <- side * side
  adj <- matrix(0L, n, n)
  idx <- function(r, c) (r - 1L) * side + c
  for (r in seq_len(side)) for (c in seq_len(side)) {
    s <- idx(r, c)
    if (r > 1L)    adj[s, idx(r - 1L, c)] <- 1L
    if (r < side)  adj[s, idx(r + 1L, c)] <- 1L
    if (c > 1L)    adj[s, idx(r, c - 1L)] <- 1L
    if (c < side)  adj[s, idx(r, c + 1L)] <- 1L
  }
  adj
}

build_data <- function(side, J = 3L, seed = 1L) {
  N   <- side * side
  adj <- grid_adj(side)
  sim <- simulate_occu_cover(
    N = N, J = J, positive = "beta", phi = 25,
    adj = adj, sigma = 0.8, alpha = 1.0, seed = seed
  )
  long <- data.frame(
    site_id = rep(seq_len(N), each = J), visit = rep(seq_len(J), times = N),
    y = as.vector(t(sim$y)),
    det_cov1 = sim$visit_data$det_cov1, pos_cov1 = sim$visit_data$pos_cov1
  )
  od <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
  cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
  y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0
  list(adj = adj, od = od, cell_dat = cell_dat, y_pos = y_pos, N = N)
}

fit_once <- function(d, diagnose_k = TRUE, k_samples = 200L,
                     n_threads_outer = 1L) {
  ctrl <- list(
    verbose = FALSE, max.iter = 200L, engine = "joint_coupled",
    sigma.grid = c(0.6, 1.2), alpha.grid = c(0, 1.0),
    adaptive.grid = FALSE, var.of.means.consistency = FALSE,
    diagnose.k = diagnose_k, k.samples = as.integer(k_samples),
    n.threads.outer = as.integer(n_threads_outer),
    progress = FALSE
  )
  suppressWarnings(tobs(
    formula   = ~ occ_cov1 + icar(graph = d$adj), data = d$cell_dat,
    family    = occu_cover("beta"),
    detection = ~ det_cov1, positive = ~ pos_cov1,
    y = d$od$y, y_pos = d$y_pos, visits = d$od$det.covs,
    method = "nested_laplace", control = ctrl
  ))
}

log_line <- function(...) {
  cat(sprintf(...)); cat("\n"); flush(stdout())
}

N_GRID <- 4L  # 2 sigma x 2 alpha

sides <- c(15L, 22L, 30L, 38L)   # 225, 484, 900, 1444 cells

log_line("# build + warmup")
datasets <- lapply(sides, function(s) {
  t0 <- proc.time()[["elapsed"]]
  d <- build_data(s)
  log_line("  built side=%d N=%d in %.1fs", s, d$N, proc.time()[["elapsed"]] - t0)
  d
})

# Warm the JIT / first-call costs on the smallest field.
invisible(fit_once(datasets[[1]], diagnose_k = FALSE))

log_line("\n# scaling: per-grid-point cost vs N (diagnose.k on vs off), serial")
res <- data.frame()
for (i in seq_along(sides)) {
  d <- datasets[[i]]
  t_on  <- system.time(f_on  <- fit_once(d, diagnose_k = TRUE))[["elapsed"]]
  t_off <- system.time(f_off <- fit_once(d, diagnose_k = FALSE))[["elapsed"]]
  row <- data.frame(
    side = sides[i], N = d$N,
    t_on = t_on, t_off = t_off,
    per_gp_on  = t_on  / N_GRID,
    per_gp_off = t_off / N_GRID,
    diag_frac  = (t_on - t_off) / t_on
  )
  res <- rbind(res, row)
  log_line("  N=%4d  t_on=%6.2fs (%.2fs/gp)  t_off=%6.2fs (%.2fs/gp)  diag=%4.0f%%",
           d$N, t_on, row$per_gp_on, t_off, row$per_gp_off, 100 * row$diag_frac)
}

log_line("\n# scaling exponents (log-log slope across consecutive sizes)")
for (i in 2:nrow(res)) {
  e_on  <- log(res$per_gp_on[i]  / res$per_gp_on[i - 1])  / log(res$N[i] / res$N[i - 1])
  e_off <- log(res$per_gp_off[i] / res$per_gp_off[i - 1]) / log(res$N[i] / res$N[i - 1])
  log_line("  N %d->%d : exponent on=%.2f off=%.2f", res$N[i - 1], res$N[i], e_on, e_off)
}

# Rprof the largest size (diagnose.k on) for an R-level breakdown -- catches
# the Pareto-k R code, post-processing vcov block, field demean, etc.
log_line("\n# Rprof breakdown at largest size (diagnose.k = TRUE)")
prof_file <- file.path("dev_notes", "_profile_101.Rprof")
Rprof(prof_file, interval = 0.01, line.profiling = TRUE)
invisible(fit_once(datasets[[length(datasets)]], diagnose_k = TRUE))
Rprof(NULL)
sp <- summaryRprof(prof_file)
log_line("## by.total (top 25)")
print(utils::head(sp$by.total, 25))
log_line("## by.self (top 25)")
print(utils::head(sp$by.self, 25))

saveRDS(res, file.path("dev_notes", "_profile_101_res.rds"))
log_line("\n# DONE")
