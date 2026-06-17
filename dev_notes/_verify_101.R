# Verify the gcol33/tulpaObs#101 fix: diagnose.k now defaults OFF on the
# occu_cover (and occu_multiscale_cover) joint-coupled paths.
#
# Correctness: the Pareto-k diagnostic is a POST-integration pass that only
# attaches `pareto_k`; it does not touch the mode-find or grid integration. So
# the default (diagnose.k = FALSE) and an explicit diagnose.k = TRUE fit must be
# BIT-IDENTICAL in betas / SDs / field. We assert that, plus: default pareto_k
# is NA, opt-in pareto_k is finite, and the default is materially faster.

suppressMessages(devtools::load_all(".", quiet = TRUE))

`%||%` <- function(a, b) if (is.null(a)) b else a

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

log_line <- function(...) { cat(sprintf(...)); cat("\n"); flush(stdout()) }

n_to <- 8L

# ---- single-block ICAR intercept field, 30x30 = 900 cells ----
side <- 30L; N <- side * side; J <- 3L
adj  <- grid_adj(side)
sim  <- simulate_occu_cover(N = N, J = J, positive = "beta", phi = 25,
                            adj = adj, sigma = 0.8, alpha = 1.0, seed = 101L)
long <- data.frame(site_id = rep(seq_len(N), each = J),
                   visit = rep(seq_len(J), times = N),
                   y = as.vector(t(sim$y)),
                   det_cov1 = sim$visit_data$det_cov1,
                   pos_cov1 = sim$visit_data$pos_cov1)
od   <- tobs_data(long, y = "y", site = "site_id", visit = "visit",
                  det.covs = c("det_cov1", "pos_cov1"))
cell_dat <- cbind(data.frame(site_id = seq_len(N)), sim$data)
y_pos <- sim$y_pos; y_pos[is.na(y_pos)] <- 0

fit_sb <- function(diag_k) {
  ctrl <- list(verbose = FALSE, max.iter = 200L, engine = "joint_coupled",
               sigma.grid = c(0.6, 1.2), alpha.grid = c(0, 1.0),
               adaptive.grid = FALSE, var.of.means.consistency = FALSE,
               n.threads.outer = n_to, progress = FALSE)
  if (!is.null(diag_k)) ctrl$diagnose.k <- diag_k
  suppressWarnings <- suppressWarnings
  suppressWarnings(tobs(
    formula = ~ occ_cov1 + icar(graph = adj), data = cell_dat,
    family = occu_cover("beta"), detection = ~ det_cov1, positive = ~ pos_cov1,
    y = od$y, y_pos = y_pos, visits = od$det.covs,
    method = "nested_laplace", control = ctrl))
}

log_line("# single-block (900 cells), %d outer threads", n_to)
t_def <- system.time(f_def <- fit_sb(NULL))[["elapsed"]]      # default -> FALSE now
t_on  <- system.time(f_on  <- fit_sb(TRUE))[["elapsed"]]      # opt-in diagnostic
log_line("  default (diagnose.k unset): %.2fs  pareto_k=%s",
         t_def, format(f_def$joint_fit$pareto_k))
log_line("  diagnose.k = TRUE        : %.2fs  pareto_k=%s",
         t_on, format(f_on$joint_fit$pareto_k))
log_line("  speedup default vs diag : %.2fx", t_on / t_def)

d_means <- max(abs(f_def$means - f_on$means))
d_sds   <- max(abs(f_def$sds   - f_on$sds))
d_field <- max(abs(f_def$spatial_field - f_on$spatial_field))
log_line("  max|d means|=%.3e  max|d sds|=%.3e  max|d field|=%.3e",
         d_means, d_sds, d_field)
stopifnot(
  "default pareto_k must be NA"        = is.na(f_def$joint_fit$pareto_k),
  "opt-in pareto_k must be finite"     = is.finite(f_on$joint_fit$pareto_k),
  "fit must be bit-identical (means)"  = d_means < 1e-9,
  "fit must be bit-identical (sds)"    = d_sds   < 1e-9,
  "fit must be bit-identical (field)"  = d_field < 1e-9,
  "default must be faster"             = t_def < t_on
)
log_line("  OK single-block")

# ---- multi-block || trend path (the issue's repro shape), 22x22 = 484 ----
side2 <- 22L; N2 <- side2 * side2
adj2  <- grid_adj(side2)
sim2  <- simulate_occu_cover(N = N2, J = J, positive = "beta", phi = 25,
                             adj = adj2, sigma = 0.8, alpha = 1.0,
                             trend = TRUE, sigma_trend = 0.7, alpha_trend = 0.9,
                             seed = 202L)
long2 <- data.frame(site_id = rep(seq_len(N2), each = J),
                    visit = rep(seq_len(J), times = N2),
                    y = as.vector(t(sim2$y)),
                    det_cov1 = sim2$visit_data$det_cov1,
                    pos_cov1 = sim2$visit_data$pos_cov1)
od2   <- tobs_data(long2, y = "y", site = "site_id", visit = "visit",
                   det.covs = c("det_cov1", "pos_cov1"))
cell_dat2 <- cbind(data.frame(site_id = seq_len(N2)), sim2$data)
y_pos2 <- sim2$y_pos; y_pos2[is.na(y_pos2)] <- 0

fit_mb <- function(diag_k) {
  # copy(spatial()) supplies the cross-arm coupling formula-natively, so the
  # alpha grids are not passed via control (mutually exclusive). Keep sigma.grid
  # small + CCD for a quick deterministic outer grid.
  ctrl <- list(verbose = FALSE, max.iter = 200L, engine = "joint_coupled",
               sigma.grid = c(0.6, 1.2), adaptive.grid = FALSE,
               var.of.means.consistency = FALSE,
               n.threads.outer = n_to, progress = FALSE)
  if (!is.null(diag_k)) ctrl$diagnose.k <- diag_k
  suppressWarnings(tobs(
    formula = ~ occ_cov1 + icar(graph = adj2) + icar(graph = adj2, weight = time),
    data = cell_dat2, family = occu_cover("beta"),
    detection = ~ det_cov1, positive = ~ pos_cov1 + copy(spatial()),
    y = od2$y, y_pos = y_pos2, visits = od2$det.covs,
    method = "nested_laplace", control = ctrl))
}

log_line("\n# multi-block || trend (484 cells), %d outer threads", n_to)
t_def2 <- system.time(f_def2 <- fit_mb(NULL))[["elapsed"]]
t_on2  <- system.time(f_on2  <- fit_mb(TRUE))[["elapsed"]]
log_line("  default: %.2fs  pareto_k=%s", t_def2, format(f_def2$joint_fit$pareto_k))
log_line("  diag=T : %.2fs  pareto_k=%s", t_on2,  format(f_on2$joint_fit$pareto_k))
log_line("  speedup: %.2fx", t_on2 / t_def2)
d_means2 <- max(abs(f_def2$means - f_on2$means))
d_field2 <- max(abs(f_def2$spatial_field - f_on2$spatial_field))
log_line("  max|d means|=%.3e  max|d field|=%.3e", d_means2, d_field2)
stopifnot(
  "mb default pareto_k must be NA"   = is.na(f_def2$joint_fit$pareto_k),
  "mb opt-in pareto_k must be finite"= is.finite(f_on2$joint_fit$pareto_k),
  "mb fit bit-identical (means)"     = d_means2 < 1e-9,
  "mb fit bit-identical (field)"     = d_field2 < 1e-9,
  "mb default faster"                = t_def2 < t_on2
)
log_line("  OK multi-block trend")
log_line("\n# DONE")
