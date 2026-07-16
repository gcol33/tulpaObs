# Probe: single-season occu() + icar() under nested_laplace, INTERIOR field.
#
# The finding note's "a separate finding worth acting on regardless": every
# "single-season is healthy" claim rested on null-field fixtures (truth
# sigma = 0). The tau grid is [0.3, 30] over 9 cells, i.e. sigma in
# [0.1826, 1.8257], so sigma = 0 is OUTSIDE it -- the marginal piles onto the
# last cell (tau = 30, 60% of the weight) and reports sigma ~= 0.22 because that
# is the smallest field the grid can express, not because it identified
# anything. A null-field fixture cannot tell "correctly shrunk" from "pinned at
# the grid ceiling".
#
# This measures single-season against a field the grid can represent in its
# INTERIOR: sigma_f = 1.0 -> tau = 1.0, which on a 9-cell log grid from 0.3 to
# 30 (ratio 1.7783: 0.3, 0.53, 0.95, 1.69, 3.0, 5.33, 9.49, 16.9, 30) sits at
# cell ~3 of 9. If single-season genuinely identifies the field precision the
# peak cell must land off BOTH boundaries.
devtools::load_all(".")  # run from the repo root

chain_adj <- function(N) {
  a <- matrix(0L, N, N)
  for (s in seq_len(N)) {
    if (s > 1L) a[s, s - 1L] <- 1L
    if (s < N)  a[s, s + 1L] <- 1L
  }
  a
}

smooth_field <- function(N, sd_target, phase) {
  f <- sin(2 * pi * (seq_len(N) / N) + phase)
  f <- f - mean(f)
  f * (sd_target / stats::sd(f))
}

sim <- function(seed, sigma_f, n_cells = 40L, reps = 6L, J = 4L, b_x = 0.5) {
  set.seed(seed)
  adj <- chain_adj(n_cells)
  f0 <- if (sigma_f > 0) smooth_field(n_cells, sigma_f, phase = 0.7) else rep(0, n_cells)
  n_sites <- n_cells * reps
  cell <- rep(seq_len(n_cells), each = reps)
  x <- as.numeric(scale(stats::rnorm(n_sites)))
  z <- stats::rbinom(n_sites, 1L,
                     plogis(stats::qlogis(0.35) + b_x * x + f0[cell]))
  y <- matrix(0L, n_sites, J)
  for (i in seq_len(n_sites)) y[i, ] <- stats::rbinom(J, 1L, z[i] * plogis(0.4))
  list(adj = adj, y = y, f0 = f0, cell = cell, data = data.frame(x = x, cell = cell))
}

fit_one <- function(s) {
  tobs(~ x + icar(graph = s$adj, group_var = "cell"),
       detection = ~ 1, data = s$data, family = occu(), y = s$y,
       method = "nested_laplace",
       control = list(verbose = FALSE, progress = FALSE))
}

n_seed <- 8L
for (sigma_f in c(0, 1.0)) {
  est <- se <- fsd <- fcor <- peak <- ncell <- pkw <- rep(NA_real_, n_seed)
  for (i in seq_len(n_seed)) {
    s <- sim(i, sigma_f)
    f <- try(fit_one(s), silent = TRUE)
    if (inherits(f, "try-error")) {
      cat("seed", i, "ERROR:", conditionMessage(attr(f, "condition")), "\n"); next
    }
    est[i] <- f$means[["psi_x"]]
    se[i]  <- if ("psi_x" %in% names(f$sds)) f$sds[["psi_x"]] else NA_real_
    ff <- as.numeric(f$spatial_field)
    fsd[i] <- stats::sd(ff)
    # spatial_field is per-cell (group_var); compare against the truth surface.
    if (length(ff) == length(s$f0) && sigma_f > 0) fcor[i] <- stats::cor(ff, s$f0)
    tg <- f$nested_laplace$occ_fit$theta_grid
    w  <- f$nested_laplace$occ_fit$weights
    if (!is.null(tg) && !is.null(w)) {
      peak[i]  <- which.max(w)
      ncell[i] <- length(w)
      pkw[i]   <- max(w) / sum(w)
    }
  }
  cat("\n=== sigma_f =", sigma_f, "(tau truth =",
      if (sigma_f > 0) 1 / sigma_f^2 else Inf, ") ===\n")
  cat("slope   median", sprintf("%.4f", stats::median(est, na.rm = TRUE)), " truth 0.5\n")
  cat("se      finite ", sum(is.finite(se)), "/", n_seed, "\n")
  cat("field sd median", sprintf("%.4f", stats::median(fsd, na.rm = TRUE)), "\n")
  if (sigma_f > 0) cat("field cor median", sprintf("%.4f", stats::median(fcor, na.rm = TRUE)),
                       " min", sprintf("%.4f", min(fcor, na.rm = TRUE)), "\n")
  cat("peak cell     ", paste(peak, collapse = ","), " of ",
      paste(unique(ncell[!is.na(ncell)]), collapse = ","), "\n")
  cat("peak weight   ", paste(sprintf("%.3f", pkw), collapse = ","), "\n")
  if (any(is.finite(se))) {
    hit <- abs(est - 0.5) <= 1.96 * se
    cat("coverage      ", sprintf("%.2f", mean(hit, na.rm = TRUE)), "\n")
  }
}
