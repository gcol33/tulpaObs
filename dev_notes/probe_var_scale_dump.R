# Dump the raw per-cell variance pieces for one seed to find why icar
# under-covers (too narrow) and bym2 over-covers (full width). Reads the engine
# fe / fev / weights via the guarded `tobs.nested.debug` attribute on the
# state posterior. Measure, don't theorize.
#
#   "/c/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/probe_var_scale_dump.R

suppressMessages(devtools::load_all("."))
options(tobs.nested.debug = TRUE)

make_grid <- function(gx = 10, gy = 10, J = 8, p_det = 0.5, seed = 1) {
  set.seed(seed)
  n <- gx * gy
  coord <- expand.grid(cx = seq_len(gx), cy = seq_len(gy))
  adj <- matrix(0, n, n)
  idx_of <- function(i, j) (j - 1) * gx + i
  for (i in seq_len(gx)) for (j in seq_len(gy)) {
    a <- idx_of(i, j)
    if (i < gx) { b <- idx_of(i + 1, j); adj[a, b] <- 1; adj[b, a] <- 1 }
    if (j < gy) { b <- idx_of(i, j + 1); adj[a, b] <- 1; adj[b, a] <- 1 }
  }
  u <- 0.9 * scale(coord$cx)[, 1] + 0.7 * scale(coord$cy)[, 1] +
       1.2 * exp(-((coord$cx - 5)^2 + (coord$cy - 5)^2) / 6)
  u <- u - mean(u)
  psi <- plogis(u)
  z <- rbinom(n, 1, psi)
  y <- matrix(0L, n, J)
  for (i in seq_len(n)) if (z[i]) y[i, ] <- rbinom(J, 1, p_det)
  list(adj = adj, y = y, psi = psi, z = z, n = n)
}

dump_one <- function(model_type, seed = 1) {
  d <- make_grid(seed = seed)
  heldout <- seq(2, d$n, by = 4)
  y <- d$y; y[heldout, ] <- NA
  adj <- d$adj
  form <- if (model_type == "icar") ~ 1 + icar(graph = adj)
          else ~ 1 + bym2(graph = adj)
  fit <- tobs(form, data = data.frame(s = seq_len(d$n)),
              family = occu(), detection = ~ 1, y = y,
              method = "nested_laplace",
              control = list(max.iter = 30L, tol = 1e-5, verbose = FALSE))
  sp  <- predict(fit, type = "state")
  eng <- attr(sp, "engine")

  cat(sprintf("\n########## %s (seed %d) ##########\n", model_type, seed))
  fe <- eng$fitted_eta; fev <- eng$fitted_eta_var; w <- eng$weights
  w <- w / sum(w)
  cat("n_grid:", nrow(fe), " N:", ncol(fe), "\n")
  cat("weights (sorted desc):",
      paste(sprintf("%.3f", sort(w, decreasing = TRUE)), collapse = " "), "\n")
  cat("effective grid size 1/sum(w^2):", sprintf("%.2f", 1 / sum(w^2)), "\n")
  cat("fitted_eta_var range:", sprintf("%.4f .. %.4f", min(fev), max(fev)), "\n")

  detected   <- which(rowSums(d$y, na.rm = TRUE) > 0 & !(seq_len(d$n) %in% heldout))
  nondet_obs <- which(rowSums(d$y, na.rm = TRUE) == 0 & !(seq_len(d$n) %in% heldout))
  decomp <- function(sites, lab) {
    if (!length(sites)) { cat(sprintf("%-12s : (none)\n", lab)); return(invisible()) }
    em_cell <- fe[, sites, drop = FALSE]            # mode per cell
    within  <- colSums(w * fev[, sites, drop = FALSE])      # E_k[var | k]
    cellmean<- colSums(w * em_cell)                          # E_k[eta]
    between <- colSums(w * sweep(em_cell, 2, cellmean)^2)    # var_k(eta mode)
    cat(sprintf("%-12s : within=%.3f  between=%.3f  total=%.3f  (sd=%.3f)\n",
                lab, mean(within), mean(between),
                mean(within + between), sqrt(mean(within + between))))
  }
  decomp(heldout,    "held-out")
  decomp(detected,   "obs detect")
  decomp(nondet_obs, "obs nondet")

  io <- heldout[1:4]
  cat("sample held-out (psi_true / psi / lo / hi / cover):\n")
  for (s in io)
    cat(sprintf("  site %3d: %.3f / %.3f / [%.3f, %.3f] %s\n",
                s, d$psi[s], sp$psi[s], sp$psi_lower[s], sp$psi_upper[s],
                d$psi[s] >= sp$psi_lower[s] & d$psi[s] <= sp$psi_upper[s]))
  invisible(NULL)
}

dump_one("icar", 1)
dump_one("bym2", 1)
