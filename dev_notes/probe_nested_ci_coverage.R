# Measure 95% CI coverage of the nested-Laplace NA-response prediction.
# psi_true = plogis(u_i); check psi_lower <= psi_true <= psi_upper at held-out
# sites (the prediction targets) and, for contrast, at observed sites (where the
# M-inflation deflates the variance). Loops seeds for a stable estimate.
#
#   "/c/Program Files/R/R-4.6.0/bin/Rscript.exe" dev_notes/probe_nested_ci_coverage.R

suppressMessages(devtools::load_all("."))

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
  list(adj = adj, y = y, psi = psi, n = n)
}

run_one <- function(model_type, seed) {
  d <- make_grid(seed = seed)
  heldout <- seq(2, d$n, by = 4)
  y <- d$y; y[heldout, ] <- NA
  adj <- d$adj
  form <- if (model_type == "icar") ~ 1 + icar(graph = adj)
          else ~ 1 + bym2(graph = adj)
  fit <- tobs(form,
              data = data.frame(s = seq_len(d$n)),
              family = occu(), detection = ~ 1, y = y,
              method = "nested_laplace",
              control = list(max.iter = 30L, tol = 1e-5, verbose = FALSE))
  sp <- predict(fit, type = "state")
  list(sp = sp, psi = d$psi, heldout = heldout)
}

summarise <- function(term, seeds = 1:6) {
  ho_cov <- c(); ob_cov <- c(); ho_w <- c(); ho_cor <- c(); ho_mae <- c()
  for (s in seeds) {
    r <- tryCatch(run_one(term, s), error = function(e) {
      cat(sprintf("  seed %d ERROR: %s\n", s, conditionMessage(e))); NULL
    })
    if (is.null(r)) next
    sp <- r$sp
    ho <- sp$heldout
    cov_i <- r$psi >= sp$psi_lower & r$psi <= sp$psi_upper
    ho_cov <- c(ho_cov, cov_i[ho])
    ob_cov <- c(ob_cov, cov_i[!ho])
    ho_w   <- c(ho_w, (sp$psi_upper - sp$psi_lower)[ho])
    ho_cor <- c(ho_cor, cor(sp$psi[ho], r$psi[ho]))
    ho_mae <- c(ho_mae, mean(abs(sp$psi[ho] - r$psi[ho])))
  }
  cat(sprintf("\n== %s ==\n", term))
  cat(sprintf("  held-out coverage : %.3f  (n=%d)\n", mean(ho_cov), length(ho_cov)))
  cat(sprintf("  observed coverage : %.3f  (n=%d)\n", mean(ob_cov), length(ob_cov)))
  cat(sprintf("  held-out CI width : %.3f (median)\n", median(ho_w)))
  cat(sprintf("  held-out cor/MAE  : %.3f / %.3f\n", mean(ho_cor), mean(ho_mae)))
}

summarise("icar")
summarise("bym2")
