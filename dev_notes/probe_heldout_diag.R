# Diagnostic: separate (a) observed-site field recovery from (b) held-out
# interpolation, using SCATTERED held-out sites (realistic) on a 2-D grid graph
# (richer neighbourhood than a chain). Also compares nested-Laplace field
# recovery against a plain NUTS spatial fit on the SAME data to localise whether
# weak recovery is the nested engine or the data.
suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE))
set.seed(11)

# 2-D grid graph (rook adjacency), gx by gy.
gx <- 9; gy <- 9; n <- gx * gy
coord <- expand.grid(cx = seq_len(gx), cy = seq_len(gy))
adj <- matrix(0, n, n)
idx_of <- function(i, j) (j - 1) * gx + i
for (i in seq_len(gx)) for (j in seq_len(gy)) {
  a <- idx_of(i, j)
  if (i < gx) { b <- idx_of(i + 1, j); adj[a, b] <- 1; adj[b, a] <- 1 }
  if (j < gy) { b <- idx_of(i, j + 1); adj[a, b] <- 1; adj[b, a] <- 1 }
}

# Smooth 2-D field: plane + bump, centred.
u_true <- 0.9 * scale(coord$cx)[, 1] + 0.7 * scale(coord$cy)[, 1] +
          1.2 * exp(-((coord$cx - 5)^2 + (coord$cy - 5)^2) / 6)
u_true <- u_true - mean(u_true)
psi_true <- plogis(0 + u_true)
z <- rbinom(n, 1, psi_true)
J <- 8; p_det <- 0.5
y_full <- matrix(0L, n, J)
for (i in seq_len(n)) if (z[i]) y_full[i, ] <- rbinom(J, 1, p_det)

# Scattered held-out: every 4th site.
heldout <- seq(2, n, by = 4)
y <- y_full
y[heldout, ] <- NA
obs <- setdiff(seq_len(n), heldout)
d <- data.frame(site = seq_len(n))

recover <- function(fit, label) {
  occ <- fit$nested_laplace$occ_fit
  modes <- occ$modes; wts <- occ$weights
  mode_bar <- as.numeric(crossprod(modes, wts))
  field_bar <- mode_bar[2:(1 + n)]
  psi <- numeric(n)
  for (i in seq_len(n)) psi[i] <- sum(plogis(modes[, 1] + modes[, 1 + i]) * wts)
  cat(sprintf("[%s] cor(field,u) obs=%.3f held=%.3f | cor(psi,true) held=%.3f | MAE held=%.3f\n",
              label,
              cor(field_bar[obs], u_true[obs]),
              cor(field_bar[heldout], u_true[heldout]),
              cor(psi[heldout], psi_true[heldout]),
              mean(abs(psi[heldout] - psi_true[heldout]))))
  invisible(psi)
}

cat("n =", n, " held-out =", length(heldout), " observed =", length(obs), "\n\n")

fit_icar <- tobs(~ 1 + icar(graph = adj), data = d, family = occu(),
                 detection = ~ 1, y = y, method = "nested_laplace",
                 control = list(max.iter = 40L, tol = 1e-5, verbose = FALSE))
recover(fit_icar, "icar")

fit_bym2 <- tobs(~ 1 + bym2(graph = adj), data = d, family = occu(),
                 detection = ~ 1, y = y, method = "nested_laplace",
                 control = list(max.iter = 40L, tol = 1e-5, verbose = FALSE))
recover(fit_bym2, "bym2")

cat("\nDONE\n")
