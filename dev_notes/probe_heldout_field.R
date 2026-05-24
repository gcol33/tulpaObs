# Probe Part A foundation: is the nested-Laplace latent field reconstructable
# as eta_i = X_i beta + field[spatial_idx(i)] (field-on-eta-scale), so that
# held-out (n_trials = 0) sites interpolate occupancy from neighbours?
#
# Strategy: simulate single-season occupancy on a chain graph with a SMOOTH
# known spatial field. Hold out a block of interior sites (all-NA histories).
# Fit icar via the nested path with those sites flagged heldout_state. Then
# reconstruct psi at held-out sites two ways and compare to the true psi:
#   (1) grid-weighted mode plug-in:  plogis(X beta_bar + field_bar[idx])
#   (2) full marginalisation:        weighted mean over grid of plogis(eta_k)
# Success = reconstructed field correlates with truth AND held-out psi tracks
# the true psi (spatial interpolation is doing real work).
suppressMessages(devtools::load_all("C:/Users/Gilles Colling/Documents/dev/tulpaObs", quiet = TRUE))
set.seed(7)

n <- 60
# Chain adjacency.
adj <- matrix(0, n, n)
for (i in seq_len(n - 1)) { adj[i, i + 1] <- 1; adj[i + 1, i] <- 1 }

# Smooth true field along the chain (sine wave), centred.
u_true <- 1.4 * sin(2 * pi * seq_len(n) / n)
u_true <- u_true - mean(u_true)
beta0  <- 0.0
psi_true <- plogis(beta0 + u_true)
z <- rbinom(n, 1, psi_true)
J <- 6; p_det <- 0.5
y <- matrix(0L, n, J)
for (i in seq_len(n)) if (z[i]) y[i, ] <- rbinom(J, 1, p_det)

# Hold out an interior block: set their histories to NA.
heldout <- 26:35
y[heldout, ] <- NA

d <- data.frame(site = seq_len(n))

fit <- tobs(~ 1 + icar(graph = adj), data = d, family = occu(),
            detection = ~ 1, y = y, method = "nested_laplace",
            control = list(max.iter = 30L, tol = 1e-4, verbose = FALSE))

occ <- fit$nested_laplace$occ_fit
p_occ <- 1L
modes <- occ$modes            # [n_grid x n_x]
wts   <- occ$weights
cat("n_x (mode width):", ncol(modes), " n_grid:", nrow(modes),
    " field comps:", ncol(modes) - p_occ, " (n =", n, ")\n")
cat("theta_grid head:\n"); print(utils::head(occ$theta_grid))

# Grid-weighted mode.
mode_bar  <- as.numeric(crossprod(modes, wts))
beta_bar  <- mode_bar[seq_len(p_occ)]
field_bar <- mode_bar[(p_occ + 1L):(p_occ + n)]    # first n comps (idx = 1:n)

# (1) plug-in psi.
eta_plug <- beta_bar + field_bar
psi_plug <- plogis(eta_plug)

# (2) marginalised psi over grid (first-n field component, no extra scale).
psi_marg <- numeric(n)
for (i in seq_len(n)) {
  eta_k <- modes[, 1] + modes[, p_occ + i]    # beta_k + field_k[i]
  psi_marg[i] <- sum(plogis(eta_k) * wts)
}

cat("\ncor(field_bar, u_true) all sites:",
    round(cor(field_bar, u_true), 3), "\n")
cat("cor(psi_plug, psi_true) HELD-OUT:",
    round(cor(psi_plug[heldout], psi_true[heldout]), 3), "\n")
cat("cor(psi_marg, psi_true) HELD-OUT:",
    round(cor(psi_marg[heldout], psi_true[heldout]), 3), "\n")
cat("\nheld-out comparison (site / true psi / plug / marg):\n")
print(round(data.frame(site = heldout, true = psi_true[heldout],
                       plug = psi_plug[heldout], marg = psi_marg[heldout]), 3))
cat("\nMAE held-out (plug):", round(mean(abs(psi_plug[heldout] - psi_true[heldout])), 3),
    " (marg):", round(mean(abs(psi_marg[heldout] - psi_true[heldout])), 3), "\n")
cat("\nDONE\n")
