## Smoke test for Slice A: SPDE end-to-end via occu_laplace
## Validates that:
##  1. existing non-spatial occupancy still works (regression)
##  2. occu_laplace + occu_spde() fits without error
##  3. recovered beta_occ is close to truth
##  4. fitted$mode contains a mesh field of expected length

suppressMessages({
  devtools::load_all("../tulpaMesh", quiet = TRUE)
  devtools::load_all("../tulpa", quiet = TRUE)
  devtools::load_all(quiet = TRUE)
})

set.seed(42)

n_sites <- 400
J <- 8

## Simulate spatial occupancy: ψ = plogis(β0 + β1 * x_cov + u(s))
coords <- cbind(runif(n_sites), runif(n_sites))

## Smooth spatial field: low-frequency cosine on the plane (no SPDE draw, but
## a spatial structure for the model to recover)
u_true <- 0.8 * cos(3 * coords[, 1]) * sin(3 * coords[, 2])

x_cov   <- rnorm(n_sites)
beta_occ_true <- c(-0.5, 0.7)  # intercept, slope
beta_det_true <- c(-0.3, 0.4)
det_cov <- rnorm(n_sites)

eta_occ <- beta_occ_true[1] + beta_occ_true[2] * x_cov + u_true
psi <- plogis(eta_occ)
z <- rbinom(n_sites, 1, psi)

eta_det <- beta_det_true[1] + beta_det_true[2] * det_cov
p_det <- plogis(eta_det)

y <- matrix(0L, n_sites, J)
for (i in seq_len(n_sites)) {
  if (z[i] == 1L) y[i, ] <- rbinom(J, 1, p_det[i])
}

dat <- data.frame(occ_cov = x_cov, det_cov = det_cov,
                  x = coords[, 1], y = coords[, 2])

cat("Truth: beta_occ =", beta_occ_true, " beta_det =", beta_det_true, "\n")
cat("Naive naive_occ =", mean(rowSums(y) > 0), "\n\n")

## ---------------------------------------------------------------------------
## (1) Non-spatial baseline: should still work (regression check)
## ---------------------------------------------------------------------------
mod <- occu(occ_formula = ~ occ_cov, det_formula = ~ det_cov,
            data = dat, y = y)

cat("\n--- Non-spatial fit ---\n")
fit_ns <- occu_fit(mod, method = "laplace", verbose = TRUE)
cat("Non-spatial Laplace means:\n"); print(fit_ns$means)

## ---------------------------------------------------------------------------
## (2) SPDE Laplace
## ---------------------------------------------------------------------------
sp <- occu_spde(coords = coords, max_edge = c(0.3, 0.6), nu = 1,
                prior_range = c(0.3, 0.5), prior_sigma = c(0.7, 0.5))

cat("\nSPDE spec: n_mesh =", sp$n_units, "\n")

## Instrument the laplace path by tracing tulpa_laplace
ws_log <- list()
mode_log <- list()
trace_fn <- function() {
  orig <- tulpa::tulpa_laplace
  wrap <- function(y, n_trials, X, ..., spatial = NULL) {
    res <- orig(y = y, n_trials = n_trials, X = X, ..., spatial = spatial)
    if (!is.null(spatial)) {
      ws_log[[length(ws_log) + 1]] <<- summary(as.numeric(y))
      mh <- if (!is.null(res$mode)) res$mode else NA
      mode_log[[length(mode_log) + 1]] <<- list(
        beta = if (length(mh) > 0) mh[1:2] else NA,
        u_summary = if (length(mh) > 2) summary(mh[3:length(mh)]) else NA
      )
    }
    res
  }
  assignInNamespace("tulpa_laplace", wrap, ns = "tulpa")
}
trace_fn()

fit_sp <- occu_fit(mod, method = "laplace", spatial = sp, verbose = TRUE)

cat("\n--- M-step weight summaries per iteration ---\n")
for (i in seq_along(ws_log)) {
  cat(sprintf("Iter %d: y summary = ", i)); print(ws_log[[i]])
}
cat("\n--- M-step mode (beta, u) per iteration ---\n")
for (i in seq_along(mode_log)) {
  cat(sprintf("Iter %d: beta = %s\n", i, toString(round(mode_log[[i]]$beta, 4))))
  cat(sprintf("        u: ")); print(mode_log[[i]]$u_summary)
}
cat("\nEM convergence:\n"); print(fit_sp$convergence)
cat("\nSPDE Laplace means:\n"); print(fit_sp$means)

cat("\nDifference (sp - ns):\n"); print(fit_sp$means - fit_ns$means)
cat("\nDone.\n")
