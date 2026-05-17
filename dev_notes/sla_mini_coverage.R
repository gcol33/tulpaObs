# sla_mini_coverage.R
#
# Minimal Tier-B-style coverage check (slimmed, 30 seeds rather than 200).
# Verifies that simplified-Laplace CIs improve coverage in regimes where
# Gaussian Laplace under-covers — typically small N and near-boundary p.
#
# This is a smoke-level coverage probe, not the production Tier B.
# Production-grade Tier B should run 200 seeds across the full N x psi x p
# grid (see dev_notes/simplified_laplace_derivation.md §8 / plan §3.6).

setwd("C:/Users/Gilles Colling/Documents/dev/tulpaObs")
suppressMessages(devtools::load_all(quiet = TRUE))

set.seed(2026)

run_one <- function(seed, n_sites, n_visits, beta_psi, beta_p, approx) {
  set.seed(seed)
  elev <- rnorm(n_sites)
  data <- data.frame(elev = elev)
  X_psi <- cbind(1, elev)
  eta_psi <- X_psi %*% beta_psi
  psi <- plogis(eta_psi)
  z <- rbinom(n_sites, 1, psi)
  p <- plogis(beta_p)
  y <- matrix(0L, n_sites, n_visits)
  for (i in seq_len(n_sites)) {
    y[i, ] <- if (z[i] == 1) rbinom(n_visits, 1, p) else 0L
  }
  if (mean(rowSums(y) > 0) == 0) return(NULL)
  fit <- tryCatch(
    tobs(
      formula   = ~ elev,
      data      = data,
      family    = occu(),
      detection = ~ 1,
      y         = y,
      approx    = approx,
      control   = list(verbose = FALSE)
    ),
    error = function(e) NULL
  )
  if (is.null(fit) || any(!is.finite(fit$means))) return(NULL)
  ci <- t(apply(fit$draws, 2, quantile, c(0.025, 0.975)))
  list(
    means = fit$means, sds = fit$sds, ci = ci,
    skew  = if (!is.null(fit$skew)) fit$skew else rep(NA_real_, length(fit$means))
  )
}

# Cases that should *favour* SLA: small N, low p, low psi.
cases <- list(
  list(n_sites = 50, n_visits = 3, beta_psi = c(-1, 0.7), beta_p = -1.2, tag = "n50_p0.23"),
  list(n_sites = 80, n_visits = 3, beta_psi = c(-0.5, 0.5), beta_p = -0.8, tag = "n80_p0.31"),
  list(n_sites = 150, n_visits = 3, beta_psi = c(0, 0.8), beta_p = -0.4, tag = "n150_p0.40")
)
seeds <- seq(101, by = 1, length.out = 30)

for (case in cases) {
  truth <- c(psi_int = case$beta_psi[1], psi_elev = case$beta_psi[2],
             p_int = case$beta_p)
  cov_g <- matrix(NA, length(seeds), length(truth))
  cov_s <- matrix(NA, length(seeds), length(truth))
  for (k in seq_along(seeds)) {
    rg <- run_one(seeds[k], case$n_sites, case$n_visits,
                  case$beta_psi, case$beta_p, "gaussian_laplace")
    rs <- run_one(seeds[k], case$n_sites, case$n_visits,
                  case$beta_psi, case$beta_p, "simplified_laplace")
    if (is.null(rg) || is.null(rs)) next
    cov_g[k, ] <- (rg$ci[, 1] <= truth) & (truth <= rg$ci[, 2])
    cov_s[k, ] <- (rs$ci[, 1] <= truth) & (truth <= rs$ci[, 2])
  }
  cov_g_rate <- colMeans(cov_g, na.rm = TRUE)
  cov_s_rate <- colMeans(cov_s, na.rm = TRUE)
  cat(sprintf("\nCase %s (%d sites, %d visits, p~%.2f):\n",
              case$tag, case$n_sites, case$n_visits, plogis(case$beta_p)))
  cat(sprintf("  truth: psi_int=%.2f psi_elev=%.2f p_int=%.2f\n",
              truth[1], truth[2], truth[3]))
  cat("  coverage (gaussian_laplace):", round(cov_g_rate, 2), "\n")
  cat("  coverage (simplified_laplace):", round(cov_s_rate, 2), "\n")
}
