# Smoke: reproduce the single-season occupancy Laplace bias and test the fix
# (weights = w_i on the detection block). Self-contained, does NOT depend on
# tulpaObs's public API (which is in mid-refactor in the working tree).
#
# Hypothesis: build_single_callbacks() encodes the detection block with all
# n_valid > 0 sites unweighted. Sites with low w_i = P(z=1|y) are treated as
# definitely occupied, biasing p̂ low and ψ̂ high.
#
# Test: replicate the M-step encoding two ways — buggy (current) and fixed
# (weights = w_i on detection) — and feed both into tulpa_em_laplace. Compare
# the modes to truth and to a quick NUTS reference (via spOccupancy as the
# external cross-check, or skip if not installed).

set.seed(42)

# -----------------------------------------------------------------------------
# Simulator (replicates simulate_occu without depending on tulpaObs)
# -----------------------------------------------------------------------------
sim_occu <- function(N = 600, J = 6, beta_occ = c(0.5, 1.2),
                     beta_det = c(0, 0.8), seed = 42) {
  set.seed(seed)
  occ_cov1 <- rnorm(N)
  det_cov1 <- rnorm(N)
  X_occ <- cbind(1, occ_cov1)
  X_det <- cbind(1, det_cov1)
  psi <- plogis(as.vector(X_occ %*% beta_occ))
  p   <- plogis(as.vector(X_det %*% beta_det))
  z   <- rbinom(N, 1, psi)
  y   <- matrix(NA_integer_, N, J)
  for (i in seq_len(N)) y[i, ] <- rbinom(J, 1, z[i] * p[i])
  list(y = y, X_occ = X_occ, X_det = X_det,
       psi = psi, p = p, z = z,
       beta_occ = beta_occ, beta_det = beta_det)
}

# -----------------------------------------------------------------------------
# Callback builder — two variants share the E-step, differ in M-step encoding
# -----------------------------------------------------------------------------
build_cb <- function(sim, fix = FALSE) {
  y <- sim$y; X_occ <- sim$X_occ; X_det <- sim$X_det
  n_sites <- nrow(y); J <- ncol(y)
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)

  n_valid <- integer(n_sites); n_det <- integer(n_sites); any_det <- logical(n_sites)
  for (i in seq_len(n_sites)) {
    valid <- y[i, ] >= 0
    n_valid[i] <- sum(valid)
    n_det[i] <- sum(y[i, valid] == 1)
    any_det[i] <- n_det[i] > 0
  }
  keep <- n_valid > 0

  extract_beta <- function(sub, p) {
    if (is.null(sub)) return(rep(0, p))
    if (!is.null(sub$mode)) return(sub$mode[seq_len(p)])
    rep(0, p)
  }

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    beta_det <- extract_beta(fits$det, p_det)
    psi <- plogis(as.vector(X_occ %*% beta_occ))
    p   <- plogis(as.vector(X_det %*% beta_det))
    w <- numeric(n_sites)
    for (i in seq_len(n_sites)) {
      if (any_det[i]) { w[i] <- 1 }
      else if (n_valid[i] == 0) { w[i] <- psi[i] }
      else {
        prod_1mp <- (1 - p[i])^n_valid[i]
        num <- psi[i] * prod_1mp
        w[i] <- num / (num + (1 - psi[i]))
      }
    }
    list(weights = w)
  }

  m_step_encode <- function(weights, ...) {
    M <- 1000L
    y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
    y_occ <- pmin(pmax(y_occ, 0L), M)
    occ_block <- list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                      family = "binomial")

    if (fix) {
      # FIXED: weight detection block by w_i.
      w_det <- weights
      w_det[any_det] <- 1   # belt-and-braces (E-step already sets this)
      keep_det <- keep & (w_det > 1e-6)
      det_block <- list(
        y         = n_det[keep_det],
        n_trials = n_valid[keep_det],
        X         = X_det[keep_det, , drop = FALSE],
        weights   = w_det[keep_det],
        family    = "binomial"
      )
    } else {
      # BUGGY: include every n_valid > 0 site unweighted.
      det_block <- list(
        y         = n_det[keep],
        n_trials = n_valid[keep],
        X         = X_det[keep, , drop = FALSE],
        family    = "binomial"
      )
    }
    list(occ = occ_block, det = det_block)
  }

  list(e_step = e_step, m_step_encode = m_step_encode)
}

run_one <- function(cb) {
  res <- tulpa::tulpa_em_laplace(
    e_step        = cb$e_step,
    m_step_encode = cb$m_step_encode,
    max_iter      = 100L, tol = 1e-6, damping = 0.3,
    correction    = "none", verbose = FALSE
  )
  list(
    beta_occ = res$fits$occ$mode[1:2],
    beta_det = res$fits$det$mode[1:2],
    n_iter   = res$n_iter,
    converged = res$converged
  )
}

# -----------------------------------------------------------------------------
# Sweep seeds, both variants
# -----------------------------------------------------------------------------
truth_occ <- c(0.5, 1.2)
truth_det <- c(0,   0.8)

cat("Truth: beta_occ = (0.500, 1.200), beta_det = (0.000, 0.800)\n\n")
cat(sprintf("%5s | %-26s | %-26s\n",
            "seed", "Laplace (current, buggy)", "Laplace (FIXED: w on det)"))
cat(sprintf("%5s-+-%-26s-+-%-26s\n",
            "-----", "--------------------------", "--------------------------"))

for (seed in c(42, 99, 123, 7, 2026)) {
  sim <- sim_occu(N = 600, J = 6, beta_occ = truth_occ, beta_det = truth_det,
                  seed = seed)
  bug   <- run_one(build_cb(sim, fix = FALSE))
  fixed <- run_one(build_cb(sim, fix = TRUE))
  cat(sprintf("%5d | psi=(%.2f, %.2f) p=(%.2f, %.2f) | psi=(%.2f, %.2f) p=(%.2f, %.2f)\n",
              seed,
              bug$beta_occ[1],   bug$beta_occ[2],   bug$beta_det[1],   bug$beta_det[2],
              fixed$beta_occ[1], fixed$beta_occ[2], fixed$beta_det[1], fixed$beta_det[2]))
}
