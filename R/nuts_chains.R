# =============================================================================
# nuts_chains.R -- Multi-chain orchestration for the NUTS engine.
#
# The C++ entry point `cpp_occu_fit()` runs a single chain. Multiple chains,
# thinning, and cross-chain convergence diagnostics are assembled here, in R,
# by running the single-chain sampler `n.chains` times with offset seeds and
# pooling the draws. This keeps the C++ kernel chain-count agnostic; the
# pooled result is shaped exactly like one `cpp_occu_fit()` return so the
# downstream unscaling / parameter-naming code is unchanged.
#
# Parallel chains (`n.threads > 1`) use a PSOCK cluster and therefore require
# tulpaObs to be installed (so workers can load it); the default
# `n.threads = 1` runs chains sequentially and works under `load_all()` too.
# =============================================================================

# Keep every `n.thin`-th row of a draw / diagnostic vector or matrix.
.tobs_thin <- function(x, n.thin) {
  if (n.thin <= 1L || is.null(x)) return(x)
  if (is.matrix(x)) {
    x[seq.int(1L, nrow(x), by = n.thin), , drop = FALSE]
  } else {
    x[seq.int(1L, length(x), by = n.thin)]
  }
}

# Run `n.chains` NUTS chains for a prepared `spec` and pool them.
.tobs_run_chains <- function(spec, n.chains = 1L, n.thin = 1L,
                             n.threads = 1L, verbose = TRUE) {
  n.chains  <- max(1L, as.integer(n.chains))
  n.thin    <- max(1L, as.integer(n.thin))
  n.threads <- max(1L, as.integer(n.threads))
  base_seed <- as.integer(spec$seed %||% 42L)

  run_one <- function(chain_id) {
    s <- spec
    s$seed    <- base_seed + chain_id - 1L
    s$verbose <- isTRUE(verbose) && chain_id == 1L   # avoid interleaved logs
    cpp_occu_fit(s)
  }

  if (n.chains == 1L) {
    chains <- list(run_one(1L))
  } else if (n.threads > 1L) {
    cl <- parallel::makeCluster(min(n.threads, n.chains))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, requireNamespace("tulpaObs", quietly = TRUE))
    parallel::clusterExport(cl, c("spec", "base_seed"), envir = environment())
    chains <- parallel::parLapply(cl, seq_len(n.chains), function(chain_id) {
      s <- spec
      s$seed    <- base_seed + chain_id - 1L
      s$verbose <- FALSE
      cpp_occu_fit(s)
    })
  } else {
    chains <- lapply(seq_len(n.chains), run_one)
  }

  # Per-chain draws (engine scale, thinned) for diagnostics; pooled draws for
  # the returned fit.
  per_chain_draws <- lapply(chains, function(ch) .tobs_thin(ch$draws, n.thin))
  pooled <- do.call(rbind, per_chain_draws)

  out <- chains[[1L]]                       # inherit colnames / extra fields
  out$draws       <- pooled
  out$means       <- colMeans(pooled)
  out$n_samples   <- nrow(pooled)
  out$n_chains    <- n.chains
  out$n_thin      <- n.thin
  out$divergent   <- unlist(lapply(chains, function(ch) .tobs_thin(ch$divergent, n.thin)))
  out$accept_prob <- unlist(lapply(chains, function(ch) .tobs_thin(ch$accept_prob, n.thin)))
  out$treedepth   <- unlist(lapply(chains, function(ch) .tobs_thin(ch$treedepth, n.thin)))
  out$epsilon     <- mean(vapply(chains, function(ch) ch$epsilon %||% NA_real_,
                                 numeric(1)), na.rm = TRUE)
  # Row -> chain map, aligned with the pooled draws (which stay in chain-major
  # order through per-process unscaling). Convergence diagnostics (Rhat / ESS)
  # are computed downstream by tulpa::diagnostics() from the named,
  # unscaled draws + this chain_id; see .tobs_fit_model().
  out$chain_id    <- rep(seq_len(n.chains), vapply(per_chain_draws, nrow, integer(1)))
  out
}


# ---------------------------------------------------------------------------
# Cross-chain convergence diagnostics for the in-tree FullGradFn samplers
# (abun / removal / distance / occu_cover / ms_occu_cover ...), which run their
# own per-chain loops in R rather than through cpp_occu_fit + .tobs_fit_model's
# tulpa::diagnostics path. Family-agnostic: every input is a list of
# per-chain draw matrices [N x P]. Single source of truth for split-R-hat /
# bulk-ESS across those paths.
# ---------------------------------------------------------------------------

# Per-parameter biased autocovariance (lags 0..n-1) via the FFT, for the
# effective-sample-size sum.
.tobs_nuts_acov <- function(x) {
  n <- length(x); x <- x - mean(x)
  nf <- 2^ceiling(log2(2 * n))
  f  <- stats::fft(c(x, rep(0, nf - n)))
  ac <- Re(stats::fft(f * Conj(f), inverse = TRUE)) / nf
  ac[seq_len(n)] / n
}

# Split-R-hat and bulk effective sample size per parameter from a list of M
# per-chain draw matrices [N x P] (Vehtari et al. 2021): split each chain in half
# (2M segments of length n), R-hat = sqrt(((n-1)/n) W + B/n) / W) with W the mean
# within-segment variance and B the between-segment variance; ESS = (2M n) / tau,
# tau = 1 + 2 sum rho_t truncated by Geyer's positive-pair rule, rho_t the
# combined autocorrelation 1 - (W - s_t)/var_plus.
.tobs_nuts_rhat_ess <- function(chains) {
  M <- length(chains); P <- ncol(chains[[1L]])
  N <- nrow(chains[[1L]]); n <- N %/% 2L
  if (n < 2L) return(list(rhat = rep(NA_real_, P), ess = rep(NA_real_, P)))
  # 2M split segments, each n draws.
  segs <- vector("list", 2L * M)
  for (m in seq_len(M)) {
    segs[[2L * m - 1L]] <- chains[[m]][seq_len(n), , drop = FALSE]
    segs[[2L * m]]      <- chains[[m]][n + seq_len(n), , drop = FALSE]
  }
  K <- length(segs)
  rhat <- numeric(P); ess <- numeric(P)
  for (p in seq_len(P)) {
    means <- vapply(segs, function(s) mean(s[, p]), 0)
    vars  <- vapply(segs, function(s) stats::var(s[, p]), 0)
    W <- mean(vars); B <- n * stats::var(means)
    var_plus <- ((n - 1) / n) * W + B / n
    rhat[p] <- if (W > 0) sqrt(var_plus / W) else NA_real_
    # combined autocorrelation rho_t (averaged segment autocovariances).
    acovs <- vapply(segs, function(s) .tobs_nuts_acov(s[, p]), numeric(n))  # n x K
    s_t <- rowMeans(acovs)                                                  # mean acov per lag
    rho <- if (var_plus > 0) 1 - (W - s_t) / var_plus else rep(0, n)
    rho[1L] <- 1
    # Geyer initial positive sequence: sum paired (rho_{2k}, rho_{2k+1}) while >0.
    tau <- 1
    t <- 2L
    while (t + 1L <= n) {
      pair <- rho[t] + rho[t + 1L]
      if (pair < 0) break
      tau <- tau + 2 * pair
      t <- t + 2L
    }
    ess[p] <- if (tau > 0) (K * n) / tau else NA_real_
  }
  list(rhat = rhat, ess = ess)
}
