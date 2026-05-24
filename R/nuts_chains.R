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
      tulpaObs:::cpp_occu_fit(s)
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
  # are computed downstream by tulpa::mcmc_diagnostics() from the named,
  # unscaled draws + this chain_id; see .tobs_fit_model().
  out$chain_id    <- rep(seq_len(n.chains), vapply(per_chain_draws, nrow, integer(1)))
  out
}
