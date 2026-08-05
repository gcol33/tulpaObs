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

# ---------------------------------------------------------------------------
# Fixed-hyper field NUTS: shared chain-run and fit-assembly tail
#
# The observation families that sample a latent field (abun, removal, distance,
# fp_occu, dyn_abun, areal or temporal) all run the same shape of sampler: the
# parameter vector is `c(coefficients, raw)` where `raw ~ N(0, I)` is the
# whitened field and the field itself is `field_load %*% raw`, with the field
# precision FIXED at its nested-Laplace posterior mean. What differs between
# them is only the C++ entry point, the coefficient names, and how the family's
# own builder turns the posterior-mean coefficients into a `tobs_fit`. The two
# helpers below carry everything either side of that.
# ---------------------------------------------------------------------------

# Run the per-chain sampler, pool the draws, and split them into the
# coefficient block (the leading `n_base` columns, whose posterior mean and
# covariance are reported) and the whitened field block, whose posterior mean
# maps back through `field_load` to the fitted field.
.tobs_nuts_field_draws <- function(run_chain, n_chains, nms, n_base, n_raw,
                                   field_load) {
  chains <- lapply(seq_len(as.integer(n_chains)), run_chain)
  draws  <- do.call(rbind, lapply(chains, `[[`, "draws"))
  colnames(draws) <- nms
  b_idx   <- seq_len(n_base)
  raw_idx <- n_base + seq_len(n_raw)
  par     <- colMeans(draws); names(par) <- nms
  list(chains     = chains,
       draws      = draws,
       nms        = nms,
       par        = par,
       cov        = stats::cov(draws[, b_idx, drop = FALSE]),
       b_idx      = b_idx,
       n_raw      = n_raw,
       field_mean = as.numeric(
         field_load %*% colMeans(draws[, raw_idx, drop = FALSE])))
}

# Attach the pooled draws, posterior moments, sampler diagnostics and the fitted
# field to the `tobs_fit` the family's own builder produced from the posterior
# mean. A non-NULL `temporal` spec routes the field to the temporal slots;
# otherwise it is the spatial field. `fl` is the field loading (its fixed `tau`
# / `rho` are recorded so the fit says which hyperparameters it sampled under).
.tobs_nuts_field_attach <- function(fit, run, log_lik, n_chains, prior_type, fl,
                                    temporal = NULL) {
  b_idx <- run$b_idx
  fit$draws <- run$draws[, b_idx, drop = FALSE]
  fit$means <- run$par[b_idx]
  fit$sds   <- sqrt(pmax(diag(run$cov), 0)); names(fit$sds) <- run$nms[b_idx]
  fit$vcov  <- run$cov
  fit$n_samples <- nrow(run$draws)
  fit$log_prob  <- rep(log_lik, nrow(run$draws))

  accept    <- unlist(lapply(run$chains, `[[`, "accept_prob"))
  divergent <- unlist(lapply(run$chains, `[[`, "divergent"))
  fit$accept_prob <- accept
  fit$divergent   <- divergent
  fit$method      <- "nuts"
  if (is.null(temporal)) {
    fit$spatial_field <- run$field_mean
  } else {
    fit$temporal <- temporal
    fit$temporal_field <- run$field_mean
  }
  fit$nuts <- list(
    accept_prob = accept, divergent = divergent,
    treedepth = as.integer(unlist(lapply(run$chains, `[[`, "treedepth"))),
    epsilon = run$chains[[1L]]$epsilon, n_chains = as.integer(n_chains),
    divergent_total = sum(divergent), tau = fl$tau, rho = fl$rho,
    n_raw = run$n_raw, prior_type = prior_type, fixed_hyper = TRUE)
  fit
}


# ---------------------------------------------------------------------------
# Non-spatial marginal NUTS: shared chain-run and fit-assembly tail
#
# The same five families (abun, removal, distance, fp_occu, dyn_abun) run an
# identically shaped sampler when there is no field: a flat coefficient vector,
# one C++ FullGradFn entry point, and a `run_chain(ch)` closure. Everything
# after that closure -- pooling the chains, naming the draw columns, the
# posterior moments, the RE tail, and the `fit$nuts` diagnostics block -- is the
# same for all of them. The two helpers below carry it, so a change to the
# reported diagnostics is made once.
# ---------------------------------------------------------------------------

# Run the per-chain sampler, pool the draws under `nms`, and collect the
# posterior moments plus the raw sampler diagnostics.
.tobs_count_nuts_run <- function(run_chain, n_chains, nms) {
  chains <- lapply(seq_len(as.integer(n_chains)), run_chain)
  draws  <- do.call(rbind, lapply(chains, `[[`, "draws"))
  colnames(draws) <- nms
  par <- colMeans(draws); names(par) <- nms
  list(chains    = chains,
       draws     = draws,
       nms       = nms,
       par       = par,
       cov       = stats::cov(draws),
       accept    = unlist(lapply(chains, `[[`, "accept_prob")),
       divergent = unlist(lapply(chains, `[[`, "divergent")),
       treedepth = as.integer(unlist(lapply(chains, `[[`, "treedepth"))),
       epsilon   = chains[[1L]]$epsilon)
}

# Replace the family builder's moment-matched draws and NA diagnostics with the
# actual NUTS posterior. `log_lik` is the data log-likelihood at the posterior
# mean (the laplace-path `logLik()` convention); `re_info` attaches the RE tail
# when the fit carries one; `extra` is the family's own `fit$nuts` entries
# (`is_nb`, `K_max`, `re_arm`, the prior scales, ...), appended after the common
# sampler diagnostics.
.tobs_count_nuts_attach <- function(fit, run, log_lik, n_chains, re_info = NULL,
                                    extra = list()) {
  n_draws <- nrow(run$draws)
  fit$draws <- run$draws
  fit <- .tobs_count_nuts_re_finish(fit, run$draws, run$par, run$cov, run$nms,
                                    re_info)
  fit$n_samples   <- n_draws
  fit$log_prob    <- rep(log_lik, n_draws)
  fit$accept_prob <- run$accept
  fit$divergent   <- run$divergent
  fit$treedepth   <- run$treedepth
  fit$epsilon     <- run$epsilon
  fit$method      <- "nuts"
  fit$nuts <- c(list(accept_prob = run$accept, divergent = run$divergent,
                     treedepth = run$treedepth, epsilon = run$epsilon,
                     n_chains = as.integer(n_chains),
                     divergent_total = sum(run$divergent)),
                extra)
  fit
}


# ---------------------------------------------------------------------------
# Shared single-intercept RE wiring for the marginal NUTS families
# (abun, removal, distance, fp_occu, dyn_abun): one grouping factor, intercept
# only, on the state OR detection arm. The C++ side (marginal_count_nuts.h)
# carries the non-centered per-site offset; these helpers resolve the RE term,
# thread it into the spec / warm-start / draw names, and surface fit$re. Slopes
# / correlated / multi-term / both-arm RE stay on the AGHQ Laplace path.
# ---------------------------------------------------------------------------

# Resolve a single intercept RE from the formula `re` terms. NULL = no RE.
# `arms` names the two predictor blocks the term may sit on (state arm first,
# detection arm second) so the count families (abun/removal: lambda/p), the
# open N-mixture and distance (lambda only) and the false-positive occupancy
# family (psi/p11) share one resolver; the returned `arm` is 0 for the first,
# 1 for the second.
.tobs_count_nuts_re_info <- function(re, model, arms = c("lambda", "p")) {
  if (is.null(re) || length(re) == 0L) return(NULL)
  if (inherits(re, "tobs_re")) re <- list(re)
  split <- .tobs_re_split_two_arms(
    re, model, arms[1L], arms[2L],
    sprintf(paste0("method = \"nuts\" with a random effect supports the RE on ",
                   "ONE arm; put it on %s OR %s, or use method = \"laplace\"."),
            arms[1L], arms[2L]))
  a1 <- split[[arms[1L]]]; a2 <- split[[arms[2L]]]
  if (length(a1) && length(a2))
    stop(sprintf(paste0("method = \"nuts\" with a random effect supports the RE ",
                        "on ONE arm; put it on %s OR %s, or use ",
                        "method = \"laplace\"."), arms[1L], arms[2L]),
         call. = FALSE)
  design <- if (length(a1)) a1 else a2
  if (length(design) != 1L || design[[1L]]$n_coefs != 1L ||
      !isTRUE(design[[1L]]$has_intercept))
    stop("method = \"nuts\" supports a single intercept random effect (1|g) on ",
         "one arm; random slopes / multiple grouping factors fit under ",
         "method = \"laplace\" (AGHQ).", call. = FALSE)
  list(arm = if (length(a1)) 0L else 1L,
       arm_tag = if (length(a1)) arms[1L] else arms[2L],
       group = as.integer(design[[1L]]$idx),
       n_groups = as.integer(design[[1L]]$n_groups),
       label = design[[1L]]$group_label %||% "g1")
}

# Merge the RE block fields into the NUTS spec list (re_arm = -1 -> no RE).
.tobs_count_nuts_re_spec <- function(spec, re_info, sigma.logr) {
  if (is.null(re_info)) { spec$re_arm <- -1L; return(spec) }
  spec$re_arm       <- re_info$arm
  spec$re_group     <- re_info$group
  spec$n_re_groups  <- re_info$n_groups
  spec$sigma_re_lsd <- sigma.logr
  spec
}

# Append z (warm-started at 0, unit metric) + log_sigma_re (log(0.5), 0.25
# metric) to the warm-start init.
.tobs_count_nuts_re_init <- function(init, lay, re_info) {
  if (is.null(re_info)) return(init)
  G <- re_info$n_groups
  init$theta0     <- c(init$theta0, rep(0, G), log(0.5))
  init$inv_metric <- c(init$inv_metric[seq_len(lay$total - G - 1L)],
                       rep(1, G), 0.25)
  init
}

# Draw-column names for the RE block (z_1..z_G + log_sigma_<arm>_<label>).
.tobs_count_nuts_re_names <- function(re_info) {
  if (is.null(re_info)) return(character(0))
  c(paste0("re_", re_info$label, "_z", seq_len(re_info$n_groups)),
    paste0("log_sigma_", re_info$arm_tag, "_", re_info$label))
}

# With an RE present, set means/sds/vcov to the full coordinate set (so they
# align with the RE draws columns and the per-process unscaler leaves the RE
# tail untouched) and attach fit$re (sigma + per-group BLUPs on the natural
# scale, b_g = sigma_re * z_g).
.tobs_count_nuts_re_finish <- function(fit, draws, par, cov, nms, re_info) {
  if (is.null(re_info)) return(fit)
  fit$means <- par
  fit$sds   <- sqrt(pmax(diag(cov), 0)); names(fit$sds) <- nms
  fit$vcov  <- cov
  ls_col  <- paste0("log_sigma_", re_info$arm_tag, "_", re_info$label)
  z_cols  <- paste0("re_", re_info$label, "_z", seq_len(re_info$n_groups))
  sig_dr  <- exp(draws[, ls_col])
  blup_dr <- sig_dr * draws[, z_cols, drop = FALSE]
  fit$re <- list(arm = re_info$arm_tag, group_label = re_info$label,
                 n_groups = re_info$n_groups,
                 sigma = mean(sig_dr), sigma_sd = stats::sd(sig_dr),
                 blup = colMeans(blup_dr), blup_sd = apply(blup_dr, 2L, stats::sd))
  fit
}


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
