# community_ploglik.R - pointwise log-likelihood for the community-occupancy
# families (ms_occu / ms_dyn_occu / ms_int_occu), so tobs_waic() / tobs_dic() /
# tobs_cpo() score them.
#
# The observation unit is the (species, site) pair: for each species the latent
# occupancy at a site is integrated out over that site's visits / seasons /
# sources using the species' own coefficients (community mean + per-species
# BLUP deviation). The community-mean posterior uncertainty comes from the
# stored community-mean pseudo-draws (object$draws); the per-species deviations
# enter at their BLUPs (the shared arms -- dynamic gamma / eps -- carry no
# deviation). This is the community generalisation of the single-family Laplace
# WAIC, which likewise scores over the fixed-effect pseudo-draws.
#
# Two marginals cover the three families:
#   * two-state (ms_occu, ms_int_occu): psi * P(data | occupied) + (1 - psi) *
#     1{no detection}. Detection is site-level per source, so a site's visits
#     collapse to per-source detection / non-detection counts. One source is the
#     ms_occu case; D sources sharing psi is ms_int_occu.
#   * HMM forward (ms_dyn_occu): reuses the single-family dynamic kernel
#     (cpp_occu_dynamic_ploglik) via a per-species model shim.
#
# The two-state R marginal is anchored to the C++ single-season kernel: for
# ms_occu (one source) it reproduces cpp_occu_single_ploglik column for column
# (asserted in test-community-diagnostics.R).

# Column layout of the community-mean draw matrix: fixed-effect betas
# concatenated in process order. Returns the 1-based column block for process k.
.tobs_community_proc_cols <- function(object, k) {
  pinfo <- object$process_info %||% object$model$process_info
  off <- if (k == 1L) 0L else sum(vapply(pinfo[seq_len(k - 1L)],
                                         function(p) p$p, integer(1)))
  off + seq_len(pinfo[[k]]$p)
}

# Per-species eta draws for a process: [n_draws x n_units]. `blup_s` (length p,
# or NULL for a shared arm) is added to the community-mean draw before the design
# multiply, so each species sees mu + b_s.
.tobs_community_eta <- function(draws, cols, X, blup_s = NULL) {
  beta <- draws[, cols, drop = FALSE]
  if (!is.null(blup_s)) beta <- sweep(beta, 2L, blup_s, "+")
  beta %*% t(X)
}

# Per-(source, species) site-level detection counts from a community detection
# array [n_sites x J_d x n_species] and its validity mask: k1 = detections,
# k0 = valid non-detections. Returns two [n_sites x n_species] integer matrices.
.tobs_community_det_counts <- function(y_arr, valid_arr) {
  d <- dim(y_arr)
  k1 <- matrix(0L, d[1L], d[3L]); k0 <- matrix(0L, d[1L], d[3L])
  for (s in seq_len(d[3L])) {
    ys <- y_arr[, , s]; vs <- valid_arr[, , s]
    k1[, s] <- rowSums(vs & ys == 1L)
    k0[, s] <- rowSums(vs & ys == 0L)
  }
  list(k1 = k1, k0 = k0)
}

# Assemble the two-state spec (psi arm + D detection sources) for ms_occu (D = 1)
# or ms_int_occu (D >= 1). `sources[[d]]` carries the draw columns, design,
# per-species BLUP, and per-(site, species) detection counts.
.tobs_community_twostate_spec <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  mt    <- model$model_type

  psi <- list(cols = .tobs_community_proc_cols(object, 1L),
              X    = model$X_psi %||% model$X_occ,
              blup = cm$blup_psi)

  if (identical(mt, "ms_occu")) {
    src <- list(list(cols = .tobs_community_proc_cols(object, 2L),
                     X    = model$X_det,
                     blup = cm$blup_p,
                     cnt  = .tobs_community_det_counts(model$y, model$valid)))
  } else {                                    # ms_int_occu
    D    <- model$n_sources
    pn   <- model$process_names
    src  <- lapply(seq_len(D), function(d) {
      list(cols = .tobs_community_proc_cols(object, 1L + d),
           X    = model$X_p[[d]],
           blup = cm[[paste0("blup_", pn[d])]],
           cnt  = .tobs_community_det_counts(model$y[[d]], model$valid[[d]]))
    })
  }
  list(psi = psi, sources = src,
       n_sites = model$n_sites, n_species = model$n_species)
}

# Two-state pointwise log-likelihood, columns species-major:
# [n_draws x (n_species * n_sites)]. Shared by ms_occu and ms_int_occu.
.tobs_ploglik_community_twostate <- function(object, draws) {
  sp <- .tobs_community_twostate_spec(object)
  n_sites <- sp$n_sites
  cols <- vector("list", sp$n_species)
  for (s in seq_len(sp$n_species)) {
    eta_psi  <- .tobs_community_eta(draws, sp$psi$cols, sp$psi$X, sp$psi$blup[s, ])
    log_psi  <- stats::plogis(eta_psi,  log.p = TRUE)   # [M x n_sites]
    log_1mps <- stats::plogis(-eta_psi, log.p = TRUE)

    logA     <- matrix(0, nrow(draws), n_sites)         # log P(data | occupied)
    det_any  <- integer(n_sites)                        # detections across sources
    for (d in seq_along(sp$sources)) {
      sd_    <- sp$sources[[d]]
      eta_p  <- .tobs_community_eta(draws, sd_$cols, sd_$X, sd_$blup[s, ])
      lp     <- stats::plogis(eta_p,  log.p = TRUE)
      l1mp   <- stats::plogis(-eta_p, log.p = TRUE)
      k1     <- sd_$cnt$k1[, s]; k0 <- sd_$cnt$k0[, s]
      logA   <- logA + sweep(lp, 2L, k1, "*") + sweep(l1mp, 2L, k0, "*")
      det_any <- det_any + k1
    }
    # A site with any detection cannot be unoccupied.
    unocc_pen <- ifelse(det_any == 0L, 0, -Inf)

    term_occ   <- log_psi + logA
    term_unocc <- sweep(log_1mps, 2L, unocc_pen, "+")
    cols[[s]]  <- .tobs_logaddexp(term_occ, term_unocc)
  }
  do.call(cbind, cols)
}

# Dynamic community: per-species shim into the single-family HMM-forward kernel
# (cpp_occu_dynamic_ploglik via .tobs_ploglik_dynamic). psi1 / p carry the
# per-species BLUP; gamma / eps are shared community coefficients.
.tobs_ploglik_community_dynamic <- function(object, draws, n.threads = 1L) {
  model <- object$model
  cm    <- object$ms_community
  n_sites <- model$n_sites; n_seasons <- model$n_seasons
  p_psi1  <- model$process_info[[1L]]$p
  p_p     <- model$process_info[[2L]]$p
  psi1_cols <- seq_len(p_psi1)
  p_cols    <- p_psi1 + seq_len(p_p)

  cols <- vector("list", model$n_species)
  for (s in seq_len(model$n_species)) {
    ys <- model$y[, , , s]
    ys[!model$valid[, , , s]] <- -1L                    # NA / invisible visits
    valid <- ys >= 0

    n_visits     <- integer(n_sites * n_seasons)
    any_detected <- logical(n_sites * n_seasons)
    for (i in seq_len(n_sites)) for (t in seq_len(n_seasons)) {
      idx <- (i - 1L) * n_seasons + (t - 1L) + 1L
      raw <- ys[i, , t]; vv <- raw >= 0
      n_visits[idx]     <- sum(vv)
      any_detected[idx] <- any(raw[vv] == 1L)
    }

    draws_s <- draws
    draws_s[, psi1_cols] <- sweep(draws[, psi1_cols, drop = FALSE], 2L,
                                  cm$blup_psi1[s, ], "+")
    draws_s[, p_cols]    <- sweep(draws[, p_cols, drop = FALSE], 2L,
                                  cm$blup_p[s, ], "+")

    shim <- list(model_type = "dynamic",
                 y = ys, n_visits = n_visits, any_detected = any_detected,
                 X_processes = list(model$X_psi1, model$X_p,
                                    model$X_gamma, model$X_eps),
                 process_info = model$process_info,
                 n_sites = n_sites, n_seasons = n_seasons)
    cols[[s]] <- .tobs_ploglik_dynamic(shim, draws_s, n.threads)
  }
  do.call(cbind, cols)
}

# Dispatch: community pointwise log-likelihood given an explicit draw matrix.
.tobs_ploglik_community <- function(object, draws, n.threads = 1L) {
  switch(object$model$model_type,
    ms_occu     = .tobs_ploglik_community_twostate(object, draws),
    ms_int_occu = .tobs_ploglik_community_twostate(object, draws),
    ms_dyn_occu = .tobs_ploglik_community_dynamic(object, draws, n.threads),
    stop("No community pointwise log-likelihood for model_type = '",
         object$model$model_type, "'.", call. = FALSE))
}

# Subsample the community-mean draws to n.draws and score (WAIC / LOO path).
.tobs_ploglik_ms_community <- function(object, n.draws = 1000L, n.threads = 1L) {
  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("Community WAIC / LOO needs the community-mean pseudo-draw matrix; ",
         "`object$draws` is missing.", call. = FALSE)
  }
  if (!is.null(n.draws) && as.integer(n.draws) < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  .tobs_ploglik_community(object, draws, n.threads)
}

# Pointwise log-likelihood at the community-mean posterior mean (DIC plug-in).
.tobs_community_loglik_at_mean <- function(object) {
  mean_draw <- matrix(colMeans(object$draws), nrow = 1L,
                      dimnames = list(NULL, colnames(object$draws)))
  as.numeric(.tobs_ploglik_community(object, mean_draw, 1L))
}
