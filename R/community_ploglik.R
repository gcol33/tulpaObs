# community_ploglik.R - pointwise log-likelihood for the community-occupancy
# families (ms_occu / ms_dyn_occu / ms_int_occu), so waic() / dic() /
# cpo() score them.
#
# The observation unit is the (species, site) pair: for each species the latent
# occupancy at a site is integrated out over that site's visits / seasons /
# sources using the species' own coefficients (community mean + per-species
# deviation). Each draw pairs a community-mean draw with per-species deviation
# draws taken jointly with it (.tobs_community_joint_draws(); the shared arms --
# dynamic gamma / eps -- carry no deviation).
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

# Per-species eta draws for a process: [n_draws x n_units]. `b_s` ([n_draws x p]
# deviation draws paired row for row with `draws`, or NULL for a shared arm) is
# added to the community-mean draw before the design multiply, so each species
# sees mu + b_s.
.tobs_community_eta <- function(draws, cols, X, b_s = NULL) {
  beta <- draws[, cols, drop = FALSE]
  if (!is.null(b_s)) beta <- beta + b_s
  beta %*% t(X)
}

# Per-species deviation draws paired row for row with community-mean draws.
#
# A species' coefficients are mu + b_s, and mu and b_s are not independent in
# the posterior: the data pin each species' own coefficients more tightly than
# the community mean, so a draw of mu moves b_s the other way. Holding b_s at
# its BLUP while mu varies spreads every species' coefficients by the full
# community-mean uncertainty, wider than that species' posterior, which
# inflates p_waic and p_D and puts WAIC and DIC apart by more than the model
# size explains (#451). The joint draw comes from what the fit carries:
#   * Laplace-EM fits carry Bf_s, the (mu, global) x b_s cross-Hessian block,
#     and Cinv_s = Cov(b_s | mu, y), so
#     b_s | mu ~ N(b_hat_s - Cinv_s Bf_s' (mu - mu_hat), Cinv_s)
#     (.tobs_sbc_community_b_draws()); the noise is drawn on a fixed stream so
#     a criterion is reproducible and the caller's stream is left untouched.
#   * NUTS fits carry the sampled chain, and mu and b_s are read off the same
#     draw (`fit$nuts$par_cols` names the community columns in it).
#   * A fit carrying neither holds b_s at its BLUP.
# `at_mean = TRUE` returns the one-row posterior mean of the same joint draw
# (the DIC plug-in): the conditional mean of b_s at the mean of mu.
# Returns list(draws = [M x K], b = list over species of [M x P]), b's columns
# laid out as `ms_community$blup_idx`.
.tobs_community_joint_draws <- function(object, n.draws = 1000L,
                                        at_mean = FALSE) {
  cm <- object$ms_community
  S  <- object$model$n_species
  bidx <- cm$blup_idx
  P  <- max(unlist(bidx))
  blup_row <- function(s) {
    b <- numeric(P)
    for (k in names(bidx)) b[bidx[[k]]] <- cm[[k]][s, ]
    b
  }
  rep_rows <- function(v, M) matrix(v, M, length(v), byrow = TRUE)

  nd <- object$nuts
  if (!is.null(nd$draws) && !is.null(nd$layout) && !is.null(nd$par_cols)) {
    mu_all <- nd$draws[, nd$par_cols, drop = FALSE]
    colnames(mu_all) <- colnames(object$draws)
    if (at_mean) {
      return(list(draws = matrix(colMeans(mu_all), 1L,
                                 dimnames = list(NULL, colnames(mu_all))),
                  b = lapply(seq_len(S), function(s) rep_rows(blup_row(s), 1L))))
    }
    M <- nrow(mu_all)
    rows <- if (is.null(n.draws) || as.integer(n.draws) >= M) seq_len(M) else
      unique(round(seq(1, M, length.out = as.integer(n.draws))))
    b <- replicate(S, matrix(0, length(rows), P), simplify = FALSE)
    for (r in seq_along(rows)) {
      B <- .ms_ocs_b_from_z(nd$draws[rows[r], ], nd$layout)
      for (s in seq_len(S)) b[[s]][r, ] <- B[s, seq_len(P)]
    }
    return(list(draws = mu_all[rows, , drop = FALSE], b = b))
  }

  draws <- object$draws
  if (at_mean) {
    draws <- matrix(colMeans(draws), 1L, dimnames = list(NULL, colnames(draws)))
  } else {
    draws <- .tobs_subset_draws(draws, n.draws)
  }
  M <- nrow(draws)
  if (is.null(cm$Bf) || is.null(cm$Cinv)) {
    return(list(draws = draws,
                b = lapply(seq_len(S), function(s) rep_rows(blup_row(s), M))))
  }
  mu_hat <- as.numeric(object$means)
  delta  <- sweep(draws, 2L, mu_hat, "-")
  b <- .tobs_with_seed(1L, lapply(seq_len(S), function(s) {
    Bf_s <- cm$Bf[[s]]; Cinv_s <- cm$Cinv[[s]]
    if (nrow(Bf_s) != ncol(draws)) {
      stop(sprintf(paste0("The cross-Hessian block of species %d has %d rows ",
                          "for %d community-mean draw columns."),
                   s, nrow(Bf_s), ncol(draws)), call. = FALSE)
    }
    if (at_mean) {
      sweep(-(delta %*% Bf_s) %*% Cinv_s, 2L, blup_row(s), "+")
    } else {
      .tobs_sbc_community_b_draws(draws, mu_hat, blup_row(s), Bf_s, Cinv_s, M)
    }
  }))
  list(draws = draws, b = b)
}

# Species s's deviation draws on one arm, `key` naming its `blup_idx` block.
.tobs_community_b_arm <- function(object, b, s, key) {
  b[[s]][, object$ms_community$blup_idx[[key]], drop = FALSE]
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
  mt    <- model$model_type

  # A shared areal field and/or latent factors on the occupancy arm enter every
  # species' psi predictor as a fixed offset: per-site for the field, per-(site,
  # species) for the factors. Both are point estimates from the block-coordinate
  # / nested-Laplace fit, so they shift each draw's eta rather than varying with
  # it. The block-coordinate path records the per-site field contribution
  # directly (a varying-coefficient field is sum_k W[i,k] F[u(i),k], not the
  # intercept field alone); the in-tree C++ plain-field path carries the
  # one-node-per-site field itself. Scoring without them would omit the model's
  # spatial / co-occurrence structure from WAIC / LOO / CPO.
  field_off  <- model$occu_field_offset %||% object$spatial_field
  factor_off <- model$occu_factor_offset
  psi <- list(cols = .tobs_community_proc_cols(object, 1L),
              X    = model$X_psi %||% model$X_occ,
              key  = "blup_psi",
              offset = function(s) {
                o <- if (is.null(field_off)) NULL else as.numeric(field_off)
                if (!is.null(factor_off)) o <- (o %||% 0) + factor_off[, s]
                o
              })

  if (identical(mt, "ms_occu")) {
    src <- list(list(cols = .tobs_community_proc_cols(object, 2L),
                     X    = model$X_det,
                     key  = "blup_p",
                     cnt  = .tobs_community_det_counts(model$y, model$valid)))
  } else {                                    # ms_int_occu
    D    <- model$n_sources
    pn   <- model$process_names
    src  <- lapply(seq_len(D), function(d) {
      list(cols = .tobs_community_proc_cols(object, 1L + d),
           X    = model$X_p[[d]],
           key  = paste0("blup_", pn[d]),
           cnt  = .tobs_community_det_counts(model$y[[d]], model$valid[[d]]))
    })
  }
  list(psi = psi, sources = src,
       n_sites = model$n_sites, n_species = model$n_species)
}

# Two-state pointwise log-likelihood, columns species-major:
# [n_draws x (n_species * n_sites)]. Shared by ms_occu and ms_int_occu.
.tobs_ploglik_community_twostate <- function(object, draws, b) {
  sp <- .tobs_community_twostate_spec(object)
  n_sites <- sp$n_sites
  cols <- vector("list", sp$n_species)
  for (s in seq_len(sp$n_species)) {
    eta_psi  <- .tobs_community_eta(draws, sp$psi$cols, sp$psi$X,
                                    .tobs_community_b_arm(object, b, s, sp$psi$key))
    psi_off  <- sp$psi$offset(s)
    if (!is.null(psi_off)) eta_psi <- sweep(eta_psi, 2L, psi_off, "+")
    log_psi  <- stats::plogis(eta_psi,  log.p = TRUE)   # [M x n_sites]
    log_1mps <- stats::plogis(-eta_psi, log.p = TRUE)

    logA     <- matrix(0, nrow(draws), n_sites)         # log P(data | occupied)
    det_any  <- integer(n_sites)                        # detections across sources
    for (d in seq_along(sp$sources)) {
      sd_    <- sp$sources[[d]]
      eta_p  <- .tobs_community_eta(draws, sd_$cols, sd_$X,
                                    .tobs_community_b_arm(object, b, s, sd_$key))
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
.tobs_ploglik_community_dynamic <- function(object, draws, b, n.threads = 1L) {
  model <- object$model
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
    draws_s[, psi1_cols] <- draws[, psi1_cols, drop = FALSE] +
      .tobs_community_b_arm(object, b, s, "blup_psi1")
    draws_s[, p_cols]    <- draws[, p_cols, drop = FALSE] +
      .tobs_community_b_arm(object, b, s, "blup_p")

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

# Community joint occupancy + cover pointwise log-likelihood, columns species-major:
# [n_draws x (n_species * n_sites)]. Per species the exact two-state occu_cover
# marginal (.occu_cover_site_ll, the same source of truth the fit's
# negative-log-posterior uses) is scored over the community-mean pseudo-draws with
# the per-species BLUP deviation plugged into each arm; the shared positive density
# handles beta / lognormal / gaussian.
.tobs_ploglik_community_occu_cover <- function(object, draws, b) {
  model <- object$model
  n_sites <- model$n_sites; n_species <- model$n_species
  mv      <- model$max_visits
  pil <- model$process_info
  p_occ <- pil[[1L]]$p; p_p <- pil[[2L]]$p; p_pos <- pil[[3L]]$p
  occ_cols <- seq_len(p_occ)
  p_cols   <- p_occ + seq_len(p_p)
  pos_cols <- p_occ + p_p + seq_len(p_pos)
  disp_col <- p_occ + p_p + p_pos + 1L
  M <- nrow(draws)
  cols <- vector("list", n_species)
  for (s in seq_len(n_species)) {
    # Per-species response slice; the shared community designs stay on model_s.
    model_s <- model
    model_s$y     <- matrix(model$y[, , s],     n_sites, mv)
    model_s$valid <- matrix(model$valid[, , s], n_sites, mv)
    model_s$y_pos <- matrix(model$y_pos[, , s], n_sites, mv)
    ll_s <- matrix(0, M, n_sites)
    b_occ <- .tobs_community_b_arm(object, b, s, "blup_occ")
    b_p   <- .tobs_community_b_arm(object, b, s, "blup_p")
    b_pos <- .tobs_community_b_arm(object, b, s, "blup_pos")
    for (m in seq_len(M)) {
      bo   <- draws[m, occ_cols] + b_occ[m, ]
      bp   <- draws[m, p_cols]   + b_p[m, ]
      bpos <- draws[m, pos_cols] + b_pos[m, ]
      eta  <- .occu_cover_eta_from_par(model_s, bo, bp, bpos)
      ll_s[m, ] <- .occu_cover_site_ll(model_s, eta$psi, eta$p_mat, eta$ep_mat,
                                       draws[m, disp_col])
    }
    cols[[s]] <- ll_s
  }
  do.call(cbind, cols)
}

# Community binned distance sampling (ms_distance): per-(species, site)
# marginal log-likelihood scored over the community-mean pseudo-draws with the
# per-species BLUP deviation plugged into each arm -- the bins are pooled by
# the same closed-form marginal the fit used (.tobs_ms_distance_engine(),
# R/ms_distance.R), so `sweep(..., value_only = TRUE)$log_lik` already returns
# the per-site value. Half-normal key only; under the hazard key the shared
# log-shape is a community `global` and rides the trailing `log_b` draw column.
.tobs_ploglik_community_ms_distance <- function(object, draws, b) {
  model <- object$model
  n_sites <- model$n_sites; n_species <- model$n_species
  eng <- .tobs_ms_distance_engine(model)
  lam_cols <- .tobs_community_proc_cols(object, 1L)
  sig_cols <- .tobs_community_proc_cols(object, 2L)
  hazard <- eng$hazard
  b_col  <- if (hazard) match("log_b", colnames(draws)) else NA_integer_
  field      <- model$distance_field_offset
  factor_off <- model$distance_factor_offset
  M <- nrow(draws)
  cols <- vector("list", n_species)
  for (s in seq_len(n_species)) {
    bl <- .tobs_community_b_arm(object, b, s, "blup_lambda")
    bs <- .tobs_community_b_arm(object, b, s, "blup_sigma")
    ll_s <- matrix(0, M, n_sites)
    for (m in seq_len(M)) {
      eta_lam <- as.numeric(model$X_processes[[1L]] %*% (draws[m, lam_cols] + bl[m, ]))
      if (!is.null(field))      eta_lam <- eta_lam + as.numeric(field)
      if (!is.null(factor_off)) eta_lam <- eta_lam + factor_off[, s]
      eta_sig <- as.numeric(model$X_processes[[2L]] %*% (draws[m, sig_cols] + bs[m, ]))
      eta_b   <- if (hazard) draws[m, b_col] else 0
      ll_s[m, ] <- eng$sweep(s, eta_lam, eta_sig, eta_b, value_only = TRUE)$log_lik
    }
    cols[[s]] <- ll_s
  }
  do.call(cbind, cols)
}

# Dispatch: community pointwise log-likelihood given community-mean draws and
# the per-species deviation draws paired with them
# (.tobs_community_joint_draws()).
.tobs_ploglik_community <- function(object, draws, b, n.threads = 1L) {
  switch(object$model$model_type,
    ms_occu       = .tobs_ploglik_community_twostate(object, draws, b),
    ms_int_occu   = .tobs_ploglik_community_twostate(object, draws, b),
    ms_dyn_occu   = .tobs_ploglik_community_dynamic(object, draws, b, n.threads),
    ms_occu_cover = .tobs_ploglik_community_occu_cover(object, draws, b),
    ms_distance   = .tobs_ploglik_community_ms_distance(object, draws, b),
    stop("No community pointwise log-likelihood for model_type = '",
         object$model$model_type, "'.", call. = FALSE))
}

# Score n.draws joint (community mean, per-species deviation) draws (WAIC / LOO
# path).
.tobs_ploglik_ms_community <- function(object, n.draws = 1000L, n.threads = 1L) {
  if (is.null(object$draws) || !is.matrix(object$draws)) {
    stop("Community WAIC / LOO needs the community-mean pseudo-draw matrix; ",
         "`object$draws` is missing.", call. = FALSE)
  }
  jd <- .tobs_community_joint_draws(object, n.draws)
  .tobs_ploglik_community(object, jd$draws, jd$b, n.threads)
}

# Pointwise log-likelihood at the community-mean posterior mean (DIC plug-in).
.tobs_community_loglik_at_mean <- function(object) {
  jd <- .tobs_community_joint_draws(object, at_mean = TRUE)
  as.numeric(.tobs_ploglik_community(object, jd$draws, jd$b, 1L))
}
