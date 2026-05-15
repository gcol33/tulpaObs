# ============================================================================
# Occupancy-specific callbacks for tulpa's EM+Laplace engine
# ============================================================================

#' Fit a tobs model with Laplace approximation (internal)
#'
#' Uses tulpa's generic EM+Laplace engine with occupancy-specific E-step and
#' M-step encoding. Supports all built-in model types. Called from
#' `.tobs_fit_model()`; not user-facing.
#'
#' @keywords internal
.tobs_laplace <- function(model, spatial = NULL, sigma_beta = 10,
                          max_iter = 50L, tol = 1e-4, damping = 0.3,
                          correction = c("auto", "mi", "gibbs", "none"),
                          n_imputations = 20L,
                          verbose = TRUE) {
  correction <- match.arg(correction)
  if (!inherits(model, "tobs_model")) stop("model must be a tobs_model object")

  .validate_spatial_laplace(spatial, model$model_type)

  callbacks <- switch(model$model_type,
    single     = build_single_callbacks(model, spatial),
    dynamic    = build_dynamic_callbacks(model, spatial),
    community  = build_community_callbacks(model, spatial),
    integrated = build_integrated_callbacks(model, spatial),
    jsdm       = build_jsdm_callbacks(model, spatial),
    stop(sprintf("Laplace not supported for model_type '%s'", model$model_type))
  )

  em_result <- tulpa::tulpa_em_laplace(
    e_step        = callbacks$e_step,
    m_step_encode = callbacks$m_step_encode,
    draw_z        = callbacks$z_draw,
    max_iter      = max_iter,
    tol           = tol,
    damping       = damping,
    correction    = correction,
    n_imputations = n_imputations,
    verbose       = verbose
  )

  build_laplace_fit(em_result, model, spatial, callbacks$p_per_submodel)
}

# ============================================================================
# Single-season callbacks
# ============================================================================
build_single_callbacks <- function(model, spatial = NULL) {
  y <- model$y
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  n_sites <- model$n_sites
  p_occ <- ncol(X_occ)
  p_det <- ncol(X_det)

  n_valid <- integer(n_sites)
  n_det <- integer(n_sites)
  any_det <- logical(n_sites)
  for (i in seq_len(n_sites)) {
    valid <- y[i, ] >= 0
    n_valid[i] <- sum(valid)
    n_det[i] <- sum(y[i, valid] == 1)
    any_det[i] <- n_det[i] > 0
  }
  keep <- n_valid > 0

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    beta_det <- extract_beta(fits$det, p_det)
    eta_occ <- as.vector(X_occ %*% beta_occ)
    sp_off <- .spatial_eta_offset(spatial, fits$occ, p_occ)
    if (length(sp_off) == n_sites) eta_occ <- eta_occ + sp_off
    psi <- plogis(eta_occ)
    p <- plogis(as.vector(X_det %*% beta_det))
    list(weights = occ_weights(psi, p, n_sites, n_valid, n_det, any_det))
  }

  m_step_encode <- function(weights, ...) {
    if (is.null(spatial)) {
      # Pseudo-binomial encoding: y = round(M*w), n_trials = M. The
      # M-inflation makes the M-step into a sharp binomial whose mode
      # equals the weighted mean, and is the historical encoding used
      # everywhere else in the package.
      M <- 1000L
      y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
      y_occ <- pmin(pmax(y_occ, 0L), M)
      occ_block <- list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                        family = "binomial")
    } else {
      # Modest pseudo-binomial encoding for SPDE: M = 4 gives some
      # fractional resolution on the weights while keeping the per-site
      # effective sample size O(1), so the SPDE prior precision is not
      # swamped by the data signal as it would be at M = 1000.
      M <- 4L
      y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
      y_occ <- pmin(pmax(y_occ, 0L), M)
      occ_block <- list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                        family = "binomial")
      occ_block <- .attach_spatial_spde(occ_block, spatial)
    }
    # Detection block: weight by w_i = P(z_i = 1 | y_i, theta). Sites that
    # the E-step thinks are likely empty (w_i ~ 0) must drop out of the
    # detection fit, otherwise they bias p_hat downward by feeding their
    # all-zero detection history as evidence about (1 - p)^J. Sites with
    # any detection have w_i = 1 (the E-step sets this).
    w_det <- weights
    w_det[any_det] <- 1
    keep_det <- keep & (w_det > 1e-6)
    list(
      occ = occ_block,
      det = list(y = n_det[keep_det], n_trials = n_valid[keep_det],
                 X = X_det[keep_det, , drop = FALSE],
                 weights = w_det[keep_det], family = "binomial")
    )
  }

  z_draw <- function(weights, ...) {
    z <- as.integer(any_det)
    z[!any_det] <- rbinom(sum(!any_det), 1, clamp_w(weights[!any_det]))
    z
  }

  hard_encode <- function(z, ...) {
    occ_sites <- which(z == 1L)
    det_keep <- occ_sites[n_valid[occ_sites] > 0]
    occ_block <- list(y = z, n_trials = rep(1L, n_sites), X = X_occ,
                      family = "binomial")
    occ_block <- .attach_spatial_spde(occ_block, spatial)
    list(
      occ = occ_block,
      det = if (length(det_keep) > 0)
        list(y = n_det[det_keep], n_trials = n_valid[det_keep],
             X = X_det[det_keep, , drop = FALSE], family = "binomial")
      else NULL
    )
  }

  init <- glm_init(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det)

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init, p_per_submodel = c(occ = p_occ, det = p_det))
}

# ============================================================================
# Dynamic occupancy callbacks
# ============================================================================
build_dynamic_callbacks <- function(model, spatial = NULL) {
  # spatial guaranteed NULL here by .validate_spatial_laplace.
  y_flat <- model$y_flat
  n_sites <- model$n_sites
  n_seasons <- model$n_seasons
  max_visits <- model$max_visits
  X_occ <- model$X_processes[[1]]  # psi1
  X_det <- model$X_processes[[2]]  # p
  X_col <- model$X_processes[[3]]  # gamma
  X_ext <- model$X_processes[[4]]  # epsilon
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)
  p_col <- ncol(X_col); p_ext <- ncol(X_ext)

  # Precompute per site-season
  nv <- model$n_visits
  ad <- model$any_detected

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    beta_det <- extract_beta(fits$det, p_det)
    beta_col <- extract_beta(fits$col, p_col)
    beta_ext <- extract_beta(fits$ext, p_ext)

    psi1 <- plogis(as.vector(X_occ %*% beta_occ))
    p <- plogis(as.vector(X_det %*% beta_det))
    gam <- plogis(as.vector(X_col %*% beta_col))
    eps <- plogis(as.vector(X_ext %*% beta_ext))

    # HMM forward pass to get P(z_it = 1 | y)
    w <- matrix(NA_real_, n_sites, n_seasons)
    for (i in seq_len(n_sites)) {
      alpha_occ <- psi1[i]; alpha_unocc <- 1 - psi1[i]
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        nv_it <- nv[idx]; det_it <- ad[idx]
        if (nv_it > 0) {
          prob_y_occ <- if (det_it) 1 else (1 - p[i])^nv_it
          prob_y_unocc <- if (det_it) 0 else 1
          post_occ <- alpha_occ * prob_y_occ
          post_unocc <- alpha_unocc * prob_y_unocc
          total <- post_occ + post_unocc
          w[i, t] <- post_occ / total
          alpha_occ <- post_occ / total
          alpha_unocc <- post_unocc / total
        } else {
          w[i, t] <- alpha_occ
        }
        if (t < n_seasons) {
          new_occ <- alpha_occ * (1 - eps[i]) + alpha_unocc * gam[i]
          new_unocc <- alpha_occ * eps[i] + alpha_unocc * (1 - gam[i])
          alpha_occ <- new_occ; alpha_unocc <- new_unocc
        }
      }
    }
    list(weights = w)
  }

  m_step_encode <- function(weights, ...) {
    w <- weights  # n_sites x n_seasons matrix
    # Occupancy: psi1 from season 1 weights
    M <- 1000L
    w1 <- w[, 1]
    y_occ <- ifelse(ad[seq(1, by = n_seasons, length.out = n_sites)], M,
                    as.integer(round(w1 * M)))
    y_occ <- pmin(pmax(y_occ, 0L), M)

    # Colonization: from transitions where z_{t-1}=0, z_t=1
    # Extinction: from transitions where z_{t-1}=1, z_t=0
    # Approximate: site-level average
    col_y <- integer(n_sites); col_n <- integer(n_sites)
    ext_y <- integer(n_sites); ext_n <- integer(n_sites)
    for (i in seq_len(n_sites)) {
      for (t in 2:n_seasons) {
        p_prev_occ <- w[i, t - 1]
        p_curr_occ <- w[i, t]
        # P(colonization event) ≈ (1-w_{t-1}) * w_t
        col_y[i] <- col_y[i] + as.integer(round((1 - p_prev_occ) * p_curr_occ * M))
        col_n[i] <- col_n[i] + as.integer(round((1 - p_prev_occ) * M))
        ext_y[i] <- ext_y[i] + as.integer(round(p_prev_occ * (1 - p_curr_occ) * M))
        ext_n[i] <- ext_n[i] + as.integer(round(p_prev_occ * M))
      }
    }
    col_n <- pmax(col_n, 1L); ext_n <- pmax(ext_n, 1L)
    col_y <- pmin(pmax(col_y, 0L), col_n)
    ext_y <- pmin(pmax(ext_y, 0L), ext_n)

    # Detection: per-(site, season) rows weighted by w[i, t] = P(z_it = 1 | y).
    # Replaces the legacy hard threshold (w > 0.5) which silently dropped
    # site-seasons in the boundary regime and double-counted detection
    # evidence for site-seasons in the high-confidence regime. X_det is
    # site-indexed in this model, so per-season rows just replicate the
    # site's covariates.
    rows_i <- integer(n_sites * n_seasons)
    det_count <- integer(n_sites * n_seasons)
    vis_count <- integer(n_sites * n_seasons)
    w_it <- numeric(n_sites * n_seasons)
    n_rows <- 0L
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        if (nv[idx] <= 0) next
        base <- (i - 1) * n_seasons * max_visits + (t - 1) * max_visits
        dc <- 0L; vc <- 0L
        for (j in seq_len(nv[idx])) {
          v <- y_flat[base + j]
          if (v >= 0) { vc <- vc + 1L; if (v == 1) dc <- dc + 1L }
        }
        if (vc == 0L) next
        w_eff <- if (dc > 0L) 1 else w[i, t]
        if (w_eff <= 1e-6) next
        n_rows <- n_rows + 1L
        rows_i[n_rows] <- i
        det_count[n_rows] <- dc
        vis_count[n_rows] <- vc
        w_it[n_rows] <- w_eff
      }
    }
    if (n_rows > 0L) {
      rows_i <- rows_i[seq_len(n_rows)]
      det_count <- det_count[seq_len(n_rows)]
      vis_count <- vis_count[seq_len(n_rows)]
      w_it <- w_it[seq_len(n_rows)]
      det_block <- list(y = det_count, n_trials = vis_count,
                        X = X_det[rows_i, , drop = FALSE],
                        weights = w_it, family = "binomial")
    } else {
      det_block <- list(y = integer(0), n_trials = integer(0),
                        X = X_det[integer(0), , drop = FALSE],
                        weights = numeric(0), family = "binomial")
    }

    list(
      occ = list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                 family = "binomial"),
      det = det_block,
      col = list(y = col_y, n_trials = col_n, X = X_col,
                 family = "binomial"),
      ext = list(y = ext_y, n_trials = ext_n, X = X_ext,
                 family = "binomial")
    )
  }

  z_draw <- function(weights, ...) {
    w <- weights
    z <- matrix(0L, n_sites, n_seasons)
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        idx <- (i - 1) * n_seasons + t
        if (ad[idx]) z[i, t] <- 1L
        else z[i, t] <- rbinom(1, 1, clamp_w(w[i, t]))
      }
    }
    z
  }

  hard_encode <- function(z, ...) {
    z1 <- z[, 1]
    # Colonization/extinction from hard transitions
    col_y <- integer(n_sites); col_n <- integer(n_sites)
    ext_y <- integer(n_sites); ext_n <- integer(n_sites)
    for (i in seq_len(n_sites)) {
      for (t in 2:n_seasons) {
        if (z[i, t - 1] == 0) { col_n[i] <- col_n[i] + 1L; if (z[i, t] == 1) col_y[i] <- col_y[i] + 1L }
        if (z[i, t - 1] == 1) { ext_n[i] <- ext_n[i] + 1L; if (z[i, t] == 0) ext_y[i] <- ext_y[i] + 1L }
      }
    }
    col_n <- pmax(col_n, 1L); ext_n <- pmax(ext_n, 1L)

    total_det <- integer(n_sites); total_vis <- integer(n_sites)
    for (i in seq_len(n_sites)) {
      for (t in seq_len(n_seasons)) {
        if (z[i, t] == 1 && nv[(i-1)*n_seasons+t] > 0) {
          base <- (i-1)*n_seasons*max_visits + (t-1)*max_visits
          for (j in seq_len(nv[(i-1)*n_seasons+t])) {
            v <- y_flat[base + j]
            if (v >= 0) { total_vis[i] <- total_vis[i]+1L; if (v==1) total_det[i] <- total_det[i]+1L }
          }
        }
      }
    }
    dk <- total_vis > 0

    list(
      occ = list(y = z1, n_trials = rep(1L, n_sites), X = X_occ,
                 family = "binomial"),
      det = if (sum(dk) > 0)
        list(y = total_det[dk], n_trials = total_vis[dk],
             X = X_det[dk,,drop=FALSE], family = "binomial")
      else NULL,
      col = list(y = col_y, n_trials = col_n, X = X_col,
                 family = "binomial"),
      ext = list(y = ext_y, n_trials = ext_n, X = X_ext,
                 family = "binomial")
    )
  }

  init <- list(
    occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)),
    det = list(beta = rep(0, p_det), se = rep(1, p_det)),
    col = list(beta = rep(0, p_col), se = rep(1, p_col)),
    ext = list(beta = rep(0, p_ext), se = rep(1, p_ext))
  )

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init,
       p_per_submodel = c(occ = p_occ, det = p_det, col = p_col, ext = p_ext))
}

# ============================================================================
# Community occupancy callbacks
# ============================================================================
build_community_callbacks <- function(model, spatial = NULL) {
  # spatial guaranteed NULL here by .validate_spatial_laplace.
  y <- model$y  # N x max_visits (expanded: site-species rows)
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  N <- model$N  # n_sites * n_species
  n_sites <- model$n_sites
  n_species <- model$n_species
  max_visits <- model$max_visits
  p_occ <- ncol(X_occ); p_det <- ncol(X_det)

  n_valid <- integer(N); n_det <- integer(N); any_det <- logical(N)
  for (i in seq_len(N)) {
    valid <- y[i, ] >= 0
    n_valid[i] <- sum(valid)
    n_det[i] <- sum(y[i, valid] == 1)
    any_det[i] <- n_det[i] > 0
  }
  keep <- n_valid > 0

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    beta_det <- extract_beta(fits$det, p_det)
    psi <- plogis(as.vector(X_occ %*% beta_occ))
    p <- plogis(as.vector(X_det %*% beta_det))
    list(weights = occ_weights(psi, p, N, n_valid, n_det, any_det))
  }

  m_step_encode <- function(weights, ...) {
    M <- 1000L
    y_occ <- ifelse(any_det, M, as.integer(round(weights * M)))
    y_occ <- pmin(pmax(y_occ, 0L), M)
    # Weight detection rows by E-step occupancy posterior. Same fix as
    # build_single_callbacks(): species-site rows with low w_i drop out
    # of the detection fit.
    w_det <- weights
    w_det[any_det] <- 1
    keep_det <- keep & (w_det > 1e-6)
    list(
      occ = list(y = y_occ, n_trials = rep(M, N), X = X_occ,
                 family = "binomial"),
      det = list(y = n_det[keep_det], n_trials = n_valid[keep_det],
                 X = X_det[keep_det, , drop = FALSE],
                 weights = w_det[keep_det], family = "binomial")
    )
  }

  z_draw <- function(weights, ...) {
    z <- as.integer(any_det)
    z[!any_det] <- rbinom(sum(!any_det), 1, clamp_w(weights[!any_det]))
    z
  }

  hard_encode <- function(z, ...) {
    occ_obs <- which(z == 1L)
    det_keep <- occ_obs[n_valid[occ_obs] > 0]
    list(
      occ = list(y = z, n_trials = rep(1L, N), X = X_occ,
                 family = "binomial"),
      det = if (length(det_keep) > 0)
        list(y = n_det[det_keep], n_trials = n_valid[det_keep],
             X = X_det[det_keep,,drop=FALSE], family = "binomial")
      else NULL
    )
  }

  init <- list(
    occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)),
    det = list(beta = rep(0, p_det), se = rep(1, p_det))
  )

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init,
       p_per_submodel = c(occ = p_occ, det = p_det))
}

# ============================================================================
# Integrated occupancy callbacks
# ============================================================================
build_integrated_callbacks <- function(model, spatial = NULL) {
  if (!is.null(spatial)) {
    stop("Integrated occupancy with SPDE is not yet plumbed in .tobs_laplace.",
         call. = FALSE)
  }
  y_sources <- model$y_sources
  site_maps <- model$site_maps
  X_occ <- model$X_processes[[1]]
  n_sites <- model$n_sites
  n_sources <- model$n_sources
  p_occ <- ncol(X_occ)

  # Per-source detection info
  src_info <- lapply(seq_len(n_sources), function(s) {
    ys <- y_sources[[s]]; ns <- nrow(ys); mv <- ncol(ys)
    nv <- integer(ns); nd <- integer(ns); ad <- logical(ns)
    for (i in seq_len(ns)) {
      valid <- ys[i, ] >= 0; nv[i] <- sum(valid)
      nd[i] <- sum(ys[i, valid] == 1); ad[i] <- nd[i] > 0
    }
    X_det <- model$X_processes[[1 + s]]
    src_rows <- site_maps[[s]] + 1L
    list(nv = nv, nd = nd, ad = ad, X_det = X_det[src_rows, , drop = FALSE],
         p_det = ncol(X_det), keep = nv > 0, src_rows = src_rows)
  })

  # Global detection status per site
  any_det_global <- logical(n_sites)
  for (s in seq_len(n_sources)) {
    for (j in seq_along(src_info[[s]]$src_rows)) {
      if (src_info[[s]]$ad[j]) any_det_global[src_info[[s]]$src_rows[j]] <- TRUE
    }
  }

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    psi <- plogis(as.vector(X_occ %*% beta_occ))
    weights <- psi  # Prior occupancy
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      beta_det <- extract_beta(fits[[paste0("det", s)]], si$p_det)
      p_s <- plogis(as.vector(si$X_det %*% beta_det))
      for (j in seq_along(si$src_rows)) {
        i <- si$src_rows[j]
        if (si$ad[j]) { weights[i] <- 1 }
        else if (si$nv[j] > 0) {
          prod_1mp <- (1 - p_s[j])^si$nv[j]
          num <- weights[i] * prod_1mp
          weights[i] <- num / (num + (1 - weights[i]) + 1e-10)
        }
      }
    }
    list(weights = weights)
  }

  m_step_encode <- function(weights, ...) {
    M <- 1000L
    y_occ <- ifelse(any_det_global, M, as.integer(round(weights * M)))
    y_occ <- pmin(pmax(y_occ, 0L), M)
    specs <- list(occ = list(y = y_occ, n_trials = rep(M, n_sites), X = X_occ,
                             family = "binomial"))
    # Per-source detection blocks: weight each row by w_i at the global
    # site mapped through src_rows. Sites where the E-step says "almost
    # certainly empty" drop out of every source's detection fit.
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      w_src <- weights[si$src_rows]
      w_src[si$ad] <- 1
      dk <- si$keep & (w_src > 1e-6)
      specs[[paste0("det", s)]] <- list(y = si$nd[dk], n_trials = si$nv[dk],
                                        X = si$X_det[dk,,drop=FALSE],
                                        weights = w_src[dk],
                                        family = "binomial")
    }
    specs
  }

  z_draw <- function(weights, ...) {
    z <- as.integer(any_det_global)
    undet <- !any_det_global
    z[undet] <- rbinom(sum(undet), 1, clamp_w(weights[undet]))
    z
  }

  hard_encode <- function(z, ...) {
    specs <- list(occ = list(y = z, n_trials = rep(1L, n_sites), X = X_occ,
                             family = "binomial"))
    for (s in seq_len(n_sources)) {
      si <- src_info[[s]]
      occ_local <- z[si$src_rows] == 1L & si$nv > 0
      if (any(occ_local)) {
        specs[[paste0("det", s)]] <- list(y = si$nd[occ_local],
                                          n_trials = si$nv[occ_local],
                                          X = si$X_det[occ_local,,drop=FALSE],
                                          family = "binomial")
      }
    }
    specs
  }

  init <- list(occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)))
  p_sub <- c(occ = p_occ)
  for (s in seq_len(n_sources)) {
    nm <- paste0("det", s)
    init[[nm]] <- list(beta = rep(0, src_info[[s]]$p_det), se = rep(1, src_info[[s]]$p_det))
    p_sub[nm] <- src_info[[s]]$p_det
  }

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init, p_per_submodel = p_sub)
}

# ============================================================================
# JSDM callbacks (no detection — simple Bernoulli)
# ============================================================================
build_jsdm_callbacks <- function(model, spatial = NULL) {
  # JSDM is N = n_sites * n_species; SPDE A_x is n_sites-indexed, so attaching
  # would require row-expansion. Deferred.
  if (!is.null(spatial)) {
    stop("JSDM with SPDE is not yet plumbed in .tobs_laplace.", call. = FALSE)
  }
  y_jsdm <- model$y_jsdm
  X_occ <- model$X_processes[[1]]
  N <- model$N
  p_occ <- ncol(X_occ)

  # No E-step needed — no latent variable (y is observed directly)
  # But we still use EM framework for consistency with species RE
  e_step <- function(fits, ...) {
    list(weights = as.numeric(y_jsdm))
  }

  m_step_encode <- function(weights, ...) {
    list(occ = list(y = as.integer(y_jsdm), n_trials = rep(1L, N), X = X_occ,
                    family = "binomial"))
  }

  z_draw <- function(weights, ...) as.integer(y_jsdm)
  hard_encode <- function(z, ...) {
    list(occ = list(y = z, n_trials = rep(1L, N), X = X_occ,
                    family = "binomial"))
  }

  init <- list(occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)))

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init, p_per_submodel = c(occ = p_occ))
}

# ============================================================================
# Shared helpers
# ============================================================================

# Validate that `spatial` (a `tobs_spatial` or NULL) can be consumed by the
# Laplace path. Slice A: SPDE on the occupancy/state submodel only, single +
# integrated + jsdm model types. Other combinations error explicitly rather
# than silently dropping the spec.
.validate_spatial_laplace <- function(spatial, model_type) {
  if (is.null(spatial)) return(invisible())
  if (!inherits(spatial, "tobs_spatial")) {
    stop("spatial must be a tobs_spatial object (from tobs_spde(), etc.)",
         call. = FALSE)
  }
  if (!identical(spatial$type, "spde")) {
    stop(sprintf(
      ".tobs_laplace currently supports spatial$type == 'spde' only (got '%s'). Use engine = 'nuts' for other spatial types.",
      spatial$type), call. = FALSE)
  }
  if (length(spatial$shared) >= 2 && isTRUE(spatial$shared[2])) {
    stop("SPDE on the detection process is not yet plumbed in .tobs_laplace; use shared = c(TRUE, FALSE).",
         call. = FALSE)
  }
  if (!isTRUE(spatial$shared[1])) {
    stop("SPDE must be attached to the occupancy/state submodel (shared[1] = TRUE).",
         call. = FALSE)
  }
  if (model_type == "community") {
    stop("Community models with SPDE are not yet plumbed in .tobs_laplace.",
         call. = FALSE)
  }
  if (model_type == "dynamic") {
    stop("Dynamic models with SPDE are not yet plumbed in .tobs_laplace; coming after single-season.",
         call. = FALSE)
  }
  invisible()
}

# Linear-predictor offset induced by the SPDE mesh field at the current fit.
# After tulpa_laplace returns mode = c(beta, u_mesh), the spatial contribution
# to eta at the observed locations is A %*% u_mesh.
.spatial_eta_offset <- function(spatial, fits_sub, p_fixed) {
  if (is.null(spatial) || is.null(fits_sub) || is.null(fits_sub$mode)) {
    return(rep(0, 0))
  }
  if (!identical(spatial$type, "spde")) return(rep(0, 0))
  mode_vec <- fits_sub$mode
  if (length(mode_vec) <= p_fixed) return(rep(0, 0))
  u <- mode_vec[(p_fixed + 1L):length(mode_vec)]
  as.numeric(spatial$tulpa_spec$A %*% u)
}

# Attach the tulpa-side spatial spec to an M-step block. The block's `spatial`
# field is forwarded as-is by tulpa_em_laplace -> tulpa_laplace.
.attach_spatial_spde <- function(block, spatial) {
  if (is.null(spatial) || !identical(spatial$type, "spde")) return(block)
  block$spatial <- spatial$tulpa_spec
  block
}

extract_beta <- function(sub, p) {
  if (is.null(sub)) return(rep(0, p))
  if (!is.null(sub$beta)) return(sub$beta)
  if (!is.null(sub$mean)) return(sub$mean)
  if (!is.null(sub$mode)) return(sub$mode[seq_len(p)])
  rep(0, p)
}

# SE for the fixed-effect block of a tulpa_laplace() fit. Reads the
# negative-log-posterior Hessian (`H_beta`, the precision matrix), inverts
# it, and returns sqrt(diag(.)) restricted to the first `p` fixed effects.
# When the inner fit had a spatial mesh field attached (`spde` / `gp`),
# tulpa_laplace skips H_beta — return NA so callers can flag the
# uncertainty as unavailable instead of carrying a placeholder.
.se_from_laplace_fit <- function(fi, p) {
  if (!is.null(fi$se)) {
    se <- as.numeric(fi$se)
    if (length(se) >= p) return(se[seq_len(p)])
  }
  H <- fi$H_beta
  if (is.null(H)) return(rep(NA_real_, p))
  cov <- tryCatch(solve(H), error = function(e) NULL)
  if (is.null(cov)) return(rep(NA_real_, p))
  d <- sqrt(pmax(diag(cov), 0))
  if (length(d) >= p) return(d[seq_len(p)])
  c(d, rep(NA_real_, p - length(d)))
}

clamp_w <- function(w) pmin(pmax(w, 0.001), 0.999)

occ_weights <- function(psi, p, N, n_valid, n_det, any_det) {
  weights <- numeric(N)
  for (i in seq_len(N)) {
    if (any_det[i]) { weights[i] <- 1 }
    else if (n_valid[i] == 0) { weights[i] <- psi[i] }
    else {
      prod_1mp <- (1 - p[i])^n_valid[i]
      num <- psi[i] * prod_1mp
      weights[i] <- num / (num + (1 - psi[i]))
    }
  }
  weights
}

glm_init <- function(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det) {
  tryCatch({
    occ_glm <- glm(any_det ~ X_occ[, -1, drop = FALSE] - 1 + X_occ[, 1], family = binomial)
    det_glm <- glm(cbind(n_det[keep], n_valid[keep] - n_det[keep]) ~
                      X_det[keep, -1, drop = FALSE] - 1 + X_det[keep, 1], family = binomial)
    list(occ = list(beta = unname(coef(occ_glm)), se = rep(1, p_occ)),
         det = list(beta = unname(coef(det_glm)), se = rep(1, p_det)))
  }, error = function(e) {
    list(occ = list(beta = rep(0, p_occ), se = rep(1, p_occ)),
         det = list(beta = rep(0, p_det), se = rep(1, p_det)))
  })
}

# Build tobs_fit from EM result
build_laplace_fit <- function(em_result, model, spatial, p_per_submodel) {
  pi_list <- model$process_info

  # Collect betas from correction (if available) or EM fits
  means <- numeric()
  sds <- numeric()
  nms <- character()

  for (k in seq_along(pi_list)) {
    pi <- pi_list[[k]]
    sub_name <- names(p_per_submodel)[k]
    if (is.null(sub_name)) sub_name <- names(p_per_submodel)[min(k, length(p_per_submodel))]

    if (is.list(em_result$pooled) && !is.null(em_result$pooled[[sub_name]])) {
      # MI/Gibbs correction pool from rubins_pool().
      cr <- em_result$pooled[[sub_name]]
      means <- c(means, cr$mean)
      sds <- c(sds, cr$se)
    } else if (!is.null(em_result$fits[[sub_name]])) {
      fi <- em_result$fits[[sub_name]]
      beta <- extract_beta(fi, pi$p)
      means <- c(means, beta)
      sds <- c(sds, .se_from_laplace_fit(fi, pi$p))
    } else {
      means <- c(means, rep(0, pi$p))
      sds <- c(sds, rep(NA_real_, pi$p))
    }
    nms <- c(nms, paste0(pi$name, "_", pi$coef_names))
  }

  names(means) <- nms
  names(sds)   <- nms
  n_params <- length(means)

  # Pseudo-draws
  n_pseudo <- 1000L
  draws <- matrix(NA_real_, n_pseudo, n_params)
  for (j in seq_len(n_params)) draws[, j] <- rnorm(n_pseudo, means[j], max(sds[j], 1e-4))
  colnames(draws) <- nms

  intercepts <- compute_intercepts(model, means)

  # When SPDE is attached to the occ submodel, the M-step mode is
  # c(beta_occ, u_mesh). Extract u_mesh so callers can inspect or project
  # the latent field to observation locations via A %*% u_mesh.
  spatial_field <- NULL
  if (!is.null(spatial) && identical(spatial$type, "spde") &&
      !is.null(em_result$fits$occ$mode)) {
    p_occ <- pi_list[[1]]$p
    mode_vec <- em_result$fits$occ$mode
    if (length(mode_vec) > p_occ) {
      spatial_field <- mode_vec[(p_occ + 1L):length(mode_vec)]
    }
  }

  structure(list(
    draws = draws, means = means, sds = sds,
    n_samples = n_pseudo, n_params = n_params,
    log_prob = rep(NA_real_, n_pseudo),
    accept_prob = rep(1, n_pseudo),
    divergent = rep(0L, n_pseudo),
    treedepth = rep(0L, n_pseudo),
    epsilon = NA_real_,
    col_names = nms, param_names = nms,
    intercepts = intercepts,
    model = model, spatial = spatial,
    spatial_field = spatial_field,
    process_info = model$process_info,
    method = "laplace",
    convergence = em_result$convergence,
    correction = em_result$correction
  ), class = c("tobs_fit", "tulpa_fit"))
}
