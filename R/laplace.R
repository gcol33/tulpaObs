# ============================================================================
# Occupancy-specific callbacks for tulpa's EM+Laplace engine
# ============================================================================

#' Fit a tobs model with Laplace approximation (internal)
#'
#' Uses tulpa's generic EM+Laplace engine ([tulpa::tulpa_em_laplace()]) with
#' occupancy-specific E-step and M-step encoding, augmented with a
#' weakly-informative quadratic prior on the fixed-effect coefficients (see
#' [occu_priors()]). The prior is attached to each M-step block as a per-block
#' `beta_prior` (see `.attach_priors_to_blocks()` in `R/occu_priors.R`), which
#' tulpa threads through every phase -- the EM iterations and the MI / Gibbs
#' correction refits alike (gcol33/tulpa#27). Supports all built-in model
#' types. Called from `.tobs_fit_model()`; not user-facing.
#'
#' @param model A `tobs_model` from `.tobs_build_model()`.
#' @param spatial Optional `tobs_spatial` spec (NULL for non-spatial).
#' @param priors Optional prior spec from [occu_priors()]. `NULL` -> use
#'   the package defaults. Pass `FALSE` (or `"none"`) to disable the prior
#'   and recover the historical unpenalised MAP behavior — the penalised
#'   objective is `Q(beta) = -log L(beta) + sum (beta_j - mu_j)^2 / (2 sd_j^2)`,
#'   so `sd_j = Inf` yields a zero penalty term.
#' @param sigma_beta Reserved for future use (NUTS-side beta prior); ignored
#'   by the EM-Laplace path.
#' @param max_iter,tol,damping EM controls.
#' @param correction Post-EM correction (`"none"`, `"mi"`, `"gibbs"`). MI /
#'   Gibbs run tulpa's post-EM Rubin-pooled correction; the fixed-effect prior
#'   (when active) threads into the correction refits, so the corrected fit is
#'   penalised the same way as the EM point estimate (gcol33/tulpa#27).
#' @param n_imputations Number of MI draws when `correction = "mi"`.
#' @param verbose Print per-iteration progress.
#' @keywords internal
.tobs_laplace <- function(model, spatial = NULL, re = NULL,
                          priors = NULL,
                          sigma_beta = 10,
                          max_iter = 50L, tol = 1e-4, damping = 0.3,
                          correction = c("auto", "mi", "gibbs", "none"),
                          n_imputations = 20L, n_gibbs = 10L, seed = NULL,
                          approx = c("gaussian_laplace", "simplified_laplace"),
                          verbose = TRUE) {
  correction <- match.arg(correction)
  approx <- match.arg(approx)
  if (!inherits(model, "tobs_model")) stop("model must be a tobs_model object")

  .validate_spatial_laplace(spatial, model$model_type)

  # Formula random effects on the deterministic path. Supported forms (iid
  # intercept, uncorrelated slopes, and correlated slopes on the occupancy
  # predictor of a single-season model) are fit via the variance-component EM
  # in R/em_laplace_re.R; everything else errors with a pointer to NUTS rather
  # than being silently dropped (gcol33/tulpaObs#11).
  if (!is.null(re)) {
    .validate_re_laplace(re, model, spatial, approx)
    em_result <- .tobs_em_laplace_re(model, re, priors = priors,
                                     max_iter = max_iter, tol = tol,
                                     damping = damping, verbose = verbose)
    re_block <- .tobs_re_param_block(em_result$re_post)
    fit <- build_laplace_fit(em_result, model, spatial,
                             c(occ = ncol(model$X_processes[[1]]),
                               det = ncol(model$X_processes[[2]])),
                             prior_spec = NULL, approx = "gaussian_laplace",
                             re_block = re_block)
    fit$re <- if (inherits(re, "tobs_re")) list(re) else re
    return(fit)
  }

  callbacks <- switch(model$model_type,
    single     = build_single_callbacks(model, spatial),
    dynamic    = build_dynamic_callbacks(model, spatial),
    community  = build_community_callbacks(model, spatial),
    integrated = build_integrated_callbacks(model, spatial),
    jsdm       = build_jsdm_callbacks(model, spatial),
    stop(sprintf("Laplace not supported for model_type '%s'", model$model_type))
  )

  prior_spec <- .resolve_occu_priors(priors)

  # Single engine for every Laplace fit: tulpa's generic EM+Laplace. The
  # fixed-effect prior is attached per M-step block as a `beta_prior` (see
  # .attach_priors_to_blocks); tulpa's block fitter applies it in every phase,
  # so a prior-aware MI/Gibbs correction comes for free (gcol33/tulpa#27).
  # Spatial fits are left unpenalised here -- the SPDE/NNGP solver carries its
  # own fixed-effect prior and tulpa_laplace() rejects `beta_prior` on the
  # spatial path (tulpaObs#5) -- so the prior is attached only when there is no
  # spatial term.
  m_step_encode <- if (is.null(spatial)) {
    function(weights, ...) {
      .attach_priors_to_blocks(callbacks$m_step_encode(weights, ...),
                               model, prior_spec)
    }
  } else {
    callbacks$m_step_encode
  }

  # MI / Gibbs draw hard z with R's RNG; seed it so the corrected fit
  # reproduces.
  if (correction %in% c("mi", "gibbs") && !is.null(seed)) {
    set.seed(as.integer(seed))
  }
  em_result <- tulpa::tulpa_em_laplace(
    e_step        = callbacks$e_step,
    m_step_encode = m_step_encode,
    draw_z        = callbacks$z_draw,
    max_iter      = max_iter,
    tol           = tol,
    damping       = damping,
    correction    = correction,
    n_imputations = n_imputations,
    n_gibbs       = n_gibbs,
    verbose       = verbose
  )

  # tulpa_em_laplace returns flat convergence fields; synthesize the nested
  # `convergence` list that build_laplace_fit() / summary() expect.
  if (is.null(em_result$convergence)) {
    em_result$convergence <- list(converged = em_result$converged,
                                  n_iter = em_result$n_iter,
                                  history = em_result$history)
  }

  fit <- build_laplace_fit(em_result, model, spatial, callbacks$p_per_submodel,
                           prior_spec = prior_spec,
                           approx = approx)
  fit$priors <- prior_spec
  # Record the seed used for a stochastic correction so the run reproduces.
  if (correction %in% c("mi", "gibbs") && !is.null(seed)) {
    fit$seed <- as.integer(seed)
  }
  fit
}

# ============================================================================
# Single-season callbacks
# ============================================================================
build_single_callbacks <- function(model, spatial = NULL) {
  y <- model$y
  X_occ <- model$X_processes[[1]]
  X_det <- model$X_processes[[2]]
  X_det_visit <- model$X_det_visit  # NULL when no visit-level covariates
  max_visits <- ncol(y)
  n_sites <- model$n_sites
  p_occ <- ncol(X_occ)
  p_det <- ncol(X_det)
  p_det_visit <- if (is.null(X_det_visit)) 0L else ncol(X_det_visit)
  p_det_total <- p_det + p_det_visit

  n_valid <- integer(n_sites)
  n_det <- integer(n_sites)
  any_det <- logical(n_sites)
  valid_mat <- matrix(FALSE, n_sites, max_visits)
  for (i in seq_len(n_sites)) {
    v <- y[i, ] >= 0
    valid_mat[i, ] <- v
    n_valid[i] <- sum(v)
    n_det[i] <- sum(y[i, v] == 1)
    any_det[i] <- n_det[i] > 0
  }
  keep <- n_valid > 0

  # Per-visit indexing for the X_det_visit path. Row r of X_det_visit
  # corresponds to (site = (r-1) %/% max_visits + 1, visit = (r-1) %% max_visits + 1).
  if (p_det_visit > 0L) {
    site_idx_all <- rep(seq_len(n_sites), each = max_visits)
    visit_idx_all <- rep(seq_len(max_visits), times = n_sites)
  }

  e_step <- function(fits, ...) {
    beta_occ <- extract_beta(fits$occ, p_occ)
    eta_occ <- as.vector(X_occ %*% beta_occ)
    sp_off <- .spatial_eta_offset(spatial, fits$occ, p_occ)
    if (length(sp_off) == n_sites) eta_occ <- eta_occ + sp_off
    psi <- plogis(eta_occ)

    if (p_det_visit == 0L) {
      beta_det <- extract_beta(fits$det, p_det)
      p_site <- plogis(as.vector(X_det %*% beta_det))
      return(list(weights = occ_weights(psi, p_site, n_sites,
                                        n_valid, n_det, any_det)))
    }

    # Visit-level path: logit(p_ij) = X_det[i,] beta_site + X_det_visit[(i-1)*J + j,] beta_visit
    beta_det <- extract_beta(fits$det, p_det_total)
    eta_site <- as.vector(X_det %*% beta_det[seq_len(p_det)])
    eta_visit_long <- as.vector(X_det_visit %*%
                                  beta_det[(p_det + 1L):p_det_total])
    # eta_visit_long is in site-major order: reshape so [i, j] = visit (i, j)
    eta_visit_mat <- matrix(eta_visit_long, n_sites, max_visits, byrow = TRUE)
    logit_p_ij <- matrix(eta_site, n_sites, max_visits) + eta_visit_mat
    logit_p_ij <- pmin(pmax(logit_p_ij, -30), 30)
    # log(1 - plogis(eta)) = -log1pexp(eta) computed stably as -pmax(eta,0) - log1p(exp(-|eta|))
    log_1mp <- -(pmax(logit_p_ij, 0) + log1p(exp(-abs(logit_p_ij))))
    log_1mp[!valid_mat] <- 0
    log_prod_1mp <- rowSums(log_1mp)
    weights <- numeric(n_sites)
    for (i in seq_len(n_sites)) {
      if (any_det[i]) {
        weights[i] <- 1
      } else if (n_valid[i] == 0L) {
        weights[i] <- psi[i]
      } else {
        num <- psi[i] * exp(log_prod_1mp[i])
        weights[i] <- num / (num + (1 - psi[i]))
      }
    }
    list(weights = weights)
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

    if (p_det_visit == 0L) {
      keep_det <- keep & (w_det > 1e-6)
      det_block <- list(y = n_det[keep_det], n_trials = n_valid[keep_det],
                        X = X_det[keep_det, , drop = FALSE],
                        weights = w_det[keep_det], family = "binomial")
    } else {
      # Per-visit detection block: one Bernoulli row per (site, valid visit)
      # whose weight is the site's posterior occupancy w_i. Combined design
      # matrix stacks site-level X_det (replicated across visits) and
      # visit-level X_det_visit (already in site-major order).
      keep_visit <- valid_mat & (w_det >= 1e-6)
      site_kept <- site_idx_all[as.vector(t(keep_visit))]
      visit_kept <- visit_idx_all[as.vector(t(keep_visit))]
      n_kept <- length(site_kept)
      if (n_kept > 0L) {
        visit_row_idx <- (site_kept - 1L) * max_visits + visit_kept
        X_combined <- cbind(
          X_det[site_kept, , drop = FALSE],
          X_det_visit[visit_row_idx, , drop = FALSE]
        )
        y_kept <- y[cbind(site_kept, visit_kept)]
        det_block <- list(y = as.integer(y_kept),
                          n_trials = rep(1L, n_kept),
                          X = X_combined,
                          weights = w_det[site_kept],
                          family = "binomial")
      } else {
        det_block <- list(y = integer(0),
                          n_trials = integer(0),
                          X = matrix(0, 0, p_det_total),
                          weights = numeric(0),
                          family = "binomial")
      }
    }
    list(occ = occ_block, det = det_block)
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
    if (p_det_visit == 0L) {
      det_block <- if (length(det_keep) > 0)
        list(y = n_det[det_keep], n_trials = n_valid[det_keep],
             X = X_det[det_keep, , drop = FALSE], family = "binomial")
      else NULL
    } else {
      keep_visit <- valid_mat & matrix(z == 1L, n_sites, max_visits)
      site_kept <- site_idx_all[as.vector(t(keep_visit))]
      visit_kept <- visit_idx_all[as.vector(t(keep_visit))]
      n_kept <- length(site_kept)
      det_block <- if (n_kept > 0L) {
        visit_row_idx <- (site_kept - 1L) * max_visits + visit_kept
        X_combined <- cbind(
          X_det[site_kept, , drop = FALSE],
          X_det_visit[visit_row_idx, , drop = FALSE]
        )
        y_kept <- y[cbind(site_kept, visit_kept)]
        list(y = as.integer(y_kept), n_trials = rep(1L, n_kept),
             X = X_combined, family = "binomial")
      } else NULL
    }
    list(occ = occ_block, det = det_block)
  }

  init <- glm_init(X_occ, X_det, any_det, n_det, n_valid, keep, p_occ, p_det)
  if (p_det_visit > 0L) {
    # Pad det init with zeros for visit-level cols so warm-start shapes
    # match the combined design when the penalized driver pulls beta_init
    # from the previous EM iteration.
    init$det$beta <- c(init$det$beta, rep(0, p_det_visit))
    init$det$se   <- c(init$det$se,   rep(1, p_det_visit))
  }

  list(e_step = e_step, m_step_encode = m_step_encode, z_draw = z_draw,
       hard_encode = hard_encode, init = init,
       p_per_submodel = c(occ = p_occ, det = p_det_total))
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
    stop("spatial must be a tobs_spatial term (from a spde() formula term)",
         call. = FALSE)
  }
  if (!identical(spatial$type, "spde")) {
    stop(sprintf(
      ".tobs_laplace currently supports spatial$type == 'spde' only (got '%s'). Use method = 'nuts' for other spatial types.",
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

# Gate the deterministic random-effect path. The variance-component EM in
# R/em_laplace_re.R fits iid intercept, uncorrelated slopes, and correlated
# slopes (a full RE covariance) on the occupancy predictor of a single-season
# model. Forms it cannot fit -- non-single families, RE + spatial, RE +
# visit-level detection, RE on / shared with detection -- error here with a
# pointer to `method = "nuts"` (which fits every RE form) rather than being
# silently dropped (gcol33/tulpaObs#11). The deterministic Laplace variance
# estimate carries the usual small-cluster (PQL) bias; NUTS is the calibrated
# route when that matters.
.validate_re_laplace <- function(re, model, spatial, approx) {
  re_list <- if (inherits(re, "tobs_re")) list(re) else re

  if (!identical(model$model_type, "single")) {
    stop(sprintf(
      "Random effects under method = 'laplace' are wired for single-season occupancy only (got model_type = '%s'). Use method = 'nuts' for random effects on this family.",
      model$model_type), call. = FALSE)
  }
  if (!is.null(spatial)) {
    stop("A random effect combined with a spatial term is not supported on the Laplace path. Use method = 'nuts'.",
         call. = FALSE)
  }
  if (!is.null(model$X_det_visit)) {
    stop("Random effects with visit-level detection covariates are not supported on the Laplace path. Use method = 'nuts'.",
         call. = FALSE)
  }
  for (r in re_list) {
    if (length(r$shared) >= 2L && isTRUE(r$shared[2]) && !isTRUE(r$shared[1])) {
      stop("A random effect on the detection predictor is not supported on the Laplace path. Use method = 'nuts'.",
           call. = FALSE)
    }
    if (length(r$shared) >= 2L && isTRUE(r$shared[2])) {
      stop("A random effect shared across occupancy and detection is not supported on the Laplace path. Use method = 'nuts'.",
           call. = FALSE)
    }
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

# Louis-corrected observed Fisher info for the occupancy fixed-effect block
# of a single-season occu fit (tulpaObs#7).
#
# Why this is needed. The inner M-step encodes the soft-imputed P(z_i = 1 | y_i)
# as a pseudo-binomial likelihood with n_trials = M (M = 1000 non-spatial,
# M = 4 spatial). The resulting inner Hessian is
#
#   H_inner = M * X' diag(psi (1 - psi)) X + P_prior
#
# i.e. the complete-data Fisher info inflated by the M trick, plus the prior
# precision. This is the wrong object for SE reporting on two counts: the M
# factor is an artefact of the M-step encoding (not data information), and
# the complete-data info ignores the missing-z variance.
#
# Louis identity for the occupancy score s_i = x_i (z_i - psi_i) gives the
# observed Fisher info at the EM stationary point:
#
#   I_obs(beta_psi) = E[-d2 log f / dbeta2 | y] - Var(s_complete | y)
#                   = X' diag(psi (1 - psi)) X - X' diag(w (1 - w)) X
#                   = X' diag(psi (1 - psi) - w (1 - w)) X
#
# where w_i = P(z_i = 1 | y_i, theta_hat) is the converged E-step weight. The
# per-site `psi(1-psi) - w(1-w)` term can be negative (the marginal log-lik
# can be locally convex at a single site), but the aggregate X' D X is PSD at
# the MLE because it equals minus the marginal log-lik Hessian at its max.
.louis_info_psi_single <- function(X_occ, beta_psi, weights,
                                   spatial = NULL, spatial_fit = NULL,
                                   prior_spec = NULL,
                                   coef_names = NULL) {
  p_psi <- length(beta_psi)
  if (p_psi == 0L) return(NULL)
  if (is.null(X_occ) || nrow(X_occ) == 0L) return(NULL)
  if (is.null(weights) || length(weights) != nrow(X_occ)) return(NULL)

  eta <- as.numeric(X_occ %*% beta_psi)
  sp_off <- .spatial_eta_offset(spatial, spatial_fit, p_psi)
  if (length(sp_off) == nrow(X_occ)) eta <- eta + sp_off
  eta <- pmin(pmax(eta, -30), 30)
  psi <- plogis(eta)

  d <- psi * (1 - psi) - weights * (1 - weights)
  I_obs <- as.matrix(crossprod(X_occ, d * X_occ))

  if (!is.null(prior_spec)) {
    if (is.null(coef_names)) coef_names <- colnames(X_occ) %||% paste0("x", seq_len(p_psi))
    pr <- .prior_for_submodel(prior_spec, "psi", coef_names)
    if (!is.null(pr)) {
      pen_prec <- ifelse(is.finite(pr$sd), 1 / (pr$sd^2), 0)
      diag(I_obs) <- diag(I_obs) + pen_prec[seq_len(p_psi)]
    }
  }
  I_obs
}

# SE vector from an observed-info matrix; returns NA of length p on failure.
.se_from_info <- function(I, p) {
  if (is.null(I)) return(rep(NA_real_, p))
  cov <- tryCatch(solve(I), error = function(e) NULL)
  if (is.null(cov)) return(rep(NA_real_, p))
  d <- sqrt(pmax(diag(cov), 0))
  if (length(d) >= p) d[seq_len(p)] else c(d, rep(NA_real_, p - length(d)))
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
build_laplace_fit <- function(em_result, model, spatial, p_per_submodel,
                              prior_spec = NULL,
                              approx = "gaussian_laplace",
                              re_block = NULL) {
  pi_list <- model$process_info

  # Collect betas from correction (if available) or EM fits
  means <- numeric()
  sds <- numeric()
  nms <- character()
  louis_psi_se <- NULL

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

      # Louis-corrected SE on the psi block of a single-season fit. The inner
      # M-step Hessian is M * I_complete + P_prior (pseudo-binomial trick);
      # I_obs = X' diag(psi(1-psi) - w(1-w)) X + P_prior is the right object
      # for SEs. See `.louis_info_psi_single` and tulpaObs#7.
      # Louis-corrected occupancy SE applies to the fixed-effect-only fit. When
      # random effects are present the occupancy block's fixed-effect SE comes
      # from the GLMM marginal precision (`H_beta`, Schur over the RE block)
      # that tulpa_laplace returns, so skip Louis on the RE path.
      use_louis <- identical(model$model_type, "single") &&
                   identical(sub_name, "occ") &&
                   !is.null(em_result$weights) &&
                   is.null(re_block)
      if (!is.null(re_block) && identical(sub_name, "occ")) {
        # Occupancy fixed-effect SE on the RE path: natural-scale observed info
        # marginalised over the random-effect block (the M-step H_beta is
        # M-inflated). Computed in .tobs_re_occ_fixed_se().
        sds <- c(sds, re_block$occ_se)
      } else if (use_louis) {
        I_obs <- .louis_info_psi_single(
          X_occ       = model$X_processes[[1]],
          beta_psi    = beta,
          weights     = em_result$weights,
          spatial     = spatial,
          spatial_fit = fi,
          prior_spec  = prior_spec,
          coef_names  = pi$coef_names
        )
        louis_psi_se <- .se_from_info(I_obs, pi$p)
        sds <- c(sds, louis_psi_se)
      } else {
        sds <- c(sds, .se_from_laplace_fit(fi, pi$p))
      }
    } else {
      means <- c(means, rep(0, pi$p))
      sds <- c(sds, rep(NA_real_, pi$p))
    }
    nms <- c(nms, paste0(pi$name, "_", pi$coef_names))
  }

  # Append visit-level detection coefficients when X_det_visit is present.
  # The detection M-step block has X of width p_det + p_det_visit; the main
  # loop above extracts only the first p_det elements (the site-level
  # detection coefs). Pull the visit-level tail and label as `p_visit_<name>`
  # so the public output matches the NUTS engine's column layout.
  if (!is.null(model$det_visit_names) && length(model$det_visit_names) > 0L) {
    p_det_visit <- length(model$det_visit_names)
    pi_p <- pi_list[[2]]  # detection process metadata
    p_det <- pi_p$p
    p_det_total <- as.integer(p_per_submodel[["det"]] %||% (p_det + p_det_visit))
    visit_idx <- (p_det + 1L):p_det_total
    visit_nms <- paste0("p_visit_", model$det_visit_names)

    if (is.list(em_result$pooled) && !is.null(em_result$pooled[["det"]])) {
      cr <- em_result$pooled[["det"]]
      visit_means <- cr$mean[visit_idx]
      visit_sds   <- cr$se[visit_idx]
    } else if (!is.null(em_result$fits[["det"]])) {
      fi_det <- em_result$fits[["det"]]
      beta_full <- extract_beta(fi_det, p_det_total)
      se_full <- .se_from_laplace_fit(fi_det, p_det_total)
      visit_means <- beta_full[visit_idx]
      visit_sds   <- se_full[visit_idx]
    } else {
      visit_means <- rep(0, p_det_visit)
      visit_sds   <- rep(NA_real_, p_det_visit)
    }

    means <- c(means, visit_means)
    sds   <- c(sds, visit_sds)
    nms   <- c(nms, visit_nms)
  }

  # Append the deterministic random-effect block (sigma hyperparameters +
  # per-group BLUPs) so the public output matches the NUTS column layout and
  # ranef() / summary() can name them (gcol33/tulpaObs#11).
  if (!is.null(re_block)) {
    means <- c(means, re_block$means)
    sds   <- c(sds, re_block$sds)
    nms   <- c(nms, re_block$names)
  }

  names(means) <- nms
  names(sds)   <- nms
  n_params <- length(means)

  # Pseudo-draws
  n_pseudo <- 1000L
  draws <- matrix(NA_real_, n_pseudo, n_params)
  for (j in seq_len(n_params)) {
    # Hyperparameters without an analytic SE (e.g. RE sigma on the
    # variance-component path) carry NA sd -> draw a near-constant column at
    # the point estimate rather than NAs.
    sd_j <- if (is.finite(sds[j])) sds[j] else 0
    draws[, j] <- rnorm(n_pseudo, means[j], max(sd_j, 1e-4))
  }
  colnames(draws) <- nms

  # Simplified-Laplace skewness correction
  # Computes gamma_j at the original observation likelihood (NOT the M-step
  # pseudo-binomial encoding — see dev_notes/simplified_laplace_derivation.md
  # §3 and dev_notes/upstream_tulpa_sla_spec.md §3 for why).
  sla_gamma <- NULL
  sla_status <- "off"
  # Simplified-Laplace skewness correction is not wired for the random-effect
  # path (the gamma derivation assumes a fixed-effect-only M-step).
  if (identical(approx, "simplified_laplace") && is.null(re_block)) {
    sla_res <- switch(model$model_type,
      single     = .sla_compute_occu_single(model, em_result,
                                            spatial = spatial,
                                            prior_spec = prior_spec),
      dynamic    = .sla_compute_dyn_occu(model, em_result,
                                         spatial = spatial,
                                         prior_spec = prior_spec),
      integrated = .sla_compute_int_occu(model, em_result,
                                         spatial = spatial,
                                         prior_spec = prior_spec),
      list(gamma = NULL, valid = FALSE,
           reason = sprintf("simplified Laplace not yet supported for model_type '%s'",
                            model$model_type))
    )
    if (isTRUE(sla_res$valid)) {
      sla_gamma  <- sla_res$gamma
      sla_status <- "simplified_laplace"
      # Align gamma names with the joint parameter ordering used in `means`
      sla_gamma <- sla_gamma[intersect(names(means), names(sla_gamma))]
      gamma_full <- setNames(rep(0, n_params), nms)
      gamma_full[names(sla_gamma)] <- sla_gamma
      sla_gamma <- gamma_full
      draws <- .sla_replace_draws(draws, means, sds, sla_gamma)
    } else {
      sla_status <- paste0("fallback_gaussian (", sla_res$reason, ")")
    }
  }

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
    skew = sla_gamma, sla_status = sla_status,
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
    re_effects = re_block$re_effects,
    convergence = em_result$convergence,
    correction = em_result$correction
  ), class = c("tobs_fit", "tulpa_fit"))
}
