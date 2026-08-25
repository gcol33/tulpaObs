# =============================================================================
# ms_occu_cover.R - community / multispecies joint occupancy-detection + cover
#
# The community version of occu_cover(): a per-species joint occupancy-cover
# model with Gaussian community hyperpriors on the per-species coefficients of
# all three arms.
#
#   z_{s,i}        ~ Bernoulli(psi_{s,i})                       (latent presence)
#   y_{s,i,j}|z=1  ~ Bernoulli(p_{s,i,j})                       (detection)
#   c_{s,i,j}|y=1  ~ f_pos(eta_pos_{s,i,j}, disp)               (positive cover)
#   logit psi_{s,i}     = X_occ_i . (mu_occ + b_occ_s)
#   logit p_{s,i,j}     = X_p_{ij} . (mu_p   + b_p_s)
#   g(cover)_{s,i,j}    = X_pos_{ij} . (mu_pos + b_pos_s)
#   b_occ_s ~ N(0, Sigma_occ), b_p_s ~ N(0, Sigma_p),
#   b_pos_s ~ N(0, Sigma_pos)                                   (community RE)
#
# The latent presence z marginalises out per species-cell in closed form (the
# same two-state mixture as occu_cover(), reused via .occu_cover_site_ll); the
# per-species coefficient deviations b_s = (b_occ_s, b_p_s, b_pos_s) are the
# random effects. The positive-arm dispersion is a shared community parameter
# (one log_disp), as in the single-species occu_cover().
#
# Fit: the shared community Laplace-EM engine (`.tobs_community_em()`,
# R/community_em.R), with the occ/p/pos arms as `arm_idx` and the shared
# log-dispersion as its one `global`. Per EM iteration the engine's arrowhead
# joint Newton finds the mode of (mu, log_disp, {b_s}) -- the per-species RE
# blocks fold out by a Schur complement -- from this family's per-species
# gradients (`.occu_cover_eta_grad`, via `.occu_cover_sp_ll`/`_sp_grad`); the
# observed-information block is the engine's own finite-difference-of-gradient
# fallback (no analytic `sp_info` here). The M-step is the closed-form
# community-covariance update Sigma_arm = mean_s[b_s b_s' + Cov(b_s | y)].
# Community-mean SEs are the marginal observed information (the Schur
# complement of the b-block, Louis 1982), read at the natural scale.
#
# This is the family wiring + the in-tree non-spatial fitter; it reuses the
# single-species occu_cover() linear-predictor builder and per-cell marginal as
# the per-species kernel (single source of truth). The community SPATIAL model --
# a shared latent field coupled across the occupancy and cover arms with
# per-species RE on all three arms -- is the reduced-rank spatial-factor fit
# reached by a front-door icar()/car_proper()/bym2() term
# (R/ms_occu_cover_spatial.R): per-species loadings on the shared field via
# Laplace-EM / NUTS, with per-species community covariances on top. The free
# per-species loading form generalises a single common field amplitude, so it
# subsumes the common-amplitude coupling of occu_cover()'s joint-coupled engine.
# Routing the community model through tulpa's joint cell-coupling engine instead
# (per-arm RE blocks integrated on the outer grid) is not viable: the joint engine
# integrates every variance component on its grid, so per-arm community RE
# variances plus the field hyperparameters exceed the grid cap -- the closed-form
# covariance M-step of the Laplace-EM is the scaling route for community variance
# components. The non-spatial dispatcher below rejects a residual structured term
# with a pointer rather than silently dropping it.
# =============================================================================


# ---------------------------------------------------------------------------
# Data binder
# ---------------------------------------------------------------------------

# Bind a community joint occupancy-cover model. `y` / `y_pos` are 3D arrays
# [n_sites x max_visits x n_species] or named lists of n_sites x max_visits
# matrices, one per species. The occupancy design X_occ is cell-level; the
# detection / cover designs are the shared site-level design plus any visit-
# level covariates (community covariates, identical across species). Detection
# is 0/1/NA; cover is meaningful only where detection == 1.
.tobs_build_ms_occu_cover <- function(occ_formula, det_formula, pos_formula,
                                      data, y, y_pos, positive, species,
                                      det_visit_formula = NULL,
                                      det_visit_data    = NULL,
                                      pos_visit_formula = NULL,
                                      pos_visit_data    = NULL) {
  to_array <- function(z, label) {
    if (is.list(z) && !is.array(z)) {
      n_sp <- length(z)
      arr <- array(NA_real_, dim = c(nrow(z[[1L]]), ncol(z[[1L]]), n_sp))
      for (s in seq_len(n_sp)) arr[, , s] <- as.matrix(z[[s]])
      attr(arr, "names_from") <- names(z)
      return(arr)
    }
    if (length(dim(z)) != 3L) {
      stop(sprintf("%s must be a 3D array [n_sites x max_visits x n_species] ",
                   label), "or a list of matrices.", call. = FALSE)
    }
    z
  }
  y     <- to_array(y,     "y")
  y_pos <- to_array(y_pos, "y_pos")
  if (!all(dim(y) == dim(y_pos))) {
    stop("y and y_pos must have identical dimensions ",
         "[n_sites x max_visits x n_species].", call. = FALSE)
  }

  n_sites    <- dim(y)[1L]
  max_visits <- dim(y)[2L]
  n_species  <- dim(y)[3L]

  species_names <- if (is.character(species)) species
                   else if (!is.null(attr(y, "names_from"))) attr(y, "names_from")
                   else paste0("sp", seq_len(n_species))
  if (length(species_names) != n_species) {
    stop(sprintf("species has %d names but y has %d species",
                 length(species_names), n_species), call. = FALSE)
  }
  .tobs_check_site_count(n_sites, nrow(data), "sites")

  # Reject structured terms on every arm (spatial sharing + per-species RE is
  # not wired; see file header).
  .occu_cover_reject_structured(occ_formula, "occupancy")
  .occu_cover_reject_structured(det_formula, "detection")
  .occu_cover_reject_structured(pos_formula, "positive cover")

  X_occ      <- stats::model.matrix(occ_formula, data)
  X_det_site <- stats::model.matrix(det_formula, data)
  X_pos_site <- stats::model.matrix(pos_formula, data)
  X_det_visit <- .tobs_build_visit_X(det_visit_formula, det_visit_data,
                                     n_sites, max_visits, arm = "detection")
  X_pos_visit <- .tobs_build_visit_X(pos_visit_formula, pos_visit_data,
                                     n_sites, max_visits, arm = "positive cover")

  det_coef_names <- colnames(X_det_site)
  pos_coef_names <- colnames(X_pos_site)
  if (!is.null(X_det_visit)) det_coef_names <- c(det_coef_names, colnames(X_det_visit))
  if (!is.null(X_pos_visit)) pos_coef_names <- c(pos_coef_names, colnames(X_pos_visit))

  # Validate + clean each species' detection / cover matrices (binary detection;
  # positive cover finite and in-support at detected visits), and stash the
  # masked, zero-filled per-species matrices.
  y_int   <- array(0L,         dim = dim(y))
  y_posn  <- array(0,          dim = dim(y))
  valid   <- array(FALSE,      dim = dim(y))
  for (s in seq_len(n_species)) {
    ys  <- matrix(as.integer(round(y[, , s])), n_sites, max_visits)
    vs  <- !is.na(ys)
    if (any(ys[vs] != 0L & ys[vs] != 1L)) {
      stop(sprintf("species '%s': y must contain only 0, 1, or NA.",
                   species_names[s]), call. = FALSE)
    }
    ys[!vs] <- 0L
    yp   <- matrix(as.numeric(y_pos[, , s]), n_sites, max_visits)
    pmsk <- vs & (ys == 1L)
    if (any(!is.finite(yp[pmsk]))) {
      stop(sprintf("species '%s': y_pos must be finite at every detected visit.",
                   species_names[s]), call. = FALSE)
    }
    if (identical(positive, "beta")) {
      if (any(pmsk & (yp <= 0 | yp >= 1))) {
        stop(sprintf("species '%s': beta cover requires 0 < y_pos < 1 at every ",
                     species_names[s]), "detected visit.", call. = FALSE)
      }
    } else if (identical(positive, "lognormal")) {
      if (any(pmsk & (yp <= 0))) {
        stop(sprintf("species '%s': lognormal cover requires y_pos > 0 at every ",
                     species_names[s]), "detected visit.", call. = FALSE)
      }
    }
    # gaussian: any finite real is in support (checked finite above).
    yp[!pmsk] <- 0
    y_int[, , s]  <- ys
    y_posn[, , s] <- yp
    valid[, , s]  <- vs
  }

  structure(list(
    model_type   = "ms_occu_cover",
    positive     = positive,
    y            = y_int,
    y_pos        = y_posn,
    valid        = valid,
    n_sites      = n_sites,
    max_visits   = max_visits,
    n_species    = n_species,
    species_names = species_names,
    X_occ        = X_occ,
    X_det_site   = X_det_site,
    X_pos_site   = X_pos_site,
    X_det_visit  = X_det_visit,
    X_pos_visit  = X_pos_visit,
    formulas     = list(occ = occ_formula, det = det_formula, pos = pos_formula),
    data         = data,
    process_info = list(
      list(name = "psi", p = ncol(X_occ),
           coef_names = colnames(X_occ), link = "logit"),
      list(name = "p",   p = length(det_coef_names),
           coef_names = det_coef_names, link = "logit"),
      list(name = "pos", p = length(pos_coef_names),
           coef_names = pos_coef_names,
           link = if (positive == "beta") "logit" else "identity")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# Per-species kernel: log-likelihood, eta gradient, coefficient gradient
# ---------------------------------------------------------------------------

# A lightweight single-species occu_cover "model" view sharing the community
# design but carrying species `s`'s detection / cover matrices. The eta / ll /
# grad helpers below read $valid / $y / $y_pos and the $X_* designs, so this
# view drives the single-species occu_cover kernel directly.
.ms_occu_cover_species_view <- function(model, s) {
  model$y     <- model$y[, , s]
  model$y_pos <- model$y_pos[, , s]
  model$valid <- model$valid[, , s]
  model
}

# Gradient of a species' marginal log-likelihood with respect to the three arm
# linear predictors and the shared log-dispersion, given the per-cell occupancy
# probability `psi`, the per-visit detection probability `p_mat`, the per-visit
# cover linear predictor `ep_mat`, and the scalar `log_disp`. The latent
# presence z is integrated out (the same det / nodet branches as
# .occu_cover_site_ll), so this is the exact marginal gradient. Returns the
# per-cell `g_psi` (length n_sites), the per-visit `g_p` / `g_pos`
# [n_sites x max_visits], and the scalar `g_ld`.
.occu_cover_eta_grad <- function(m, psi, p_mat, ep_mat, log_disp) {
  cl <- .tobs_clamp_eta
  valid <- m$valid; y <- m$y; y_pos <- m$y_pos
  n_sites <- m$n_sites; max_visits <- m$max_visits

  log_1mp <- ifelse(valid, log(1 - p_mat), 0)
  prod1mp <- exp(rowSums(log_1mp))
  any_det <- rowSums(y * valid) > 0
  A <- psi * prod1mp                  # occupied-but-undetected mass per cell
  B <- 1 - psi                        # unoccupied mass
  L <- A + B                          # no-detection marginal per cell

  # ---- occupancy predictor ----
  g_psi <- numeric(n_sites)
  g_psi[any_det] <- 1 - psi[any_det]
  nd <- !any_det
  g_psi[nd] <- psi[nd] * (1 - psi[nd]) * (prod1mp[nd] - 1) / L[nd]

  # ---- detection predictor ----
  g_p <- matrix(0, n_sites, max_visits)
  det_rows <- which(any_det)
  if (length(det_rows)) {
    pm <- p_mat[det_rows, , drop = FALSE]
    ym <- y[det_rows, , drop = FALSE]
    vm <- valid[det_rows, , drop = FALSE]
    gd <- ifelse(ym == 1L, 1 - pm, -pm)
    gd[!vm] <- 0
    g_p[det_rows, ] <- gd
  }
  nd_rows <- which(nd)
  if (length(nd_rows)) {
    pm <- p_mat[nd_rows, , drop = FALSE]
    vm <- valid[nd_rows, , drop = FALSE]
    gn <- -(A[nd_rows] / L[nd_rows]) * pm
    gn[!vm] <- 0
    g_p[nd_rows, ] <- gn
  }

  # ---- cover predictor + shared dispersion ----
  g_pos <- matrix(0, n_sites, max_visits)
  g_ld  <- 0
  pos_mask <- valid & (y == 1L)
  if (any(pos_mask)) {
    yp_safe <- y_pos; yp_safe[!pos_mask] <- 0.5   # keep logs finite off the mask
    if (identical(m$positive, "beta")) {
      phi <- exp(log_disp)
      mu  <- stats::plogis(cl(ep_mat))
      a   <- mu * phi; b <- (1 - mu) * phi
      dmu <- phi * (-digamma(a) + digamma(b) + log(yp_safe) - log(1 - yp_safe))
      gpe <- dmu * mu * (1 - mu)
      gpe[!pos_mask] <- 0
      g_pos <- gpe
      dld <- phi * (digamma(phi) - digamma(a) * mu - digamma(b) * (1 - mu) +
                    mu * log(yp_safe) + (1 - mu) * log(1 - yp_safe))
      g_ld <- sum(dld[pos_mask])
    } else if (identical(m$positive, "gaussian")) {
      # Identity-Gaussian arm: mu = eta, residual on the raw response, no log
      # Jacobian.
      sigma <- exp(log_disp)
      r   <- (yp_safe - ep_mat) / sigma
      gpe <- r / sigma                            # (c - eta) / sigma^2
      gpe[!pos_mask] <- 0
      g_pos <- gpe
      g_ld <- sum((r^2 - 1)[pos_mask])
    } else {
      sigma <- exp(log_disp)
      r   <- (log(yp_safe) - ep_mat) / sigma
      gpe <- r / sigma                            # (log c - eta) / sigma^2
      gpe[!pos_mask] <- 0
      g_pos <- gpe
      g_ld <- sum((r^2 - 1)[pos_mask])
    }
  }
  list(g_psi = g_psi, g_p = g_p, g_pos = g_pos, g_ld = g_ld)
}

# Chain the eta-gradient to the packed coefficient gradient
# c(beta_occ, beta_p, beta_pos, log_disp), following the fitter's site-level +
# visit-level packing on the detection and cover arms.
.occu_cover_coef_grad <- function(m, eg) {
  g_occ <- as.numeric(crossprod(m$X_occ, eg$g_psi))

  g_det <- as.numeric(crossprod(m$X_det_site, rowSums(eg$g_p)))
  if (!is.null(m$X_det_visit)) {
    g_det <- c(g_det, as.numeric(crossprod(m$X_det_visit, as.vector(t(eg$g_p)))))
  }
  g_pos <- as.numeric(crossprod(m$X_pos_site, rowSums(eg$g_pos)))
  if (!is.null(m$X_pos_visit)) {
    g_pos <- c(g_pos, as.numeric(crossprod(m$X_pos_visit, as.vector(t(eg$g_pos)))))
  }
  c(g_occ, g_det, g_pos, eg$g_ld)
}

# Per-species marginal log-likelihood at an arm-split coefficient triple.
.occu_cover_sp_ll <- function(m, bo, bp, bpos, log_disp) {
  eta <- .occu_cover_eta_from_par(m, bo, bp, bpos)
  sum(.occu_cover_site_ll(m, eta$psi, eta$p_mat, eta$ep_mat, log_disp))
}

# Per-species packed coefficient gradient at an arm-split coefficient triple.
.occu_cover_sp_grad <- function(m, bo, bp, bpos, log_disp) {
  eta <- .occu_cover_eta_from_par(m, bo, bp, bpos)
  eg  <- .occu_cover_eta_grad(m, eta$psi, eta$p_mat, eta$ep_mat, log_disp)
  .occu_cover_coef_grad(m, eg)
}


# ---------------------------------------------------------------------------
# Laplace-EM fitter
# ---------------------------------------------------------------------------

# Gauss-Hermite nodes/weights (physicists', weight exp(-x^2)) by Golub-Welsch:
# eigendecomposition of the symmetric Jacobi matrix of the Hermite recurrence.
# Returns nodes `x` and weights `w` with sum(w) = sqrt(pi).
.ms_gh_quad <- function(n) {
  if (n <= 1L) return(list(x = 0, w = sqrt(pi)))
  i <- seq_len(n - 1L)
  b <- sqrt(i / 2)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1L)] <- b
  J[cbind(i + 1L, i)] <- b
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(x = e$values[ord], w = (sqrt(pi) * e$vectors[1L, ]^2)[ord])
}

# AGHQ debias of the community covariance (the single-arm #47 fix generalised
# to the joint community RE) lives in R/community_em.R as
# `.tobs_cem_aghq_sigma()`/`.tobs_cem_reproject_cinv()`, run inside
# `.tobs_community_em()` for every consumer of the shared engine, this family
# included (#269). `arm_idx` here is `list(occ=, p=, pos=)`, so the shared
# function's arm-agnostic `lapply()` returns exactly the `list(occ=, p=,
# pos=)` this family's callers expect.

# Fit the community joint occupancy-cover model via the shared community
# Laplace-EM engine (`.tobs_community_em()`, R/community_em.R): the RE arms
# (occ, p, pos) live in `theta`, the shared cover log-dispersion is the one
# `global`. `model` is the bound ms_occu_cover model. Returns a `tobs_fit`
# (via build_ms_occu_cover_fit).
.tobs_fit_ms_occu_cover <- function(model,
                                    priors     = NULL,
                                    max.iter   = 200L,
                                    tol        = 1e-4,
                                    sigma.beta = 5,
                                    verbose    = TRUE,
                                    ...) {
  dots <- list(...)
  newton.max <- as.integer(dots$newton.max %||% 30L)

  pi_list <- model$process_info
  P_occ <- pi_list[[1L]]$p
  P_p   <- pi_list[[2L]]$p
  P_pos <- pi_list[[3L]]$p
  P     <- P_occ + P_p + P_pos
  S     <- model$n_species

  occ_idx <- seq_len(P_occ)
  p_idx   <- P_occ + seq_len(P_p)
  pos_idx <- P_occ + P_p + seq_len(P_pos)
  arm_idx <- list(occ = occ_idx, p = p_idx, pos = pos_idx)

  views <- lapply(seq_len(S), function(s) .ms_occu_cover_species_view(model, s))

  split_theta <- function(theta) list(bo = theta[occ_idx], bp = theta[p_idx],
                                       bpos = theta[pos_idx])
  sp_ll <- function(s, theta, global) {
    th <- split_theta(theta)
    .occu_cover_sp_ll(views[[s]], th$bo, th$bp, th$bpos, global)
  }
  sp_grad <- function(s, theta, global) {
    th <- split_theta(theta)
    .occu_cover_sp_grad(views[[s]], th$bo, th$bp, th$bpos, global)
  }
  # No analytic sp_info: leaves the engine's own FD-of-sp_grad path (byte-
  # identical to this family's former bespoke copy of the same finite
  # difference).

  # ---- warm start ----
  is_beta  <- identical(model$positive, "beta")
  is_gauss <- identical(model$positive, "gaussian")
  any_det_all <- mean(vapply(views, function(v) mean(rowSums(v$y * v$valid) > 0),
                             numeric(1)))
  pos_vals <- unlist(lapply(views, function(v) v$y_pos[v$valid & v$y == 1L]))
  mu <- numeric(P)
  mu[occ_idx][1L] <- stats::qlogis(min(max(any_det_all, 1e-3), 1 - 1e-3))
  ld <- if (is_beta) log(10) else log(0.4)
  if (length(pos_vals) > 0L) {
    if (is_beta) {
      mu[pos_idx][1L] <- stats::qlogis(min(max(mean(pos_vals), 1e-3), 1 - 1e-3))
    } else if (is_gauss) {
      mu[pos_idx][1L] <- mean(pos_vals)
      ld <- log(stats::sd(pos_vals) + 0.1)
    } else {
      mu[pos_idx][1L] <- mean(log(pos_vals))
      ld <- log(stats::sd(log(pos_vals)) + 0.1)
    }
  }

  # AGHQ variance-component debias is ON by default for this family (unlike
  # ms_occu(): the joint 3-arm occ/p/pos community RE is where the shared
  # engine's re_aghq / re_aghq_maxdim / init_b / init_Sigma hooks were
  # developed and first validated). Disable with control$re.aghq = FALSE, set
  # nodes with control$n.quad; larger RE dims keep the EM covariance + the
  # attenuation flag.
  re_aghq  <- !isFALSE(dots$re.aghq)
  aghq_nq  <- as.integer(dots$n.quad %||% .tobs_n_quad("ms_occu_cover"))
  aghq_cap <- as.integer(dots$re.aghq.maxdim %||% 4L)

  fit <- .tobs_community_em(
    S = S, P = P, arm_idx = arm_idx,
    sp_ll = sp_ll, sp_grad = sp_grad,
    init_mu = mu, init_global = ld,
    penalize_global = FALSE, sigma_beta = sigma.beta, priors = priors,
    sigma_init = 0.3, max_iter = as.integer(max.iter), tol = as.numeric(tol),
    newton_max = newton.max, verbose = isTRUE(verbose),
    re_aghq = re_aghq, n_quad = aghq_nq, re_aghq_maxdim = aghq_cap
  )

  build_ms_occu_cover_fit(model, fit$mu, fit$global, fit$b_list, fit$Sigma,
                          fit$Cinv, fit$Bf, fit$Vf, arm_idx, F_val = fit$logML,
                          converged = fit$converged, n_iter = fit$n_iter,
                          debias_method = fit$debias_method)
}

# ---------------------------------------------------------------------------
# Wrap the EM output into a tobs_fit
# ---------------------------------------------------------------------------

build_ms_occu_cover_fit <- function(model, mu, ld, b_list, Sigma, Cinv_list,
                                    Bf_list, Vf, arm_idx, F_val, converged,
                                    n_iter, debias_method = "none") {
  pi_list <- model$process_info
  P_occ <- pi_list[[1L]]$p
  P_p   <- pi_list[[2L]]$p
  P_pos <- pi_list[[3L]]$p
  P     <- P_occ + P_p + P_pos
  is_beta <- identical(model$positive, "beta")

  beta_names <- c(
    paste0("psi_", pi_list[[1L]]$coef_names),
    paste0("p_",   pi_list[[2L]]$coef_names),
    paste0("pos_", pi_list[[3L]]$coef_names)
  )
  disp_name <- if (is_beta) "log_phi" else "log_sigma_pos"
  par_names <- c(beta_names, disp_name)

  # Per-species community structure (mu + BLUP deviations) per arm.
  B <- do.call(rbind, b_list)                 # S x P
  arm_block <- function(arm) {
    idx <- arm_idx[[arm]]
    blup <- B[, idx, drop = FALSE]
    coef <- sweep(blup, 2L, mu[idx], "+")
    rownames(blup) <- rownames(coef) <- model$species_names
    list(blup = blup, coef = coef)
  }
  occ_b <- arm_block("occ"); p_b <- arm_block("p"); pos_b <- arm_block("pos")
  colnames(occ_b$blup) <- colnames(occ_b$coef) <- pi_list[[1L]]$coef_names
  colnames(p_b$blup)   <- colnames(p_b$coef)   <- pi_list[[2L]]$coef_names
  colnames(pos_b$blup) <- colnames(pos_b$coef) <- pi_list[[3L]]$coef_names

  Sigma_occ <- Sigma$occ; Sigma_p <- Sigma$p; Sigma_pos <- Sigma$pos
  dimnames(Sigma_occ) <- list(pi_list[[1L]]$coef_names, pi_list[[1L]]$coef_names)
  dimnames(Sigma_p)   <- list(pi_list[[2L]]$coef_names, pi_list[[2L]]$coef_names)
  dimnames(Sigma_pos) <- list(pi_list[[3L]]$coef_names, pi_list[[3L]]$coef_names)

  ms_community <- list(
      Sigma_occ = Sigma_occ, Sigma_p = Sigma_p, Sigma_pos = Sigma_pos,
      sd_occ = sqrt(pmax(diag(Sigma_occ), 0)),
      sd_p   = sqrt(pmax(diag(Sigma_p),   0)),
      sd_pos = sqrt(pmax(diag(Sigma_pos), 0)),
      coef_occ = occ_b$coef, coef_p = p_b$coef, coef_pos = pos_b$coef,
      blup_occ = occ_b$blup, blup_p = p_b$blup, blup_pos = pos_b$blup,
      # Per-species posterior covariance Cov(b_s|y) (Louis 1982, from the
      # community EM's own Newton solve, conditional on the converged
      # community mean) -- what a per-species-coefficient consumer (SBC's
      # "rank a fixed species set" design, a calibrated per-species CI) needs
      # beyond the point BLUP; not previously exposed on the fit object.
      # Covers the full b_s vector across all three arms (occ + p + pos). Bf
      # = the (mu,log_disp)-b_s cross-Hessian block from the same Newton
      # solve: mu/log_disp and b_s are NOT independent in the posterior, and
      # Bf is what lets a consumer draw them jointly instead -- see
      # .tobs_sbc_community_b_draws (R/sbc.R). NULL on a NUTS fit (no Newton
      # solve to read it from).
      Cinv = Cinv_list, Bf = Bf_list,
      # The community-MEAN estimates (coef / vcov / confint) are unbiased. The
      # community VARIANCE components (Sigma_occ/Sigma_p/Sigma_pos and their
      # sd_*) carry Laplace small-cluster attenuation at small per-species n.
      # When `debias_method == "aghq"` they have been debiased by adaptive
      # Gauss-Hermite quadrature of the exact per-species RE posterior
      # (generalising the single-arm #47 fix); otherwise (large RE dim /
      # re.aghq = FALSE) they remain the attenuated EM lower bound.
      var_attenuation = if (identical(debias_method, "aghq")) list(
        affects = character(0),
        means_affected = FALSE,
        source = "laplace_small_cluster",
        debias = "aghq",
        note = paste0(
          "Community variance components debiased by adaptive Gauss-Hermite ",
          "quadrature of the exact per-species RE posterior; ",
          "community means are unaffected.")
      ) else list(
        affects = c("Sigma_occ", "Sigma_p", "Sigma_pos",
                    "sd_occ", "sd_p", "sd_pos"),
        means_affected = FALSE,
        source = "laplace_small_cluster",
        debias = "none",
        note = paste0(
          "Community variance components carry Laplace small-cluster ",
          "attenuation (reported as a lower bound); community means are ",
          "unaffected. AGHQ debias applies for small RE dim; ",
          "this fit kept the EM variance (re.aghq = FALSE or RE dim too large).")
      )
  )

  .tobs_cem_finalize_fit(
    means = c(mu, ld), V = Vf, par_names = par_names,
    model = model, process_info = pi_list, N = sum(model$valid),
    log_prob_val = F_val, converged = converged, n_iter = n_iter,
    extra = list(positive = model$positive, ms_community = ms_community)
  )
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "ms_occu_cover")
# ---------------------------------------------------------------------------

# Per-species BLUP deviations, long form: one row per (species, arm, term).
.tobs_ranef_ms_occu_cover <- function(object) {
  .tobs_ranef_ms_long(object$ms_community,
                      c(psi = "blup_occ", p = "blup_p", pos = "blup_pos"))
}

# Per-species posterior-mean linear predictors: site-level occupancy psi
# [n_sites x n_species], site-level detection p [n_sites x n_species], and the
# expected positive cover on the response scale [n_sites x n_species].
# Positive-arm cover on the response scale: beta -> plogis(eta), gaussian ->
# eta, lognormal -> exp(eta + sigma^2 / 2) at the fit's own dispersion. One
# reader, so fitted() and predict(newdata = ) cannot drift apart.
.tobs_ms_cover_response <- function(eta, object) {
  pos <- object$model$positive
  if (identical(pos, "beta")) return(stats::plogis(eta))
  if (identical(pos, "gaussian")) return(eta)
  exp(eta + exp(object$means[[length(object$means)]])^2 / 2)
}

.tobs_fitted_ms_occu_cover <- function(object) {
  model <- object$model
  cm    <- object$ms_community
  X_occ      <- model$X_occ
  X_det_site <- model$X_det_site
  X_pos_site <- model$X_pos_site
  p_det_site <- ncol(X_det_site)
  p_pos_site <- ncol(X_pos_site)
  psi <- stats::plogis(X_occ %*% t(cm$coef_occ))
  p   <- stats::plogis(X_det_site %*%
                       t(cm$coef_p[, seq_len(p_det_site), drop = FALSE]))
  eta_pos <- X_pos_site %*% t(cm$coef_pos[, seq_len(p_pos_site), drop = FALSE])
  cover <- .tobs_ms_cover_response(eta_pos, object)
  dimnames(psi) <- dimnames(p) <- dimnames(cover) <-
    list(NULL, model$species_names)
  list(psi = psi, p = p, cover = cover)
}

# residuals() for the community occupancy-cover families (the spatial-factor
# fit aliases onto this): the state-level residual of the fitted per-species
# occupancy against the ever-detected indicator, the same surface the community
# occupancy families report. The cover arm is observed only where the species
# was detected, so its residual is not a per-site series and is not reported
# here; `det` is NULL for the same reason it is on the community occupancy
# families (no smoothed per-visit state is stored).
.tobs_residuals_ms_occu_cover <- function(object, type) {
  psi <- fitted(object)$psi
  list(occ = .tobs_resid_binary(.tobs_community_ever_detected(object$model),
                                psi, type),
       det = NULL)
}

# Draw community joint occupancy-cover data under the fitted per-species
# coefficients, at the observed visit pattern. Returns 3D arrays matching the
# input y / y_pos.
.tobs_simulate_ms_occu_cover <- function(object, nsim = 1) {
  model <- object$model
  cm    <- object$ms_community
  n_sites <- model$n_sites; max_visits <- model$max_visits
  n_species <- model$n_species
  pos_code <- .occu_cover_pos_code(model$positive)
  disp <- exp(object$means[[length(object$means)]])
  cl <- .tobs_clamp_eta

  # Per-species predictors (community means, deterministic); the z + detection +
  # cover draws run in cpp_simulate_ms_occu_cover from R's RNG stream in the
  # former order (byte-identical).
  psi <- matrix(0, n_sites, n_species)
  p_mat <- array(0, c(n_sites, max_visits, n_species))
  ep_mat <- array(0, c(n_sites, max_visits, n_species))
  for (s in seq_len(n_species)) {
    eta <- .occu_cover_eta_from_par(model, cm$coef_occ[s, ], cm$coef_p[s, ],
                                    cm$coef_pos[s, ])
    psi[, s] <- eta$psi; p_mat[, , s] <- eta$p_mat; ep_mat[, , s] <- eta$ep_mat
  }
  res <- cpp_simulate_ms_occu_cover(psi, as.numeric(p_mat), as.numeric(ep_mat),
    as.integer(model$valid), as.numeric(disp), pos_code,
    n_sites, max_visits, n_species, as.integer(nsim))
  res <- lapply(res, function(r) {
    dn <- list(NULL, NULL, model$species_names)
    dimnames(r$y) <- dn; dimnames(r$y_pos) <- dn; r
  })
  if (nsim == 1L) res[[1]] else res
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate community (multispecies) joint occupancy-cover data
#'
#' Per-species joint occupancy-detection + cover model with Gaussian community
#' hyperpriors on the per-species coefficients of all three arms:
#' `beta_occ_s ~ N(mu_occ, diag(sd_occ^2))`,
#' `beta_p_s ~ N(mu_p, diag(sd_p^2))`,
#' `beta_pos_s ~ N(mu_pos, diag(sd_pos^2))`, then per cell
#' `z_{s,i} ~ Bernoulli(psi_{s,i})`, per visit
#' `y_{s,i,j} | z = 1 ~ Bernoulli(p_{s,i,j})`, and per detected visit
#' `c_{s,i,j} ~ f_pos(eta_pos_{s,i,j}, disp)`. The site and visit covariates are
#' shared across species (community covariates). The returned `y` / `y_pos` are
#' 3D arrays `[n_sites x J x n_species]` suitable for [tobs()] with
#' [ms_occu_cover()] (`y_pos` is `NA` where not detected).
#'
#' @param n_species Number of species (default 12).
#' @param N Number of sites / cells (default 120).
#' @param J Number of replicate visits (default 5).
#' @param n_occ_covs,n_det_covs,n_pos_covs Number of covariates on each arm
#'   (drawn IID standard normal, shared across species).
#' @param mu_occ,mu_p,mu_pos Community-mean coefficient vectors
#'   `c(intercept, slopes...)` on each arm's link scale. Defaults pick
#'   weakly-informative values.
#' @param sd_occ,sd_p,sd_pos Per-coefficient community SD on each arm (length 1,
#'   recycled, or one per coefficient). Default 0.5 / 0.4 / 0.4.
#' @param positive `"lognormal"` (default) or `"beta"`.
#' @param phi Beta precision when `positive = "beta"` (default 30).
#' @param sigma_pos Lognormal residual SD when `positive = "lognormal"`
#'   (default 0.4).
#' @param seed Optional integer seed.
#' @return A list with `y` (3D detection array), `y_pos` (3D cover array, `NA`
#'   where not detected), `data` (per-cell covariate frame), `visit_data`
#'   (per-visit covariate frame, `N*J` rows in site-major order), `species`
#'   (species names), and `truth` (community means / SDs, per-species
#'   coefficients, the dispersion, and the latent state).
#' @export
simulate_ms_occu_cover <- function(n_species = 12, N = 120, J = 5,
                                   n_occ_covs = 1, n_det_covs = 1, n_pos_covs = 1,
                                   mu_occ = NULL, mu_p = NULL, mu_pos = NULL,
                                   sd_occ = 0.5, sd_p = 0.4, sd_pos = 0.4,
                                   positive = c("lognormal", "beta", "gaussian"),
                                   phi = 30, sigma_pos = 0.4, seed = NULL) {
  positive <- match.arg(positive)
  if (!is.null(seed)) set.seed(seed)
  N <- as.integer(N); J <- as.integer(J)
  is_beta  <- identical(positive, "beta")
  is_gauss <- identical(positive, "gaussian")

  if (is.null(mu_occ)) mu_occ <- c(stats::qlogis(0.4), rep(0.6, n_occ_covs))
  if (is.null(mu_p))   mu_p   <- c(0.0, rep(-0.3, n_det_covs))
  if (is.null(mu_pos)) {
    pos_int <- if (is_beta) stats::qlogis(0.3) else if (is_gauss) 2 else log(0.1)
    mu_pos <- c(pos_int, rep(0.4, n_pos_covs))
  }
  P_occ <- length(mu_occ); P_p <- length(mu_p); P_pos <- length(mu_pos)
  rec <- function(sd, p) if (length(sd) == 1L) rep(sd, p) else sd
  sd_occ <- rec(sd_occ, P_occ); sd_p <- rec(sd_p, P_p); sd_pos <- rec(sd_pos, P_pos)

  # Per-species coefficients = community mean + Gaussian deviation.
  draw_beta <- function(mu, sd) {
    matrix(stats::rnorm(n_species * length(mu), 0, rep(sd, each = n_species)),
           n_species, length(mu)) +
      matrix(mu, n_species, length(mu), byrow = TRUE)
  }
  beta_occ <- draw_beta(mu_occ, sd_occ)
  beta_p   <- draw_beta(mu_p,   sd_p)
  beta_pos <- draw_beta(mu_pos, sd_pos)

  # Shared (community) covariates.
  make_covs <- function(n, prefix) {
    if (n <= 0L) return(NULL)
    df <- as.data.frame(matrix(stats::rnorm(N * n), N, n))
    names(df) <- paste0(prefix, seq_len(n)); df
  }
  occ_covs <- make_covs(n_occ_covs, "occ_cov")
  data <- if (is.null(occ_covs)) data.frame(row.names = seq_len(N)) else occ_covs
  X_occ <- if (is.null(occ_covs)) stats::model.matrix(~ 1, data.frame(row.names = seq_len(N)))
           else stats::model.matrix(~ ., occ_covs)

  det_covs <- if (n_det_covs <= 0L) NULL
              else { df <- as.data.frame(matrix(stats::rnorm(N * J * n_det_covs),
                                                N * J, n_det_covs))
                     names(df) <- paste0("det_cov", seq_len(n_det_covs)); df }
  pos_covs <- if (n_pos_covs <= 0L) NULL
              else { df <- as.data.frame(matrix(stats::rnorm(N * J * n_pos_covs),
                                                N * J, n_pos_covs))
                     names(df) <- paste0("pos_cov", seq_len(n_pos_covs)); df }
  visit_data <- do.call(cbind, Filter(Negate(is.null), list(det_covs, pos_covs)))
  if (is.null(visit_data)) visit_data <- data.frame(row.names = seq_len(N * J))
  X_p   <- if (is.null(det_covs)) stats::model.matrix(~ 1, data.frame(row.names = seq_len(N * J)))
           else stats::model.matrix(~ ., det_covs)
  X_pos <- if (is.null(pos_covs)) stats::model.matrix(~ 1, data.frame(row.names = seq_len(N * J)))
           else stats::model.matrix(~ ., pos_covs)

  species_names <- paste0("sp", seq_len(n_species))
  y     <- array(NA_integer_, dim = c(N, J, n_species),
                 dimnames = list(NULL, NULL, species_names))
  y_pos <- array(NA_real_,    dim = c(N, J, n_species),
                 dimnames = list(NULL, NULL, species_names))
  z_all <- matrix(NA_integer_, N, n_species)

  for (s in seq_len(n_species)) {
    psi <- stats::plogis(as.numeric(X_occ %*% beta_occ[s, ]))
    eta_p   <- as.numeric(X_p   %*% beta_p[s, ])
    eta_pos <- as.numeric(X_pos %*% beta_pos[s, ])
    z <- stats::rbinom(N, 1L, psi); z_all[, s] <- z
    for (i in seq_len(N)) {
      if (z[i] == 0L) { y[i, , s] <- 0L; next }
      for (j in seq_len(J)) {
        idx <- (i - 1L) * J + j
        d <- stats::rbinom(1L, 1L, stats::plogis(eta_p[idx]))
        y[i, j, s] <- d
        if (d == 1L) {
          y_pos[i, j, s] <- if (is_beta) {
            mu <- stats::plogis(eta_pos[idx]); stats::rbeta(1L, mu * phi, (1 - mu) * phi)
          } else if (is_gauss) {
            stats::rnorm(1L, eta_pos[idx], sigma_pos)
          } else exp(stats::rnorm(1L, eta_pos[idx], sigma_pos))
        }
      }
    }
  }

  list(
    y = y, y_pos = y_pos, data = data, visit_data = visit_data,
    species = species_names,
    truth = list(
      mu_occ = mu_occ, mu_p = mu_p, mu_pos = mu_pos,
      sd_occ = sd_occ, sd_p = sd_p, sd_pos = sd_pos,
      beta_occ = beta_occ, beta_p = beta_p, beta_pos = beta_pos,
      z = z_all, positive = positive,
      phi       = if (is_beta) phi else NA_real_,
      sigma_pos = if (is_beta) NA_real_ else sigma_pos,
      disp      = if (is_beta) phi else sigma_pos
    )
  )
}
