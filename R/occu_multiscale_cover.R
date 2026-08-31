# =============================================================================
# occu_multiscale_cover.R
# - three-level occupancy + cover hurdle.
#
# Cell-level occupancy (psi) gates plot-level availability (theta), which gates
# per-visit detection (p) and the cover hurdle (pos). The cell is the areal
# field node; plots are the y-rows; visits are the y-columns. Both z (over
# cells) and a (over plots) marginalize in closed form, so the joint marginal
# log-likelihood is exact (the cell-coupling spec in
# src/cell_coupling_occu_multiscale_cover.h). Spatial-only: the fit is the
# joint nested-Laplace cell-coupling path (.tobs_fit_occu_multiscale_cover_*).
#
# Builder, dispatcher, and simulator. The joint fitter lives in
# occu_multiscale_cover_joint.R.
# =============================================================================


# Build the multiscale model object. `y` / `y_pos` are [n_plots x max_visits];
# `plot_cell` maps each plot (row) to its areal cell (1..n_cells). The psi
# design is per-CELL (cell-level covariates aggregated to the first plot of each
# cell), theta is per-PLOT, p / pos are per-visit (plot-level block + optional
# visit-varying). Detection NAs mask visits out; cover is read only where y == 1.
.tobs_build_occu_multiscale_cover <- function(occ_formula, theta_formula,
                                              det_formula, pos_formula,
                                              data, y, y_pos, plot_cell, n_cells,
                                              positive,
                                              det_visit_formula = NULL,
                                              det_visit_data    = NULL,
                                              pos_visit_formula = NULL,
                                              pos_visit_data    = NULL) {
  if (!is.matrix(y) || !is.matrix(y_pos)) {
    stop("y and y_pos must be matrices (n_plots x max_visits).", call. = FALSE)
  }
  if (!all(dim(y) == dim(y_pos))) {
    stop("y and y_pos must have identical dimensions.", call. = FALSE)
  }
  .tobs_check_site_count(nrow(y), nrow(data), "rows (plots)")
  n_plots    <- nrow(y)
  max_visits <- ncol(y)
  plot_cell  <- as.integer(plot_cell)
  if (length(plot_cell) != n_plots || anyNA(plot_cell) ||
      min(plot_cell) < 1L || max(plot_cell) > n_cells) {
    stop(sprintf(paste0(
      "plot_cell must be an integer cell index in 1..%d, one per plot ",
      "(%d plots)."), n_cells, n_plots), call. = FALSE)
  }

  y_int <- matrix(as.integer(y), n_plots, max_visits)
  valid <- !is.na(y_int)
  if (any(y_int[valid] != 0L & y_int[valid] != 1L)) {
    stop("y must contain only 0, 1, or NA (binary detection per visit).",
         call. = FALSE)
  }
  y_int[!valid] <- 0L

  # Cover values through the same validator the two-level family uses: a
  # detected visit whose cover was not recorded is carried as NA and drops only
  # its cover factor (cover missing at random); undetected visits zero-fill.
  y_pos_num <- .occu_cover_validate_pos_values(
    matrix(as.numeric(y_pos), n_plots, max_visits),
    valid & (y_int == 1L), positive)

  # The state-process formula's FE part is cell-level. Build it over plots, then
  # aggregate to one row per cell (the first plot of each cell). Cells with no
  # plot get the intercept + zeroed covariates (prior-mean prediction node).
  X_occ_plot <- stats::model.matrix(occ_formula, data)
  first_plot <- match(seq_len(n_cells), plot_cell)
  X_psi <- matrix(0, n_cells, ncol(X_occ_plot),
                  dimnames = list(NULL, colnames(X_occ_plot)))
  int_col <- match("(Intercept)", colnames(X_occ_plot))
  if (!is.na(int_col)) X_psi[, int_col] <- 1
  present <- !is.na(first_plot)
  X_psi[present, ] <- X_occ_plot[first_plot[present], , drop = FALSE]

  X_theta    <- stats::model.matrix(theta_formula, data)
  X_p_site   <- stats::model.matrix(det_formula, data)
  X_pos_site <- stats::model.matrix(pos_formula, data)
  X_p_visit  <- .tobs_build_visit_X(det_visit_formula, det_visit_data,
                                    n_plots, max_visits, arm = "detection")
  X_pos_visit <- .tobs_build_visit_X(pos_visit_formula, pos_visit_data,
                                     n_plots, max_visits, arm = "positive cover")

  det_coef_names <- colnames(X_p_site)
  pos_coef_names <- colnames(X_pos_site)
  if (!is.null(X_p_visit))   det_coef_names <- c(det_coef_names, colnames(X_p_visit))
  if (!is.null(X_pos_visit)) pos_coef_names <- c(pos_coef_names, colnames(X_pos_visit))

  structure(list(
    model_type  = "occu_multiscale_cover",
    positive    = positive,
    y           = y_int,
    y_pos       = y_pos_num,
    valid       = valid,
    n_cells     = n_cells,
    n_plots     = n_plots,
    n_sites     = n_plots,   # generic accessors (nobs) read plot-level n
    max_visits  = max_visits,
    plot_cell   = plot_cell,
    X_psi       = X_psi,
    X_theta     = X_theta,
    X_p_site    = X_p_site,
    X_pos_site  = X_pos_site,
    X_p_visit   = X_p_visit,
    X_pos_visit = X_pos_visit,
    formulas    = list(psi = occ_formula, theta = theta_formula,
                       p = det_formula, pos = pos_formula),
    data        = data,
    process_info = list(
      list(name = "psi",   p = ncol(X_psi),
           coef_names = colnames(X_psi),   link = "logit"),
      list(name = "theta", p = ncol(X_theta),
           coef_names = colnames(X_theta), link = "logit"),
      list(name = "p",     p = length(det_coef_names),
           coef_names = det_coef_names,    link = "logit"),
      list(name = "pos",   p = length(pos_coef_names),
           coef_names = pos_coef_names,
           link = if (positive == "beta") "logit" else "identity")
    )
  ), class = "tobs_model")
}


# ---------------------------------------------------------------------------
# S3 helpers (routed from methods.R by model_type == "occu_multiscale_cover")
# ---------------------------------------------------------------------------

# Arm coefficients, split by process block, from a flat coefficient vector in
# the SAME layout as object$means / one row of object$draws -- so this reads
# either the posterior mean or one posterior draw.
.occu_mscale_cover_split_beta <- function(beta, pi_list) {
  pp  <- vapply(pi_list, function(x) as.integer(x$p), integer(1))
  off <- cumsum(c(0L, pp))
  list(psi   = beta[off[1] + seq_len(pp[1])],
       theta = beta[off[2] + seq_len(pp[2])],
       p     = beta[off[3] + seq_len(pp[3])],
       pos   = beta[off[4] + seq_len(pp[4])])
}

# The per-unit hierarchy surface at an arbitrary coefficient vector `beta` (the
# posterior mean for fitted(), or one posterior draw for predict()'s interval).
# psi is per cell (the areal field is added with coefficient 1); theta / p /
# cover are per plot (site-level designs; the visit-level part of p / cover is
# constant across visits at the posterior mean and so reads off the site
# design). The cover arm carries the shared field with the copy coefficient
# `alpha`. Lognormal cover reports the conditional MEAN
# E[cover | detected] = exp(eta + sigma_pos^2 / 2); beta cover reports the
# conditional mean plogis(eta). `field` and `alpha` are handed in separately
# since a draw does not carry its own field realization (see
# .tobs_predict_occu_multiscale_cover()).
.occu_mscale_cover_surface_at <- function(model, beta, field, alpha, sigma_pos) {
  b   <- .occu_mscale_cover_split_beta(beta, model$process_info)
  pc  <- model$plot_cell

  eta_psi  <- as.numeric(model$X_psi %*% b$psi) + field                  # per cell
  psi_cell <- stats::plogis(eta_psi)
  theta    <- stats::plogis(as.numeric(model$X_theta %*% b$theta))       # per plot
  pdet     <- stats::plogis(as.numeric(model$X_p_site %*% b$p))          # per plot
  eta_pos  <- as.numeric(model$X_pos_site %*% b$pos) + alpha * field[pc] # per plot
  cover <- if (identical(model$positive, "beta")) {
    stats::plogis(eta_pos)
  } else if (identical(model$positive, "gaussian")) {
    eta_pos
  } else {
    exp(eta_pos + 0.5 * sigma_pos^2)
  }
  list(psi = psi_cell, theta = theta, p = pdet, cover = cover)
}

# fitted(): per-unit posterior-mean predictions for each level of the
# hierarchy. `p_marginal` is the per-plot marginal detection probability
# psi_cell * theta * p (a single visit).
.tobs_fitted_occu_multiscale_cover <- function(object) {
  model <- object$model
  field <- as.numeric(object$spatial_field %||% rep(0, model$n_cells))   # per cell
  alpha <- unname(object$means["alpha"]); if (!is.finite(alpha)) alpha <- 0
  sigma_pos <- .occu_mscale_cover_sigma_pos(object$means)
  s  <- .occu_mscale_cover_surface_at(model, object$means, field, alpha, sigma_pos)
  pc <- model$plot_cell
  c(s, list(field = field, p_marginal = s$psi[pc] * s$theta * s$p))
}

# Lognormal-cover residual SD at a named coefficient vector (object$means, or
# one row of object$draws), robust to the dispersion naming: the spatial path
# reports `phi_pos` (already sigma_pos on the natural scale), the non-spatial
# path a log-scale `log_sigma_pos`. Returns 0 if neither (beta arm /
# unavailable), giving the conditional median exp(eta).
.occu_mscale_cover_sigma_pos <- function(m) {
  if ("phi_pos" %in% names(m) && is.finite(m[["phi_pos"]])) return(unname(m[["phi_pos"]]))
  if ("log_sigma_pos" %in% names(m) && is.finite(m[["log_sigma_pos"]]))
    return(exp(unname(m[["log_sigma_pos"]])))
  0
}

# residuals() for occu_multiscale_cover(): the CELL-level residual of the
# fitted psi against the ever-detected indicator (any plot in the cell, at
# any visit) -- occupancy is a cell-level state that gates plot availability,
# not something observed at the plot. Availability / detection / cover are
# each observed at a different unit (plot, visit, detected-visit), so none of
# them is a per-cell series and only the state residual is reported; `det` is
# NULL for the same reason it is on occu_cover().
.tobs_residuals_occu_multiscale_cover <- function(object, type) {
  model <- object$model
  y_masked <- model$y
  y_masked[!model$valid] <- NA_integer_
  plot_det <- apply(y_masked, 1, function(row) as.integer(any(row == 1L, na.rm = TRUE)))
  pc <- model$plot_cell
  cell_det <- vapply(seq_len(model$n_cells),
                     function(c) as.integer(any(plot_det[pc == c] == 1L)),
                     integer(1))
  psi <- .tobs_fitted_occu_multiscale_cover(object)$psi
  list(occ = .tobs_resid_binary(cell_det, psi, type), det = NULL)
}

# predict(): in-sample posterior arm predictions, draw-based (mean/sd/interval
# via .occu_cover_summ(), R/occu_cover_predict.R -- the occu_cover / cover
# joint predictor's own summary table). Each arm's coefficient draws come from
# `object$draws` (a Gaussian resample of the arm coefficients on every engine,
# .rmvn(n, means, V) on the joint route), sliced by the same
# process_info-ordered offsets .occu_mscale_cover_split_beta() reads at the
# posterior mean. The shared field / copy amplitude `alpha` are point
# estimates -- not drawn per iteration on any of the three engines -- so they
# enter as a fixed offset, the same way .tobs_fitted_occu_multiscale_cover()
# and .tobs_count_arm_predict() (R/field_offset.R) treat a fitted field. The
# areal field is tied to the cell graph, so prediction at new covariates /
# cells (X.0 / newdata) is not supported (no field at an unseen cell),
# matching ms_occu_cover_spatial().
.tobs_predict_occu_multiscale_cover <- function(object, type = "state",
                                                level = 0.95) {
  type <- switch(type,
                 state = "psi", occupancy = "psi", psi = "psi",
                 availability = "theta", theta = "theta",
                 detection = "p", p = "p",
                 cover = "cover", cover_cond = "cover",
                 stop(sprintf(paste0("predict() type '%s' is not supported for ",
                      "occu_multiscale_cover(). Use \"state\"/\"psi\", ",
                      "\"availability\"/\"theta\", \"detection\"/\"p\", or ",
                      "\"cover\"."), type), call. = FALSE))
  model   <- object$model
  draws   <- object$draws
  n_draws <- nrow(draws)
  pi_list <- model$process_info
  pp  <- vapply(pi_list, function(x) as.integer(x$p), integer(1))
  off <- cumsum(c(0L, pp))
  arm_cols <- function(k) off[k] + seq_len(pp[k])
  field <- as.numeric(object$spatial_field %||% rep(0, model$n_cells))
  alpha <- if ("alpha" %in% colnames(draws)) draws[, "alpha"]
           else rep(0, n_draws)
  pc    <- model$plot_cell

  if (identical(type, "psi")) {
    eta <- draws[, arm_cols(1L), drop = FALSE] %*% t(model$X_psi)
    eta <- sweep(eta, 2L, field, "+")
    return(.occu_cover_summ(t(stats::plogis(eta)), seq_len(model$n_cells), level))
  }
  if (identical(type, "theta")) {
    eta <- draws[, arm_cols(2L), drop = FALSE] %*% t(model$X_theta)
    return(.occu_cover_summ(t(stats::plogis(eta)), pc, level))
  }
  if (identical(type, "p")) {
    eta <- draws[, arm_cols(3L), drop = FALSE] %*% t(model$X_p_site)
    return(.occu_cover_summ(t(stats::plogis(eta)), pc, level))
  }
  # cover
  eta <- draws[, arm_cols(4L), drop = FALSE] %*% t(model$X_pos_site)
  eta <- eta + outer(alpha, field[pc])
  mat <- if (identical(model$positive, "beta")) {
    stats::plogis(eta)
  } else if (identical(model$positive, "gaussian")) {
    eta
  } else {
    sigma_pos <- .occu_mscale_cover_sigma_pos(object$means)
    exp(eta + 0.5 * sigma_pos^2)
  }
  .occu_cover_summ(t(mat), pc, level)
}


# ---------------------------------------------------------------------------
# Non-spatial Laplace path (iid cells, no areal field). The exact three-level
# marginal is optimised directly: z (cells) and a (plots) marginalize in closed
# form, the cover hurdle multiplies the detected-visit branch. Same closed-form
# structure as the joint-coupled engine with the field fixed at 0.
# ---------------------------------------------------------------------------

# Per-arm visit-level eta matrix [n_plots x max_visits]: the site predictor
# broadcast across a plot's visits, plus the optional visit-varying part (whose
# rows are ordered (plot - 1) * max_visits + visit, matching .tobs_build_visit_X).
.occu_ms_eta_visit <- function(Xs, bs, Xv, bv, n_plots, J) {
  eta <- matrix(as.numeric(Xs %*% bs), n_plots, J)
  if (!is.null(Xv) && length(bv) > 0L) {
    ev <- as.numeric(Xv %*% bv)
    eta <- eta + matrix(ev, n_plots, J, byrow = TRUE)
  }
  eta
}

# Exact three-level marginal log-likelihood at the packed parameter vector
# par = (beta_psi, beta_theta, beta_p[site,visit], beta_pos[site,visit], log_disp).
# `idx` carries each block's coordinate indices and the per-arm site/visit split.
# `per_cell = TRUE` returns the length-n_cells per-cell log-likelihood vector
# (the pointwise unit WAIC / LOO score), otherwise their sum.
.occu_mscale_cover_nonspatial_ll <- function(par, model, idx, per_cell = FALSE) {
  J        <- model$max_visits
  n_plots  <- model$n_plots
  n_cells  <- model$n_cells
  pc       <- model$plot_cell
  valid    <- model$valid
  y        <- model$y
  ypos     <- model$y_pos
  is_beta  <- identical(model$positive, "beta")
  is_gauss <- identical(model$positive, "gaussian")

  eta_psi   <- as.numeric(model$X_psi %*% par[idx$psi])           # n_cells
  eta_theta <- as.numeric(model$X_theta %*% par[idx$theta])       # n_plots
  eta_p   <- .occu_ms_eta_visit(model$X_p_site,   par[idx$p_site],
                                model$X_p_visit,  par[idx$p_visit],   n_plots, J)
  eta_pos <- .occu_ms_eta_visit(model$X_pos_site, par[idx$pos_site],
                                model$X_pos_visit, par[idx$pos_visit], n_plots, J)

  clp   <- function(x) pmin(pmax(x, 1e-12), 1 - 1e-12)
  psi   <- clp(stats::plogis(eta_psi))
  theta <- clp(stats::plogis(eta_theta))
  p_mat <- clp(stats::plogis(eta_p))

  log_p   <- ifelse(valid, log(p_mat),     0)
  log_1mp <- ifelse(valid, log(1 - p_mat), 0)
  det_v   <- valid & (y == 1L)
  log_hdet <- ifelse(valid, ifelse(y == 1L, log_p, log_1mp), 0)

  # Cover density at detected visits with an observed cover: a detected visit
  # carrying NA cover keeps its detection factor and drops only this one.
  cov_v    <- det_v & is.finite(ypos)
  cover_lf <- matrix(0, n_plots, J)
  if (any(cov_v)) {
    ep <- eta_pos[cov_v]; cv <- ypos[cov_v]
    # One boundary policy for the positive arm, shared with the coupling
    # kernel and the samplers: `.tobs_log_safe` on every log, no clamp on the
    # RESPONSE. Nudging the cover to 1e-9 evaluated a density at a value the
    # data does not carry; the shared policy keeps the density finite at a
    # cover of exactly 0 or 1 and leaves the response alone.
    cover_lf[cov_v] <- .occu_cover_pos_logdens(
      cv, ep, exp(par[idx$disp]),
      if (is_beta) "beta" else if (is_gauss) "gaussian" else "lognormal")
  }

  sum_hdet  <- rowSums(log_hdet)
  sum_1mp   <- rowSums(log_1mp)
  sum_cover <- rowSums(cover_lf)
  det_plot  <- rowSums(det_v) > 0

  log_theta   <- log(theta); log_1mtheta <- log(1 - theta)
  # Plot log-prob given z = 1: detected -> a = 1 forced; else marginalize a.
  log_pj_det <- log_theta + sum_hdet + sum_cover
  ln_a <- log_theta + sum_1mp; ln_b <- log_1mtheta
  m_p  <- pmax(ln_a, ln_b)
  log_pj_nodet <- m_p + log(exp(ln_a - m_p) + exp(ln_b - m_p))
  log_pj <- ifelse(det_plot, log_pj_det, log_pj_nodet)

  # Aggregate plots to their cells.
  sum_logpj_cell <- numeric(n_cells)
  agg <- rowsum(log_pj, group = pc, reorder = FALSE)
  sum_logpj_cell[as.integer(rownames(agg))] <- agg[, 1L]
  det_cell <- logical(n_cells)
  aggd <- rowsum(as.numeric(det_plot), group = pc, reorder = FALSE)
  det_cell[as.integer(rownames(aggd))] <- aggd[, 1L] > 0

  log_psi <- log(psi); log_1mpsi <- log(1 - psi)
  ll_det <- log_psi + sum_logpj_cell
  lc_a <- log_psi + sum_logpj_cell; lc_b <- log_1mpsi
  m_c  <- pmax(lc_a, lc_b)
  ll_nodet <- m_c + log(exp(lc_a - m_c) + exp(lc_b - m_c))
  ll_cell <- ifelse(det_cell, ll_det, ll_nodet)
  if (per_cell) ll_cell else sum(ll_cell)
}


# Per-cell pointwise log-likelihood over a posterior draw matrix (the WAIC / LOO
# unit is the cell, the top-level marginalised observation). Reuses the exact
# three-level marginal LL with `per_cell = TRUE`, so the pointwise score and the
# fitted log-likelihood share one source of truth. Returns an [n_draws x n_cells]
# matrix. Used for both the NUTS draws (calibrated WAIC) and the Laplace
# pseudo-draws.
# Batched [n_draws x n_cells] per-cell pointwise log-likelihood via the C++
# kernel (cpp_occu_mscale_cover_ploglik), parallel over draws. Mirrors
# .occu_mscale_cover_nonspatial_ll (the R oracle) draw for draw. `draws` is the
# [S x total] parameter matrix.
.occu_mscale_cover_ploglik_core <- function(model, draws, n_threads = 1L) {
  idx <- .tobs_occu_mscale_cover_nuts_layout(model)
  off_w <- function(v) if (length(v) == 0L) c(0L, 0L)
                       else c(as.integer(v[1L]) - 1L, length(v))
  pw <- off_w(idx$psi); tw <- off_w(idx$theta)
  psw <- off_w(idx$p_site); pvw <- off_w(idx$p_visit)
  posw <- off_w(idx$pos_site); posvw <- off_w(idx$pos_visit)
  np <- model$n_plots; J <- model$max_visits
  empty_v <- function(m) if (is.null(m)) matrix(0, np * J, 0L) else m
  y <- model$y; storage.mode(y) <- "integer"
  cpp_occu_mscale_cover_ploglik(
    draws, model$X_psi, model$X_theta, model$X_p_site, empty_v(model$X_p_visit),
    model$X_pos_site, empty_v(model$X_pos_visit),
    y, model$y_pos, model$valid, as.integer(model$plot_cell),
    .occu_cover_pos_code(model$positive),
    pw[1], pw[2], tw[1], tw[2], psw[1], psw[2], pvw[1], pvw[2],
    posw[1], posw[2], posvw[1], posvw[2], as.integer(idx$disp) - 1L,
    max(1L, as.integer(n_threads)))
}

.tobs_ploglik_occu_multiscale_cover <- function(object, n.draws = 1000L,
                                                n.threads = 1L) {
  model <- object$model
  draws <- object$draws
  if (is.null(draws) || !is.matrix(draws)) {
    stop("Pointwise log-likelihood needs a posterior draw matrix; ",
         "`object$draws` is missing or not a matrix.", call. = FALSE)
  }
  if (!is.null(n.draws) && n.draws < nrow(draws)) {
    draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
  }
  .occu_mscale_cover_ploglik_core(model, draws, n.threads)
}

# Non-spatial Laplace fit: BFGS over the exact marginal (numeric gradient),
# observed-information vcov from the marginal Hessian. Returns a tobs_fit shaped
# like the joint-coupled output (minus the field / hyperparameters).
.tobs_fit_occu_multiscale_cover_laplace <- function(model, priors = NULL,
                                                    max.iter = 300L, tol = 1e-7,
                                                    verbose = TRUE, sigma.beta = 5,
                                                    ...) {
  pi_list <- model$process_info
  p_psi   <- pi_list[[1L]]$p; p_theta <- pi_list[[2L]]$p
  p_p     <- pi_list[[3L]]$p; p_pos   <- pi_list[[4L]]$p
  ps_p    <- ncol(model$X_p_site);   ps_pos <- ncol(model$X_pos_site)
  off     <- cumsum(c(0L, p_psi, p_theta, p_p, p_pos))
  idx <- list(
    psi   = off[1] + seq_len(p_psi),
    theta = off[2] + seq_len(p_theta),
    p     = off[3] + seq_len(p_p),
    pos   = off[4] + seq_len(p_pos),
    disp  = off[5] + 1L)
  idx$p_site   <- idx$p[seq_len(ps_p)]
  idx$p_visit  <- if (p_p > ps_p)   idx$p[(ps_p + 1L):p_p]     else integer(0)
  idx$pos_site <- idx$pos[seq_len(ps_pos)]
  idx$pos_visit<- if (p_pos > ps_pos) idx$pos[(ps_pos + 1L):p_pos] else integer(0)
  n_par <- off[5] + 1L
  is_beta <- identical(model$positive, "beta")

  par_names <- c(
    paste0("psi_",   pi_list[[1L]]$coef_names),
    paste0("theta_", pi_list[[2L]]$coef_names),
    paste0("p_",     pi_list[[3L]]$coef_names),
    paste0("pos_",   pi_list[[4L]]$coef_names),
    if (is_beta) "log_phi" else "log_sigma_pos")

  # Warm starts: occupancy / availability / detection intercepts from the
  # empirical any-detection rate; cover intercept + dispersion from the detected
  # positive values, mirroring .tobs_fit_occu_cover.
  start <- numeric(n_par)
  any_det_plot <- rowSums(model$y * model$valid) > 0
  rate <- min(max(mean(any_det_plot), 1e-3), 1 - 1e-3)
  start[idx$psi[1]]   <- stats::qlogis(rate)
  start[idx$theta[1]] <- stats::qlogis(0.7)
  start[idx$p[1]]     <- 0
  is_gauss <- identical(model$positive, "gaussian")
  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  pos_vals <- pos_vals[is.finite(pos_vals)]
  if (length(pos_vals) > 0L) {
    if (is_beta) {
      start[idx$pos[1]] <- stats::qlogis(min(max(mean(pos_vals), 1e-3), 1 - 1e-3))
      start[idx$disp]   <- log(10)
    } else if (is_gauss) {
      start[idx$pos[1]] <- mean(pos_vals)
      start[idx$disp]   <- log(stats::sd(pos_vals) + 0.1)
    } else {
      start[idx$pos[1]] <- mean(log(pos_vals))
      start[idx$disp]   <- log(stats::sd(log(pos_vals)) + 0.1)
    }
  } else {
    start[idx$disp] <- if (is_beta) log(10) else log(0.4)
  }

  # Weakly-informative Gaussian prior on the betas (dispersion stays flat),
  # matching the occu_cover convention. priors = FALSE disables it.
  pprec <- numeric(n_par)
  if (!isFALSE(priors)) {
    beta_idx <- c(idx$psi, idx$theta, idx$p, idx$pos)
    pprec[beta_idx] <- 1 / sigma.beta^2
  }

  neg_pen <- function(par) {
    ll <- .occu_mscale_cover_nonspatial_ll(par, model, idx)
    -(ll - 0.5 * sum(pprec * par^2))
  }
  .prog <- tulpa:::.tulpa_iter_progress("occu-ms-cover-laplace",
                                        as.integer(max.iter), unit = "iter")
  neg_pen_p <- function(par) { .prog$tick(); neg_pen(par) }
  opt <- stats::optim(start, neg_pen_p, method = "BFGS",
                      control = list(maxit = as.integer(max.iter), reltol = tol))
  .prog$finish()
  par <- opt$par
  H   <- stats::optimHess(par, neg_pen)
  V   <- tryCatch(solve(H), error = function(e)
                  matrix(NA_real_, n_par, n_par))
  dimnames(V) <- list(par_names, par_names)
  means <- par; names(means) <- par_names
  sds   <- sqrt(pmax(diag(V), 0)); names(sds) <- par_names

  n_draws <- 1000L
  draws <- .rmvn(n_draws, means, V); colnames(draws) <- par_names
  ll_val <- .occu_mscale_cover_nonspatial_ll(par, model, idx)

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = n_par,
    log_prob     = rep(ll_val, n_draws),
    log_lik      = ll_val,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = NULL,
    spatial_field = NULL,
    method       = "laplace",
    positive     = model$positive,
    convergence  = list(converged = opt$convergence == 0L,
                        n_iter = unname(opt$counts[[1L]]))
  )), class = c("tobs_fit", "tulpa_fit"))
}


# ---------------------------------------------------------------------------
# Dispatcher (wired into tobs.R's switch). Cells are declared via an
# icar()/bym2() term's group_var; method = "nested_laplace" fits the shared
# areal field, method = "laplace" the non-spatial (iid-cell) marginal.
# ---------------------------------------------------------------------------
.dispatch_occu_multiscale_cover <- function(formula, data, family, detection,
                                            y, visits, engine, priors, control,
                                            approx = "gaussian_laplace",
                                            correction = "none", ...) {
  dots <- list(...)
  if (is.null(detection)) {
    stop("occu_multiscale_cover() requires a `detection` formula (the ",
         "per-visit detection arm p).", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu_multiscale_cover() requires `y` (an n_plots x max_visits ",
         "detection-history matrix).", call. = FALSE)
  }
  if (is.null(dots$y_pos)) {
    stop("occu_multiscale_cover() requires `y_pos` (n_plots x max_visits; ",
         "values used only where y == 1).", call. = FALSE)
  }
  # `alpha.grid` states the copy axis's nodes, `alpha.n` how many nodes the
  # engine's own axis is read at; one block takes one of them.
  .tobs_check_alpha_control(control, "occu_multiscale_cover()")

  # Identifiability surfaced. The availability (theta) and detection (p) levels
  # separate only with replicate visits WITHIN a plot. With at most one observed
  # visit per plot (single releves) the fit identifies cell occupancy (psi) and
  # the product theta * p, reducing to occu_cover(). Flag the reduction so the
  # per-level theta / p coefficients are not over-interpreted.
  if (is.matrix(y) && max(rowSums(!is.na(y)), 0L) < 2L) {
    message("occu_multiscale_cover(): the data carry no within-plot ",
            "replication (at most one visit per plot), so availability (theta) ",
            "and detection (p) are not separately identified -- the fit ",
            "identifies cell occupancy (psi) and the product theta * p. Add ",
            "within-plot replicate visits to separate the two levels, or use ",
            "occu_cover() if a two-level model is intended.")
  }

  theta_formula <- dots$availability %||% ~ 1
  pos_formula   <- dots$positive %||% detection
  is_nuts       <- identical(engine, "nuts")
  non_spatial   <- identical(engine, "laplace") || is_nuts

  spatial_info <- .occu_cover_spatial_fields(formula, data)
  if (is.null(spatial_info)) {
    stop("occu_multiscale_cover() declares cells via an areal term: the state ",
         "formula must carry icar(graph = adj, group_var = \"<cell>\") naming ",
         "the per-plot cell column (the graph is used for the field under ",
         "method = \"nested_laplace\" and ignored under method = \"laplace\").",
         call. = FALSE)
  }
  if (isTRUE(spatial_info$correlated)) {
    stop(paste0(
      "occu_multiscale_cover(): a correlated spatial bar (`|`, free-Sigma MCAR) ",
      "is wired on occu_cover()/cover(), not the multiscale path. Use the ",
      "INDEPENDENT spelling `||` (or the two-term areal form) here."),
      call. = FALSE)
  }
  gv <- spatial_info$group_var
  if (is.null(gv)) {
    stop("occu_multiscale_cover() requires group_var on the icar()/bym2() ",
         "term naming the integer cell column (plot -> cell), e.g. ",
         "icar(graph = adj, group_var = \"cell\").", call. = FALSE)
  }
  if (!gv %in% names(data)) {
    stop(sprintf("occu_multiscale_cover() group_var '%s' is not a column of data.",
                 gv), call. = FALSE)
  }

  # Detection / availability / cover arms carry no structured terms.
  .occu_cover_reject_structured(detection,     "detection")
  .occu_cover_reject_structured(theta_formula, "availability")
  .occu_cover_reject_structured(pos_formula,   "positive cover")

  fields  <- spatial_info$fields
  n_cells <- nrow(fields[[1L]]$graph)
  plot_cell <- as.integer(data[[gv]])

  # method = "nuts" samples the EXACT non-spatial three-level marginal (iid
  # cells, field fixed at 0); the cell-declaring icar()/bym2() term supplies the
  # plot -> cell map only, its graph is ignored. The shared areal field (and any
  # coupled SVC / trend field) is grid-integrated under method = "nested_laplace"
  # and has no sampled-field route here -- gate it with a pointer.
  if (is_nuts && length(fields) > 1L) {
    stop("occu_multiscale_cover() method = \"nuts\" is the non-spatial path ",
         "(iid cells, no areal field); a coupled SVC / trend field is not ",
         "sampled. Use method = \"nested_laplace\" for the shared / trend ",
         "field, or keep a single cell-declaring icar()/bym2() term for the ",
         "non-spatial NUTS fit.", call. = FALSE)
  }

  vd_det <- .normalize_visits(visits, detection,
                              n_sites = nrow(y), max_visits = ncol(y))
  vd_pos <- .normalize_visits(visits, pos_formula,
                              n_sites = nrow(y), max_visits = ncol(y))

  model <- .tobs_build_occu_multiscale_cover(
    occ_formula       = spatial_info$fe,
    theta_formula     = theta_formula,
    det_formula       = vd_det$det_formula,
    pos_formula       = vd_pos$det_formula,
    data              = data,
    y                 = y,
    y_pos             = dots$y_pos,
    plot_cell         = plot_cell,
    n_cells           = n_cells,
    positive          = family$params$positive,
    det_visit_formula = vd_det$det_visit_formula,
    det_visit_data    = vd_det$visits,
    pos_visit_formula = vd_pos$det_visit_formula,
    pos_visit_data    = vd_pos$visits
  )

  control[["engine"]] <- NULL
  if (is_nuts) {
    fit_args <- c(list(model = model, priors = priors), control)
    return(do.call(.tobs_fit_occu_multiscale_cover_nuts, fit_args))
  }
  if (non_spatial) {
    fit_args <- c(list(model = model, priors = priors), control)
    return(do.call(.tobs_fit_occu_multiscale_cover_laplace, fit_args))
  }
  fit_args <- c(list(model = model, fields = fields, priors = priors), control)
  do.call(.tobs_fit_occu_multiscale_cover_joint, fit_args)
}


# ---------------------------------------------------------------------------
# Simulator
# ---------------------------------------------------------------------------

#' Simulate a three-level occupancy + cover hurdle data set
#'
#' Draws cell occupancy `z_c ~ Bernoulli(psi_c)`, plot availability
#' `a_cj | z_c=1 ~ Bernoulli(theta_cj)`, visit detection
#' `y_cjv | a_cj=1 ~ Bernoulli(p_cjv)`, and the cover hurdle
#' `cover_cjv | y_cjv=1 ~ f_pos`, with a shared areal ICAR field on the
#' occupancy (`sigma`) and cover (`alpha * sigma`) arms. The data-generating
#' process for [occu_multiscale_cover()] recovery tests.
#'
#' @param n_cells number of areal cells (graph nodes).
#' @param plots_per_cell plots (the availability units) per cell.
#' @param visits_per_plot detection replicate visits per plot.
#' @param beta_psi,beta_theta,beta_p,beta_pos length-2 `c(intercept, slope)`
#'   coefficients for the four arms (one covariate each).
#' @param positive `"lognormal"` or `"beta"` cover arm.
#' @param phi cover dispersion (lognormal log-scale SD, or beta precision).
#' @param adj `n_cells x n_cells` adjacency; default a 1-D chain.
#' @param sigma areal-field marginal SD on the occupancy arm.
#' @param alpha cover-arm field scaling (`alpha * sigma` on the cover arm).
#' @param trend add a second, spatially-varying-coefficient (SVC) areal field
#'   weighted per cell by a cell-level covariate `tcov` (added to `data`):
#'   `tcov_c * sigma_trend * f_trend[c]` on occupancy and
#'   `tcov_c * alpha_trend * sigma_trend * f_trend[c]` on cover. Default `FALSE`.
#' @param sigma_trend,alpha_trend trend-field marginal SD and its cover-arm
#'   scaling (used only when `trend = TRUE`).
#' @param seed optional RNG seed.
#' @return A list with `y`, `y_pos` (`[n_plots x visits_per_plot]`), the
#'   plot-level `data` (cell id `cell`, covariates), `adj`, and `truth`.
#' @export
simulate_occu_multiscale_cover <- function(n_cells = 60L,
                                           plots_per_cell = 4L,
                                           visits_per_plot = 2L,
                                           beta_psi   = c(0.4, 0.6),
                                           beta_theta = c(0.2, 0.5),
                                           beta_p     = c(0.0, 0.5),
                                           beta_pos   = c(log(0.10), -0.4),
                                           positive   = c("lognormal", "beta",
                                                           "gaussian"),
                                           phi        = 0.35,
                                           adj        = NULL,
                                           sigma      = 0.7,
                                           alpha      = 1.0,
                                           trend      = FALSE,
                                           sigma_trend = 0.7,
                                           alpha_trend = 1.0,
                                           seed       = NULL) {
  positive <- match.arg(positive)
  if (!is.null(seed)) set.seed(seed)
  n_cells <- as.integer(n_cells)
  J <- as.integer(visits_per_plot)

  if (is.null(adj)) {
    adj <- matrix(0, n_cells, n_cells)
    for (i in seq_len(n_cells - 1L)) { adj[i, i + 1L] <- 1; adj[i + 1L, i] <- 1 }
  }

  # ICAR field, Sorbye-Rue scaled (geo-mean marginal variance 1).
  Q       <- .occu_cover_icar_Q(adj)
  scale_q <- .occu_cover_icar_scale(adj)
  draw_field <- function() as.numeric(.tobs_draw_icar_unit(Q, scale_q))
  f <- draw_field()

  # Optional spatially-varying trend: a SECOND ICAR field f_trend, weighted per
  # cell by a cell-level covariate `tcov`, adding tcov_c * sigma_trend * f_trend
  # to occupancy and tcov_c * alpha_trend * sigma_trend * f_trend to cover.
  f_trend <- NULL; tcov <- NULL
  if (isTRUE(trend)) {
    f_trend <- draw_field()
    tcov    <- stats::rnorm(n_cells)
  }

  # Plot -> cell map (balanced).
  plot_cell <- rep(seq_len(n_cells), each = plots_per_cell)
  n_plots   <- length(plot_cell)

  # Covariates: cell-level (psi), plot-level (theta, p, pos).
  x_cell  <- stats::rnorm(n_cells)
  x_plot  <- stats::rnorm(n_plots)
  x_pdet  <- stats::rnorm(n_plots)
  x_cov   <- stats::rnorm(n_plots)

  # Cell occupancy (intercept field + optional cell-level SVC trend).
  eta_psi <- beta_psi[1L] + beta_psi[2L] * x_cell + sigma * f
  if (isTRUE(trend)) eta_psi <- eta_psi + tcov * sigma_trend * f_trend
  psi     <- stats::plogis(eta_psi)
  z       <- stats::rbinom(n_cells, 1L, psi)

  # Plot availability (only meaningful where the cell is occupied).
  eta_theta <- beta_theta[1L] + beta_theta[2L] * x_plot
  theta     <- stats::plogis(eta_theta)
  a         <- stats::rbinom(n_plots, 1L, theta) * z[plot_cell]

  # Detection + cover per visit.
  y     <- matrix(0L, n_plots, J)
  y_pos <- matrix(NA_real_, n_plots, J)
  eta_p   <- beta_p[1L]   + beta_p[2L]   * x_pdet
  p_plot  <- stats::plogis(eta_p)
  eta_pos <- beta_pos[1L] + beta_pos[2L] * x_cov + alpha * sigma * f[plot_cell]
  if (isTRUE(trend)) {
    eta_pos <- eta_pos +
      tcov[plot_cell] * alpha_trend * sigma_trend * f_trend[plot_cell]
  }
  for (i in seq_len(n_plots)) {
    if (a[i] == 1L) {
      yi <- stats::rbinom(J, 1L, p_plot[i])
      y[i, ] <- yi
      det <- which(yi == 1L)
      if (length(det)) {
        if (positive == "beta") {
          mu <- stats::plogis(eta_pos[i])
          y_pos[i, det] <- stats::rbeta(length(det), mu * phi, (1 - mu) * phi)
        } else if (positive == "gaussian") {
          y_pos[i, det] <- stats::rnorm(length(det), eta_pos[i], phi)
        } else {
          y_pos[i, det] <- stats::rlnorm(length(det), eta_pos[i], phi)
        }
      }
    }
  }

  data <- data.frame(cell = plot_cell, x_cell = x_cell[plot_cell],
                     x_plot = x_plot, x_pdet = x_pdet, x_cov = x_cov)
  if (isTRUE(trend)) data$tcov <- tcov[plot_cell]

  list(
    y     = y,
    y_pos = y_pos,
    data  = data,
    adj   = adj,
    truth = list(beta_psi = beta_psi, beta_theta = beta_theta,
                 beta_p = beta_p, beta_pos = beta_pos,
                 positive = positive, phi = phi,
                 sigma = sigma, alpha = alpha,
                 f = f, psi = psi, z = z, theta = theta, a = a,
                 plot_cell = plot_cell,
                 trend = isTRUE(trend), sigma_trend = if (isTRUE(trend)) sigma_trend else NA_real_,
                 alpha_trend = if (isTRUE(trend)) alpha_trend else NA_real_,
                 f_trend = f_trend, tcov = tcov)
  )
}
