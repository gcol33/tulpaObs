# ---------------------------------------------------------------------------
# Pointwise log-likelihood (WAIC / PSIS-LOO) -- gcol33/tulpaObs#26
# ---------------------------------------------------------------------------

# Pointwise log-likelihood [n_draws x n_sites] for an occu_cover() fit: the
# per-site marginal log-likelihood (latent occupancy state integrated out)
# evaluated at each posterior draw. The spatial nested-Laplace fit samples the
# grid-integrated joint (betas + shared field) via the joint substrate; the
# non-spatial Laplace fit reuses the stored coefficient draws (no field). Both
# feed the same per-draw site-likelihood accumulation.
# Per-draw arm coefficients + dispersion + shared-field contributions for an
# occu_cover() fit. Joint path: sample the grid-integrated posterior and
# accumulate each block's field on the occupancy / cover arms. Non-spatial path:
# read the stored coefficient draws (no field). Returns a list consumed by both
# the pointwise log-likelihood and the posterior-mean plug-in.
.tobs_occu_cover_components <- function(object, n.draws = 1000L) {
  model   <- object$model
  pi_list <- model$process_info
  p_occ   <- pi_list[[1L]]$p
  p_det   <- pi_list[[2L]]$p
  p_pos   <- pi_list[[3L]]$p
  n_sites <- model$n_sites

  if (!is.null(.tobs_joint_fit(object))) {
    bundle <- .tobs_joint_draws(object, n = n.draws)
    S      <- bundle$n
    b_occ  <- bundle$b$occ
    b_det  <- bundle$b$det
    b_pos  <- bundle$b$pos
    disp   <- bundle$disp
    # Each site reads its field node (cell) via site_cell; the field draws
    # blk$z carry one column per node, so blk$z[, units] broadcasts the shared
    # cell field across that cell's sites. Identity when no group_var.
    units  <- model$site_cell %||% seq_len(n_sites)
    wfun   <- function(nm) {
      if (!nm %in% names(model$data)) {
        stop("occu_cover WAIC: trend-field weight column '", nm,
             "' is not in the fitted data.", call. = FALSE)
      }
      as.numeric(model$data[[nm]])
    }
    field_occ <- matrix(0, n_sites, S)
    field_pos <- matrix(0, n_sites, S)
    for (blk in bundle$blocks) {
      z_unit <- blk$z[, units, drop = FALSE]   # [S x n_sites]
      w <- if (is.null(blk$weight)) rep(1, n_sites) else wfun(blk$weight)
      field_occ <- field_occ + t(z_unit * blk$amp_occ) * w
      field_pos <- field_pos + t(z_unit * blk$amp_pos) * w
    }
  } else {
    draws <- object$draws
    if (is.null(draws) || !is.matrix(draws)) {
      stop("occu_cover WAIC: the fit carries no posterior draw matrix.",
           call. = FALSE)
    }
    if (!is.null(n.draws) && n.draws < nrow(draws)) {
      draws <- draws[seq_len(as.integer(n.draws)), , drop = FALSE]
    }
    b_occ <- draws[, seq_len(p_occ), drop = FALSE]
    b_det <- draws[, p_occ + seq_len(p_det), drop = FALSE]
    b_pos <- draws[, p_occ + p_det + seq_len(p_pos), drop = FALSE]
    disp  <- exp(draws[, ncol(draws)])         # trailing log-dispersion column
    S     <- nrow(draws)
    fld <- .tobs_occu_cover_v3_field(object, n_sites, S)
    field_occ <- fld$field_occ
    field_pos <- fld$field_pos
  }
  list(b_occ = b_occ, b_det = b_det, b_pos = b_pos, disp = disp,
       field_occ = field_occ, field_pos = field_pos)
}

# Per-cell field draws for the v3 nested-Laplace occu_cover spatial path, which
# stores no joint object to sample but DOES carry the per-cell field posterior in
# `field_table` (z_mean / z_sd) plus the copy coefficient `alpha` (and, for a
# spatially-varying trend, `trend_field_table` / `alpha_trend` weighted by the
# trend covariate). The per-site pointwise log-likelihood for site i depends only
# on the field at that site's own cell, so per-cell marginal Gaussian draws
# N(z_mean, z_sd) give the exact per-observation marginal -- the cross-cell
# covariance never enters a single site's term. Returns the [n_sites x S]
# occupancy / cover field contributions (zero matrices when no field table is
# present, recovering the non-spatial behaviour).
.tobs_occu_cover_v3_field <- function(object, n_sites, S) {
  zero <- matrix(0, n_sites, S)
  ft <- object$field_table
  if (is.null(ft) || is.null(ft$z_mean)) {
    return(list(field_occ = zero, field_pos = zero))
  }
  model    <- object$model
  site_cell <- model$site_cell %||% seq_len(n_sites)
  # Per-cell field draws (one column per draw), mapped to sites by site_cell.
  draw_cell_field <- function(tab) {
    n_cell <- nrow(tab)
    z <- matrix(stats::rnorm(n_cell * S, tab$z_mean, tab$z_sd), n_cell, S)
    z[site_cell, , drop = FALSE]                # [n_sites x S]
  }
  alpha <- object$means[["alpha"]] %||% object$spatial$alpha_mean %||% 0
  f_occ <- draw_cell_field(ft)
  field_occ <- f_occ
  field_pos <- alpha * f_occ

  # Spatially-varying trend field, weighted per cell by its trend covariate.
  tt <- object$trend_field_table
  if (!is.null(tt) && !is.null(tt$z_mean)) {
    alpha_tr <- object$means[["alpha_trend"]] %||% 0
    wname <- object$trend_weight %||% object$trend_weights[[1L]]
    w <- if (!is.null(wname) && wname %in% names(model$data)) {
      as.numeric(model$data[[wname]])
    } else rep(1, n_sites)
    f_tr <- draw_cell_field(tt) * w             # [n_sites x S], row-scaled
    field_occ <- field_occ + f_tr
    field_pos <- field_pos + alpha_tr * f_tr
  }
  list(field_occ = field_occ, field_pos = field_pos)
}

.tobs_ploglik_occu_cover <- function(object, n.draws = 1000L, n.threads = 1L) {
  c0 <- .tobs_occu_cover_components(object, n.draws)
  .occu_cover_ploglik_core(object$model, c0$b_occ, c0$b_det, c0$b_pos,
                           c0$disp, c0$field_occ, c0$field_pos,
                           n_threads = n.threads)
}

# Pointwise log-likelihood at the posterior mean (length n_sites): the per-site
# marginal evaluated at the mean arm coefficients, mean dispersion, and mean
# shared field -- the same core driven by a one-draw mean.
.tobs_occu_cover_loglik_at_mean <- function(object, n.draws = 1000L) {
  c0 <- .tobs_occu_cover_components(object, n.draws)
  as.numeric(.occu_cover_ploglik_core(
    object$model,
    matrix(colMeans(c0$b_occ), nrow = 1L),
    matrix(colMeans(c0$b_det), nrow = 1L),
    matrix(colMeans(c0$b_pos), nrow = 1L),
    mean(c0$disp),
    matrix(rowMeans(c0$field_occ), ncol = 1L),
    matrix(rowMeans(c0$field_pos), ncol = 1L)
  ))
}

# Accumulate the [S x n_sites] pointwise log-likelihood from per-draw arm
# coefficients (`b_occ` / `b_det` / `b_pos`, each [S x p_arm]), per-draw
# dispersion `disp` [S], and per-arm shared-field contributions `field_occ` /
# `field_pos` [n_sites x S] (zero matrices for a non-spatial fit). The detection
# and cover arms split into site-level and visit-level coefficient blocks exactly
# as the fitter packs them.
# Per-draw linear-predictor blocks for an occu_cover() fit: the occupancy
# `eta_psi_all` [n_sites x S] (with the shared field), and the site-level +
# optional visit-level detection / cover blocks the fitter packs. Split out so
# the pointwise log-likelihood, the posterior predictive check, and the PIT all
# build the per-draw psi / p / cover predictors from one source.
.occu_cover_eta_components <- function(model, b_occ, b_det, b_pos,
                                       field_occ, field_pos) {
  p_det_site <- ncol(model$X_det_site)
  p_pos_site <- ncol(model$X_pos_site)
  has_det_visit <- !is.null(model$X_det_visit)
  has_pos_visit <- !is.null(model$X_pos_visit)
  list(
    eta_psi_all = tcrossprod(model$X_occ, b_occ) + field_occ,
    eta_p_site_all = tcrossprod(model$X_det_site,
                                b_det[, seq_len(p_det_site), drop = FALSE]),
    eta_pos_site_all = tcrossprod(model$X_pos_site,
                                  b_pos[, seq_len(p_pos_site), drop = FALSE]) +
                       field_pos,
    eta_p_visit_all = if (has_det_visit) {
      tcrossprod(model$X_det_visit,
                 b_det[, p_det_site + seq_len(ncol(model$X_det_visit)),
                       drop = FALSE])
    } else NULL,
    eta_pos_visit_all = if (has_pos_visit) {
      tcrossprod(model$X_pos_visit,
                 b_pos[, p_pos_site + seq_len(ncol(model$X_pos_visit)),
                       drop = FALSE])
    } else NULL,
    has_det_visit = has_det_visit, has_pos_visit = has_pos_visit
  )
}

# Draw `d`'s occupancy predictor and the [n_sites x max_visits] detection /
# cover linear predictors, folding the visit-level block in site-major order.
.occu_cover_draw_eta <- function(comp, d, n_sites, max_visits) {
  p_eta <- matrix(comp$eta_p_site_all[, d], n_sites, max_visits)
  if (comp$has_det_visit) {
    p_eta <- p_eta + matrix(comp$eta_p_visit_all[, d], n_sites, max_visits,
                            byrow = TRUE)
  }
  ep_mat <- matrix(comp$eta_pos_site_all[, d], n_sites, max_visits)
  if (comp$has_pos_visit) {
    ep_mat <- ep_mat + matrix(comp$eta_pos_visit_all[, d], n_sites, max_visits,
                              byrow = TRUE)
  }
  list(psi_eta = comp$eta_psi_all[, d], p_eta = p_eta, ep_mat = ep_mat)
}

# Compact (ragged) counterpart of .occu_cover_draw_eta: the detection / cover
# predictors are length-V (one per valid visit), built by broadcasting the
# site-level block onto each visit's site and adding the per-visit block --
# never a padded [n_sites x max_visits] matrix.
.occu_cover_draw_eta_ragged <- function(comp, d, site_of_visit) {
  p_eta <- comp$eta_p_site_all[site_of_visit, d]
  if (comp$has_det_visit) p_eta <- p_eta + comp$eta_p_visit_all[, d]
  ep <- comp$eta_pos_site_all[site_of_visit, d]
  if (comp$has_pos_visit) ep <- ep + comp$eta_pos_visit_all[, d]
  list(psi_eta = comp$eta_psi_all[, d], p_eta = p_eta, ep = ep)
}

# Compact (ragged) counterpart of .occu_cover_site_ll: the per-site occu-cover
# marginal from length-V per-visit predictors grouped by `model$site_of_visit`,
# returned as a length-n_sites vector. Algebraically identical to the dense
# version (same MacKenzie mixture, same per-visit Beta/lognormal cover term);
# rowsum() accumulates each site's visits instead of a matrix rowSums. Scoped to
# cover_aggregate = "none" (the only mode the ragged path builds).
.occu_cover_site_ll_ragged <- function(model, psi, p_vec, ep_vec, log_disp) {
  site    <- model$site_of_visit
  y       <- model$y_det_visit
  n_sites <- model$n_sites
  g       <- factor(site, levels = seq_len(n_sites))

  log_p   <- log(p_vec)
  log_1mp <- log(1 - p_vec)
  log_h_det <- ifelse(y == 1L, log_p, log_1mp)

  is_beta <- identical(model$positive %||% "lognormal", "beta")
  pos_ld  <- numeric(length(y))
  det     <- y == 1L
  if (any(det)) {
    pos_ld[det] <- .occu_cover_pos_logdens(model$y_pos_visit[det], ep_vec[det],
                                           exp(log_disp), is_beta)
  }

  sum_log_h  <- as.numeric(rowsum(log_h_det, g))
  sum_log1mp <- as.numeric(rowsum(log_1mp,   g))
  cover_term <- as.numeric(rowsum(pos_ld,    g))
  n_det      <- as.numeric(rowsum(as.numeric(y), g))

  log_psi   <- log(pmax(psi, 1e-300))
  log_1mpsi <- log(pmax(1 - psi, 1e-300))
  det_ll    <- log_psi + sum_log_h + cover_term
  nodet_ll  <- .tobs_logsumexp2(log_psi + sum_log1mp, log_1mpsi)
  ifelse(n_det > 0, det_ll, nodet_ll)
}

# Available system RAM in bytes, or NA if it cannot be read without a hard
# dependency. Linux (the LiSC server) exposes /proc/meminfo MemAvailable; the
# optional `ps` package covers the other platforms. NA -> the caller falls back
# to a fixed default budget.
.occu_cover_free_ram <- function() {
  v <- tryCatch({
    if (file.exists("/proc/meminfo")) {
      ln <- grep("^MemAvailable:", readLines("/proc/meminfo", n = 60L), value = TRUE)
      if (length(ln)) as.numeric(sub("\\D+(\\d+).*", "\\1", ln[1L])) * 1024 else NA_real_
    } else if (requireNamespace("ps", quietly = TRUE)) {
      as.numeric(ps::ps_system_memory()[["avail"]])
    } else NA_real_
  }, error = function(e) NA_real_)
  if (length(v) != 1L || !is.finite(v)) NA_real_ else v
}

# Draw-chunk size for the WAIC pointwise log-likelihood. The heavy transient is
# the two [n_plots x chunk] per-visit eta matrices (detection + cover), so bound
# `2 * n_plots * chunk * 8` bytes to a fraction of free RAM (a 4 GB default when
# RAM cannot be probed). WAIC is a sum over draws, so chunking is exact. Returns
# a chunk in 1..n_draws.
.occu_cover_waic_chunk <- function(n_plots, n_draws, frac = 0.4) {
  free   <- .occu_cover_free_ram()
  budget <- if (is.finite(free)) frac * free else 4e9
  per_draw <- 2 * max(as.numeric(n_plots), 1) * 8
  max(1L, min(as.integer(n_draws), as.integer(budget / per_draw)))
}

# Flatten a dense (padded [n_sites x max_visits]) no-aggregation occu_cover model
# to the ragged one-row-per-valid-visit form the C++ pointwise kernel consumes.
# The dense visit designs are site-major with n_sites * max_visits rows (cell
# (i, v) at row (i - 1) * max_visits + v), which is exactly the column-major
# position of that cell in t(valid); so the valid-cell selector indexes both the
# response and the visit designs. Cells are enumerated site-major, visit
# ascending -- the same order the dense rowSums accumulates -- so the kernel's
# per-site sums are byte-identical to .occu_cover_site_ll on the dense grid.
.occu_cover_dense_ragged <- function(model) {
  mv  <- model$max_visits
  sel <- which(t(model$valid))               # site-major, visit-ascending
  v_idx <- ((sel - 1L) %% mv) + 1L
  s_idx <- ((sel - 1L) %/% mv) + 1L
  cell  <- cbind(s_idx, v_idx)
  list(
    site_of_visit = as.integer(s_idx),
    y_det_visit   = as.integer(model$y[cell]),
    y_pos_visit   = as.numeric(model$y_pos[cell]),
    X_det_visit   = if (!is.null(model$X_det_visit))
                      model$X_det_visit[sel, , drop = FALSE] else NULL,
    X_pos_visit   = if (!is.null(model$X_pos_visit))
                      model$X_pos_visit[sel, , drop = FALSE] else NULL,
    V = length(sel)
  )
}

# `chunk` (draws per block) defaults to a memory-adaptive size; the [n_plots x
# chunk] eta matrices are the WAIC's memory peak, and processing draws in blocks
# keeps that bounded while the returned [S x n_sites] pointwise log-likelihood is
# byte-identical to the unchunked result (each draw's row depends only on that
# draw).
.occu_cover_ploglik_core <- function(model, b_occ, b_det, b_pos, disp,
                                     field_occ, field_pos, chunk = NULL,
                                     n_threads = 1L) {
  n_sites    <- model$n_sites
  max_visits <- model$max_visits
  S  <- nrow(b_occ)
  cl <- .tobs_clamp_eta
  is_ragged <- isTRUE(model$ragged)
  mode <- model$cover_aggregate %||% "none"

  # No-aggregation path (every ragged fit, and the dense grid without cover
  # aggregation): the C++ kernel mirrors .occu_cover_site_ll_ragged draw for
  # draw and parallelises over draws, with no draw-chunking (each draw's
  # per-visit predictors live in thread-private scratch, not the [V x n_draws]
  # transient .occu_cover_waic_chunk bounds). A dense grid is flattened to the
  # same one-row-per-valid-visit form (.occu_cover_dense_ragged), summed in the
  # same visit order as the dense rowSums, so both feed one kernel.
  if (identical(mode, "none")) {
    rg <- if (is_ragged) {
      list(site_of_visit = as.integer(model$site_of_visit),
           y_det_visit   = as.integer(model$y_det_visit),
           y_pos_visit   = as.numeric(model$y_pos_visit),
           X_det_visit   = model$X_det_visit, X_pos_visit = model$X_pos_visit,
           V = length(model$site_of_visit))
    } else .occu_cover_dense_ragged(model)
    empty_v <- function(m) if (is.null(m)) matrix(0, rg$V, 0L) else m
    return(cpp_occu_cover_ploglik_ragged(
      X_occ = model$X_occ, X_det_site = model$X_det_site,
      X_pos_site = model$X_pos_site,
      X_det_visit = empty_v(rg$X_det_visit),
      X_pos_visit = empty_v(rg$X_pos_visit),
      site_of_visit = rg$site_of_visit,
      y_det_visit   = rg$y_det_visit,
      y_pos_visit   = rg$y_pos_visit,
      b_occ = b_occ, b_det = b_det, b_pos = b_pos, disp = disp,
      field_occ = field_occ, field_pos = field_pos,
      is_beta   = identical(model$positive %||% "lognormal", "beta"),
      eta_bound = .TOBS_ETA_BOUND,
      n_threads = max(1L, as.integer(n_threads))))
  }

  # Aggregated (mean / median / latent) paths stay in R, draw-chunked
  # to bound the [n_plots x chunk] per-visit eta transient. Detected-unit cover
  # values are draw-invariant, so resolve them once (gcol33/tulpaObs#34).
  units <- if (identical(mode, "none")) NULL else .occu_cover_unit_cover(model)

  n_plots <- max(
    if (!is.null(model$X_det_visit)) nrow(model$X_det_visit) else 0L,
    if (!is.null(model$X_pos_visit)) nrow(model$X_pos_visit) else 0L,
    n_sites)
  if (is.null(chunk)) chunk <- .occu_cover_waic_chunk(n_plots, S)

  ll <- matrix(0, S, n_sites)
  for (st in seq.int(1L, S, by = chunk)) {
    idx  <- st:min(st + chunk - 1L, S)
    comp <- .occu_cover_eta_components(model,
              b_occ[idx, , drop = FALSE], b_det[idx, , drop = FALSE],
              b_pos[idx, , drop = FALSE],
              field_occ[, idx, drop = FALSE], field_pos[, idx, drop = FALSE])
    for (j in seq_along(idx)) {
      d <- idx[j]
      de    <- .occu_cover_draw_eta(comp, j, n_sites, max_visits)
      psi   <- stats::plogis(cl(de$psi_eta))
      p_mat <- stats::plogis(cl(de$p_eta))
      ll[d, ] <- .occu_cover_site_ll(model, psi, p_mat, de$ep_mat,
                                     log(disp[d]), units = units)
    }
  }
  ll
}


# ---------------------------------------------------------------------------
# Posterior predictive check + PIT (gcol33/tulpaObs#27)
# ---------------------------------------------------------------------------

# Posterior-predictive cover discrepancy for one draw, at the granularity the
# fitter optimised (gcol33/tulpaObs#34). Returns the observed and replicated
# positive-part discrepancy (`obs` / `rep`) the PPC adds to the detection term.
# `ep_mat` is the draw's [n_sites x max_visits] cover predictor; `disp` its
# dispersion (beta precision / lognormal residual SD, or the cover-RE SD under
# "latent"); `units` the per-unit detected covers from .occu_cover_unit_cover().
# Mode branches mirror .occu_cover_cover_term: "none" scores the per-visit cells;
# "mean" / "median" one aggregated cover per detected unit at the unit predictor
# and dispersion the fit held; "latent" the per-unit covers replicated through
# the shared cover RE (u_i ~ N(0, sigma_u^2)) at the within-unit residual disp2.
.occu_cover_ppc_cover <- function(model, ep_mat, disp, units, is_beta, stat_fn,
                                  cl) {
  draw_pos <- function(eta, d) {
    if (is_beta) {
      mu <- stats::plogis(cl(eta))
      stats::rbeta(length(eta), mu * d, (1 - mu) * d)
    } else {
      exp(stats::rnorm(length(eta), eta, d))
    }
  }
  mean_pos <- function(eta, d) {
    if (is_beta) stats::plogis(cl(eta)) else exp(cl(eta) + d^2 / 2)
  }
  mode <- model$cover_aggregate %||% "none"

  if (identical(mode, "none")) {
    pos_mask <- model$valid & (model$y == 1L)
    if (!any(pos_mask)) return(c(obs = 0, rep = 0))
    Epos <- mean_pos(ep_mat, disp)
    yrep <- matrix(draw_pos(as.vector(ep_mat), disp), nrow(ep_mat), ncol(ep_mat))
    return(c(obs = stat_fn(model$y_pos[pos_mask], Epos[pos_mask]),
             rep = stat_fn(yrep[pos_mask],        Epos[pos_mask])))
  }

  ps <- units$pos_site
  if (length(ps) == 0L) return(c(obs = 0, rep = 0))
  eta <- ep_mat[ps, 1L]                     # unit-level cover predictor

  if (identical(mode, "latent")) {
    disp2   <- model$cover_latent_disp2
    m       <- lengths(units$vals)
    unit_of <- rep(seq_along(ps), m)
    v_all   <- unlist(units$vals, use.names = FALSE)
    eta_all <- eta[unit_of]
    u_all   <- stats::rnorm(length(ps), 0, disp)[unit_of]
    e_all   <- mean_pos(eta_all, disp2)
    return(c(obs = stat_fn(v_all, e_all),
             rep = stat_fn(draw_pos(eta_all + u_all, disp2), e_all)))
  }

  aggfun <- if (identical(mode, "median")) stats::median else mean
  yv   <- vapply(units$vals, function(v) as.numeric(aggfun(v)), numeric(1))
  Epos <- mean_pos(eta, disp)
  c(obs = stat_fn(yv, Epos), rep = stat_fn(draw_pos(eta, disp), Epos))
}

# Posterior predictive check for an occu_cover() fit. Per draw, the latent
# occupancy z is sampled from its full conditional given the detection history
# (the spOccupancy construction), detection replicates from Bernoulli(z p), and
# the cover replicate is built at the granularity the fit used (per-visit for
# cover_aggregate = "none", one aggregated cover per detected unit for "mean" /
# "median", and the shared cover-RE marginal for "latent"; gcol33/tulpaObs#34).
# The discrepancy is a Freeman-Tukey (or chi-squared) sum over the detection
# cells plus the positive-part term, returning a Bayesian p-value.
.tobs_ppc_occu_cover <- function(object,
                                 fit.stat = c("freeman-tukey", "chi-squared"),
                                 n.samples = 500) {
  fit.stat <- match.arg(fit.stat)
  model    <- object$model
  positive <- model$positive %||% "lognormal"
  is_beta  <- identical(positive, "beta")
  c0   <- .tobs_occu_cover_components(object, n.samples)
  comp <- .occu_cover_eta_components(model, c0$b_occ, c0$b_det, c0$b_pos,
                                     c0$field_occ, c0$field_pos)
  disp <- c0$disp
  S <- nrow(c0$b_occ)
  n_sites <- model$n_sites; max_visits <- model$max_visits
  y <- model$y; valid <- model$valid
  cl <- .tobs_clamp_eta
  stat_fn <- if (fit.stat == "freeman-tukey") {
    function(o, e) sum((sqrt(o) - sqrt(e))^2, na.rm = TRUE)
  } else {
    function(o, e) sum((o - e)^2 / (e + 1e-10), na.rm = TRUE)
  }
  any_det <- rowSums(y * valid, na.rm = TRUE) > 0
  n_valid <- rowSums(valid)
  mode <- model$cover_aggregate %||% "none"

  # No-aggregation path: the per-draw simulation (latent z, detection replicate,
  # cover replicate) runs in cpp_occu_cover_ppc, which draws from R's RNG stream
  # in the SAME order as the former R loop, so under a fixed seed the discrepancy
  # is byte-identical. Build the per-draw predictors (deterministic) here.
  if (identical(mode, "none")) {
    psi_all <- matrix(0, n_sites, S)
    p_all   <- matrix(0, n_sites, S * max_visits)
    ep_all  <- matrix(0, n_sites, S * max_visits)
    for (s in seq_len(S)) {
      de   <- .occu_cover_draw_eta(comp, s, n_sites, max_visits)
      cols <- (s - 1L) * max_visits + seq_len(max_visits)
      psi_all[, s]   <- stats::plogis(cl(de$psi_eta))
      p_all[, cols]  <- stats::plogis(cl(de$p_eta))
      ep_all[, cols] <- de$ep_mat
    }
    vint <- valid; storage.mode(vint) <- "integer"
    yint <- y;     storage.mode(yint) <- "integer"
    r <- cpp_occu_cover_ppc(psi_all, p_all, ep_all, yint, model$y_pos, vint,
                            as.integer(any_det), as.integer(n_valid), disp,
                            is_beta, identical(fit.stat, "freeman-tukey"))
    return(list(fit.y = r$fit.y, fit.y.rep = r$fit.y.rep,
                bayesian.p = mean(r$fit.y.rep > r$fit.y)))
  }

  # Aggregated (mean / median / latent) cover discrepancy: the detection
  # replicate plus the aggregated / latent cover replicate run in
  # cpp_occu_cover_ppc_agg with matched RNG order (byte-identical). The observed
  # aggregates / detected covers are draw-invariant and gathered once.
  units <- .occu_cover_unit_cover(model)
  psi_all <- matrix(0, n_sites, S)
  p_all   <- matrix(0, n_sites, S * max_visits)
  ep_all  <- matrix(0, n_sites, S * max_visits)
  for (s in seq_len(S)) {
    de   <- .occu_cover_draw_eta(comp, s, n_sites, max_visits)
    cols <- (s - 1L) * max_visits + seq_len(max_visits)
    psi_all[, s]   <- stats::plogis(cl(de$psi_eta))
    p_all[, cols]  <- stats::plogis(cl(de$p_eta))
    ep_all[, cols] <- de$ep_mat
  }
  vint <- valid; storage.mode(vint) <- "integer"
  yint <- y;     storage.mode(yint) <- "integer"
  ps0  <- as.integer(units$pos_site - 1L)
  if (identical(mode, "latent")) {
    mode_code <- 2L
    vals_flat <- as.numeric(unlist(units$vals, use.names = FALSE))
    unit_off  <- as.integer(c(0L, cumsum(lengths(units$vals))))
    yv <- numeric(0); disp2 <- model$cover_latent_disp2 %||% 0
  } else {
    mode_code <- 1L
    aggfun <- if (identical(mode, "median")) stats::median else mean
    yv <- vapply(units$vals, function(v) as.numeric(aggfun(v)), numeric(1))
    vals_flat <- numeric(0); unit_off <- 0L; disp2 <- 0
  }
  r <- cpp_occu_cover_ppc_agg(psi_all, p_all, ep_all, yint, vint,
    as.integer(any_det), as.integer(n_valid), disp, mode_code, ps0,
    as.numeric(yv), vals_flat, unit_off, disp2, is_beta,
    identical(fit.stat, "freeman-tukey"))
  list(fit.y = r$fit.y, fit.y.rep = r$fit.y.rep,
       bayesian.p = mean(r$fit.y.rep > r$fit.y))
}

# Per-draw CDF limits for the occu_cover() per-site detection summary
# (any-detection vs all-zero), marginalized over the latent occupancy state with
# the shared field folded in per site. Returns the [S x n_sites] lower / upper
# CDF limits of the ordered detected / non-detected outcome -- the randomized-PIT
# building block shared by the posterior PIT and the LOO-PIT.
.occu_cover_pit_cdf_limits <- function(object, n.samples) {
  model <- object$model
  c0    <- .tobs_occu_cover_components(object, n.samples)
  valid <- model$valid; y <- model$y
  any_det <- as.integer(rowSums(y * valid, na.rm = TRUE) > 0)
  vint <- valid; storage.mode(vint) <- "integer"
  empty_v <- function(m) if (is.null(m))
    matrix(0, model$n_sites * model$max_visits, 0L) else m
  # The per-draw detection-summary CDF limits are deterministic; the former R
  # loop now runs in cpp_occu_cover_cdf_limits, parallel over draws.
  cpp_occu_cover_cdf_limits(model$X_occ, model$X_det_site,
                            empty_v(model$X_det_visit), c0$b_occ, c0$b_det,
                            c0$field_occ, vint, any_det, 1L)
}

# Randomized PIT for an occu_cover() fit, on the per-site detection summary
# (any-detection vs all-zero) marginalized over the latent occupancy state, with
# the shared field projected per site. The detected / non-detected outcome is the
# ordered event; the left and right CDF limits feed the engine's randomized PIT.
.tobs_pit_occu_cover <- function(object, n.samples = 250) {
  lim <- .occu_cover_pit_cdf_limits(object, n.samples)
  tulpa::tulpa_pit(lim$cdf_upper, cdf_lower = lim$cdf_lower)
}

# Leave-one-out randomized PIT for an occu_cover() fit -- the INLA `cpo$pit`
# analogue. Per site, the draws are reweighted by their PSIS leave-one-out
# importance weights (so the predictive distribution excludes that site's own
# contribution), then the LOO-weighted CDF limits feed the randomized PIT. The
# loglik matrix supplying the weights is the field-folded pointwise log-
# likelihood, so the LOO-PIT inherits the full-model predictor. A site whose PSIS
# weights are degenerate (k unavailable) falls back to the equal-weight posterior
# CDF for that site.
.tobs_loo_pit_occu_cover <- function(object, n.samples = 1000L, ll = NULL) {
  lim <- .occu_cover_pit_cdf_limits(object, n.samples)
  if (is.null(ll)) ll <- .tobs_ploglik_occu_cover(object, n.samples)
  .tobs_loo_pit_from_limits(ll, lim$cdf_lower, lim$cdf_upper)
}

# LOO-weighted randomized PIT from a pointwise log-likelihood matrix `ll`
# [S x N] and the per-draw CDF limits `Fl` / `Fu` [S x N] of the ordered outcome.
# For each observation the PSIS leave-one-out weights w_is (from tulpa_psis on
# -ll[, i]) reweight the per-draw CDF limits to the LOO predictive limits
# Fl_loo_i = sum_s w_is Fl[s, i], Fu_loo_i likewise, then a single uniform draw
# placed between them gives the randomized LOO-PIT. The single source of truth
# behind every family's LOO-PIT so the construction matches INLA's cpo$pit.
.tobs_loo_pit_from_limits <- function(ll, Fl, Fu) {
  # Per-observation PSIS leave-one-out weighting of the CDF limits + a uniform
  # jitter, batched in tulpa's cpp_psis_loo_pit (PSIS columns parallel, the runif
  # in index order), so it is byte-identical to the former per-column R loop.
  tail_len <- getFromNamespace(".psis_tail_len", "tulpa")(nrow(ll), NULL)
  getFromNamespace("cpp_psis_loo_pit", "tulpa")(
    ll, Fl, Fu, as.integer(tail_len), 1L)
}


