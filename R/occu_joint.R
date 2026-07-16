# =============================================================================
# occu_joint.R - single-arm occupancy nested-Laplace path for the
# standalone occu() varying-coefficient (SVC) spatial bar (gcol33/tulpaObs#81).
#
# This is occu_cover()'s joint direct-grid engine (occu_cover_joint.R)
# with the cover (positive) arm removed: the occupancy (psi) + detection (p)
# arms only, the shared cell-indexed intercept + trend areal field on the
# occupancy arm, no coupling alpha. The field hyperparameters integrate on a
# direct outer (sigma, sigma_trend) grid, so the fit cannot oscillate the way
# the EM fixed-point path (.tobs_em_nested_laplace) does at EVA scale.
#
# Routes through tulpa_nested_laplace_joint(cell_coupling = "occu_only") with a
# 2-arm responses list:
#   * psi arm: one row per occupancy unit (site), no observed data; the
#              cell-coupling spec writes its grad + neg-hess from the per-site
#              occupancy mixture
#   * p   arm: one row per valid visit; spatial_idx is the 0-sentinel so the
#              field never enters the detection predictor
#
# The occu_only cell-coupling spec (registered from tulpaObs's .onLoad) reuses
# the SAME occu_det_psi_p_block / occu_nodet_block the occu_cover specs use, so
# the occupancy / detection derivatives are single-source.
#
# Field assembly (CSR, demeaning, multi-block trend copy), the law-of-total-
# covariance posterior moments, the rmvn draws and the hyperparameter
# marginalisation reuse the occu_cover_joint helpers; only the arm
# roster (2 arms, no pos / no alpha) and the per-arm coefficient unpacking
# differ.
# =============================================================================


# Assemble the two-arm `responses` list (psi + p) the joint engine consumes for
# a standalone occu() fit. Compacts visit-level rows to valid visits only. The
# psi arm carries one row per site (cell_obs_map = site, the occupancy unit the
# mixture groups by); the p arm carries one row per valid visit. The field rides
# the occupancy arm only -- the p arm's spatial_idx is the 0-sentinel that the
# joint engine reads as "no field on this row" (l_b > 0 gate in
# nested_laplace_joint_multi.h's INDEXED_SINGLE scatter), the occupancy-only
# analogue of the occu_cover p arm's field_coef = 0.
.occu_joint_arms <- function(model, n_cells, site_cell) {
  n_sites    <- model$n_sites
  max_visits <- model$max_visits

  # The standalone occu() `single` model encodes NA visits as y < 0; the valid
  # mask is y >= 0. X_occ / X_det are X_processes[[1]] / [[2]] (autoscaled
  # upstream by .tobs_fit_model); X_det_visit is the optional visit-level design.
  y_mat       <- model$y
  valid_mat   <- y_mat >= 0
  X_occ       <- model$X_processes[[1L]]
  X_det_site  <- model$X_processes[[2L]]
  X_det_visit <- model$X_det_visit

  valid_flat <- as.logical(t(valid_mat))           # site-major: site 1 visits 1..J
  y_flat     <- as.numeric(t(y_mat))
  site_flat  <- rep(seq_len(n_sites), each = max_visits)

  keep           <- which(valid_flat)
  n_visits_valid <- length(keep)
  if (n_visits_valid == 0L) {
    stop("occu joint: no valid visits in the data.", call. = FALSE)
  }
  y_det_visit   <- as.integer(y_flat[keep])
  site_of_visit <- as.integer(site_flat[keep])
  cell_of_visit <- as.integer(site_cell[site_of_visit])

  # Detection design on each valid visit: the site-level row broadcast to the
  # visit, plus any visit-level columns (sliced to valid visits, site-major).
  X_p <- X_det_site[site_of_visit, , drop = FALSE]
  if (!is.null(X_det_visit)) {
    X_p <- cbind(X_p, X_det_visit[keep, , drop = FALSE])
  }

  arm_psi <- list(
    y            = rep(0, n_sites),
    n_trials     = rep(0L, n_sites),
    X            = X_occ,
    spatial_idx  = as.integer(site_cell),
    family       = "binomial",
    phi          = 1.0,
    coupled      = TRUE,
    cell_obs_map = seq_len(n_sites)
  )
  arm_p <- list(
    y            = as.numeric(y_det_visit),
    n_trials     = rep(1L, n_visits_valid),
    X            = X_p,
    spatial_idx  = rep(0L, n_visits_valid),
    family       = "binomial",
    phi          = 1.0,
    coupled      = TRUE,
    cell_obs_map = site_of_visit
  )

  list(responses      = list(psi = arm_psi, p = arm_p),
       site_of_visit  = site_of_visit,
       cell_of_visit  = cell_of_visit,
       n_visits_valid = n_visits_valid)
}


# Resolve per-arm weakly-informative fixed-effect priors for the two coupled
# arms, returning list(psi=, p=) of list(mean, prec) (NULL per arm -> the weak
# engine default). The detection (p) intercept prior keeps the coupled occupancy
# mixture off the psi = 1 boundary at weak detection -- the same load-bearing
# default the occu_cover coupled arms carry (the logit-scale score vanishes as
# psi -> 1, so an unpenalised intercept runs away). `priors = FALSE` / "none"
# disables both. A supplied occu_priors() / list overrides the matching arm(s).
.occu_joint_arm_priors <- function(priors, responses) {
  if (identical(priors, FALSE) || identical(priors, "none")) {
    return(list(psi = NULL, p = NULL))
  }
  to_prec <- function(pr) {
    if (is.null(pr) || length(pr$sd) == 0L) return(NULL)
    list(mean = as.numeric(pr$mean), prec = pmax(1 / pr$sd^2, 1e-4))
  }
  occ_spec <- .resolve_occu_priors(priors)        # NULL / list / occu_priors
  list(
    psi = to_prec(.prior_for_submodel(occ_spec, "psi", colnames(responses$psi$X))),
    p   = to_prec(.prior_for_submodel(occ_spec, "p",   colnames(responses$p$X)))
  )
}


# Single-arm occupancy joint-coupled fitter. Calls tulpa_nested_laplace_joint()
# with the 2-arm responses and the occu_only cell-coupling spec, integrating the
# field hyperparameter(s) on the outer (sigma [, sigma_trend]) grid, then unpacks
# the integrated posterior into a tobs_fit. The shaping mirrors
# .occu_cover_jc_postprocess (cover arm dropped) so methods.R / the joint
# substrate read it without per-engine branching.
#
# `model` is the standalone occu() `single` model (autoscaled upstream by
# .tobs_fit_model). `fields` is the intercept-first ordered list of `tobs_spatial`
# specs from .tobs_resolve_occu_spatial_fields(): the unweighted intercept field
# first, then any weighted SVC (trend) fields, all on one areal graph.
.tobs_fit_occu_joint <- function(model, fields,
                                         priors  = NULL,
                                         max.iter = 200L,
                                         tol      = 1e-6,
                                         verbose  = TRUE,
                                         ...) {
  if (!inherits(model, "tobs_model") ||
      !identical(model$model_type, "single")) {
    stop("occu joint engine fits a single-season occu() model.",
         call. = FALSE)
  }
  adj <- fields[[1L]]$graph

  # Field nodes (cells) and occupancy units (sites). With an areal group_var many
  # sites share one cell field node; site_cell maps each site to its node. The
  # fields carry the node-index column as `group_var`; resolve it against the
  # model data the single model retained.
  n_cells <- nrow(adj)
  n_sites <- model$n_sites
  gv <- fields[[1L]]$group_var
  if (!is.null(gv)) {
    if (is.null(model$data) || !gv %in% names(model$data)) {
      stop(sprintf(paste0("occu joint: group_var '%s' is not a column of the ",
                          "model data."), gv), call. = FALSE)
    }
    site_cell <- as.integer(model$data[[gv]])
    if (length(site_cell) != n_sites || anyNA(site_cell) ||
        min(site_cell) < 1L || max(site_cell) > n_cells) {
      stop(sprintf(paste0(
        "occu joint: group_var '%s' must be an integer cell index in ",
        "1..%d, one per site (%d sites)."), gv, n_cells, n_sites),
        call. = FALSE)
    }
  } else {
    if (n_cells != n_sites) {
      stop(sprintf(paste0(
        "occu joint: the areal field has %d nodes but the model has %d ",
        "sites; map sites to cells with group_var on the spatial term, or pass ",
        "one node per site."), n_cells, n_sites), call. = FALSE)
    }
    site_cell <- seq_len(n_sites)
  }
  model$site_cell <- site_cell
  model$n_cells   <- n_cells

  dots <- list(...)

  # Coupled trend (SVC) fields: each is a per-cell-weighted areal field that
  # contributes weight_i * sigma_trend * z[cell_i] on the occupancy predictor.
  # They arrive as the weighted entries of `fields` (each carrying a resolved
  # per-site `$weight` and a `$weight_label`). The intercept field is fields[[1]];
  # with at least one trend field the fit takes the multi-block path.
  weighted_fields <- fields[-1L]
  coupled_trends <- lapply(weighted_fields, function(f) {
    w <- as.numeric(f$weight)
    if (length(w) != n_sites || any(!is.finite(w))) {
      stop(sprintf(paste0(
        "occu joint: trend field weight must be a finite per-site ",
        "numeric vector of length %d."), n_sites), call. = FALSE)
    }
    list(weight = w, weight_label = f$weight_label %||% "trend")
  })
  n_trend   <- length(coupled_trends)
  has_trend <- n_trend > 0L

  # The single-arm occupancy field has no cover arm to copy onto, so each ICAR
  # block grids on its precision tau directly (the non-copy multi-block path,
  # axis `b<k>.tau`). The user-facing `sigma.grid` is the field SD (the units the
  # occu_cover / cover paths grid on); convert it to tau = 1 / sigma^2 so the SD
  # knob is honoured, and the postprocess reports sigma = 1 / sqrt(tau) as a
  # derived quantity marginalized over the grid. The 4-point default keeps the
  # base outer grid at 4^n_fields cells (16 for an intercept + one trend field),
  # comparable to the occu_cover deliverable's 4-point sigma grid; adaptive
  # refinement densifies it where the posterior piles up at an edge.
  sigma_grid <- dots$sigma.grid %||% exp(seq(log(0.15), log(3), length.out = 4))
  tau_grid   <- 1.0 / (as.numeric(sigma_grid)^2)

  arms_out      <- .occu_joint_arms(model, n_cells, site_cell)
  responses     <- arms_out$responses
  cell_of_visit <- arms_out$cell_of_visit
  n_v           <- arms_out$n_visits_valid

  arm_priors <- .occu_joint_arm_priors(priors, responses)
  for (nm in c("psi", "p")) {
    ap <- arm_priors[[nm]]
    if (!is.null(ap)) {
      responses[[nm]]$beta_prior_mean <- ap$mean
      responses[[nm]]$beta_prior_prec <- ap$prec
    }
  }

  csr <- .occu_cover_adj_to_csr(adj)
  icar_template <- function(extra = list()) {
    c(list(
        type            = "icar",
        n_spatial_units = csr$n_spatial_units,
        adj_row_ptr     = csr$adj_row_ptr,
        adj_col_idx     = csr$adj_col_idx,
        n_neighbors     = csr$n_neighbors,
        tau_grid        = tau_grid
      ), extra)
  }

  if (has_trend) {
    # Multi-block path: the intercept ICAR block plus one ICAR block per coupled
    # trend field, all on the same graph. The field rides the occupancy (psi) arm
    # only; the p arm carries the 0-sentinel node index, so its field weight is
    # irrelevant. There is no cover arm and no copy -- each block is an
    # independent areal field whose amplitude (sigma / sigma_trend) integrates on
    # the outer grid. Per-block svc_weight injects the per-row field weight on the
    # psi (per-site) arm; the p arm weight is a placeholder of length n_v.
    spatial_idx_arms <- list(as.integer(site_cell), rep(0L, n_v))
    make_block <- function(weight_site) {
      w_psi <- if (is.null(weight_site)) rep(1.0, n_sites)
               else as.numeric(weight_site)
      icar_template(list(
        spatial_idx = spatial_idx_arms,
        svc_weight  = list(w_psi, rep(1.0, n_v))
      ))
    }
    prior_arg <- c(
      list(make_block(NULL)),
      lapply(coupled_trends, function(tf) make_block(tf$weight))
    )
  } else {
    # Single shared intercept field on the occupancy arm; the multi-block driver
    # with an explicit per-arm spatial_idx (the p arm's 0-sentinel excludes the
    # field there). Routing through the multi-block driver -- rather than the
    # single-block backend, which derives the per-arm field map from a pos arm's
    # field_coef this fit does not have -- keeps the assembly identical to the
    # trend path with n_trend = 0.
    prior_arg <- list(icar_template(list(
      spatial_idx = list(as.integer(site_cell), rep(0L, n_v)))))
  }

  fit_call <- list(
    responses     = responses,
    prior         = prior_arg,
    cell_coupling = "occu_only",
    control = list(
      max_iter        = as.integer(max.iter),
      tol             = as.numeric(tol),
      n_threads       = as.integer(dots$n.threads %||% 1L),
      store_Q         = TRUE,
      hessian         = dots$hessian %||% "lm",
      inner_refresh   = as.integer(dots$inner.refresh %||% 1L),
      # The outer grid runs across `n.threads.outer` threads (each cell an
      # independent inner Laplace solve). At EVA scale the per-cell inner solve on
      # the full areal field dominates, so the outer grid is the parallelism that
      # makes the fit land in occu_cover-like wall time; default it to the machine
      # parallelism (leaving a few cores) rather than serial, since a standalone
      # occu() SVC fit on real data is the EVA-scale workload this path targets.
      n_threads_outer = as.integer(
        dots$n.threads.outer %||% max(1L, parallel::detectCores() - 4L)),
      force_sparse    = isTRUE(dots$force.sparse),
      adaptive_grid             = dots$adaptive.grid             %||% TRUE,
      adaptive_grid_edge_thresh = dots$adaptive.grid.edge.thresh %||% 0.02,
      adaptive_grid_max_passes  = dots$adaptive.grid.max.passes  %||% 1L,
      var_of_means_consistency  = dots$var.of.means.consistency  %||% TRUE,
      var_of_means_tolerance    = dots$var.of.means.tolerance    %||% 0.7,
      # Outer Pareto-k-hat accuracy diagnostic defaults OFF on this path. At EVA
      # scale it importance-samples the hyperparameter posterior with `k_samples`
      # extra inner solves on the full areal field, which dominates the runtime
      # (the engine note records ~50x at EVA scale) -- and the single-arm
      # occupancy field posterior is well-behaved, so the diagnostic is an opt-in,
      # not a default cost. Set control$diagnose.k = TRUE to compute it.
      diagnose_k = dots$diagnose.k %||% FALSE,
      # diagnose.draws is the precision knob (k.samples is the legacy alias); the
      # outer Pareto-k is scored ONCE over this many importance draws.
      k_samples = as.integer(dots$diagnose.draws %||% dots$k.samples %||% 500L),
      # Bootstrap outer Pareto-k uncertainty (gcol33/tulpa#127): SE / 95% CI /
      # band_confident from resampling the raw log-ratios (NO new solves). Raise
      # diagnose.draws, not k.bootstrap, for a tighter k. k.tail.points (NULL =
      # automatic PSIS rule) is an expert control; k.conf.bands the band boundaries.
      k_bootstrap   = as.integer(dots$k.bootstrap %||% 1000L),
      k_tail_points = if (is.null(dots$k.tail.points)) NULL else as.integer(dots$k.tail.points),
      k_conf_bands  = dots$k.conf.bands %||% c(0.5, 0.7),
      checkpoint = dots$checkpoint,
      integration = dots$integration,
      progress          = dots[["progress"]] %||% TRUE,
      progress.every    = dots$progress.every,
      progress.throttle = dots$progress.throttle,
      progress.file     = dots$progress.file
    )
  )

  ctx <- list(adj = adj, pi_list = model$process_info, n_cells = n_cells,
              has_trend = has_trend, n_trend = n_trend,
              coupled_trends = coupled_trends, model = model,
              n_threads = as.integer(
                dots$n.threads.outer %||% max(1L, parallel::detectCores() - 4L)))

  fit <- do.call(tulpa::tulpa_nested_laplace_joint, fit_call)
  .occu_jc_postprocess(fit, ctx)
}


# Post-process an occu single-arm joint-coupled engine fit into a tobs_fit. The
# shaping mirrors .occu_cover_jc_postprocess with the pos (cover) arm removed:
# two arms (psi, p), no alpha copy, the field amplitude(s) reported as
# sigma / sigma_trend and the cell-indexed fields as spatial_field / trend_field.
.occu_jc_postprocess <- function(fit, ctx) {
  adj            <- ctx$adj
  pi_list        <- ctx$pi_list
  n_cells        <- ctx$n_cells
  has_trend      <- ctx$has_trend
  n_trend        <- ctx$n_trend
  coupled_trends <- ctx$coupled_trends
  model          <- ctx$model

  layout <- fit$arm_layout
  p_psi  <- layout$p[1L]
  p_p    <- layout$p[2L]
  bpsi_idx <- layout$beta_start[1L] + seq_len(p_psi)
  bp_idx   <- layout$beta_start[2L] + seq_len(p_p)

  ok_cells <- which(is.finite(fit$log_marginal))
  if (length(ok_cells) == 0L) {
    stop("occu joint: inner Newton failed at every grid cell. ",
         "Bump control$max.iter or tighten control$tol.", call. = FALSE)
  }
  if (length(ok_cells) < length(fit$log_marginal)) {
    n_bad <- length(fit$log_marginal) - length(ok_cells)
    warning(sprintf(paste0(
      "occu joint: dropping %d / %d outer-grid cell(s) ",
      "whose inner Newton did not converge."),
      n_bad, length(fit$log_marginal)), call. = FALSE)
  }
  w_raw <- exp(fit$log_marginal[ok_cells] - max(fit$log_marginal[ok_cells]))
  w     <- w_raw / sum(w_raw)

  if (!any(is.finite(fit$weights) & fit$weights > 0)) {
    w_full <- numeric(length(fit$log_marginal))
    w_full[ok_cells] <- w
    fit$weights <- w_full
  }

  modes <- fit$modes[ok_cells, , drop = FALSE]
  beta_psi_m <- as.numeric(crossprod(w, modes[, bpsi_idx, drop = FALSE]))
  beta_p_m   <- as.numeric(crossprod(w, modes[, bp_idx,   drop = FALSE]))

  p_beta        <- p_psi + p_p
  field_starts0 <- layout$field_starts %||% layout$phi_start
  n_fields      <- length(field_starts0)
  field_idx     <- as.integer(unlist(lapply(field_starts0,
                                            function(s0) s0 + seq_len(n_cells))))
  idx_joint <- c(bpsi_idx, bp_idx, field_idx)
  blocks    <- .joint_inner_vcov_block(fit, idx_joint, n_dense = p_psi + p_p,
                                       n_threads = ctx$n_threads)

  if (is.null(blocks)) {
    # Older tulpa without stored per-grid Q: marginal-only diagonal (var-of-means
    # plus the diagonal inner-Laplace variance), no betas+field cross-covariance.
    beta_idx_all <- c(bpsi_idx, bp_idx)
    inner_var    <- .joint_inner_var(fit, beta_idx_all)
    total_var <- function(modes_block, mean_vec, iv_block) {
      vom <- as.numeric(crossprod(w, modes_block^2)) - mean_vec^2
      mov <- if (is.null(iv_block)) {
        0
      } else {
        iv_k <- iv_block[ok_cells, , drop = FALSE]
        iv_k[!is.finite(iv_k)] <- 0
        as.numeric(crossprod(w, iv_k))
      }
      pmax(vom + mov, 0)
    }
    iv_psi <- if (is.null(inner_var)) NULL
              else inner_var[, seq_len(p_psi), drop = FALSE]
    iv_p   <- if (is.null(inner_var)) NULL
              else inner_var[, p_psi + seq_len(p_p), drop = FALSE]
    sds_beta <- c(
      sqrt(total_var(modes[, bpsi_idx, drop = FALSE], beta_psi_m, iv_psi)),
      sqrt(total_var(modes[, bp_idx,   drop = FALSE], beta_p_m,   iv_p))
    )
    beta_block <- diag(sds_beta^2, nrow = p_beta)

    field_modes    <- modes[, field_idx, drop = FALSE]
    field_at_cell  <- as.numeric(crossprod(w, field_modes))
    field_var      <- as.numeric(crossprod(w, field_modes^2)) - field_at_cell^2
    field_demeaned <- .occu_cover_demean_fields(field_at_cell, n_cells, n_fields)
    Vj <- NULL
  } else {
    p_joint     <- length(idx_joint)
    modes_joint <- modes[, idx_joint, drop = FALSE]
    mbar_joint  <- as.numeric(crossprod(w, modes_joint))
    Vj <- matrix(0, p_joint, p_joint)
    for (kk in seq_along(ok_cells)) {
      dk     <- modes_joint[kk, ] - mbar_joint
      Ck     <- blocks[[ ok_cells[kk] ]]
      within <- if (is.null(Ck) || anyNA(Ck)) matrix(0, p_joint, p_joint)
                else as.matrix(Ck)
      Vj <- Vj + w[kk] * (within + tcrossprod(dk))
    }
    Vj <- (Vj + t(Vj)) / 2
    diag_Vj    <- diag(Vj)
    sds_beta   <- sqrt(pmax(diag_Vj[seq_len(p_beta)], 0))
    beta_block <- Vj[seq_len(p_beta), seq_len(p_beta), drop = FALSE]

    n_field_cols   <- n_fields * n_cells
    field_at_cell  <- mbar_joint[p_beta + seq_len(n_field_cols)]
    field_var      <- diag_Vj[p_beta + seq_len(n_field_cols)]
    field_demeaned <- .occu_cover_demean_fields(field_at_cell, n_cells, n_fields)
  }
  field_sd <- sqrt(pmax(field_var, 0))

  sd_psi <- sds_beta[seq_len(p_psi)]
  sd_p   <- sds_beta[p_psi + seq_len(p_p)]
  means  <- c(beta_psi_m, beta_p_m)
  sds    <- c(sd_psi,     sd_p)
  # The p arm carries the site-level detection design followed by any visit-level
  # columns (X_p = cbind(X_det_site, X_det_visit)), so the p block is wider than
  # `pi_list[[2]]$coef_names` (site-level only) when a visit detection design is
  # present. Name the visit tail `p_visit_<name>`, matching the laplace / NUTS
  # output layout; the per-process autoscaler only scales the site-level head
  # (it never sees X_det_visit), so the unscale that walks `pi_list[[2]]$p` site
  # coefficients is unaffected by the wider name vector.
  psi_names <- paste0("psi_", pi_list[[1L]]$coef_names)
  p_names   <- paste0("p_", pi_list[[2L]]$coef_names)
  dv_names  <- model$det_visit_names
  if (!is.null(dv_names) && length(dv_names) > 0L) {
    p_names <- c(p_names, paste0("p_visit_", dv_names))
  }
  par_names <- c(psi_names, p_names)

  # Hyperparameters: surface sigma (and sigma_trend per trend field) from the
  # joint posterior moments on the filtered grid (the marginalize-derived-
  # quantities rule -- a grid-weighted summary, not a plug-in at the mode).
  tg_full   <- fit$theta_grid
  tg_ok     <- tg_full[ok_cells, , drop = FALSE]
  tg_names  <- colnames(tg_full)
  hyper_means <- numeric(0)
  hyper_sds   <- numeric(0)
  hyper_vals  <- list()
  hyper_names <- character(0)
  # Each ICAR block grids on its precision tau (axis `b<k>.tau`); the reported
  # field SD is the derived quantity sigma = 1 / sqrt(tau), summarized as a
  # grid-weighted mean / SD of the per-cell sigma values (the marginalize-derived-
  # quantities rule -- not 1 / sqrt(mean(tau))). The per-cell sigma values feed
  # the law-of-total-covariance block below so sigma's covariance with the betas
  # and the other field SDs is carried, not dropped.
  pick_sigma <- function(public, tau_col) {
    j <- match(tau_col, tg_names)
    if (is.na(j)) return(invisible(NULL))
    tau  <- as.numeric(tg_ok[, j])
    vals <- 1.0 / sqrt(tau)
    m <- sum(w * vals)
    v <- sum(w * vals^2) - m^2
    hyper_means[[public]] <<- m
    hyper_sds  [[public]] <<- sqrt(max(v, 0))
    hyper_vals [[public]] <<- vals
    hyper_names <<- c(hyper_names, public)
  }
  pick_sigma("sigma", "b1.tau")
  for (j in seq_len(n_trend)) {
    suffix <- if (n_trend == 1L) "" else as.character(j)
    pick_sigma(paste0("sigma_trend", suffix), sprintf("b%d.tau", j + 1L))
  }
  if (length(hyper_names) > 0L) {
    means <- c(means, unlist(hyper_means)[hyper_names])
    sds   <- c(sds,   unlist(hyper_sds)[hyper_names])
    par_names <- c(par_names, hyper_names)
  }
  names(means) <- par_names
  names(sds)   <- par_names

  # Parameter-surface vcov by the law of total covariance over the outer grid.
  # The betas carry within + between (beta_block); the hyperparameters ARE the
  # grid coordinates, so within a cell their variance / cross-covariance with the
  # betas is zero, leaving only the between (mode-deviation) term, computed
  # jointly so the cross- and hyper-hyper covariance are retained.
  n_par <- length(means)
  V <- matrix(0, n_par, n_par)
  V[seq_len(p_beta), seq_len(p_beta)] <- beta_block
  if (length(hyper_names) > 0L) {
    hyper_idx <- p_beta + seq_along(hyper_names)
    H    <- do.call(cbind, hyper_vals[hyper_names])
    H_dm <- sweep(H, 2L, unlist(hyper_means)[hyper_names], "-")
    if (!is.null(Vj)) {
      beta_modes <- modes[, c(bpsi_idx, bp_idx), drop = FALSE]
      B_dm   <- sweep(beta_modes, 2L, means[seq_len(p_beta)], "-")
      cross  <- crossprod(B_dm * w, H_dm)
      hyhy   <- crossprod(H_dm * w, H_dm)
      hyhy   <- (hyhy + t(hyhy)) / 2
      V[seq_len(p_beta), hyper_idx] <- cross
      V[hyper_idx, seq_len(p_beta)] <- t(cross)
      V[hyper_idx, hyper_idx]       <- hyhy
    } else {
      diag(V)[hyper_idx] <- sds[hyper_idx]^2
    }
  }
  dimnames(V) <- list(par_names, par_names)

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  field_block <- function(b) {
    idx <- (b - 1L) * n_cells + seq_len(n_cells)
    list(mean = field_demeaned[idx], sd = field_sd[idx])
  }
  field_z_table <- function(blk) {
    data.frame(cell = seq_len(n_cells), z_mean = blk$mean, z_sd = blk$sd,
               z_lower = blk$mean - 1.96 * blk$sd,
               z_upper = blk$mean + 1.96 * blk$sd)
  }
  fblocks <- lapply(seq_len(n_fields), field_block)

  field_intercept <- fblocks[[1L]]$mean
  field_table     <- field_z_table(fblocks[[1L]])

  trend_blocks <- if (n_fields >= 2L) fblocks[-1L] else list()
  trend_means  <- lapply(trend_blocks, function(b) b$mean)
  trend_tables <- lapply(trend_blocks, field_z_table)
  trend_labels <- vapply(coupled_trends, function(tf) tf$weight_label,
                         character(1))
  if (length(trend_labels) == length(trend_means)) {
    names(trend_means)  <- trend_labels
    names(trend_tables) <- trend_labels
  }
  field_trend       <- if (length(trend_means))  trend_means[[1L]]  else NULL
  trend_field_table <- if (length(trend_tables)) trend_tables[[1L]] else NULL

  # Joint betas+field posterior for field-aware predict (the occupancy psi at a
  # new cell marginalizes the betas + field). Fields are stacked in block order
  # (intercept then trend fields), demeaned to the sum-to-zero convention.
  field_par_names <- unlist(c(
    list(paste0("field_", seq_len(n_cells))),
    lapply(seq_len(n_fields - 1L), function(j) {
      suffix <- if (n_fields - 1L == 1L) "" else as.character(j)
      paste0("trend_field", suffix, "_", seq_len(n_cells))
    })
  ))
  joint_par_names <- c(psi_names, p_names, field_par_names)
  joint_means <- c(beta_psi_m, beta_p_m, field_demeaned)
  names(joint_means) <- joint_par_names
  if (!is.null(Vj)) dimnames(Vj) <- list(joint_par_names, joint_par_names)

  log_lik_val <- sum(w * fit$log_marginal[ok_cells])

  spatial_summary <- list(
    type = "icar", graph = adj,
    sigma_mean = unname(hyper_means["sigma"]),
    alpha_mean = NULL,
    sigma_trend_mean = if (has_trend)
      unname(hyper_means[if (n_trend == 1L) "sigma_trend" else "sigma_trend1"])
      else NULL,
    alpha_trend_mean = NULL)

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(log_lik_val, n_draws),
    log_lik      = log_lik_val,
    N            = sum(model$y >= 0)),
    .tobs_na_nuts_diagnostics(n_draws),
    .tobs_promote_pareto_k(fit),
    list(
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = spatial_summary,
    spatial_field = field_intercept,
    trend_field   = field_trend,
    trend_fields  = if (length(trend_means))  trend_means  else NULL,
    field_table  = field_table,
    trend_field_table  = trend_field_table,
    trend_field_tables = if (length(trend_tables)) trend_tables else NULL,
    trend_weight  = if (has_trend) trend_labels[[1L]] else NULL,
    trend_weights = if (has_trend) trend_labels        else NULL,
    joint_par_names = joint_par_names,
    joint_means     = joint_means,
    joint_vcov      = Vj,
    method       = "joint",
    occu_only_joint = TRUE,
    joint_fit    = fit,
    convergence  = list(converged = TRUE, n_iter = NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}
