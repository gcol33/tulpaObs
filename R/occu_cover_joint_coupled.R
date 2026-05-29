# =============================================================================
# occu_cover_joint_coupled.R - joint nested-Laplace path for occu_cover()
# driven by the cell-coupling registry (gcol33/tulpa#32 Layer B.2 consumer).
#
# Routes through tulpa_nested_laplace_joint(cell_coupling =
# "occu_cover_lognormal") with a 3-arm responses list:
#   * psi arm: one row per cell, no observed data; spec writes its grad
#              + neg-hess from the cell-level occupancy mixture
#   * p   arm: one row per valid visit; field_coef = 0 (no shared field)
#   * pos arm: one row per valid visit; field_coef = list(name = "alpha")
#              so the cover arm's field amplitude is alpha * sigma on the
#              outer (sigma, alpha) grid
#
# The cell-coupling spec (registered from tulpaObs's .onLoad) writes every
# closed-form derivative -- the per-arm scatter in the joint engine is
# skipped for `coupled = TRUE` arms and the per-cell density (det / nodet
# branches) drives the inner Newton.
#
# Spatial path only. Non-spatial occu_cover stays on method = "laplace";
# the spatial v3 nested-Laplace path (.tobs_fit_occu_cover_nested) stays
# as the default under method = "nested_laplace". This file is reached
# via control$engine = "joint_coupled" -- same pattern as v2_joint.
# =============================================================================


# Convert a dense adjacency matrix to the CSR (row_ptr, col_idx, n_neighbors)
# layout the joint engine's ICAR block consumes. Indices are 0-based, sorted
# within each row. Errors on isolated nodes (no neighbours), matching
# .occu_cover_icar_Q.
.occu_cover_adj_to_csr <- function(adj) {
  n <- nrow(adj)
  nbr <- lapply(seq_len(n), function(i) sort(which(adj[i, ] != 0)) - 1L)
  n_neighbors <- vapply(nbr, length, integer(1))
  if (any(n_neighbors == 0L)) {
    isolated <- which(n_neighbors == 0L)
    stop(sprintf("ICAR graph has %d isolated node(s): %s. ",
                 length(isolated),
                 paste(utils::head(isolated, 5L), collapse = ", ")),
         "Drop them or connect them before fitting.", call. = FALSE)
  }
  list(n_spatial_units = n,
       adj_row_ptr     = as.integer(c(0L, cumsum(n_neighbors))),
       adj_col_idx     = as.integer(unlist(nbr)),
       n_neighbors     = as.integer(n_neighbors))
}


# Assemble the three-arm `responses` list the joint engine consumes.
# Compacts visit-level rows to valid visits only (matching the v3
# nested-Laplace mask). Visit ordering within a cell follows site-major
# layout: valid visits 1..max_visits for cell c are emitted in increasing
# visit index, so the p and pos arms see identical (cell, visit) pairings
# in identical row order -- the spec reads them positionally inside each
# cell.
.occu_cover_build_joint_coupled_arms <- function(model, sigma_pos_init,
                                                  alpha_grid) {
  n_cells    <- model$n_sites
  max_visits <- model$max_visits

  valid_flat  <- as.logical(t(model$valid))   # site-major: site 1 visits 1..J, ...
  y_flat      <- as.numeric(t(model$y))
  y_pos_flat  <- as.numeric(t(model$y_pos))
  cell_flat   <- rep(seq_len(n_cells), each = max_visits)
  visit_flat  <- rep(seq_len(max_visits), times = n_cells)

  keep             <- which(valid_flat)
  n_visits_valid   <- length(keep)
  y_det_visit      <- as.integer(y_flat[keep])
  y_pos_visit      <- y_pos_flat[keep]
  cell_idx_visit   <- as.integer(cell_flat[keep])

  if (n_visits_valid == 0L) {
    stop("occu_cover joint_coupled: no valid visits in the data.",
         call. = FALSE)
  }

  # Site-level + visit-level fixed-effect design on each observation arm.
  X_p_site <- model$X_det_site
  X_p <- X_p_site[cell_idx_visit, , drop = FALSE]
  if (!is.null(model$X_det_visit)) {
    X_p <- cbind(X_p, model$X_det_visit[keep, , drop = FALSE])
  }
  X_pos_site <- model$X_pos_site
  X_pos <- X_pos_site[cell_idx_visit, , drop = FALSE]
  if (!is.null(model$X_pos_visit)) {
    X_pos <- cbind(X_pos, model$X_pos_visit[keep, , drop = FALSE])
  }

  # psi arm: one row per cell. y / n_trials / family are placeholders -- the
  # per-obs scatter is skipped for coupled = TRUE and the cell-coupling spec
  # writes every derivative from the cell-level occupancy mixture.
  arm_psi <- list(
    y            = rep(0, n_cells),
    n_trials     = rep(0L, n_cells),
    X            = model$X_occ,
    spatial_idx  = seq_len(n_cells),
    family       = "binomial",
    phi          = 1.0,
    coupled      = TRUE,
    cell_obs_map = seq_len(n_cells)
  )

  # p arm: one row per valid visit. field_coef = 0 -> no field scatter, so
  # spatial_idx is ignored; pass 0L to satisfy the length check.
  arm_p <- list(
    y            = as.numeric(y_det_visit),
    n_trials     = rep(1L, n_visits_valid),
    X            = X_p,
    spatial_idx  = rep(0L, n_visits_valid),
    family       = "binomial",
    phi          = 1.0,
    field_coef   = 0,
    coupled      = TRUE,
    cell_obs_map = cell_idx_visit
  )

  # pos arm: one row per valid visit. field_coef carries the alpha axis;
  # spatial_idx maps each visit to its cell so alpha * sigma * z[cell]
  # contributes the cover-arm field. phi is the lognormal SD on the log
  # scale; the spec reads y_cell.phi(2) = sigma_pos.
  arm_pos <- list(
    y            = y_pos_visit,
    n_trials     = rep(1L, n_visits_valid),
    X            = X_pos,
    spatial_idx  = cell_idx_visit,
    family       = "lognormal",
    phi          = sigma_pos_init,
    field_coef   = list(name = "alpha", grid = alpha_grid),
    coupled      = TRUE,
    cell_obs_map = cell_idx_visit
  )

  list(psi = arm_psi, p = arm_p, pos = arm_pos)
}


# Per-arm coef names for unpacking joint modes into a (means, sds) table.
.occu_cover_joint_coupled_coef_names <- function(model) {
  pi_list <- model$process_info
  list(
    psi = pi_list[[1L]]$coef_names,
    p   = pi_list[[2L]]$coef_names,
    pos = pi_list[[3L]]$coef_names
  )
}


# Joint-coupled fitter. Calls tulpa_nested_laplace_joint() with the
# 3-arm responses and the occu_cover_lognormal cell-coupling spec, then
# unpacks the integrated posterior into a tobs_fit shaped to match
# .tobs_fit_occu_cover_nested's output (so methods.R / generic accessors
# work without per-engine branching).
#
# Hyperparam grid: outer axes are (sigma, alpha), Cartesian product
# defaulting to the engine's 5-point sigma grid and the engine's
# 6-point alpha grid (incl. 0). sigma_pos is fixed pre-fit at the
# empirical SD of log(y_pos) at detected visits; pass
# control$phi.grid.pos to integrate over it as a phi_grid axis on the
# pos arm.
.tobs_fit_occu_cover_joint_coupled <- function(model, adj,
                                                priors    = NULL,
                                                max.iter  = 200L,
                                                tol       = 1e-6,
                                                verbose   = TRUE,
                                                sigma.beta = 5,
                                                ...) {
  if (!identical(model$positive, "lognormal")) {
    stop("occu_cover() joint_coupled engine supports positive = ",
         "\"lognormal\" only (the occu_cover_lognormal cell-coupling spec is ",
         "lognormal-specific). Beta positive arm is a follow-up.",
         call. = FALSE)
  }

  pi_list <- model$process_info
  n_cells <- model$n_sites
  if (nrow(adj) != n_cells) {
    stop(sprintf("Spatial graph has %d nodes but data has %d cells.",
                 nrow(adj), n_cells), call. = FALSE)
  }

  dots <- list(...)

  # Pre-fit sigma_pos at the empirical log-SD of detected cover values.
  # Matches the v3 nested-Laplace warm start for log_disp.
  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  sigma_pos_init <- if (length(pos_vals) > 0L) {
    max(stats::sd(log(pos_vals)), 0.05) + 0.05
  } else {
    0.4
  }

  alpha_grid <- dots$alpha.grid %||%
                c(0, exp(seq(log(0.1), log(3), length.out = 5)))
  sigma_grid <- dots$sigma.grid %||%
                exp(seq(log(0.1), log(3), length.out = 5))

  responses <- .occu_cover_build_joint_coupled_arms(
    model           = model,
    sigma_pos_init  = sigma_pos_init,
    alpha_grid      = alpha_grid
  )

  csr <- .occu_cover_adj_to_csr(adj)
  prior_block <- list(
    type            = "icar",
    n_spatial_units = csr$n_spatial_units,
    adj_row_ptr     = csr$adj_row_ptr,
    adj_col_idx     = csr$adj_col_idx,
    n_neighbors     = csr$n_neighbors,
    sigma_grid      = sigma_grid
  )

  # Optional sigma_pos integration via a phi_grid axis on the pos arm.
  phi_grid_pos <- dots$phi.grid.pos
  phi_grid_arg <- if (!is.null(phi_grid_pos))
                    list(pos = as.numeric(phi_grid_pos))
                  else NULL

  fit <- tulpa::tulpa_nested_laplace_joint(
    responses     = responses,
    prior         = prior_block,
    phi_grid      = phi_grid_arg,
    cell_coupling = "occu_cover_lognormal",
    control = list(
      max_iter  = as.integer(max.iter),
      tol       = as.numeric(tol),
      n_threads = as.integer(dots$n.threads %||% 1L),
      store_Q   = TRUE,
      # Adaptive-grid refinement defaults OFF for the first joint_coupled
      # release. Edge-detection compares per-cell log_marginal across the
      # outer grid; a degenerate hyperpoint (e.g. sigma at the lower bound
      # combined with alpha = 0) can land non-finite log_marginal that
      # propagates NA into the refine helper's threshold comparison. Pass
      # `control$adaptive.grid = TRUE` to opt in once the per-cell log_marginal
      # stays finite across the whole grid.
      adaptive_grid             = dots$adaptive.grid             %||% FALSE,
      adaptive_grid_edge_thresh = dots$adaptive.grid.edge.thresh %||% 0.02,
      adaptive_grid_max_passes  = dots$adaptive.grid.max.passes  %||% 1L
    )
  )

  # Unpack per-arm posterior means + SDs from the joint modes.
  # arm_layout$beta_start[k] is the 0-based offset of arm k's betas in the
  # latent vector; arm_layout$p[k] is the count.
  layout <- fit$arm_layout
  p_psi  <- layout$p[1L]
  p_p    <- layout$p[2L]
  p_pos  <- layout$p[3L]
  bpsi_idx <- layout$beta_start[1L] + seq_len(p_psi)
  bp_idx   <- layout$beta_start[2L] + seq_len(p_p)
  bpos_idx <- layout$beta_start[3L] + seq_len(p_pos)

  # Drop outer-grid cells whose inner Newton did not converge (NaN
  # log_marginal). Their modes hold NaN that would poison every weighted
  # sum below if left in. Zero-mass cells stay represented in fit$weights;
  # we just route around them when computing posterior moments.
  ok_cells <- which(is.finite(fit$log_marginal))
  if (length(ok_cells) == 0L) {
    stop("occu_cover joint_coupled: inner Newton failed at every grid cell. ",
         "Bump control$max.iter or tighten control$tol.", call. = FALSE)
  }
  if (length(ok_cells) < length(fit$log_marginal)) {
    n_bad <- length(fit$log_marginal) - length(ok_cells)
    warning(sprintf(
      "occu_cover joint_coupled: dropping %d / %d outer-grid cell(s) ",
      n_bad, length(fit$log_marginal)),
      "whose inner Newton did not converge.", call. = FALSE)
  }
  w_raw <- exp(fit$log_marginal[ok_cells] - max(fit$log_marginal[ok_cells]))
  w     <- w_raw / sum(w_raw)
  modes <- fit$modes[ok_cells, , drop = FALSE]
  beta_psi_m  <- as.numeric(crossprod(w, modes[, bpsi_idx, drop = FALSE]))
  beta_p_m    <- as.numeric(crossprod(w, modes[, bp_idx,   drop = FALSE]))
  beta_pos_m  <- as.numeric(crossprod(w, modes[, bpos_idx, drop = FALSE]))

  # Per-coefficient SD via var-of-means across the outer grid (Rubin's
  # between-imputation variance). Underestimates the full posterior SD by
  # the within-cell Laplace contribution (the law-of-total-variance inner
  # term); the inner term requires diag(Q^-1) at each grid cell and Q has
  # the ICAR field's sum-to-zero null direction so a naive Cholesky-then-
  # solve blows up on the fixed-effect intercepts that share variance with
  # the field. Inner-term integration is a follow-up (constrained Sigma-
  # column solve like sla_cover_hurdle_joint.R).
  outer_sd <- function(modes_block, mean_vec) {
    outer_part <- as.numeric(crossprod(w, modes_block^2)) - mean_vec^2
    sqrt(pmax(outer_part, 0))
  }
  sd_psi <- outer_sd(modes[, bpsi_idx, drop = FALSE], beta_psi_m)
  sd_p   <- outer_sd(modes[, bp_idx,   drop = FALSE], beta_p_m)
  sd_pos <- outer_sd(modes[, bpos_idx, drop = FALSE], beta_pos_m)

  means <- c(beta_psi_m, beta_p_m, beta_pos_m)
  sds   <- c(sd_psi,    sd_p,     sd_pos)

  par_names <- c(
    paste0("psi_", pi_list[[1L]]$coef_names),
    paste0("p_",   pi_list[[2L]]$coef_names),
    paste0("pos_", pi_list[[3L]]$coef_names)
  )

  # Hyperparams: surface sigma, alpha (and optionally sigma_pos when on the
  # outer grid) from the joint posterior moments. Recompute on the filtered
  # grid -- fit$theta_mean uses the engine's full weights, which are NaN
  # when any cell's log_marginal is non-finite.
  tg_full   <- fit$theta_grid
  tg_ok     <- tg_full[ok_cells, , drop = FALSE]
  tg_names  <- colnames(tg_full)
  hyper_means <- numeric(0)
  hyper_sds   <- numeric(0)
  hyper_names <- character(0)
  pick <- function(name) {
    j <- match(name, tg_names)
    if (is.na(j)) return(invisible(NULL))
    vals <- as.numeric(tg_ok[, j])
    m <- sum(w * vals)
    v <- sum(w * vals^2) - m^2
    hyper_means[[name]] <<- m
    hyper_sds  [[name]] <<- sqrt(max(v, 0))
    hyper_names <<- c(hyper_names, name)
  }
  for (nm in c("sigma", "alpha", "phi_pos")) pick(nm)
  if (length(hyper_names) > 0L) {
    means <- c(means, unlist(hyper_means)[hyper_names])
    sds   <- c(sds,   unlist(hyper_sds)[hyper_names])
    par_names <- c(par_names, hyper_names)
  }

  names(means) <- par_names
  names(sds)   <- par_names
  V_dummy <- diag(sds^2)
  dimnames(V_dummy) <- list(par_names, par_names)

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V_dummy)
  colnames(draws) <- par_names

  # Field summary at the posterior-mean grid mode (simple weighted average
  # across grid cells, like nested fitter's z_final).
  field_at_cell <- as.numeric(crossprod(w, modes[, layout$phi_start +
                                                   seq_len(n_cells),
                                                 drop = FALSE]))
  field_at_cell <- field_at_cell - mean(field_at_cell)
  field_sd      <- sqrt(pmax(
    as.numeric(crossprod(w, modes[, layout$phi_start + seq_len(n_cells),
                                  drop = FALSE]^2)) - field_at_cell^2,
    0))

  field_table <- data.frame(
    cell    = seq_len(n_cells),
    z_mean  = field_at_cell,
    z_sd    = field_sd,
    z_lower = field_at_cell - 1.96 * field_sd,
    z_upper = field_at_cell + 1.96 * field_sd
  )

  structure(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V_dummy,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(0, n_draws),
    log_lik      = sum(w * fit$log_marginal[ok_cells]),
    N            = sum(model$valid),
    accept_prob  = rep(1, n_draws),
    divergent    = rep(0L, n_draws),
    treedepth    = rep(0L, n_draws),
    epsilon      = NA_real_,
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = list(type = "icar", graph = adj,
                        sigma_mean = unname(hyper_means["sigma"]),
                        alpha_mean = unname(hyper_means["alpha"])),
    spatial_field = field_at_cell,
    field_table  = field_table,
    method       = "joint_coupled",
    positive     = model$positive,
    joint_fit    = fit,
    convergence  = list(converged = TRUE, n_iter = NA_integer_)
  ), class = c("tobs_fit", "tulpa_fit"))
}
