# =============================================================================
# occu_multiscale_cover_joint_coupled.R - joint nested-Laplace path for the
# three-level occupancy + cover hurdle (gcol33/tulpaObs#29).
#
# Four-arm generalisation of occu_cover's joint_coupled path: a cell-level
# occupancy gate (psi), a plot-level availability gate (theta), per-visit
# detection (p) and the cover hurdle (pos). Drives
# tulpa_nested_laplace_joint(cell_coupling = "occu_multiscale_cover_*") with:
#   * psi   arm: one row per cell,           field_coef = 1 (shared field)
#   * theta arm: one row per plot,           field_coef = 0
#   * p     arm: one row per valid visit,    field_coef = 0
#   * pos   arm: one row per valid visit,    field_coef = list(name = "alpha")
#
# The cell-coupling spec (registered per-fit, carrying this fit's per-cell plot
# structure) writes every closed-form derivative; the per-arm scatter is
# skipped for coupled = TRUE arms and the nested per-cell density (branch A /
# branch B) drives the inner Newton. Both z (over cells) and a (over plots)
# marginalize in closed form, so the joint marginal log-likelihood is exact.
#
# Single shared intercept field only (the (sigma, alpha) grid); spatially
# varying trend fields are an occu_cover feature not yet carried here. Reuses
# occu_cover's CSR / field-demean / inner-vcov / rmvn helpers.
# =============================================================================


# Assemble the four-arm `responses` list. Visit rows are compacted to valid
# visits and laid PLOT-MAJOR within each cell (plot 1's visits, plot 2's, ...),
# the order build_cell_rows_from_arms reproduces from input row order so the
# spec's per-cell `cell_plot_sizes` partition aligns with arm 2/3's flat rows.
.occu_ms_cover_build_arms <- function(model, sigma_pos_init, alpha_grid,
                                      positive = "lognormal") {
  n_cells   <- model$n_cells
  n_plots   <- model$n_plots
  J         <- model$max_visits
  plot_cell <- model$plot_cell

  valid_cols <- lapply(seq_len(n_plots), function(i) which(model$valid[i, ]))
  nvis       <- vapply(valid_cols, length, integer(1))
  if (any(nvis == 0L)) {
    stop("occu_multiscale_cover: every plot must have >= 1 valid visit; ",
         sprintf("plot(s) %s have none.",
                 paste(utils::head(which(nvis == 0L), 5L), collapse = ", ")),
         call. = FALSE)
  }

  # Plot-major visit emission (plot 1, then plot 2, ...).
  vis_plot   <- rep(seq_len(n_plots), nvis)
  cell_of_v  <- plot_cell[vis_plot]
  flat_idx   <- unlist(lapply(seq_len(n_plots),
                              function(i) (i - 1L) * J + valid_cols[[i]]))
  y_det_v    <- unlist(lapply(seq_len(n_plots),
                              function(i) model$y[i, valid_cols[[i]]]))
  y_pos_v    <- unlist(lapply(seq_len(n_plots),
                              function(i) model$y_pos[i, valid_cols[[i]]]))
  n_v        <- length(vis_plot)

  # Per-arm visit-level designs (plot-level block + optional visit-varying).
  X_p   <- model$X_p_site[vis_plot, , drop = FALSE]
  if (!is.null(model$X_p_visit))   X_p   <- cbind(X_p,   model$X_p_visit[flat_idx, , drop = FALSE])
  X_pos <- model$X_pos_site[vis_plot, , drop = FALSE]
  if (!is.null(model$X_pos_visit)) X_pos <- cbind(X_pos, model$X_pos_visit[flat_idx, , drop = FALSE])

  arm_psi <- list(
    y            = rep(0, n_cells),
    n_trials     = rep(0L, n_cells),
    X            = model$X_psi,
    spatial_idx  = seq_len(n_cells),
    family       = "binomial",
    phi          = 1.0,
    coupled      = TRUE,
    cell_obs_map = seq_len(n_cells)
  )
  arm_theta <- list(
    y            = rep(0, n_plots),
    n_trials     = rep(0L, n_plots),
    X            = model$X_theta,
    spatial_idx  = rep(0L, n_plots),
    family       = "binomial",
    phi          = 1.0,
    field_coef   = 0,
    coupled      = TRUE,
    cell_obs_map = as.integer(plot_cell)
  )
  arm_p <- list(
    y            = as.numeric(y_det_v),
    n_trials     = rep(1L, n_v),
    X            = X_p,
    spatial_idx  = rep(0L, n_v),
    family       = "binomial",
    phi          = 1.0,
    field_coef   = 0,
    coupled      = TRUE,
    cell_obs_map = as.integer(cell_of_v)
  )
  arm_pos <- list(
    y            = y_pos_v,
    n_trials     = rep(1L, n_v),
    X            = X_pos,
    spatial_idx  = as.integer(cell_of_v),
    family       = positive,
    phi          = sigma_pos_init,
    field_coef   = list(name = "alpha", grid = alpha_grid),
    coupled      = TRUE,
    cell_obs_map = as.integer(cell_of_v)
  )

  # Per-cell plot structure for the spec: plots grouped cell-major, ascending
  # plot index within cell (order(plot_cell) is a stable radix sort on the
  # integer cell id), with each plot's valid-visit count.
  n_plots_per_cell <- tabulate(plot_cell, nbins = n_cells)
  plot_sizes_flat  <- nvis[order(plot_cell)]

  list(responses        = list(psi = arm_psi, theta = arm_theta,
                               p = arm_p, pos = arm_pos),
       n_plots_per_cell = as.integer(n_plots_per_cell),
       plot_sizes_flat  = as.integer(plot_sizes_flat),
       n_visits_valid   = n_v)
}


# Per-arm weakly-informative fixed-effect priors for the four coupled arms.
# psi / theta / p carry the occu_priors() defaults (the detection-arm intercept
# prior in particular keeps the coupled occupancy mixture off the psi = 1
# boundary); pos carries cover_priors() only when supplied. theta reuses the
# detection-arm prior shape (a logit-scale availability gate).
.occu_ms_cover_arm_priors <- function(priors, responses) {
  if (identical(priors, FALSE) || identical(priors, "none")) {
    return(list(psi = NULL, theta = NULL, p = NULL, pos = NULL))
  }
  to_prec <- function(pr) {
    if (is.null(pr) || length(pr$sd) == 0L) return(NULL)
    list(mean = as.numeric(pr$mean), prec = pmax(1 / pr$sd^2, 1e-4))
  }
  occ_spec <- if (inherits(priors, "occu_priors")) priors
              else if (inherits(priors, "cover_priors")) occu_priors()
              else .resolve_occu_priors(priors)
  cover_spec <- if (inherits(priors, "cover_priors")) priors else NULL

  list(
    psi   = to_prec(.prior_for_submodel(occ_spec, "psi", colnames(responses$psi$X))),
    theta = to_prec(.prior_for_submodel(occ_spec, "p",   colnames(responses$theta$X))),
    p     = to_prec(.prior_for_submodel(occ_spec, "p",   colnames(responses$p$X))),
    pos   = to_prec(.cover_arm_prior(cover_spec, "pos",  colnames(responses$pos$X)))
  )
}


# Joint-coupled fitter for the three-level family. Registers the stateful
# multiscale cell-coupling spec carrying this fit's per-cell plot structure,
# calls tulpa_nested_laplace_joint() with the four arms, then unpacks the
# integrated posterior into a tobs_fit shaped like the occu_cover joint fit.
.tobs_fit_occu_multiscale_cover_joint_coupled <- function(model, fields,
                                                          priors    = NULL,
                                                          max.iter  = 200L,
                                                          tol       = 1e-6,
                                                          verbose   = TRUE,
                                                          ...) {
  adj <- fields[[1L]]$graph
  is_beta <- identical(model$positive, "beta")
  is_lnrm <- identical(model$positive, "lognormal")
  if (!is_beta && !is_lnrm) {
    stop("occu_multiscale_cover() supports positive = \"lognormal\" or ",
         "\"beta\".", call. = FALSE)
  }
  if (length(fields) > 1L) {
    stop("occu_multiscale_cover() couples a single shared intercept field; ",
         "spatially varying trend fields are not yet supported.", call. = FALSE)
  }

  pi_list <- model$process_info
  n_cells <- nrow(adj)
  if (model$n_cells != n_cells) {
    stop(sprintf(paste0(
      "occu_multiscale_cover: %d cells in the model but the graph has %d ",
      "nodes."), model$n_cells, n_cells), call. = FALSE)
  }

  dots <- list(...)

  # Pre-fit the pos-arm dispersion at the empirical sample value at detected
  # visits (lognormal: SD of log y_pos; beta: moment-matched precision).
  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  phi_pos_init <- if (is_beta) {
    if (length(pos_vals) >= 2L) {
      mu_hat  <- mean(pos_vals); var_hat <- max(stats::var(pos_vals), 1e-6)
      max((mu_hat * (1 - mu_hat)) / var_hat - 1, 1)
    } else 10
  } else {
    if (length(pos_vals) > 0L) max(stats::sd(log(pos_vals)), 0.05) + 0.05 else 0.4
  }

  alpha_grid <- dots$alpha.grid %||% c(0, exp(seq(log(0.1), log(3), length.out = 5)))
  sigma_grid <- dots$sigma.grid %||% exp(seq(log(0.1), log(3), length.out = 5))

  arms_out  <- .occu_ms_cover_build_arms(
    model = model, sigma_pos_init = phi_pos_init,
    alpha_grid = alpha_grid, positive = model$positive)
  responses <- arms_out$responses

  # Attach per-arm fixed-effect priors.
  arm_priors <- .occu_ms_cover_arm_priors(priors, responses)
  for (nm in c("psi", "theta", "p", "pos")) {
    ap <- arm_priors[[nm]]
    if (!is.null(ap)) {
      responses[[nm]]$beta_prior_mean <- ap$mean
      responses[[nm]]$beta_prior_prec <- ap$prec
    }
  }

  # Register the stateful spec carrying this fit's per-cell plot structure
  # (last-writer-wins under a fixed name; the previous fit's spec is released).
  spec_name <- cpp_register_occu_multiscale_cover_coupling(
    positive         = model$positive,
    n_plots_per_cell = arms_out$n_plots_per_cell,
    plot_sizes_flat  = arms_out$plot_sizes_flat)

  csr <- .occu_cover_adj_to_csr(adj)
  prior_arg <- list(
    type            = "icar",
    n_spatial_units = csr$n_spatial_units,
    adj_row_ptr     = csr$adj_row_ptr,
    adj_col_idx     = csr$adj_col_idx,
    n_neighbors     = csr$n_neighbors,
    sigma_grid      = sigma_grid
  )

  phi_grid_pos <- dots$phi.grid.pos
  phi_grid_arg <- if (!is.null(phi_grid_pos)) list(pos = as.numeric(phi_grid_pos)) else NULL

  fit_call <- list(
    responses     = responses,
    prior         = prior_arg,
    phi_grid      = phi_grid_arg,
    cell_coupling = spec_name,
    control = list(
      max_iter        = as.integer(max.iter),
      tol             = as.numeric(tol),
      n_threads       = as.integer(dots$n.threads %||% 1L),
      store_Q         = TRUE,
      hessian         = dots$hessian %||% (if (is_beta) "fisher" else "lm"),
      inner_refresh   = as.integer(dots$inner.refresh %||% 1L),
      n_threads_outer = as.integer(dots$n.threads.outer %||% 1L),
      force_sparse    = isTRUE(dots$force.sparse),
      adaptive_grid             = dots$adaptive.grid             %||% TRUE,
      adaptive_grid_edge_thresh = dots$adaptive.grid.edge.thresh %||% 0.02,
      adaptive_grid_max_passes  = dots$adaptive.grid.max.passes  %||% 1L,
      diagnose_k = dots$diagnose.k %||% TRUE,
      k_samples  = as.integer(dots$k.samples %||% 200L),
      checkpoint = dots$checkpoint
    )
  )

  fit <- do.call(tulpa::tulpa_nested_laplace_joint, fit_call)

  # Unpack per-arm posterior means + SDs from the joint modes (4 arms).
  layout <- fit$arm_layout
  p_psi   <- layout$p[1L]; p_theta <- layout$p[2L]
  p_p     <- layout$p[3L]; p_pos   <- layout$p[4L]
  bpsi_idx   <- layout$beta_start[1L] + seq_len(p_psi)
  btheta_idx <- layout$beta_start[2L] + seq_len(p_theta)
  bp_idx     <- layout$beta_start[3L] + seq_len(p_p)
  bpos_idx   <- layout$beta_start[4L] + seq_len(p_pos)

  ok_cells <- which(is.finite(fit$log_marginal))
  if (length(ok_cells) == 0L) {
    stop("occu_multiscale_cover: inner Newton failed at every grid cell. ",
         "Bump control$max.iter or tighten control$tol.", call. = FALSE)
  }
  if (length(ok_cells) < length(fit$log_marginal)) {
    n_bad <- length(fit$log_marginal) - length(ok_cells)
    warning(sprintf(
      "occu_multiscale_cover: dropping %d / %d outer-grid cell(s) whose inner ",
      n_bad, length(fit$log_marginal)),
      "Newton did not converge.", call. = FALSE)
  }
  w_raw <- exp(fit$log_marginal[ok_cells] - max(fit$log_marginal[ok_cells]))
  w     <- w_raw / sum(w_raw)
  modes <- fit$modes[ok_cells, , drop = FALSE]
  beta_psi_m   <- as.numeric(crossprod(w, modes[, bpsi_idx,   drop = FALSE]))
  beta_theta_m <- as.numeric(crossprod(w, modes[, btheta_idx, drop = FALSE]))
  beta_p_m     <- as.numeric(crossprod(w, modes[, bp_idx,     drop = FALSE]))
  beta_pos_m   <- as.numeric(crossprod(w, modes[, bpos_idx,   drop = FALSE]))

  p_beta        <- p_psi + p_theta + p_p + p_pos
  field_starts0 <- layout$field_starts %||% layout$phi_start
  field_idx     <- as.integer(field_starts0[1L] + seq_len(n_cells))
  idx_joint     <- c(bpsi_idx, btheta_idx, bp_idx, bpos_idx, field_idx)
  blocks        <- .joint_inner_vcov_block(fit, idx_joint)

  if (is.null(blocks)) {
    # Older tulpa without stored per-grid Q: marginal-only diagonal fallback.
    beta_idx_all <- c(bpsi_idx, btheta_idx, bp_idx, bpos_idx)
    inner_var    <- .joint_inner_var(fit, beta_idx_all)
    total_var <- function(modes_block, mean_vec, iv_block) {
      vom <- as.numeric(crossprod(w, modes_block^2)) - mean_vec^2
      mov <- if (is.null(iv_block)) 0 else {
        iv_k <- iv_block[ok_cells, , drop = FALSE]
        iv_k[!is.finite(iv_k)] <- 0
        as.numeric(crossprod(w, iv_k))
      }
      pmax(vom + mov, 0)
    }
    cuts <- cumsum(c(0L, p_psi, p_theta, p_p, p_pos))
    iv_of <- function(k) if (is.null(inner_var)) NULL
                         else inner_var[, (cuts[k] + 1L):cuts[k + 1L], drop = FALSE]
    sds_beta <- c(
      sqrt(total_var(modes[, bpsi_idx,   drop = FALSE], beta_psi_m,   iv_of(1L))),
      sqrt(total_var(modes[, btheta_idx, drop = FALSE], beta_theta_m, iv_of(2L))),
      sqrt(total_var(modes[, bp_idx,     drop = FALSE], beta_p_m,     iv_of(3L))),
      sqrt(total_var(modes[, bpos_idx,   drop = FALSE], beta_pos_m,   iv_of(4L)))
    )
    beta_block    <- diag(sds_beta^2, nrow = p_beta)
    field_modes   <- modes[, field_idx, drop = FALSE]
    field_at_cell <- as.numeric(crossprod(w, field_modes))
    field_var     <- as.numeric(crossprod(w, field_modes^2)) - field_at_cell^2
    field_demeaned <- field_at_cell - mean(field_at_cell)
    Vj <- NULL
  } else {
    p_joint     <- length(idx_joint)
    modes_joint <- modes[, idx_joint, drop = FALSE]
    mbar_joint  <- as.numeric(crossprod(w, modes_joint))
    Vj <- matrix(0, p_joint, p_joint)
    for (kk in seq_along(ok_cells)) {
      dk     <- modes_joint[kk, ] - mbar_joint
      Ck     <- blocks[[ ok_cells[kk] ]]
      within <- if (is.null(Ck) || anyNA(Ck)) matrix(0, p_joint, p_joint) else as.matrix(Ck)
      Vj <- Vj + w[kk] * (within + tcrossprod(dk))
    }
    Vj <- (Vj + t(Vj)) / 2
    diag_Vj    <- diag(Vj)
    sds_beta   <- sqrt(pmax(diag_Vj[seq_len(p_beta)], 0))
    beta_block <- Vj[seq_len(p_beta), seq_len(p_beta), drop = FALSE]
    field_at_cell  <- mbar_joint[p_beta + seq_len(n_cells)]
    field_var      <- diag_Vj[p_beta + seq_len(n_cells)]
    field_demeaned <- field_at_cell - mean(field_at_cell)
  }
  field_sd <- sqrt(pmax(field_var, 0))

  cuts <- cumsum(c(0L, p_psi, p_theta, p_p, p_pos))
  sd_psi   <- sds_beta[(cuts[1] + 1L):cuts[2]]
  sd_theta <- sds_beta[(cuts[2] + 1L):cuts[3]]
  sd_p     <- sds_beta[(cuts[3] + 1L):cuts[4]]
  sd_pos   <- sds_beta[(cuts[4] + 1L):cuts[5]]

  means <- c(beta_psi_m, beta_theta_m, beta_p_m, beta_pos_m)
  sds   <- c(sd_psi,     sd_theta,     sd_p,     sd_pos)
  par_names <- c(
    paste0("psi_",   pi_list[[1L]]$coef_names),
    paste0("theta_", pi_list[[2L]]$coef_names),
    paste0("p_",     pi_list[[3L]]$coef_names),
    paste0("pos_",   pi_list[[4L]]$coef_names)
  )

  # Hyperparameters from the joint posterior moments.
  tg_full  <- fit$theta_grid
  tg_ok    <- tg_full[ok_cells, , drop = FALSE]
  tg_names <- colnames(tg_full)
  hyper_means <- numeric(0); hyper_sds <- numeric(0); hyper_names <- character(0)
  pick <- function(name) {
    j <- match(name, tg_names)
    if (is.na(j)) return(invisible(NULL))
    vals <- as.numeric(tg_ok[, j])
    m <- sum(w * vals); v <- sum(w * vals^2) - m^2
    hyper_means[[name]] <<- m; hyper_sds[[name]] <<- sqrt(max(v, 0))
    hyper_names <<- c(hyper_names, name)
  }
  for (nm in c("sigma", "alpha", "phi_pos")) pick(nm)
  if (length(hyper_names) > 0L) {
    means <- c(means, unlist(hyper_means)[hyper_names])
    sds   <- c(sds,   unlist(hyper_sds)[hyper_names])
    par_names <- c(par_names, hyper_names)
  }
  names(means) <- par_names; names(sds) <- par_names

  n_par <- length(means)
  V <- matrix(0, n_par, n_par)
  V[seq_len(p_beta), seq_len(p_beta)] <- beta_block
  if (length(hyper_names) > 0L) {
    hyper_idx <- p_beta + seq_along(hyper_names)
    diag(V)[hyper_idx] <- sds[hyper_idx]^2
  }
  dimnames(V) <- list(par_names, par_names)

  n_draws <- 1000L
  draws <- .occu_cover_rmvn(n_draws, means, V)
  colnames(draws) <- par_names

  field_table <- data.frame(
    cell = seq_len(n_cells), z_mean = field_demeaned, z_sd = field_sd,
    z_lower = field_demeaned - 1.96 * field_sd,
    z_upper = field_demeaned + 1.96 * field_sd)

  log_lik_val <- sum(w * fit$log_marginal[ok_cells])
  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(log_lik_val, n_draws),
    log_lik      = log_lik_val,
    N            = sum(model$valid)),
    .tobs_na_nuts_diagnostics(n_draws),
    list(
    col_names    = par_names,
    param_names  = par_names,
    process_info = pi_list,
    model        = model,
    spatial      = list(type = "icar", graph = adj,
                        sigma_mean = unname(hyper_means["sigma"]),
                        alpha_mean = unname(hyper_means["alpha"])),
    spatial_field = field_demeaned,
    field_table  = field_table,
    method       = "joint_coupled",
    positive     = model$positive,
    joint_fit    = fit,
    convergence  = list(converged = TRUE, n_iter = NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}
