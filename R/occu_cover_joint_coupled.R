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
                                                  alpha_grid,
                                                  positive = "lognormal") {
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
  # contributes the cover-arm field. `phi` is the pos-arm dispersion --
  # the lognormal SD on the log scale for `positive = "lognormal"` and the
  # beta precision for `positive = "beta"` (the spec reads y_cell.phi(2)
  # and interprets it per its policy). `family` is unused for coupled arms
  # (per-obs scatter + per-obs log-lik are both skipped); we tag it with
  # the positive family so the responses list reads as intended.
  arm_pos <- list(
    y            = y_pos_visit,
    n_trials     = rep(1L, n_visits_valid),
    X            = X_pos,
    spatial_idx  = cell_idx_visit,
    family       = positive,
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
  is_beta <- identical(model$positive, "beta")
  is_lnrm <- identical(model$positive, "lognormal")
  if (!is_beta && !is_lnrm) {
    stop("occu_cover() joint_coupled engine supports positive = ",
         "\"lognormal\" or \"beta\".", call. = FALSE)
  }
  spec_name <- if (is_beta) "occu_cover_beta" else "occu_cover_lognormal"

  pi_list <- model$process_info
  n_cells <- model$n_sites
  if (nrow(adj) != n_cells) {
    stop(sprintf("Spatial graph has %d nodes but data has %d cells.",
                 nrow(adj), n_cells), call. = FALSE)
  }

  dots <- list(...)

  # Pre-fit the pos-arm dispersion at the empirical sample value at detected
  # visits. For lognormal, that's the SD of log(y_pos) (matching the v3
  # nested-Laplace warm start for log_disp); for beta, a moment-matched
  # precision from the sample mean and variance of y_pos in (0, 1).
  pos_vals <- model$y_pos[model$valid & model$y == 1L]
  phi_pos_init <- if (is_beta) {
    if (length(pos_vals) >= 2L) {
      mu_hat   <- mean(pos_vals)
      var_hat  <- max(stats::var(pos_vals), 1e-6)
      max((mu_hat * (1 - mu_hat)) / var_hat - 1, 1)
    } else {
      10
    }
  } else {
    if (length(pos_vals) > 0L) {
      max(stats::sd(log(pos_vals)), 0.05) + 0.05
    } else {
      0.4
    }
  }
  sigma_pos_init <- phi_pos_init  # passed through as pos-arm phi

  alpha_grid <- dots$alpha.grid %||%
                c(0, exp(seq(log(0.1), log(3), length.out = 5)))
  sigma_grid <- dots$sigma.grid %||%
                exp(seq(log(0.1), log(3), length.out = 5))

  responses <- .occu_cover_build_joint_coupled_arms(
    model           = model,
    sigma_pos_init  = sigma_pos_init,
    alpha_grid      = alpha_grid,
    positive        = model$positive
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
    cell_coupling = spec_name,
    control = list(
      max_iter  = as.integer(max.iter),
      tol       = as.numeric(tol),
      n_threads = as.integer(dots$n.threads %||% 1L),
      store_Q   = TRUE,
      hessian   = dots$hessian %||% "lm",
      # Adaptive-grid refinement defaults ON. Non-convergent inner Newton
      # cells (degenerate sigma + small non-zero alpha hyperpoints) drop to
      # -Inf log_marginal under the engine's NaN-safe edge-score path
      # (tulpa/R/hyper_grid_refine.R::.hyper_axis_edge_scores), so the
      # refinement walks finite mass only and never trips on a missing-value
      # threshold compare.
      adaptive_grid             = dots$adaptive.grid             %||% TRUE,
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

  # Joint (betas + field) posterior covariance by the law of total covariance
  # over the outer hyperparameter grid:
  #
  #   Cov(x | y) = sum_k w_k [ Cov(x | y, theta_k) + (m_k - mbar)(m_k - mbar)' ]
  #
  # where x = (beta_psi, beta_p, beta_pos, field) stacked, Cov(x | y, theta_k)
  # is the inner-Laplace covariance at grid cell k (the ICAR sum-to-zero
  # constrained sub-block of Q_k^-1, so the intercept covariance is
  # data-identified rather than collapsing along the (intercept, mean(phi))
  # near-null direction of the improper field prior), and (m_k - mbar) is the
  # between-grid mode deviation. Carrying the field block (not just the betas)
  # means downstream derived quantities (delta_p, delta_cover) can marginalize
  # the joint betas+field posterior instead of a marginal-only diagonal. The
  # constrained per-grid block comes from `.joint_inner_vcov_block()`
  # (family_cover_hurdle.R) -- the dense-block analogue of `.joint_inner_var()`,
  # same constraint correction, one source of truth across families.
  p_beta    <- p_psi + p_p + p_pos
  field_idx <- layout$phi_start + seq_len(n_cells)
  idx_joint <- c(bpsi_idx, bp_idx, bpos_idx, field_idx)
  blocks    <- .joint_inner_vcov_block(fit, idx_joint)

  if (is.null(blocks)) {
    # Older tulpa without stored per-grid Q: fall back to the marginal-only
    # diagonal (var-of-means plus the diagonal inner-Laplace variance) so the
    # fit still completes, with no betas+field cross-covariance.
    beta_idx_all <- c(bpsi_idx, bp_idx, bpos_idx)
    inner_var    <- .joint_inner_var(fit, beta_idx_all)
    total_var <- function(modes_block, mean_vec, iv_block) {
      vom <- as.numeric(crossprod(w, modes_block^2)) - mean_vec^2
      mov <- if (is.null(iv_block)) {
        0
      } else {
        iv_k <- iv_block[ok_cells, , drop = FALSE]
        iv_k[!is.finite(iv_k)] <- 0  # rank-deficient Q -> var-of-means only
        as.numeric(crossprod(w, iv_k))
      }
      pmax(vom + mov, 0)
    }
    iv_psi <- if (is.null(inner_var)) NULL
              else inner_var[, seq_len(p_psi), drop = FALSE]
    iv_p   <- if (is.null(inner_var)) NULL
              else inner_var[, p_psi + seq_len(p_p), drop = FALSE]
    iv_pos <- if (is.null(inner_var)) NULL
              else inner_var[, p_psi + p_p + seq_len(p_pos), drop = FALSE]
    sds_beta <- c(
      sqrt(total_var(modes[, bpsi_idx, drop = FALSE], beta_psi_m, iv_psi)),
      sqrt(total_var(modes[, bp_idx,   drop = FALSE], beta_p_m,   iv_p)),
      sqrt(total_var(modes[, bpos_idx, drop = FALSE], beta_pos_m, iv_pos))
    )
    beta_block <- diag(sds_beta^2, nrow = p_beta)

    field_modes <- modes[, field_idx, drop = FALSE]
    field_at_cell <- as.numeric(crossprod(w, field_modes))
    field_var <- as.numeric(crossprod(w, field_modes^2)) - field_at_cell^2
    field_demeaned <- field_at_cell - mean(field_at_cell)
    Vj <- NULL  # no joint covariance available
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
    Vj <- (Vj + t(Vj)) / 2  # symmetrize off floating-point constraint residuals
    diag_Vj    <- diag(Vj)
    sds_beta   <- sqrt(pmax(diag_Vj[seq_len(p_beta)], 0))
    beta_block <- Vj[seq_len(p_beta), seq_len(p_beta), drop = FALSE]

    # Field summary uses the full (within + between) variance, demeaned to the
    # sum-to-zero convention the field-block covariance already sits under.
    field_at_cell  <- mbar_joint[p_beta + seq_len(n_cells)]
    field_var      <- diag_Vj[p_beta + seq_len(n_cells)]
    field_demeaned <- field_at_cell - mean(field_at_cell)
  }
  field_sd <- sqrt(pmax(field_var, 0))

  sd_psi <- sds_beta[seq_len(p_psi)]
  sd_p   <- sds_beta[p_psi + seq_len(p_p)]
  sd_pos <- sds_beta[p_psi + p_p + seq_len(p_pos)]

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

  # Parameter-surface vcov: betas correlated (the beta block of the joint
  # covariance), hyperparameters diagonal (grid-summarised, no cross-covariance
  # with the betas). Block-diagonal across the two.
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
    cell    = seq_len(n_cells),
    z_mean  = field_demeaned,
    z_sd    = field_sd,
    z_lower = field_demeaned - 1.96 * field_sd,
    z_upper = field_demeaned + 1.96 * field_sd
  )

  # Joint betas+field posterior for downstream derived-quantity prediction
  # (delta_p / delta_cover marginalized over the full correlated posterior).
  # `joint_means` carries the field in the same demeaned convention as
  # `spatial_field`; `joint_vcov` is the law-of-total-covariance Vj (NULL on
  # the older-tulpa diagonal fallback).
  joint_par_names <- c(
    paste0("psi_",   pi_list[[1L]]$coef_names),
    paste0("p_",     pi_list[[2L]]$coef_names),
    paste0("pos_",   pi_list[[3L]]$coef_names),
    paste0("field_", seq_len(n_cells))
  )
  joint_means <- c(beta_psi_m, beta_p_m, beta_pos_m, field_demeaned)
  names(joint_means) <- joint_par_names
  if (!is.null(Vj)) dimnames(Vj) <- list(joint_par_names, joint_par_names)

  structure(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
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
    spatial_field = field_demeaned,
    field_table  = field_table,
    joint_par_names = joint_par_names,
    joint_means     = joint_means,
    joint_vcov      = Vj,
    method       = "joint_coupled",
    positive     = model$positive,
    joint_fit    = fit,
    convergence  = list(converged = TRUE, n_iter = NA_integer_)
  ), class = c("tobs_fit", "tulpa_fit"))
}
