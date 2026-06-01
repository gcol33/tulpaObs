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


# Demean each field block independently to the sum-to-zero convention each
# field's covariance already sits under. `vals` is the stacked per-cell field
# means (n_fields blocks of n_cells columns, in block order).
.occu_cover_demean_fields <- function(vals, n_cells, n_fields) {
  out <- numeric(length(vals))
  for (b in seq_len(n_fields)) {
    idx <- (b - 1L) * n_cells + seq_len(n_cells)
    out[idx] <- vals[idx] - mean(vals[idx])
  }
  out
}


# Resolve the optional trend spec. `trend` is `control$trend`:
#   * NULL / FALSE -> no trend field (single shared-intercept field).
#   * list(weight = "<col>") -> a per-cell numeric column of the cell data
#     weighting the trend field.
# Returns NULL when no trend, or a list with `time_cell` (length n_cells) and
# `weight` (the column name) otherwise.
.occu_cover_resolve_trend <- function(trend, model) {
  if (is.null(trend) || isFALSE(trend)) return(NULL)
  if (!is.list(trend) || is.null(trend$weight)) {
    stop("occu_cover() trend spec must be a list naming the weighting ",
         "covariate, e.g. control = list(trend = list(weight = \"time\")).",
         call. = FALSE)
  }
  wcol <- trend$weight
  if (!is.character(wcol) || length(wcol) != 1L) {
    stop("occu_cover() trend$weight must be a single column name.",
         call. = FALSE)
  }
  data <- model$data
  if (is.null(data) || !wcol %in% names(data)) {
    stop(sprintf(paste0(
      "occu_cover() trend$weight = '%s' is not a column of the cell data. ",
      "Supply it as a per-cell covariate in the data frame passed to tobs()."),
      wcol), call. = FALSE)
  }
  time_cell <- as.numeric(data[[wcol]])
  if (length(time_cell) != model$n_sites || any(!is.finite(time_cell))) {
    stop(sprintf(paste0(
      "occu_cover() trend$weight = '%s' must be a finite per-cell numeric ",
      "vector of length %d."), wcol, model$n_sites), call. = FALSE)
  }
  list(time_cell = time_cell, weight = wcol)
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
                                                  positive = "lognormal",
                                                  multi = FALSE,
                                                  n_cells = NULL,
                                                  site_cell = NULL) {
  n_sites    <- model$n_sites
  max_visits <- model$max_visits
  if (is.null(site_cell)) site_cell <- seq_len(n_sites)
  if (is.null(n_cells))   n_cells   <- max(site_cell)

  valid_flat  <- as.logical(t(model$valid))   # site-major: site 1 visits 1..J, ...
  y_flat      <- as.numeric(t(model$y))
  y_pos_flat  <- as.numeric(t(model$y_pos))
  site_flat   <- rep(seq_len(n_sites), each = max_visits)

  keep             <- which(valid_flat)
  n_visits_valid   <- length(keep)
  y_det_visit      <- as.integer(y_flat[keep])
  y_pos_visit      <- y_pos_flat[keep]
  # Each valid visit carries its occupancy unit (site, for the mixture grouping
  # via cell_obs_map) and its field node (cell, for the shared field via
  # spatial_idx). The two coincide when site_cell is the identity.
  site_of_visit    <- as.integer(site_flat[keep])
  cell_of_visit    <- as.integer(site_cell[site_of_visit])

  if (n_visits_valid == 0L) {
    stop("occu_cover joint_coupled: no valid visits in the data.",
         call. = FALSE)
  }

  # Site-level + visit-level fixed-effect design on each observation arm.
  X_p_site <- model$X_det_site
  X_p <- X_p_site[site_of_visit, , drop = FALSE]
  if (!is.null(model$X_det_visit)) {
    X_p <- cbind(X_p, model$X_det_visit[keep, , drop = FALSE])
  }
  X_pos_site <- model$X_pos_site
  X_pos <- X_pos_site[site_of_visit, , drop = FALSE]
  if (!is.null(model$X_pos_visit)) {
    X_pos <- cbind(X_pos, model$X_pos_visit[keep, , drop = FALSE])
  }

  # psi arm: one row per site (occupancy unit). spatial_idx maps the site to its
  # field node (cell); cell_obs_map indexes the occupancy unit itself. y /
  # n_trials / family are placeholders -- the per-obs scatter is skipped for
  # coupled = TRUE and the cell-coupling spec writes every derivative from the
  # per-site occupancy mixture.
  arm_psi <- list(
    y            = rep(0, n_sites),
    n_trials     = rep(0L, n_sites),
    X            = model$X_occ,
    spatial_idx  = as.integer(site_cell),
    family       = "binomial",
    phi          = 1.0,
    coupled      = TRUE,
    cell_obs_map = seq_len(n_sites)
  )

  # p arm: one row per valid visit. field_coef = 0 excludes the detection
  # predictor from every shared field (the per-arm field_coef multiplies the
  # field amplitude on EVERY block, so one scalar decouples the p arm from
  # both the intercept field and the trend field). spatial_idx is then a
  # placeholder.
  arm_p <- list(
    y            = as.numeric(y_det_visit),
    n_trials     = rep(1L, n_visits_valid),
    X            = X_p,
    spatial_idx  = rep(0L, n_visits_valid),
    family       = "binomial",
    phi          = 1.0,
    field_coef   = 0,
    coupled      = TRUE,
    cell_obs_map = site_of_visit
  )

  # pos arm: one row per valid visit. spatial_idx maps each visit to its cell.
  # `phi` is the pos-arm dispersion -- the lognormal SD on the log scale for
  # `positive = "lognormal"` and the beta precision for `positive = "beta"`
  # (the spec reads y_cell.phi(2) and interprets it per its policy). `family`
  # is unused for coupled arms (per-obs scatter + per-obs log-lik are both
  # skipped); we tag it with the positive family so the responses list reads
  # as intended.
  #
  # Coupling onto the field differs by path:
  #   * single-block (multi = FALSE): field_coef carries the alpha axis on the
  #     one shared block.
  #   * multi-block (multi = TRUE): the per-block copy spec carries each alpha
  #     axis, so the pos arm carries NO field_coef (the engine rejects copy +
  #     field_coef together).
  arm_pos <- list(
    y            = y_pos_visit,
    n_trials     = rep(1L, n_visits_valid),
    X            = X_pos,
    spatial_idx  = cell_of_visit,
    family       = positive,
    phi          = sigma_pos_init,
    coupled      = TRUE,
    cell_obs_map = site_of_visit
  )
  if (!multi) {
    arm_pos$field_coef <- list(name = "alpha", grid = alpha_grid)
  }

  list(responses      = list(psi = arm_psi, p = arm_p, pos = arm_pos),
       site_of_visit  = site_of_visit,
       cell_of_visit  = cell_of_visit,
       n_visits_valid = n_visits_valid)
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
.tobs_fit_occu_cover_joint_coupled <- function(model, fields,
                                                priors    = NULL,
                                                max.iter  = 200L,
                                                tol       = 1e-6,
                                                verbose   = TRUE,
                                                sigma.beta = 5,
                                                ...) {
  # `fields` is the coupled-field list from .occu_cover_spatial_fields(): the
  # unweighted intercept field first, then any weighted SVC fields. They share
  # one areal graph, so the base graph drives the (single) CSR.
  adj <- fields[[1L]]$graph
  is_beta <- identical(model$positive, "beta")
  is_lnrm <- identical(model$positive, "lognormal")
  if (!is_beta && !is_lnrm) {
    stop("occu_cover() joint_coupled engine supports positive = ",
         "\"lognormal\" or \"beta\".", call. = FALSE)
  }
  spec_name <- if (is_beta) "occu_cover_beta" else "occu_cover_lognormal"

  pi_list <- model$process_info
  # Field nodes (cells) and occupancy units (sites) are distinct under
  # group_var: many sites can share one cell field node. site_cell maps each
  # site to its node; absent, the two coincide 1:1.
  n_cells   <- nrow(adj)
  n_sites   <- model$n_sites
  site_cell <- model$site_cell %||% seq_len(n_sites)
  if (length(site_cell) != n_sites || max(site_cell) > n_cells ||
      min(site_cell) < 1L) {
    stop(sprintf(paste0(
      "occu_cover joint_coupled: site_cell must map %d sites into 1..%d ",
      "graph nodes."), n_sites, n_cells), call. = FALSE)
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

  # Coupled trend (SVC) fields: each is a per-cell-weighted areal field that
  # contributes weight_i * sigma_trend * z[cell_i] on occupancy and
  # weight_i * alpha_trend * sigma_trend * z[cell_i] on cover. They arrive
  # either as weighted areal terms in the formula (`fields[-1]`, each carrying a
  # resolved per-cell `$weight`) or via the back-compat `control$trend =
  # list(weight = "<col>")`. The two routes are mutually exclusive. With at
  # least one trend field the fit takes the multi-block copy path; absent, the
  # single shared-intercept field.
  coupled_trends <- lapply(fields[-1L], function(f) {
    list(weight = f$weight, weight_label = f$weight_label %||% "trend")
  })
  trend_spec <- .occu_cover_resolve_trend(dots$trend, model)
  if (!is.null(trend_spec)) {
    if (length(coupled_trends) > 0L) {
      stop("occu_cover(): give the trend field EITHER as a weighted areal term ",
           "in the formula (icar(graph = adj, weight = col)) OR via ",
           "control$trend, not both.", call. = FALSE)
    }
    coupled_trends <- list(list(weight = trend_spec$time_cell,
                                weight_label = trend_spec$weight))
  }
  n_trend   <- length(coupled_trends)
  has_trend <- n_trend > 0L

  arms_out <- .occu_cover_build_joint_coupled_arms(
    model           = model,
    sigma_pos_init  = sigma_pos_init,
    alpha_grid      = alpha_grid,
    positive        = model$positive,
    multi           = has_trend,
    n_cells         = n_cells,
    site_cell       = site_cell
  )
  responses      <- arms_out$responses
  site_of_visit  <- arms_out$site_of_visit
  cell_of_visit  <- arms_out$cell_of_visit
  n_v            <- arms_out$n_visits_valid

  csr <- .occu_cover_adj_to_csr(adj)
  icar_template <- function(extra = list()) {
    c(list(
        type            = "icar",
        n_spatial_units = csr$n_spatial_units,
        adj_row_ptr     = csr$adj_row_ptr,
        adj_col_idx     = csr$adj_col_idx,
        n_neighbors     = csr$n_neighbors,
        sigma_grid      = sigma_grid
      ), extra)
  }

  # Optional sigma_pos integration via a phi_grid axis on the pos arm.
  phi_grid_pos <- dots$phi.grid.pos
  phi_grid_arg <- if (!is.null(phi_grid_pos))
                    list(pos = as.numeric(phi_grid_pos))
                  else NULL

  if (has_trend) {
    # Multi-block path: the intercept ICAR block plus one ICAR block per coupled
    # trend field, all on the same graph and each copied onto the pos arm with
    # its own alpha axis. The p arm is excluded from every field via its
    # field_coef = 0. Per-block svc_weight injects the per-row field weight on
    # the psi (per-cell) and pos (per-visit) arms; the p-arm weight is
    # irrelevant (field_coef = 0 already zeroes the p field).
    # Field node per arm row: psi rows are sites (-> site_cell), p / pos rows are
    # visits (-> cell_of_visit). The SVC weight is per occupancy unit (site) on
    # the psi arm and per visit's site on the pos arm; the field it multiplies
    # is the per-cell node addressed by spatial_idx.
    spatial_idx_arms <- list(as.integer(site_cell), cell_of_visit, cell_of_visit)
    make_block <- function(weight_site) {
      w_psi   <- if (is.null(weight_site)) rep(1.0, n_sites)
                 else as.numeric(weight_site)
      w_visit <- w_psi[site_of_visit]
      icar_template(list(
        spatial_idx = spatial_idx_arms,
        svc_weight  = list(w_psi, rep(1.0, n_v), w_visit)
      ))
    }
    alpha_grid_trend <- dots$alpha.grid.trend %||% alpha_grid
    prior_arg <- c(
      list(make_block(NULL)),
      lapply(coupled_trends, function(tf) make_block(tf$weight))
    )
    copy_arg <- c(
      list(list(arm = "pos", block = 1L, alpha_grid = alpha_grid)),
      lapply(seq_len(n_trend), function(j)
        list(arm = "pos", block = j + 1L, alpha_grid = alpha_grid_trend))
    )
  } else {
    prior_arg <- icar_template()
    copy_arg  <- NULL
  }

  fit_call <- list(
    responses     = responses,
    prior         = prior_arg,
    phi_grid      = phi_grid_arg,
    cell_coupling = spec_name,
    control = list(
      max_iter  = as.integer(max.iter),
      tol       = as.numeric(tol),
      n_threads = as.integer(dots$n.threads %||% 1L),
      store_Q   = TRUE,
      # Inner-Newton curvature (gcol33/tulpa#46). The beta positive arm's
      # observed mixture Hessian is indefinite away from the mode, so observed-
      # curvature Newton steps stall and the inner Newton hits max.iter in every
      # grid cell (non-convergence -- the dominant cost). Expected/Fisher
      # information is PSD by construction, so the steps are well-conditioned and
      # the inner Newton converges in ~12 steps instead. The reported SEs,
      # log_det and grid weights are unchanged: the final mode-pass always
      # re-factorizes with the observed Hessian; the curvature mode only steers
      # the path to the mode. The lognormal arm is exactly quadratic (one inner
      # step), so observed curvature is already optimal -> keep "lm".
      hessian   = dots$hessian %||% (if (is_beta) "fisher" else "lm"),
      # Cholesky factor reuse (Shamanskii / chord) is exposed but defaults off:
      # at realistic field sizes the sparse factorization is milliseconds, so
      # reuse buys nothing once the iteration count is fixed by the curvature.
      inner_refresh = as.integer(dots$inner.refresh %||% 1L),
      # Outer-grid parallelism (gcol33/tulpa#46, lever 2). The cover hurdle's
      # large spatial field takes the sparse inner-solve path, whose outer grid
      # now runs across `n.threads.outer` threads (per-thread Hessian builder /
      # scratch / specs). Defaults to serial; set it for the full-field fits.
      # `force.sparse` forces the sparse path on small fields (testing / the
      # parallel and factor-reuse paths live there).
      n_threads_outer = as.integer(dots$n.threads.outer %||% 1L),
      force_sparse    = isTRUE(dots$force.sparse),
      # Adaptive-grid refinement defaults ON. Non-convergent inner Newton
      # cells (degenerate sigma + small non-zero alpha hyperpoints) drop to
      # -Inf log_marginal under the engine's NaN-safe edge-score path
      # (tulpa/R/hyper_grid_refine.R::.hyper_axis_edge_scores), so the
      # refinement walks finite mass only and never trips on a missing-value
      # threshold compare.
      adaptive_grid             = dots$adaptive.grid             %||% TRUE,
      adaptive_grid_edge_thresh = dots$adaptive.grid.edge.thresh %||% 0.02,
      adaptive_grid_max_passes  = dots$adaptive.grid.max.passes  %||% 1L,
      # Outer Pareto-k accuracy diagnostic. The engine draws `k_samples`
      # hyperparameter points and re-solves the inner Laplace at each; on a large
      # field the Gaussian proposal lands many draws at extreme sigma where the
      # inner Newton stalls to max.iter, so the diagnostic can cost far more than
      # the grid integration itself (~50x at EVA scale). Forward the knobs so a
      # production fit can disable it (`control$diagnose.k = FALSE`) or shrink the
      # sample; small fits keep the engine default.
      diagnose_k = dots$diagnose.k %||% TRUE,
      k_samples  = as.integer(dots$k.samples %||% 200L),
      # Grid-cell checkpoint/resume (gcol33/tulpa#50). An EVA-scale occu_cover
      # fit runs for hours; `control$checkpoint = list(path =, resume =)` makes
      # the outer grid append each completed cell to `path` and a resume run
      # load the finished cells and solve only the rest, so a killed/rebooted
      # fit resumes instead of restarting. Forwarded verbatim to the engine.
      checkpoint = dots$checkpoint
    )
  )
  if (!is.null(copy_arg)) fit_call$copy <- copy_arg

  fit <- do.call(tulpa::tulpa_nested_laplace_joint, fit_call)

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
  # One latent field block per coupled spatial field. The multi-block layout
  # reports every ICAR/structured field's 0-based offset in `field_starts`;
  # the single-block joint layout reports the one field's offset in `phi_start`.
  field_starts0 <- layout$field_starts %||% layout$phi_start
  n_fields      <- length(field_starts0)
  field_idx     <- as.integer(unlist(lapply(field_starts0,
                                            function(s0) s0 + seq_len(n_cells))))
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

    field_modes   <- modes[, field_idx, drop = FALSE]
    field_at_cell <- as.numeric(crossprod(w, field_modes))
    field_var     <- as.numeric(crossprod(w, field_modes^2)) - field_at_cell^2
    field_demeaned <- .occu_cover_demean_fields(field_at_cell, n_cells, n_fields)
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
    # sum-to-zero convention the field-block covariance already sits under. One
    # block of n_cells columns per coupled field, in block order.
    n_field_cols   <- n_fields * n_cells
    field_at_cell  <- mbar_joint[p_beta + seq_len(n_field_cols)]
    field_var      <- diag_Vj[p_beta + seq_len(n_field_cols)]
    field_demeaned <- .occu_cover_demean_fields(field_at_cell, n_cells, n_fields)
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
  # pick2 reads a multi-block axis column (`b<k>.<axis>`) under a public name.
  pick2 <- function(public, col) {
    j <- match(col, tg_names)
    if (is.na(j)) return(invisible(NULL))
    vals <- as.numeric(tg_ok[, j])
    m <- sum(w * vals)
    v <- sum(w * vals^2) - m^2
    hyper_means[[public]] <<- m
    hyper_sds  [[public]] <<- sqrt(max(v, 0))
    hyper_names <<- c(hyper_names, public)
  }
  if (has_trend) {
    # Multi-block: block 1 is the intercept field, blocks 2.. the trend fields.
    # A single trend field keeps the bare sigma_trend/alpha_trend names; several
    # are indexed (sigma_trend1, alpha_trend1, ...).
    pick2("sigma", "b1.sigma")
    pick2("alpha", "b1.alpha")
    for (j in seq_len(n_trend)) {
      suffix <- if (n_trend == 1L) "" else as.character(j)
      pick2(paste0("sigma_trend", suffix), sprintf("b%d.sigma", j + 1L))
      pick2(paste0("alpha_trend", suffix), sprintf("b%d.alpha", j + 1L))
    }
    pick("phi_pos")
  } else {
    for (nm in c("sigma", "alpha", "phi_pos")) pick(nm)
  }
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

  # Split the stacked per-field summaries into one block of n_cells per coupled
  # field. Field 1 is the intercept field (the back-compat `spatial_field`);
  # fields 2.. are the spatially-varying trend fields, in block order.
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
  # Back-compat single-trend accessors (the first trend field).
  field_trend       <- if (length(trend_means))  trend_means[[1L]]  else NULL
  trend_field_table <- if (length(trend_tables)) trend_tables[[1L]] else NULL

  # Joint betas+field posterior for downstream derived-quantity prediction
  # (delta_p / delta_cover marginalized over the full correlated posterior).
  # `joint_means` carries every field in the same demeaned convention as
  # `spatial_field`; `joint_vcov` is the law-of-total-covariance Vj (NULL on
  # the older-tulpa diagonal fallback). Fields are stacked in block order
  # (intercept then trend fields).
  field_par_names <- unlist(c(
    list(paste0("field_", seq_len(n_cells))),
    lapply(seq_len(n_fields - 1L), function(j) {
      suffix <- if (n_fields - 1L == 1L) "" else as.character(j)
      paste0("trend_field", suffix, "_", seq_len(n_cells))
    })
  ))
  joint_par_names <- c(
    paste0("psi_",   pi_list[[1L]]$coef_names),
    paste0("p_",     pi_list[[2L]]$coef_names),
    paste0("pos_",   pi_list[[3L]]$coef_names),
    field_par_names
  )
  joint_means <- c(beta_psi_m, beta_p_m, beta_pos_m, field_demeaned)
  names(joint_means) <- joint_par_names
  if (!is.null(Vj)) dimnames(Vj) <- list(joint_par_names, joint_par_names)

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
                        alpha_mean = unname(hyper_means["alpha"]),
                        sigma_trend_mean = if (has_trend)
                          unname(hyper_means[if (n_trend == 1L) "sigma_trend"
                                             else "sigma_trend1"]) else NULL,
                        alpha_trend_mean = if (has_trend)
                          unname(hyper_means[if (n_trend == 1L) "alpha_trend"
                                             else "alpha_trend1"]) else NULL),
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
    method       = "joint_coupled",
    positive     = model$positive,
    joint_fit    = fit,
    convergence  = list(converged = TRUE, n_iter = NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}
