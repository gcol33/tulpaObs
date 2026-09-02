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
.tobs_fit_occu_cover_joint <- function(model, fields,
                                                priors    = NULL,
                                                re_spec   = NULL,
                                                correlated = FALSE,
                                                pos_armspec = NULL,
                                                det_armspec = NULL,
                                                max.iter  = 200L,
                                                tol       = 1e-6,
                                                verbose   = TRUE,
                                                sigma.beta = 5,
                                                .batch_collect = FALSE,
                                                ...) {
  # `fields` is the coupled-field list from .occu_cover_spatial_fields(): the
  # unweighted intercept field first, then any weighted SVC fields. They share
  # one areal graph, so the base graph drives the (single) CSR.
  adj <- fields[[1L]]$graph
  is_beta  <- identical(model$positive, "beta")
  is_lnrm  <- identical(model$positive, "lognormal")
  is_gauss <- identical(model$positive, "gaussian")
  if (!is_beta && !is_lnrm && !is_gauss) {
    stop("occu_cover() joint engine supports positive = ",
         "\"lognormal\", \"beta\", or \"gaussian\".", call. = FALSE)
  }
  # Cover-arm granularity. "mean" / "median" route through the `_agg` spec (one
  # log f_pos at the per-unit mean / median); "latent" routes through the
  # stateful `_latent` spec (a per-unit cover RE integrated out, one marginal
  # per unit); "none" is the per-visit spec.
  cover_aggregate <- model$cover_aggregate %||% "none"
  is_latent  <- identical(cover_aggregate, "latent")
  if (is_latent && is_gauss) {
    stop("occu_cover(response = \"gaussian\") has no latent cover-aggregate ",
         "variant; use cover_aggregate = \"none\" / \"mean\" / \"median\".",
         call. = FALSE)
  }
  aggregated <- !identical(cover_aggregate, "none") && !is_latent
  spec_base  <- if (is_beta) "occu_cover_beta"
                else if (is_gauss) "occu_cover_gaussian"
                else "occu_cover_lognormal"
  spec_name  <- if (is_latent) paste0(spec_base, "_latent")
                else if (aggregated) paste0(spec_base, "_agg")
                else spec_base

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
      "occu_cover joint: site_cell must map %d sites into 1..%d ",
      "graph nodes."), n_sites, n_cells), call. = FALSE)
  }

  dots <- list(...)

  # Per-group random intercepts. The occupancy-arm RE (`re_spec`, one code per
  # site) and the observation-arm REs (model$re_det / model$re_pos, one code per
  # detection / positive-cover row) each join the fit as an `iid` prior block
  # whose per-group latent rides ONE arm. Any RE block forces the multi-block
  # driver (the field amplitude becomes an explicit copy spec). Not composed with
  # the cover-latent RE (the latent spec carries its own per-unit cover RE) nor
  # with the batched fused path (one species at a time).
  has_re     <- !is.null(re_spec)            # occupancy (psi) arm
  has_re_det <- !is.null(model$re_det)       # detection (p) arm
  has_re_pos <- !is.null(model$re_pos)       # positive-cover arm
  has_any_re <- has_re || has_re_det || has_re_pos
  # Arm-specific cover field: an independent, non-copied ICAR block on the cover
  # (pos) arm alone, composed with the shared occupancy field. Forces the
  # multi-block driver (like a trend field / RE block). Not composed with the
  # latent cover RE, the correlated MCAR (gated at parse), or the batched fused
  # path (one species at a time, no extra block). Arm-specific fields carry the
  # detection (p) arm as well as the cover (pos) arm (each an independent,
  # non-copied ICAR block on that arm alone). Both force the multi-block driver.
  has_pos_armspec <- !is.null(pos_armspec)
  has_det_armspec <- !is.null(det_armspec)
  has_armspec     <- has_pos_armspec || has_det_armspec
  if (has_any_re && is_latent) {
    stop("occu_cover(): a per-group RE and cover_aggregate = \"latent\" cannot ",
         "be combined (the latent path carries its own per-unit cover RE).",
         call. = FALSE)
  }
  if (has_any_re && isTRUE(.batch_collect)) {
    stop("occu_cover(): the batched fused path does not carry a per-group RE ",
         "block.", call. = FALSE)
  }
  if (has_armspec && is_latent) {
    stop("occu_cover(): an arm-specific field (to = \"positive\" / \"detection\") ",
         "does not compose with cover_aggregate = \"latent\".", call. = FALSE)
  }
  if (has_armspec && isTRUE(.batch_collect)) {
    stop("occu_cover(): the batched fused path does not carry an arm-specific ",
         "field.", call. = FALSE)
  }

  # Correlated (`|`) free-Sigma MCAR field: one coupled block over the bar's
  # intercept + coefficient fields, copied onto the cover arm with one amplitude
  # alpha. Scoped to the standard (non-latent, unbatched) path; the latent cover
  # RE and the fused batch driver are not composed with it.
  if (correlated) {
    if (is_latent) {
      stop("occu_cover(): a correlated spatial bar (`|`, free-Sigma MCAR) does ",
           "not compose with cover_aggregate = \"latent\".", call. = FALSE)
    }
    if (isTRUE(.batch_collect)) {
      stop("occu_cover(): the batched fused path does not carry a correlated ",
           "MCAR field.", call. = FALSE)
    }
    if (has_any_re) {
      stop("occu_cover(): a per-group RE does not compose with a correlated ",
           "spatial bar (`|`, free-Sigma MCAR) on the joint engine.",
           call. = FALSE)
    }
  }

  # Pre-fit the pos-arm dispersion(s). The non-latent paths carry a single
  # dispersion on the pos arm's phi slot; the latent path carries the integrated
  # cover-latent SD (sigma_u) there and holds a SECOND, within-unit dispersion
  # (disp2_fixed) fixed in the stateful spec.
  disp2_fixed <- NULL
  if (is_latent) {
    # The within-unit dispersion (disp2) is FIXED and captured in the spec;
    # sigma_u (the integrated cover-latent SD) rides the pos arm's phi axis.
    # Pre-fit disp2 from the WITHIN-unit spread and seed sigma_u from the
    # BETWEEN-unit spread: Var(log y) = sigma_eps^2 + sigma_u^2, so pre-fitting
    # disp2 at the total spread would swallow sigma_u and leave it unidentified.
    det_mat   <- model$valid & (model$y == 1L) & is.finite(model$y_pos)
    det_sites <- which(rowSums(det_mat) > 0L)
    site_vals <- lapply(det_sites, function(i)
                        as.numeric(model$y_pos[i, det_mat[i, ]]))
    has_2 <- length(det_sites) >= 2L
    if (is_beta) {
      all_v   <- unlist(site_vals)
      mu_hat  <- mean(all_v)
      var_hat <- max(stats::var(all_v), 1e-6)
      disp2_fixed  <- max((mu_hat * (1 - mu_hat)) / var_hat - 1, 1)
      site_mu      <- vapply(site_vals, function(v)
        stats::qlogis(min(max(mean(v), 1e-3), 1 - 1e-3)), numeric(1))
      sigma_u_init <- if (has_2) max(stats::sd(site_mu), 0.1) else 0.5
    } else {
      logvals    <- lapply(site_vals, log)
      site_means <- vapply(logvals, mean, numeric(1))
      m_per      <- lengths(logvals)
      within_ss  <- sum(vapply(logvals, function(lv)
        if (length(lv) >= 2L) sum((lv - mean(lv))^2) else 0, numeric(1)))
      within_df  <- sum(pmax(m_per - 1L, 0L))
      disp2_fixed  <- if (within_df > 0L) max(sqrt(within_ss / within_df), 0.05)
                      else max(stats::sd(unlist(logvals)), 0.05)
      between_var  <- if (has_2)
                        stats::var(site_means) - disp2_fixed^2 / mean(m_per)
                      else NA_real_
      sigma_u_init <- if (is.finite(between_var))
                        max(sqrt(max(between_var, 1e-4)), 0.1) else 0.5
    }
    sigma_pos_init <- sigma_u_init
  } else {
    # Pre-fit the single pos-arm dispersion at the empirical sample value of the
    # cover observations the arm actually models. For lognormal, the SD of
    # log(y_pos); for beta, a moment-matched precision. Under mean / median
    # aggregation the modelled observation is the per-unit mean / median, so the
    # dispersion is pre-fit on those aggregated values.
    pos_vals <- if (aggregated) {
      aggfun  <- if (identical(cover_aggregate, "median")) stats::median else mean
      det_mat <- model$valid & (model$y == 1L) & is.finite(model$y_pos)
      sw      <- which(rowSums(det_mat) > 0L)
      vapply(sw, function(i) as.numeric(aggfun(model$y_pos[i, det_mat[i, ]])),
             numeric(1))
    } else {
      pv <- model$y_pos[model$valid & model$y == 1L]
      pv[is.finite(pv)]
    }
    phi_pos_init <- if (is_beta) {
      if (length(pos_vals) >= 2L) {
        mu_hat   <- mean(pos_vals)
        var_hat  <- max(stats::var(pos_vals), 1e-6)
        max((mu_hat * (1 - mu_hat)) / var_hat - 1, 1)
      } else {
        10
      }
    } else if (is_gauss) {
      # Identity-Gaussian arm: sigma_pos is the SD of the raw response (which
      # may be negative), not of log(y_pos).
      if (length(pos_vals) > 0L) max(stats::sd(pos_vals), 0.05) + 0.05 else 0.4
    } else {
      if (length(pos_vals) > 0L) {
        max(stats::sd(log(pos_vals)), 0.05) + 0.05
      } else {
        0.4
      }
    }
    sigma_pos_init <- phi_pos_init  # passed through as pos-arm phi
  }

  # A defaulted axis arrives already marked from the helper and reaches its block
  # unsorted, so the mark survives; the pos-arm amplitude below re-derives
  # through `sort()` and is marked again there.
  alpha_axis <- .tobs_alpha_axis_base(dots)
  sigma_grid <- dots$sigma.grid %||% .tobs_default_sigma_grid()

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

  arms_out <- .occu_cover_build_joint_arms(
    model           = model,
    sigma_pos_init  = sigma_pos_init,
    alpha_axis      = alpha_axis,
    positive        = model$positive,
    multi           = has_trend || has_any_re || has_armspec,
    n_cells         = n_cells,
    site_cell       = site_cell,
    cover_aggregate = cover_aggregate,
    det_field       = has_det_armspec,
    # Detection-pattern compression (exact): on for the single-species path,
    # off for the batched fused solve (per-species detection differs, so the
    # nodet rows are not exchangeable across species). getOption escape hatch
    # (default on) lets a fit force the uncompressed build for an equivalence check.
    compress_nodet  = !isTRUE(.batch_collect) &&
      isTRUE(getOption("tulpaObs.compress_nodet", TRUE))
  )
  responses      <- arms_out$responses
  site_of_visit  <- arms_out$site_of_visit
  cell_of_visit  <- arms_out$cell_of_visit
  n_v            <- arms_out$n_visits_valid
  pos_site       <- arms_out$pos_site
  n_pos_rows     <- arms_out$n_pos_rows
  pos_cover_values <- arms_out$pos_cover_values
  re_det_terms   <- arms_out$re_det_terms
  re_pos_terms   <- arms_out$re_pos_terms

  # Attach the per-arm fixed-effect priors. These reach tulpa's joint engine as
  # per-arm `beta_prior_mean` / `beta_prior_prec` on each response and replace
  # the engine's uniform weak default in add_per_arm_beta_re_priors().
  arm_priors <- .occu_cover_coupled_arm_priors(priors, responses)
  for (nm in c("psi", "p", "pos")) {
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
        sigma_grid      = sigma_grid
      ), extra)
  }

  # Per-group RE blocks. Each random-intercept term (one per arm for #56/#102;
  # several per arm for crossed / nested #103) contributes one `iid` prior block
  # whose per-group latent rides THAT arm only: obs_idx is the 3-element (psi, p,
  # pos) list of per-row group codes, with the targeted arm carrying the codes
  # and the other two zeroed (0 = no RE for that row, the engine's scatter skip).
  # Each block's SD integrates on the outer grid over its `sigma_grid`. They
  # trail the field block(s), so the field copy indices stay valid; the blocks
  # carry no copy (each rides its own arm). `re_descs` records each block's arm +
  # grouping metadata in prior order so the postprocess maps each block back to
  # its `b<k>.sigma` axis and BLUP columns. Each per-group RE block's SD axis.
  # Declared as ours whenever the matching `control$re.sigma.grid*` is unset, so
  # the engine may recentre it; `re_auto()` pairs each call site with its own
  # knob.
  re_grid_default <- exp(seq(log(0.05), log(2), length.out = 6L))
  re_auto <- function(user) .tobs_mark_auto(user %||% re_grid_default,
                                            is.null(user))
  re_blocks <- list()
  re_descs  <- list()   # one descriptor per RE TERM (may span several blocks)
  zero_psi <- rep(0L, n_sites)
  zero_p   <- rep(0L, n_v)
  zero_pos <- rep(0L, n_pos_rows)
  # obs_idx: per-row group codes with only the targeted arm carrying codes (the
  # other two zeroed: 0 = no RE for that row, the engine's scatter skip).
  obs_for <- function(arm, codes) switch(arm,
    psi = list(as.integer(codes), zero_p, zero_pos),
    p   = list(zero_psi, as.integer(codes), zero_pos),
    pos = list(zero_psi, zero_p, as.integer(codes)))
  # Per-arm design-weight list: the targeted arm gets `w` (a coefficient's design
  # column); the other two get unit weights of the right length (unused under
  # their 0 obs_idx, but length-matched for the engine's per-arm validation).
  wt_for <- function(arm, w) {
    base <- list(rep(1.0, n_sites), rep(1.0, n_v), rep(1.0, n_pos_rows))
    base[[match(arm, c("psi", "p", "pos"))]] <- as.numeric(w)
    base
  }
  # Emit the prior block(s) for one RE term and record ONE descriptor spanning
  # them:
  #   * intercept (n_coefs == 1, !correlated): one scalar `iid` block.
  #   * uncorrelated slope (!correlated, Z present): one weighted `iid` block per
  #     coefficient -- svc_weight = that coefficient's design column (the
  # intercept column is all-ones, so its block is the scalar iid).
  #   * correlated slope: one multivariate `miid` block over a free cross-coef
  # Sigma -- field_weight = the design columns.
  # The descriptor records the [block_start, n_blocks] run, in emission order, so
  # the postprocess pulls the right latent columns; the blocks trail the field
  # block(s), so the field copy indices stay valid.
  add_re_term <- function(arm, codes, grid, Z, n_groups, var, levels,
                          n_coefs, coef_names, correlated, logchol_grid = NULL,
                          coef_scales = NULL) {
    obs_idx <- obs_for(arm, codes)
    b0 <- length(re_blocks)
    # `as.numeric()` drops the auto-grid marker the caller applied, so re-apply
    # it to the vector that actually reaches the block.
    grid <- .tobs_num_auto(grid)
    if (!isTRUE(correlated)) {
      for (cc in seq_len(n_coefs)) {
        blk <- list(type = "iid", n_units = as.integer(n_groups),
                    sigma_grid = grid, obs_idx = obs_idx)
        if (!is.null(Z)) blk$svc_weight <- wt_for(arm, Z[, cc])
        re_blocks[[length(re_blocks) + 1L]] <<- blk
      }
    } else {
      field_weight <- lapply(seq_len(n_coefs),
                             function(cc) wt_for(arm, Z[, cc]))
      blk <- list(type = "miid", n_groups = as.integer(n_groups),
                  n_fields = as.integer(n_coefs), obs_idx = obs_idx,
                  field_weight = field_weight)
      # A coarse free-Sigma grid (or the user's) so the block composes with the
      # shared field + copy under the engine's outer-grid cap.
      lc <- logchol_grid %||% .occu_cover_miid_logchol_grid(n_coefs)
      if (!is.null(lc)) blk$logchol_grid <- as.matrix(lc)
      re_blocks[[length(re_blocks) + 1L]] <<- blk
    }
    re_descs[[length(re_descs) + 1L]] <<- list(
      arm = arm, var = var, levels = levels,
      n_groups = as.integer(n_groups), n_coefs = as.integer(n_coefs),
      coef_names = coef_names, correlated = isTRUE(correlated),
      has_intercept = identical(coef_names[[1L]], "(Intercept)"),
      coef_scales = if (is.null(coef_scales)) rep(1, n_coefs)
                    else as.numeric(coef_scales),
      block_start = b0 + 1L, n_blocks = length(re_blocks) - b0)
  }
  if (has_re) {
    add_re_term("psi", re_spec$group_idx, re_auto(dots$re.sigma.grid),
                NULL, re_spec$n_groups, re_spec$var %||% NA_character_,
                re_spec$levels, 1L, "(Intercept)", FALSE)
  }
  if (has_re_det) {
    for (d in re_det_terms) {
      add_re_term("p", d$codes, re_auto(dots$re.sigma.grid.p),
                  d$Z, d$n_groups, d$var, d$levels, d$n_coefs, d$coef_names,
                  d$correlated, logchol_grid = dots$re.logchol.grid.p,
                  coef_scales = d$coef_scales)
    }
  }
  if (has_re_pos) {
    for (d in re_pos_terms) {
      add_re_term("pos", d$codes, re_auto(dots$re.sigma.grid.pos),
                  d$Z, d$n_groups, d$var, d$levels, d$n_coefs, d$coef_names,
                  d$correlated, logchol_grid = dots$re.logchol.grid.pos,
                  coef_scales = d$coef_scales)
    }
  }

  # Arm-specific cover field blocks. Each field column of the `to = "positive"` bar
  # becomes ONE non-copied ICAR block scattering on the cover (pos) arm alone: the
  # psi + detection rows carry the 0-sentinel node so they skip it
  # (nested_laplace_joint_multi.h's `l_b > 0` guard), and it has no copy entry --
  # its amplitude is its OWN sigma (b<k>.sigma), decoupled from the occupancy
  # field's alpha copy. A trend (non-intercept) column carries its per-cell weight
  # on the pos rows. These trail the occupancy field blocks so the copy indices
  # (which name occupancy blocks only) stay valid; the RE blocks trail them.
  # `pos_field_specs` records each field block's arm + weight column so the
  # postprocess and the draw substrate map the blocks back to (occ vs pos)
  # amplitudes without re-deriving the layout. The non-copied ICAR block uses the
  # single-arm precision parameterization (axis b<k>.tau, sigma = 1/sqrt(tau)),
  # like the cover() arm-specific path -- the copy reparameterization (b<k>.sigma +
  # b<k>.alpha) applies only to a copied field. So the amplitude grid is passed as
  # tau = 1 / sigma^2. `sort()` / `as.numeric()` drop the auto-grid marker, so it
  # is re-applied on the translated tau vector; the source vector's own marker is
  # the provenance, which covers both the explicit pos-field grid and the shared
  # sigma grid it falls back to.
  pos_armspec_sigma_grid <- dots$sigma.grid.pos.field %||% sigma_grid
  pos_armspec_tau_grid   <- .tobs_mark_auto(
    sort(1.0 / as.numeric(pos_armspec_sigma_grid)^2),
    tulpa::is_auto_grid(pos_armspec_sigma_grid))

  # One arm-specific field -> ICAR block(s) that scatter on ONE arm's rows: the
  # node index lands in that arm's slot of the 3-slot spatial_idx (psi = site rows
  # [1], detection = visit rows [2], cover = pos rows [3]) and the other two carry
  # the 0-sentinel node so they skip it. Same shape for the cover (pos) and
  # detection (p) arms; only the target slot and the row->site map differ, so both
  # go through one builder (no per-arm copy of the block logic).
  arm_field_blocks <- function(af, arm) {
    idx_site <- as.integer(af$idx_obs)
    if (length(idx_site) != n_sites)
      stop(sprintf(paste0(
        "occu_cover(): the arm-specific %s field node index has %d values but ",
        "there are %d sites."),
        if (identical(arm, "pos")) "cover" else "detection",
        length(idx_site), n_sites), call. = FALSE)
    slot    <- if (identical(arm, "pos")) 3L else 2L
    row_map <- if (identical(arm, "pos")) pos_site else site_of_visit
    node    <- idx_site[row_map]
    zeros_i <- list(rep(0L, n_sites), rep(0L, n_v), rep(0L, n_pos_rows))
    zeros_w <- list(rep(0.0, n_sites), rep(0.0, n_v), rep(0.0, n_pos_rows))
    blocks <- list(); specs <- list()
    for (f in af$fields) {
      sidx <- zeros_i; sidx[[slot]] <- as.integer(node)
      blk <- list(
        type            = "icar",
        n_spatial_units = csr$n_spatial_units,
        adj_row_ptr     = csr$adj_row_ptr,
        adj_col_idx     = csr$adj_col_idx,
        n_neighbors     = csr$n_neighbors,
        tau_grid        = pos_armspec_tau_grid,
        spatial_idx     = sidx)
      if (!isTRUE(f$is_intercept)) {
        wt <- zeros_w; wt[[slot]] <- as.numeric(f$weight)[row_map]
        blk$svc_weight <- wt
      }
      blocks[[length(blocks) + 1L]] <- blk
      specs[[length(specs) + 1L]] <- list(
        arm = arm,
        weight = if (isTRUE(f$is_intercept)) NULL else f$column_name,
        is_intercept = isTRUE(f$is_intercept),
        column_name = f$column_name)
    }
    list(blocks = blocks, specs = specs)
  }

  pos_armspec_blocks <- list()
  pos_field_specs    <- list()
  for (as_arm in list(list(af = pos_armspec, arm = "pos"),
                      list(af = det_armspec, arm = "p"))) {
    if (is.null(as_arm$af)) next
    built <- arm_field_blocks(as_arm$af, as_arm$arm)
    pos_armspec_blocks <- c(pos_armspec_blocks, built$blocks)
    pos_field_specs    <- c(pos_field_specs, built$specs)
  }

  # Field-block descriptors in emitted (prior) order: the shared occupancy
  # intercept field, its coupled trend fields, then the arm-specific cover fields.
  # Consumed by the postprocess (per-block sigma naming, occ-vs-pos partition) and
  # the draw substrate (per-block occ / pos amplitude).
  field_specs <- c(
    list(list(arm = "shared", weight = NULL, is_intercept = TRUE)),
    lapply(coupled_trends, function(tf) list(
      arm = "shared", weight = tf$weight_label, is_intercept = FALSE)),
    pos_field_specs)

  # Pos-arm phi axis on the outer grid. For the latent path the pos arm's phi IS
  # sigma_u (the cover-latent SD), integrated over `sigma.u.grid` (default a
  # log-spaced grid around the between-unit init); the within-unit dispersion is
  # fixed in the spec. Otherwise the phi slot is sigma_pos and the optional
  # `phi.grid.pos` integrates it.
  if (is_latent) {
    su_grid <- dots$sigma.u.grid %||%
               (sigma_u_init * exp(seq(log(0.4), log(2.5), length.out = 4L)))
    phi_grid_arg <- list(pos = as.numeric(su_grid))
  } else {
    phi_grid_pos <- dots$phi.grid.pos
    phi_grid_arg <- if (!is.null(phi_grid_pos))
                      list(pos = as.numeric(phi_grid_pos))
                    else NULL
  }

  # Register the stateful latent spec for THIS fit: it captures the per-unit
  # detected cover values (indexed by pos-arm row, the order the builder emits)
  # and the fixed within-unit dispersion. Last-writer-wins under the fixed name;
  # the joint driver holds the resolved shared_ptr for the duration of the fit.
  if (is_latent) {
    n_quad_latent <- as.integer(dots$n.quad %||% .tobs_n_quad(
      if (is_beta) "cover_latent_beta" else "cover_latent_lognormal"))
    if (is_beta) {
      cpp_register_occu_cover_beta_latent_coupling(
        pos_cover_values, disp2_fixed, n_quad_latent)
    } else {
      cpp_register_occu_cover_lognormal_latent_coupling(
        pos_cover_values, disp2_fixed, n_quad_latent)
    }
  }

  if (correlated) {
    # Correlated free-Sigma MCAR: ONE coupled block over the bar's intercept +
    # coefficient fields, sharing a free cross-covariance Sigma (x) Q^-1 on the
    # occupancy arm (the within-arm relationship among the fields, integrated
    # over the outer CCD in log-Cholesky coords), copied onto the cover arm with
    # one amplitude alpha (the cross-arm transfer). The p (detection) arm carries
    # no field: its 0-sentinel cell index skips the block (mcar_block_factory's
    # `cell < 1 => skip`). Per (field, arm) weights mirror the trend path -- the
    # intercept is all-ones, each coefficient its per-site design column, and the
    # cover arm slices both by `pos_site`.
    field_weight_site <- c(
      list(rep(1.0, n_sites)),
      lapply(coupled_trends, function(tf) as.numeric(tf$weight))
    )
    p_mcar <- length(field_weight_site)
    pos_field_node <- as.integer(site_cell[pos_site])
    mcar_block <- list(
      type            = "mcar",
      n_spatial_units = csr$n_spatial_units,
      n_fields        = as.integer(p_mcar),
      adj_row_ptr     = csr$adj_row_ptr,
      adj_col_idx     = csr$adj_col_idx,
      n_neighbors     = csr$n_neighbors,
      spatial_idx     = list(as.integer(site_cell), rep(0L, n_v),
                             pos_field_node),
      field_weight    = lapply(field_weight_site, function(w)
        list(as.numeric(w), rep(1.0, n_v), as.numeric(w[pos_site])))
    )
    prior_arg <- list(mcar_block)
    copy_arg  <- .tobs_alpha_copy_spec("pos", 1L, alpha_axis)
    # The MCAR block carries p(p+1)/2 + 1 latent axes (log-Cholesky Sigma +
    # alpha), so the outer grid uses the mode-centred CCD by default rather than
    # a dense tensor (the same recipe the cover-hurdle MCAR path uses).
    if (is.null(dots$integration)) dots$integration <- "ccd"
  } else if (has_trend) {
    # Multi-block path: the intercept ICAR block plus one ICAR block per coupled
    # trend field, all on the same graph and each copied onto the pos arm with
    # its own alpha axis. The p arm is excluded from every field via its
    # field_coef = 0. Per-block svc_weight injects the per-row field weight on
    # the psi (per-cell) and pos (per-visit) arms; the p-arm weight is
    # irrelevant (field_coef = 0 already zeroes the p field).
    # Field node per arm row: psi rows are sites (-> site_cell), p rows are
    # visits (-> cell_of_visit), pos rows are either visits (per-visit cover) or
    # aggregated occupancy units (cell-aggregated cover); `pos_site` is the site
    # behind each pos row either way, so its field node is site_cell[pos_site]
    # and its SVC weight is w_psi[pos_site]. Under per-visit cover pos_site ==
    # site_of_visit, so this reduces to the previous cell_of_visit / w_visit.
    pos_field_node <- as.integer(site_cell[pos_site])
    # When the detection arm carries its own non-copied block -- an RE or an
    # arm-specific field -- its field_coef is 1 so that block scatters, so the
    # shared field must be skipped on detection by the 0-node sentinel rather
    # than by field_coef = 0.
    det_field_node <- if (has_re_det || has_det_armspec) rep(0L, n_v)
                      else cell_of_visit
    spatial_idx_arms <- list(as.integer(site_cell), det_field_node, pos_field_node)
    make_block <- function(weight_site) {
      w_psi <- if (is.null(weight_site)) rep(1.0, n_sites)
               else as.numeric(weight_site)
      w_pos <- w_psi[pos_site]
      icar_template(list(
        spatial_idx = spatial_idx_arms,
        svc_weight  = list(w_psi, rep(1.0, n_v), w_pos)
      ))
    }
    alpha_axis_trend <- .tobs_alpha_axis_trend(dots, alpha_axis)
    prior_arg <- c(
      list(make_block(NULL)),
      lapply(coupled_trends, function(tf) make_block(tf$weight))
    )
    copy_arg <- c(
      list(.tobs_alpha_copy_spec("pos", 1L, alpha_axis)),
      lapply(seq_len(n_trend), function(j)
        .tobs_alpha_copy_spec("pos", j + 1L, alpha_axis_trend))
    )
    # The arm-specific cover fields (non-copied, pos arm only) trail the occupancy
    # field blocks, then the RE blocks; the copy indices above name occupancy
    # blocks only, so they stay valid.
    prior_arg <- c(prior_arg, pos_armspec_blocks, re_blocks)
  } else if (has_any_re || has_armspec) {
    # Single shared occupancy field + per-group REs and/or arm-specific cover
    # fields: the multi-block driver with the occupancy field as block 1 (alpha
    # copy onto cover), then the non-copied arm-specific cover block(s), then the
    # iid RE block(s) -- each rides its own arm with no copy.
    field_block <- icar_template(list(
      spatial_idx = lapply(responses, function(a) as.integer(a$spatial_idx))))
    prior_arg <- c(list(field_block), pos_armspec_blocks, re_blocks)
    copy_arg  <- .tobs_alpha_copy_spec("pos", 1L, alpha_axis)
  } else if (isTRUE(.batch_collect)) {
    # Single-field, batched fused path: run the MULTI-block driver so the alpha
    # axis is an explicit copy spec and the per-arm field-node map is an explicit
    # `spatial_idx` (the single-block backend derives both from the pos arm's
    # field_coef; the multi-block driver needs them spelled out). A multi-block
    # fit with this copy is bit-identical to the single-block fit at a shared
    # grid (dev_notes/_probe_mb_vs_sb_occucover.R).
    prior_arg <- icar_template(list(
      spatial_idx = lapply(responses, function(a) as.integer(a$spatial_idx))))
    copy_arg  <- .tobs_alpha_copy_spec("pos", 1L, alpha_axis)
  } else {
    prior_arg <- icar_template()
    copy_arg  <- NULL
  }

  fit_call <- list(
    responses     = responses,
    prior         = prior_arg,
    phi_grid      = phi_grid_arg,
    cell_coupling = spec_name,
    control = c(list(
      max_iter  = as.integer(max.iter),
      tol       = as.numeric(tol),
      n_threads = as.integer(dots$n.threads %||% 1L),
      store_Q   = TRUE,
      # Inner-Newton curvature. The beta positive arm's observed mixture Hessian
      # is indefinite away from the mode, so observed- curvature Newton steps
      # stall and the inner Newton hits max.iter in every grid cell
      # (non-convergence -- the dominant cost). Expected/Fisher information is
      # PSD by construction, so the steps are well-conditioned and the inner
      # Newton converges in ~12 steps instead. The reported SEs, log_det and grid
      # weights are unchanged: the final mode-pass always re-factorizes with the
      # observed Hessian; the curvature mode only steers the path to the mode.
      # The lognormal arm is exactly quadratic (one inner step), so observed
      # curvature is already optimal -> keep "lm".
      hessian   = dots$hessian %||% (if (is_beta) "fisher" else "lm"),
      # Cholesky factor reuse (Shamanskii / chord) is exposed but defaults off for
      # the grid fit. Reuse also makes the off-factor scatter `grad_only` (skipping
      # the beta Hessian fill, the dominant per-iteration cost --
      # dev_notes/_profile_pareto_k.R), so it is NOT just a factorize saving; the
      # diagnostic re-solves enable it (refresh 4) because they need only the
      # converged log-marginal. The grid fit keeps refresh 1 by default since its
      # SEs/log-det use the true per-iteration curvature; raise via
      # `control$inner.refresh` if a grid fit is scatter-bound and SEs allow it.
      inner_refresh = as.integer(dots$inner.refresh %||% 1L),
      # Outer-grid parallelism (lever 2). The cover hurdle's large spatial
      # field takes the sparse inner-solve path, whose outer grid now runs
      # across `n.threads.outer` threads (per-thread Hessian builder / scratch
      # / specs). Defaults to serial; set it for the full-field fits.
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
      # Var-of-means consistency pass (tulpa engine, defaults ON in the joint
      # path) refines a sharply peaked axis post-integration -- independent of
      # adaptive_grid. Exposed so a fit can request a genuinely fixed outer grid
      # (`adaptive.grid = FALSE` AND `var.of.means.consistency = FALSE`), which
      # the fused batch driver requires for per-species bit-identity.
      var_of_means_consistency  = dots$var.of.means.consistency  %||% TRUE,
      var_of_means_min_ess      = dots$var.of.means.min.ess,
      # Cheap-pass screening of the outer grid. Every cell of a copy fit's
      # grid costs a full inner Newton on the areal field, while the posterior
      # mass sits on a handful of them, so `prune = TRUE` first sweeps the
      # lattice with a short warm-started Newton per cell and full-solves only
      # the cells whose screened weight clears `prune.tol`. The engine's safety
      # gate re-solves the full grid whenever the cheap ranking disagrees with
      # the full-solve argmax, so a pruned posterior is never silently wrong --
      # a screen that does not pay costs time, not accuracy. NULL leaves the
      # engine default (off).
      #
      # `[[` (exact) not `$`: `prune` is a unique prefix of `prune.tol`, so
      # `dots$prune` reads the TOLERANCE on a fit that sets only the tolerance.
      prune     = dots[["prune"]],
      prune_tol = dots[["prune.tol"]],
      # Outer-grid placement. `auto.recenter` moves a default axis onto the
      # hyperparameter mode and refits, which costs a second full solve of the
      # grid; `recenter.pilot` detects that placement on a THINNED grid instead,
      # so the full grid is solved once, at the placed axes. The joint path
      # recentres on the whole grid's collapsed-edge regime rather than a
      # per-axis rail, so `auto.recenter` takes TRUE / FALSE here and the
      # per-axis policy names the standalone path accepts are refused.
      auto_recenter   = dots[["auto.recenter"]],
      recenter_pilot  = dots[["recenter.pilot"]],
      # Prior on the cross-arm copy scale. The engine defaults to an exponential
      # continuum ("exponential") plus a point mass at alpha = 0 carrying half
      # the prior. `copy.slab = "flat"` makes the continuum flat in log alpha
      # over the span `alpha.grid` declares, and `copy.atom.mass = 0` drops the
      # point mass; NULL leaves both at the engine defaults.
      copy_slab       = dots$copy.slab,
      copy_atom_mass  = dots$copy.atom.mass,
      # Outer Pareto-k-hat accuracy diagnostic defaults OFF. It draws `k_samples`
      # extra hyperparameter points and re-solves the inner Laplace at each on the
      # full areal field, so it dominates the runtime -- ~200 re-solves vs the
      # grid's ~30-70 (measured 84-98% of wall time across field sizes). Per-phase
      # profiling (dev_notes/_profile_pareto_k.R) shows the binding per-solve cost
      # is the per-Newton-iteration Hessian/gradient SCATTER (the beta arm's
      # per-observation digamma/trigamma fill, 73-83%), NOT the sparse Cholesky
      # factorize (a flat ~0.5 ms, 8-12%, not super-linear up to ~1100 cells). cut
      # the diagnostic 2-4x (Shamanskii reuse + loosened inner tol + near-neighbour
      # batch order) with the k-hat byte-stable, but it stays OFF by default: it
      # reports k-hat only -- it does not move the betas / SDs / field -- so it is
      # an opt-in validation pass, matching the occu_joint path. Set
      # control$diagnose.k = TRUE to compute it (control$diagnose.draws sizes the
      # importance batch).
      diagnose_k = dots$diagnose.k %||% FALSE,
      # diagnose.draws is the diagnostic's precision knob (k.samples is the legacy
      # alias). The outer Pareto-k is scored ONCE over this many importance draws.
      k_samples = as.integer(dots$diagnose.draws %||% dots$k.samples %||% 500L),
      # Bootstrap outer Pareto-k uncertainty. The k-hat's sampling uncertainty is
      # bootstrapped from its raw importance log-ratios (k.bootstrap replicates, NO
      # new solves): reports the SE, 95% CI and the band_confident flag. A tighter
      # k needs more actual tail ratios -- raise diagnose.draws, NOT k.bootstrap.
      # k.tail.points (NULL = automatic PSIS rule) is an expert tail-threshold
      # control; k.conf.bands the reliability-band boundaries.
      k_bootstrap   = as.integer(dots$k.bootstrap %||% 1000L),
      k_tail_points = if (is.null(dots$k.tail.points)) NULL else as.integer(dots$k.tail.points),
      k_conf_bands  = dots$k.conf.bands %||% c(0.5, 0.7),
      # Diagnostic parallelism. When `diagnose.k = TRUE` the `k.samples` importance
      # re-solves are independent and run after the grid (every core free), each
      # solved single-threaded, so widening their outer pool is a bit-identical
      # wall-clock speedup. NULL (default) follows the fit's own thread grant
      # (`n.threads.outer` / inner `n.threads`); "auto" grabs the performance
      # cores; an integer pins the width. Forwarded verbatim.
      k_threads  = dots$k.threads,
      # Grid-cell checkpoint/resume. An EVA-scale occu_cover fit runs for
      # hours; `control$checkpoint = list(path =, resume =)` makes the outer
      # grid append each completed cell to `path` and a resume run load the
      # finished cells and solve only the rest, so a killed/rebooted fit
      # resumes instead of restarting. Forwarded verbatim to the engine.
      checkpoint = dots$checkpoint,
      # Outer-grid node layout. "ccd" places a central composite design over
      # the >= 3 latent axes (intercept + trend sigma/alpha) and crosses the
      # pos-arm phi tensor on top; "grid" forces the dense tensor. Forwarded
      # so a two-field trend fit can request CCD from the consumer side; NULL
      # falls through to the engine default.
      integration = dots$integration,
      # Outer-grid progress + ETA. Two channels, like the cover() hurdle, both
      # ON by default: `progress` gates the Rcout console progress bar -- ON by
      # default (NOT tied to `verbose`), set dots$progress = FALSE to silence
      # it; `progress.file` writes the ETA to disk and is emitted whenever it
      # is non-empty, independent of `progress`/`verbose` -- the channel that
      # survives a detached Start-Process stdout buffer. An explicit dotted key
      # overrides. `[[` (exact) not `$`: `dots$progress` prefix-matches
      # `progress.file`.
      progress          = dots[["progress"]] %||% TRUE,
      progress.every    = dots$progress.every,
      progress.throttle = dots$progress.throttle,
      progress.file     = dots$progress.file
    ),
    .tobs_adaptive_grid_control(dots))
  )
  if (!is.null(copy_arg)) fit_call$copy <- copy_arg
  # Regularizing hyperpriors on the outer grid, forwarded from control to the
  # joint driver's `prior_sigma` / `prior_alpha` / `prior_phi`. Each is a
  # list(<family>, <params>), e.g. list("pc.prec", c(U, alpha)) for a Penalized
  # Complexity prior (Simpson et al. 2017) on the spatial field SD, or
  # list("half_normal", scale). NULL (default) leaves the flat hyperprior. The
  # PC prior shrinks the field-SD upper tail toward the no-spatial base model
  # unless the data identifies a larger amplitude, so a weakly-identified field
  # is not driven to an inflated SD that widens every per-cell interval.
  if (!is.null(dots[["prior.sigma"]])) fit_call$prior_sigma <- dots[["prior.sigma"]]
  if (!is.null(dots[["prior.alpha"]])) fit_call$prior_alpha <- dots[["prior.alpha"]]
  if (!is.null(dots[["prior.phi"]]))   fit_call$prior_phi   <- dots[["prior.phi"]]

  ctx <- list(adj = adj, is_latent = is_latent, pi_list = pi_list,
              n_cells = n_cells,
              disp2_fixed   = if (is_latent) disp2_fixed   else NULL,
              n_quad_latent = if (is_latent) n_quad_latent else NULL,
              sigma_pos_init = sigma_pos_init, has_trend = has_trend,
              n_trend = n_trend, coupled_trends = coupled_trends, model = model,
              re_spec = re_spec,
              # Per-block RE descriptors in emitted (prior) order: each carries
              # the arm, grouping var + levels, group count, and coefficient
              # shape, so the postprocess maps each block back to its BLUP
              # columns, sigma axis, and per-arm summary.
              re_descs = re_descs,
              mcar = correlated,
              n_fields_mcar = if (correlated) 1L + n_trend else NULL,
              # Arm-specific cover field: `field_specs` labels every field block
              # (shared occupancy vs pos-arm), `n_occ_fields` is the occupancy field
              # count (intercept + trends) so the postprocess partitions the
              # trailing pos-arm blocks; `pos_field_specs` carries each pos block's
              # weight column for reporting.
              field_specs = field_specs,
              n_occ_fields = 1L + n_trend,
              has_pos_armspec = has_armspec,
              pos_field_specs = pos_field_specs,
              n_threads = as.integer(dots$n.threads.outer %||% 1L))

  # Batched fused path: return the assembled call + context instead of fitting,
  # so .tobs_fit_occu_cover_batch_fused can run B species through one fused
  # multi-block solve and post-process each with the shared ctx. Eligibility (no
  # pos-arm phi axis -> no latent / phi.grid.pos) is judged by the caller from
  # `fit_call$phi_grid` + `is_latent`.
  if (isTRUE(.batch_collect)) {
    return(structure(
      list(fit_call = fit_call, ctx = ctx, sigma_pos_init = sigma_pos_init,
           is_latent = is_latent, spec_name = spec_name, has_trend = has_trend),
      class = "occu_cover_jc_prep"))
  }

  fit <- do.call(tulpa::tulpa_nested_laplace_joint, fit_call)

  .occu_cover_jc_postprocess(fit, ctx)
}

