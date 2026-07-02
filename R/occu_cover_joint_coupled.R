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
                                                  site_cell = NULL,
                                                  cover_aggregate = "none",
                                                  det_field = FALSE) {
  n_sites    <- model$n_sites
  max_visits <- model$max_visits
  if (is.null(site_cell)) site_cell <- seq_len(n_sites)
  if (is.null(n_cells))   n_cells   <- max(site_cell)

  # Compact the visit-level data to one row per VALID visit, site-major (site 1's
  # visits in increasing visit index, then site 2's, ...). The compact (ragged)
  # model already carries exactly these rows in that order, so it is read
  # directly; the dense model flattens its [n_sites x max_visits] grid and keeps
  # the valid cells. Both yield identical `(y_det_visit, y_pos_visit,
  # site_of_visit)`, so the arms below are identical -- the only difference is
  # whether the padded grid was ever materialised.
  if (isTRUE(model$ragged)) {
    site_of_visit  <- model$site_of_visit
    y_det_visit    <- model$y_det_visit
    y_pos_visit    <- model$y_pos_visit
    n_visits_valid <- length(site_of_visit)
    keep           <- seq_len(n_visits_valid)   # X_*_visit already compacted
    cell_of_visit  <- as.integer(site_cell[site_of_visit])
  } else {
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
  }

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
  # field_coef gates EVERY latent block on this arm (it is the per-arm arm_scale
  # multiplier), so it stays 0 to decouple the detection predictor from the
  # shared field -- UNLESS the detection arm carries its own non-copied block: a
  # random effect (gcol33/tulpaObs#102) or an arm-specific spatial field
  # (gcol33/tulpa#140). The shared field is decoupled from detection by its
  # `spatial_idx = 0` sentinel either way (the engine skips a 0 node before any
  # field indexing), so when the detection arm carries its own block field_coef
  # is 1 so that block scatters onto the detection rows; the shared field's
  # detection `spatial_idx` is forced to the same sentinel.
  det_field_coef <- if (!is.null(model$re_det) || isTRUE(det_field)) 1.0 else 0
  arm_p <- list(
    y            = as.numeric(y_det_visit),
    n_trials     = rep(1L, n_visits_valid),
    X            = X_p,
    spatial_idx  = rep(0L, n_visits_valid),
    family       = "binomial",
    phi          = 1.0,
    field_coef   = det_field_coef,
    coupled      = TRUE,
    cell_obs_map = site_of_visit
  )

  # pos arm rows. `cover_aggregate = "none"` (per-visit) keeps one row per valid
  # detected visit, aligned with the p arm so the cell-coupling spec reads them
  # positionally. "mean" / "median" (tulpaObs#33) collapse the cover arm to ONE
  # row per occupancy unit (site) that has any detection, carrying the
  # mean / median cover over that site's detected visits and the site-level
  # positive design; the `_agg` cell-coupling spec evaluates the cover density
  # once per cell so the cover arm contributes at the cell scale rather than the
  # per-visit scale (otherwise a cell with many detected plots drives the shared
  # field far more than the single occupancy observation for that cell).
  pos_cover_values <- NULL
  if (identical(cover_aggregate, "none")) {
    pos_site  <- site_of_visit
    pos_cell  <- cell_of_visit
    y_pos_arm <- y_pos_visit
    X_pos_arm <- X_pos
  } else {
    # "mean" / "median" / "latent" all carry one pos row per detected occupancy
    # unit with the site-level positive design; they differ only in the cover
    # response. "latent" (per-unit cover RE integrated out) keeps every detected
    # visit's cover in `pos_cover_values` for the stateful latent spec and uses
    # the per-site mean only as the arm's placeholder y (the spec reads its own
    # captured data, not y(2, j)).
    units    <- .occu_cover_unit_cover(model)
    pos_site <- units$pos_site
    if (length(pos_site) == 0L) {
      stop("occu_cover joint_coupled: no detected visits to aggregate cover ",
           "over.", call. = FALSE)
    }
    if (identical(cover_aggregate, "latent")) {
      pos_cover_values <- units$vals
      y_pos_arm <- vapply(units$vals, mean, numeric(1))   # placeholder
    } else {
      aggfun <- if (identical(cover_aggregate, "median")) stats::median else mean
      y_pos_arm <- vapply(units$vals, function(v) as.numeric(aggfun(v)),
                          numeric(1))
    }
    X_pos_arm <- X_pos_site[pos_site, , drop = FALSE]   # site-level design only
    pos_cell  <- as.integer(site_cell[pos_site])
  }
  n_pos_rows <- length(pos_site)

  # Observation-arm RE terms (gcol33/tulpaObs#102, #103), aligned to the arm rows
  # by the same `keep` the detection arm uses. model$re_det / model$re_pos are
  # per-term LISTS (one entry per crossed / nested / slope term); each term's
  # site-major group codes -- and, for a random slope, its per-row design `Z`
  # (intercept + covariate columns) -- are subset by `keep`. The detection arm is
  # one row per valid visit; per-visit cover (the only mode an obs-arm RE
  # supports) is the same row set, so both subset by `keep`.
  keep_re_term <- function(d) c(d, list(codes = as.integer(d$codes_flat[keep]),
    Z = if (!is.null(d$Z)) d$Z[keep, , drop = FALSE] else NULL))
  re_det_terms <- if (!is.null(model$re_det)) lapply(model$re_det, keep_re_term)
                  else NULL
  re_pos_terms <- if (!is.null(model$re_pos) && identical(cover_aggregate, "none"))
                  lapply(model$re_pos, keep_re_term) else NULL

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
    y            = as.numeric(y_pos_arm),
    n_trials     = rep(1L, n_pos_rows),
    X            = X_pos_arm,
    spatial_idx  = pos_cell,
    family       = positive,
    phi          = sigma_pos_init,
    coupled      = TRUE,
    cell_obs_map = as.integer(pos_site)
  )
  if (!multi) {
    arm_pos$field_coef <- list(name = "alpha", grid = alpha_grid)
  }

  list(responses      = list(psi = arm_psi, p = arm_p, pos = arm_pos),
       site_of_visit  = site_of_visit,
       cell_of_visit  = cell_of_visit,
       n_visits_valid = n_visits_valid,
       pos_site       = as.integer(pos_site),
       n_pos_rows     = n_pos_rows,
       pos_cover_values = pos_cover_values,
       re_det_terms   = re_det_terms,
       re_pos_terms   = re_pos_terms)
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


# Resolve per-arm weakly-informative fixed-effect priors for the coupled arms,
# returning list(psi=, p=, pos=) of list(mean, prec) (NULL per arm -> the weak
# engine default). All three arms carry weakly-informative defaults; the
# intercept priors are load-bearing in the shared-field path:
#   * The detection (p) intercept prior keeps the coupled occupancy mixture off
#     the psi = 1 boundary at weak detection (the logit-scale score vanishes as
#     psi -> 1, so an unpenalised intercept runs away).
#   * The cover (pos) intercept prior keeps the cover intercept off the
#     field-level confound (tulpaObs#32). The cover arm sees the shared field
#     only at detected visits, so its intercept trades off against the field
#     level over those cells -- a direction the sum-to-zero field constraint
#     does not pin when low-occupancy regions carry no cover. Left at the
#     engine's flat 1e-4 default that intercept floats to a huge posterior SD
#     (occupancy stays tight: it is regularised and observes every cell), which
#     blows up predict()'s conditional cover via Jensen. The weakly-informative
#     cover_priors() default (pos_intercept sd 3) bounds it without biasing the
#     data-identified mode.
# `priors = FALSE` / "none" disables all three. A supplied occu_priors() /
# cover_priors() / list overrides the matching arm(s); the cover arm still gets
# the cover_priors() default unless a cover_priors object narrows it. Precisions
# are 1 / sd^2, floored at the engine's own weak default (1e-4) so an Inf-sd
# bucket reproduces the pre-existing weak ridge rather than dropping the diagonal.
.occu_cover_coupled_arm_priors <- function(priors, responses) {
  if (identical(priors, FALSE) || identical(priors, "none")) {
    return(list(psi = NULL, p = NULL, pos = NULL))
  }
  to_prec <- function(pr) {
    if (is.null(pr) || length(pr$sd) == 0L) return(NULL)
    list(mean = as.numeric(pr$mean), prec = pmax(1 / pr$sd^2, 1e-4))
  }
  occ_spec <- if (inherits(priors, "occu_priors")) priors
              else if (inherits(priors, "cover_priors")) occu_priors()
              else .resolve_occu_priors(priors)        # NULL / list -> defaults
  cover_spec <- if (inherits(priors, "cover_priors")) priors else cover_priors()

  list(
    psi = to_prec(.prior_for_submodel(occ_spec, "psi", colnames(responses$psi$X))),
    p   = to_prec(.prior_for_submodel(occ_spec, "p",   colnames(responses$p$X))),
    pos = to_prec(.cover_arm_prior(cover_spec, "pos",  colnames(responses$pos$X)))
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
  is_beta <- identical(model$positive, "beta")
  is_lnrm <- identical(model$positive, "lognormal")
  if (!is_beta && !is_lnrm) {
    stop("occu_cover() joint_coupled engine supports positive = ",
         "\"lognormal\" or \"beta\".", call. = FALSE)
  }
  # Cover-arm granularity. "mean" / "median" (tulpaObs#33) route through the
  # `_agg` spec (one log f_pos at the per-unit mean / median); "latent" routes
  # through the stateful `_latent` spec (a per-unit cover RE integrated out, one
  # marginal per unit); "none" is the per-visit spec.
  cover_aggregate <- model$cover_aggregate %||% "none"
  is_latent  <- identical(cover_aggregate, "latent")
  aggregated <- !identical(cover_aggregate, "none") && !is_latent
  spec_base  <- if (is_beta) "occu_cover_beta" else "occu_cover_lognormal"
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
      "occu_cover joint_coupled: site_cell must map %d sites into 1..%d ",
      "graph nodes."), n_sites, n_cells), call. = FALSE)
  }

  dots <- list(...)

  # Per-group random intercepts (gcol33/tulpaObs#56, #102). The occupancy-arm RE
  # (`re_spec`, one code per site) and the observation-arm REs (model$re_det /
  # model$re_pos, one code per detection / positive-cover row) each join the fit
  # as an `iid` prior block whose per-group latent rides ONE arm. Any RE block
  # forces the multi-block driver (the field amplitude becomes an explicit copy
  # spec). Not composed with the cover-latent RE (the latent spec carries its own
  # per-unit cover RE) nor with the batched fused path (one species at a time).
  has_re     <- !is.null(re_spec)            # occupancy (psi) arm
  has_re_det <- !is.null(model$re_det)       # detection (p) arm
  has_re_pos <- !is.null(model$re_pos)       # positive-cover arm
  has_any_re <- has_re || has_re_det || has_re_pos
  # Arm-specific cover field (gcol33/tulpaObs#110): an independent, non-copied
  # ICAR block on the cover (pos) arm alone, composed with the shared occupancy
  # field. Forces the multi-block driver (like a trend field / RE block). Not
  # composed with the latent cover RE, the correlated MCAR (gated at parse), or
  # the batched fused path (one species at a time, no extra block).
  # Arm-specific fields carry the detection (p) arm as well as the cover (pos) arm
  # (each an independent, non-copied ICAR block on that arm alone). Both force the
  # multi-block driver.
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

  # Correlated (`|`) free-Sigma MCAR field (gcol33/tulpaObs#63): one coupled
  # block over the bar's intercept + coefficient fields, copied onto the cover
  # arm with one amplitude alpha. Scoped to the standard (non-latent, unbatched)
  # path; the latent cover RE and the fused batch driver are not composed with it.
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
    det_mat   <- model$valid & (model$y == 1L)
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
      det_mat <- model$valid & (model$y == 1L)
      sw      <- which(rowSums(det_mat) > 0L)
      vapply(sw, function(i) as.numeric(aggfun(model$y_pos[i, det_mat[i, ]])),
             numeric(1))
    } else {
      model$y_pos[model$valid & model$y == 1L]
    }
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
  }

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
    multi           = has_trend || has_any_re || has_armspec,
    n_cells         = n_cells,
    site_cell       = site_cell,
    cover_aggregate = cover_aggregate,
    det_field       = has_det_armspec
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

  # Per-group RE blocks (gcol33/tulpaObs#56, #102, #103). Each random-intercept
  # term (one per arm for #56/#102; several per arm for crossed / nested #103)
  # contributes one `iid` prior block whose per-group latent rides THAT arm only:
  # obs_idx is the 3-element (psi, p, pos) list of per-row group codes, with the
  # targeted arm carrying the codes and the other two zeroed (0 = no RE for that
  # row, the engine's scatter skip). Each block's SD integrates on the outer grid
  # over its `sigma_grid`. They trail the field block(s), so the field copy
  # indices stay valid; the blocks carry no copy (each rides its own arm).
  # `re_descs` records each block's arm + grouping metadata in prior order so the
  # postprocess maps each block back to its `b<k>.sigma` axis and BLUP columns.
  re_grid_default <- exp(seq(log(0.05), log(2), length.out = 6L))
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
  # them (gcol33/tulpaObs#102, #103):
  #   * intercept (n_coefs == 1, !correlated): one scalar `iid` block.
  #   * uncorrelated slope (!correlated, Z present): one weighted `iid` block per
  #     coefficient -- svc_weight = that coefficient's design column (the
  #     intercept column is all-ones, so its block is the scalar iid; tulpa#114).
  #   * correlated slope: one multivariate `miid` block over a free cross-coef
  #     Sigma -- field_weight = the design columns (tulpa#114).
  # The descriptor records the [block_start, n_blocks] run, in emission order, so
  # the postprocess pulls the right latent columns; the blocks trail the field
  # block(s), so the field copy indices stay valid.
  add_re_term <- function(arm, codes, grid, Z, n_groups, var, levels,
                          n_coefs, coef_names, correlated, logchol_grid = NULL,
                          coef_scales = NULL) {
    obs_idx <- obs_for(arm, codes)
    b0 <- length(re_blocks)
    if (!isTRUE(correlated)) {
      for (cc in seq_len(n_coefs)) {
        blk <- list(type = "iid", n_units = as.integer(n_groups),
                    sigma_grid = as.numeric(grid), obs_idx = obs_idx)
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
      # shared field + copy under the engine's outer-grid cap (gcol33/tulpa#114).
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
    add_re_term("psi", re_spec$group_idx, dots$re.sigma.grid %||% re_grid_default,
                NULL, re_spec$n_groups, re_spec$var %||% NA_character_,
                re_spec$levels, 1L, "(Intercept)", FALSE)
  }
  if (has_re_det) {
    for (d in re_det_terms) {
      add_re_term("p", d$codes, dots$re.sigma.grid.p %||% re_grid_default,
                  d$Z, d$n_groups, d$var, d$levels, d$n_coefs, d$coef_names,
                  d$correlated, logchol_grid = dots$re.logchol.grid.p,
                  coef_scales = d$coef_scales)
    }
  }
  if (has_re_pos) {
    for (d in re_pos_terms) {
      add_re_term("pos", d$codes, dots$re.sigma.grid.pos %||% re_grid_default,
                  d$Z, d$n_groups, d$var, d$levels, d$n_coefs, d$coef_names,
                  d$correlated, logchol_grid = dots$re.logchol.grid.pos,
                  coef_scales = d$coef_scales)
    }
  }

  # Arm-specific cover field blocks (gcol33/tulpaObs#110). Each field column of the
  # `to = "positive"` bar becomes ONE non-copied ICAR block scattering on the
  # cover (pos) arm alone: the psi + detection rows carry the 0-sentinel node so
  # they skip it (nested_laplace_joint_multi.h's `l_b > 0` guard), and it has no
  # copy entry -- its amplitude is its OWN sigma (b<k>.sigma), decoupled from the
  # occupancy field's alpha copy. A trend (non-intercept) column carries its
  # per-cell weight on the pos rows. These trail the occupancy field blocks so the
  # copy indices (which name occupancy blocks only) stay valid; the RE blocks
  # trail them. `pos_field_specs` records each field block's arm + weight column so
  # the postprocess and the draw substrate map the blocks back to (occ vs pos)
  # amplitudes without re-deriving the layout.
  # The non-copied ICAR block uses the single-arm precision parameterization
  # (axis b<k>.tau, sigma = 1/sqrt(tau)), like the cover() arm-specific path -- the
  # copy reparameterization (b<k>.sigma + b<k>.alpha) applies only to a copied
  # field. So the amplitude grid is passed as tau = 1 / sigma^2.
  pos_armspec_sigma_grid <- dots$sigma.grid.pos.field %||% sigma_grid
  pos_armspec_tau_grid   <- sort(1.0 / as.numeric(pos_armspec_sigma_grid)^2)

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
    n_quad_latent <- as.integer(dots$n.quad %||% (if (is_beta) 15L else 1L))
    if (is_beta) {
      cpp_register_occu_cover_beta_latent_coupling(
        pos_cover_values, disp2_fixed, n_quad_latent)
    } else {
      cpp_register_occu_cover_lognormal_latent_coupling(
        pos_cover_values, disp2_fixed, n_quad_latent)
    }
  }

  if (correlated) {
    # Correlated free-Sigma MCAR (gcol33/tulpaObs#63): ONE coupled block over the
    # bar's intercept + coefficient fields, sharing a free cross-covariance
    # Sigma (x) Q^-1 on the occupancy arm (the within-arm relationship among the
    # fields, integrated over the outer CCD in log-Cholesky coords), copied onto
    # the cover arm with one amplitude alpha (the cross-arm transfer). The p
    # (detection) arm carries no field: its 0-sentinel cell index skips the block
    # (mcar_block_factory's `cell < 1 => skip`). Per (field, arm) weights mirror
    # the trend path -- the intercept is all-ones, each coefficient its per-site
    # design column, and the cover arm slices both by `pos_site`.
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
    copy_arg  <- list(arm = "pos", block = 1L, alpha_grid = alpha_grid)
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
    # When the detection arm carries its own non-copied block -- an RE
    # (gcol33/tulpaObs#102) or an arm-specific field (gcol33/tulpa#140) -- its
    # field_coef is 1 so that block scatters, so the shared field must be skipped
    # on detection by the 0-node sentinel rather than by field_coef = 0.
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
    # The arm-specific cover fields (non-copied, pos arm only) trail the occupancy
    # field blocks, then the RE blocks; the copy indices above name occupancy
    # blocks only, so they stay valid (gcol33/tulpaObs#110).
    prior_arg <- c(prior_arg, pos_armspec_blocks, re_blocks)
  } else if (has_any_re || has_armspec) {
    # Single shared occupancy field + per-group REs and/or arm-specific cover
    # fields: the multi-block driver with the occupancy field as block 1 (alpha
    # copy onto cover), then the non-copied arm-specific cover block(s), then the
    # iid RE block(s) -- each rides its own arm with no copy (gcol33/tulpaObs#110).
    field_block <- icar_template(list(
      spatial_idx = lapply(responses, function(a) as.integer(a$spatial_idx))))
    prior_arg <- c(list(field_block), pos_armspec_blocks, re_blocks)
    copy_arg  <- list(arm = "pos", block = 1L, alpha_grid = alpha_grid)
  } else if (isTRUE(.batch_collect)) {
    # Single-field, batched fused path: run the MULTI-block driver so the alpha
    # axis is an explicit copy spec and the per-arm field-node map is an explicit
    # `spatial_idx` (the single-block backend derives both from the pos arm's
    # field_coef; the multi-block driver needs them spelled out). A multi-block
    # fit with this copy is bit-identical to the single-block fit at a shared
    # grid (dev_notes/_probe_mb_vs_sb_occucover.R).
    prior_arg <- icar_template(list(
      spatial_idx = lapply(responses, function(a) as.integer(a$spatial_idx))))
    copy_arg  <- list(arm = "pos", block = 1L, alpha_grid = alpha_grid)
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
      # Cholesky factor reuse (Shamanskii / chord) is exposed but defaults off
      # for the grid fit. Reuse also makes the off-factor scatter `grad_only`
      # (skipping the beta Hessian fill, the dominant per-iteration cost --
      # dev_notes/_profile_pareto_k.R), so it is NOT just a factorize saving; the
      # tulpa#118 diagnostic re-solves enable it (refresh 4) because they need only
      # the converged log-marginal. The grid fit keeps refresh 1 by default since
      # its SEs/log-det use the true per-iteration curvature; raise via
      # `control$inner.refresh` if a grid fit is scatter-bound and SEs allow it.
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
      # Var-of-means consistency pass (tulpa engine, defaults ON in the joint
      # path) refines a sharply peaked axis post-integration -- independent of
      # adaptive_grid. Exposed so a fit can request a genuinely fixed outer grid
      # (`adaptive.grid = FALSE` AND `var.of.means.consistency = FALSE`), which
      # the fused batch driver requires for per-species bit-identity
      # (gcol33/tulpa#69, gcol33/tulpaObs#58).
      var_of_means_consistency  = dots$var.of.means.consistency  %||% TRUE,
      var_of_means_tolerance    = dots$var.of.means.tolerance    %||% 0.7,
      # Outer Pareto-k-hat accuracy diagnostic defaults OFF (gcol33/tulpaObs#101).
      # It draws `k_samples` extra hyperparameter points and re-solves the inner
      # Laplace at each on the full areal field, so it dominates the runtime --
      # ~200 re-solves vs the grid's ~30-70 (measured 84-98% of wall time across
      # field sizes). Per-phase profiling (dev_notes/_profile_pareto_k.R) shows the
      # binding per-solve cost is the per-Newton-iteration Hessian/gradient SCATTER
      # (the beta arm's per-observation digamma/trigamma fill, 73-83%), NOT the
      # sparse Cholesky factorize (a flat ~0.5 ms, 8-12%, not super-linear up to
      # ~1100 cells). tulpa#118 cut the diagnostic 2-4x (Shamanskii reuse + loosened
      # inner tol + near-neighbour batch order) with the k-hat byte-stable, but it
      # stays OFF by default: it reports k-hat only -- it does not move the betas /
      # SDs / field -- so it is an opt-in validation pass, matching the
      # occu_joint_coupled path. Set control$diagnose.k = TRUE to compute it
      # (control$diagnose.draws sizes the importance batch).
      diagnose_k = dots$diagnose.k %||% FALSE,
      # diagnose.draws is the diagnostic's precision knob (k.samples is the legacy
      # alias). The outer Pareto-k is scored ONCE over this many importance draws.
      diagnose_draws = as.integer(dots$diagnose.draws %||% dots$k.samples %||% 500L),
      # Bootstrap outer Pareto-k uncertainty (gcol33/tulpa#127). The k-hat's
      # sampling uncertainty is bootstrapped from its raw importance log-ratios
      # (k.bootstrap replicates, NO new solves): reports the SE, 95% CI and the
      # band_confident flag. A tighter k needs more actual tail ratios -- raise
      # diagnose.draws, NOT k.bootstrap. k.tail.points (NULL = automatic PSIS rule)
      # is an expert tail-threshold control; k.conf.bands the reliability-band
      # boundaries.
      k_bootstrap   = as.integer(dots$k.bootstrap %||% 1000L),
      k_tail_points = if (is.null(dots$k.tail.points)) NULL else as.integer(dots$k.tail.points),
      k_conf_bands  = dots$k.conf.bands %||% c(0.5, 0.7),
      # Diagnostic parallelism (gcol33/tulpa#117). When `diagnose.k = TRUE` the
      # `k.samples` importance re-solves are independent and run after the grid
      # (every core free), each solved single-threaded, so widening their outer
      # pool is a bit-identical wall-clock speedup. NULL (default) follows the
      # fit's own thread grant (`n.threads.outer` / inner `n.threads`); "auto"
      # grabs the performance cores; an integer pins the width. Forwarded verbatim.
      k_threads  = dots$k.threads,
      # Grid-cell checkpoint/resume (gcol33/tulpa#50). An EVA-scale occu_cover
      # fit runs for hours; `control$checkpoint = list(path =, resume =)` makes
      # the outer grid append each completed cell to `path` and a resume run
      # load the finished cells and solve only the rest, so a killed/rebooted
      # fit resumes instead of restarting. Forwarded verbatim to the engine.
      checkpoint = dots$checkpoint,
      # Outer-grid node layout (gcol33/tulpa#61, tulpaObs#31). "ccd" places a
      # central composite design over the >= 3 latent axes (intercept + trend
      # sigma/alpha) and crosses the pos-arm phi tensor on top; "grid" forces
      # the dense tensor. Forwarded so a two-field trend fit can request CCD
      # from the consumer side; NULL falls through to the engine default.
      integration = dots$integration,
      # Outer-grid progress + ETA (gcol33/tulpa#45, tulpaObs#43). Two channels,
      # like the cover() hurdle, both ON by default: `progress` gates the Rcout
      # console progress bar -- ON by default (NOT tied to `verbose`), set
      # dots$progress = FALSE to silence it; `progress.file` writes the ETA to
      # disk and is emitted whenever it is non-empty, independent of
      # `progress`/`verbose` -- the channel that survives a detached
      # Start-Process stdout buffer. An explicit dotted key overrides.
      # `[[` (exact) not `$`: `dots$progress` prefix-matches `progress.file`.
      progress          = dots[["progress"]] %||% TRUE,
      progress.every    = dots$progress.every,
      progress.throttle = dots$progress.throttle,
      progress.file     = dots$progress.file
    )
  )
  if (!is.null(copy_arg)) fit_call$copy <- copy_arg

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
              # columns, sigma axis, and per-arm summary (gcol33/tulpaObs#103).
              re_descs = re_descs,
              mcar = correlated,
              n_fields_mcar = if (correlated) 1L + n_trend else NULL,
              # Arm-specific cover field (gcol33/tulpaObs#110): `field_specs` labels
              # every field block (shared occupancy vs pos-arm), `n_occ_fields` is
              # the occupancy field count (intercept + trends) so the postprocess
              # partitions the trailing pos-arm blocks; `pos_field_specs` carries
              # each pos block's weight column for reporting.
              field_specs = field_specs,
              n_occ_fields = 1L + n_trend,
              has_pos_armspec = has_armspec,
              pos_field_specs = pos_field_specs,
              n_threads = as.integer(dots$n.threads.outer %||% 1L))

  # Batched fused path (gcol33/tulpa#66): return the assembled call + context
  # instead of fitting, so .tobs_fit_occu_cover_batch_fused can run B species
  # through one fused multi-block solve and post-process each with the shared
  # ctx. Eligibility (no pos-arm phi axis -> no latent / phi.grid.pos) is judged
  # by the caller from `fit_call$phi_grid` + `is_latent`.
  if (isTRUE(.batch_collect)) {
    return(structure(
      list(fit_call = fit_call, ctx = ctx, sigma_pos_init = sigma_pos_init,
           is_latent = is_latent, spec_name = spec_name, has_trend = has_trend),
      class = "occu_cover_jc_prep"))
  }

  fit <- do.call(tulpa::tulpa_nested_laplace_joint, fit_call)

  .occu_cover_jc_postprocess(fit, ctx)
}

# Compact-but-PRINCIPLED outer-grid for a correlated-slope `miid` block's free
# Sigma, in the engine's column-major lower-triangular log-Cholesky coordinates
# (gcol33/tulpa#114). A p = 2 (intercept + one slope) block -- the common
# `(1 + x | g)` -- gets a (sigma_0, sigma_1, rho) tensor sized to compose with the
# shared field + copy amplitude under the engine's outer-grid cap. The grid is a
# coarsened version of the engine's `.mcar_default_logchol_grid`: SYMMETRIC
# correlation nodes that include 0 and reach strong +/- (so the marginal
# correlation is not forced into a lop-sided range), and log-spaced SD nodes
# spanning small to large. The slope covariate is standardized
# (.occu_cover_obs_re_design), so a fixed SD bracket is meaningful for any
# covariate scale. Users widen it via `control$re.logchol.grid.p` /
# `re.logchol.grid.pos`. For p != 2 (multi-slope correlated, rare) return NULL so
# the engine fills its own default; that design is grid-heavy and usually needs
# an explicit coarse grid.
.occu_cover_miid_logchol_grid <- function(p, sig_grid = NULL, rho_grid = NULL) {
  if (!identical(as.integer(p), 2L)) return(NULL)
  sig_grid <- sig_grid %||% c(0.35, 0.8, 1.6)
  rho_grid <- rho_grid %||% c(-0.7, -0.3, 0, 0.3, 0.7)
  g <- expand.grid(s1 = sig_grid, s2 = sig_grid, rho = rho_grid,
                   KEEP.OUT.ATTRS = FALSE)
  out <- cbind(L11 = log(g$s1),
               L21 = g$rho * g$s2,
               L22 = log(g$s2 * sqrt(1 - g$rho^2)))
  colnames(out) <- c("L11", "L21", "L22")
  as.matrix(out)
}

# Public hyperparameter name for each RE block, aligned with the descriptor
# list (gcol33/tulpaObs#103). A lone term on an arm keeps the legacy bare name
# (sigma_re / sigma_re_p / sigma_re_pos for psi / detection / positive cover);
# crossed / nested terms sharing an arm are disambiguated by the grouping var
# (sigma_re_p_<var>), so every block's variance gets a distinct, stable name.
.occu_cover_re_sigma_names <- function(re_descs) {
  if (!length(re_descs)) return(character(0))
  base <- c(psi = "sigma_re", p = "sigma_re_p", pos = "sigma_re_pos")
  arms <- vapply(re_descs, `[[`, character(1), "arm")
  counts <- table(arms)
  vapply(seq_along(re_descs), function(i) {
    arm <- re_descs[[i]]$arm
    nm  <- base[[arm]]
    if (counts[[arm]] > 1L) {
      var <- re_descs[[i]]$var
      nm  <- paste0(nm, "_", make.names(if (is.na(var)) as.character(i) else var))
    }
    nm
  }, character(1))
}

# Post-process an occu_cover joint-coupled engine fit into a tobs_fit. `fit` is
# the tulpa_nested_laplace_joint return (single-species) or a per-species slice
# of a batched fused fit assembled to the same shape (gcol33/tulpa#66); `ctx`
# carries the Part-A context the shaping needs. Marginalisation here is a
# weighted sum over outer-grid cells (order-invariant), so a fused fixed-grid
# slice and an adaptive single-species fit shape identically given the same
# cells.
.occu_cover_jc_postprocess <- function(fit, ctx) {
  adj            <- ctx$adj
  is_latent      <- ctx$is_latent
  pi_list        <- ctx$pi_list
  n_cells        <- ctx$n_cells
  disp2_fixed    <- ctx$disp2_fixed
  n_quad_latent  <- ctx$n_quad_latent
  sigma_pos_init <- ctx$sigma_pos_init
  has_trend      <- ctx$has_trend
  n_trend        <- ctx$n_trend
  coupled_trends <- ctx$coupled_trends
  model          <- ctx$model

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

  # Reconcile the engine's grid weights with the ok-cell weights when the engine
  # left none usable. tulpa_posterior_draws() (predict / WAIC grid sampling) reads
  # fit$weights; when some cells carry a non-finite log_marginal (a corner of the
  # grid where the inner Newton -- e.g. the beta latent spec's Gauss-Hermite arm
  # -- did not converge) an unguarded weight normalization upstream collapses
  # fit$weights to all NaN (gcol33/tulpa#65), so the sampler finds no
  # positive-weight cell. Fall back to the same pure softmax over
  # finite-log_marginal cells the reported posterior moments use, so predict() /
  # WAIC stay consistent with the point estimates. Untouched when the engine
  # weights are already usable (every finite-grid fit).
  if (!any(is.finite(fit$weights) & fit$weights > 0)) {
    w_full <- numeric(length(fit$log_marginal))
    w_full[ok_cells] <- w
    fit$weights <- w_full
  }

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
  mcar      <- isTRUE(ctx$mcar)
  # One latent field block per coupled spatial field. The multi-block layout
  # reports every ICAR/structured field's 0-based offset in `field_starts`;
  # the single-block joint layout reports the one field's offset in `phi_start`.
  # The correlated MCAR field is a SINGLE block of p sub-fields laid out
  # contiguously (slot a*n_cells + cell for sub-field a), so the decode treats it
  # as p sub-fields over one contiguous run -- every per-field summary below
  # (demeaning, z-tables, block slicing) is shared with the independent path.
  field_starts0 <- layout$field_starts %||% layout$phi_start
  if (mcar) {
    n_fields  <- as.integer(ctx$n_fields_mcar)
    field_idx <- as.integer(field_starts0[[1L]] + seq_len(n_fields * n_cells))
  } else {
    n_fields  <- length(field_starts0)
    field_idx <- as.integer(unlist(lapply(field_starts0,
                                          function(s0) s0 + seq_len(n_cells))))
  }
  idx_joint <- c(bpsi_idx, bp_idx, bpos_idx, field_idx)
  blocks    <- .joint_inner_vcov_block(fit, idx_joint, n_dense = p_beta,
                                       n_threads = ctx$n_threads)

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

  # Per-group RE BLUPs (gcol33/tulpaObs#56, #102, #103). The RE blocks trail the
  # n_fields field blocks, so term i's blocks sit at layout positions
  # n_fields + block_start .. (+ n_blocks - 1); each block's latent is a
  # contiguous run in `modes`. The per-(group, coefficient) posterior-mean offset
  # is the grid-weighted mean of those columns (centred per coefficient), the SD
  # their grid-weighted posterior SD. A random intercept / uncorrelated slope has
  # one `iid` block per coefficient (`block c` -> column `c`); a correlated slope
  # is one `miid` block whose latent is coefficient-major (field a, group g at
  # (a-1)*n_groups + g), reshaped to [n_groups x n_coefs]. `latent_idx` is the
  # coefficient-major column run (the predict draws reshape it identically). The
  # per-coefficient variance / correlation is filled from the hyper axes below.
  re_descs <- ctx$re_descs %||% list()
  re_sig_names <- .occu_cover_re_sigma_names(re_descs)
  re_terms <- vector("list", length(re_descs))
  if (length(re_descs) > 0L) {
    bstart <- layout$block_start
    bsize  <- layout$block_size
    for (i in seq_along(re_descs)) {
      d   <- re_descs[[i]]
      nc  <- d$n_coefs; ng <- d$n_groups
      lay <- n_fields + d$block_start + seq_len(d$n_blocks) - 1L
      if (is.null(bstart) || length(bstart) < max(lay)) next
      B_mean <- matrix(0, ng, nc); B_sd <- matrix(0, ng, nc); lat <- integer(0)
      grid_moments <- function(cols) {
        u_mod <- modes[, cols, drop = FALSE]
        u_hat <- as.numeric(crossprod(w, u_mod))
        list(mean = u_hat, var = as.numeric(crossprod(w, u_mod^2)) - u_hat^2)
      }
      if (!isTRUE(d$correlated)) {
        for (cc in seq_len(nc)) {
          cols <- bstart[lay[cc]] + seq_len(bsize[lay[cc]])   # length n_groups
          mom  <- grid_moments(cols)
          B_mean[, cc] <- mom$mean - mean(mom$mean)
          B_sd[, cc]   <- sqrt(pmax(mom$var, 0))
          lat <- c(lat, cols)
        }
      } else {
        cols <- bstart[lay[1L]] + seq_len(bsize[lay[1L]])     # n_coefs*n_groups
        mom  <- grid_moments(cols)
        Hm <- matrix(mom$mean, ng, nc); Hs <- matrix(sqrt(pmax(mom$var, 0)), ng, nc)
        for (cc in seq_len(nc)) B_mean[, cc] <- Hm[, cc] - mean(Hm[, cc])
        B_sd[] <- Hs
        lat <- cols
      }
      # Back-transform a slope coefficient's BLUP from the standardized covariate
      # the fit ran on to its natural units (`b_raw = b_std / scale`); the
      # intercept's scale is 1. The per-coefficient SD is rescaled below the same
      # way; correlation is scale-free.
      sc <- d$coef_scales %||% rep(1, nc)
      if (any(sc != 1)) {
        B_mean <- sweep(B_mean, 2L, sc, "/")
        B_sd   <- sweep(B_sd,   2L, sc, "/")
      }
      re_terms[[i]] <- c(d, list(blup_mat = B_mean, blup_sd_mat = B_sd,
                                 prior_pos = lay, latent_idx = as.integer(lat)))
    }
    re_terms <- Filter(Negate(is.null), re_terms)
  }

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
  hyper_vals  <- list()
  hyper_names <- character(0)
  pick <- function(name, public = name) {
    j <- match(name, tg_names)
    if (is.na(j)) return(invisible(NULL))
    vals <- as.numeric(tg_ok[, j])
    m <- sum(w * vals)
    v <- sum(w * vals^2) - m^2
    hyper_means[[public]] <<- m
    hyper_sds  [[public]] <<- sqrt(max(v, 0))
    hyper_vals [[public]] <<- vals
    hyper_names <<- c(hyper_names, public)
  }
  # On the latent path the pos arm's phi axis IS the cover-latent SD; surface it
  # as `sigma_u` rather than the engine's generic `phi_pos`.
  phi_pos_public <- if (is_latent) "sigma_u" else "phi_pos"
  # pick2 reads a multi-block axis column (`b<k>.<axis>`) under a public name.
  pick2 <- function(public, col) {
    j <- match(col, tg_names)
    if (is.na(j)) return(invisible(NULL))
    vals <- as.numeric(tg_ok[, j])
    m <- sum(w * vals)
    v <- sum(w * vals^2) - m^2
    hyper_means[[public]] <<- m
    hyper_sds  [[public]] <<- sqrt(max(v, 0))
    hyper_vals [[public]] <<- vals
    hyper_names <<- c(hyper_names, public)
  }
  # put_derived stores the grid-weighted moments of a DERIVED per-cell quantity
  # (e.g. a sigma / correlation reconstructed from log-Cholesky axes), so it is
  # marginalized over the joint posterior rather than plugged in at the mode.
  put_derived <- function(public, vals) {
    vals <- as.numeric(vals)
    mn <- sum(w * vals); vv <- sum(w * vals^2) - mn^2
    hyper_means[[public]] <<- mn
    hyper_sds  [[public]] <<- sqrt(max(vv, 0))
    hyper_vals [[public]] <<- vals
    hyper_names <<- c(hyper_names, public)
  }
  # Clean a coefficient name for use in a hyperparameter name: `(Intercept)` ->
  # `intercept`, other punctuation collapsed to `_`.
  .re_coef_tag <- function(x) {
    x <- gsub("\\(Intercept\\)", "intercept", x)
    gsub("(^_|_$)", "", gsub("[^A-Za-z0-9]+", "_", x))
  }
  # A per-group RE block also forces the multi-block driver (its iid block is an
  # extra prior block), so the field axes carry the `b<k>.` names even without a
  # trend field. Each RE block's variance is its `b<n_fields+i>.sigma` axis.
  has_re     <- !is.null(ctx$re_spec)
  has_any_re <- length(re_descs) > 0L
  # Arm-specific cover field blocks (gcol33/tulpaObs#110) trail the occupancy
  # field blocks: fields 1..n_occ_fields are the shared occupancy intercept +
  # trends (copied to cover with alpha), the rest are the non-copied pos-arm
  # fields. Each reports its own sigma (b<k>.sigma) with no alpha copy axis.
  has_pos_armspec <- isTRUE(ctx$has_pos_armspec)
  n_occ_fields    <- ctx$n_occ_fields %||% n_fields
  pos_field_specs <- ctx$pos_field_specs %||% list()
  if (mcar) {
    # Free-Sigma MCAR hyperparameters (gcol33/tulpaObs#63). Reconstruct Sigma per
    # outer-grid cell from the log-Cholesky axes b1.L<i><j>, derive each field SD
    # (sigma_mcar<a>, a = 1 intercept .. p) and each cross-correlation
    # (rho_mcar_<a><b>), and carry the per-cell values so their grid-weighted
    # moments AND their between-cell covariance with the betas / other hypers are
    # marginalized over the joint posterior, not plugged in at the mode (the
    # marginalize-derived-quantities rule). The whole correlated field's copy
    # amplitude onto the cover arm is the alpha_mcar grid axis.
    p_f  <- n_fields
    m_lc <- p_f * (p_f + 1L) / 2L
    axis_nm <- character(m_lc); tt <- 1L
    for (j in seq_len(p_f)) for (i in j:p_f) {
      axis_nm[tt] <- sprintf("b1.L%d%d", i, j); tt <- tt + 1L
    }
    lc_cols <- match(axis_nm, tg_names)
    n_ok    <- length(ok_cells)
    sd_mat  <- matrix(NA_real_, n_ok, p_f)
    n_rho   <- p_f * (p_f - 1L) / 2L
    rho_mat <- matrix(NA_real_, n_ok, n_rho)
    for (k in seq_len(n_ok)) {
      L   <- .cover_mcar_logchol_to_L(as.numeric(tg_ok[k, lc_cols]), p_f)
      Sig <- L %*% t(L)
      sds_k <- sqrt(pmax(diag(Sig), 0))
      sd_mat[k, ] <- sds_k
      cc <- 1L
      for (a in seq_len(p_f - 1L)) for (b in (a + 1L):p_f) {
        rho_mat[k, cc] <- Sig[a, b] / max(sds_k[a] * sds_k[b], 1e-12); cc <- cc + 1L
      }
    }
    for (a in seq_len(p_f)) put_derived(sprintf("sigma_mcar%d", a), sd_mat[, a])
    cc <- 1L
    for (a in seq_len(p_f - 1L)) for (b in (a + 1L):p_f) {
      put_derived(sprintf("rho_mcar_%d%d", a, b), rho_mat[, cc]); cc <- cc + 1L
    }
    pick2("alpha_mcar", "b1.alpha")
    pick("phi_pos", phi_pos_public)
  } else if (has_trend || has_any_re || has_pos_armspec) {
    # Multi-block: block 1 is the intercept field, blocks 2.. the trend fields,
    # then the RE block(s). A single trend field keeps the bare
    # sigma_trend/alpha_trend names; several are indexed (sigma_trend1, ...).
    # Each RE block's SD is reported by `re_sig_names[i]` (sigma_re / sigma_re_p /
    # sigma_re_pos for a lone term on an arm; suffixed by grouping var for crossed
    # / nested terms sharing an arm).
    pick2("sigma", "b1.sigma")
    pick2("alpha", "b1.alpha")
    for (j in seq_len(n_trend)) {
      suffix <- if (n_trend == 1L) "" else as.character(j)
      pick2(paste0("sigma_trend", suffix), sprintf("b%d.sigma", j + 1L))
      pick2(paste0("alpha_trend", suffix), sprintf("b%d.alpha", j + 1L))
    }
    # Arm-specific cover fields (gcol33/tulpaObs#110): blocks n_occ_fields+1 ..
    # n_fields, each a NON-copied ICAR with its own precision axis (b<k>.tau,
    # sigma = 1/sqrt(tau)) and NO alpha copy. Report the grid-weighted marginal SD
    # (marginalize-derived-quantities). A lone intercept field keeps the bare
    # `sigma_pos_field`; a covariate column is suffixed by its name.
    for (j in seq_along(pos_field_specs)) {
      spec  <- pos_field_specs[[j]]
      blk_k <- n_occ_fields + j
      # sigma_pos_field for the cover arm, sigma_p_field for the detection arm.
      base_nm <- paste0("sigma_", if (identical(spec$arm, "p")) "p" else "pos",
                        "_field")
      nm <- if (isTRUE(spec$is_intercept)) base_nm
            else paste0(base_nm, "_", .re_coef_tag(spec$column_name))
      tau_col <- sprintf("b%d.tau", blk_k)
      sig_col <- sprintf("b%d.sigma", blk_k)
      if (sig_col %in% tg_names) {
        pick2(nm, sig_col)
      } else if (tau_col %in% tg_names) {
        put_derived(nm, 1.0 / sqrt(as.numeric(tg_ok[, match(tau_col, tg_names)])))
      }
    }
    # Per-term RE variance components. An intercept / uncorrelated-slope term
    # has one `b<P>.sigma` axis per coefficient; a correlated-slope term has one
    # `miid` block whose log-Cholesky axes b<P>.L<ij> reconstruct a free Sigma,
    # marginalized to per-coefficient SDs + cross-correlations over the grid. The
    # per-coefficient SD (and correlation) are stored back on `re_terms` so the
    # fit summary keeps the structured covariance, not just the scalar names.
    # Back-transform a slope coefficient's reported SD from the standardized
    # covariate the fit ran on to its natural units (`sigma_raw = sigma_std /
    # scale`); the intercept scale is 1. Divides the stored hyper mean / SD / per-
    # cell values so the fit summary AND its vcov are on the natural scale.
    rescale_hyper <- function(nm, s) {
      if (s == 1 || is.null(hyper_means[[nm]])) return(invisible())
      hyper_means[[nm]] <<- hyper_means[[nm]] / s
      hyper_sds  [[nm]] <<- hyper_sds  [[nm]] / s
      hyper_vals [[nm]] <<- hyper_vals [[nm]] / s
    }
    for (i in seq_along(re_terms)) {
      trm  <- re_terms[[i]]
      base <- re_sig_names[i]
      nc   <- trm$n_coefs; cn <- trm$coef_names; PP <- trm$prior_pos
      sc   <- trm$coef_scales %||% rep(1, nc)
      tags <- vapply(cn, .re_coef_tag, character(1))
      sigma_vec <- stats::setNames(rep(NA_real_, nc), cn)
      cor_mat   <- NULL
      if (!isTRUE(trm$correlated)) {
        for (cc in seq_len(nc)) {
          nm <- if (nc == 1L) base else paste0(base, "_", tags[cc])
          pick2(nm, sprintf("b%d.sigma", PP[cc]))
          rescale_hyper(nm, sc[cc])
          sigma_vec[cc] <- hyper_means[[nm]] %||% NA_real_
        }
      } else {
        p_f <- nc; P1 <- PP[1L]
        axis_nm <- character(p_f * (p_f + 1L) / 2L); tt <- 1L
        for (jj in seq_len(p_f)) for (ii in jj:p_f) {
          axis_nm[tt] <- sprintf("b%d.L%d%d", P1, ii, jj); tt <- tt + 1L
        }
        lc_cols <- match(axis_nm, tg_names)
        n_ok2 <- length(ok_cells)
        sd_mat2 <- matrix(NA_real_, n_ok2, p_f)
        n_rho2  <- p_f * (p_f - 1L) / 2L
        rho_mat2 <- matrix(NA_real_, n_ok2, max(n_rho2, 1L))
        for (k in seq_len(n_ok2)) {
          L   <- .cover_mcar_logchol_to_L(as.numeric(tg_ok[k, lc_cols]), p_f)
          Sig <- L %*% t(L); sds_k <- sqrt(pmax(diag(Sig), 0))
          sd_mat2[k, ] <- sds_k; rr <- 1L
          for (a in seq_len(p_f - 1L)) for (b in (a + 1L):p_f) {
            rho_mat2[k, rr] <- Sig[a, b] / max(sds_k[a] * sds_k[b], 1e-12); rr <- rr + 1L
          }
        }
        cor_mat <- diag(nc)
        for (cc in seq_len(nc)) {
          nm <- paste0(base, "_", tags[cc])
          put_derived(nm, sd_mat2[, cc] / sc[cc])   # natural-scale slope SD
          sigma_vec[cc] <- hyper_means[[nm]]
        }
        cbase <- sub("^sigma", "cor", base); rr <- 1L
        for (a in seq_len(nc - 1L)) for (b in (a + 1L):nc) {
          nm <- paste0(cbase, "_", tags[a], "_", tags[b]); put_derived(nm, rho_mat2[, rr])
          cor_mat[a, b] <- cor_mat[b, a] <- hyper_means[[nm]]; rr <- rr + 1L
        }
      }
      re_terms[[i]]$sigma <- sigma_vec
      re_terms[[i]]$cor   <- cor_mat
    }
    pick("phi_pos", phi_pos_public)
  } else {
    pick("sigma"); pick("alpha"); pick("phi_pos", phi_pos_public)
  }
  if (length(hyper_names) > 0L) {
    means <- c(means, unlist(hyper_means)[hyper_names])
    sds   <- c(sds,   unlist(hyper_sds)[hyper_names])
    par_names <- c(par_names, hyper_names)
  }

  names(means) <- par_names
  names(sds)   <- par_names

  # Parameter-surface vcov by the law of total covariance over the outer grid.
  # The betas carry within + between (`beta_block`); the hyperparameters ARE the
  # grid coordinates, so within a cell their variance -- and their covariance
  # with the betas -- is exactly zero, leaving only the between (variance- and
  # covariance-of-modes) term. That between term is computed jointly over
  # (betas, hyperparameters) so the cross-covariance and the hyper-hyper
  # covariance are both retained rather than dropped to a diagonal block.
  n_par <- length(means)
  V <- matrix(0, n_par, n_par)
  V[seq_len(p_beta), seq_len(p_beta)] <- beta_block
  if (length(hyper_names) > 0L) {
    hyper_idx <- p_beta + seq_along(hyper_names)
    H <- do.call(cbind, hyper_vals[hyper_names])              # n_ok x n_hyper
    H_dm <- sweep(H, 2L, unlist(hyper_means)[hyper_names], "-")
    if (!is.null(Vj)) {
      # Joint per-grid covariance available: betas carry the real within+between
      # block, so the exact between cross- and hyper-covariance keep V PSD.
      beta_modes <- modes[, c(bpsi_idx, bp_idx, bpos_idx), drop = FALSE]
      B_dm   <- sweep(beta_modes, 2L, means[seq_len(p_beta)], "-")
      cross  <- crossprod(B_dm * w, H_dm)                     # p_beta x n_hyper
      hyhy   <- crossprod(H_dm * w, H_dm)                     # n_hyper x n_hyper
      hyhy   <- (hyhy + t(hyhy)) / 2
      V[seq_len(p_beta), hyper_idx] <- cross
      V[hyper_idx, seq_len(p_beta)] <- t(cross)
      V[hyper_idx, hyper_idx]       <- hyhy
    } else {
      # Fallback (older tulpa, no per-grid block): betas are a marginal-only
      # diagonal, so a beta-hyper cross block could break PSD. Keep the hyper
      # block diagonal, matching the betas' marginal treatment.
      diag(V)[hyper_idx] <- sds[hyper_idx]^2
    }
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

  # Occupancy fields are blocks 1..n_occ_fields (intercept then coupled trends);
  # the arm-specific cover fields (gcol33/tulpaObs#110) are the trailing blocks.
  field_intercept <- fblocks[[1L]]$mean
  field_table     <- field_z_table(fblocks[[1L]])

  trend_blocks <- if (n_occ_fields >= 2L) fblocks[2:n_occ_fields] else list()
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

  # Arm-specific cover field posteriors (gcol33/tulpaObs#110): the per-cell z
  # tables for the independent cover-arm field(s), surfaced separately from the
  # occupancy fields so the user can inspect the cover trend map. The first is the
  # intercept field (bare `pos_field` / `pos_field_table`); a covariate column is
  # keyed by its name.
  pos_field_blocks <- if (has_pos_armspec && n_fields > n_occ_fields)
                        fblocks[(n_occ_fields + 1L):n_fields] else list()
  pos_field_means  <- lapply(pos_field_blocks, function(b) b$mean)
  pos_field_tables <- lapply(pos_field_blocks, field_z_table)
  if (length(pos_field_specs) == length(pos_field_means)) {
    pos_field_labels <- vapply(pos_field_specs, function(s)
      if (isTRUE(s$is_intercept)) "(Intercept)" else s$column_name, character(1))
    names(pos_field_means)  <- pos_field_labels
    names(pos_field_tables) <- pos_field_labels
  }
  pos_field       <- if (length(pos_field_means))  pos_field_means[[1L]]  else NULL
  pos_field_table <- if (length(pos_field_tables)) pos_field_tables[[1L]] else NULL

  # Joint betas+field posterior for downstream derived-quantity prediction
  # (delta_p / delta_cover marginalized over the full correlated posterior).
  # `joint_means` carries every field in the same demeaned convention as
  # `spatial_field`; `joint_vcov` is the law-of-total-covariance Vj (NULL on
  # the older-tulpa diagonal fallback). Fields are stacked in block order
  # (occupancy intercept, occupancy trends, then arm-specific cover fields).
  n_occ_trend <- max(n_occ_fields - 1L, 0L)
  field_par_names <- unlist(c(
    list(paste0("field_", seq_len(n_cells))),
    lapply(seq_len(n_occ_trend), function(j) {
      suffix <- if (n_occ_trend == 1L) "" else as.character(j)
      paste0("trend_field", suffix, "_", seq_len(n_cells))
    }),
    lapply(seq_len(n_fields - n_occ_fields), function(j) {
      suffix <- if (n_fields - n_occ_fields == 1L) "" else as.character(j)
      paste0("pos_field", suffix, "_", seq_len(n_cells))
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

  # Record the fixed within-unit dispersion (sigma_eps / beta precision) so the
  # latent-path predict can reconstruct the marginal cover (it pairs with the
  # integrated sigma_u reported in `means`).
  if (is_latent) {
    model$cover_latent_disp2 <- disp2_fixed
    model$cover_latent_nquad <- n_quad_latent
  }
  # Record the pos-arm dispersion the spec held fixed (sigma_pos for non-latent;
  # the latent path integrates sigma_u on the grid instead). The pointwise
  # log-likelihood reads it to score the cover term at the fitted dispersion
  # rather than a bare unit default (gcol33/tulpaObs#34).
  if (!is_latent) model$cover_pos_disp <- sigma_pos_init

  # Spatial summary. The correlated MCAR field reports its per-field SDs
  # (sigma_mcar, intercept first) and cross-correlations (rho_mcar) alongside the
  # single copy amplitude (alpha_mcar); sigma_mean / alpha_mean stay populated
  # (intercept-field SD + copy) so the shared print / summary / predict layer
  # reads it without branching on the field structure.
  if (mcar) {
    rho_nms <- grep("^rho_mcar_", names(hyper_means), value = TRUE)
    spatial_summary <- list(
      type = "mcar", graph = adj,
      sigma_mean = unname(hyper_means["sigma_mcar1"]),
      alpha_mean = unname(hyper_means["alpha_mcar"]),
      sigma_trend_mean = if (n_fields >= 2L)
        unname(hyper_means["sigma_mcar2"]) else NULL,
      alpha_trend_mean = NULL,
      sigma_mcar = vapply(seq_len(n_fields), function(a)
        unname(hyper_means[[sprintf("sigma_mcar%d", a)]]), numeric(1)),
      rho_mcar   = if (length(rho_nms))
        unname(unlist(hyper_means[rho_nms])) else numeric(0),
      rho_mcar_names = rho_nms,
      alpha_mcar = unname(hyper_means["alpha_mcar"]))
  } else {
    spatial_summary <- list(
      type = "icar", graph = adj,
      sigma_mean = unname(hyper_means["sigma"]),
      alpha_mean = unname(hyper_means["alpha"]),
      sigma_trend_mean = if (has_trend)
        unname(hyper_means[if (n_trend == 1L) "sigma_trend"
                           else "sigma_trend1"]) else NULL,
      alpha_trend_mean = if (has_trend)
        unname(hyper_means[if (n_trend == 1L) "alpha_trend"
                           else "alpha_trend1"]) else NULL)
  }

  structure(c(list(
    draws        = draws,
    means        = means,
    sds          = sds,
    vcov         = V,
    n_samples    = n_draws,
    n_params     = length(means),
    log_prob     = rep(log_lik_val, n_draws),
    log_lik      = log_lik_val,
    N            = if (isTRUE(model$ragged)) model$n_visits_valid else sum(model$valid)),
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
    trend_weight  = if (has_trend) trend_labels[[1L]]
                    else if (length(pos_field_specs)) {
                      cols <- Filter(Negate(is.null),
                                     lapply(pos_field_specs, `[[`, "column_name"))
                      nonint <- Filter(function(s) !isTRUE(s$is_intercept),
                                       pos_field_specs)
                      if (length(nonint)) nonint[[1L]]$column_name else NULL
                    } else NULL,
    trend_weights = if (has_trend) trend_labels        else NULL,
    # Arm-specific cover field (gcol33/tulpaObs#110): `field_specs` labels every
    # field block (shared occupancy vs pos arm) + its weight column so the draw
    # substrate maps each block to (occ, pos) amplitudes; the pos-field tables are
    # the independent cover field posterior for user inspection.
    field_specs      = ctx$field_specs,
    has_pos_armspec  = has_pos_armspec,
    pos_field        = pos_field,
    pos_field_table  = pos_field_table,
    pos_fields       = if (length(pos_field_means))  pos_field_means  else NULL,
    pos_field_tables = if (length(pos_field_tables)) pos_field_tables else NULL,
    joint_par_names = joint_par_names,
    joint_means     = joint_means,
    joint_vcov      = Vj,
    method       = "joint_coupled",
    positive     = model$positive,
    # Per-term RE summaries (gcol33/tulpaObs#56, #102, #103): a flat list, one
    # entry per RE block, each carrying its arm, grouping var + observed levels,
    # variance component, centred per-group BLUPs + SDs, group count, and latent
    # column indices (for the marginalized predict draws). A lone term on an arm
    # is keyed by the arm ("psi" / "p" / "pos"); crossed / nested terms sharing an
    # arm are keyed "<arm>:<var>". ranef() stacks them; predict() sums every term
    # on the predicted arm.
    re           = if (length(re_terms)) {
                     arms <- vapply(re_terms, `[[`, character(1), "arm")
                     counts <- table(arms)
                     keys <- vapply(re_terms, function(t)
                       if (counts[[t$arm]] > 1L) paste0(t$arm, ":", t$var)
                       else t$arm, character(1))
                     stats::setNames(lapply(re_terms, function(t) {
                       # Intercept term -> per-group BLUP vector (back-compat);
                       # slope term -> [n_groups x n_coefs] matrix with coef-named
                       # columns, plus the slope covariate names predict() weights
                       # rows by. `sigma` is the per-coefficient SD; `cor` the free
                       # cross-coefficient correlation matrix (NULL when scalar).
                       nc <- t$n_coefs
                       blup <- if (nc == 1L) as.numeric(t$blup_mat[, 1L]) else {
                         m <- t$blup_mat; colnames(m) <- t$coef_names; m }
                       blup_sd <- if (nc == 1L) as.numeric(t$blup_sd_mat[, 1L]) else {
                         m <- t$blup_sd_mat; colnames(m) <- t$coef_names; m }
                       covnms <- if (isTRUE(t$has_intercept)) t$coef_names[-1L]
                                 else t$coef_names
                       list(arm        = t$arm,
                            var        = if (is.na(t$var)) NULL else t$var,
                            sigma      = t$sigma,
                            cor        = t$cor,
                            blup       = blup,
                            blup_sd    = blup_sd,
                            n_groups   = t$n_groups,
                            n_coefs    = nc,
                            coef_names = t$coef_names,
                            covariate_names = covnms,
                            coef_scales = t$coef_scales %||% rep(1, nc),
                            has_intercept   = isTRUE(t$has_intercept),
                            correlated = t$correlated,
                            levels     = t$levels,
                            latent_idx = t$latent_idx)
                     }), keys)
                   } else NULL,
    joint_fit    = fit,
    convergence  = list(converged = TRUE, n_iter = NA_integer_)
  )), class = c("tobs_fit", "tulpa_fit"))
}
