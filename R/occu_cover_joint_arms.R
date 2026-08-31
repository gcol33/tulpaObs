# =============================================================================
# occu_cover_joint.R
# - joint nested-Laplace path for occu_cover() driven by the cell-coupling
# registry ( consumer).
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
# via control$engine = "joint" -- same pattern as v2_joint.
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
# Detection-pattern compression (exact sufficient statistics). The all-undetected
# visits within a site enter the occupancy mixture only through
# prod_v (1 - p_v) = prod_u (1 - p_u)^{w_u}, so visits that share a detection
# design row (identical eta) collapse to one row of multiplicity w_u carried in
# the arm's n_trials slot; the cell-coupling kernel honours the weight exactly.
# Detected visits stay individual -- each carries its own per-visit cover, so the
# cover arm and its alignment with the detection arm are untouched. Returns the
# representative-visit indices (`sel`, one per compressed row, site-major) and the
# integer weight per row. Base R only (no data.table dependency); the grouping key
# is (site, full detection design row) compared bitwise, so it is exact.
.occu_cover_compress_nodet_visits <- function(site_of_visit, X_p, y_det) {
  n   <- length(site_of_visit)
  det <- y_det > 0.5
  ndi <- which(!det)                                   # non-detected visit rows
  if (length(ndi) == 0L)
    return(list(sel = seq_len(n), weight = rep(1L, n)))

  key <- cbind(site_of_visit[ndi], X_p[ndi, , drop = FALSE])   # (site, det design)
  o   <- do.call(order, lapply(seq_len(ncol(key)), function(j) key[, j]))
  ks  <- key[o, , drop = FALSE]
  if (nrow(ks) == 1L) {
    newgrp <- TRUE
  } else {
    diff_prev <- ks[-1L, , drop = FALSE] != ks[-nrow(ks), , drop = FALSE]
    newgrp    <- c(TRUE, rowSums(diff_prev) > 0L)       # exact bitwise change
  }
  grp      <- cumsum(newgrp)                            # group id in sorted order
  rep_sort <- which(!duplicated(grp))                  # first sorted row per group
  rep_orig <- ndi[o[rep_sort]]                          # representative, original index
  wt_grp   <- as.integer(tabulate(grp))                # group multiplicities

  sel    <- c(which(det), rep_orig)                    # detected (w=1) + nodet groups
  weight <- c(rep(1L, sum(det)), wt_grp)
  ord    <- order(site_of_visit[sel], sel)             # site-major, stable
  list(sel = sel[ord], weight = weight[ord])
}

.occu_cover_build_joint_arms <- function(model, sigma_pos_init,
                                                  alpha_axis,
                                                  positive = "lognormal",
                                                  multi = FALSE,
                                                  n_cells = NULL,
                                                  site_cell = NULL,
                                                  cover_aggregate = "none",
                                                  det_field = FALSE,
                                                  compress_nodet = FALSE) {
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
    stop("occu_cover joint: no valid visits in the data.",
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

  # Detection-pattern compression: collapse each site's all-undetected visits that
  # share a detection design row to one weighted row (weight in the p arm's
  # n_trials). Gated to the plain single-species per-visit path -- an obs-arm RE or
  # an arm-specific detection field makes eta depend on more than X_p (so identical
  # X_p rows are no longer exchangeable), and per-visit cover ("none") is required
  # so detected visits stay individual. Off (compress_nodet = FALSE) leaves every
  # arm byte-identical to the uncompressed build.
  do_compress <- isTRUE(compress_nodet) && is.null(model$re_det) &&
    is.null(model$re_pos) && !isTRUE(det_field) &&
    identical(cover_aggregate, "none")
  if (do_compress) {
    cmp           <- .occu_cover_compress_nodet_visits(site_of_visit, X_p,
                                                       y_det_visit)
    sel           <- cmp$sel
    visit_weight  <- cmp$weight
    site_of_visit <- site_of_visit[sel]
    cell_of_visit <- cell_of_visit[sel]
    X_p           <- X_p[sel, , drop = FALSE]
    X_pos         <- X_pos[sel, , drop = FALSE]
    y_det_visit   <- y_det_visit[sel]
    y_pos_visit   <- y_pos_visit[sel]
    n_visits_valid <- length(sel)
  } else {
    visit_weight  <- rep(1L, n_visits_valid)
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
  # field amplitude on EVERY block, so one scalar decouples the p arm from both
  # the intercept field and the trend field). spatial_idx is then a placeholder.
  # field_coef gates EVERY latent block on this arm (it is the per-arm arm_scale
  # multiplier), so it stays 0 to decouple the detection predictor from the
  # shared field -- UNLESS the detection arm carries its own non-copied block: a
  # random effect or an arm-specific spatial field. The shared field is
  # decoupled from detection by its `spatial_idx = 0` sentinel either way (the
  # engine skips a 0 node before any field indexing), so when the detection arm
  # carries its own block field_coef is 1 so that block scatters onto the
  # detection rows; the shared field's detection `spatial_idx` is forced to the
  # same sentinel.
  det_field_coef <- if (!is.null(model$re_det) || isTRUE(det_field)) 1.0 else 0
  arm_p <- list(
    y            = as.numeric(y_det_visit),
    n_trials     = as.integer(visit_weight),
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
  # positionally. "mean" / "median" collapse the cover arm to ONE row per
  # occupancy unit (site) that has any detection, carrying the mean / median
  # cover over that site's detected visits and the site-level positive design;
  # the `_agg` cell-coupling spec evaluates the cover density once per cell so
  # the cover arm contributes at the cell scale rather than the per-visit scale
  # (otherwise a cell with many detected plots drives the shared field far more
  # than the single occupancy observation for that cell).
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
      stop("occu_cover joint: no detected visits to aggregate cover ",
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

  # Observation-arm RE terms, aligned to the arm rows by the same `keep` the
  # detection arm uses. model$re_det / model$re_pos are per-term LISTS (one entry
  # per crossed / nested / slope term); each term's site-major group codes --
  # and, for a random slope, its per-row design `Z` (intercept + covariate
  # columns) -- are subset by `keep`. The detection arm is one row per valid
  # visit; per-visit cover (the only mode an obs-arm RE supports) is the same row
  # set, so both subset by `keep`.
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
    arm_pos$field_coef <- .tobs_alpha_field_coef(alpha_axis)
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


# Resolve per-arm weakly-informative fixed-effect priors for the coupled arms,
# returning list(psi=, p=, pos=) of list(mean, prec) (NULL per arm -> the weak
# engine default). All three arms carry weakly-informative defaults; the
# intercept priors are load-bearing in the shared-field path:
#   * The detection (p) intercept prior keeps the coupled occupancy mixture off
#     the psi = 1 boundary at weak detection (the logit-scale score vanishes as
#     psi -> 1, so an unpenalised intercept runs away).
#   * The cover (pos) intercept prior keeps the cover intercept off the
# field-level confound. The cover arm sees the shared field
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


