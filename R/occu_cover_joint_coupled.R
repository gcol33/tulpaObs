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
                                                  cover_aggregate = "none") {
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
       pos_cover_values = pos_cover_values)
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
    multi           = has_trend,
    n_cells         = n_cells,
    site_cell       = site_cell,
    cover_aggregate = cover_aggregate
  )
  responses      <- arms_out$responses
  site_of_visit  <- arms_out$site_of_visit
  cell_of_visit  <- arms_out$cell_of_visit
  n_v            <- arms_out$n_visits_valid
  pos_site       <- arms_out$pos_site
  n_pos_rows     <- arms_out$n_pos_rows
  pos_cover_values <- arms_out$pos_cover_values

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

  if (has_trend) {
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
    spatial_idx_arms <- list(as.integer(site_cell), cell_of_visit, pos_field_node)
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
      checkpoint = dots$checkpoint,
      # Outer-grid node layout (gcol33/tulpa#61, tulpaObs#31). "ccd" places a
      # central composite design over the >= 3 latent axes (intercept + trend
      # sigma/alpha) and crosses the pos-arm phi tensor on top; "grid" forces
      # the dense tensor. Forwarded so a two-field trend fit can request CCD
      # from the consumer side; NULL falls through to the engine default.
      integration = dots$integration
    )
  )
  if (!is.null(copy_arg)) fit_call$copy <- copy_arg

  ctx <- list(adj = adj, is_latent = is_latent, pi_list = pi_list,
              n_cells = n_cells,
              disp2_fixed   = if (is_latent) disp2_fixed   else NULL,
              n_quad_latent = if (is_latent) n_quad_latent else NULL,
              sigma_pos_init = sigma_pos_init, has_trend = has_trend,
              n_trend = n_trend, coupled_trends = coupled_trends, model = model)

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
  pick <- function(name, public = name) {
    j <- match(name, tg_names)
    if (is.na(j)) return(invisible(NULL))
    vals <- as.numeric(tg_ok[, j])
    m <- sum(w * vals)
    v <- sum(w * vals^2) - m^2
    hyper_means[[public]] <<- m
    hyper_sds  [[public]] <<- sqrt(max(v, 0))
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
