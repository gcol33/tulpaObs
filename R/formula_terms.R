# =============================================================================
# formula_terms.R — Structured terms for tobs() formulas
#
# tulpaObs specifies spatial / temporal / random-effect / latent structure
# *inside* the process formulas, the way lme4 (`(1|g)`), mgcv (`s()`), brms
# (`gp()`), and INLA (`f()`) do. Each special term is a constructor that
# returns a classed `tobs_*` spec. The parser (`.tobs_parse_formula`) detects
# these by name against the registry below and evaluates them with the model
# `data` columns in scope, so bare symbols (`gp(lon, lat)`, `re(observer)`)
# resolve to columns exactly like `s()` inside `gam()`.
#
# The term constructors are NOT exported. They are only meaningful inside a
# tobs() formula and are resolved through the registry, never called by the
# user directly. Process membership (which linear predictor an effect enters)
# is determined by *which process formula* the term appears in — there is no
# `shared` argument. A single realization is shared across processes with
# `copy("id")`.
# =============================================================================


# ---------------------------------------------------------------------------
# Internal: stamp class + id onto a spec list
# ---------------------------------------------------------------------------
.tobs_term <- function(spec, class, id = NULL, label = NULL) {
  spec$id    <- id
  spec$label <- label
  structure(spec, class = class)
}

# Collect a coords matrix from either an explicit `coords =` matrix or two or
# more bare coordinate columns passed positionally (resolved to vectors).
.tobs_collect_coords <- function(dots, coords, term) {
  if (!is.null(coords)) {
    if (!is.matrix(coords) || ncol(coords) != 2L) {
      stop(sprintf("%s(): `coords` must be a matrix with 2 columns.", term),
           call. = FALSE)
    }
    return(coords)
  }
  if (length(dots) < 2L) {
    stop(sprintf(
      "%s(): supply two coordinate columns, e.g. %s(lon, lat), or coords = M.",
      term, term), call. = FALSE)
  }
  m <- tryCatch(do.call(cbind, dots[1:2]),
                error = function(e) stop(sprintf(
                  "%s(): could not bind coordinate columns: %s",
                  term, conditionMessage(e)), call. = FALSE))
  if (ncol(m) != 2L) {
    stop(sprintf("%s(): coordinates must resolve to 2 columns.", term),
         call. = FALSE)
  }
  m
}

# Resolve a grouping/time vector to 1-based integer codes plus a level count.
.tobs_index_codes <- function(x, term, arg) {
  if (is.factor(x) || is.character(x)) {
    f <- as.factor(x)
    return(list(idx = as.integer(f), n = nlevels(f)))
  }
  if (is.numeric(x)) {
    f <- as.factor(x)
    return(list(idx = as.integer(f), n = nlevels(f)))
  }
  stop(sprintf("%s(): `%s` must be a factor, character, or numeric column.",
               term, arg), call. = FALSE)
}

# Resolve an optional per-node SVC weight for an areal term. A non-NULL `weight`
# is a numeric column with one value per graph node; it turns the field into a
# spatially-varying coefficient whose contribution to a predictor row is
# `weight_i * amplitude * z[node_i]` (the areal analogue of INLA's
# f(node, weight, model = ...)). NULL is an unweighted (intercept) field. Only
# the occu_cover joint spatial path consumes the weight; the other areal
# consumers reject a weighted term via .tobs_reject_weighted_spatial().
.tobs_resolve_field_weight <- function(weight, n_units, term, per_obs = FALSE) {
  if (is.null(weight)) return(NULL)
  w <- tryCatch(as.numeric(weight),
                error = function(e) stop(sprintf(
                  "%s(): `weight` must be a numeric column.", term),
                  call. = FALSE))
  # With group_var the weight is per observation (one per site / data row),
  # not per graph node; its length is validated downstream against the site
  # count. Without group_var it is a per-node SVC covariate.
  if (!per_obs && length(w) != n_units) {
    stop(sprintf(
      "%s(): `weight` has length %d but the graph has %d node(s).",
      term, length(w), n_units), call. = FALSE)
  }
  if (any(!is.finite(w))) {
    stop(sprintf("%s(): `weight` has non-finite entries.", term), call. = FALSE)
  }
  w
}


# ---------------------------------------------------------------------------
# Spatial terms (areal)
# ---------------------------------------------------------------------------

# Areal terms accept an optional `group_var`: the name of a column that maps
# each observation to an areal unit (a graph node), for graphs defined over
# regions rather than over observations directly. When NULL the graph is over
# the observations 1:1. `group_var` is a deferred column-name string (matching
# tulpa's spatial_*() convention) — it is resolved by the engine, not here.
.tobs_check_graph <- function(graph, term) {
  if (!is.matrix(graph)) {
    stop(sprintf("%s(): `graph` must be an adjacency matrix.", term), call. = FALSE)
  }
  if (!isSymmetric(graph)) {
    stop(sprintf("%s(): `graph` must be symmetric.", term), call. = FALSE)
  }
}

# icar(graph)              — intrinsic CAR over an adjacency graph
#
# `weight` (optional) is a per-node numeric column that makes this a
# spatially-varying coefficient (a weighted field, `weight_i * z[node_i]`)
# instead of a plain intercept field; see .tobs_resolve_field_weight().
.tobs_term_icar <- function(graph, group_var = NULL, weight = NULL, id = NULL) {
  .tobs_check_graph(graph, "icar")
  csr <- adjacency_to_csr(graph)
  wlabel <- if (is.null(weight)) NULL else deparse(substitute(weight))
  weight <- .tobs_resolve_field_weight(weight, nrow(graph), "icar",
                                       per_obs = !is.null(group_var))
  .tobs_term(list(
    type = "icar", n_units = nrow(graph), graph = graph, group_var = group_var,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, weight = weight, weight_label = wlabel
  ), class = "tobs_spatial", id = id, label = "icar")
}

# bym2(graph, scale_factor) — BYM2 reparameterization of ICAR + IID
.tobs_term_bym2 <- function(graph, scale_factor = NULL, group_var = NULL,
                            weight = NULL, id = NULL) {
  .tobs_check_graph(graph, "bym2")
  csr <- adjacency_to_csr(graph)
  if (is.null(scale_factor)) scale_factor <- compute_bym2_scale(graph)
  wlabel <- if (is.null(weight)) NULL else deparse(substitute(weight))
  weight <- .tobs_resolve_field_weight(weight, nrow(graph), "bym2",
                                       per_obs = !is.null(group_var))
  .tobs_term(list(
    type = "bym2", n_units = nrow(graph), graph = graph, group_var = group_var,
    adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
    n_neighbors = csr$n_neighbors, scale_factor = scale_factor,
    weight = weight, weight_label = wlabel
  ), class = "tobs_spatial", id = id, label = "bym2")
}

# car(graph) / car_proper(graph) — (improper / proper) CAR areal fields.
# Both carry the graph for the cover-hurdle nested-Laplace engine, which
# precomputes its own adjacency via tulpa::spatial_car(). car_proper also
# carries the CSR adjacency so it can drive the occupancy / N-mixture
# multi-block nested-Laplace kernel (cpp_nested_laplace_multi reads the CSR
# directly when assembling the proper-CAR precision Q = tau (D - rho W)).
.tobs_term_car <- function(graph, group_var = NULL, id = NULL) {
  .tobs_check_graph(graph, "car")
  .tobs_term(list(type = "car", n_units = nrow(graph), graph = graph,
                  group_var = group_var),
             class = "tobs_spatial", id = id, label = "car")
}

.tobs_term_car_proper <- function(graph, group_var = NULL, weight = NULL,
                                  id = NULL) {
  .tobs_check_graph(graph, "car_proper")
  csr <- adjacency_to_csr(graph)
  wlabel <- if (is.null(weight)) NULL else deparse(substitute(weight))
  weight <- .tobs_resolve_field_weight(weight, nrow(graph), "car_proper")
  .tobs_term(list(type = "car_proper", n_units = nrow(graph), graph = graph,
                  group_var = group_var,
                  adj_row_ptr = csr$row_ptr, adj_col_idx = csr$col_idx,
                  n_neighbors = csr$n_neighbors,
                  weight = weight, weight_label = wlabel),
             class = "tobs_spatial", id = id, label = "car_proper")
}


# ---------------------------------------------------------------------------
# Spatial terms (continuous, NNGP / SPDE)
# ---------------------------------------------------------------------------

# gp(lon, lat, ...)        — NNGP-approximated Gaussian process
.tobs_term_gp <- function(..., coords = NULL, cov = "exponential", nu = 1.5,
                          nn = 15, id = NULL,
                          sigma2_prior_U = 1.0, sigma2_prior_alpha = 0.01,
                          phi_prior_lower = 0.01, phi_prior_upper = 10.0) {
  coords <- .tobs_collect_coords(list(...), coords, "gp")
  cov <- match.arg(cov, c("exponential", "matern", "gaussian", "spherical"))
  n <- nrow(coords); nn <- min(nn, n - 1L)
  nngp <- compute_nngp_neighbors(coords, nn)
  .tobs_term(list(
    type = "gp", n_obs = n, nn = nn,
    coords = as.vector(t(coords)),
    nn_idx = as.vector(t(nngp$nn_idx)),
    nn_dist = as.vector(t(nngp$nn_dist)),
    nn_neighbor_dist = as.vector(nngp$nn_neighbor_dist),
    nn_order = nngp$nn_order, nn_order_inv = nngp$nn_order_inv,
    cov_type = cov, nu = nu,
    sigma2_prior_U = sigma2_prior_U, sigma2_prior_alpha = sigma2_prior_alpha,
    phi_prior_lower = phi_prior_lower, phi_prior_upper = phi_prior_upper
  ), class = "tobs_spatial", id = id, label = "gp")
}

# multiscale_gp(lon, lat, ...) — two-scale (local + regional) NNGP
.tobs_term_multiscale_gp <- function(..., coords = NULL, cov = "exponential",
                                     nu = 1.5, nn_local = 15, nn_regional = 15,
                                     id = NULL,
                                     range_local_lower = 0.01, range_local_upper = 10.0,
                                     range_regional_lower = 0.01, range_regional_upper = 100.0,
                                     sigma2_local_prior_U = 1.0, sigma2_local_prior_alpha = 0.01,
                                     sigma2_regional_prior_U = 1.0, sigma2_regional_prior_alpha = 0.01) {
  coords <- .tobs_collect_coords(list(...), coords, "multiscale_gp")
  cov <- match.arg(cov, c("exponential", "matern", "gaussian", "spherical"))
  n <- nrow(coords)
  nn_local <- min(nn_local, n - 1L); nn_regional <- min(nn_regional, n - 1L)
  loc <- compute_nngp_neighbors(coords, nn_local)
  reg <- compute_nngp_neighbors(coords, nn_regional)
  .tobs_term(list(
    type = "multiscale_gp", n_obs = n, coords = as.vector(t(coords)),
    nn_local = nn_local,
    nn_idx_local = as.vector(t(loc$nn_idx)),
    nn_dist_local = as.vector(t(loc$nn_dist)),
    nn_neighbor_dist_local = as.vector(loc$nn_neighbor_dist),
    nn_order_local = loc$nn_order, nn_order_inv_local = loc$nn_order_inv,
    nn_regional = nn_regional,
    nn_idx_regional = as.vector(t(reg$nn_idx)),
    nn_dist_regional = as.vector(t(reg$nn_dist)),
    nn_neighbor_dist_regional = as.vector(reg$nn_neighbor_dist),
    nn_order_regional = reg$nn_order, nn_order_inv_regional = reg$nn_order_inv,
    cov_type = cov, nu = nu,
    range_local_lower = range_local_lower, range_local_upper = range_local_upper,
    range_regional_lower = range_regional_lower, range_regional_upper = range_regional_upper,
    sigma2_local_prior_U = sigma2_local_prior_U, sigma2_local_prior_alpha = sigma2_local_prior_alpha,
    sigma2_regional_prior_U = sigma2_regional_prior_U, sigma2_regional_prior_alpha = sigma2_regional_prior_alpha
  ), class = "tobs_spatial", id = id, label = "multiscale_gp")
}

# spde(lon, lat, ...)      — continuous Matern field via triangular mesh
.tobs_term_spde <- function(..., coords = NULL, mesh = NULL, max_edge = NULL,
                            cutoff = 0, nu = 1, id = NULL,
                            prior_range = c(0.5, 0.5), prior_sigma = c(1, 0.5)) {
  if (is.null(coords) && is.null(mesh)) {
    coords <- .tobs_collect_coords(list(...), coords, "spde")
  }
  tulpa_spec <- tulpa::spatial_spde(
    coords = coords, data = NULL, mesh = mesh,
    max_edge = max_edge, cutoff = cutoff, nu = nu,
    prior_range = prior_range, prior_sigma = prior_sigma
  )
  .tobs_term(list(
    type = "spde", tulpa_spec = tulpa_spec, n_units = tulpa_spec$n_mesh,
    nu = nu, prior_range = prior_range, prior_sigma = prior_sigma
  ), class = "tobs_spatial", id = id, label = "spde")
}


# ---------------------------------------------------------------------------
# Random effects, temporal structure, SVC, latent factors
# ---------------------------------------------------------------------------

# re(group, ...)           — grouped random effect (intercept / slope / iid)
#
# `covariate` may be a single column (one random slope) or several stacked
# columns (a multi-slope block): pass a numeric matrix (e.g. `cbind(x, z)`),
# several column names, or a one-sided formula. Every column becomes a slope
# sharing the group's correlated covariance. `intercept = FALSE` drops the
# implicit group intercept, giving a slope-only block (lme4 `(0 + x | g)`).
.tobs_term_re <- function(group, type = c("intercept", "slope", "iid"),
                          covariate = NULL, model = "iid",
                          correlated = TRUE, intercept = TRUE,
                          sigma_scale = 1, id = NULL) {
  type  <- match.arg(type)
  model <- match.arg(model, c("iid", "ar1", "rw1", "rw2"))
  if (type == "slope" && is.null(covariate)) {
    stop("re(): `covariate` must be given for a random slope.", call. = FALSE)
  }
  if (!isTRUE(intercept) && type != "slope") {
    stop("re(): `intercept = FALSE` is only meaningful for a random slope.",
         call. = FALSE)
  }
  # A one-sided formula spells out several slope columns by name; normalise it
  # to the term labels so build_re_spec() resolves them against the data.
  if (inherits(covariate, "formula")) {
    covariate <- attr(stats::terms(covariate), "term.labels")
  }
  codes <- .tobs_index_codes(group, "re", "group")
  .tobs_term(list(
    group_idx = codes$idx, n_groups = codes$n,
    type = type, covariate = covariate, model = model,
    correlated = correlated, intercept = isTRUE(intercept),
    sigma_scale = sigma_scale
  ), class = "tobs_re", id = id, label = "re")
}

# temporal(time, ...)      — AR1 / RW1 / RW2 / IID temporal field
.tobs_term_temporal <- function(time, type = c("ar1", "rw1", "rw2", "iid"),
                                group = NULL, cyclic = FALSE,
                                tau_shape = 1, tau_rate = 0.01, id = NULL) {
  type <- match.arg(type)
  tcodes <- .tobs_index_codes(time, "temporal", "time")
  spec <- list(
    type = type, time_idx = tcodes$idx, n_times = tcodes$n,
    cyclic = cyclic, tau_shape = tau_shape, tau_rate = tau_rate
  )
  if (!is.null(group)) {
    gcodes <- .tobs_index_codes(group, "temporal", "group")
    spec$group_idx <- gcodes$idx
    spec$n_groups  <- gcodes$n
  }
  .tobs_term(spec, class = "tobs_temporal", id = id, label = "temporal")
}

# svc(lon, lat, indices)   — spatially varying coefficients on design columns
.tobs_term_svc <- function(..., indices, coords = NULL, cov = "exponential",
                           nn = 15, id = NULL, sigma2_prior_scale = 1,
                           phi_prior_lower = 0.01, phi_prior_upper = 10) {
  coords <- .tobs_collect_coords(list(...), coords, "svc")
  cov <- match.arg(cov, c("exponential", "matern", "gaussian"))
  n <- nrow(coords); nn <- min(nn, n - 1L)
  nngp <- compute_nngp_neighbors(coords, nn)
  .tobs_term(list(
    indices = as.integer(indices), n_svc = length(indices),
    n_obs = n, nn = nn, coords = as.vector(t(coords)),
    nn_idx = as.vector(t(nngp$nn_idx)), nn_dist = as.vector(t(nngp$nn_dist)),
    nn_order = nngp$nn_order, nn_order_inv = nngp$nn_order_inv,
    cov_type = cov, sigma2_prior_scale = sigma2_prior_scale,
    phi_prior_lower = phi_prior_lower, phi_prior_upper = phi_prior_upper
  ), class = "tobs_svc", id = id, label = "svc")
}

# latent(n_factors, ...)   — latent factors for community models
.tobs_term_latent <- function(n_factors, constraint = 0,
                              sigma_prior_rate = 1, id = NULL) {
  .tobs_term(list(
    n_factors = as.integer(n_factors), constraint = as.integer(constraint),
    sigma_prior_rate = sigma_prior_rate
  ), class = "tobs_latent", id = id, label = "latent")
}

# copy(id)                 — share one realization of a named term across
#                            another process's linear predictor
.tobs_term_copy <- function(id, scale = NULL) {
  if (!is.character(id) || length(id) != 1L) {
    stop("copy(): `id` must be a single string naming a term's `id =`.",
         call. = FALSE)
  }
  .tobs_term(list(ref = id, scale = scale),
             class = "tobs_copy", id = id, label = "copy")
}


# ---------------------------------------------------------------------------
# Spatial umbrella
# ---------------------------------------------------------------------------

# Spatial model names, in registry order. Single source of truth for the
# `spatial()` umbrella's `model =` choices; an areal/continuous term added to
# the registry below is exposed through `spatial()` by listing it here too.
.tobs_spatial_models <- c("icar", "bym2", "car", "car_proper",
                          "gp", "multiscale_gp", "spde")

# spatial(..., model = ...) — single-verb umbrella over the areal (icar / bym2
# / car / car_proper) and continuous (gp / multiscale_gp / spde) spatial terms,
# mirroring temporal()'s one-verb-plus-`type=` surface and INLA's
# f(i, model = ...). Dispatches to the specific constructor via the registry,
# forwarding `...` (bare coord columns, `graph =`, and any per-model arguments)
# and `id` unchanged, so `spatial(graph = adj, model = "bym2")` is identical to
# `bym2(graph = adj)` and `spatial(lon, lat, model = "spde")` to
# `spde(lon, lat)`. Only spatial models dispatch here; re() / temporal() /
# svc() / latent() keep their own verbs.
#
# A varying-coefficient bar (`spatial(~ 1 + w || node, graph = adj, to = ...)`)
# is a second surface: an lme4-style coefficient formula whose intercept column
# is the unweighted shared field and whose covariate columns are weight-scaled
# coefficient fields, all over the graph node index `node`, copied onto the arms
# named by `to`. The bar is expanded later (where the model data is in scope),
# so here it is captured into a `tobs_spatial_bar` spec; see
# .tobs_spatial_bar_spec().
#
# Named arguments are validated against the target constructor's formals before
# dispatch. The areal terms have no `...`, so R already rejects an unknown named
# arg, but the continuous terms (gp / multiscale_gp / spde) take coords through
# `...` and would otherwise *silently* absorb a typo'd or wrong-model named
# argument (`spatial(lon, lat, model = "gp", graph = adj)`) as if it were a
# coordinate. Checking names here closes that gap and names `spatial()`/`model`
# in the error. Positional arguments (the bare coordinate columns) carry no
# name and pass through untouched.
.tobs_term_spatial <- function(..., model = .tobs_spatial_models, id = NULL) {
  dots <- list(...)

  # Varying-coefficient bar form: the first positional argument is a one-sided
  # coefficient formula carrying a `|` / `||` grouping bar. Capture it together
  # with `graph`, `model`, `to` for later (data-aware) expansion; the bar's node
  # index is the areal group_var. The first dot is positional when it is unnamed
  # (no names at all, or an empty-string name in slot 1).
  first_unnamed <- length(dots) >= 1L &&
    (is.null(names(dots)) || !nzchar(names(dots)[[1L]]))
  if (first_unnamed && inherits(dots[[1L]], "formula") &&
      tulpa::tulpa_is_spatial_bar(dots[[1L]])) {
    return(.tobs_spatial_bar_spec(dots[[1L]], dots[-1L], model = model, id = id))
  }

  model <- match.arg(model)
  ctor  <- .tobs_terms[[model]]

  named <- names(list(...))
  named <- named[nzchar(named)]
  if (length(named)) {
    fmls <- setdiff(names(formals(ctor)), "...")
    bad  <- setdiff(named, fmls)
    if (length(bad)) {
      stop(sprintf(
        "spatial(model = \"%s\"): unknown argument%s %s. %s takes: %s.",
        model, if (length(bad) > 1L) "s" else "",
        paste0("`", bad, "`", collapse = ", "), model,
        paste0("`", fmls, "`", collapse = ", ")),
        call. = FALSE)
    }
  }

  ctor(..., id = id)
}

# Arm labels of the cover hurdle: presence (the y > 0 Bernoulli arm) and
# positive (the y | y > 0 arm). `to =` validates against this set; summary() and
# the coefficient output print the same labels (formula label == output label).
.tobs_cover_arms <- c("presence", "positive")

# Capture a `spatial(~ 1 + w || node, graph = adj, to = ...)` varying-coefficient
# bar into a deferred spec. The bar's left-hand side (intercept + covariate
# columns) and node index are stored verbatim; expansion against the model data
# happens in .tobs_expand_spatial_bar(), where each design column becomes either
# the unweighted intercept areal field (intercept column) or a weight-scaled
# trend areal field (covariate column), all on the graph node index. `to` names
# the arms that share this one latent field (presence-anchored, copied to the
# other arm with an estimated coupling alpha); it defaults to all arms so the
# common shared call can omit it. The areal `model` is inherited from the
# umbrella's `model =` (defaults to icar, matching the bare spatial() default).
.tobs_spatial_bar_spec <- function(bar_formula, rest, model = .tobs_spatial_models,
                                   id = NULL) {
  graph <- rest$graph
  to    <- rest$to
  # Only `graph`, `to`, and (optionally) `model` are accepted on the bar form;
  # the node index is the bar RHS, weights are the bar LHS, so group_var/weight
  # named args would be redundant. Reject anything else with a clear pointer.
  known <- c("graph", "to")
  extra <- setdiff(names(rest)[nzchar(names(rest))], known)
  if (length(extra)) {
    stop(sprintf(paste0(
      "spatial(<bar>, ...): unexpected argument%s %s. The bar form takes ",
      "`graph` and `to` only; the node index is the bar right-hand side and ",
      "the per-coefficient weights are the bar left-hand side."),
      if (length(extra) > 1L) "s" else "",
      paste0("`", extra, "`", collapse = ", ")), call. = FALSE)
  }
  if (is.null(graph)) {
    stop("spatial(<bar>, graph = ): a `graph` (adjacency matrix) is required.",
         call. = FALSE)
  }
  .tobs_check_graph(graph, "spatial")

  # `model` may arrive as the full default vector (umbrella default) or a single
  # name; collapse to one areal model. Only areal models carry a graph.
  model <- match.arg(model)
  if (!model %in% c("icar", "bym2", "car", "car_proper")) {
    stop(sprintf(paste0(
      "spatial(<bar>, model = \"%s\"): a varying-coefficient bar needs an areal ",
      "model (icar / bym2 / car / car_proper)."), model), call. = FALSE)
  }

  if (is.null(to)) {
    to <- .tobs_cover_arms
  }
  if (!is.character(to) || length(to) < 1L) {
    stop("spatial(<bar>, to = ): `to` must be a character vector of arm labels.",
         call. = FALSE)
  }
  bad_to <- setdiff(to, .tobs_cover_arms)
  if (length(bad_to)) {
    stop(sprintf(paste0(
      "spatial(<bar>, to = ): unknown arm label%s %s. Valid arms: %s."),
      if (length(bad_to) > 1L) "s" else "",
      paste0("\"", bad_to, "\"", collapse = ", "),
      paste0("\"", .tobs_cover_arms, "\"", collapse = ", ")),
      call. = FALSE)
  }
  to <- unique(to)

  # `|` => correlated (MCAR) fields, `||` => independent fields. Read the bar's
  # top operator off the AST (data-free), so the correlated path can be gated
  # before the model data is in scope.
  correlated <- identical(bar_formula[[2L]][[1L]], as.name("|"))

  .tobs_term(list(
    type        = model,
    is_bar      = TRUE,
    bar_formula = bar_formula,
    correlated  = correlated,
    graph       = graph,
    to          = to
  ), class = "tobs_spatial", id = id, label = model)
}

# Expand a captured varying-coefficient bar spec against the model data into the
# pair the existing weighted-areal-term machinery (gcol33/tulpaObs#59) consumes:
# one unweighted intercept areal field plus, per bar covariate column, a
# weight-scaled trend areal field, all on the same graph keyed by the bar's node
# index (the areal group_var). Each is a plain `tobs_spatial` term identical to
# what `icar(graph, group_var = node)` / `icar(graph, weight = col,
# group_var = node)` would produce, so the bar form desugars to exactly the
# two-term coupled cover path with no engine change. Returns a list of
# `tobs_spatial` terms in column order (intercept first).
# Validate a bar's node-index column against the graph dimension (the bar RHS is
# the graph node index, the old group_var). Shared by the independent expansion
# (.tobs_expand_spatial_bar) and the arm-specific field builder
# (.tobs_armspecific_bar_fields), so one source of truth for the check.
.tobs_validate_bar_node <- function(node, graph, data) {
  if (is.null(data[[node]])) {
    stop(sprintf(paste0(
      "spatial(<bar>): node-index column \"%s\" not found in the data."), node),
      call. = FALSE)
  }
  n_nodes <- nrow(graph)
  lv <- unique(stats::na.omit(data[[node]]))
  if (is.numeric(data[[node]]) && !is.factor(data[[node]])) {
    rng <- range(lv)
    if (rng[1L] < 1L || rng[2L] > n_nodes) {
      stop(sprintf(paste0(
        "spatial(<bar>): node index \"%s\" has values in [%g, %g] but the graph ",
        "has %d node(s); the node index must be a 1..%d code into the graph."),
        node, rng[1L], rng[2L], n_nodes, n_nodes), call. = FALSE)
    }
  } else if (length(lv) > n_nodes) {
    stop(sprintf(paste0(
      "spatial(<bar>): node index \"%s\" has %d distinct levels but the graph ",
      "has %d node(s)."), node, length(lv), n_nodes), call. = FALSE)
  }
  invisible(node)
}

.tobs_expand_spatial_bar <- function(spec, data) {
  specs <- tulpa::tulpa_bar_field_specs(spec$bar_formula, data)
  node  <- attr(specs, "node")
  ctor  <- .tobs_terms[[spec$type]]

  # Validate the node index against the graph dimension (the bar RHS is the
  # graph node index, the old group_var).
  .tobs_validate_bar_node(node, spec$graph, data)

  out <- vector("list", length(specs))
  for (i in seq_along(specs)) {
    col <- specs[[i]]
    if (isTRUE(col$is_intercept)) {
      term <- ctor(graph = spec$graph, group_var = node, id = spec$id)
    } else {
      term <- ctor(graph = spec$graph, group_var = node, id = spec$id)
      term$weight       <- as.numeric(col$weight)
      term$weight_label <- col$column_name
    }
    out[[i]] <- term
  }
  out
}

# Expand a captured arm-specific (single-arm `to`) spatial bar
# (gcol33/tulpaObs#65) against the model data into the per-field design columns
# the cover-hurdle joint driver places on ONE arm with no cross-arm copy. Unlike
# `.tobs_expand_spatial_bar` (which desugars a shared field to the existing
# copied two-block machinery), this keeps the bar as a single self-describing
# spec: the areal `type`, `graph`, the node-index column and per-obs node codes,
# the single target `arm`, and the per-field design weights (the intercept's are
# all-ones, a covariate column's is its per-row value). The fitter builds one
# non-copied areal block per field, restricted to `arm` via a 0-sentinel
# spatial_idx on the other arm. `data_obs` is the NA-dropped data.
.tobs_armspecific_bar_fields <- function(spec, data_obs) {
  arm <- spec$to
  if (length(arm) != 1L) {
    stop("internal: .tobs_armspecific_bar_fields expects a single-arm `to`.",
         call. = FALSE)
  }
  if (!spec$type %in% c("icar", "car", "car_proper")) {
    stop(sprintf(paste0(
      "spatial(<bar>, to = \"%s\"): an arm-specific field uses an intrinsic ",
      "areal model (icar / car / car_proper); model = \"%s\" is not supported ",
      "(the bym2 phi+theta mix is deferred for separate per-arm fields)."),
      arm, spec$type), call. = FALSE)
  }
  specs <- tulpa::tulpa_bar_field_specs(spec$bar_formula, data_obs)
  node  <- attr(specs, "node")
  .tobs_validate_bar_node(node, spec$graph, data_obs)

  n_obs <- nrow(data_obs)
  fields <- lapply(specs, function(col) {
    weight <- if (isTRUE(col$is_intercept) || is.null(col$weight)) {
      rep(1.0, n_obs)
    } else {
      as.numeric(col$weight)
    }
    list(is_intercept = isTRUE(col$is_intercept),
         column_name  = col$column_name,
         weight       = weight)
  })

  list(
    arm     = arm,
    type    = spec$type,
    graph   = spec$graph,
    node    = node,
    idx_obs = as.integer(data_obs[[node]]),
    fields  = fields
  )
}


# ---------------------------------------------------------------------------
# Registry — name -> constructor. Adding a component is one entry here.
# ---------------------------------------------------------------------------
.tobs_terms <- list(
  icar          = .tobs_term_icar,
  bym2          = .tobs_term_bym2,
  car           = .tobs_term_car,
  car_proper    = .tobs_term_car_proper,
  gp            = .tobs_term_gp,
  multiscale_gp = .tobs_term_multiscale_gp,
  spde          = .tobs_term_spde,
  spatial       = .tobs_term_spatial,
  re            = .tobs_term_re,
  temporal      = .tobs_term_temporal,
  svc           = .tobs_term_svc,
  latent        = .tobs_term_latent,
  copy          = .tobs_term_copy
)

# Names of the registered special terms (used by the parser to detect them).
.tobs_term_names <- function() names(.tobs_terms)


# A varying-coefficient areal field -- a per-node SVC weight
# (icar/bym2/car_proper `weight =`), the multi-field intercept + SVC container,
# or a `spatial(~ ... || node)` bar -- is supported on the occu_cover() joint
# spatial path and on the standalone occu() nested-Laplace path. Every other
# areal consumer treats the field as a plain intercept field; rather than
# silently dropping the weighting, error with a pointer to the supported paths.
.tobs_reject_weighted_spatial <- function(spec, context) {
  if (!inherits(spec, "tobs_spatial")) return(invisible(spec))
  is_svc <- !is.null(spec$weight) || isTRUE(spec$is_multifield) ||
            isTRUE(spec$is_bar)
  if (is_svc) {
    stop(sprintf(paste0(
      "%s: a spatially-varying coefficient (a weighted areal term, a multi-",
      "field intercept + trend structure, or a `spatial(~ ... || node)` bar) ",
      "is supported on the occu_cover() joint spatial path and on the ",
      "standalone occu() nested-Laplace path (method = \"nested_laplace\"). ",
      "Drop the varying-coefficient field here."), context),
      call. = FALSE)
  }
  invisible(spec)
}


# Convert an areal `tobs_spatial` term into the `tulpa_spatial` spec the
# cover-hurdle nested-Laplace path consumes. The areal terms retain the raw
# adjacency `graph` (and an optional `group_var` mapping observations to graph
# nodes), so the tulpa-side spec is rebuilt from them. SPDE terms already
# carry a `tulpa_spec`; continuous (gp / multiscale_gp) terms are not areal
# and are not supported by the cover engine. tulpa has no spatial_icar(); an
# intrinsic ICAR is the improper CAR, so icar maps to spatial_car().
.tobs_term_to_tulpa_spatial <- function(spec) {
  .tobs_reject_weighted_spatial(spec, "cover() spatial")
  if (identical(spec$type, "spde")) return(spec$tulpa_spec)
  if (is.null(spec$graph)) {
    stop(sprintf(
      "cover() spatial supports areal terms icar()/bym2()/car()/car_proper() ",
      "(got '%s').", spec$type), call. = FALSE)
  }
  level_args <- if (!is.null(spec$group_var)) {
    list(level = "group", group_var = spec$group_var)
  } else {
    list(level = "obs")
  }
  ctor <- switch(spec$type,
    bym2       = tulpa::spatial_bym2,
    car        = tulpa::spatial_car,
    car_proper = tulpa::spatial_car_proper,
    icar       = tulpa::spatial_car,   # intrinsic CAR
    stop(sprintf("cover() spatial does not support term type '%s'.", spec$type),
         call. = FALSE)
  )
  do.call(ctor, c(list(spec$graph), level_args))
}


# ---------------------------------------------------------------------------
# Print methods for the specs (used when echoing a parsed model)
# ---------------------------------------------------------------------------

#' @export
print.tobs_spatial <- function(x, ...) {
  n <- if (!is.null(x$n_units)) x$n_units else x$n_obs
  cat(sprintf("tobs spatial term: %s (%d units)%s\n", x$type, n,
              if (!is.null(x$id)) sprintf(" [id: %s]", x$id) else ""))
  if (identical(x$type, "spde")) {
    cat(sprintf("  Matern nu=%g, mesh=%d nodes\n", x$nu, x$n_units))
  }
  if (!is.null(x$weight)) {
    cat(sprintf("  Spatially-varying coefficient, weight: %s\n",
                x$weight_label %||% "<numeric>"))
  }
  invisible(x)
}

#' @export
print.tobs_re <- function(x, ...) {
  cat(sprintf("tobs re term: %s (%s, %d groups)%s\n", x$type, x$model,
              x$n_groups, if (!is.null(x$id)) sprintf(" [id: %s]", x$id) else ""))
  if (!is.null(x$covariate)) {
    nm <- if (is.character(x$covariate)) x$covariate
          else if (is.matrix(x$covariate)) colnames(x$covariate)
          else NULL
    label <- if (length(nm)) paste(nm, collapse = ", ")
             else if (is.matrix(x$covariate)) sprintf("%d slopes", ncol(x$covariate))
             else "1 slope"
    cat(sprintf("  Covariate: %s%s\n", label,
                if (isFALSE(x$intercept)) " (no intercept)" else ""))
  }
  invisible(x)
}

#' @export
print.tobs_temporal <- function(x, ...) {
  cat(sprintf("tobs temporal term: %s (%d times)%s\n", x$type, x$n_times,
              if (!is.null(x$id)) sprintf(" [id: %s]", x$id) else ""))
  if (isTRUE(x$cyclic)) cat("  Cyclic: yes\n")
  invisible(x)
}

#' @export
print.tobs_copy <- function(x, ...) {
  cat(sprintf("tobs copy term: -> %s%s\n", x$ref,
              if (!is.null(x$scale)) sprintf(" (scaled)") else ""))
  invisible(x)
}
