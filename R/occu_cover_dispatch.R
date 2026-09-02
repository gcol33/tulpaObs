# ---------------------------------------------------------------------------
# Formula-native cross-arm coupling: share(), the INLA `copy =` idiom
#
# The occurrence arm carries a spatial field; the positive (cover) arm declares
# it carries a scaled copy of that field with a share() selector, the DAG edge
# u_occ -> cover placed in the formula:
#
#   share(spatial(), alpha = grid(g))   the unique occurrence spatial effect,
#                                       one amplitude g over every block
#   share(spatial(cell_idx), ...)       disambiguate by grouping variable when
#                                       the occurrence arm has several spatials
#   share(spatial(), terms = list(intercept = grid(g0), time.sc = grid(g1)))
#                                       a per-block amplitude
#   share("occ_space", ...)             explicit-name reference (lower-level),
#                                       requires spatial(..., name = "occ_space")
#
# No name is needed in the common case: spatial() selects the occurrence arm's
# spatial effect structurally. The engine still reads the coupling amplitude axes
# off `control$alpha.grid` / `control$alpha.grid.trend`, so the formula share() is
# translated into those axes here and the downstream fit is unchanged:
#
#   whole-field amplitude g            -> alpha.grid = g, alpha.grid.trend = g
#   terms = list(intercept=, time.sc=) -> alpha.grid = ., alpha.grid.trend = .
#
# `alpha = grid(g)` integrates over g; a scalar fixes it (a length-1 grid). The
# block layout follows .occu_cover_spatial_fields(): the unweighted intercept
# field is block 1, weighted trend field(s) block 2+. Decoupling an arm is
# structural -- write spatial() (an own field) or omit share() -- not a magic
# alpha of 0; 0 is only ever one value you could place in a grid.
# ---------------------------------------------------------------------------

# Parse share() terms off the positive formula, returning the stripped
# fixed-effects positive formula plus the list of tobs_copy specs. The copy
# special is the only structured term allowed on the positive arm; any other is
# left in place for .occu_cover_reject_structured() to reject.
.occu_cover_extract_pos_copies <- function(pos_formula) {
  if (is.null(pos_formula)) return(list(formula = pos_formula, copies = list()))
  tt   <- stats::terms(pos_formula, keep.order = TRUE)
  labs <- attr(tt, "term.labels")

  # Only the share() calls are evaluated. Parsing the WHOLE formula here
  # evaluates every term against `data = NULL`, so a structured term naming a
  # data column (`re(cell)`, `temporal(year)`) died as "object cell not found"
  # before the arm rejector could name the unsupported feature.
  share_labs <- character(0); keep <- character(0)
  for (lab in labs) {
    e    <- tryCatch(str2lang(lab), error = function(...) NULL)
    head <- if (is.call(e) && is.symbol(e[[1L]])) as.character(e[[1L]])
            else NA_character_
    .tobs_check_retired_term(head)
    if (identical(head, "share")) {
      share_labs <- c(share_labs, lab)
    } else if (!is.na(head) && head %in% c("|", "||")) {
      # terms() strips the parentheses off a bar; reformulate() would rebuild
      # `... + 1 | g`, which re-parses as `(... + 1) | g` and is no longer a bar.
      keep <- c(keep, sprintf("(%s)", lab))
    } else {
      keep <- c(keep, lab)
    }
  }

  copies <- list()
  if (length(share_labs)) {
    sf <- stats::reformulate(termlabels = share_labs, intercept = FALSE)
    environment(sf) <- environment(pos_formula)
    copies <- Filter(function(t) inherits(t, "tobs_copy"),
                     .tobs_parse_formula(sf, data = NULL)$terms)
  }

  fe <- stats::reformulate(
    termlabels = if (length(keep)) keep else "1",
    intercept  = as.logical(attr(tt, "intercept")))
  environment(fe) <- environment(pos_formula)
  list(formula = fe, copies = copies)
}

# Spatial-field constructors that declare a NEW latent field (unlike share(), which
# reuses a named one). A term with one of these heads is a field; placement in an
# arm's formula puts the field on that arm.
.occu_cover_field_ctors <- c("spatial", "icar", "bym2", "car", "car_proper")

# Placement is the canonical way to put a field on an arm: a spatial-field term
# written in the detection or positive formula declares a field ON that arm. Pull
# such terms off their arm formula and carry them as (call, arm) pairs, so the
# arm-generic resolver (.occu_cover_spatial_fields) can evaluate each field spec
# and tag it with its arm alongside the occurrence formula's own fields. share()
# and RE terms are handled separately and left untouched. Returns the unchanged
# occurrence formula, the stripped detection / positive formulas, and the lifted
# arm-field calls.
.occu_cover_lift_arm_fields <- function(occ_formula, det_formula, pos_formula) {
  arm_fields <- list()

  strip_arm <- function(arm_formula, arm) {
    if (is.null(arm_formula)) return(arm_formula)
    tt   <- stats::terms(arm_formula, keep.order = TRUE)
    labs <- attr(tt, "term.labels")
    keep <- character(0)
    for (lab in labs) {
      e    <- tryCatch(str2lang(lab), error = function(...) NULL)
      head <- if (is.call(e) && is.symbol(e[[1L]])) as.character(e[[1L]]) else NA_character_
      .tobs_check_retired_term(head)
      if (!is.na(head) && head %in% .occu_cover_field_ctors) {
        # The arm is fixed by placement; keep the field call unevaluated and
        # carry its arm, to tag the evaluated spec later.
        arm_fields[[length(arm_fields) + 1L]] <<- list(call = e, arm = arm)
      } else if (!is.na(head) && head %in% c("|", "||")) {
        # An lme4 RE bar. terms() strips the parentheses off `(1 | g)` down to the
        # label `1 | g`; reformulate() would rebuild it as `... + 1 | g`, which R
        # re-parses as `(... + 1) | g` -- no longer an RE bar. Restore the parens
        # so the downstream RE parse (.occu_cover_obs_re_parse) still sees a bar.
        keep <- c(keep, sprintf("(%s)", lab))
      } else {
        keep <- c(keep, lab)
      }
    }
    fe <- stats::reformulate(
      termlabels = if (length(keep)) keep else "1",
      intercept  = as.logical(attr(tt, "intercept")))
    environment(fe) <- environment(arm_formula)
    fe
  }

  det2 <- strip_arm(det_formula, "detection")
  pos2 <- strip_arm(pos_formula, "positive")
  list(occ = occ_formula, det = det2, pos = pos2, arm_fields = arm_fields)
}

# Map the positive arm's share() specs onto the coupling-amplitude grids the
# joint fitter reads (control$alpha.grid for the intercept block,
# control$alpha.grid.trend for the trend block). `spatial_info` carries the
# resolved fields (block 1 = intercept, block 2+ = weighted trend), each with a
# `field_name` and a `component` label. Returns the updated control list.
#
# The amplitude axes come ENTIRELY from share(): a block a share() names takes the
# axis it states (or the engine's own axis when it states none), and a block no
# share() names is pinned at alpha = 0 -- the field rides occupancy alone. That
# holds however the occurrence formula spells the field. `name =` is not the
# signal: the selector is type-carrying (share(spatial())), so a field needs no
# name to be copied and carrying one does not couple it.
#
# The one exception is a fit that sets control$alpha.grid[.trend] itself, which
# is the lower-level spelling of the same axes: those grids stand as given, and
# are refused alongside a share().
#
# Everything a copy states about its coefficient is carried here: the amplitude
# axis (stated nodes, or a resolution for the engine's own axis) and the
# hyperprior on the coefficient itself, which lands on control$prior.alpha.
.occu_cover_apply_copy_coupling <- function(copies, spatial_info, control) {
  has_control_alpha <- any(c("alpha.grid", "alpha.grid.trend") %in% names(control))
  if (has_control_alpha && length(copies) > 0L) {
    stop("occu_cover(): set the cross-arm coupling with share() in the positive ",
         "formula OR control$alpha.grid[.trend], not both.", call. = FALSE)
  }
  # `control$alpha.n[.trend]` is a RESOLUTION for the engine's own axis, so it
  # composes with a share() that names no amplitude -- that copy asks for the
  # default axis, and the resolution says how finely to read it. A share() that
  # STATES one (nodes, or grid(n = )) is the same request spelled twice, and is
  # refused.
  .tobs_check_alpha_copy(
    any(vapply(copies, .tobs_copy_states_amplitude, logical(1))),
    control, "occu_cover()")
  # Resolved before the early returns below: a prior written on a share() that
  # selects nothing is a mistake wherever it is written, and the conflict with
  # the control spelling does not depend on there being a field to copy.
  alpha_prior <- .tobs_copy_prior_resolve(copies, control, "occu_cover()")
  # control$alpha.grid is the low-level amplitude knob: when set (and no share())
  # the engine reads the grids as given.
  if (has_control_alpha) return(control)

  if (is.null(spatial_info)) {
    if (length(copies) > 0L) {
      stop("occu_cover(): share() needs a spatial field on the occurrence ",
           "formula, e.g. spatial(~ 1 || cell, graph = adj).", call. = FALSE)
    }
    return(control)
  }

  # Coupling is formula-native and explicit: a share() carries the occurrence
  # spatial field onto the cover arm with the amplitude it names; a block with no
  # share() is decoupled (alpha pinned 0), the field rides occupancy only. There
  # is no implicit default coupling.
  #
  # Component labels of the resolved field blocks. Block 1 is the intercept
  # field; blocks 2+ are weighted trend fields. "trend" is an alias for the
  # single trend block, the column name (e.g. "time.sc") names it explicitly.
  has_trend  <- length(spatial_info$fields) > 1L
  components <- vapply(spatial_info$fields,
                       function(f) f$component %||% NA_character_, character(1))
  node       <- spatial_info$group_var

  # Default to decoupled: every block pinned at alpha = 0. A share() then sets the
  # amplitude axis on the block(s) it names -- to the nodes it states, or to NULL
  # for a share() naming no amplitude, which leaves the block on the engine's own
  # axis (at `control$alpha.n`'s resolution when one is asked for). NULL is also
  # what removes the key below, so a bare share() reaches the fitter exactly as an
  # unset knob does.
  amp_int   <- .tobs_copy_amp(grid = 0)
  amp_trend <- if (has_trend) .tobs_copy_amp(grid = 0) else NULL

  # Apply one (component, amplitude) assignment, returning the canonical block
  # role ("intercept" / "trend") it resolved to. `comp = NULL` is the whole
  # field. `cp_label` names the share() in any error.
  apply_component <- function(comp, amp, cp_label) {
    if (is.null(comp)) {
      amp_int <<- amp
      if (has_trend) amp_trend <<- amp
      return(c("intercept", if (has_trend) "trend"))
    }
    if (identical(comp, "intercept")) {
      amp_int <<- amp
      return("intercept")
    }
    if (identical(comp, "trend") ||
        (has_trend && comp %in% stats::na.omit(components[-1L]))) {
      if (!has_trend) {
        stop(sprintf(paste0(
          "%s: the spatial field has no trend component (it is a single ",
          "intercept field)."), cp_label), call. = FALSE)
      }
      amp_trend <<- amp
      return("trend")
    }
    avail <- paste0("\"", stats::na.omit(components), "\"", collapse = ", ")
    stop(sprintf(paste0(
      "%s: unknown field component \"%s\". Available component(s): %s, or the ",
      "whole field."), cp_label, comp, avail), call. = FALSE)
  }

  for (cp in copies) {
    cp_label <- sprintf("share(%s)", cp$id)

    # The positive arm couples the occurrence spatial field through a selector
    # (spatial() / spatial(<grouping_var>)); a string reference is not a coupling
    # selector here.
    if (is.null(cp$selector_type)) {
      stop(sprintf(paste0(
        "%s: select the occurrence spatial field with share(spatial()) or ",
        "share(spatial(%s)), not a string."), cp_label, node %||% "<grouping_var>"),
        call. = FALSE)
    }
    if (!is.null(cp$selector_group)) {
      if (is.null(node)) {
        stop(sprintf(paste0(
          "%s: the occurrence spatial field has no named grouping variable; ",
          "use share(spatial())."), cp_label), call. = FALSE)
      }
      if (!identical(cp$selector_group, node)) {
        stop(sprintf(paste0(
          "%s: no spatial effect grouped on \"%s\"; the occurrence spatial ",
          "field is on \"%s\". Use share(spatial(%s)) or share(spatial())."),
          cp_label, cp$selector_group, node, node), call. = FALSE)
      }
    }

    # Per-component grids (terms = list(...)) must address every field block, so
    # no block is silently left at alpha = 0. The whole-field / single-component
    # forms set the blocks they name; any unnamed block stays decoupled.
    if (!is.null(cp$copy_terms)) {
      covered <- character(0)
      for (k in names(cp$copy_terms)) {
        covered <- union(covered,
                         apply_component(k, .tobs_copy_amp_of(cp$copy_terms[[k]]),
                                         cp_label))
      }
      required <- c("intercept", if (has_trend) "trend")
      missing_blocks <- setdiff(required, covered)
      if (length(missing_blocks)) {
        blocks <- paste0("\"", stats::na.omit(components), "\"", collapse = ", ")
        stop(sprintf(paste0(
          "%s: terms = must give an amplitude for every field block; %s left ",
          "unaddressed. Field blocks: %s."), cp_label,
          paste0("\"", missing_blocks, "\"", collapse = ", "), blocks),
          call. = FALSE)
      }
    } else {
      apply_component(cp$component,
                      .tobs_copy_amp_of(list(grid = cp$alpha_grid,
                                             n = cp$alpha_n,
                                             integrate = cp$alpha_integrate)),
                      cp_label)
    }
  }

  control <- .tobs_alpha_control_set(control, amp_int, "alpha.grid", "alpha.n")
  if (has_trend) {
    control <- .tobs_alpha_control_set(control, amp_trend,
                                       "alpha.grid.trend", "alpha.n.trend")
  }
  # A decoupled block is the single node 0, and carries no coefficient to
  # regularize; the prior reaches the engine only if some block is copied, and
  # only one block may be, since the engine bakes it on the first
  # (`.tobs_check_alpha_prior()`).
  copied <- c(!.tobs_copy_amp_decoupled(amp_int),
              if (has_trend) !.tobs_copy_amp_decoupled(amp_trend))
  .tobs_check_alpha_prior(alpha_prior, sum(copied), "occu_cover()")
  control[["prior.alpha"]] <- if (sum(copied) == 0L) NULL else alpha_prior
  control
}


# Map the positive arm's share() spec(s) onto the sampled field's coupling
# amplitude on the NUTS spatial path, including the no-copy case, which
# pins the amplitude at 0.
#
# A bare areal term on the psi formula loads on the OCCURRENCE arm alone under
# both engines: a term's process is the formula it sits in, and share() is how a
# field is also placed on the cover arm. Running the translation unconditionally
# is what makes that true here -- with no share() the shared translation pins
# alpha at 0, exactly as the grid-integrated route already did, instead of
# leaving the default amplitude axis in place for the warm fit to estimate and
# the sampler to then integrate over ( made that axis visible work rather than a
# silent default).
#
# `control$alpha.grid` is a live knob here, not a deterministic-backend-only
# one: the sampler's warm nested-Laplace fit integrates that axis, and the axis
# span becomes the support of the sampled alpha's flat prior
# (`.occu_cover_nuts_hyper_bounds()`), so `alpha = grid(c(...))` bounds the
# amplitude and a scalar `alpha =` collapses the axis to one node and pins it.
# The same share() therefore means the same thing under both engines.
#
# The translation is the shared `.occu_cover_apply_copy_coupling()`, so the
# selector rules, the per-component addressing and every error message are one
# implementation. It reads two things off the resolved field description -- the
# per-block roster (one entry per block, each with its `component` role) and the
# grouping variable a `share(spatial(<var>))` selector must match. The sampled
# path resolves that roster itself (`.occu_cover_spatial_fields()` cannot supply
# it: it resolves the grid-integrated path's fields and rejects `car_proper()`,
# the sampled field's primary kind), so the roster here is the intercept field
# plus one entry per varying-coefficient field, and a per-component
# `share(spatial(), terms = list(...))` addresses them by the same names it does
# under nested_laplace.
.occu_cover_nuts_copy_control <- function(copies, nuts_sp, control) {
  field_view <- list(
    fields    = c(list(list(component = "intercept")),
                  lapply(nuts_sp$trends %||% list(),
                         function(tf) list(component = tf$label))),
    group_var = nuts_sp$group_var)
  .occu_cover_apply_copy_coupling(copies, field_view, control)
}


# ---------------------------------------------------------------------------
# Dispatcher (wired into tobs.R's switch)
# ---------------------------------------------------------------------------

.dispatch_occu_cover <- function(formula, data, family, detection, y, visits,
                                  engine, priors, control,
                                  approx = "gaussian_laplace",
                                  correction = "none", ...) {
  dots <- list(...)

  if (is.null(detection)) {
    stop("occu_cover() requires a `detection` formula.", call. = FALSE)
  }
  if (is.null(y)) {
    stop("occu_cover() requires `y` (N x J detection-history matrix).",
         call. = FALSE)
  }
  if (is.null(dots$y_pos)) {
    stop("occu_cover() requires `y_pos` (N x J positive-cover matrix; ",
         "values used only where y == 1).", call. = FALSE)
  }
  # `alpha.grid` states the copy axis's nodes, `alpha.n` how many nodes the
  # engine's own axis is read at; one block takes one of them.
  .tobs_check_alpha_control(control, "occu_cover()")

  # Compact (ragged) input: `y` is a tobs_ragged carrier (tobs_data(compact =
  # TRUE)). It feeds the joint nested-Laplace engine one valid visit at a time,
  # with no padded grid and so no per-site visit cap. `y_pos` may be the paired
  # ragged carrier or a bare length-V vector aligned to `y`. Scoped to the
  # spatial joint path (the only consumer of the compacted arms); every other
  # configuration that would read the dense `y` / `valid` grid is gated off below
  # with a clear error rather than a silent dense rebuild.
  ragged <- inherits(y, "tobs_ragged")
  if (ragged) {
    y_pos_arg <- dots$y_pos
    if (inherits(y_pos_arg, "tobs_ragged")) {
      if (!identical(y_pos_arg$site, y$site) || !identical(y_pos_arg$visit, y$visit))
        stop("occu_cover(): compact `y` and `y_pos` are not aligned (different ",
             "site / visit order). Build both with tobs_data(compact = TRUE) on ",
             "the same df / site / visit.", call. = FALSE)
      y_pos_values <- y_pos_arg$values
    } else {
      y_pos_values <- y_pos_arg
    }
  }

  pos_formula <- dots$positive

  # Placement -> arm: a spatial-field term written in the positive (or detection)
  # formula is pulled off that arm formula and carried as a (call, arm) pair, so
  # the arm-generic spatial resolver evaluates and tags every arm's fields
  # together. Done before the RE parse and the positive-defaults-to-detection
  # fallback, so each downstream step sees a field-free arm formula. share() stays
  # on its arm's formula (it is a reference, not a new field). An arm is chosen by
  # placement (write the field in that arm's formula); a field is shared across
  # arms with share().
  lifted         <- .occu_cover_lift_arm_fields(formula, detection, pos_formula)
  formula        <- lifted$occ
  detection      <- lifted$det
  pos_formula    <- lifted$pos
  arm_fields     <- lifted$arm_fields

  if (is.null(pos_formula)) pos_formula <- detection

  # Observation-arm random intercept: a `(1 | g)` / `re(g)` on the detection or
  # positive-cover formula adds a random effect on that arm -- an iid latent
  # block on the joint nested-Laplace fit, a non-centered block with its own
  # sampled SD under method = "nuts". Parse it off FIRST -- before the share()
  # extraction and every design build -- so each downstream consumer sees a
  # clean fixed-effects formula; the grouping is resolved to per-visit codes
  # once the model (and its `valid` mask) is built. The plain `laplace` route is
  # a coefficient-marginal fit with no latent block at all, so it errors rather
  # than silently dropping the RE.
  det_re_parse <- .occu_cover_obs_re_parse(detection,   "detection")
  pos_re_parse <- .occu_cover_obs_re_parse(pos_formula, "positive cover")
  has_obs_re   <- !is.null(det_re_parse) || !is.null(pos_re_parse)
  if (has_obs_re && !engine %in% c("nested_laplace", "nuts")) {
    stop("occu_cover(): a random effect on the detection / positive-cover arm ",
         "needs method = \"nested_laplace\" (the joint nested-Laplace engine ",
         "carries the RE as a latent block) or method = \"nuts\" (which samples ",
         "it, group SD included); got method = \"", engine, "\".",
         call. = FALSE)
  }
  if (ragged && identical(engine, "nuts")) {
    stop("occu_cover(): compact (ragged) input feeds the joint nested-Laplace ",
         "engine; method = \"nuts\" reads the dense detection grid. Use ",
         "method = \"nested_laplace\", or build the data densely for NUTS.",
         call. = FALSE)
  }
  # Random intercepts (crossed, nested) AND random slopes are supported on the
  # detection / positive-cover arms: an intercept rides one `iid` block per term,
  # an uncorrelated slope one weighted `iid` block per coefficient, a correlated
  # slope one multivariate free-Sigma `miid` block. The slope blocks need tulpa's
  # joint engine >= 0.0.39, which the DESCRIPTION Imports floor enforces.
  if (!is.null(det_re_parse)) detection   <- det_re_parse$fe
  if (!is.null(pos_re_parse)) pos_formula <- pos_re_parse$fe

  # Formula-native cross-arm coupling: pull any share() term off the positive
  # formula and keep the stripped fixed-effects design. The share() specs are
  # translated into the engine's coupling-amplitude grids once the occupancy
  # field blocks are resolved (below). Stripping here lets every downstream
  # consumer (NUTS branch, design build, structured-term rejection) see a clean
  # fixed-effects positive formula.
  pos_copy   <- .occu_cover_extract_pos_copies(pos_formula)
  pos_copies <- pos_copy$copies
  pos_formula <- pos_copy$formula

  # The engine gate on share(residual = ) is read HERE, before the sampled-field
  # branch returns: the deviation is a grid-integrated latent block, so a fit
  # asking for one under another engine has to be told, not quietly given a fit
  # without it. The remaining gates need the resolved field roster and run below.
  .occu_cover_residual_check(.occu_cover_residual_spec(pos_copies),
                             list(), engine)

  # Spatial NUTS path: a single areal term on the psi formula -- icar() / bym2() /
  # car_proper(), or the equivalent single-column bar spatial(~ 1 || cell, graph =
  # adj) -- under method = "nuts" samples a non-centered coupled field, and its
  # hyperparameters, jointly with the coefficient marginal (rather than
  # grid-integrating them, as nested_laplace does). Detected separately from the
  # grid-integrated fields below, because the sampled field is a NUTS-only
  # structure.
  if (identical(engine, "nuts")) {
    nuts_sp <- .occu_cover_nuts_spatial_term(formula, data)
    if (!is.null(nuts_sp)) {
      if (has_obs_re) {
        stop("occu_cover(): method = \"nuts\" samples an observation-arm random ",
             "effect on the NON-SPATIAL target; composed with the coupled areal ",
             "field it needs the grid-integrated method = \"nested_laplace\" ",
             "(which carries both as latent blocks).",
             call. = FALSE)
      }
      .occu_cover_reject_structured(detection,   "detection")
      .occu_cover_reject_structured(pos_formula, "positive cover")
      # The share() specs were stripped off the positive formula above, so they
      # have to be translated here rather than at the shared call below this
      # branch. Ahead of the design build, so a share() the field cannot satisfy
      # is refused before any fitting work.
      control <- .occu_cover_nuts_copy_control(pos_copies, nuts_sp, control)
      vd_det  <- .normalize_visits(visits, detection,
                                   n_sites = nrow(y), max_visits = ncol(y))
      vd_pos  <- .normalize_visits(visits, pos_formula,
                                   n_sites = nrow(y), max_visits = ncol(y))
      model_sp <- .tobs_build_occu_cover(
        occ_formula = nuts_sp$fe, det_formula = vd_det$det_formula,
        pos_formula = vd_pos$det_formula, data = data, y = y,
        y_pos = dots$y_pos, positive = family$params$positive,
        det_visit_formula = vd_det$det_visit_formula,
        det_visit_data    = vd_det$visits,
        pos_visit_formula = vd_pos$det_visit_formula,
        pos_visit_data    = vd_pos$visits)
      model_sp$cover_aggregate <- "none"
      # Resolve the site -> field-node map (group_var lets sites > cells).
      sp_graph <- nuts_sp$spatial$graph
      n_cells_f <- nrow(sp_graph)
      gv <- nuts_sp$group_var
      if (!is.null(gv)) {
        if (!gv %in% names(data))
          stop(sprintf("occu_cover() group_var '%s' is not a column of data.",
                       gv), call. = FALSE)
        site_cell <- as.integer(data[[gv]])
        if (length(site_cell) != model_sp$n_sites || anyNA(site_cell) ||
            min(site_cell) < 1L || max(site_cell) > n_cells_f)
          stop(sprintf(paste0(
            "occu_cover() group_var '%s' must be an integer cell index in ",
            "1..%d, one per site (%d sites)."), gv, n_cells_f, model_sp$n_sites),
            call. = FALSE)
      } else {
        if (model_sp$n_sites != n_cells_f)
          stop(sprintf(paste0(
            "occu_cover() NUTS spatial: %d sites but the graph has %d nodes. ",
            "Map sites to cells with group_var on the car_proper() term, or ",
            "match the site count to the graph."),
            model_sp$n_sites, n_cells_f), call. = FALSE)
        site_cell <- seq_len(model_sp$n_sites)
      }
      model_sp$site_cell <- site_cell
      model_sp$n_cells   <- n_cells_f
      return(do.call(.tobs_fit_occu_cover_nuts_spatial,
                     c(list(model = model_sp, spatial = nuts_sp$spatial,
                            trends = nuts_sp$trends %||% list(),
                            priors = priors), control)))
    }
  }

  # Detect the coupled spatial field(s) on the psi formula. The spatial path is
  # the joint nested-Laplace engine (shared field(s) across the psi and cover
  # arms); the non-spatial path is plain Laplace on the exact two-state
  # marginal. A weighted areal term adds a second coupled (SVC) field.
  spatial_info <- .occu_cover_spatial_fields(formula, data, arm_fields)
  has_spatial  <- !is.null(spatial_info)

  # Translate the positive arm's share() spec(s) into the engine coupling grids
  # now that the occupancy field blocks are resolved. Every fit carrying a
  # spatial field gets control$alpha.grid[.trend] set here -- to the amplitude
  # each share() names, and to 0 (decoupled) for a block none names -- unless the
  # fit stated those grids itself, in which case they stand. A share() with no
  # spatial field to select is an error surfaced just above.
  if (length(pos_copies) > 0L && !has_spatial) {
    stop("occu_cover(): share() on the positive arm needs a spatial field on the ",
         "occurrence formula (e.g. spatial(~ 1 || cell, graph = adj, name = ",
         "\"occ_space\")) under method = \"nested_laplace\".", call. = FALSE)
  }
  control <- .occu_cover_apply_copy_coupling(pos_copies, spatial_info, control)

  # A share() may also ask for an arm-specific DEVIATION beside the copy
  # (`residual =`), which the cover arm carries as its own latent block. Resolved
  # here, next to the amplitude it belongs with, and gated against the
  # configurations the arm-specific block does not carry; the block itself is
  # built below, once the site -> node map is known.
  residual_spec <- .occu_cover_residual_spec(pos_copies)
  if (!is.null(residual_spec)) {
    if (!has_spatial) {
      stop("occu_cover(): share(residual = ) is a deviation from the shared ",
           "field, so it needs one on the occurrence formula (e.g. ",
           "spatial(~ 1 || cell, graph = adj)).", call. = FALSE)
    }
    .occu_cover_residual_check(residual_spec, spatial_info, engine)
  }

  # Resolve cover aggregation. NULL (unset) -> "mean" on the shared-field spatial
  # path (so the cover arm contributes at the cell scale and does not outweigh
  # occupancy on the shared field), "none" (per-visit) on the non-spatial path
  # (no shared field to over-weight). `agg_explicit` records whether the user set
  # it: an explicit mean / median on an unsupported configuration errors, whereas
  # the bare default quietly falls back to per-visit cover so a plain visit-level
  # fit keeps working.
  #
  # Aggregated cover is a per-cell observation, so its positive design must be
  # cell-level (resolved from the cell `data`). A `positive` formula that
  # references a visit-level covariate (a name carried in `visits`) is a
  # per-visit design and cannot be aggregated: an explicit request errors, the
  # bare default falls back to per-visit cover.
  visit_cov_names <- if (is.null(visits)) character(0)
                     else if (is.data.frame(visits) || is.list(visits)) names(visits)
                     else character(0)
  pos_is_visit_level <- length(intersect(all.vars(pos_formula),
                                          visit_cov_names)) > 0L

  agg_explicit    <- !is.null(family$params$cover_aggregate)
  cover_aggregate <- family$params$cover_aggregate %||%
                     (if (has_spatial) "mean" else "none")
  if (!has_spatial && cover_aggregate != "none") {
    stop(sprintf(paste0(
      "occu_cover(cover_aggregate = \"%s\") aggregates the cover arm on the ",
      "shared-field spatial path (method = \"nested_laplace\"); the non-spatial ",
      "laplace fit uses per-visit cover (cover_aggregate = \"none\")."),
      cover_aggregate), call. = FALSE)
  }
  if (cover_aggregate != "none" && pos_is_visit_level) {
    if (agg_explicit) {
      stop(sprintf(paste0(
        "occu_cover() cell-aggregated cover (cover_aggregate = \"%s\") needs a ",
        "cell-level positive design, but the `positive` formula references the ",
        "visit-level covariate(s) %s (carried in `visits`). Use a cell-level ",
        "positive covariate (a column of `data`), or cover_aggregate = \"none\" ",
        "for per-visit cover."), cover_aggregate,
        paste(intersect(all.vars(pos_formula), visit_cov_names),
              collapse = ", ")), call. = FALSE)
    }
    cover_aggregate <- "none"
  }
  # An arm-specific cover field is scored per detected visit (its node/weight
  # index the pos-arm visit rows), so it needs per-visit cover; an explicit
  # aggregation errors, the bare default falls back to per-visit.
  if (has_spatial &&
      (!is.null(spatial_info$pos_armspec) || !is.null(residual_spec)) &&
      cover_aggregate != "none") {
    if (agg_explicit) {
      stop(sprintf(paste0(
        "occu_cover(): an arm-specific cover field (to = \"positive\" / ",
        "share(residual = )) uses ",
        "per-visit cover (cover_aggregate = \"none\"); it cannot map onto ",
        "cell-aggregated cover rows. Got cover_aggregate = \"%s\"."),
        cover_aggregate), call. = FALSE)
    }
    cover_aggregate <- "none"
  }

  if (has_spatial && engine == "laplace") {
    stop("occu_cover() found a spatial term (icar/bym2) in the psi formula ",
         "but method = \"laplace\" is non-spatial. Use method = ",
         "\"nested_laplace\" for the spatial path.", call. = FALSE)
  }
  if (has_spatial && engine == "nuts") {
    # A car_proper() term would already have routed to the spatial NUTS fitter
    # above; reaching here means an intrinsic icar()/bym2() field, whose flat
    # field-mean direction needs the grid-integrated nested-Laplace path (or the
    # sampled-field community route) -- it is not the fixed-hyper NUTS structure.
    stop("occu_cover() with method = \"nuts\" samples a FIXED-HYPER proper-CAR ",
         "shared field; the intrinsic icar()/bym2() field on the psi formula ",
         "has a flat field-mean direction needing a sum-to-zero reparameterisation ",
         "for NUTS -- use method = \"nested_laplace\" (the shared coupled field is ",
         "grid-integrated), car_proper() for the NUTS shared field, or ",
         "ms_occu_cover() + icar() for a sampled shared field.",
         call. = FALSE)
  }
  if (!has_spatial && engine == "nested_laplace") {
    stop("occu_cover() with method = \"nested_laplace\" requires a spatial ",
         "term (icar() or bym2()) on the psi formula.", call. = FALSE)
  }
  if (ragged && !has_spatial) {
    stop("occu_cover(): compact (ragged) input is implemented for the joint ",
         "nested-Laplace path (a spatial term on the occurrence formula). For a ",
         "non-spatial fit, build the data densely (tobs_data() without ",
         "compact = TRUE).", call. = FALSE)
  }
  if (ragged && cover_aggregate != "none") {
    stop("occu_cover(): compact (ragged) input uses per-visit cover ",
         "(cover_aggregate = \"none\"); the cell-aggregated cover path reads the ",
         "dense detection grid. Pass cover_aggregate = \"none\", or build densely.",
         call. = FALSE)
  }

  fe_formula <- if (has_spatial) spatial_info$fe else formula

  # Detection / cover arms never carry a spatial term (the shared field is
  # on the latent state z, not on the observation process). Other structured
  # terms (re, temporal, ...) are not supported on any arm in v1/v2.
  .occu_cover_reject_structured(detection,   "detection")
  .occu_cover_reject_structured(pos_formula, "positive cover")

  # Visit-design normalization + model build. The compact (ragged) path builds
  # the visit designs on the V valid rows directly (no padded grid); the dense
  # path flattens the [n_sites x max_visits] grid. Both converge to the same
  # `tobs_model` shape (the compact one carries `ragged = TRUE` and pre-compacted
  # visit-level fields), so the joint-coupled arm builder downstream is the same.
  if (ragged) {
    n_v    <- y$n_visits
    vd_det <- .normalize_visits_ragged(visits, detection, n_visits_valid = n_v)
    vd_pos <- .normalize_visits_ragged(visits, pos_formula, n_visits_valid = n_v)
    model  <- .tobs_build_occu_cover_ragged(
      occ_formula       = fe_formula,
      det_formula       = vd_det$det_formula,
      pos_formula       = vd_pos$det_formula,
      data              = data,
      y_ragged          = y,
      y_pos_values      = y_pos_values,
      positive          = family$params$positive,
      det_visit_formula = vd_det$det_visit_formula,
      det_visit_data    = vd_det$visits,
      pos_visit_formula = vd_pos$det_visit_formula,
      pos_visit_data    = vd_pos$visits)
  } else {
    vd_det <- .normalize_visits(visits, detection,
                                n_sites = nrow(y), max_visits = ncol(y))
    # Positive design. Per-visit cover reads the visit-level positive formula from
    # `visits`; cell-aggregated cover reads a cell-level positive design directly
    # from `data` (one value per occupancy unit) and carries no visit-level term.
    if (cover_aggregate == "none") {
      vd_pos            <- .normalize_visits(visits, pos_formula,
                                             n_sites = nrow(y), max_visits = ncol(y))
      pos_site_formula  <- vd_pos$det_formula
      pos_visit_formula <- vd_pos$det_visit_formula
      pos_visit_data    <- vd_pos$visits
    } else {
      # Cell-aggregated cover: the positive design is cell-level, so there is no
      # per-visit positive frame (and no per-visit cover RE, gated below).
      vd_pos            <- NULL
      pos_site_formula  <- pos_formula
      pos_visit_formula <- NULL
      pos_visit_data    <- NULL
    }

    model <- .tobs_build_occu_cover(
      occ_formula      = fe_formula,
      det_formula      = vd_det$det_formula,
      pos_formula      = pos_site_formula,
      data             = data,
      y                = y,
      y_pos            = dots$y_pos,
      positive         = family$params$positive,
      det_visit_formula = vd_det$det_visit_formula,
      det_visit_data    = vd_det$visits,
      pos_visit_formula = pos_visit_formula,
      pos_visit_data    = pos_visit_data
    )
  }

  model$cover_aggregate <- cover_aggregate

  # Observation-arm random effects: resolve each detection / positive-cover RE
  # term to its per-(site, visit) design (group codes + slope weights) now that
  # the model carries its `valid` mask. Both hosting engines read model$re_det /
  # model$re_pos from here -- the joint nested-Laplace builder subsets each term's
  # codes by the same `keep` it uses for the arm rows, the sampler scatters them
  # on the padded grid.
  model <- .occu_cover_attach_obs_re(model, det_re_parse, pos_re_parse, data,
                                     vd_det$visits, vd_pos$visits)

  if (has_spatial) {
    fields      <- spatial_info$fields
    base_graph  <- fields[[1L]]$graph

    # Resolve the site -> field-node map. With group_var the occupancy units
    # (sites, one per row of `data` / `y`) map onto fewer field nodes (cells),
    # so the same cell field is shared across that cell's sites (e.g. cell-year
    # sites sharing one cell). Without group_var the two coincide 1:1.
    n_cells_field <- nrow(base_graph)
    gv <- spatial_info$group_var
    if (!is.null(gv)) {
      if (!gv %in% names(data)) {
        stop(sprintf("occu_cover() group_var '%s' is not a column of data.", gv),
             call. = FALSE)
      }
      site_cell <- as.integer(data[[gv]])
      if (length(site_cell) != model$n_sites || anyNA(site_cell) ||
          min(site_cell) < 1L || max(site_cell) > n_cells_field) {
        stop(sprintf(paste0(
          "occu_cover() group_var '%s' must be an integer cell index in 1..%d, ",
          "one per site (%d sites)."), gv, n_cells_field, model$n_sites),
          call. = FALSE)
      }
    } else {
      if (model$n_sites != n_cells_field) {
        stop(sprintf(paste0(
          "occu_cover() spatial: %d sites but the graph has %d nodes. Map sites ",
          "to cells with group_var = \"<col>\" on the icar()/bym2() term (e.g. ",
          "site = cell-year), or match the site count to the graph."),
          model$n_sites, n_cells_field), call. = FALSE)
      }
      site_cell <- seq_len(model$n_sites)
    }
    model$site_cell <- site_cell
    model$n_cells   <- n_cells_field

    # share(residual = "full"): the deviation spans the node set, which is what
    # the arm-specific block already is, so it is built here -- where the site ->
    # node map is resolved -- and handed to the fitter through the same slot a
    # placed positive-arm field fills. One block builder, two spellings. A
    # rank-r deviation is a different block and is prepared at the fit call.
    if (!is.null(residual_spec) && identical(residual_spec$rank, "full")) {
      spatial_info$armspec[["pos"]] <- .occu_cover_residual_armspec(
        residual_spec, site_cell, model$n_sites)
      spatial_info$pos_armspec <- spatial_info$armspec[["pos"]]
    }

    # Optional per-group random intercept on the occupancy arm, layered on the
    # shared field. The grouping is per occupancy unit (one code per site / data
    # row); validate its length and carry it to the fitter. It also lands on the
    # model, so the per-site offset the criteria add to the occupancy predictor
    # is built from the same codes the fit ran on -- the occupancy counterpart of
    # model$re_det / model$re_pos.
    re_spec <- spatial_info$re
    if (!is.null(re_spec)) {
      if (length(re_spec$group_idx) != model$n_sites) {
        stop(sprintf(paste0(
          "occu_cover() spatial + RE: the random-effect grouping has %d codes ",
          "but there are %d occupancy units (sites)."),
          length(re_spec$group_idx), model$n_sites), call. = FALSE)
      }
      model$re_psi <- re_spec
    }

    # joint (3-arm nested-Laplace via tulpa's cell_coupling spec) is the
    # default: outer-grid integration over (sigma, alpha [, sigma_trend,
    # alpha_trend]) with inner Newton driven by the occu_cover_{lognormal,beta}
    # cell-coupling spec. 150-300x faster than v3 at N=100 and reliably completes
    # at N=200+ where v3 trips on a missing-value compare in its outer BFGS. v3
    # pure-R nested-Laplace and v2's joint Laplace stay reachable via
    # control$engine = "v3_nested" / "v2_joint" as debug escape hatches; both
    # take only the single intercept field.
    correlated <- isTRUE(spatial_info$correlated)
    # Default 3-arm nested-Laplace fitter (coupling lives in the positive
    # formula via share()); "v2_joint" / "v3_nested" are the single-field escape
    # hatches handled below.
    engine_pick <- control[["engine"]] %||% "joint"
    control[["engine"]] <- NULL
    if (correlated && engine_pick %in% c("v2_joint", "v3_nested")) {
      stop(sprintf(paste0(
        "occu_cover(): a correlated spatial bar (`|`, free-Sigma MCAR) needs ",
        "the default joint engine; the \"%s\" escape hatch couples a ",
        "single shared field only."), engine_pick), call. = FALSE)
    }
    if (engine_pick %in% c("v2_joint", "v3_nested")) {
      if (length(spatial_info$armspec)) {
        stop(sprintf(paste0(
          "occu_cover() an arm-specific field (to = \"positive\" / ",
          "\"detection\") needs the default joint engine; the \"%s\" ",
          "escape hatch couples a single shared field only."), engine_pick),
          call. = FALSE)
      }
      if (!is.null(re_spec) || !is.null(model$re_det) || !is.null(model$re_pos)) {
        stop(sprintf(paste0(
          "occu_cover() per-group RE needs the default joint engine; ",
          "the \"%s\" escape hatch has no RE block."),
          engine_pick), call. = FALSE)
      }
      # The v2/v3 escape hatches model per-visit cover only; cell-aggregated
      # cover is a joint feature. An explicit request errors; the bare
      # default falls back to per-visit on these engines.
      if (cover_aggregate != "none") {
        if (agg_explicit) {
          stop(sprintf(paste0(
            "occu_cover() cell-aggregated cover (cover_aggregate = \"%s\") is ",
            "wired on the default joint engine; the \"%s\" escape hatch ",
            "models per-visit cover only."), cover_aggregate, engine_pick),
            call. = FALSE)
        }
        model$cover_aggregate <- "none"
      }
      if (length(fields) > 1L) {
        stop(sprintf(paste0(
          "occu_cover() engine \"%s\" couples a single shared field; ",
          "weighted SVC field(s) need the default joint engine."),
          engine_pick), call. = FALSE)
      }
      if (!is.null(gv)) {
        stop(sprintf(paste0(
          "occu_cover() engine \"%s\" binds the field 1:1 to sites and does ",
          "not support group_var; use the default joint engine."),
          engine_pick), call. = FALSE)
      }
      fit_args <- c(list(model = model, adj = base_graph, priors = priors),
                    control)
      fitter <- if (engine_pick == "v2_joint") .tobs_fit_occu_cover_spatial
                else .tobs_fit_occu_cover_nested
      return(do.call(fitter, fit_args))
    }
    fit_args <- c(list(model = model, fields = fields, priors = priors,
                       re_spec = re_spec, correlated = correlated,
                       pos_armspec = spatial_info$armspec[["pos"]],
                       det_armspec = spatial_info$armspec[["p"]]),
                  control)
    # A rank-r deviation is anchored on a warm fit of the model without it: that
    # fit supplies the shared field the basis is orthogonalized against and the
    # scale the basis coefficients are pinned at. Prepared from the SAME
    # `fit_args`, so the warm fit and the fit it warms differ in the deviation
    # and in nothing else.
    if (!is.null(residual_spec)) {
      prep <- .occu_cover_residual_prepare(residual_spec, fit_args, base_graph)
      if (!is.null(prep)) {
        fit_args$pos_residual <- list(basis = prep$basis, sigma = prep$sigma)
        residual_spec$basis   <- prep$basis
        residual_spec$sigma   <- prep$sigma
      }
    }
    return(.occu_cover_residual_attach(
      do.call(.tobs_fit_occu_cover_joint, fit_args), residual_spec))
  }

  # Non-spatial NUTS: sample the exact two-state coefficient marginal (the
  # in-tree FullGradFn), warm-started at the Laplace mode. Other non-spatial
  # routes (only "laplace" here) fit the direct Laplace optim.
  if (identical(engine, "nuts")) {
    return(do.call(.tobs_fit_occu_cover_nuts,
                   c(list(model = model, priors = priors), control)))
  }

  fit_args <- list(model = model, method = engine, priors = priors)
  fit_args <- c(fit_args, control)
  do.call(.tobs_fit_occu_cover, fit_args)
}


